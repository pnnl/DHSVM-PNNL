/*
 * SUMMARY:      RouteSubSurface.c - Route subsurface flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bill Perkins
 * ORG:          PNNL
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Route subsurface flow
 * DESCRIP-END.
 * FUNCTIONS:    RouteSubSurface()
 * COMMENTS:
 * $Id: RouteSubSurface.c,v3.1.2 2013/08/18 ning Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "soilmoisture.h"
#include "slopeaspect.h"
#include "DHSVMChannel.h"

#ifndef MIN_GRAD
#define MIN_GRAD .3		/* minimum slope for flow to channel */
#endif


/*****************************************************************************
  RouteSubSurface()

  Sources: 
  Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed 
      hydrology-vegetation model for complex terrain, Water Resour. Res.,
      30(6), 1665-1679, 1994.

  Quinn, P., K. Beven, P. Chevallier, and O. Planchon, The prediction of 
      hillslope flow paths for distributed hydrological modelling using 
      digital terrain models, Hydrological Processes, 5, 59-79, 1991.

  This routine follows Wigmosta et al. [1994] in calculating the subsurface
  flow.  The local gradient is based on the local hydraulic head, consisting 
  of the height of the pixel surface minus the depth of the water table 
  below the water surface.  This has the disadvantage that the local gradients
  have to be recalculated for each pixel for each timestep.  In Wigmosta et
  al. [1994] the local gradient is taken to be equal to the slope of the land
  surface, which is a reasonable assunption for mountainous areas.  For the 
  flat boreal forest landscape it is probably better to use the slope
  of the water surface.

  Set the gradient with pixels that are outside tha basin to zero.  This 
  ensures that there is no flux of water across the basin boundary.  In the 
  current implementation water can only leave the basin as surface flow.  
  This may not be entirely realistic, and should be analyzed further.  
  One consequence of this could be that the soil in the basin is more 
  saturated than it would be if subsurface flow out of the basin would
  be allowed.

  The surrounding grid cells are numbered in the following way

                |-----| DX

          0-----1-----2  ---
	  |\    |    /|   |
          | \   |   / |   |
          |  \  |  /  |   | DY
          |   \ | /   |   |
          |    \|/    |   |
          7-----*-----3  ---
          |    /|\    |
          |   / | \   |
          |  /  |  \  |
          | /   |   \ |
          |/    |    \|
          6-----5-----4

  For the current implementation it is assumed that the resolution is the 
  same in both the X and the Y direction.  If this is not the case an error
  message is generated and the program exits.  The reason is that the 
  formulation for the flow width in the diagonal direction changes if the
  grid is not square.  The method for defining the flow widths in the case
  of square grids is taken from Quinn et al [1991]

  Update Jan 2004 COD
  When Gradient = WATERTABLE, the watertable was used to route the
  surface water. This was because of the common use of TopoMap.Dir and 
  TopoMap.TotalDir. These are now for surface routing (always) and subsurface 
  routing (when Gradient = TOPOGRAPHY). Subsurface routing directions 
  and FlowGrad (SubDir, SubTotalDir, SubFlowGrad) for Gradient = WATERTABLE 
  are now determined locally here (in RouteSubsurface.c.)

  WORK IN PROGRESS
*****************************************************************************/
void RouteSubSurface(int Dt, MAPSIZE *Map, TOPOPIX **TopoMap,
		     VEGTABLE *VType, VEGPIX **VegMap,
		     ROADSTRUCT **Network, SOILTABLE *SType,
		     SOILPIX **SoilMap, CHANNEL *ChannelData,
		     TIMESTRUCT *Time, OPTIONSTRUCT *Options, 
		     char *DumpPath, int MaxStreamID, SNOWPIX **SnowMap)
{
  const char *Routine = "RouteSubSurface";
  int x;			/* counter */
  int y;			/* counter */
  int i, j;	        /* counters */
  float BankHeight;
  float *Adjust;
  float fract_used;
  float depth;
  float OutFlow;
  float water_out_road;
  float Transmissivity;
  float AvailableWater;
  int k;
  float **SubFlowGrad;	        /* Magnitude of subsurface flow gradient slope * width */
  unsigned char ***SubDir;      /* Fraction of flux moving in each direction*/ 
  unsigned int **SubTotalDir;	/* Sum of Dir array */

  int count, totalcount;
  float mgrid, sat;
  char buffer[32];
  char satoutfile[100];         /* Character arrays to hold file name. */ 
  FILE *fs;                     /* File pointer. */

  /*****************************************************************************
   Allocate memory 
  ****************************************************************************/
  
  if (!(SubFlowGrad = (float **)calloc(Map->NY, sizeof(float *))))
    ReportError((char *) Routine, 1);
  for(i=0; i<Map->NY; i++) {
    if (!(SubFlowGrad[i] = (float *)calloc(Map->NX, sizeof(float))))
      ReportError((char *) Routine, 1);
  }
  
  if (!((SubDir) = (unsigned char ***) calloc(Map->NY, sizeof(unsigned char **))))
    ReportError((char *) Routine, 1);
  for (i=0; i<Map->NY; i++) {
    if (!((SubDir)[i] = (unsigned char **) calloc(Map->NX, sizeof(unsigned char*))))
	  ReportError((char *) Routine, 1);
	for (j=0; j<Map->NX; j++) {
      if (!(SubDir[i][j] = (unsigned char *)calloc(NDIRS, sizeof(unsigned char ))))
		ReportError((char *) Routine, 1);
    }
  }

  if (!(SubTotalDir = (unsigned int **)calloc(Map->NY, sizeof(unsigned int *))))
    ReportError((char *) Routine, 1);
  for (i=0; i<Map->NY; i++) {
    if (!(SubTotalDir[i] = (unsigned int *)calloc(Map->NX, sizeof(unsigned int))))
      ReportError((char *) Routine, 1);
  }
  
  /* reset the saturated subsurface flow to zero */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
		SoilMap[y][x].SatFlow = 0;
		SoilMap[y][x].RoadInt = 0;
      }
    }
  }

  if (Options->FlowGradient == WATERTABLE)
    HeadSlopeAspect(Map, TopoMap, SoilMap, SubFlowGrad, SubDir, SubTotalDir);

  /* next sweep through all the grid cells, calculate the amount of
     flow in each direction, and divide the flow over the surrounding
     pixels */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
		if (Options->FlowGradient == TOPOGRAPHY){
		  SubTotalDir[y][x] = TopoMap[y][x].TotalDir;
	      SubFlowGrad[y][x] = TopoMap[y][x].FlowGrad;
		  for (k = 0; k < NDIRS; k++) 
	        SubDir[y][x][k] = TopoMap[y][x].Dir[k];
		}
		BankHeight = (Network[y][x].BankHeight > SoilMap[y][x].Depth) ?
	    SoilMap[y][x].Depth : Network[y][x].BankHeight;
	    Adjust = Network[y][x].Adjust;
	    fract_used = 0.0f;
		water_out_road = 0.0;
		
		if (!channel_grid_has_channel(ChannelData->stream_map, x, y)) {
	      for (k = 0; k < NDIRS; k++) {
			fract_used += (float) SubDir[y][x][k];
		  }
		  if (SubTotalDir[y][x] > 0)
	        fract_used /= (float) SubTotalDir[y][x];
		  else
	        fract_used = 0.;
		  
		  /* only bother calculating subsurface flow if water table is above bedrock */
		  if (SoilMap[y][x].TableDepth < SoilMap[y][x].Depth) {
	        depth = ((SoilMap[y][x].TableDepth > BankHeight) ?
				SoilMap[y][x].TableDepth : BankHeight);
			
			Transmissivity = CalcTransmissivity(SoilMap[y][x].Depth, depth,
				 SoilMap[y][x].KsLat,
				 SType[SoilMap[y][x].Soil - 1].KsLatExp,
                 SType[SoilMap[y][x].Soil - 1].DepthThresh);
			
			OutFlow = 
				(Transmissivity * fract_used * SubFlowGrad[y][x] * Dt) / (Map->DX * Map->DY);
			
			/* check whether enough water is available for redistribution */
			AvailableWater =
				CalcAvailableWater(VType[VegMap[y][x].Veg - 1].NSoilLayers,
				 SoilMap[y][x].Depth, VType[VegMap[y][x].Veg - 1].RootDepth,
				 SoilMap[y][x].Porosity, SoilMap[y][x].FCap,
				 SoilMap[y][x].TableDepth, Adjust);
			OutFlow = (OutFlow > AvailableWater) ? AvailableWater : OutFlow;
		  }
		  else {
	        depth = SoilMap[y][x].Depth;
	        OutFlow = 0.0f;
		  }
		  
		  /* compute road interception if water table is above road cut */
		  if (SoilMap[y][x].TableDepth < BankHeight &&
			  channel_grid_has_channel(ChannelData->road_map, x, y)) {
		    if (SubTotalDir[y][x] > 0)
	          fract_used = ((float) Network[y][x].fraction /
			    (float)SubTotalDir[y][x]);
			else
	          fract_used = 0.;
			Transmissivity =
				 CalcTransmissivity(BankHeight, SoilMap[y][x].TableDepth,
				 SoilMap[y][x].KsLat,
				 SType[SoilMap[y][x].Soil - 1].KsLatExp,
                 SType[SoilMap[y][x].Soil - 1].DepthThresh);
			
			water_out_road = (Transmissivity * fract_used *
			      SubFlowGrad[y][x] * Dt) / (Map->DX * Map->DY);
			
			AvailableWater =
				CalcAvailableWater(VType[VegMap[y][x].Veg - 1].NSoilLayers,
				 BankHeight, VType[VegMap[y][x].Veg - 1].RootDepth,
				 SoilMap[y][x].Porosity,
				 SoilMap[y][x].FCap,
				 SoilMap[y][x].TableDepth, Adjust);
			
			water_out_road = 
				(water_out_road > AvailableWater) ? AvailableWater : water_out_road;
			
			/* increase lateral inflow to road channel */
			SoilMap[y][x].RoadInt = water_out_road;
			channel_grid_inc_inflow(ChannelData->road_map, x, y,
				    water_out_road * Map->DX * Map->DY);
		  }
		  /* Subsurface Component - Decrease water change by outwater */
		  SoilMap[y][x].SatFlow -= OutFlow + water_out_road;
		  
		  /* Assign the water to appropriate surrounding pixels */
		  if (SubTotalDir[y][x] > 0)
	        OutFlow /= (float) SubTotalDir[y][x];
		  else
	        OutFlow = 0.;
		  
		  for (k = 0; k < NDIRS; k++) {
	        int nx = xdirection[k] + x;
	        int ny = ydirection[k] + y;
	        if (valid_cell(Map, nx, ny)) {
	          SoilMap[ny][nx].SatFlow += OutFlow * SubDir[y][x][k];
			}
		  }
		}
	    else {			/* cell has a stream channel */
	      if (SoilMap[y][x].TableDepth < BankHeight &&
	        channel_grid_has_channel(ChannelData->stream_map, x, y)) {
			float gradient = 4.0 * (BankHeight - SoilMap[y][x].TableDepth);
			if (gradient < 0.0)
	          gradient = 0.0;
			Transmissivity =
				CalcTransmissivity(BankHeight, SoilMap[y][x].TableDepth,
				 SoilMap[y][x].KsLat,
				 SType[SoilMap[y][x].Soil - 1].KsLatExp,
                 SType[SoilMap[y][x].Soil - 1].DepthThresh);

			OutFlow = (Transmissivity * gradient * Dt) / (Map->DX * Map->DY);
			
			/* check whether enough water is available for redistribution */
			AvailableWater = 
				 CalcAvailableWater(VType[VegMap[y][x].Veg - 1].NSoilLayers,
				 BankHeight, VType[VegMap[y][x].Veg - 1].RootDepth,
				 SoilMap[y][x].Porosity,
				 SoilMap[y][x].FCap,
				 SoilMap[y][x].TableDepth, Adjust);
			
			OutFlow = (OutFlow > AvailableWater) ? AvailableWater : OutFlow;
			
			/* remove water going to channel from the grid cell */
			SoilMap[y][x].SatFlow -= OutFlow;
			
			/* contribute to channel segment lateral inflow */
			channel_grid_inc_inflow(ChannelData->stream_map, x, y,
				    OutFlow * Map->DX * Map->DY);
			
			SoilMap[y][x].ChannelInt += OutFlow;
		  }
		}
      }
    }
  }

 for(i=0; i<Map->NY; i++) { 
    free(SubTotalDir[i]);
    free(SubFlowGrad[i]);
    for(j=0; j<Map->NX; j++){
      free(SubDir[i][j]);
    }
    free(SubDir[i]);
  }
  free(SubDir);
  free(SubTotalDir);
  free(SubFlowGrad);

  /**********************************************************************/
  /* Dump saturation extent file to screen.
     Saturation extent is based on the number of pixels with a water table 
     that is at least MTHRESH of soil depth. */ 
  
  count =0;
  totalcount = 0;
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	     mgrid = (SoilMap[y][x].Depth - SoilMap[y][x].TableDepth)/SoilMap[y][x].Depth;
	     if (mgrid > MTHRESH) 
		   count += 1;
		 totalcount += 1;
      }
    }
  }
 
  sat = 100.*((float)count/(float)totalcount);
  
  sprintf(satoutfile, "%ssaturation_extent.txt", DumpPath);
  
  if((fs = fopen(satoutfile,"a")) == NULL){
    printf("Cannot open saturation extent output file.\n");
    exit(0);
  }
  
  SPrintDate(&(Time->Current), buffer);
  fprintf(fs, "%-20s %.4f \n", buffer, sat); 
  fclose(fs);    
}

