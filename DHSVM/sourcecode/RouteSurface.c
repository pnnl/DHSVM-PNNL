/*
 * SUMMARY:      RouteSurface.c - Route surface flow
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Route surface flow
 * DESCRIP-END.
 * FUNCTIONS:    RouteSurface()
 * Modification: Changes are made to exclude the impervious channel cell (with 
                 a non-zero impervious fraction) from surface routing. In the original
				 code, some impervious channel cells are routed to themselves causing 
				 overestimated runoff in those cells (Ning, 2013).
              
 * $Id: RouteSurface.c, v3.1.2  2013/3/21   Ning Exp $  
 */
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "slopeaspect.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
/*****************************************************************************
  RouteSurface()

  If the watertable calculated in WaterTableDepth() was negative, then water is
  ponding on the surface.  At the moment no ponding of water is allowed in
  DHSVM, and the "excess" water is routed to the outlet one pixel per time step
  However, if the pixel contains an impervious fraction, then the surface water
  is immediately routed to the nearest downslope pixel that contains a channel.
  The net effect is that all pixels that have an impervious area are directly
  connected (over the coarse of a single time step) to the channel network, this
  assumption is likely to be true for small urban basins, and perhaps even for
  large rural basins with some urban development

  If Overland Routing = KINEMATIC, then "excess" water is routed to the outlet 
  using a infinite difference approximation to the kinematic wave solution of 
  the Saint-Venant equations.

*****************************************************************************/
void RouteSurface(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
		  SOILPIX ** SoilMap, OPTIONSTRUCT *Options,
		  UNITHYDR ** UnitHydrograph,
		  UNITHYDRINFO * HydrographInfo, float *Hydrograph,
		  DUMPSTRUCT *Dump, VEGPIX ** VegMap, VEGTABLE * VType, 
		  SOILTABLE *SType, CHANNEL *ChannelData, SEDPIX **SedMap,
		  PRECIPPIX **PrecipMap, SEDTABLE *SedType,
		  float Tair, float Rh, float *SedDiams)
{
  const char *Routine = "RouteSurface";
  int Lag;			/* Lag time for hydrograph */
  int Step;
  float StreamFlow;
  int TravelTime;
  int WaveLength;
  TIMESTRUCT NextTime;
  TIMESTRUCT VariableTime;
  int i, j, x, y, n, k;         /* Counters */
  float **Runon;                /* (m3/s) */

  /*************************** Kinematic wave routing**************************************** */
  float knviscosity;           /* kinematic viscosity JSL */  
  double slope, alpha, beta;   /* Slope is manning's slope; 
								alpha is channel parameter including wetted perimeter,
								manning's n, and manning's slope.  
								Beta is 3/5 */
  double outflow;              /* Outflow of water from a pixel during a sub-time step (m3/s)
							   outflow is not entirely true for channel cells*/
  double sedoutflow;           /* Outflow used for sediment routing purposes (m3/s) */
  float VariableDT;            /* Maximum stable time step (s) */  
  float **SedIn, SedOut;       /* (m3/m3) */  
  float DR;                    /* Potential erosion due to leaf drip */ 
  float DS;                    /* Median particle diameter (m) */
  float Cd;                    /* Drag coefficient */
  float vs, vs_last;           /* Settling velocity (m/s)*/
  float Rn;                    /* Particle Reynolds number */
  float h;                     /* Water depth (m) */
  float term1, term2, term3;
  float streampower;           /* Unit streampower from KINEROS (M/s) */
  float TC;                    /* Transport capacity (m3/m3) */
  float Fw;                    /* Water depth correction factor */
  float floweff;               /* Flow efficiency similar to Morgan */
  int sedbin;                  /* Particle bin that erosion is added to */
  /* Check to see if calculations for surface erosion should be done */
  if (Options->SurfaceErosion) {
     if ((SedIn = (float **) calloc(Map->NY, sizeof(float *))) == NULL) 
        ReportError((char *) Routine, 1);
     for (y = 0; y < Map->NY; y++) {
        if ((SedIn[y] = (float *) calloc(Map->NX, sizeof(float))) == NULL) {
	    ReportError((char *) Routine, 1);
        }
     }
     /* Initiazlied Variables */
     for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	   SedIn[y][x] = 0;
	}
     }
  }

  /* Allocate memory for Runon Matrix */
  if (Options->HasNetwork)  {
     if ((Runon = (float **) calloc(Map->NY, sizeof(float *))) == NULL) 
	 ReportError((char *) Routine, 1);
     for (y = 0; y < Map->NY; y++) {
	 if ((Runon[y] = (float *) calloc(Map->NX, sizeof(float))) == NULL) {
	    ReportError((char *) Routine, 1);
	 }
     }
     /* Initialize Runon variables */
     for (y = 0; y < Map->NY; y++) {
	for (x = 0; x < Map->NX; x++) {
	   Runon[y][x] = 0.;
	}
     }
     /* Option->Routing = false when routing = conventional */
     if(!Options->Routing) {
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	    if (INBASIN(TopoMap[y][x].Mask)) {
		SoilMap[y][x].Runoff = SoilMap[y][x].IExcess; 
		SoilMap[y][x].IExcess = 0;
		SoilMap[y][x].DetentionIn = 0;
	    }
         }
	} 
	for (y = 0; y < Map->NY; y++) {
	  for (x = 0; x < Map->NX; x++) {
	     if (INBASIN(TopoMap[y][x].Mask)) {
		 if (!channel_grid_has_channel(ChannelData->stream_map, x, y)) {
		    if (VType[VegMap[y][x].Veg - 1].ImpervFrac > 0.0) {
			/* Calculate the outflow from impervious portion of urban cell straight to nearest channel cell */		
			SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess += 
				(1 - VType[VegMap[y][x].Veg - 1].DetentionFrac) * 
				VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].Runoff;
			/* Retained water in detention storage */
			SoilMap[y][x].DetentionIn = VType[VegMap[y][x].Veg - 1].DetentionFrac * 
				VType[VegMap[y][x].Veg - 1].ImpervFrac * SoilMap[y][x].Runoff;		
			/* Retained water in Detention storage routed to channel */   
			SoilMap[y][x].DetentionStorage += SoilMap[y][x].DetentionIn;
			SoilMap[y][x].DetentionOut = SoilMap[y][x].DetentionStorage * VType[VegMap[y][x].Veg - 1].DetentionDecay;
			SoilMap[TopoMap[y][x].drains_y][TopoMap[y][x].drains_x].IExcess += SoilMap[y][x].DetentionOut;
			SoilMap[y][x].DetentionStorage -= SoilMap[y][x].DetentionOut;
			if (SoilMap[y][x].DetentionStorage < 0.0) 
			   SoilMap[y][x].DetentionStorage = 0.0;
			/* Route the runoff from pervious portion of urban cell to the neighboring cell */       
			for (n = 0; n < NDIRS; n++) {
			   int xn = x + xdirection[n];
			   int yn = y + ydirection[n];
			   if (valid_cell(Map, xn, yn)) {
				SoilMap[yn][xn].IExcess += (1 - VType[VegMap[y][x].Veg - 1].ImpervFrac) * SoilMap[y][x].Runoff 
					*((float) TopoMap[y][x].Dir[n] /(float) TopoMap[y][x].TotalDir);
			   }
			}
		    }
		    else {
			for (n = 0; n < NDIRS; n++) {
			   int xn = x + xdirection[n];
			   int yn = y + ydirection[n];
			   if (valid_cell(Map, xn, yn)) {
				SoilMap[yn][xn].IExcess += SoilMap[y][x].Runoff *((float)TopoMap[y][x].Dir[n]/(float)TopoMap[y][x].TotalDir);
			   }
			}
		    }
		 }
		 else if (channel_grid_has_channel(ChannelData->stream_map, x, y)){
			SoilMap[y][x].IExcess += SoilMap[y][x].Runoff;
               }
	     }
	  }
      }
}/* end if Options->routing = conventional */
/***********************************************************************************************************************/ 
else {/* Begin code for kinematic wave routing. */ 
	NextTime = *Time;
	VariableTime = *Time;
      
	/* Holds the value of the next DHSVM time step. */ 
	IncreaseTime(&NextTime);
      
    /* Use the Courant condition to find the maximum stable time step (in seconds). Must be an even increment of Dt. */
    VariableDT = FindDT(SoilMap, Map, Time, TopoMap, SType); 
      
	/*Numcells = cells in the basin*/
    for (k = 0; k < Map->NumCells; k++) {
		y = Map->OrderedCells[k].y;
		x = Map->OrderedCells[k].x;
		SoilMap[y][x].Runoff = 0.;		  
		if(Options->SurfaceErosion){
			SedMap[y][x].SedFluxOut = 0.;
			SedMap[y][x].Erosion = 0.;
		}
	}
	
	/* estimate kinematic viscosity through interpolation JSL */
    knviscosity=viscosity(Tair, Rh);
    /* converting units to m2/sec */
    knviscosity /= 1000. * 1000.;
	
	/* Must loop through surface routing multiple times within one DHSVM  model time step. */
	while (Before(&(VariableTime.Current), &(NextTime.Current))) {
		/* Loop thru all of the cells in descending order of elevation */
		for (k = (Map->NumCells)-1; k >-1;  k--) {
			y = Map->OrderedCells[k].y;
			x = Map->OrderedCells[k].x;
			outflow = SoilMap[y][x].startRunoff;   
			slope = TopoMap[y][x].Slope;
			if (slope == 0) slope=0.0001;
			else if (slope < 0) {
				printf("negative slope in RouteSurface.c\n");
				exit(0);
			}
			beta = 3./5.;
			alpha = pow(SType[SoilMap[y][x].Soil-1].Manning* pow((double)Map->DX,2./3.)/sqrt(slope),beta);
			
			/* Calculate discharge (m3/s) from the grid cell using an explicit 
			  finite difference solution of the linear kinematic wave. */
			if(Runon[y][x] > 0.0001 || outflow > 0.0001) {
				outflow = ((VariableDT/Map->DX)*Runon[y][x] + alpha*beta*outflow * pow((outflow+Runon[y][x])/2.0,beta-1.) +
					SoilMap[y][x].IExcess*Map->DX*VariableDT/Time->Dt)/ ((VariableDT/Map->DX) + alpha*beta*pow((outflow+
					  Runon[y][x])/2.0, beta-1.));
			}
			else if(SoilMap[y][x].IExcess > 0.0)
				outflow = SoilMap[y][x].IExcess*Map->DX*Map->DY/Time->Dt; 
			else
				outflow = 0.0; 
			if(outflow < 0.0) 
				outflow = 0.0; 
			
			/* Save flow depth and outflow for sediment routing */
			sedoutflow = outflow;
			h = SoilMap[y][x].IExcess;
			if (channel_grid_has_channel(ChannelData->stream_map, x, y)  
				  || (channel_grid_has_channel(ChannelData->road_map, x, y) 
				  && !channel_grid_has_sink(ChannelData->road_map, x, y))) {
			    /*  Recalculate for pixels with channels for sediment erosion  */
				if(Runon[y][x] > 0.0001 || outflow > 0.0001) {
					sedoutflow = ((VariableDT/Map->DX)*Runon[y][x] + alpha*beta*outflow * pow((outflow+Runon[y][x])/2.0,beta-1.) +
						(SoilMap[y][x].IExcessSed)*Map->DX*VariableDT/Time->Dt)/((VariableDT/Map->DX) + alpha*beta*pow((outflow+
						Runon[y][x])/2.0, beta-1.));
			     }
			     else if(SoilMap[y][x].IExcessSed > 0.0)
						  sedoutflow = SoilMap[y][x].IExcessSed*Map->DX*Map->DY/Time->Dt; 
				 else
					 sedoutflow = 0.0; 
				 if(sedoutflow < 0.0) 
					 sedoutflow = 0.0; 
				 outflow = 0.0;
				 h = SoilMap[y][x].IExcessSed;
				 if(sedoutflow > (SoilMap[y][x].IExcessSed*(Map->DX*Map->DY)/Time->Dt + Runon[y][x])) 
					 sedoutflow = SoilMap[y][x].IExcessSed*(Map->DX*Map->DY)/Time->Dt + (Runon[y][x]);
				 SoilMap[y][x].IExcessSed += (Runon[y][x] - sedoutflow)* VariableDT/(Map->DX*Map->DY);
			}	  
			  /*Make sure calculated outflow doesn't exceed available water, and update surface water storage */
			  if(outflow > (SoilMap[y][x].IExcess*(Map->DX*Map->DY)/Time->Dt + Runon[y][x])) 
				  outflow = SoilMap[y][x].IExcess*(Map->DX*Map->DY)/Time->Dt + (Runon[y][x]);
			  
			  SoilMap[y][x].IExcess += (Runon[y][x] - outflow)* VariableDT/(Map->DX*Map->DY);
			  
			  /*************************************************************/
			  /* PERFORM HILLSLOPE SEDIMENT ROUTING.                       */
			  /*************************************************************/
			  SedOut = 0.;   
			  if(Options->SurfaceErosion) {
				  /* Only perform sediment routing if there is outflow, and and for erodibility sil types */
				  if((sedoutflow > 0.) && (SedType[SoilMap[y][x].Soil-1].KIndex > -999.)){
					  DS = SedType[SoilMap[y][x].Soil-1].d50 * (float)MMTOM;

					  /* calculate unit streampower = u*S (m/s)  */
					  streampower = (sedoutflow/Map->DX/h)*slope;

					  /* avoid dividing by zero */ 
					  if (h <= 0.) streampower = 0.;
					  
					  /* Only continue sediment routing if there is depth greater 
					  than 0.001m  and streampower is greater than critical streampower */
					  
					  if((h > 0.001) && (streampower > SETTLECRIT)){
					  /* First find potential erosion due to rainfall Morgan et al. (1998). 
					  Momentum squared of the precip is determined in MassEnergyBalance.c 
					  Converting from (kg/m2*s) to (m3/s*m) */
						  if (h <= PrecipMap[y][x].Dm)
							  Fw = 1.;
						  else
							  Fw = exp(1 - (h/PrecipMap[y][x].Dm));
						  
						  /* If there is an understory, it is assumed to cover the entire grid cell. DR is in kg/m^2*s) */
                                            /* Modified by Ning (3/21/2013, 4:04 pm) */
						  if (VType[VegMap[y][x].Veg - 1].OverStory == TRUE) {
							  if (VType[VegMap[y][x].Veg - 1].UnderStory == FALSE)
								  DR = SedType[SoilMap[y][x].Soil-1].KIndex * Fw * PrecipMap[y][x].MomentSq; 
							  if (VType[VegMap[y][x].Veg - 1].UnderStory == TRUE)
								  DR = SedType[SoilMap[y][x].Soil-1].KIndex * Fw * 
                                                             (1 - VType[VegMap[y][x].Veg - 1].Fract[1])*PrecipMap[y][x].MomentSq;
						  }
						  else if (VType[VegMap[y][x].Veg - 1].UnderStory == TRUE)
							  /* There is no Overstory, then (1-VType->Fract[0]) is the fraction of understory */	
							  DR = SedType[SoilMap[y][x].Soil-1].KIndex * Fw * 
                                                         (1-VType[VegMap[y][x].Veg - 1].Fract[0])*PrecipMap[y][x].MomentSq; 
						  /* no vegetation */
						  else
							  DR = 0.;
						  
						  /* converting units to m3 m-1 s-1*/
						  DR = DR/PARTDENSITY * Map->DX;
						  
						  /* Calculate settling velocity iteratively (kineros) initial guess */
						  vs = sqrt((4./3.) * G * ((PARTDENSITY/WATER_DENSITY) - 1.)*DS);
						  vs_last = 999.;
						  
						  while (fabs(vs_last - vs) > (0.0001 * vs_last)) {
							  vs_last = vs;
							  Rn = (vs * DS) / knviscosity; 
							  Cd = (24./Rn) + (3./(pow((double)Rn, 0.5))) + 0.34;
							  vs = sqrt((4./3.) * G * ((PARTDENSITY/WATER_DENSITY) - 1.)*(DS/Cd));
						  }
						  floweff = 0.79*exp(-0.6*SedType[SoilMap[y][x].Soil-1].Cohesion.mean); 

						 /* calculate transport capacity (eq. 7 kineros) */
						 TC = 0.05/(DS*pow((PARTDENSITY/WATER_DENSITY-1.),2.))*sqrt(slope*h/G)*(streampower-SETTLECRIT);
						 
						 /* Calculate sediment mass balance. */
						 term1 = (TIMEWEIGHT/Map->DX);
						 term2 = alpha/(2.*VariableDT);
						 term3 = (1.-TIMEWEIGHT)/Map->DX;
						 
						 SedOut = (SedIn[y][x]*(term1*Runon[y][x]-term2*pow((double)Runon[y][x], beta)) +
							 SedMap[y][x].OldSedOut*(term2*pow((double)SoilMap[y][x].startRunoff, beta) -
						  term3*SoilMap[y][x].startRunoff) + SedMap[y][x].OldSedIn*(term2*pow((double)SoilMap[y][x].startRunon, beta) + 
						 term3*SoilMap[y][x].startRunon) + DR + floweff*Map->DY*vs*TC)/
						 (term2*pow(sedoutflow, beta) + term1*sedoutflow + floweff*Map->DY*vs);
						 
						 if(SedOut >= TC) 
							 SedOut = TC;
						 SedMap[y][x].OldSedOut = SedOut;
					     SedMap[y][x].OldSedIn = SedIn[y][x];
						 SedMap[y][x].SedFluxOut += (SedOut*sedoutflow*VariableDT);  /* total sediment (m3) */
						 SedMap[y][x].Erosion += (SedIn[y][x]*Runon[y][x] - SedOut*sedoutflow)*VariableDT/(Map->DX*Map->DY)*1000.;  /* total depth of erosion (mm) */
					  } /* end if((h > DS) && (streampower > SETTLECRIT){ */
					  else {
						  SedMap[y][x].OldSedOut = 0.;
						  SedMap[y][x].OldSedIn = 0.;
						  SedOut = 0.;
					  }	    
				  } /*   end if((sedoutflow > 0.) */
				  else {
					  SedMap[y][x].OldSedOut = 0.;
					  SedMap[y][x].OldSedIn = 0.;
					  SedOut = 0.;
				  }	     
			  } /* end of if Options->SurfaceErosion */
			  
			  /* Save sub-timestep runoff for q(i)(t-1) and q(i-1)(t-1) of next time step. */
			  SoilMap[y][x].startRunoff = sedoutflow;
	          SoilMap[y][x].startRunon = Runon[y][x];
			  
			  /******************************* Calculate total runoff in m/dhsvm timestep*************************************************/
	          /* This is serving the purpose of holding this value for Cournat condition calculation. It does not trully repreent
	          the Runoff, because if there is a channel there is no outflow. Instead, IExcess is updated based on Runon in the same manner
	          of the original DHSVM */
			  SoilMap[y][x].Runoff += sedoutflow*VariableDT/(Map->DX*Map->DY); 
			  
			  /* Sediment from pixels with channels goes into the channel.  This assumes that all surface erosion is of the smallest particle sizes (first size, index 0) */
	          /* Note that stream_map and road_map are indexed by [x][y], unlike the other "map"-type variables. */
			  
			  if((Options->SurfaceErosion)&&(SedOut > 0.)){
				  
				  /* Determine which particle bin sediment gets added to SedDiams */
				  sedbin = 0;
				  if (SedType[SoilMap[y][x].Soil-1].d50 > SedDiams[NSEDSIZES-1])
					  sedbin = NSEDSIZES-1;
				  else {
					  for (j=0; j < NSEDSIZES; j++){
						  if (SedType[SoilMap[y][x].Soil-1].d50 <= SedDiams[j]){
							  sedbin = j - 1;
							  break;
						  }
					  }
					  if (sedbin < 0) sedbin = 0;
				  }
				  
				  if (channel_grid_has_channel(ChannelData->stream_map, x, y)) {
				      /* Converting SedOut from m3/m3 to kg for channel routing */
					  ChannelData->stream_map[x][y]->channel->sediment.overlandinflow[sedbin] += 
						  (SedOut*sedoutflow*VariableDT*PARTDENSITY)/(Map->DX*Map->DY);
					  SedOut = 0.;
				  }
				  
				  if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
					  /* Converting SedOut from m3/m3 to kg for channel routing */
 					  ChannelData->road_map[x][y]->channel->sediment.overlandinflow[sedbin] += 
 					  (SedOut*sedoutflow*VariableDT*PARTDENSITY)/(Map->DX*Map->DY); 
 	                  SedOut = 0.; 
				  } 
			  }	  
			  
			  /********************************************************************************************************/
			  /* Redistribute surface water to downslope pixels. */
			  if(outflow > 0.) {  
				  for (n = 0; n < NDIRS; n++) {
					  int xn = x + xdirection[n];
					  int yn = y + ydirection[n];
					  
					  /* If a channel cell runoff does not go to downslope pixels. */
					  if (valid_cell(Map, xn, yn)) {
						  if (INBASIN(TopoMap[yn][xn].Mask)) {
							  Runon[yn][xn] += outflow * ((float) TopoMap[y][x].Dir[n]/(float) TopoMap[y][x].TotalDir);
							  
							  /* No need to distribute sediment if there isn't any*/
							  if((Options->SurfaceErosion)&&(SedOut > 0.)){ 	
								  SedIn[yn][xn] += SedOut * ((float) TopoMap[y][x].Dir[n] /(float) TopoMap[y][x].TotalDir);
							  }
						  }
					  }

				  } /* end loop thru possible flow directions */
			  }
			  
			  /* Initialize runon for next timestep. */
			  Runon[y][x] = 0.0;
			  /* Initialize SedIn for next timestep. */
			  if(Options->SurfaceErosion) 
				  SedIn[y][x] = 0.0;
} /* end loop thru ordered basin cells */

/* Increases time by VariableDT. */


IncreaseVariableTime(&VariableTime, VariableDT, &NextTime);

/*************************************************************/
} /* End of internal time step loop. */
}/* End of code added for kinematic wave routing. */
    
if(Options->SurfaceErosion) {
	for (y = 0; y < Map->NY; y++) {
		free(SedIn[y]);
	}
	free(SedIn);
}
    
for (y = 0; y < Map->NY; y++) {
	free(Runon[y]);
}
free(Runon);
}
  
/* MAKE SURE THIS WORKS WITH A TIMESTEP IN SECONDS */
else {			/* No network, so use unit hydrograph method */
	for (y = 0; y < Map->NY; y++) {
		for (x = 0; x < Map->NX; x++) {
			if (INBASIN(TopoMap[y][x].Mask)) {
				TravelTime = (int) TopoMap[y][x].Travel;
				if (TravelTime != 0) {
					WaveLength = HydrographInfo->WaveLength[TravelTime - 1];
					for (Step = 0; Step < WaveLength; Step++) {
						Lag = UnitHydrograph[TravelTime - 1][Step].TimeStep;
						Hydrograph[Lag] += SoilMap[y][x].Runoff * UnitHydrograph[TravelTime - 1][Step].Fraction;

					}
					SoilMap[y][x].Runoff = 0.0;
				}
			}
          }
    }
    
    StreamFlow = 0.0;
    for (i = 0; i < Time->Dt; i++)
      StreamFlow += (Hydrograph[i] * Map->DX * Map->DY) / Time->Dt;
    
    /* Advance Hydrograph */
    for (i = 0; i < Time->Dt; i++) {
      for (j = 0; j < HydrographInfo->TotalWaveLength - 1; j++) {
	Hydrograph[j] = Hydrograph[j + 1];
      }

    }
    
    /* Set the last elements of the hydrograph to zero */
    for (i = 0; i < Time->Dt; i++)
      Hydrograph[HydrographInfo->TotalWaveLength - (i + 1)] = 0.0;
    
    PrintDate(&(Time->Current), Dump->Stream.FilePtr);
    fprintf(Dump->Stream.FilePtr, " %g\n", StreamFlow);
  }
}

/*****************************************************************************
  FindDT()
  Find the variable time step that will satisfy the courant condition for stability 
  in overland flow routing.
*****************************************************************************/
float FindDT(SOILPIX **SoilMap, MAPSIZE *Map, TIMESTRUCT *Time, 
	     TOPOPIX **TopoMap, SOILTABLE *SType)
{
  int x, y;
  /* JSL: slope is manning's slope; alpha is channel parameter including wetted perimeter, 
     manning's n, and manning's slope.  Beta is 3/5 */
  float slope;
  double beta, alpha;
  double Ck;
  float DT, minDT;
  float numinc;
  float maxRunoff;
  
  maxRunoff = -99.;
  minDT = 36000.;
  
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
		  if (SoilMap[y][x].Runoff >0.0){
			  slope = TopoMap[y][x].Slope;
			  if (slope <= 0) slope = 0.0001;
			  beta = 3./5.;
			  alpha = pow((double)SType[SoilMap[y][x].Soil-1].Manning *pow((double)Map->DX,(double)(2./3.))/sqrt(slope), (double)beta);
	  
	  /* Calculate flow velocity from discharge  using manning's equation. */
	  Ck = 1./(alpha*beta*pow((double)SoilMap[y][x].Runoff, beta -1.));
	  
	  if(SoilMap[y][x].Runoff > maxRunoff) {
	    maxRunoff = SoilMap[y][x].Runoff;
	  }

	  
	  if(Map->DX/Ck < minDT)
	    minDT = Map->DX/Ck;
	  }
      }
    }
  }
  /* Find the time step that divides evenly into Time->DT */
  
  numinc = (float) ceil((double)Time->Dt/minDT);
  DT = Time->Dt/numinc;
  
  if(DT > Time->Dt)
    DT = (float) Time->Dt;
  
  return DT;
}
/*****************************************************************************
  SedimentFlag()
  To determine when the surface erosion and kinematic routing will be 
  performed
*****************************************************************************/
void SedimentFlag(OPTIONSTRUCT *Options,  TIMESTRUCT * Time)
{
  int i;             /* Counter */
  int flag = 0;
  int oldrouting;
  double StartSed;      /* Start of the current or next surface erosion 
			calculation period (Julian day) */
  double EndSed;        /* End of the current or next surface erosion 
			calculation period (Julian day) */
  
  StartSed = 0.; 
  EndSed = 0.;      
			

  if ((Options->ErosionPeriod) && (Time->Current.Julian==Time->Start.Julian)){
    Options->OldSedFlag = 1;
  }
  
  if ((Options->Routing) && (Time->Current.Julian==Time->Start.Julian)){
    Options->OldRouteFlag = 1;
  }
  
  oldrouting = Options->Routing;

   /* Determine next period for calculating surface erosion */
  if (Time->NSETotalSteps == 0){
    StartSed = Time->Start.Julian; 
    EndSed = Time->End.Julian; 
  }
  else{
    for (i = 0; i < Time->NSETotalSteps; i++){
      if((Time->Current.Julian <= Time->EndSed[i].Julian) && flag < 1){
	StartSed = Time->StartSed[i].Julian; 
	EndSed = Time->EndSed[i].Julian; 
	flag = 1;
      }
    }
  }
  
  if((Options->Sediment) && (Options->ErosionPeriod)) {
     if ((Time->Current.Julian >= StartSed) && (Time->Current.Julian <= EndSed)){
      Options->SurfaceErosion=TRUE;
    }
    else Options->SurfaceErosion=FALSE;
  }

  if ((Options->OldSedFlag != Options->SurfaceErosion) && 
      (Time->Current.Julian != Time->Start.Julian)) {
    if (Options->SurfaceErosion)
      printf("Beginning surface erosion model calculations.\n");
    else
      printf("Ending surface erosion model calculations.\n");
  }
  Options->OldSedFlag=Options->SurfaceErosion;
  
  if(Options->SurfaceErosion){
    Options->Routing=TRUE;
    if (oldrouting!=Options->Routing)
      printf("Turning on kinematic routing calculations.\n");
  }
  
  else if (!Options->OldRouteFlag){
    Options->Routing=FALSE;
    if (oldrouting!=Options->Routing)
      printf("Turning off kinematic routing calculations.\n");
  }
}



































































