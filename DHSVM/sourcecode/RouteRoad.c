/*
 * SUMMARY:      RouteRoad.c - Route surface flow and erosion for Forest Roads
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Colleen O. Doten
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       dhsvm@hydro.washington.edu
 * ORIG-DATE:    Nov-03
 * DESCRIPTION:  Route surface flow and erosion for Forest Roads
 * DESCRIP-END.
 * FUNCTIONS:    RouteRoad()
 * COMMENTS:
 * $Id: RouteRoad.c,v 1.9 2004/08/18 01:01:32 colleen Exp $     
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "channel_grid.h"
#include "channel.h"

/*****************************************************************************
  Function name: RouteRoad()

  Purpose      : Calculate routing of water and sediment across the road 
                 surface using a four-point finite difference solution of the 
                 kinematic wave approximation of the Saint-Venant equations

  Required     :
  
  Returns      : void

  Modifies     :
 
  Comments     : This can only be run if the sediment model is run. Therefore
                 there are checks for the sediment model being run.
 
  Sources: 
   Smith, R.E., D.C. Goodrich, and C.L. Unkrich, Simulation of selected events on 
   the Catsop catchment by KINEROS, 1999, A report for the GCTE conference on 
   catchment scale erosion models, Catena, 36, pp. 457-475.

   Smith, R.E., D.C. Goodrich, and C.L. Unkrich, 1995, KINERSO - A Kinematic 
   Runoff and Erosion Model, In: Singh, V.J. (Ed.), Computer Models of Watershed
   Hydrology, Chapter 20, Water Resources Publications, 697-732.
   
  These two papers and the code that can be downloaded from www.tucson.ars.ag.gov/kineros
  have a different formulation of the equation for hydraulic erosion. The equation used 
  here, is eh = cg(Cmx-Cs)A. cg = vs/h*CH (variables defined below).

*****************************************************************************/
void RouteRoad(MAPSIZE * Map, TIMESTRUCT * Time, TOPOPIX ** TopoMap,
	       SOILPIX ** SoilMap, ROADSTRUCT ** Network, SOILTABLE * SType,
	       CHANNEL * ChannelData, PRECIPPIX ** PrecipMap, SEDPIX **SedMap,
	       float Tair, float Rh, float *SedDiams) 
{
  const char *Routine = "RouteRoad";
  int i,j,x,y;                   /* Counters */
  float dx, dy;                  /* Road grid cell dimensions (m)*/
  float cells;                   /* Number of road grid cells in a basin grid 
				    cell */
  float roadwater;               /* Depth (m) of water on the road surface */
  float knviscosity;             /* kinematic viscosity (mm2/s) */ 
  double slope;                   /* Slope of road surface the water travels 
				    over (m/m)*/    
  float *Runon;                  /* Water on moving downslope along the road (m) */
  double beta, alpha;             /* For Mannings. alpha is channel parameter 
				    including wetted perimeter, n, and slope. */
  double outflow;                 /* Flow out of a road grid cell in the 
				    subtimestep (m3/s) */
  double check; 

  TIMESTRUCT NextTime;
  TIMESTRUCT VariableTime;
  float VariableDT;              /* Maximum stable time step (s) */  

  float ES;                      /* Rain splash erosion (m2/s)*/
  float ch;                      	/* damping effectiveness of surface water */
  float h;                       /* Water depth (m) */
  float k;                       	/* Reduction factor due to surface water depth */
  float DS;                      /* Median particle size (m) */
  float Cd;                      /* Drag coefficient */
  float vs, vs_last;             /* Settling velocity (m/s)*/
  float cg;                      /* Transfer rate coefficient (1/s)*/
  float CH;                      /* inversely related to soil cohesion or other restrictions on
				    soil entrainment by flowing water; = 1 during deposition */
  float Rn;                      /* Particle Reynolds number */
  float streampower;             /* Streampower (m/s) */
  float Cmx;                     /* Transport capacity (m3/m3) */
  float SedOut;                  /* Current local sediment concentration (m3/m3) */
  float *SedIn;                   /* Current inflow sediment concentration (m3/m3) */
  float term1, term2, term3;
  int sedbin;                   /* Particle bin that road erosion is added to */

  if ((Runon = (float *) calloc(CELLFACTOR, sizeof(float))) == NULL)
    ReportError((char *) Routine, 1);
  
  if ((SedIn = (float *) calloc(CELLFACTOR, sizeof(float))) == NULL)
    ReportError((char *) Routine, 1);

  NextTime = *Time;
  
  /* Holds the value of the next DHSVM time step. */ 
  IncreaseTime(&NextTime);

  knviscosity = viscosity(Tair, Rh);
 
  /* Since Network.IExcess stays in the cell it is generated in,
     route each basin grid cell, with a road, separately */
  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      if (INBASIN(TopoMap[y][x].Mask)) {
	if (channel_grid_has_channel(ChannelData->road_map, x, y)) {
	  
	   VariableTime = *Time;

	  /* Discretizing road into a grid for finite difference solution.
	     This assume the cells are oriented with the direction of flow. */
	  
	  dx = Network[y][x].FlowLength/(float) CELLFACTOR;
	  dy = dx; /* road grid cells are square */
	  cells = Network[y][x].RoadArea/(dx*dy);

	  slope = Network[y][x].FlowSlope;
	  if (slope == 0) slope=0.0001;
	  else if (slope < 0) {
	    printf("RouteRoad.c: negative slope\n");
	    exit(0);
	  }
	  
	  beta = 3./5.;
	  alpha = pow(Network[y][x].RoadClass->friction_road*pow((double)dx,2./3.)/sqrt(slope),beta);

	  /* Evenly distribute water over road surface. */
	  roadwater = (Network[y][x].IExcess * Map->DX * Map->DY)/
	    (Network[y][x].RoadArea);
	   
	  for (i = 0; i < CELLFACTOR; i++){
	    if(Network[y][x].h[i] < 0){
	      printf ("RouteRoad: Negative Network.IExcess(%e)\n",Network[y][x].h[i]);
	      exit(0);
	    }
	    Network[y][x].h[i] += roadwater;
	  }
	

	  /* Perform sediment calculations that only need to be
	     performed once for the coarse grid cell */
	  Network[y][x].Erosion = 0.;
	  SedMap[y][x].RoadSed = 0.;
	  DS = Network[y][x].RoadClass->d50_road * MMTOM;
	  
	  /* Calculate settling velocity iteratively 
	     initial guess */
	  vs = sqrt((4./3.) * G * ((PARTDENSITY/WATER_DENSITY) - 1.)*DS);
	  vs_last = 999.;
	  
	  while (fabs(vs_last - vs) > 0.0001 * vs_last) {
	    vs_last = vs;
	    Rn = (vs * DS * 1000. * 1000.) / knviscosity; 
	    Cd = (24./Rn) + (3./(pow((double)Rn, 0.5))) + 0.34;
	    vs = sqrt((4./3.) * G * ((PARTDENSITY/WATER_DENSITY) - 1.)*(DS/Cd));
	  }
		  
	  /* Use the Courant condition to find the maximum stable time step 
	     (in seconds). Must be an even increment of Dt. */
	  VariableDT = FindDTRoad(Network, Time, y, x, dx, beta, alpha);  
	  
	  /* Must loop through road segment routing multiple times within 
	     one DHSVM model time step. */
	  while (Before(&(VariableTime.Current), &(NextTime.Current))) {
	    
	    /* Loop through road grid cells starting at crown or road edge*/
	    for (i = 0; i < CELLFACTOR; i++){
	      
	      outflow = Network[y][x].startRunoff[i]; 
	      
	      /* Calculate discharge from the road segment using an explicit
		 finite difference solution of the linear kinematic wave. */
	      if(Runon[i] > 0.0001 || outflow > 0.0001) {
		outflow = ((VariableDT/dx)*Runon[i] + alpha*beta*outflow * 
			   pow((outflow+Runon[i])/2.,beta-1.) +
			   Network[y][x].h[i]*dx*VariableDT/Time->Dt)/
		  ((VariableDT/dx) + alpha*beta*pow((outflow+Runon[i])/2., beta-1.));
		
	      }
	      else if(Network[y][x].h[i] > 0.0)
		outflow = Network[y][x].h[i]*dy*dx/(float)Time->Dt;
	      else
		outflow = 0.0;
	      
	      if(outflow < 0.0) outflow = 0.0;
	      
	      h = Network[y][x].h[i];
	      
	      /*Update surface water storage.  Make sure calculated outflow 
		doesn't exceed available water.  Otherwise, update surface 
		water storage. */
	      if(outflow > (Network[y][x].h[i]*dy*dx)/(float)Time->Dt + Runon[i]){ 
		outflow = (Network[y][x].h[i]*dy*dx)/(float)Time->Dt + Runon[i];
	      }
	      	      
	      /* Accounting for rounding errors */
	      check = Network[y][x].h[i] + ((Runon[i]-outflow)*VariableDT/(dy*dx));

	      if ((check < .0000001) && (check > -.0000001))
		Network[y][x].h[i] = 0.;
	      else
		Network[y][x].h[i] += (Runon[i]-outflow)*VariableDT/(dy*dx);
	 
	      /*************************************************************/
	      /* PERFORM ROAD SEDIMENT ROUTING.                            */
	      /*************************************************************/
	      /* Don't perform sediment calculations if the road is paved, there is no outflow
		 or if the flow depth is less than the median particle size */
	      SedOut = 0.;

	      if((Network[y][x].RoadClass->erodibility_coeff < 999999 ) && 
		 (outflow > 0.) && (h > DS)){

		/* Calculating rainsplash erosion */
		ch = 656.; /* based on KINSED.for from KINEROS2 code */
		k = exp(-ch*h); 

		/* PrecipMap[y][x].RainFall/Time->Dt is the rainfall intensity (m/s)*/
		ES = Network[y][x].RoadClass->erodibility_coeff * k * 
		  pow((PrecipMap[y][x].RainFall/(float)Time->Dt), 2); 

		if (ES < 0.)
		  ES = 0.; 
		
		/* Calculating hydraulic erosion */
		if (Network[y][x].OldSedOut[i] <  Network[y][x].OldSedIn[i])
		  CH = 1.; /* upper limit, indicates deposition */
		else
		  CH = Network[y][x].RoadClass->erodibility_coeff_overland;
		
		cg = CH * vs/h; 
			
		/* Calculate unit streampower = u*S (m/s) */
		streampower = (outflow /(h * dx)) * slope;

		if (streampower > 0.0004){
		  Cmx = 0.05/(DS*pow((PARTDENSITY/WATER_DENSITY-1.),2.))*pow((slope*h/G),0.5) * 
		    (streampower - 0.0004);

		  /* Calculate sediment mass balance. */
		  term1 = (TIMEWEIGHT/dx);
		  term2 = alpha/(2.*VariableDT);
		  term3 = (1.-TIMEWEIGHT)/dx;
				  
		  SedOut = (SedIn[i]*(term1*Runon[i]-term2*pow((double)Runon[i], beta)) +
			    Network[y][x].OldSedOut[i]*(term2*pow((double)Network[y][x].startRunoff[i], beta) -
							term3*Network[y][x].startRunoff[i]) +
			    Network[y][x].OldSedIn[i]*(term2*pow((double)Network[y][x].startRunon[i], beta) + 
						       term3*Network[y][x].startRunon[i]) + ES + 
			    (cg*Cmx*alpha*pow(outflow, beta)))/
		    (term2*pow(outflow, beta) + term1*outflow + cg*alpha*pow(outflow,beta));
		  
		  if(SedOut >= Cmx) /* then deposition */
		    SedOut = Cmx;
		  
		  if ((Cmx > 1) || (SedIn[i] > 1) || (SedOut > 1)){
		    printf("RouteRoad: Invalid results Cmx(%e) SedIn(%e) SedOut(%e)\n",
			   Cmx, SedIn[i], SedOut);
		    printf("DS %f slope %f h %e outflow %e dx %f streampow %e\n", 
			   DS, slope, h, outflow, dx, streampower);
		  }

		}
	      } /* end  if((outflow > 0.) && (h > DS)){ */
	      else 
		SedOut = 0.0;

	      Network[y][x].OldSedOut[i] = SedOut;
	      Network[y][x].OldSedIn[i] = SedIn[i];
	      /* total depth of erosion (m) over entire grid cell */
	      Network[y][x].Erosion += (((SedIn[i]*Runon[i] - SedOut*outflow)*VariableDT)/
					(Map->DX * Map->DY))*(cells/(float)CELLFACTOR); 	      
	      	      
	      /* Save sub-timestep runoff for q(i)(t-1) and q(i-1)(t-1) of next time step. */
	      Network[y][x].startRunoff[i] = outflow;
	      Network[y][x].startRunon[i] = Runon[i];
	      
	      /* Redistribute surface water and sediment to downslope pixel. */
	      if(outflow > 0.){
		if(i < (CELLFACTOR-1)){
		  Runon[i+1] += outflow;
		  SedIn[i+1] += SedOut;
		}
		/* If last pixel send outflow off road */
		else {

		  /* Determine which particle bin sediment gets added to */
		  sedbin = 0;
		  if (Network[y][x].RoadClass->d50_road > SedDiams[NSEDSIZES-1])
		    sedbin = NSEDSIZES-1;
		  else {
		    for (j=0; j < NSEDSIZES; j++){
		      if (Network[y][x].RoadClass->d50_road <= SedDiams[j]){
			sedbin = j - 1;
			break;
		      }
		    }
		    if (sedbin < 0) sedbin = 0;
		  }

		  /* Multiple results by the number of cells in a row */
		  if(Network[y][x].RoadClass->crown == CHAN_OUTSLOPED){
		    SoilMap[y][x].IExcess += ((outflow*VariableDT)/(Map->DX * Map->DY))*
		      (cells/(float)CELLFACTOR);
		    
		    SedMap[y][x].RoadSed += ((SedOut*outflow*VariableDT)/
					     (Map->DX * Map->DY))*(cells/(float)CELLFACTOR);
		  }
		  /* If the road is crowned, then the same amount of outflow goes to 
		     the ditch and off the road edge into the same pixel. This is similar to 
		     culvert flow. 0.5 accounts for 1/2 the cells on are one side of the
		     crown. */
		  else if(Network[y][x].RoadClass->crown == CHAN_CROWNED){
		    channel_grid_inc_inflow(ChannelData->road_map, x, y, 
					    ((outflow*VariableDT)/
					     (Map->DX * Map->DY))*0.5*(cells/(float)CELLFACTOR));
		    
		    SoilMap[y][x].RoadInt += ((outflow*VariableDT)/(Map->DX * Map->DY))*0.5*
		      (cells/(float)CELLFACTOR);
		
		    SoilMap[y][x].IExcess += ((outflow*VariableDT)/(Map->DX * Map->DY))*0.5*
		      (cells/(float)CELLFACTOR);

		    /* Converting SedOut from m3/m3 to kg for channel routing */
 		    ChannelData->road_map[x][y]->channel->sediment.overroadinflow[sedbin] += 
 		      SedOut*outflow*VariableDT*PARTDENSITY*0.5*(cells/(float)CELLFACTOR); 

		  		    
 		    SedMap[y][x].RoadSed += ((SedOut*outflow*VariableDT)/ 
 					     (Map->DX * Map->DY))*0.5*(cells/(float)CELLFACTOR); 

		  }
		  else { /* INSLOPED and all goes to ditch */
		    channel_grid_inc_inflow(ChannelData->road_map, x, y, 
					    ((outflow*VariableDT)/
					     (Map->DX * Map->DY))*(cells/(float)CELLFACTOR));
		    
		    SoilMap[y][x].RoadInt += ((outflow*VariableDT)/(Map->DX * Map->DY))*
		      (cells/(float)CELLFACTOR);
		    
		    /* Converting SedOut from m3/m3 to kg for channel routing */
		    ChannelData->road_map[x][y]->channel->sediment.overroadinflow[sedbin] +=
		      SedOut*outflow*VariableDT*PARTDENSITY*(cells/(float)CELLFACTOR);

		  }
		}
	      }
	      /* Initialize for next timestep. */
	      Runon[i] = 0.0;
	      SedIn[i] = 0.0;
	    }
	    /* Increases time by VariableDT. */
	    IncreaseVariableTime(&VariableTime, VariableDT, &NextTime);
	    
	  } /* End of internal time step loop. */
	  /* Initialize for next DHSVM time step */
	  Network[y][x].IExcess = 0.0;
	}
      }
    }
  }/* End loop through basin grid cells */
  free(Runon);
  free(SedIn);
}

/*****************************************************************************
  FindDTRoad()
  Find the variable time step that will satisfy the courant condition for stability 
  in overland flow routing.
*****************************************************************************/

float FindDTRoad(ROADSTRUCT **Network, TIMESTRUCT *Time, int y, int x, 
	     float dx, float beta, float alpha)
{
  int i;           /* counter */
  float Ck;        /* Flow velocity based on discharge using Manning's equation. */
  float DT, minDT; /* seconds */
  float numinc;
  float runoff, maxRunoff;
  
  maxRunoff = -99.;
  minDT = Time->Dt; 
  
  for (i = 0; i < CELLFACTOR; i++){
    
    runoff = Network[y][x].startRunoff[i];
    if (runoff == 0) runoff = 0.000000001;
    
    Ck = 1./(alpha*beta*pow((double)runoff, beta -1.)); 
    
    if(runoff > maxRunoff)
      maxRunoff = runoff;
    
    if(dx/Ck < minDT)
      minDT = dx/Ck;
    
  }
  
  /* Find the time step that divides evenly into Time->DT */
  numinc = (float) ceil((double)Time->Dt/minDT);
  DT = Time->Dt/numinc;
  
  if(DT > Time->Dt)
    DT = (float) Time->Dt;
  
  return DT;
}

