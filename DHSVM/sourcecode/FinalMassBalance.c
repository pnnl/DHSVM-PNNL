/*
 * SUMMARY:      MassBalance.c - calculate basin-wide mass balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Battelle - Pacific Northwest National Laboratory
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Oct-96
 * DESCRIPTION:  Calculate water and sediment mass balance errors
 *               
 * DESCRIP-END.
 * FUNCTIONS:    FinalMassBalance()
 * COMMENTS:
 * $Id: FinalMassBalance.c,v 1.18 2004/08/18 01:01:28 colleen Exp $
 * Modification made on 2013/1/11
 * $Id: FinalMassBalance.c, v3.1.2 Ning Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"

/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  Only the runoff is calculated as a total volume instead
  of an average.  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().

  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void FinalMassBalance(FILES * Out, AGGREGATED * Total, WATERBALANCE * Mass,
		      OPTIONSTRUCT * Options, float roadarearatio)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float MassError;		/* mass balance error m  */
  float Input;
  float MWMMassError;
  float SedInput, SedOutput, SedMassError;

  NewWaterStorage = Total->Soil.IExcess + Total->Road.IExcess + 
    Total->CanopyWater + Total->SoilWater +
    Total->Snow.Swq + Total->Soil.SatFlow + Total->Soil.DetentionStorage;

  Output = Mass->CumChannelInt + ( Mass->CumRoadInt  -
    Mass->CumCulvertReturnFlow ) + Mass->CumET;

  Input = Mass->CumPrecipIn + Mass->CumSnowVaporFlux -
    Mass->CumCulvertReturnFlow;

  MassError = (NewWaterStorage - Mass->StartWaterStorage) +
   Output - Input;

  /* Print the runoff final balance results to the screen */
  fprintf(stderr, "\n  ********************************               Depth");
  fprintf(stderr, "\n  Runoff Final Mass Balance                        mm"); 
  fprintf(stderr, "\n  ********************************        ------------"); 
  fprintf(stderr, "\n  Total Inflow ...................        %12.3f", Input*1000);
  fprintf(stderr, "\n      Precip/Inflow ..............        %12.3f", Mass->CumPrecipIn*1000); 
  fprintf(stderr, "\n      SnowVaporFlux ..............        %12.3f", Mass->CumSnowVaporFlux*1000);
  fprintf(stderr, "\n  Total Outflow ..................        %12.3f", Output*1000);
  fprintf(stderr, "\n      ET .........................        %12.3f", Mass->CumET*1000);   
  fprintf(stderr, "\n      ChannelInt .................        %12.3f", Mass->CumChannelInt*1000); 
  fprintf(stderr, "\n      RoadInt ....................        %12.3f", (Mass->CumRoadInt  -
                                                                        Mass->CumCulvertReturnFlow) * 1000);
  fprintf(stderr, "\n  Storage Change .................        %12.3f", (NewWaterStorage - Mass->StartWaterStorage)*1000);
  fprintf(stderr, "\n      Initial Storage ............        %12.3f", Mass->StartWaterStorage*1000);
  fprintf(stderr, "\n      Final Storage ..............        %12.3f", NewWaterStorage*1000);
  fprintf(stderr, "\n          Final SWQ ..............        %12.3f", Total->Snow.Swq*1000);
  fprintf(stderr, "\n          Final Soil Moisture ....        %12.3f", (Total->SoilWater + Total->Soil.SatFlow)*1000);
  fprintf(stderr, "\n          Final Surface ..........        %12.3f", (Total->Soil.IExcess  + 
						                               Total->CanopyWater + Total->Soil.DetentionStorage)*1000);
  fprintf(stderr, "\n          Final Road Surface .....        %12.3f\n", Total->Road.IExcess*1000);
  fprintf(stderr, "\n  Mass added to glacier ..........        %12.3f\n", Total->Snow.Glacier*1000);
  fprintf(stderr, "  ******************************************************");
  fprintf(stderr, "\n  Mass Error (mm).................        %12.3f", MassError*1000);
  
    /* Print the runoff final balance results to the output file named final.mass.balance */
  fprintf(Out->FilePtr, "\n  ********************************               Depth");
  fprintf(Out->FilePtr, "\n  Runoff Final Mass Balance                        mm"); 
  fprintf(Out->FilePtr, "\n  ********************************        ------------"); 
  fprintf(Out->FilePtr, "\n  Total Inflow ...................        %12.3f", Input*1000);
  fprintf(Out->FilePtr, "\n      Precip/Inflow ..............        %12.3f", Mass->CumPrecipIn*1000); 
  fprintf(Out->FilePtr, "\n      SnowVaporFlux ..............        %12.3f", Mass->CumSnowVaporFlux*1000);
  fprintf(Out->FilePtr, "\n  Total Outflow ..................        %12.3f", Output*1000);
  fprintf(Out->FilePtr, "\n      ET .........................        %12.3f", Mass->CumET*1000);   
  fprintf(Out->FilePtr, "\n      ChannelInt .................        %12.3f", Mass->CumChannelInt*1000); 
  fprintf(Out->FilePtr, "\n      RoadInt ....................        %12.3f", (Mass->CumRoadInt  -
                                                                        Mass->CumCulvertReturnFlow) * 1000);
  fprintf(Out->FilePtr, "\n  Storage Change .................        %12.3f", (NewWaterStorage - Mass->StartWaterStorage)*1000);
  fprintf(Out->FilePtr, "\n      Initial Storage ............        %12.3f", Mass->StartWaterStorage*1000);
  fprintf(Out->FilePtr, "\n      Final Storage ..............        %12.3f", NewWaterStorage*1000);
  fprintf(Out->FilePtr, "\n          Final SWQ ..............        %12.3f", Total->Snow.Swq*1000);
  fprintf(Out->FilePtr, "\n          Final Soil Moisture ....        %12.3f", (Total->SoilWater + Total->Soil.SatFlow)*1000);
  fprintf(Out->FilePtr, "\n          Final Surface ..........        %12.3f", (Total->Soil.IExcess  + 
						                               Total->CanopyWater + Total->Soil.DetentionStorage)*1000);
  fprintf(Out->FilePtr, "\n          Final Road Surface .....        %12.3f\n", Total->Road.IExcess*1000);
  fprintf(Out->FilePtr, "\n  Mass added to glacier ..........        %12.3f\n", Total->Snow.Glacier*1000);
  fprintf(Out->FilePtr, "  ******************************************************");
  fprintf(Out->FilePtr, "\n  Mass Error (mm).................        %12.3f", MassError*1000);
	  	  
  if(Options->Sediment){
    fprintf(stderr, "\nFinal Sediment Mass Balance\n");   
    if (Options->MassWaste){     
      MWMMassError = Mass->CumMassWasting - Mass->CumSedimentToChannel - 
	Mass->CumMassDeposition;
      /* These are totals, not averages */
      fprintf(stderr, " \nTotal Mass Wasting\n");
      fprintf(stderr, " MassWasted (m3): %.2e\n", Mass->CumMassWasting);
      fprintf(stderr, " SedimentToChannel (m3): %.2e\n", 
	      Mass->CumSedimentToChannel);
      fprintf(stderr, " MassDepostion (m3): %.2e\n", Mass->CumMassDeposition);
      fprintf(stderr, " Mass Error (m3): %e\n", MWMMassError);
    }
    /* These are averages */
    if (Options->InitSedFlag){
      fprintf(stderr, " \nAverage Surface Erosion\n");
      fprintf(stderr, " Surface Erosion (mm): %.2e\n", 
	      Mass->CumSedimentErosion);
      fprintf(stderr, " Surface Erosion (kg/hectare): %.2e\n", 
	      Mass->CumSedimentErosion * PARTDENSITY * (float)MMTOM * 10000.);
    }
    
    if (Options->RoadRouting){
      fprintf(stderr, " \nBasin Average Road Surface Erosion\n");
      fprintf(stderr, " Road Surface Erosion (mm): %.2e\n", 
	      Mass->CumRoadErosion * 1000.);
      fprintf(stderr, " Road Surface Erosion (kg/hectare): %.2e\n", 
	      Mass->CumRoadErosion * PARTDENSITY * 10000.);
      fprintf(stderr, "Road Sediment to Hillslope (mm): %.2e\n", 
	      Mass->CumRoadSedHill * 1000.);

      /* roadarearatio is used to convert basin average road erosion to 
	 road erosion averaged over the road surface area */

      fprintf(stderr, " \nAverage Road Surface Erosion\n");
      fprintf(stderr, " Road Surface Erosion (mm): %.2e\n", 
	      Mass->CumRoadErosion/roadarearatio * 1000.);
      fprintf(stderr, " Road Surface Erosion (kg/hectare): %.2e\n", 
	      Mass->CumRoadErosion/roadarearatio * PARTDENSITY * 10000.);
      fprintf(stderr, "Road Sediment to Hillslope (mm): %.2e\n", 
	      Mass->CumRoadSedHill/roadarearatio * 1000.);
    }
    /* These are totals */
    SedInput = Mass->CumDebrisInflow + 
      (Mass->CumSedOverlandInflow - Mass->CumCulvertSedToChannel) + 
      Mass->CumSedOverroadInflow;

    /* NOTE: CulvertSedToChannel + CulvertReturnSedFlow = CulvertSedFlow */
    SedOutput = Mass->CumSedimentOutflow - 
      (Mass->CumCulvertSedToChannel + Mass->CumCulvertReturnSedFlow) ; 
    
    SedMassError = (Total->ChannelSedimentStorage + 
		    Total->ChannelSuspendedSediment - 
		    Mass->StartChannelSedimentStorage) + 
      SedOutput - SedInput;  
    
    fprintf(stderr, " \nChannel Erosion");
    fprintf(stderr, " \nInflow %.2e:\n", SedInput); 
    fprintf(stderr, " DebrisInflow (kg): %e\n", Mass->CumDebrisInflow); 
    fprintf(stderr, " OverlandInflow (kg): %e\n", 
	    Mass->CumSedOverlandInflow);
    fprintf(stderr, " OverroadInflow (kg): %e\n", 
	     Mass->CumSedOverroadInflow);
    
    if (Options->ChannelRouting){
      fprintf(stderr, " \nOutflow %.2e:\n",  SedOutput);
      fprintf(stderr, " SedimentOutflow (kg): %e\n", Mass->CumSedimentOutflow);
      fprintf(stderr, " CulvertReturnSedFlow (kg): %e\n", 
	      Mass->CumCulvertReturnSedFlow); 
      fprintf(stderr, " CulvertSedToChannel (kg): %e\n", 
	      Mass->CumCulvertSedToChannel);    
      
      fprintf(stderr, " \nStorage:\n");
      fprintf(stderr, " Initial Storage (kg): %e\n", 
	      Mass->StartChannelSedimentStorage); 
      fprintf(stderr, " End of Run Storage (kg): %e\n", 
	      Total->ChannelSedimentStorage + Total->ChannelSuspendedSediment); 
      fprintf(stderr, " \tFinal Bed Storage (kg): %e\n", 
	      Total->ChannelSedimentStorage); 
      fprintf(stderr, " \tFinal Suspended Sediment (kg): %e\n", 
	      Total->ChannelSuspendedSediment); 
      fprintf(stderr, " \nMass Error (kg): %e\n", SedMassError); 
    }
  }
}

