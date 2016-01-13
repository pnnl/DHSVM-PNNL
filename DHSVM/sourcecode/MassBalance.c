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
 * FUNCTIONS:    MassBalance()
 * COMMENTS:
 * Modification made on 2012/12/31
 * $Id: MassBalance.c, v 4.0 Ning Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "constants.h"
#include "Calendar.h"

/*****************************************************************************
  Aggregate()
  
  Calculate the average values for the different fluxes and state variables
  over the basin.  Only the runoff and some of the sediment variables (as 
  noted) are calculated as a totals (i.e. runoff is total volume) instead
  of an average.  In the current implementation the local radiation
  elements are not stored for the entire area.  Therefore these components
  are aggregated in AggregateRadiation() inside MassEnergyBalance().

  The aggregated values are set to zero in the function RestAggregate,
  which is executed at the beginning of each time step.
*****************************************************************************/
void MassBalance(DATE *Current, FILES *Out, FILES *SedOut, AGGREGATED *Total,
		 WATERBALANCE *Mass, OPTIONSTRUCT * Options)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float Input;
  float MassError;		/* mass balance error m  */
  float MWMMassError;            /* mass wasting mass balance error m3  */
  float SedInput, SedOutput, SedMassError;  /* sediment mass balance variables 
					       for channel network */

  NewWaterStorage = Total->Soil.IExcess + Total->Road.IExcess + 
    Total->CanopyWater + Total->SoilWater +
    Total->Snow.Swq + Total->Soil.SatFlow + Total->Soil.DetentionStorage;

  Output = Total->ChannelInt + Total->RoadInt + Total->Evap.ETot;
  Input = Total->Precip.Precip + Total->Snow.VaporMassFlux +
    Total->Snow.CanopyVaporMassFlux + Total->CulvertReturnFlow;

  MassError = (NewWaterStorage - Mass->OldWaterStorage) + Output -
    Total->Precip.Precip - Total->Snow.VaporMassFlux -
    Total->Snow.CanopyVaporMassFlux - Total->CulvertReturnFlow;

  /* update */
  Mass->OldWaterStorage = NewWaterStorage;
  Mass->CumPrecipIn += Total->Precip.Precip;
  Mass->CumIExcess += Total->Soil.IExcess;
  Mass->CumChannelInt += Total->ChannelInt;
  Mass->CumRoadInt += Total->RoadInt;
  Mass->CumET += Total->Evap.ETot;
  Mass->CumSnowVaporFlux += Total->Snow.VaporMassFlux +
    Total->Snow.CanopyVaporMassFlux;
  Mass->CumCulvertReturnFlow += Total->CulvertReturnFlow;
  Mass->CumCulvertToChannel += Total->CulvertToChannel;
  Mass->CumRunoffToChannel += Total->RunoffToChannel;
  
  PrintDate(Current, Out->FilePtr);
  fprintf(Out->FilePtr, " %7.4f  %7.4f  %6.3f  %8.4f  %.2e  \
%.2e  %5.2f  %5.2f  %7.4f  %7.4f  %7.4f  %6.3f  %.2e  %5.2f  %.2e  %7.3f \n",
	  Total->Soil.IExcess, Total->CanopyWater, Total->SoilWater, 
	  Total->Snow.Swq, Total->Soil.SatFlow, 
	  Total->ChannelInt,  Total->RoadInt, Total->CulvertReturnFlow, Total->Evap.ETot,  
	  Total->Precip.Precip, Total->Snow.VaporMassFlux, 
	  Total->Snow.CanopyVaporMassFlux, Mass->OldWaterStorage, Total->CulvertToChannel,
	  Total->RunoffToChannel, MassError);

  if(Options->Sediment){
    /* Calculate sediment mass errors */
    MWMMassError = Total->Fine.MassWasting -Total->Fine.SedimentToChannel - 
      Total->Fine.MassDeposition ;
    
    SedInput = Total->DebrisInflow + 
      (Total->SedimentOverlandInflow - Total->CulvertSedToChannel) + 
      Total->SedimentOverroadInflow ;
    
    SedOutput = Total->SedimentOutflow  - 
      (Total->CulvertSedToChannel + Total->CulvertReturnSedFlow); 
    
    SedMassError = (Total->ChannelSedimentStorage + 
		    Total->ChannelSuspendedSediment - 
		    Mass->LastChannelSedimentStorage) + 
      SedOutput - SedInput;
    
    /* update */
    
    /* Mass Wasting - tota,l m3 */
    Mass->CumMassWasting += Total->Fine.MassWasting;
    Mass->CumSedimentToChannel += Total->Fine.SedimentToChannel;
    Mass->CumMassDeposition += Total->Fine.MassDeposition;
    
    /* Surface Erosion - ave, mm */
    Mass->CumSedimentErosion += Total->Sediment.Erosion;
    
    /* Road Erosion - ave, m */
    Mass->CumRoadErosion += Total->Road.Erosion;
    Mass->CumRoadSedHill += Total->Sediment.RoadSed;
    
    /* Channel Erosion - total, kg */
    Mass->CumDebrisInflow += Total->DebrisInflow;
    
    Mass->CumSedOverlandInflow += Total->SedimentOverlandInflow;
    Mass->CumSedOverroadInflow += Total->SedimentOverroadInflow;
    
    Mass->CumCulvertSedToChannel += Total->CulvertSedToChannel;
    Mass->CumCulvertReturnSedFlow += Total->CulvertReturnSedFlow;
    Mass->CumSedimentOutflow += Total->SedimentOutflow;
    
    Mass->LastChannelSedimentStorage = Total->ChannelSedimentStorage + 
      Total->ChannelSuspendedSediment;   
    
    PrintDate(Current, SedOut->FilePtr);
    
    fprintf(SedOut->FilePtr, " %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g \n", 
	    Total->Fine.MassWasting, Total->Fine.SedimentToChannel, 
	    Total->Fine.MassDeposition, MWMMassError, Total->Sediment.Erosion,
	    Total->Road.Erosion,Total->Sediment.RoadSed, 
	    Total->DebrisInflow, Total->SedimentOverlandInflow, 
	    Total->SedimentOverroadInflow, Total->SedimentOutflow, 
	    Total->CulvertReturnSedFlow, Total->CulvertSedToChannel,
	    Mass->LastChannelSedimentStorage, SedMassError);
  } 
}

/*        1. Total mass wasted (m3) */
/*        2. Total mass wasted delivered to channel (m3) */
/*        3. Total mass deposition (m3) */
/*        4. Total mass wasting mass error (m3) */
/*        5. Ave hillslope erosion (mm) */
/*        6. Ave road erosion (m) */
/*        7. Ave road erosion delivered to hillslope (m)*/                      
/*        8. Total debris inflow (kg) */
/*        9. Total overland inflow (kg) */
/*       10. Total overroad inflow (kg) */
/*       11. Total sediment outflow (kg) */
/*       12. Total culvert return sediment flow (kg) */
/*       13. Total culvert sediment to channel (kg) */
/*       14. Total amount of sediment stored in channels (kg) */
/*       15. Total channel erosion mass balance error for the current time step (kg) */
  
