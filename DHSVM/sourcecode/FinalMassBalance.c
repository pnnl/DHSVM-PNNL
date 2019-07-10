/*
 * SUMMARY:      MassBalance.c - calculate basin-wide mass balance
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Battelle - Pacific Northwest National Laboratory
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Oct-96
 * DESCRIPTION:  Calculate water mass balance errors
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
void FinalMassBalance(FILES *Out, AGGREGATED *Total, WATERBALANCE *Mass)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float MassError;		/* mass balance error m  */
  float Input;

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
  fprintf(stderr, "\n  Total Inflow ...................        %.3f", Input*1000);
  fprintf(stderr, "\n      Precip/Inflow ..............        %.3f", Mass->CumPrecipIn*1000); 
  fprintf(stderr, "\n      SnowVaporFlux ..............        %.3f", Mass->CumSnowVaporFlux*1000);
  fprintf(stderr, "\n  Total Outflow ..................        %.3f", Output*1000);
  fprintf(stderr, "\n      ET .........................        %.3f", Mass->CumET*1000);   
  fprintf(stderr, "\n      ChannelInt .................        %.3f", Mass->CumChannelInt*1000); 
  fprintf(stderr, "\n      RoadInt ....................        %.3f", (Mass->CumRoadInt  -
                                                                        Mass->CumCulvertReturnFlow) * 1000);
  fprintf(stderr, "\n  Storage Change .................        %.3f", (NewWaterStorage - Mass->StartWaterStorage)*1000);
  fprintf(stderr, "\n      Initial Storage ............        %.3f", Mass->StartWaterStorage*1000);
  fprintf(stderr, "\n      Final Storage ..............        %.3f", NewWaterStorage*1000);
  fprintf(stderr, "\n          Final SWQ ..............        %.3f", Total->Snow.Swq*1000);
  fprintf(stderr, "\n          Final Soil Moisture ....        %.3f", (Total->SoilWater + Total->Soil.SatFlow)*1000);
  fprintf(stderr, "\n          Final Surface ..........        %.3f", (Total->Soil.IExcess  + 
						                               Total->CanopyWater + Total->Soil.DetentionStorage)*1000);
  fprintf(stderr, "\n          Final Road Surface .....        %.3f\n", Total->Road.IExcess*1000);
  fprintf(stderr, "\n  Mass added to glacier ..........        %.3f\n", Total->Snow.Glacier*1000);
  fprintf(stderr, "  ******************************************************");
  fprintf(stderr, "\n  Mass Error (mm).................        %.3f", MassError*1000);
  
    /* Print the runoff final balance results to the output file named final.mass.balance */
  fprintf(Out->FilePtr, "\n  ********************************               Depth");
  fprintf(Out->FilePtr, "\n  Runoff Final Mass Balance                        mm"); 
  fprintf(Out->FilePtr, "\n  ********************************        ------------"); 
  fprintf(Out->FilePtr, "\n  Total Inflow ...................        %.3f", Input*1000);
  fprintf(Out->FilePtr, "\n      Precip/Inflow ..............        %.3f", Mass->CumPrecipIn*1000); 
  fprintf(Out->FilePtr, "\n      SnowVaporFlux ..............        %.3f", Mass->CumSnowVaporFlux*1000);
  fprintf(Out->FilePtr, "\n  Total Outflow ..................        %.3f", Output*1000);
  fprintf(Out->FilePtr, "\n      ET .........................        %.3f", Mass->CumET*1000);   
  fprintf(Out->FilePtr, "\n      ChannelInt .................        %.3f", Mass->CumChannelInt*1000); 
  fprintf(Out->FilePtr, "\n      RoadInt ....................        %.3f", (Mass->CumRoadInt  -
                                                                        Mass->CumCulvertReturnFlow) * 1000);
  fprintf(Out->FilePtr, "\n  Storage Change .................        %.3f", (NewWaterStorage - Mass->StartWaterStorage)*1000);
  fprintf(Out->FilePtr, "\n      Initial Storage ............        %.3f", Mass->StartWaterStorage*1000);
  fprintf(Out->FilePtr, "\n      Final Storage ..............        %.3f", NewWaterStorage*1000);
  fprintf(Out->FilePtr, "\n          Final SWQ ..............        %.3f", Total->Snow.Swq*1000);
  fprintf(Out->FilePtr, "\n          Final Soil Moisture ....        %.3f", (Total->SoilWater + Total->Soil.SatFlow)*1000);
  fprintf(Out->FilePtr, "\n          Final Surface ..........        %.3f", (Total->Soil.IExcess  + 
						                               Total->CanopyWater + Total->Soil.DetentionStorage)*1000);
  fprintf(Out->FilePtr, "\n          Final Road Surface .....        %.3f\n", Total->Road.IExcess*1000);
  fprintf(Out->FilePtr, "\n  Mass added to glacier ..........        %.3f\n", Total->Snow.Glacier*1000);
  fprintf(Out->FilePtr, "  ******************************************************");
  fprintf(Out->FilePtr, "\n  Mass Error (mm).................        %.3f", MassError*1000);
  
     /* error check: negative soil moisture and surface ponding */
  if (Total->SoilWater + Total->Soil.SatFlow < 0) {
    fprintf(stderr,
      "FINAL MASS BALANCE ERROR:  Negative soil moisture %.3f\n", (Total->SoilWater + Total->Soil.SatFlow) * 1000);
    fprintf(Out->FilePtr,
      "FINAL MASS BALANCE ERROR:  Negative soil moisture %.3f\n", (Total->SoilWater + Total->Soil.SatFlow) * 1000);
  }
  if ((Total->Soil.IExcess + Total->CanopyWater + Total->Soil.DetentionStorage)/ Input > 0.1) {
    fprintf(stderr, "FINAL MASS BALANCE ERROR:  TOO MUCH SURFACE WATER PONDING %.3f\n", 
      (Total->Soil.IExcess + Total->CanopyWater + Total->Soil.DetentionStorage) * 1000);
    fprintf(Out->FilePtr, "FINAL MASS BALANCE ERROR:  TOO MUCH SURFACE WATER PONDING %.3f\n",
      (Total->Soil.IExcess + Total->CanopyWater + Total->Soil.DetentionStorage) * 1000);
  }
	  	 
}

