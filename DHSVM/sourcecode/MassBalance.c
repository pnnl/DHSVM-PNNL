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
   MassBalance()

   Calculate the average values for the different fluxes and state variables
   over the basin.
   In the current implementation the local radiation
   elements are not stored for the entire area.  Therefore these components
   are aggregated in AggregateRadiation() inside MassEnergyBalance().

   The aggregated values are set to zero in the function RestAggregate,
   which is executed at the beginning of each time step.
 *****************************************************************************/
void MassBalance(DATE *Current, DATE *Start, FILES *Out, AGGREGATED *Total, WATERBALANCE *Mass)
{
  float NewWaterStorage;	/* water storage at the end of the time step */
  float Output;			/* total water flux leaving the basin;  */
  float Input;
  float MassError;		/* mass balance error m  */

  NewWaterStorage = Total->Soil.IExcess + Total->Road.IExcess +
    Total->CanopyWater + Total->SoilWater +
    Total->Snow.Swq + Total->Soil.SatFlow + Total->Soil.DetentionStorage +
    Total->Snow.Iwq + Total->Snow.IceRemoved;

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

  if (IsEqualTime(Current, Start)) {
    fprintf(Out->FilePtr, "         Date        ");
    fprintf(Out->FilePtr, " Precip(m) ");
    fprintf(Out->FilePtr, " Snow(m) ");
    fprintf(Out->FilePtr, " IExcess(m) ");
    fprintf(Out->FilePtr, " Swq   Melt Iwq GlMelt IceRemoved ");
    fprintf(Out->FilePtr, " TotalET ");   /* total evapotranspiration*/
    fprintf(Out->FilePtr, " CanopyInt ");   /* canopy intercepted rain + snow*/
    fprintf(Out->FilePtr, " TotSoilMoist ");
    fprintf(Out->FilePtr, " SatFlow ");
    fprintf(Out->FilePtr, " SnowVaporFlux ");
    fprintf(Out->FilePtr, " ChannelInt RoadInt CulvertInt"),
      fprintf(Out->FilePtr, " PixelShortIn PixelNetShort NetShort.Layer1 NetShort.Layer2 PixelNetRadiation Tair Error");
    fprintf(Out->FilePtr, "\n");
  }
  PrintDate(Current, Out->FilePtr);
  fprintf(Out->FilePtr, " %g  %g  %g  %g  %g  %g  %g  %g  %g  %g %g\
      %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g  %g \n",
    Total->Precip.Precip, Total->Precip.SnowFall, Total->Soil.IExcess,
    Total->Snow.Swq, Total->Snow.Melt, Total->Snow.Iwq, Total->Snow.GlMelt, Total->Snow.IceRemoved, Total->Evap.ETot,
    Total->CanopyWater, Total->SoilWater, Total->Soil.SatFlow, Total->Snow.VaporMassFlux,
    Total->ChannelInt, Total->RoadInt, Total->CulvertToChannel,
    Total->Rad.BeamIn + Total->Rad.DiffuseIn, Total->Rad.PixelNetShort,
    Total->Rad.NetShort[0], Total->Rad.NetShort[1], Total->NetRad, Total->Rad.Tair, MassError);
}
