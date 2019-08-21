/*
 * SUMMARY:      channel_complt.c - Calculate mass and energy balance
 * USAGE:        Part of DHSVM-RBM
 *               This function is used to output files for John's RBM model
 *
 * AUTHOR:       Nathalie Voisin
 * ORG:          University of Washington, Department of Civil Engineering
 * ORIG-DATE:    Aug-10
 * DESCRIPTION:  Calculate mass and energy balance at each pixel
 * DESCRIP-END.
 * FUNCTIONS:    channel_save_outflow_text_cplmt()
                 channel_save_outflow_cplmt()
 * Modification 
 * $Id: channel_complt.c, v 3.2  2013/04/23   Ning Exp $    
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "errorhandler.h"
//#include "channel.h"
#include "functions.h"
#include "constants.h"
#include "tableio.h"
#include "settings.h"

/* -------------------------------------------------------------
   ---------------------- Channel Functions --------------------
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   channel_save_outflow_cplmt
   Saves the channel outflow using a text string as the time field
   for John's RBM model
   ------------------------------------------------------------- */
int
channel_save_outflow_text_cplmt
(TIMESTRUCT *Time, char *tstring, Channel * net, CHANNEL * netfile , int flag)
{
  int err = 0;
  int Dt;
  FILE * out, *out2, *out9, *out10, *out11, *out13, *out14, *out15;

  Dt = Time->Dt;

  out = netfile->streamoutflow;
  out2 = netfile->streaminflow;
  out9 = netfile->streamVP;
  out10 = netfile->streamWND;
  out11 = netfile->streamATP;
  out13 = netfile->streamNLW;
  out14 = netfile->streamNSW;
  out15 = netfile->streamMelt;                              

  /* print the start date and end date. Note that the true start date is a time step behind the 
  user specified date when the model outputs data */
  if (flag == 1) {
    PrintRBMStartDate(Dt, &(Time->Current), out);
    fprintf(out, " ");
    PrintDate(&(Time->End), out);
    fprintf(out, " %d", Dt/3600);
    fprintf(out, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out2);
    fprintf(out2, " ");
    PrintDate(&(Time->End), out2);
    fprintf(out2, " %d", Dt/3600);
    fprintf(out2, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out9);
    fprintf(out9, " ");
    PrintDate(&(Time->End), out9);
    fprintf(out9, " %d", Dt/3600);
    fprintf(out9, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out10);
    fprintf(out10, " ");
    PrintDate(&(Time->End), out10);
    fprintf(out10, " %d", Dt/3600);
    fprintf(out10, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out11);
    fprintf(out11, " ");
    PrintDate(&(Time->End), out11);
    fprintf(out11, " %d", Dt/3600);
    fprintf(out11, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out13);
    fprintf(out13, " ");
    PrintDate(&(Time->End), out13);
    fprintf(out13, " %d", Dt/3600);
    fprintf(out13, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out14);
    fprintf(out14, " ");
    PrintDate(&(Time->End), out14);
    fprintf(out14, " %d", Dt/3600);
    fprintf(out14, "\n");
    PrintRBMStartDate(Dt, &(Time->Current), out15);
    fprintf(out15, " ");
    PrintDate(&(Time->End), out15);
    fprintf(out15, " %d", Dt / 3600);
    fprintf(out15, "\n");                       
  }

  if (flag == 1) {
    fprintf(out, "Date ");
    fprintf(out2, "Date ");
    fprintf(out9, "Date ");
    fprintf(out10, "Date ");
    fprintf(out11, "Date ");
    fprintf(out13, "Date ");
    fprintf(out14, "Date ");
    fprintf(out15, "Date ");                                       

    for (; net != NULL; net = net->next) {
	    fprintf(out, "%d ", net->id);
      fprintf(out2, "%d ", net->id);
      fprintf(out9, "%d ", net->id);
      fprintf(out10, "%d ", net->id);
      fprintf(out11, "%d ", net->id);
      fprintf(out13, "%d ", net->id);
      fprintf(out14, "%d ", net->id);
      fprintf(out15, "%d ", net->id);                                  
    }
    fprintf(out, "\n");
    fprintf(out2, "\n");
    fprintf(out9, "\n");
    fprintf(out10, "\n");
    fprintf(out11, "\n");
    fprintf(out13, "\n");
    fprintf(out14, "\n");
    fprintf(out15, "\n");                         
  }

  Time->Current.JDay = DayOfYear(Time->Current.Year, Time->Current.Month, Time->Current.Day);
  Time->Start.JDay = DayOfYear(Time->Start.Year, Time->Start.Month, Time->Start.Day);

  if ((Time->Current.JDay>=Time->Start.JDay+1) || 
	  (Time->Current.Year>Time->Start.Year)) {
  if (fprintf(out, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_outflow: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out2, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_inflow: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out9, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_ActualVaporPressure: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out10, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_Wind: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out11, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_AirTemp: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out13, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_NetLW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out14, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_NetSW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out15, "%s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR, "channel_save_NetSW: write error:%s", strerror(errno));
    err++;
  }                                                

  for (; net != NULL; net = net->next) {
    if (fprintf(out, "%.6f ", net->outflow/Dt) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_outflow: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out2, "%.6f ", net->inflow/Dt) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_inflow: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out9, "%.2f ", net->VP ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ActualVaporPressure: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out10, "%.2f ", net->WND ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_Wind: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out11, "%.2f ", net->ATP ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_AirTemp: write error:%s", strerror(errno));
      err++;
    }
	  if (fprintf(out13, "%.2f ", net->NLW) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_NetLW: write error:%s", strerror(errno));
      err++;
    }
	  if (fprintf(out14, "%.2f ", net->NSW) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_NetSW: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out15, "%.6f ", net->melt/Dt) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_Melt: write error:%s", strerror(errno));
      err++;
    }                                                       
  }
  fprintf(out, "\n");
  fprintf(out2, "\n");
  fprintf(out9, "\n");
  fprintf(out10, "\n");
  fprintf(out11, "\n");
  fprintf(out13, "\n");
  fprintf(out14, "\n");
  fprintf(out15, "\n");                       
  }

  return (err);
}


