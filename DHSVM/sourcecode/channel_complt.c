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
  FILE * out, *out2, *out4, *out3, *out5,*out9, *out10, *out11, 
	  *out13, *out14, *out6, *out7;

  Dt = Time->Dt;

  out = netfile->streamoutflow;
  out2 = netfile->streaminflow;
  out4 = netfile->streamISW;
  out5 = netfile->streamILW;
  out6 = netfile->streamBeam;
  out7 = netfile->streamDiffuse;
  out9 = netfile->streamVP;
  out10 = netfile->streamWND;
  out11 = netfile->streamATP;
  out13 = netfile->streamNLW;
  out14 = netfile->streamNSW;
  out3 = netfile->streamSkyView;

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
  PrintRBMStartDate(Dt, &(Time->Current), out3);
  fprintf(out3, " ");
  PrintDate(&(Time->End), out3);
  fprintf(out3, " %d", Dt/3600);
  fprintf(out3, "\n");
  PrintRBMStartDate(Dt, &(Time->Current), out4);
  fprintf(out4, " ");
  PrintDate(&(Time->End), out4);
  fprintf(out4, " %d", Dt/3600);
  fprintf(out4, "\n");
  PrintRBMStartDate(Dt, &(Time->Current), out5);
  fprintf(out5, " ");
  PrintDate(&(Time->End), out5);
  fprintf(out5, " %d", Dt/3600);
  fprintf(out5, "\n");
  PrintRBMStartDate(Dt, &(Time->Current), out6);
  fprintf(out6, " ");
  PrintDate(&(Time->End), out6);
  fprintf(out6, " %d", Dt/3600);
  fprintf(out6, "\n");
  PrintRBMStartDate(Dt, &(Time->Current), out7);
  fprintf(out7, " ");
  PrintDate(&(Time->End), out7);
  fprintf(out7, " %d", Dt/3600);
  fprintf(out7, "\n");
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
  }

  if (flag == 1) {
    fprintf(out, "                    ");
    fprintf(out2, "                    ");
	fprintf(out3, "                    ");
    fprintf(out4, "                    ");
    fprintf(out5, "                    ");
	fprintf(out6, "                    ");
    fprintf(out7, "                    ");
    fprintf(out9, "                    ");
    fprintf(out10, "                    ");
    fprintf(out11, "                    ");
    fprintf(out13, "                    ");
    fprintf(out14, "                    ");

    for (; net != NULL; net = net->next) {
	  fprintf(out, "%12d ", net->id);
      fprintf(out2, "%12d ", net->id);
	  fprintf(out3, "%12d ", net->id);
      fprintf(out4, "%8d ", net->id);
      fprintf(out5, "%8d ", net->id);

	  fprintf(out6, "%8d ", net->id);
      fprintf(out7, "%8d ", net->id);

      fprintf(out9, "%9d ", net->id);
      fprintf(out10, "%8d ", net->id);
      fprintf(out11, "%5d ", net->id);

	  fprintf(out13, "%8d ", net->id);
	  fprintf(out14, "%8d ", net->id);
    }
    fprintf(out, "\n");
    fprintf(out2, "\n");
	fprintf(out3, "\n");
    fprintf(out4, "\n");
    fprintf(out5, "\n");
	
	fprintf(out6, "\n");
    fprintf(out7, "\n");

    fprintf(out9, "\n");
    fprintf(out10, "\n");
    fprintf(out11, "\n");

	fprintf(out13, "\n");
	fprintf(out14, "\n");
  }

  Time->Current.JDay = DayOfYear(Time->Current.Year, Time->Current.Month, Time->Current.Day);
  Time->Start.JDay = DayOfYear(Time->Start.Year, Time->Start.Month, Time->Start.Day);

  if ((Time->Current.JDay>=Time->Start.JDay+1) || 
	  (Time->Current.Year>Time->Start.Year)) {
  if (fprintf(out, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_outflow: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out2, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_inflow: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out3, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_inflow: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out4, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_ISW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out5, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_ILW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out9, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_ActualVaporPressure: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out10, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_Wind: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out11, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_AirTemp: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out13, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_NetLW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out14, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_NetSW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out6, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_ILW: write error:%s", strerror(errno));
    err++;
  }
  if (fprintf(out7, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,"channel_save_ILW: write error:%s", strerror(errno));
    err++;
  }

  for (; net != NULL; net = net->next) {
    if (fprintf(out, "%12.4f ", net->outflow/Dt) == EOF) {
	  error_handler(ERRHDL_ERROR, "channel_save_outflow: write error:%s", strerror(errno));
	  err++;
	}
    if (fprintf(out2, "%12.4f ", net->inflow/Dt) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_inflow: write error:%s", strerror(errno));
      err++;
	}
    if (fprintf(out4, "%8.2f ", net->ISW ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ISW: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out5, "%8.2f ", net->ILW ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ILW: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out9, "%9.2f ", net->VP ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ActualVaporPressure: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out10, "%8.2f ", net->WND ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_Wind: write error:%s", strerror(errno));
      err++;
    }
    if (fprintf(out11, "%5.2f ", net->ATP ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_AirTemp: write error:%s", strerror(errno));
      err++;
    }
	if (fprintf(out13, "%8.2f ", net->NLW) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_NetLW: write error:%s", strerror(errno));
      err++;
    }
	if (fprintf(out14, "%8.2f ", net->NSW) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_NetSW: write error:%s", strerror(errno));
      err++;
    }
	if (fprintf(out6, "%8.2f ", net->Beam ) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ISW: write error:%s", strerror(errno));
      err++;
    }
	if (fprintf(out7, "%8.2f ", net->Diffuse) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ILW: write error:%s", strerror(errno));
      err++;
    }
	if (fprintf(out3, "%8.2f ", net->skyview) == EOF) {
      error_handler(ERRHDL_ERROR, "channel_save_ILW: write error:%s", strerror(errno));
      err++;
    }
  }
  fprintf(out, "\n");
  fprintf(out2, "\n");
  fprintf(out3, "\n");
  fprintf(out4, "\n");
  fprintf(out5, "\n");

  fprintf(out6, "\n");
  fprintf(out7, "\n");

  fprintf(out9, "\n");
  fprintf(out10, "\n");
  fprintf(out11, "\n");

  fprintf(out13, "\n");
  fprintf(out14, "\n");
  }

  return (err);
}


