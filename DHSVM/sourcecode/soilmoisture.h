/*
 * SUMMARY:      soilmoisture.h - Header file for soil moisture calculations
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Mark Wigmosta
 * ORG:          Batelle Pacific Northwest Laboratories
 * E-MAIL:       ms_wigmosta@pnl.gov
 * ORIG-DATE:    Jul-1996
 * DESCRIPTION:  Header file for soil moisture calculations, including
 *               corrections for road and channel cuts
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: soilmoisture.h,v 1.7 2004/02/04 03:08:56 colleen Exp $     
 */

#ifndef SOILMOISTURE_H
#define SOILMOISTURE_H

#define NO_CUT -10

void AdjustStorage(int NSoilLayers, float TotalDepth, float *RootDepth,
		   float Area, float DX, float DY, float BankHeight, 
		   float *PercArea, float *Adjust, int *CutBankZone);

float CalcAvailableWater(int NRootLayers, float TotalDepth, float *RootDepth,
			 float *Moisture, float *FCap, float TableDepth,
			 float *Adjust);

float CalcTotalWater(int NSoilLayers, float TotalDepth, float *RootDepth,
		     float *Moist, float *Adjust);

void CutBankGeometry(int i, float RootDepth, float TopZone, float BankHeight,
		     float Area, float DX, float DY, float *PercArea, 
		     float *Adjust, int *CutBankZone);

void UnsaturatedFlow(int Dt, float DX, float DY, float Infiltration, 
		     float RoadbedInfiltration, float SatFlow, int NSoilLayers, 
		     float TotalDepth, float Area, float *RootDepth, float *Ks, 
		     float *PoreDist, float *Porosity, float *FCap, 
		     float *Perc, float *PercArea, float *Adjust, 
		     int CutBankZone, float BankHeight, float *TableDepth, 
		     float *Runoff, float *Moist, int RoadRouteOption,
		     int InfiltOption, float *RoadIExcess);

float WaterTableDepth(int NRootLayers, float TotalDepth, float *RootDepth,
		      float *Porosity, float *FCap, float *Adjust,
		      float *Moist);

#endif
