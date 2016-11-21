/*
 * SUMMARY:      varid.h - header for VarID.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Sat Jan 30 13:38:22 1999
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:
 * $Id: varid.h,v 1.4 2003/07/01 21:26:33 olivier Exp $     
 */

#ifndef VARID_H
#define VARID_H

#ifndef ENDOFLIST
#define ENDOFLIST -1
#endif

void GetVarAttr(MAPDUMP * DMap);
void GetVarFileLabel(int ID, char *FileLabel);
void GetVarFileName(int ID, int Layer, unsigned char Resolution,
		    char *FileName);
void GetVarFormat(int ID, char *Format);
void GetVarLongName(int ID, int Layer, char *LongName);
int GetVarNLayers(int ID, int MaxSoilLayers, int MaxVegLayers);
void GetVarName(int ID, int Layer, char *Name);
void GetVarNumberType(int ID, int *NumberType);
void GetVarUnits(int ID, char *Units);
unsigned char IsMultiLayer(int ID);
unsigned char IsValidID(int ID);

#endif
