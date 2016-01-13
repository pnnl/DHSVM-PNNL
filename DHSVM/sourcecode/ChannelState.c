/*
 * SUMMARY:      ChannelState.c - Read and Store the channel state
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Tue Jan  5 16:52:42 1999
 * DESCRIPTION:  Store the state of the channel.  The channel state file
                 contains two columns.  The first column contains the unique
		 channel ID's, the second the storage in the segment in m3.
		 This is the content of the fields  
		   _channel_rec_
		     SegmentID id
		     float storage
 * DESCRIP-END.
 * FUNCTIONS:    ReadChannelState()
                 StoreChannelState()
 * COMMENTS:
 * $Id: ChannelState.c,v 1.5 2006/10/03 22:50:22 nathalie Exp $     
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h> 
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "sizeofnt.h"
#include "channel.h"

typedef struct _RECORDSTRUCT {
  SegmentID id;
  float storage;
} RECORDSTRUCT;

int CompareRecord(const void *record1, const void *record2);
int CompareRecordID(const void *key, const void *record);

/*****************************************************************************
  ReadChannelState()

  Read the state of the channel from a previous run.  Currently just read an
  ASCII file, with the unique channel IDs in the first column and the amount
  of storage in the second column (m3)
*****************************************************************************/
void ReadChannelState(char *Path, DATE * Now, Channel * Head)
{
  char InFileName[BUFSIZ + 1] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *InFile = NULL;
  int i = 0;
  int NLines = 0;
  int max_seg = 0;
  RECORDSTRUCT *Match = NULL;
  RECORDSTRUCT *Record = NULL;

  /* Re-create the storage file name and open it */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day,
	  Now->Year, Now->Hour, Now->Min, Now->Sec);
  sprintf(InFileName, "%sChannel.State.%s", Path, Str);
  OpenFile(&InFile, InFileName, "r", TRUE);
  NLines = CountLines(InFile);
  rewind(InFile);

  /* Allocate memory and read the file */
  Record = (RECORDSTRUCT *) calloc(NLines, sizeof(RECORDSTRUCT));
  if (Record == NULL)
    ReportError("ReadChannelState", 1);
  for (i = 0; i < NLines; i++)
    fscanf(InFile, "%hu %f", &(Record[i].id), &(Record[i].storage));
  qsort(Record, NLines, sizeof(RECORDSTRUCT), CompareRecord);

  /* Assign the storages to the correct IDs */
  Current = Head;
  while (Current) {
    Match = bsearch(&(Current->id), Record, NLines, sizeof(RECORDSTRUCT),
		    CompareRecordID);
    if (Current->id > max_seg)
      max_seg = Current->id;
    if (Match == NULL)
      ReportError("ReadChannelState", 55);
    Current->storage = Match->storage;
    Current = Current->next;
  }

  /* Clean up */
  if (Record)
    free(Record);
  fclose(InFile);
}

/*****************************************************************************
  StoreChannelState()

  Store the current state of the channel, i.e. the storage in each channel 
  segment.

*****************************************************************************/
void StoreChannelState(char *Path, DATE * Now, Channel * Head)
{
  char OutFileName[BUFSIZ + 1] = "";
  char Str[BUFSIZ + 1] = "";
  Channel *Current = NULL;
  FILE *OutFile = NULL;

  printf("storing channel state \n");
  /* Create storage file */
  sprintf(Str, "%02d.%02d.%04d.%02d.%02d.%02d", Now->Month, Now->Day,
	  Now->Year, Now->Hour, Now->Min, Now->Sec);
  sprintf(OutFileName, "%sChannel.State.%s", Path, Str);
  OpenFile(&OutFile, OutFileName, "w", TRUE);

  /* Store data */
  Current = Head;
  while (Current) {
    fprintf(OutFile, "%12hu ", Current->id);
    fprintf(OutFile, "%12g\n", Current->storage);
    Current = Current->next;
  }

  /* Close file */
  fclose(OutFile);
}

/*****************************************************************************
  CompareRecord()

  Compare two RECORDSTRUCT elements for qsort
*****************************************************************************/
int CompareRecord(const void *record1, const void *record2)
{
  RECORDSTRUCT *x = NULL;
  RECORDSTRUCT *y = NULL;

  x = (RECORDSTRUCT *) record1;
  y = (RECORDSTRUCT *) record2;

  return (int) x->id - y->id;
}

/*****************************************************************************
  CompareRecordID()

  Compare RECORDSTRUCT element with an ID to see if the RECORDSTRUCT has the
  right ID for bsearch
*****************************************************************************/
int CompareRecordID(const void *key, const void *record)
{
  SegmentID *x = NULL;
  RECORDSTRUCT *y = NULL;

  x = (SegmentID *) key;
  y = (RECORDSTRUCT *) record;

  return (int) (*x - y->id);
}
