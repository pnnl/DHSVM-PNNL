/* -------------------------------------------------------------
   file: channel_grid.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created January  5, 1996 by  William A Perkins
   $Id: channel_grid.c,v 1.12 2004/10/07 20:51:08 jlanini Exp $
   ------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif
#include <errno.h>
#include <string.h>

#include "channel_grid.h"
#include "tableio.h"
#include "errorhandler.h"
#include "settings.h"
#include "data.h"
#include "DHSVMChannel.h"
#include "constants.h"

/* -------------------------------------------------------------
   local function prototype
   ------------------------------------------------------------- */
static ChannelMapRec *alloc_channel_map_record(void);
static ChannelMapPtr **channel_grid_create_map(int cols, int rows);
Channel *Find_First_Segment(ChannelMapPtr ** map, int col, int row, float SlopeAspect, 
			    char *Continue);
Channel *Find_Next_Segment(ChannelMapPtr ** map, int curr_col, int curr_row, int next_col, 
			   int next_row, int CurrentID, int NextID, char *Continue, 
			   float *SedimentToChannel, float *SedimentMass);
char channel_grid_has_intersection(ChannelMapPtr ** map, int Currid, int Nextid, int row, 
				   int col, int Flag);

/* -------------------------------------------------------------
   local module variables
   ------------------------------------------------------------- */
static int channel_grid_cols = 0;
static int channel_grid_rows = 0;
static char channel_grid_initialized = FALSE;

/* -------------------------------------------------------------
   RouteDebrisFlow
   ------------------------------------------------------------- */
void RouteDebrisFlow(float *SedimentToChannel, int y, int x, float SlopeAspect, 
		     CHANNEL *ChannelData, MAPSIZE *Map)
{
  Channel *CurrentSeg;
  char Continue, Match;
  int SearchRadius, Flag;
  int i, j, curr_inti, curr_intj,next_inti,next_intj;
  int MaxRadius = 200;
  float SedimentMass[NSEDSIZES];

  /* Initialize */
  for(i = 0; i < NSEDSIZES; i++){
    SedimentMass[i] = 0.;
  }

  /* Find pointer to segment where debris flow enters channel network.
     If multiple segments exist in the current gridcell, the debris flow
     enters the channel with aspect closest to the aspect of the debris flow. */

  CurrentSeg = Find_First_Segment(ChannelData->stream_map, x, y, SlopeAspect, &Continue);

  /* Debris flow hits head wall; all sediment is deposited in channel. */
  if(Continue == FALSE) {
    CurrentSeg->sediment.tempvol += *SedimentToChannel;
    *SedimentToChannel = 0.0;
  }

  // fprintf(stdout, "Debris flow enters network at segment %d, row=%d, col=%d\n", CurrentSeg->id, y, x);
  
  /* Continue until basin mouth, or debris flow stops. */
  while (CurrentSeg->outlet != NULL && Continue ) {

   // fprintf(stdout, "Next segment equals %d\n", CurrentSeg->outlet->id);

    /* Find row and column of next intersection; first check current cell. */
    Flag = 0;
    Match = FALSE;
    Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
    					  CurrentSeg->outlet->id, y, x, Flag);

    if(Match==TRUE) {
      curr_inti = y;
      curr_intj = x;
      next_inti = y;
      next_intj = x;
    }
      
    /* Find row and column of next intersection; search radially outward. */
    SearchRadius = 1;
    while(!Match && Flag!=2) {
      for(i=y-SearchRadius; i<= y+SearchRadius; i++) {
	 for(j=x-SearchRadius; j<= x+SearchRadius; j++) {
	   if(i>=0 && j>=0 && i< Map->NY && j<Map->NX) {
	     if(i == y-SearchRadius || i == y+SearchRadius) {
	       Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
						     CurrentSeg->outlet->id, i, j, Flag);
	       //fprintf(stdout, "Looking for intersection in row=%d, col=%d\n",i,j);
	     }
	     else if(j == x-SearchRadius || j ==  x+SearchRadius) {
	       Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
						     CurrentSeg->outlet->id, i, j, Flag);
	       //fprintf(stdout, "Looking for intersection in row=%d, col=%d\n",i,j);
	     }
	   }
	   if(Match == TRUE && Flag==0) {
	     curr_inti = i;
	     curr_intj = j;
	     next_inti = i;
	     next_intj = j;
	     Flag = 1;
	   }
	 }
      }
      if(Flag==1)
	Match = TRUE;
      SearchRadius += 1;
      if(SearchRadius>MaxRadius)  /*if cannot find intersection, just find
				    nearest cell that contains the next segment*/
	Flag=2;      
    }

    /*Repeat the search if MaxRadius was exceeded */
    SearchRadius = 1;
    if( Flag == 2 ) { 
      /*first find the nearest cell in next segment*/
      while(!Match) {
	for(i=y-SearchRadius; i<= y+SearchRadius; i++) {
	  for(j=x-SearchRadius; j<= x+SearchRadius; j++) {
	    if(i>=0 && j>=0 && i< Map->NY && j<Map->NX) {
	      if(i == y-SearchRadius || i == y+SearchRadius) {
		Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
						      CurrentSeg->outlet->id, i, j, Flag);
		//fprintf(stdout, "Looking for intersection in row=%d, col=%d\n",i,j);
	      }
	      else if(j == x-SearchRadius || j ==  x+SearchRadius) {
		Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
						      CurrentSeg->outlet->id, i, j, Flag);
		//fprintf(stdout, "Looking for intersection in row=%d, col=%d\n",i,j);
	      }
	    }
	    if(Match == TRUE && Flag==2) {
	      next_inti = i;
	      next_intj = j;
	      Flag = 1;
	    }
	  }
	}
	if(Flag==1)
	  Match = TRUE;
	SearchRadius += 1;
      }
      /*now find the nearest cell to that cell in the current segment*/
      Flag=3;
      SearchRadius = 1;
      Match = FALSE;
      while(!Match) {
	for(i=next_inti-SearchRadius; i<= next_inti+SearchRadius; i++) {
	  for(j=next_intj-SearchRadius; j<= next_intj+SearchRadius; j++) {
	    if(i>=0 && j>=0 && i< Map->NY && j<Map->NX) {
	      if(i == next_inti-SearchRadius || i == next_inti+SearchRadius) {
		Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
						      CurrentSeg->outlet->id, i, j, Flag);
		//fprintf(stdout, "Looking for intersection in row=%d, col=%d\n",i,j);
	      }
	      else if(j == next_intj-SearchRadius || j ==  next_intj+SearchRadius) {
		Match = channel_grid_has_intersection(ChannelData->stream_map, CurrentSeg->id, 
						      CurrentSeg->outlet->id, i, j, Flag);
		//fprintf(stdout, "Looking for intersection in row=%d, col=%d\n",i,j);
	      }
	    }
	    if(Match == TRUE && Flag==3) {
	      curr_inti = i;
	      curr_intj = j;
	      Flag = 1;
	    }
	  }
	}
	if(Flag==1)
	  Match = TRUE;
	SearchRadius += 1;
      }
    }

//  fprintf(stdout, "Intersection found at row %d, col %d for current seg, at row %d, col %d for next segment\n", curr_inti, curr_intj, next_inti, next_intj);
    
    /* Now have location of intersection, check channel aspect. */
    CurrentSeg = Find_Next_Segment(ChannelData->stream_map, curr_intj, curr_inti, 
				   next_intj, next_inti, CurrentSeg->id, 
				   CurrentSeg->outlet->id, &Continue, SedimentToChannel,
				   SedimentMass);
  }
   /*  fprintf(stdout, "Debris flow stopped, Segment = %d\n", CurrentSeg->id); */
}

/* -------------------------------------------------------------
   Find_Next_Segment
   ------------------------------------------------------------- */

Channel *Find_Next_Segment(ChannelMapPtr ** map, int curr_col, int curr_row, 
			   int next_col, int next_row,int CurrentID, int NextID, 
			   char *Continue, float *SedimentToChannel, float *SedimentMass)
{
  ChannelMapPtr curr_cell = map[curr_col][curr_row];
  ChannelMapPtr next_cell = map[next_col][next_row];
  float test;
  int i;
  float CurrentAspect, NextAspect;
  Channel *CurrPtr;
  Channel *NextPtr;
 

  while (curr_cell != NULL) {
    if(curr_cell->channel->id == CurrentID) {
      CurrentAspect = curr_cell->aspect;
      CurrPtr = curr_cell->channel;
      NextPtr = curr_cell->channel->outlet;
    }
    curr_cell = curr_cell->next;
  }
  while (next_cell != NULL) {
    if(next_cell->channel->id == NextID) {
      NextAspect = next_cell->aspect;
    }
    next_cell = next_cell->next;
  }

  test = fabs(CurrentAspect - NextAspect);
  if(test > PI) {
    if(CurrentAspect < NextAspect)
      test = fabs(CurrentAspect - (NextAspect - 2.*PI));
    else
      test = fabs(NextAspect - (CurrentAspect - 2.*PI));
  }
  if(test < 0. || test > PI) {
    printf("Problem in Find_Next_Segment\n");
    exit(0);
  }

  if(test <= 70.*PI/180. && NextPtr->slope > 0.061) {    
    *Continue = TRUE;
    *SedimentToChannel += CurrPtr->sediment.tempvol;
    CurrPtr->sediment.tempvol = 0.0;
    for(i = 0; i < NSEDSIZES; i++){
       SedimentMass[i] += CurrPtr->sediment.tempmass[i];
       CurrPtr->sediment.tempmass[i] = 0.;
     }
   //  printf(" Movin sediment from %d\n",CurrPtr->id);
  }
  else {
    *Continue = FALSE;
    *SedimentToChannel += CurrPtr->sediment.tempvol;
    CurrPtr->sediment.tempvol = 0.;
    for(i = 0; i < NSEDSIZES; i++){
       SedimentMass[i] += CurrPtr->sediment.tempmass[i];
       CurrPtr->sediment.tempmass[i] = 0.;
     }
  //   printf("2Movin sediment from %d\n",CurrPtr->id);
    if(test > 70.*PI/180.) { 
      NextPtr->sediment.tempvol += *SedimentToChannel/2.;
      CurrPtr->sediment.tempvol = *SedimentToChannel/2.;
      *SedimentToChannel = 0.0;
      for(i = 0; i < NSEDSIZES; i++){
 	NextPtr->sediment.tempmass[i] += SedimentMass[i]/2.;
 	CurrPtr->sediment.tempmass[i] = SedimentMass[i]/2.;
 	SedimentMass[i] = 0.;
      }
  //     printf(" Dumpin sediment in %d %d \n",NextPtr->id,CurrPtr->id);
    }
    else {
      NextPtr->sediment.tempvol += *SedimentToChannel;
      *SedimentToChannel = 0.0;
      for(i = 0; i < NSEDSIZES; i++){
 	NextPtr->sediment.tempmass[i] += SedimentMass[i];
 	SedimentMass[i] = 0.;
      }
  //     printf("2Dumpin sediment in %d \n",NextPtr->id);
    }
  }

  return (NextPtr);
}

/* -------------------------------------------------------------
   Find_First_Segment
   ------------------------------------------------------------- */

Channel *Find_First_Segment(ChannelMapPtr ** map, int col, int row, float SlopeAspect, char *Continue)
{
  ChannelMapPtr cell = map[col][row];
  float test;
  float DeltaAspect;
  Channel *Ptr;
  
  DeltaAspect = 2.*PI;

  while (cell != NULL) {
    test = fabs(SlopeAspect - cell->aspect);

    if(test > PI) {
      if(SlopeAspect < cell->aspect)
	test = fabs(SlopeAspect - (cell->aspect - 2.*PI));
      else
	test = fabs(cell->aspect - (SlopeAspect - 2.*PI));
    }
    if(test < 0. || test > PI) {
      printf("Problem in Find_First_Segment\n");
      exit(0);
    }
    
    if( test < DeltaAspect) {
      Ptr = cell->channel;
      DeltaAspect = test;
    }
    cell = cell->next;
  }
  if(DeltaAspect <= 70.*PI/180.)
    *Continue = TRUE;
  else
    *Continue = FALSE;

  return (Ptr);
}

/* -------------------------------------------------------------
   alloc_channel_map_record
   ------------------------------------------------------------- */
static ChannelMapRec *alloc_channel_map_record(void)
{
  ChannelMapRec *p;
  if ((p = (ChannelMapRec *) malloc(sizeof(ChannelMapRec))) == NULL) {
    error_handler(ERRHDL_FATAL,
		  "alloc_channel_map_record: %s", strerror(errno));
  }
  p->length = 0.0;
  p->aspect = 0.0;
  p->sink = FALSE;
  p->channel = NULL;
  p->next = NULL;
  return (p);
}

/* -------------------------------------------------------------
   channel_grid_create_map
   ------------------------------------------------------------- */
static ChannelMapPtr **channel_grid_create_map(int cols, int rows)
{
  ChannelMapPtr **map;
  int row, col;
  ChannelMapPtr *junk;

  if ((map = (ChannelMapPtr **) malloc(cols * sizeof(ChannelMapPtr *))) == NULL) {
    error_handler(ERRHDL_FATAL, "channel_grid_create_map: malloc failed: %s",
		  strerror(errno));
  }
  if ((junk =
       (ChannelMapPtr *) malloc(rows * cols * sizeof(ChannelMapPtr))) == NULL) {
    free(map);
    error_handler(ERRHDL_FATAL,
		  "channel_grid_create_map: malloc failed: %s",
		  strerror(errno));
  }

  for (col = 0; col < cols; col++) {
    if (col == 0) {
      map[col] = junk;
    }
    else {
      map[col] = &(map[0][col * rows]);
    }
    for (row = 0; row < rows; row++) {
      map[col][row] = NULL;
    }
  }
  return (map);
}

/* -------------------------------------------------------------
   free_channel_map_record
   ------------------------------------------------------------- */
static void free_channel_map_record(ChannelMapRec * cell)
{
  if (cell->next != NULL) {
    free_channel_map_record(cell->next);
  }
  free(cell);
}

/* -------------------------------------------------------------
   channel_grid_free_map
   ------------------------------------------------------------- */
void channel_grid_free_map(ChannelMapPtr ** map)
{
  int c, r;
  for (c = 0; c < channel_grid_cols; c++) {
    for (r = 0; r < channel_grid_rows; r++) {
      if (map[c][r] != NULL) {
	free_channel_map_record(map[c][r]);
      }
    }
  }
  free(map[0]);
  free(map);
}

/* -------------------------------------------------------------
   ------------------- Input Functions -------------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   channel_grid_read_map
   ------------------------------------------------------------- */
ChannelMapPtr **channel_grid_read_map(Channel * net, const char *file,
				      SOILPIX ** SoilMap)
{
  ChannelMapPtr **map;
  static const int fields = 8;
  static char *sink_words[2] = {
    "SINK", "\n"
  };
  static TableField map_fields[8] = {
    {"Column", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
    {"Row", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
    {"Segment ID", TABLE_INTEGER, TRUE, FALSE, {0.0}, "", NULL},
    {"Segment Length", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Cut Height", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Cut Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Segment Azimuth", TABLE_REAL, FALSE, FALSE, {0.0}, "", NULL},
    {"Sink?", TABLE_WORD, FALSE, FALSE, {0.0}, "", sink_words}
  };
  int done, err = 0;

  if (!channel_grid_initialized) {
    error_handler(ERRHDL_ERROR,
		  "channel_grid_read_map: channel_grid module not initialized");
    return NULL;
  }

  error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
		  "channel.grid_read_map: unable to read file \"%s\"", file);
    return NULL;
  }

  map = channel_grid_create_map(channel_grid_cols, channel_grid_rows);

  done = FALSE;
  while (!done) {
    int i;
    int row = 0, col = 0;
    int rec_err = 0;
    ChannelMapPtr cell;

    done = (table_get_fields(fields, map_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++) {
	if (map_fields[i].read)
	  break;
      }
      if (i >= fields)
	continue;
    }

    if (map_fields[0].read) {
      if (map_fields[0].value.integer < 0 ||
	  map_fields[0].value.integer >= channel_grid_cols) {
	rec_err++;
      }
      else {
	col = map_fields[0].value.integer;
      }
    }
    else {
      rec_err++;
    }
    if (map_fields[1].read) {
      if (map_fields[1].value.integer < 0 ||
	  map_fields[1].value.integer >= channel_grid_rows) {
	rec_err++;
      }
      else {
	row = map_fields[1].value.integer;
      }
    }
    else {
      rec_err++;
    }

    if (rec_err) {
      error_handler(ERRHDL_ERROR,
		    "%s: line %d: bad coordinates", file, table_lineno());
      err++;
      continue;
    }

    if (map[col][row] != NULL) {
      cell = map[col][row];
      while (cell->next != NULL)
	cell = cell->next;
      cell->next = alloc_channel_map_record();
      cell = cell->next;
    }
    else {
      map[col][row] = alloc_channel_map_record();
      cell = map[col][row];
    }

    for (i = 2; i < fields; i++) {
      if (map_fields[i].read) {
	switch (i) {
	case 2:
	  if ((cell->channel =
	       channel_find_segment(net,
				    map_fields[i].value.integer)) == NULL) {
	    error_handler(ERRHDL_ERROR,
			  "%s, line %d: unable to locate segment %d", file,
			  table_lineno(), map_fields[i].value.integer);
	    err++;
	  }
	  break;
	case 3:
	  cell->length = map_fields[i].value.real;
	  if (cell->length < 0.0) {
	    error_handler(ERRHDL_ERROR,
			  "%s, line %d: bad length", file, table_lineno());
	    err++;
	  }
	  break;
	case 4:
	  cell->cut_height = map_fields[i].value.real;
	  if (cell->cut_height > SoilMap[row][col].Depth) {
	    printf("warning overriding cut depths with 0.95 soil depth \n");
	    cell->cut_height = SoilMap[row][col].Depth*0.95;
	  }
	  if (cell->cut_height < 0.0
	      || cell->cut_height > SoilMap[row][col].Depth) {
	    error_handler(ERRHDL_ERROR, "%s, line %d: bad cut_depth", file,
			  table_lineno());
	    err++;
	  }
	  break;
	case 5:
	  cell->cut_width = map_fields[i].value.real;
	  if (cell->cut_width < 0.0) {
	    error_handler(ERRHDL_ERROR,
			  "%s, line %d: bad cut_width", file, table_lineno());
	    err++;
	  }
	  break;
	case 6:
	  /* road aspect is read in degrees and
	     stored in radians */
	  cell->aspect = map_fields[i].value.real * PI / 180.0;
	  break;
	case 7:
	  cell->sink = TRUE;
	  break;
	default:
	  error_handler(ERRHDL_FATAL,
			"channel_grid_read_map: this should not happen");
	  break;
	}
      }
    }

  }

  table_errors += err;
  error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: %s: %d errors, %d warnings",
		file, table_errors, table_warnings);

  table_close();

  error_handler(ERRHDL_STATUS,
		"channel_grid_read_map: done reading file \"%s\"", file);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
		  "channel_grid_read_map: %s: too many errors", file);
    channel_grid_free_map(map);
    map = NULL;
  }

  return (map);
}

/* -------------------------------------------------------------
   ---------------------- Query Functions ---------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   channel_grid_has_channel
   ------------------------------------------------------------- */
int channel_grid_has_channel(ChannelMapPtr ** map, int col, int row)
{
  if (map != NULL)
    return (map[col][row] != NULL);
  else
    return FALSE;
}

/* -------------------------------------------------------------
   channel_grid_has_sink
   ------------------------------------------------------------- */
int channel_grid_has_sink(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  char test = FALSE;

  while (cell != NULL) {
    test = (test || cell->sink);
    cell = cell->next;
  }
  return (test);
}

/* -------------------------------------------------------------
   channel_grid_has_intersection
   ------------------------------------------------------------- */
char channel_grid_has_intersection(ChannelMapPtr ** map, int Currid, int Nextid,  int row, int col, int Flag)
{
  ChannelMapPtr cell = map[col][row];
  char Next = FALSE;
  char Current = FALSE;
  char Intersection = FALSE;
  
  if( Flag < 2) { /* search for intersecting cells given Currid and Nextid
		      occur in the cell of intersection */
    while (cell != NULL) {
      if(cell->channel->id == Currid)
	Current = TRUE;
      if(cell->channel->id == Nextid)
	Next = TRUE;
      cell = cell->next;
    }
    
    if(Current && Next)
      Intersection = TRUE;
  }

  else if(Flag == 2) { /* if Flag == 2 , only search for the nearest cell that has id of next
	    segment */
    while (cell != NULL) {
      if(cell->channel->id == Nextid)
	Next = TRUE;
      cell = cell->next;
    }
    
    if(Next)
      Intersection = TRUE;
  }
  else { /* if Flag == 3 , only search for the nearest cell that has id of current
	    segment */
    while (cell != NULL) {
      if(cell->channel->id == Currid)
	Current = TRUE;
      cell = cell->next;
    }
    
    if(Current)
      Intersection = TRUE;
  }

  return (Intersection);
}

/* -------------------------------------------------------------
   channel_grid_cell_length
   returns the total length of channel(s) in the cell.
   ------------------------------------------------------------- */
double channel_grid_cell_length(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double len = 0.0;

  while (cell != NULL) {
    len += cell->length;
    cell = cell->next;
  }
  return len;
}

/* -------------------------------------------------------------
   channel_grid_cell_width
   returns a length-weighted average of the channel widths in the cell
   ------------------------------------------------------------- */
double channel_grid_cell_width(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double len = channel_grid_cell_length(map, col, row);
  double width = 0.0;

  if (len > 0.0) {
    while (cell != NULL) {
      width += cell->cut_width * cell->length;
      cell = cell->next;
    }
    width /= len;
  }

  return width;
}

/* -------------------------------------------------------------
   channel_grid_cell_bankheight
   ------------------------------------------------------------- */
double channel_grid_cell_bankht(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double len = channel_grid_cell_length(map, col, row);
  double height = 0.0;

  if (len > 0.0) {
    while (cell != NULL) {
      height += cell->cut_height * cell->length;
      cell = cell->next;
    }
    height /= len;
  }
  return (height);
}

/* -------------------------------------------------------------
   channel_grid_inc_inflow
   Given a flow (or actually mass), this function increases the inflow
   any channel(s) in the cell in proportion to their length within the
   cell.
   ------------------------------------------------------------- */
void channel_grid_inc_inflow(ChannelMapPtr ** map, int col, int row, float mass)
{
  ChannelMapPtr cell = map[col][row];
  float len = channel_grid_cell_length(map, col, row);

  /* 
     if (mass > 0 && len <= 0.0) {
     error_handler(ERRHDL_ERROR,
     "channel_grid_inc_inflow: attempt to add flow in cell with no channels! (col=%d, row=%d)", 
     col, row);
     return;
     }
   */

  while (cell != NULL) {
    cell->channel->lateral_inflow += mass * cell->length / len;
    cell = cell->next;
  }
}

/* -------------------------------------------------------------
   channel_grid_outflow
   If the channel(s) within the cell are marked as ``sinks'', this
   function totals the mass from the channels(s) and returns the total
   mass.
   ------------------------------------------------------------- */
double channel_grid_outflow(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  double mass = 0.0;

  while (cell != NULL) {
    if (cell->sink) {
      mass += cell->channel->outflow;
    }
    cell = cell->next;
  }
  return mass;
}

/* -------------------------------------------------------------
   channel_grid_sed_outflow
   If the channel(s) within the cell are marked as ``sinks'', this
   function totals the mass from the channels(s) and returns the total
   mass.
   ------------------------------------------------------------- */
double channel_grid_sed_outflow(ChannelMapPtr ** map, int col, int row, int i)
{
  ChannelMapPtr cell = map[col][row];
  double mass = 0.0;

  while (cell != NULL) {
    if (cell->sink) {
      mass += cell->channel->sediment.outflow[i];
    }
    cell = cell->next;
  }
  return mass;
}

/* -------------------------------------------------------------
   channel_grid_flowlength
   returns the flowlength along a road surface in a channel
   if there is more than one road in a grid cell, the road
   with the greatest surface area is used to calculate the 
   flowlength.
   This can result in a flolen that is greater than the 
   horizontal length of the road in the cell. 
  ------------------------------------------------------------- */
double channel_grid_flowlength(ChannelMapPtr ** map, int col, int row, float floslope)
{
  ChannelMapPtr cell = map[col][row];
  double flolen = 0.0;
  double area;
  double maxarea = 0.0;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    if(area > maxarea){
      flolen = cell->cut_width * (floslope/ROADCROWN)*sqrt(1+pow(ROADCROWN,2));
      maxarea = area;
    }
    if(flolen < cell->cut_width)
      flolen = cell->cut_width;
    /* If crowned, only one half goes to ditch. */ 
    if (cell->channel->class2->crown == CHAN_CROWNED) 
      flolen *= 0.5;
    
    cell = cell->next;
  }
  return flolen;
}

/* -------------------------------------------------------------
   channel_grid_flowslope
   returns the flowlength along a road surface in a grid cell
   if there is more than one road in a grid cell, the road
   with the greatest surface area is used to calculate the 
   flowslope
   ------------------------------------------------------------- */

double channel_grid_flowslope(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  float floslope = 0.0;
  double area;
  double maxarea = 0.0;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    if(area > maxarea){
      floslope = sqrt(pow(ROADCROWN, 2) + pow((cell->channel->slope),2));
      maxarea = area;
    }
    cell = cell->next;
  }
  return floslope;
}

/* -------------------------------------------------------------
   channel_grid_class
   returns the erodibility coeffienct of the road surface in a
   gird cell. if there is more than one road in a grid cell, the road
   with the greatest surface area is used
   ------------------------------------------------------------- */

ChannelClass* channel_grid_class(ChannelMapPtr ** map, int col, int row)
{
  ChannelMapPtr cell = map[col][row];
  ChannelClass *pntr;
  double area;
  double maxarea = 0.0;

  while (cell != NULL) {
    area = cell->length * cell->cut_width;
    if(area > maxarea){
      pntr = cell->channel->class2;
      maxarea = area;
    }
    cell = cell->next;
  }
  return pntr;
}

/* -------------------------------------------------------------
   ---------------------- Module Functions ---------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   channel_grid_init
   ------------------------------------------------------------- */
void channel_grid_init(int cols, int rows)
{
  channel_grid_cols = cols;
  channel_grid_rows = rows;
  channel_grid_initialized = 1;
}

/* -------------------------------------------------------------
   channel_grid_done
   ------------------------------------------------------------- */
void channel_grid_done(void)
{
  /* ? */
}

#ifdef TEST_MAIN
/* -------------------------------------------------------------
   interpolate
   ------------------------------------------------------------- */
static float interpolate(int n, float *x, float *y, float x0)
{
  int i;
  if (x0 <= x[0]) {
    return ((x0 - x[0]) / (x[1] - x[0]) * (y[1] - y[0]) + y[0]);
  }
  for (i = 0; i < n - 1; i++) {
    if (x0 < x[i + 1]) {
      return ((x0 - x[i]) / (x[i + 1] - x[i]) * (y[i + 1] - y[i]) + y[i]);
    }
  }
  return ((x0 - x[i - 1]) / (x[i] - x[i - 1]) * (y[i] - y[i - 1]) + y[i]);
}

/* -------------------------------------------------------------
   Main Program
   ------------------------------------------------------------- */
int main(int argc, char **argv)
{
  static const int columns = 5;
  static const int rows = 10;

  int r, c;

  ChannelClass *class;
  Channel *simple = NULL, *current;
  ChannelMapPtr **map = NULL;

  static int interval = 3600;	/* seconds */
  static float timestep = 1.0;	/* hour */
  static float endtime = 144.0;
#define TIMES 6
  static float bndflow[TIMES] = { 0.0, 0.0, 300.0, 300.0, 0.0, 0.0 };
  static float bndtime[TIMES] = { 0.0, 12.0, 36.0, 48.0, 60.0, 1000.0 };
  float time;

  /* module initialization */

  error_handler_init(argv[0], NULL, ERRHDL_DEBUG);
  channel_init();
  channel_grid_init(columns, rows);

  /* read channel classes */

  if ((class = channel_read_classes("example_classes.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_classes.dat: trouble reading file");
  }

  /* read a network */

  if ((simple = channel_read_network("example_network.dat", class)) == NULL) {
    error_handler(ERRHDL_FATAL, "example_network.dat: trouble reading file");
  }

  /* read channel map */

  if ((map = channel_grid_read_map(simple, "example_map.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_map.dat: trouble reading file");
  }

  /* check channel_grid_read_map */

  printf("channel_grid_read_map check:\n");
  for (r = rows - 1; r >= 0; r--) {
    for (c = 0; c < columns; c++) {
      ChannelMapPtr cell = map[c][r];
      int count;
      for (count = 0; cell != NULL; cell = cell->next) {
	count++;
      }
      printf("%3d", count);
    }
    printf("\n");
  }
  printf("\n");

  /* check channel_grid_cell_length */

  printf("channel_grid_cell_length check:\n");
  for (r = rows - 1; r >= 0; r--) {
    for (c = 0; c < columns; c++) {
      printf(" %8.2g", channel_grid_cell_length(map, c, r));
    }
    printf("\n");
  }
  printf("\n");

  /* use routing example to test channel_grid_inc_inflow and
     channel_grid_outflow */

  /* initialize flows */

  for (current = simple; current != NULL; current = current->next) {
    current->inflow = bndflow[0] * timestep;
    current->outflow = bndflow[0] * timestep;
  }

  /* time loop */

  for (time = 0.0; time <= endtime; time += timestep) {
    float inflow = interpolate(TIMES, bndtime, bndflow, time) * interval;
    float outflow;

    channel_step_initialize_network(simple);
    channel_grid_inc_inflow(map, 2, 0, inflow);
    (void) channel_route_network(simple, interval);
    outflow = channel_grid_outflow(map, 2, 6);
    channel_save_outflow(time * interval, simple, stdout);
    printf("outflow: %8.3g\n", outflow);
  }

  /* deallocate memory */

  channel_grid_free_map(map);
  channel_free_network(simple);
  channel_free_classes(class);

  /* module shutdown */

  channel_grid_done();
  channel_done();
  error_handler_done();
  exit(0);

}
#endif




