/* -------------------------------------------------------------
   file: channel.c
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created October 24, 1995 by  William A Perkins
   $Id: channel.c,v 1.14 2004/10/07 20:51:08 jlanini Exp $
   ------------------------------------------------------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "errorhandler.h"
#include "channel.h"
#include "constants.h"
#include "tableio.h"
#include "settings.h"

/* for test msw */
#define TEST_MAIN 0
/* end test */

/* -------------------------------------------------------------
   -------------- ChannelClass Functions -----------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   alloc_channel_class
   ------------------------------------------------------------- */
static ChannelClass *alloc_channel_class(void)
{
  ChannelClass *p;

  if ((p = (ChannelClass *) malloc(sizeof(ChannelClass))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_class: malloc failed: %s",
		  strerror(errno));
    return NULL;
  }

  p->id = 0;
  p->width = 0.0;
  p->bank_height = 0.0;
  p->friction = 0.0;
  p->infiltration = 0.0;
  p->crown = CHAN_OUTSLOPED;
  p->erodibility_coeff = 0.0;
  p->erodibility_coeff_overland = 0.0;
  p->d50_road = 0.0;
  p->friction_road = 0.0;
  p->next = (ChannelClass *) NULL;

  return p;
}

/* -------------------------------------------------------------
   channel_free_classes
   ------------------------------------------------------------- */
void channel_free_classes(ChannelClass * head)
{
  if (head->next != NULL) {
    channel_free_classes(head->next);
  }
  free(head);
}

/* -------------------------------------------------------------
   find_channel_class
   ------------------------------------------------------------- */
static ChannelClass *find_channel_class(ChannelClass * list, ClassID id)
{
  while (list != (ChannelClass *) NULL) {
    if (list->id == id)
      break;
    list = list->next;
  }
  return list;
}

/* -------------------------------------------------------------
   channel_read_classes
   This function opens and reads the specified file and returns a
   linked list of ChannelClass structs.  If anything goes wrong, NULL
   is returned and any ChannelClass structs are destroyed.
   ------------------------------------------------------------- */
ChannelClass *channel_read_classes(const char *file, int ChanType, int Sediment)
{
  ChannelClass *head = NULL, *current = NULL;
  static const int fields = 10;
  int done;
  int err = 0;
  static char *crown_words[4] = {
    "OUTSLOPED", "CROWNED", "INSLOPED", NULL
  };

  static TableField class_fields[10] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Channel Width", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Bank (stream) or Cut Height (road)", TABLE_REAL, TRUE, FALSE, {0.0}, "",
     NULL},
    {"Friction Coefficient (Manning's n)", TABLE_REAL, TRUE, FALSE, {0.0}, "",
     NULL},
    {"Maximum Road Infiltration Rate (m/s)", TABLE_REAL, FALSE, FALSE, {0.0}, 
     "", NULL},
    {"Road Crown Type", TABLE_WORD, FALSE, FALSE, {0}, "", crown_words},
    {"Road Erodibility Coefficient", TABLE_REAL, FALSE, FALSE, {0.0}, "", NULL},
    {"Road Erodibility Overland Coefficient", TABLE_REAL, FALSE, FALSE, {0.0}, "", NULL},
    {"Road d50", TABLE_REAL, FALSE, FALSE, {0.0}, "", NULL},
    {"Road Friction Coefficient (Manning's n)", TABLE_REAL, FALSE, FALSE, {0.0}, 
     "", NULL}
  };

  // Extra fields are required if we're dealing with a road
  if (ChanType == road_class) {
    // max infiltration rate and crown type are required
    class_fields[4].required = TRUE;
    class_fields[5].required = TRUE;
    // Road erodibility coefficient and d50 are required if Sediment == TRUE,
    if (Sediment) {
      class_fields[6].required = TRUE;
      class_fields[7].required = TRUE;
      class_fields[8].required = TRUE;
      class_fields[9].required = TRUE;
    }
  }

  error_handler(ERRHDL_STATUS,
		"channel_read_classes: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
		  "channel_read_classes: unable to open file \"%s\": %s",
		  file, strerror(errno));
    return NULL;
  }

  done = FALSE;
  while (!done) {
    int i;

    done = (table_get_fields(fields, class_fields) < 0);
    if (done) {
      for (i = 0; i < fields; i++) {
	if (class_fields[i].read)
	  break;
      }
      if (i >= fields)
	continue;
    }

    if (head == NULL) {
      head = alloc_channel_class();
      current = head;
    }
    else {
      current->next = alloc_channel_class();
      current = current->next;
    }

    for (i = 0; i < fields; i++) {
      if (class_fields[i].read) {
	switch (i) {
	case 0:
	  current->id = class_fields[i].value.integer;
	  if (current->id <= 0) {
	    error_handler(ERRHDL_ERROR,
			  "%s: class %d: class id invalid", file, current->id);
	    err++;
	  }
	  break;
	  // We do not want these 3 values to be 0.0 EVER
	  // Exit if this is the case.
	case 1:
		if (class_fields[i].value.real>0.0)
			current->width = class_fields[i].value.real;
		else error_handler(ERRHDL_FATAL,
			  "channel_read_classes: %s: width cannot be 0.0", file);
	  break;
	case 2:
		if (class_fields[i].value.real>0.0)
			current->bank_height = class_fields[i].value.real;
		else error_handler(ERRHDL_FATAL,
			  "channel_read_classes: %s: bank cannot be 0.0", file);
	  break;
	case 3:
		if (class_fields[i].value.real>0.0)
			current->friction = class_fields[i].value.real;
		else error_handler(ERRHDL_FATAL,
			  "channel_read_classes: %s: friction cannot be 0.0", file);
	  break;
	case 4:
	  current->infiltration = class_fields[i].value.real;
          break;
	case 5:
	  switch (class_fields[i].value.integer) {
	  case -1:
	    error_handler(ERRHDL_ERROR,
			  "channel_read_classes: %s: unknown road crown type: %s",
			  file, class_fields[i].field);
	    err++;
	    break;
	  case 0:
	    current->crown = CHAN_OUTSLOPED;
	    break;
	  case 1:
	    current->crown = CHAN_CROWNED;
	    break;
	  case 2:
	    current->crown = CHAN_INSLOPED;
	    break;
	  default:
	    error_handler(ERRHDL_FATAL,
			  "channel_read_classes: this should not happen");
	  }
	  break;
	case 6:
	  current->erodibility_coeff = class_fields[i].value.real;
          break;
	case 7:
	  current->erodibility_coeff_overland = class_fields[i].value.real;
          break;
	case 8:
	  current->d50_road = class_fields[i].value.real;
          break;
	case 9:
	  current->friction_road = class_fields[i].value.real;
          break;
	default:
	  error_handler(ERRHDL_FATAL,
			"channel_read_classes: this should not happen either");
	}
      }
    }
  }
  error_handler(ERRHDL_STATUS,
		"channel_read_classes: %s: %d errors, %d warnings",
		file, table_errors, table_warnings);

  table_close();

  error_handler(ERRHDL_STATUS,
		"channel_read_classes: done reading file \"%s\"", file);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
		  "channel_read_classes: %s: too many errors", file);
    channel_free_classes(head);
    head = NULL;
  }

  return (head);
}

/* -------------------------------------------------------------
   ---------------------- Channel Functions --------------------
   ------------------------------------------------------------- */

/* -------------------------------------------------------------
   alloc_channel_segment
   ------------------------------------------------------------- */
static Channel *alloc_channel_segment(void)
{
  Channel *seg;
  if ((seg = (Channel *) malloc(sizeof(Channel))) == NULL) {
    error_handler(ERRHDL_ERROR, "alloc_channel_segment: malloc failed: %s",
		  strerror(errno));
    return NULL;
  }
  seg->id = 0;
  seg->order = 0;
  seg->record_name = NULL;
  seg->record = FALSE;
  seg->length = 0.0;
  seg->slope = 0.0;
  seg->class2 = NULL;
  seg->lateral_inflow = 0.0;
  seg->last_inflow = 0.0;
  seg->last_outflow = 0.0;
  seg->inflow = 0.0;
  seg->outflow = 0.0;
  seg->storage= 0.0;
  seg->last_lateral_inflow = 0.0;
  seg->outlet = NULL;
  seg->next = NULL;

  return seg;
}

/* -------------------------------------------------------------
   channel_find_segment
   A simple linear search of the channel network to find a segment
   with the given id
   ------------------------------------------------------------- */
Channel *channel_find_segment(Channel * head, SegmentID id)
{
  for (; head != NULL; head = head->next) {
    if (head->id == id)
      break;
  }
  if (head == NULL) {
    error_handler(ERRHDL_WARNING,
		  "channel_find_segment: unable to find segment %d", id);
  }
  else {
    error_handler(ERRHDL_DEBUG, "channel_find_segment: found segment %d", id);
  }

  return head;
}

/* -------------------------------------------------------------
   initialize_sediment_mass
   ------------------------------------------------------------- */
void initialize_sediment_mass(Channel * head, float **InitialSegmentSedimentm)
{
  int i;
  
  for (; head != NULL; head = head->next) {
    for(i=0;i<NSEDSIZES;i++) {
      InitialSegmentSedimentm[head->id][i] += head->sediment.mass[i]; 
    }
  }
}

/* -------------------------------------------------------------
   initialize_sediment_array
   ------------------------------------------------------------- */
void initialize_sediment_array(Channel * head, float *InitialSegmentSediment,
			       float **InitialSegmentSedimentm)
{
  int i;
  
  for (; head != NULL; head = head->next) {
    InitialSegmentSediment[head->id] += head->sediment.tempvol; 
    for(i=0;i<NSEDSIZES;i++) {
      InitialSegmentSedimentm[head->id][i] += head->sediment.tempmass[i]; 
    }
  }
}
/* -------------------------------------------------------------
   update_sediment_array
   ------------------------------------------------------------- */
void count_sediment_mass(Channel * head, float *InitialSegmentSediment)
{
  int i;
  float junk=0;

  for (; head != NULL; head = head->next){
   /*  head->sediment.tempvol = InitialSegmentSediment[head->id]; */
    for(i=0;i<NSEDSIZES;i++) {
     /*  head->sediment.tempmass[i] = head->sediment.mass[i]; */
       junk+=head->sediment.mass[i];
    }
  }
  printf("%f\n", junk);
}

/* -------------------------------------------------------------
   count_sediment_mass
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   update_sediment_array
   ------------------------------------------------------------- */
void update_sediment_array(Channel * head, float *InitialSegmentSediment, float **InitialSegmentSedimentm)
{
  int i;

  for (; head != NULL; head = head->next){
    head->sediment.tempvol = InitialSegmentSediment[head->id];
    for(i=0;i<NSEDSIZES;i++) {
      head->sediment.tempmass[i] = InitialSegmentSedimentm[head->id][i];
    }
  }
}

/* -------------------------------------------------------------
   update_sediment_array
   ------------------------------------------------------------- */
void update_sediment_mass(Channel * head, float **SegmentSedimentm, int massitertemp)
{
  int i;

  
  for (; head != NULL; head = head->next){
    for(i=0;i<NSEDSIZES;i++) { 
      head->sediment.mass[i] = SegmentSedimentm[head->id][i]
	/(float)massitertemp;
	} 
  }
}

/* -------------------------------------------------------------
   sed_vol_to_distrib_mass
   ------------------------------------------------------------- */
  /* sediment volume inflow distribute by sediment size, convert to mass */
void sed_vol_to_distrib_mass(Channel * head, float *volumearray)
{
  int i;
  float bulkporosity;
  bulkporosity = 0.245+0.14*pow(DEBRISd50,-0.21); /* Komura, 1961 relation */
  for (; head != NULL; head = head->next) {
    for(i=0;i<NSEDSIZES;i++) {
      head->sediment.debrisinflow[i] = 
	volumearray[head->id]*(1.-bulkporosity)*PARTDENSITY*(1./(float)NSEDSIZES);
    }
  }
}
/* -------------------------------------------------------------
   channel_routing_parameters
   ------------------------------------------------------------- */
void channel_routing_parameters(Channel * network, int deltat)
{
  /*   float ck; */
  float y;
  Channel *segment;

  for (segment = network; segment != NULL; segment = segment->next) {

    y = segment->class2->bank_height * 0.75;

    /* must keep 0 <= X <= 0.5 */

    /* Next line commented out by Bart Nijssen Wed Feb  3 18:57:09 1999 */
    /*     segment->X = 0.5*(1 - 0.6*y/(segment->slope*segment->length)); */

    /* must keep 0 <= X <= 0.5 */
    /* below portion commented out by westrick, shallow slopes and short segs */
    /* will yield very small y values */
    /*    if( segment->X < 0.0  ) { */
    /*  segment->X = 0.; */
    /*  y = (5.0 / 3.0 )*segment->slope*segment->length; */
    /* y which forces X = 0 */
    /*  error_handler(ERRHDL_WARNING, */
    /*                "channel_routing_parameters: segment %d: X = %f", */
    /*                segment->id, segment->X); */
    /* } */

    /* negative outflow can result unless 
       K < DT/2X so don't let K exceed that value */

    /* Next part commented out by Bart Nijssen Wed Feb  3 18:56:02 1999, since
       it is currently not being used.  Note that in the next few lines deltat
       is still in hours, NOT in seconds */
    /*     ck = 6000*sqrt(segment->slope)*pow(y,
       2.0/3.0)/segment->class->friction; */
    /*     segment->K = segment->length/ck;      */

    /*     if(segment->K > (deltat/(2.0*segment->X) ) ) { */
    /*       segment->K = deltat / (2.0*segment->X); */
    /*       error_handler(ERRHDL_WARNING, */
    /*                     "channel_routing_parameters: segment %d: K = %f", */
    /*                     segment->id, segment->K); */
    /*     } */

    /*  for new routing scheme */
    segment->K = sqrt(segment->slope) * pow((double)y, 2.0 / 3.0) /
      (segment->class2->friction * segment->length);
    segment->X = exp(-segment->K * deltat);
  }

  return;
}

/* -------------------------------------------------------------
   channel_read_network
   ------------------------------------------------------------- */
Channel *channel_read_network(const char *file, ChannelClass * class_list, int *MaxID)
{
  Channel *head = NULL, *current = NULL;
  int err = 0;
  int done;
  static const int fields = 8;
  static char *save_words[2] = {
    "SAVE", "\0"
  };
  static TableField chan_fields[8] = {
    {"ID", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Order", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Slope", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Length", TABLE_REAL, TRUE, FALSE, {0.0}, "", NULL},
    {"Class", TABLE_INTEGER, TRUE, FALSE, {0}, "", NULL},
    {"Outlet ID", TABLE_INTEGER, FALSE, FALSE, {0}, "", NULL},
    {"Save Flag", TABLE_WORD, FALSE, FALSE, {0}, "", save_words},
    {"Save Name", TABLE_STRING, FALSE, FALSE, {0.0}, "", NULL}
  };

  error_handler(ERRHDL_STATUS,
		"channel_read_network: reading file \"%s\"", file);

  if (table_open(file) != 0) {
    error_handler(ERRHDL_ERROR,
		  "channel_read_network: unable to open file \"%s\": %s",
		  file, strerror(errno));
    return NULL;
  }

  *MaxID = 0;
  done = FALSE;
  while (!done) {
    int i;

    done = (table_get_fields(fields, chan_fields) < 0);

    if (done) {
      for (i = 0; i < fields; i++) {
	if (chan_fields[i].read)
	  break;
      }
      if (i >= fields)
	continue;
    }

    if (head == NULL) {
      head = alloc_channel_segment();
      current = head;
    }
    else {
      current->next = alloc_channel_segment();
      current = current->next;
    }

    for (i = 0; i < fields; i++) {
      if (chan_fields[i].read) {
	switch (i) {
	case 0:
	  current->id = chan_fields[i].value.integer;
	  if(current->id > *MaxID) *MaxID = current->id;
	  if (current->id <= 0) {
	    error_handler(ERRHDL_ERROR,
			  "%s: segment %d: channel id invalid",
			  file, current->id);
	    err++;
	  }
	  break;
	case 1:
	  if (chan_fields[i].value.integer > 0) {
	    current->order = chan_fields[i].value.integer;
	  }
	  else {
	    error_handler(ERRHDL_ERROR,
			  "%s: segment %d: channel order (%d) invalid",
			  file, current->id, chan_fields[i].value.integer);
	    err++;
	  }
	  break;
	case 2:
	  if (chan_fields[i].value.real > 0) {
	    current->slope = chan_fields[i].value.real;
	  }
	  else {
	    error_handler(ERRHDL_ERROR,
			  "%s: segment %d: channel slope (%f) invalid",
			  file, current->id, chan_fields[i].value.real);
	    err++;
	  }
	  break;
	case 3:
	  if (chan_fields[i].value.real > 0) {
	    current->length = chan_fields[i].value.real;
	  }
	  else {
	    error_handler(ERRHDL_ERROR,
			  "%s: segment %d: channel length (%f) invalid",
			  file, current->id, chan_fields[i].value.real);
	    err++;
	  }
	  break;
	case 4:
	  if ((current->class2 =
	       find_channel_class(class_list, chan_fields[i].value.integer)
	      ) == NULL) {
	    error_handler(ERRHDL_ERROR,
			  "%s: segment %d: channel class %d not found",
			  file, current->id, chan_fields[i].value.integer);
	    err++;
	  }
	  break;
	case 5:
	  current->outlet = (Channel *) chan_fields[i].value.integer;
	  break;
	case 6:
	  current->record = TRUE;
	  break;
	case 7:
	  current->record_name = (char *) strdup(chan_fields[i].field);
	  break;
	default:
	  error_handler(ERRHDL_FATAL,
			"channel_read_network: what is this field %d?", i);
	  break;
	}
      }
    }
  }

  table_close();

  /* find segment outlet segments, if
     specified */

  for (current = head; current != NULL; current = current->next) {
    int outid = (int) current->outlet;

    if (outid != 0) {
      current->outlet = channel_find_segment(head, outid);
      if (current->outlet == NULL) {
	error_handler(ERRHDL_ERROR,
		      "%s: cannot find outlet (%d) for segment %d",
		      file, outid, current->id);
	err++;
      }
    }
  }

  table_errors += err;

  error_handler(ERRHDL_STATUS,
		"channel_read_network: %s: %d errors, %d warnings",
		file, table_errors, table_warnings);

  if (table_errors) {
    error_handler(ERRHDL_ERROR,
		  "channel_read_network: %s: too many errors", file);
    channel_free_network(head);
    head = NULL;
  }

  return (head);
}

/* -------------------------------------------------------------
   channel_route_segment
   ------------------------------------------------------------- */
static int channel_route_segment(Channel * segment, int deltat)
{
  float K = segment->K;
  float X = segment->X;
  /*   float C1, C2, C3, C4; */
  float inflow, last_inflow;
  float outflow, last_outflow, lateral_inflow, last_lateral_inflow;
  float storage;
  int err = 0;

  /* change masses to rates */

  last_inflow = segment->last_inflow / deltat;
  inflow = segment->inflow / deltat;
  last_outflow = segment->last_outflow / deltat;
  lateral_inflow = segment->lateral_inflow / deltat;
  last_lateral_inflow = segment->last_lateral_inflow / deltat;

  /* The following lines are currently not used and have therefore been
     commented out, Bart Nijssen, Wed Feb  3 19:06:33 1999.  Note that in this
     sections the units for C1, C2, C3, and C4 need to be checked before using
     the code.  The units of deltat have been changed from hours to seconds, but 
     that might not yet be reflected correctly in the next couple of lines */
  /*   C1 = deltat - 2*K*X; */
  /*   C2 = deltat + 2*K*X; */
  /*   C3 = 2*K*(1 - X) - deltat; */
  /*   C4 = lateral_inflow*deltat; */

  /*   outflow = ( C1*inflow + C2*last_inflow + C3*last_outflow )/  */
  /*     (2*K*(1 - X) + deltat); */

  /*   storage = K*(X*inflow + (1 - X)*outflow); */

  /* new code */
  storage = ((inflow + lateral_inflow) / K) +
    (segment->storage - (inflow + lateral_inflow) / K) * X;
  if (storage < 0.0)
    storage = 0.0;
  outflow = (inflow + lateral_inflow) - (storage - segment->storage) / deltat;

/*  FOR TEST  
outflow = inflow + lateral_inflow;
storage = 0.0;   
*/

  segment->outflow = outflow * deltat;
  segment->storage = storage;

  if (segment->outlet != NULL)
    segment->outlet->inflow += segment->outflow;

  return (err);
}

/* -------------------------------------------------------------
   channel_route_network
   ------------------------------------------------------------- */
int channel_route_network(Channel * net, int deltat)
{
  int order;
  int order_count;
  int err = 0;
  Channel *current;

  for (order = 1;; order += 1) {
    order_count = 0;
    current = net;
    while (current != NULL) {
      if (current->order == order) {
	err += channel_route_segment(current, deltat);
	order_count += 1;
      }
      current = current->next;
    }
    if (order_count == 0)
      break;
  }
  return (err);
}

/* -------------------------------------------------------------
   channel_step_initialize_network
   ------------------------------------------------------------- */
int channel_step_initialize_network(Channel * net)
{
  if (net != NULL) {
    net->last_inflow = net->inflow;
    net->inflow = 0.0;
    net->lateral_inflow = 0.0;
    net->last_outflow = net->outflow;
    net->last_storage = net->storage;
    net->last_lateral_inflow = net->lateral_inflow;
    channel_step_initialize_network(net->next);
  }
  return (0);
}

/* -------------------------------------------------------------
   channel_step_initialize_sednetwork
   currently this isn't used because variables need initialization
   at different times
   ------------------------------------------------------------- */
int channel_step_initialize_sednetwork(Channel * net)
{
  int i;

  if (net != NULL) {
    for(i = 0; i < NSEDSIZES;i++) {
      net->sediment.debrisinflow[i] = 0.;
      net->sediment.overlandinflow[i] = 0.;
      net->sediment.overroadinflow[i] = 0.;
      net->sediment.inflow[i] = 0.;
    }
    channel_step_initialize_sednetwork(net->next);
  }
  return (0);
}

/* -------------------------------------------------------------
   channel_save_outflow
   This routine saves the channel output
   ------------------------------------------------------------- */
int channel_save_outflow(double time, Channel * net, FILE * out, FILE * out2)
{
  char buffer[16];
  sprintf(buffer, "%12.5g", time);
  return (channel_save_outflow_text(buffer, net, out, out2, 0));
}

/* -------------------------------------------------------------
   channel_save_outflow_wtext
   Saves the channel outflow using a text string as the time field
   ------------------------------------------------------------- */
int
channel_save_outflow_text(char *tstring, Channel * net, FILE * out,
			  FILE * out2, int flag)
{
  int err = 0;
  float total_outflow = 0.0;
  float total_lateral_inflow = 0.0;
  float total_storage = 0.0;
  float total_storage_change = 0.0;
  float total_error = 0.0;

  if (flag == 1) {
    fprintf(out2, "DATE ");
    for (; net != NULL; net = net->next) {
      total_lateral_inflow += net->lateral_inflow;
      if (net->outlet == NULL) {
	total_outflow += net->outflow;
      }
      if (net->record)
	fprintf(out2, "%s ", net->record_name);
    }
    fprintf(out2, "\n");
  }

  //tsstring = date in the form of 01.01.1915-00:00:00 
  if (fprintf(out2, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,
		  "channel_save_outflow: write error:%s", strerror(errno));
    err++;
  }

  for (; net != NULL; net = net->next) {
    total_lateral_inflow += net->lateral_inflow;
    if (net->outlet == NULL) {
      total_outflow += net->outflow;
    }
    total_storage += net->storage;
    total_storage_change += net->storage - net->last_storage;

    if (net->record) {
      if (fprintf(out, "%15s %10d %12.5g %12.5g %12.5g %12.5g",
		  tstring, net->id, net->inflow, net->lateral_inflow,
		  net->outflow, net->storage - net->last_storage) == EOF) {
	error_handler(ERRHDL_ERROR,
		      "channel_save_outflow: write error:%s", strerror(errno));
	err++;
      }
      if (fprintf(out2, "%12.5g ", net->outflow) == EOF) {
	error_handler(ERRHDL_ERROR,
		      "channel_save_outflow: write error:%s", strerror(errno));
	err++;
      }
      if (net->record_name != NULL) {
	if (fprintf(out, "   \"%s\"\n", net->record_name) == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_outflow: write error:%s",
			strerror(errno));
	  err++;
	}

      }
      else {
	if (fprintf(out, "\n") == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_outflow: write error:%s",
			strerror(errno));
	  err++;
	}
      }
    }
  }
  total_error = total_storage_change - total_lateral_inflow + total_outflow;
  if (fprintf(out, "%15s %10d %12.5g %12.5g %12.5g %12.5g %12.5g \"Totals\"\n",
	      tstring, 0, total_lateral_inflow,
	      total_outflow, total_storage,
	      total_storage_change, total_error) == EOF) {
    error_handler(ERRHDL_ERROR,
		  "channel_save_outflow: write error:%s", strerror(errno));
    err++;
  }
  fprintf(out2, "\n");

  return (err);
}

/* -------------------------------------------------------------
   channel_free_network
   ------------------------------------------------------------- */
void channel_free_network(Channel * net)
{
  if (net->next != NULL) {
    channel_free_network(net->next);
  }
  free(net);
}

/* -------------------------------------------------------------
   channel_init
   ------------------------------------------------------------- */
void channel_init(void)
{
  /* do nothing */
  return;
}

/* -------------------------------------------------------------
   channel_done
   ------------------------------------------------------------- */
void channel_done(void)
{
  /* do nothing */
  return;
}

#if TEST_MAIN

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
  static int interval = 3600;	/* timestep in seconds */
  static int timestep;
  static int endtime = 144;
#define TIMES 6
  static float bndflow[TIMES] = { 0.0, 0.0, 300.0, 300.0, 0.0, 0.0 };
  static float bndtime[TIMES] = { 0.0, 12.0, 36.0, 48.0, 60.0, 1000.0 };

  float time;
  ChannelClass *class;
  Channel *simple = NULL, *current, *tail;

  error_handler_init(argv[0], NULL, ERRHDL_ERROR);
  channel_init();

  /* read classes */

  if ((class = channel_read_classes("example_classes.dat")) == NULL) {
    error_handler(ERRHDL_FATAL, "example_classes.dat: trouble reading file");
  }

  /* read a network */

  if ((simple = channel_read_network("example_network.dat", class)) == NULL) {
    error_handler(ERRHDL_FATAL, "example_network.dat: trouble reading file");
  }

  /* initialize flows */

  for (current = simple; current != NULL; current = current->next) {
    current->inflow = bndflow[0];
    current->outflow = bndflow[0];
    current->outlet = current->next;
    tail = current;
  }

  /* time loop */

  for (timestep = 0; timestep <= endtime; timestep++) {
    float inflow = interpolate(TIMES, bndtime, bndflow, timestep) * interval;
    float outflow;

    channel_step_initialize_network(simple);
    simple->inflow = inflow;
    (void) channel_route_network(simple, interval);
    outflow = tail->outflow / interval;
    channel_save_outflow(timestep * interval, simple, stdout);
  }

  channel_free_network(simple);
  channel_free_classes(class);
  channel_done();
  error_handler_done();
  exit(0);
}
#endif

/* -------------------------------------------------------------
   channel_save_sed_outflow_text
   Saves the channel sediment outflow using a text string as the time field
   ------------------------------------------------------------- */
int
channel_save_sed_outflow_text(char *tstring, Channel * net, FILE * out,
			  FILE * out2, int flag)
{
  int err = 0;

  /* print header line first time through */
  if (flag == 1) {
    fprintf(out2, "DATE ");
    for (; net != NULL; net = net->next) {
      if (net->record)
	fprintf(out2, "%s ", net->record_name);
    }
    fprintf(out2, "\n");
  }

  if (fprintf(out2, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,
		  "channel_save_sed_outflow: write error:%s", strerror(errno));
    err++;
  }
  
  for (; net != NULL; net = net->next) {
    if (net->record) {
      
      if (fprintf(out, "%15s %10d %12.5g %12.5g",
		  tstring, net->id, net->sediment.totalmass, net->sediment.outflowconc) == EOF) {
      	error_handler(ERRHDL_ERROR,
		      "channel_save_sed_outflow: write error:%s", strerror(errno));
	err++;
      }
      if (fprintf(out2, "%12.5g ", net->sediment.outflowconc) == EOF) {
	error_handler(ERRHDL_ERROR,
		      "channel_save_sed_outflow: write error:%s", strerror(errno));
	err++;
      }
      if (net->record_name != NULL) {
	if (fprintf(out, "   \"%s\"\n", net->record_name) == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_sed_outflow: write error:%s",
			strerror(errno));
	  err++;
	}
	
      }
      else {
	if (fprintf(out, "\n") == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_outflow: write error:%s",
			strerror(errno));
	  err++;
	}
      }
    }
  }
  fprintf(out2, "\n");

  return (err);
}

/* -------------------------------------------------------------
   channel_save_sed_inflow_text
   Saves the channel sediment lateral inflows using a text string as the time field
   ------------------------------------------------------------- */
int
channel_save_sed_inflow_text(char *tstring, Channel * net, FILE * out,
			     float *SedDiams, int flag)
{
  int err = 0;
  int i, j;
  int count = 0;
  
  /* print header line first time through */
  if (flag == 1) {
    fprintf(out, "DATE ");
    for (; net != NULL; net = net->next) {
      if (net->record){
	for(i = 0; i < NSEDSIZES;i++) { 
	  fprintf(out, "%s ", net->record_name);
	}
	count++;
      }
    }
    fprintf(out, "\n");
    
    fprintf(out, "SEDDIAMS ");
    for (j = 0 ; j < count; j++){
  	for(i = 0; i < NSEDSIZES;i++) { 
	  fprintf(out, "%f ", SedDiams[i]);
	}
    }
    fprintf(out, "\n");
  } 
  
  if (fprintf(out, "%15s ", tstring) == EOF) {
    error_handler(ERRHDL_ERROR,
		  "channel_save_sed_inflow: write error:%s", strerror(errno));
    err++;
  }
  
  for (; net != NULL; net = net->next) {
    if (net->record) {
      for(i = 0; i < NSEDSIZES;i++) { 
	if (fprintf(out, "%12.5g ", net->sediment.debrisinflow[i] + net->sediment.overlandinflow[i] +
		    net->sediment.overroadinflow[i]) == EOF) {
	  error_handler(ERRHDL_ERROR,
			"channel_save_sed_inflow: write error:%s", strerror(errno));
	  err++;
	}
      }
    }
  }
  fprintf(out, "\n");
  
  return (err);
}

