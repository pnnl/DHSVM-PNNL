/* -------------------------------------------------------------
   file: channel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created October 24, 1995 by  William A Perkins
   $Id: channel.h,v3.1.2 2014/1/2 ning Exp $
   ------------------------------------------------------------- */

#ifndef _channel_h_
#define _channel_h_

typedef unsigned short int SegmentID, ClassID;

/* -------------------------------------------------------------
   struct ChannelClass
   ------------------------------------------------------------- */
typedef enum {
  CHAN_OUTSLOPED, CHAN_CROWNED, CHAN_INSLOPED
} ChannelCrownType;

typedef struct _channel_class_rec_ {
  ClassID id;			/* unique identifier */

  float width;			/* ``channel'' width */
  float bank_height;		/* bank height for streams (or cut height for roads) */
  float friction;		        /* Manning's n for the channel*/
  /* The following variables are only used when the channel class is a road network. */
  float infiltration;		/* infiltration through ditch surface - roads only
				   Note: this may not be what you think it is, so be
				   sure to read the documentation before you use
				   it.  It is ONLY used for road networks and if
				   the option ROAD INFILTRATION is set to TRUE. */
  ChannelCrownType crown;	/* crown type - roads only */
  struct _channel_class_rec_ *next;

} ChannelClass;

/* -------------------------------------------------------------
   structure for storing the channel riparian vegetation information
   (used only if stream temperature is set true)
   ------------------------------------------------------------- */
typedef struct {
	float TREEHEIGHT;
	float BUFFERWIDTH;          /* The width of canopy buffer */
	float OvhCoeff;             /* A percentage of tree height thats used to represent overhanging canopy */
	float ExtnCoeff[12];        /* Monthly Extinction coefficient */
	float Extn;
	float CanopyBankDist;       /* Distance from bank to canopy */
    float StreamWidth;          /* segment width used in riparian shading module */
} CHANRVEG;

/* -------------------------------------------------------------
   struct Channel
   This is the basic unit of channel information.
   ------------------------------------------------------------- */
struct _channel_rec_ {
  SegmentID id;
  unsigned order;		/* determines computation order */
  char *record_name;	/* The name this segment is to have in the output, if output is recorded */
  char record;			/* TRUE if outflow values are to be saved by channel_save_outflow */
  float length;			/* Parameters */
  float slope;
  float K;              /* Travel time constant, a function of slope */
  float X;              /* Weighting factor (0~1), exponential function of K */
  ChannelClass *class2;	/* ChannelClass identifier */

  /* necessary routing terms */
  float lateral_inflow;	/* cubic meters */
  float last_inflow;	/* cubic meters */
  float last_outflow;	/* cubic meters */
  float last_storage;	/* cubic meters */
  float inflow;			/* cubic meters */
  float outflow;		/* cubic meters */
  float storage;		/* cubic meters */
  float last_lateral_inflow;
  /* Added for John's RBM model */
  float ATP;	        /* Avg air temp (C) */
  float ISW;            /* Incident incoming shortwave radiation (W/m2) */
  float Beam;           /* Incident incoming beam shortwave radiation (W/m2) */
  float Diffuse;        /* Incident incoming diffuse shortwave radiation (W/m2) */
  float NSW;            /* Net shortwave (W/m2) */
  float ILW;            /* Incident incoming longwave radiation (W/m2) */
  float NLW;            /* Net incoming longwave radiation (W/m2) */
  float VP;             /* Actual vapor pressure (pa) */
  float WND;	        /* Wind (m/s) */
  float azimuth;        /* segment azimuth (degrees) */
  float skyview;
  int Ncells;	        /* Number of grid cells crossed by the segment*/

  CHANRVEG rveg;        /* riparian veg sub-structure */

  struct _channel_rec_ *outlet;	/* NULL if does not drain to another segment */
  struct _channel_rec_ *next;
};
typedef struct _channel_rec_ Channel, *ChannelPtr;

/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

				/* ChannelClass */

ChannelClass *channel_read_classes(const char *file, int ChanType);
void channel_free_classes(ChannelClass *head);

				/* Channel */

Channel *channel_read_network(const char *file, ChannelClass * class_list, int *MaxID);
int channel_read_rveg_param(Channel *net, const char *file, int *MaxID);
void channel_routing_parameters(Channel *net, int deltat);
Channel *channel_find_segment(Channel *net, SegmentID id);
int channel_step_initialize_network(Channel *net);
int channel_incr_lat_inflow(Channel *segment, float linflow);
int channel_route_network(Channel *net, int deltat);
int channel_save_outflow(double time, Channel * net, FILE *file, FILE *file2);
int channel_save_outflow_text(char *tstring, Channel *net, FILE *out,
			      FILE *out2, int flag);
void channel_free_network(Channel *net);

				/* Module */

void channel_init(void);
void channel_done(void);

#endif
