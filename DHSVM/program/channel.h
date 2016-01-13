/* -------------------------------------------------------------
   file: channel.h
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Battelle Memorial Institute
   Pacific Northwest Laboratory
   ------------------------------------------------------------- */
/* -------------------------------------------------------------
   Created October 24, 1995 by  William A Perkins
   $Id: channel.h,v 1.12 2004/10/07 20:51:08 jlanini Exp $
   ------------------------------------------------------------- */

#ifndef _channel_h_
#define _channel_h_

#define NSEDSIZES 3 /* number of sediment sizes used for transport */

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
  float erodibility_coeff;	/* erodibility coefficient - roads only
				   Note: this is only used when SEDIMENT is set to TRUE. */
  float erodibility_coeff_overland;	/* erodibility coefficient overland - roads only
				   Note: this is only used when SEDIMENT is set to TRUE. */
  float d50_road;                /* mean diameter (mm) for the road segment - roads only
				   Note: this is only used when SEDIMENT is set to TRUE. */
  float friction_road;           /* Manning's n for the road surface - roads only
				    Note: this is only used when SEDIMENT is set to TRUE. */
  struct _channel_class_rec_ *next;

} ChannelClass;

/* -------------------------------------------------------------
   structure for storing the channel sediment information 
   ------------------------------------------------------------- */
typedef struct {
  float mass[NSEDSIZES];             /* stored mass, kg, of sediment in channel by D */
  float debrisinflow[NSEDSIZES];     /* inflow of sediment from mass wasting, kg, by D */
  float overlandinflow[NSEDSIZES];   /* inflow of sediment from overland
					flow erosion, kg, by D */
  float overroadinflow[NSEDSIZES];   /* inflow of sediment from over road
					flow erosion, kg, by D */
  float inflow[NSEDSIZES];           /* inflow from upstream reach, kg, by D */
  float inflowrate[NSEDSIZES];           /* inflow from upstream reach, kg/s, by D */
  float last_inflowrate[NSEDSIZES]; 
  float outflow[NSEDSIZES];          /* outflow to downstream reach, kg, by D */
  float last_outflow[NSEDSIZES];          /* outflow to downstream reach, kg, by D */
  float outflowrate[NSEDSIZES];          /* outflow to downstream reach, kg/s, by D */
  float last_outflowrate[NSEDSIZES]; 
  float tempvol;           /*volume of debris inflow - total. Temporary space */
  float tempmass[NSEDSIZES];   /* channel storage, in kg, that can move due to 
				  debris flows. Temporary space */
  float totalmass;         /* total sediment mass in reach */
  float outflowconc;       /* outflow concentration in ppm */
} CHANSED;

/* -------------------------------------------------------------
   struct Channel
   This is the basic unit of channel information.
   ------------------------------------------------------------- */
struct _channel_rec_ {
  SegmentID id;

  unsigned order;		/* determines computation order */
  char *record_name;		/* The name this segment is to have in
				   the output, if output is recorded */
  char record;			/* TRUE if outflow values are to be
				   saved by channel_save_outflow */

  float length;			/* Parameters */
  float slope;
  float K;
  float X;

  ChannelClass *class2;		/* ChannelClass identifier */

  /* necessary routing terms */

  float lateral_inflow;		/* cubic meters */
  float last_inflow;		/* cubic meters */
  float last_outflow;		/* cubic meters */
  float last_storage;		/* cubic meters */
  float inflow;			/* cubic meters */
  float outflow;		/* cubic meters */
  float storage;		/* cubic meters */
  float last_lateral_inflow;	/* cubic meters */

  CHANSED sediment;            /* sediment sub-structure */

  struct _channel_rec_ *outlet;	/* NULL if does not drain to another segment */

  struct _channel_rec_ *next;
};
typedef struct _channel_rec_ Channel, *ChannelPtr;

/* -------------------------------------------------------------
   externally available routines
   ------------------------------------------------------------- */

				/* ChannelClass */

ChannelClass *channel_read_classes(const char *file, int ChanType, int Sediment);
void channel_free_classes(ChannelClass * head);

				/* Channel */

Channel *channel_read_network(const char *file, ChannelClass * class_list, int *MaxID);
void channel_routing_parameters(Channel * net, int deltat);
Channel *channel_find_segment(Channel * net, SegmentID id);
int channel_step_initialize_network(Channel * net);
int channel_step_initialize_sednetwork(Channel * net);
int channel_incr_lat_inflow(Channel * segment, float linflow);
int channel_route_network(Channel * net, int deltat);
int channel_save_outflow(double time, Channel * net, FILE * file, FILE * file2);
int channel_save_outflow_text(char *tstring, Channel * net, FILE * out,
			      FILE * out2, int flag);
int channel_save_sed_outflow_text(char *tstring, Channel * net, FILE * out,
			      FILE * out2, int flag);
int channel_save_sed_inflow_text(char *tstring, Channel * net, FILE * out,
				 float *SedDiams, int flag);
void channel_free_network(Channel * net);

				/* Module */

void channel_init(void);
void channel_done(void);
void initialize_sediment_array(Channel * head, float *InitialSegmentSediment, 
			       float **InitialSegmentSedimentm);
void initialize_sediment_mass(Channel * head, float **InitialSegmentSedimentm);
void count_sediment_mass(Channel * head, float *InitialSegmentSediment);
void update_sediment_array(Channel * head, float *InitialSegmentSediment, float **InitialSegmentSedimentm);
void update_sediment_mass(Channel * head, float **SegmentSedimentm, int massitertemp);
void sed_vol_to_distrib_mass(Channel * head, float *volumearray);

#endif
