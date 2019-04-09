/*
 * SUMMARY:      ParallelChannel.c
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    February 2017
 * DESCRIPTION: Routines for maintaining channel network state in a
 * parallel environment.
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2019-04-09 10:35:45 d3g096
 * COMMENTS:
 *
 *    All processes have a copy of the channel network.  All processes
 *    accumulate lateral inflow into individual segmests of the
 *    channel network.  The root process collects inflow from all
 *    processes and routes the channel netork.  The root process also
 *    writes channel network output as necessary.  Channel routing
 *    results are then distributed to all other processes.
 * 
 *    Let's not worry about the RBM fields for now, but they are
 *    probably simple accumulations too.
 */

#include <ga.h>
#include <stdlib.h>

#include "array_alloc.h"
#include "channel.h"
#include "ParallelDHSVM.h"
#include "ParallelChannel.h"

enum _channel_state_slot {
  LateralInflow = 0,
  Inflow,
  Outflow,
  Storage,
  ATP,
  ISW,
  Beam,
  Diffuse,
  NSW,
  ILW,
  VP,
  WND,
  Azimuth,
  SkyView,
  NCells
};
typedef enum _channel_state_slot ChannelStateIdx;
static const int NChannelState = NCells + 1;

/******************************************************************************/
/*                                ChannelStateGA                              */
/******************************************************************************/
int
ChannelStateGA(Channel *net)
{
  int ga, ndim, dims[GA_MAX_DIM];
  int nsegment;
  Channel *current;
  
  nsegment = 0;
  current = net;
  while (current != NULL) {
    nsegment++;
    current = current->next;
  }

  if (nsegment <= 0) return 0;

  ndim = 2;
  dims[0] = nsegment;
  dims[1] = NChannelState;
  ga = NGA_Create(C_FLOAT, ndim, &dims[0], "Channel State", NULL);
  /* GA_Print_distribution(ga); */
  return ga;
}

/******************************************************************************/
/*                           ChannelGatherCellCount                           */
/******************************************************************************/
void
ChannelGatherCellCount(Channel *net, int ga)
{
  static float one = 1.0;
  int idx, nsegment;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float *ncells, value;
  Channel *current;

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next);
  nsegment = idx;
  
  /* this an all-reduce of channel network's lateral inflow */

  ncells = (float *)calloc(nsegment, sizeof(float));
  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    ncells[idx] = current->Ncells;
  }

  lo[0] = 0;
  lo[1] = NCells;
  hi[0] = nsegment - 1;
  hi[1] = NCells;
  ld[0] = 1;
  ld[1] = 1;
  NGA_Zero_patch(ga, lo, hi);
  NGA_Acc(ga, lo, hi, &ncells[0], ld, &one);
  ParallelBarrier();

  NGA_Get(ga, lo, hi, &ncells[0], ld);

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    current->Ncells = ncells[idx];
  }

  free(ncells);
}

/******************************************************************************/
/*                    ChannelGatherLateralInflow                              */
/******************************************************************************/
void
ChannelGatherLateralInflow(Channel *net, int ga)
{
  static float one = 1.0;
  int idx, nsegment;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float *lflow, value;
  Channel *current;

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next);
  nsegment = idx;
  
  /* this an all-reduce of channel network's lateral inflow */

  lflow = (float *)calloc(nsegment, sizeof(float));
  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    lflow[idx] = current->lateral_inflow;
  }

#if 1

  lo[0] = 0;
  lo[1] = LateralInflow;
  hi[0] = nsegment - 1;
  hi[1] = LateralInflow;
  ld[0] = 1;
  ld[1] = 1;
  NGA_Zero_patch(ga, lo, hi);
  NGA_Acc(ga, lo, hi, &lflow[0], ld, &one);
  ParallelBarrier();

  /* GA_Print(ga); */

  NGA_Get(ga, lo, hi, &lflow[0], ld);

#else

  GA_Fgop(&lflow[0], nsegment, "+");

#endif
  

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    current->lateral_inflow = lflow[idx];
  }

  free(lflow);
}


/******************************************************************************/
/*                           ChannelGatherHeatBudget                          */
/******************************************************************************/
void
ChannelGatherHeatBudget(ChannelPtr net, int ga)
{
  static float one = 1.0;
  int idx, nsegment, nfield, f;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float **tmp, value;
  Channel *current;

  nfield = SkyView - ATP + 1;

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next);
  nsegment = idx;
  
  tmp = calloc_2D_float(nfield, nsegment);
  
  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    f = 0;
    tmp[f++][idx] = current->ATP;
    tmp[f++][idx] = current->ISW;
    tmp[f++][idx] = current->Beam;
    tmp[f++][idx] = current->Diffuse;
    tmp[f++][idx] = current->NSW;
    tmp[f++][idx] = current->ILW;
    tmp[f++][idx] = current->NLW;
    tmp[f++][idx] = current->VP;
    tmp[f++][idx] = current->WND;
    tmp[f++][idx] = current->azimuth;
    tmp[f++][idx] = current->skyview;
  }

  
  lo[0] = 0;
  lo[1] = ATP;
  hi[0] = nsegment - 1;
  hi[1] = SkyView;
  ld[0] = nfield;
  ld[1] = nsegment;
  NGA_Zero_patch(ga, lo, hi);
  NGA_Acc(ga, lo, hi, &tmp[0][0], ld, &one);

  ParallelBarrier();

  NGA_Get(ga, lo, hi, &tmp[0][0], ld);

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    f = 0;
    current->ATP      = tmp[f++][idx];
    current->ISW      = tmp[f++][idx];
    current->Beam     = tmp[f++][idx];
    current->Diffuse  = tmp[f++][idx];
    current->NSW      = tmp[f++][idx];
    current->ILW      = tmp[f++][idx];
    current->NLW      = tmp[f++][idx];
    current->VP       = tmp[f++][idx];
    current->WND      = tmp[f++][idx];
    current->azimuth  = tmp[f++][idx];
    current->skyview  = tmp[f++][idx];
  }

  free_2D_float(tmp);
}

/******************************************************************************/
/*                            ChannelDistributeState                          */
/* Currently Unused                                                           */
/******************************************************************************/
void
ChannelDistributeState(Channel *net, int ga)
{
  int gatype, ndim, dims[GA_MAX_DIM];
  int idx;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float value[NChannelState];
  Channel *current;

  ld[0] = 1;
  ld[1] = 1;

  NGA_Inquire(ga, &gatype, &ndim, &dims[0]);

  /* collect state from root process (which presumably did the
     routing) and put it in the channel state GA */
  if (ParallelRank() == 0){
    for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
      lo[0] = idx;
      lo[1] = LateralInflow;
      hi[0] = idx;
      hi[1] = Storage;
      value[LateralInflow] = current->lateral_inflow;
      value[Inflow] = current->inflow;
      value[Outflow] = current->outflow;
      value[Storage] = current->storage;
      NGA_Put(ga, lo, hi, &value[0], ld);
    }
  }

  ParallelBarrier();

  /* get the channel state from the GA and put it in the local copy of
     the channel network */

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    lo[0] = idx;
    lo[1] = LateralInflow;
    hi[0] = idx;
    hi[1] = Storage;
    NGA_Get(ga, lo, hi, &value[0], ld);
    current->lateral_inflow = value[LateralInflow];
    current->inflow = value[Inflow];
    current->outflow = value[Outflow];
    current->storage = value[Storage];
  }
  ParallelBarrier();

  
  
}
