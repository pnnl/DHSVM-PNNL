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
 * LAST CHANGE: 2017-02-15 10:42:35 d3g096
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

#include "channel.h"
#include "ParallelDHSVM.h"
#include "ParallelChannel.h"

enum _channel_state_slot {
  LateralInflow = 0,
  Inflow,
  Outflow,
  Storage
};
typedef enum _channel_state_slot ChannelStateIdx;
static const int NChannelState = Storage + 1;

/******************************************************************************/
/*                                ChannelStateGA                              */
/******************************************************************************/
int
ChannelStateGA(Channel *net)
{
  int ga, pg, ndim, dims[GA_MAX_DIM];
  int nsegment;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM];
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
}

/******************************************************************************/
/*                    ChannelGatherLateralInflow                              */
/******************************************************************************/
void
ChannelGatherLateralInflow(Channel *net, int ga)
{
  static const one = 1.0;
  int gatype, ndim, dims[GA_MAX_DIM];
  int idx, nsegment;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float value;
  Channel *current;

  ld[0] = 1;
  ld[1] = 1;

  NGA_Inquire(ga, &gatype, &ndim, &dims[0]);
  nsegment = dims[0];

  /* put zero in the lateral_inflow slot in all segments */

  value = 0.0;
  lo[0] = 0;
  lo[1] = LateralInflow;
  hi[0] = nsegment - 1;
  hi[1] = LateralInflow;
  NGA_Fill_patch(ga, lo, hi, &value);

  /* go thru the network and accumulate lateral inflow into the GA */

  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    lo[0] = idx;
    lo[1] = LateralInflow;
    hi[0] = idx;
    hi[1] = LateralInflow;
    value = current->lateral_inflow;
    NGA_Acc(ga, lo, hi, &value, ld, &one);
  }
  GA_Sync();

  /* get the result back in the local network */
  for (idx = 0, current = net; current != NULL; ++idx, current = current->next) {
    lo[0] = idx;
    lo[1] = LateralInflow;
    hi[0] = idx;
    hi[1] = LateralInflow;
    NGA_Get(ga, lo, hi, &value, ld);
    current->lateral_inflow = value;
  }
  GA_Sync();
}

/******************************************************************************/
/*                            ChannelDistributeState                          */
/******************************************************************************/
void
ChannelDistributeState(Channel *net, int ga)
{
  int gatype, ndim, dims[GA_MAX_DIM];
  int idx, nsegment;
  int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM];
  float value[NChannelState];
  Channel *current;

  ld[0] = 1;
  ld[1] = 1;

  NGA_Inquire(ga, &gatype, &ndim, &dims[0]);
  nsegment = dims[0];

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

  GA_Sync();

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
  GA_Sync();


  
}
