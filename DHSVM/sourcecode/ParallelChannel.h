/*
 * SUMMARY:      ParallelChannel.h
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    February 2017
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2017-02-15 09:47:48 d3g096
 * COMMENTS:
 */

#ifndef _ParallelChannel_h_
#define _ParallelChannel_h_

int ChannelStateGA(Channel *net);
void ChannelGatherLateralInflow(Channel *net, int ga);
void ChannelDistributeState(Channel *net, int ga);


#endif

