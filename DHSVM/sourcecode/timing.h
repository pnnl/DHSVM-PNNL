/*
 * SUMMARY:      timing.h
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       William A. Perkins
 * ORG:          Pacific NW National Laboratory
 * E-MAIL:       william.perkins@pnnl.gov
 * ORIG-DATE:    May 2018
 * DESCRIPTION:  
 *
 * DESCRIP-END.cd
 * FUNCTIONS:    
 * LAST CHANGE: 2018-05-09 08:12:44 d3g096
 * COMMENTS:
 */

#ifndef _timing_h_
#define _timing_h_

#ifndef TIMING_MAX_LEVEL
#define TIMING_MAX_LEVEL 1
#endif

#if defined(GPTL_TIMING)

#include "gptl.h"

#define TIMING_INIT()                           \
  GPTLsetoption(GPTLabort_on_error, 1);         \
  GPTLsetoption(GPTLsync_mpi, 1);               \
  GPTLsetoption(GPTLpercent, 1);                \
  GPTLsetoption(GPTLoverhead, 1);               \
  GPTLinitialize();
  
#define TIMING_TASK_START(name, level) if (level <= TIMING_MAX_LEVEL) GPTLstart(name);
#define TIMING_TASK_END(name, level) if (level <= TIMING_MAX_LEVEL) GPTLstop(name);
#define TIMING_DONE(me) GPTLpr(me)

#else

#define TIMING_INIT() /* */
#define TIMING_TASK_START(name, level) /* */
#define TIMING_TASK_END(name, level) /* */
#define TIMING_DONE(me) /* */

#endif

#endif

