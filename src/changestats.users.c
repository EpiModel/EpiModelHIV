/*  File src/changestats.users.c in package ergm.userterms, part of the Statnet suite
 *  of packages for network analysis, http://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  http://statnet.org/attribution
 *
 *  Copyright 2003-2013 Statnet Commons
 */
#include "changestats.users.h"

/* absdiffnodemix */

CHANGESTAT_FN(d_absdiffnodemix) {
  double change; Vertex t, h; int i, nnodes, nstats, statnum; 
  ZERO_ALL_CHANGESTATS(i);
  FOR_EACH_TOGGLE(i) {
    t = TAIL(i); h = HEAD(i);
    nnodes = INPUT_PARAM[0];
    nstats = INPUT_PARAM[1];
    change = fabs(INPUT_PARAM[t+1] - INPUT_PARAM[h+1]);
    for (statnum = 0; statnum < INPUT_PARAM[1]; statnum++) 
      {    
        if ((INPUT_PARAM[nnodes+t+1] == INPUT_PARAM[2*nnodes+statnum+2] && 
            INPUT_PARAM[nnodes+h+1] == INPUT_PARAM[2*nnodes+nstats+statnum+2]) ||
            (INPUT_PARAM[nnodes+t+1] == INPUT_PARAM[2*nnodes+nstats+statnum+2] && 
            INPUT_PARAM[nnodes+h+1] == INPUT_PARAM[2*nnodes+statnum+2]))
          {CHANGE_STAT[statnum] += IS_OUTEDGE(t,h) ? -change : change;}
      }
    TOGGLE_IF_MORE_TO_COME(i);
  }
  UNDO_PREVIOUS_TOGGLES(i);
}
