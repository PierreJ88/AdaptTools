#ifndef __MS_CALLS_H
#define __MS_CALLS_H

/* data structure */
typedef struct _MSst MSst;

/* prototypes */
MSst *MS_init(int dim,int ver);
int   MS_stop(MSst *msst);
int   MS_mshmet(MSst *msst);


#endif