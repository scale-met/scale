#ifndef SCALE_H
#define SCALE_H

#include "scale_log.h"
#include "scale_openmp.h"

#ifndef VECTLEN
#define VECTLEN 8
#endif

#ifndef CACHELINESIZE
#define CACHELINESIZE 64
#endif


#define LSIZE CACHELINESIZE/RP


#endif
