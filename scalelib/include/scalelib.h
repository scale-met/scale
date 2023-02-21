#ifndef SCALE_H
#define SCALE_H

#include "scale_log.h"
#include "scale_openmp.h"

#ifndef VECTLEN
#ifdef _OPENACC
#define VECTLEN 32
#else
#define VECTLEN 8
#endif
#endif

#ifndef CACHELINESIZE
#define CACHELINESIZE 128
#endif


#endif
