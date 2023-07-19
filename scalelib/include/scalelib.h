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

#ifdef _OPENACC
#define LSIZE 1
#else
#ifdef SINGLE
#define LSIZE (CACHELINESIZE / 8)
#else
#define LSIZE (CACHELINESIZE / 16)
#endif
#endif


#endif
