#ifndef PROFILER_H_
#define PROFILER_H_

#if defined(USE_SCALASCA)
#include "epik_user.inc"
#define PROFILE_REGION_DECLARE(rid)  EPIK_USER_REG(rid, "rid")
#define PROFILE_REGION_BEGIN(rid) EPIK_USER_START(rid)
#define PROFILE_REGION_END(rid) EPIK_USER_END(rid)
#else
#define PROFILE_REGION_DECLARE(id)  
#define PROFILE_REGION_BEGIN(id) 
#define PROFILE_REGION_END(id) 
#endif

#endif 
