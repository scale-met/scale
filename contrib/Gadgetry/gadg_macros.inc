#ifndef IG_GADG_MACROS_INC
#define IG_GADG_MACROS_INC
/*
 ! Auther:: Akio Kawano
 ! Copyright 2006 Akio Kawano. All rights reserved.
 */
#define GADG_IS_ODD(x)    btest((x),0)
#define GADG_IS_EVEN(x)   (.not.btest((x),0))
#endif /* IG_GADG_MACROS_INC */
