#ifndef IG_GADG_CONFIG_INC
#define IG_GADG_CONFIG_INC
 
# include "gadg_config_platform.inc"

! Sets default values.

# ifndef PP_CONF_VEC_REG_LEN
#   define PP_CONF_VEC_REG_LEN     256
# endif

# ifndef PP_CONF_USE_FOR_ABORT
#   define PP_CONF_USE_FOR_ABORT
# endif

# ifndef PP_CONF_ABORT
#   define PP_CONF_ABORT           stop
# endif

#endif /* IG_GADG_CONFIG_INC */
