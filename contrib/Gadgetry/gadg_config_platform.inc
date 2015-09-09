/*
! Configuration file for C-preprocessor.
! You can modify the macro-definitions to fit your platform
!
!*** Switches
!
! PP_CONF_VECTORIZE     : Outputs vectorized code if defined.
! PP_CONF_USE_MPI       : Define if MPI is used.
! PP_CONF_64BIT_CPU     : Define if 64bit-CPU's are used.
!
!*** CONSTANTS
!
! PP_CONF_VEC_REG_LEN   : Vector register length. default is '256'.
! PP_CONF_ABORT         : Directive for abnormal termination.
!                         Default is 'STOP'.
!                         This option is ignored if CONF_USE_MPI is defined.
! PP_CONF_USE_FOR_ABORT : Define this as 'use A_MODULE'.
!                         If CONF_ABORT needs to use A_MODULE.

!*** Switches for compiler dependent optimal codes
! PP_COMPILER_DEPENDENT_A      in "gadg_algorithm_accum.inc"

*/
!#define PP_CONF_USE_MPI

# if defined(ES) || defined(SX)
    ! Fortran/ES, Fortran/SX
#   define PP_CONF_VECTORIZE
#   define PP_CONF_VEC_REG_LEN     256
#   define PP_CONF_ABORT           call abterm(-1)
#   define PP_CONF_64BIT_CPU
#   define PP_COMPILER_DEPENDENT_A

# elif defined(__GNUC__)
    ! g95 compiler, gfortran
#   define PP_CONF_ABORT           call abort

# elif defined(__INTEL_COMPILER)
    ! Intel fortran
#   define PP_CONF_USE_FOR_ABORT   use ifport
#   define PP_CONF_ABORT           call abort
# ifndef __i386
#   define PP_CONF_64BIT_CPU
# endif

# else
!!COMPILE_ERROR "Unknown platform"
# endif

/* EOF */
