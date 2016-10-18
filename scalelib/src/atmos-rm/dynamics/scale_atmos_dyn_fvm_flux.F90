
!-------------------------------------------------------------------------------
!> module scale_atmos_dyn_fvm_flux
!!
!! @par Description
!!          FVM flux scheme
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-04-18 (S.Nishizawa) [new]
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_dyn_fvm_flux
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_index
  use scale_tracer
  use scale_process
#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_FVM_flux_setup

  abstract interface
     subroutine valueW( &
          valW, &
          mflx, val, GSQRT, &
          CDZ )
       use scale_precision
       use scale_grid_index
       implicit none
       real(RP), intent(out) :: valW (KA)
       real(RP), intent(in)  :: mflx (KA)
       real(RP), intent(in)  :: val  (KA)
       real(RP), intent(in)  :: GSQRT(KA)
       real(RP), intent(in)  :: CDZ(KA)
     end subroutine valueW
     subroutine flux_phi( &
          flux, &
          mflx, val, GSQRT, &
          num_diff, &
          CDZ, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_grid_index
       implicit none
       real(RP), intent(out) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mflx    (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_phi
     subroutine flux_mom( &
          flux, &
          mom, val, DENS, &
          GSQRT, MAPF, &
          num_diff, &
          CDZ, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_grid_index
       implicit none
       real(RP), intent(out) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: MAPF    (   IA,JA,2)
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_mom
     subroutine flux_z( &
          flux, &
          mom, val, DENS, &
          GSQRT, J33G, &
          num_diff, &
          CDZ, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_grid_index
       implicit none
       real(RP), intent(out) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: J33G
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_z
     subroutine flux_wz( &
          flux, &
          mom, val, DENS, &
          GSQRT, J33G, &
          num_diff, &
          CDZ, FDZ, &
          dtrk, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_grid_index
       implicit none
       real(RP), intent(out) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: J33G
       real(RP), intent(in)  :: num_diff(KA,IA,JA)
       real(RP), intent(in)  :: CDZ(KA)
       real(RP), intent(in)  :: FDZ(KA-1)
       real(RP), intent(in)  :: dtrk
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_wz
     subroutine flux_j( &
          flux, &
          mom, val, DENS, &
          GSQRT, JG, MAPF, &
          CDZ, &
          IIS, IIE, JJS, JJE )
       use scale_precision
       use scale_grid_index
       implicit none
       real(RP), intent(out) :: flux    (KA,IA,JA)
       real(RP), intent(in)  :: mom     (KA,IA,JA)
       real(RP), intent(in)  :: val     (KA,IA,JA)
       real(RP), intent(in)  :: DENS    (KA,IA,JA)
       real(RP), intent(in)  :: GSQRT   (KA,IA,JA)
       real(RP), intent(in)  :: JG      (KA,IA,JA)
       real(RP), intent(in)  :: MAPF    (   IA,JA,2)
       real(RP), intent(in)  :: CDZ(KA)
       integer,  intent(in)  :: IIS, IIE, JJS, JJE
     end subroutine flux_j
  end interface

  procedure(valueW), pointer :: ATMOS_DYN_FVM_flux_valueW_Z => NULL()
  public :: ATMOS_DYN_FVM_flux_valueW_Z


  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxZ_XYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYZ
  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxZ_XYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYZ_tracer

  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxX_XYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYZ
  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxX_XYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYZ_tracer

  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxY_XYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYZ
  procedure(flux_phi), pointer :: ATMOS_DYN_FVM_fluxY_XYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYZ_tracer



  procedure(flux_wz), pointer :: ATMOS_DYN_FVM_fluxZ_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYW
  procedure(flux_wz), pointer :: ATMOS_DYN_FVM_fluxZ_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XYW_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYW
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XYW_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYW
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XYW_tracer


  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_XYW
  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_XYW_tracer

  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_XYW => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_XYW
  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_XYW_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_XYW_tracer


  procedure(flux_z), pointer :: ATMOS_DYN_FVM_fluxZ_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_UYZ
  procedure(flux_z), pointer :: ATMOS_DYN_FVM_fluxZ_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_UYZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxX_UYZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_UYZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxY_UYZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_UYZ_tracer


  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_UYZ
  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_UYZ_tracer

  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_UYZ => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_UYZ
  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_UYZ_tracer


  procedure(flux_z), pointer :: ATMOS_DYN_FVM_fluxZ_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XVZ
  procedure(flux_z), pointer :: ATMOS_DYN_FVM_fluxZ_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxZ_XVZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XVZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxX_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxX_XVZ_tracer

  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XVZ
  procedure(flux_mom), pointer :: ATMOS_DYN_FVM_fluxY_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxY_XVZ_tracer


  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_XVZ
  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ13_XVZ_tracer

  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_XVZ => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_XVZ
  procedure(flux_j), pointer :: ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => NULL()
  public :: ATMOS_DYN_FVM_fluxJ23_XVZ_tracer



  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  
  !-----------------------------------------------------------------------------
  !> setup
  subroutine ATMOS_DYN_FVM_flux_setup( &
       scheme, scheme_tracer )
    use scale_process, only: &
         PRC_MPIstop

   use scale_atmos_dyn_fvm_flux_ud1wrap, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud1wrap, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud1wrap, &
      ATMOS_DYN_FVM_fluxX_XYW_ud1wrap, &
      ATMOS_DYN_FVM_fluxY_XYW_ud1wrap, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud1wrap, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud1wrap, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud1wrap, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud1wrap

   use scale_atmos_dyn_fvm_flux_cd2, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd2, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd2, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd2, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd2, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd2, &
      ATMOS_DYN_FVM_fluxX_XYW_cd2, &
      ATMOS_DYN_FVM_fluxY_XYW_cd2, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd2, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd2, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd2, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd2, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd2

   use scale_atmos_dyn_fvm_flux_ud3, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud3, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud3, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud3, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud3, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud3, &
      ATMOS_DYN_FVM_fluxX_XYW_ud3, &
      ATMOS_DYN_FVM_fluxY_XYW_ud3, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud3, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud3, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud3, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud3, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud3

   use scale_atmos_dyn_fvm_flux_ud3Koren1993, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud3Koren1993, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud3Koren1993

   use scale_atmos_dyn_fvm_flux_cd4, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd4, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd4, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd4, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd4, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd4, &
      ATMOS_DYN_FVM_fluxX_XYW_cd4, &
      ATMOS_DYN_FVM_fluxY_XYW_cd4, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd4, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd4, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd4, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd4, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd4

   use scale_atmos_dyn_fvm_flux_ud5, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_ud5, &
      ATMOS_DYN_FVM_fluxZ_XYZ_ud5, &
      ATMOS_DYN_FVM_fluxX_XYZ_ud5, &
      ATMOS_DYN_FVM_fluxY_XYZ_ud5, &
      ATMOS_DYN_FVM_fluxZ_XYW_ud5, &
      ATMOS_DYN_FVM_fluxX_XYW_ud5, &
      ATMOS_DYN_FVM_fluxY_XYW_ud5, &
      ATMOS_DYN_FVM_fluxJ13_XYW_ud5, &
      ATMOS_DYN_FVM_fluxJ23_XYW_ud5, &
      ATMOS_DYN_FVM_fluxZ_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxX_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxY_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_ud5, &
      ATMOS_DYN_FVM_fluxZ_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxX_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxY_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_ud5, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_ud5

   use scale_atmos_dyn_fvm_flux_cd6, only: &
      ATMOS_DYN_FVM_flux_valueW_Z_cd6, &
      ATMOS_DYN_FVM_fluxZ_XYZ_cd6, &
      ATMOS_DYN_FVM_fluxX_XYZ_cd6, &
      ATMOS_DYN_FVM_fluxY_XYZ_cd6, &
      ATMOS_DYN_FVM_fluxZ_XYW_cd6, &
      ATMOS_DYN_FVM_fluxX_XYW_cd6, &
      ATMOS_DYN_FVM_fluxY_XYW_cd6, &
      ATMOS_DYN_FVM_fluxJ13_XYW_cd6, &
      ATMOS_DYN_FVM_fluxJ23_XYW_cd6, &
      ATMOS_DYN_FVM_fluxZ_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxX_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxY_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxJ13_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxJ23_UYZ_cd6, &
      ATMOS_DYN_FVM_fluxZ_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxX_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxY_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxJ13_XVZ_cd6, &
      ATMOS_DYN_FVM_fluxJ23_XVZ_cd6

    implicit none
    character(len=*), intent(in) :: scheme
    character(len=*), intent(in) :: scheme_tracer

    select case ( scheme )

    case ( "UD1" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud1 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud1wrap

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud1wrap

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud1wrap

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud1wrap

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_ud1wrap



      if ( IHALO < 1 ) then
         write(*,*) 'xxx IHALO must be >= ', 1
         call PRC_MPIstop
      end if
      if ( JHALO < 1 ) then
         write(*,*) 'xxx JHALO must be >= ', 1
         call PRC_MPIstop
      end if


    case ( "CD2" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the cd2 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd2

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd2

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd2

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd2

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_cd2

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd2

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd2

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_cd2

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_cd2

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_cd2

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd2

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd2

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_cd2

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_cd2

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_cd2

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd2

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd2

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_cd2

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_cd2



      if ( IHALO < 1 ) then
         write(*,*) 'xxx IHALO must be >= ', 1
         call PRC_MPIstop
      end if
      if ( JHALO < 1 ) then
         write(*,*) 'xxx JHALO must be >= ', 1
         call PRC_MPIstop
      end if


    case ( "UD3" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud3 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud3

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud3

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud3

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud3

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_ud3

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud3

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud3

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_ud3

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_ud3

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_ud3

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud3

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud3

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_ud3

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_ud3

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_ud3

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud3

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud3

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_ud3

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_ud3



      if ( IHALO < 2 ) then
         write(*,*) 'xxx IHALO must be >= ', 2
         call PRC_MPIstop
      end if
      if ( JHALO < 2 ) then
         write(*,*) 'xxx JHALO must be >= ', 2
         call PRC_MPIstop
      end if


    case ( "UD3KOREN1993" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud3Koren1993 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_ud3Koren1993



      if ( IHALO < 2 ) then
         write(*,*) 'xxx IHALO must be >= ', 2
         call PRC_MPIstop
      end if
      if ( JHALO < 2 ) then
         write(*,*) 'xxx JHALO must be >= ', 2
         call PRC_MPIstop
      end if


    case ( "CD4" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the cd4 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd4

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd4

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd4

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd4

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_cd4

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd4

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd4

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_cd4

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_cd4

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_cd4

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd4

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd4

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_cd4

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_cd4

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_cd4

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd4

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd4

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_cd4

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_cd4



      if ( IHALO < 2 ) then
         write(*,*) 'xxx IHALO must be >= ', 2
         call PRC_MPIstop
      end if
      if ( JHALO < 2 ) then
         write(*,*) 'xxx JHALO must be >= ', 2
         call PRC_MPIstop
      end if


    case ( "UD5" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud5 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_ud5

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_ud5

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_ud5

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_ud5

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_ud5

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_ud5

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_ud5

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_ud5

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_ud5

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_ud5

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_ud5

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_ud5

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_ud5

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_ud5

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_ud5

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_ud5

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_ud5

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_ud5

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_ud5



      if ( IHALO < 3 ) then
         write(*,*) 'xxx IHALO must be >= ', 3
         call PRC_MPIstop
      end if
      if ( JHALO < 3 ) then
         write(*,*) 'xxx JHALO must be >= ', 3
         call PRC_MPIstop
      end if


    case ( "CD6" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the cd6 scheme is used for flux calculation'

      ATMOS_DYN_FVM_flux_valueW_Z => ATMOS_DYN_FVM_flux_valueW_Z_cd6

      ATMOS_DYN_FVM_fluxZ_XYZ => ATMOS_DYN_FVM_fluxZ_XYZ_cd6

      ATMOS_DYN_FVM_fluxX_XYZ => ATMOS_DYN_FVM_fluxX_XYZ_cd6

      ATMOS_DYN_FVM_fluxY_XYZ => ATMOS_DYN_FVM_fluxY_XYZ_cd6

      ATMOS_DYN_FVM_fluxZ_XYW => ATMOS_DYN_FVM_fluxZ_XYW_cd6

      ATMOS_DYN_FVM_fluxX_XYW => ATMOS_DYN_FVM_fluxX_XYW_cd6

      ATMOS_DYN_FVM_fluxY_XYW => ATMOS_DYN_FVM_fluxY_XYW_cd6

      ATMOS_DYN_FVM_fluxJ13_XYW => ATMOS_DYN_FVM_fluxJ13_XYW_cd6

      ATMOS_DYN_FVM_fluxJ23_XYW => ATMOS_DYN_FVM_fluxJ23_XYW_cd6

      ATMOS_DYN_FVM_fluxZ_UYZ => ATMOS_DYN_FVM_fluxZ_UYZ_cd6

      ATMOS_DYN_FVM_fluxX_UYZ => ATMOS_DYN_FVM_fluxX_UYZ_cd6

      ATMOS_DYN_FVM_fluxY_UYZ => ATMOS_DYN_FVM_fluxY_UYZ_cd6

      ATMOS_DYN_FVM_fluxJ13_UYZ => ATMOS_DYN_FVM_fluxJ13_UYZ_cd6

      ATMOS_DYN_FVM_fluxJ23_UYZ => ATMOS_DYN_FVM_fluxJ23_UYZ_cd6

      ATMOS_DYN_FVM_fluxZ_XVZ => ATMOS_DYN_FVM_fluxZ_XVZ_cd6

      ATMOS_DYN_FVM_fluxX_XVZ => ATMOS_DYN_FVM_fluxX_XVZ_cd6

      ATMOS_DYN_FVM_fluxY_XVZ => ATMOS_DYN_FVM_fluxY_XVZ_cd6

      ATMOS_DYN_FVM_fluxJ13_XVZ => ATMOS_DYN_FVM_fluxJ13_XVZ_cd6

      ATMOS_DYN_FVM_fluxJ23_XVZ => ATMOS_DYN_FVM_fluxJ23_XVZ_cd6



      if ( IHALO < 3 ) then
         write(*,*) 'xxx IHALO must be >= ', 3
         call PRC_MPIstop
      end if
      if ( JHALO < 3 ) then
         write(*,*) 'xxx JHALO must be >= ', 3
         call PRC_MPIstop
      end if


    case default
       write(*,*) 'xxx scheme is invalid: ', scheme
       call PRC_MPIstop
    end select

    select case ( scheme_tracer )

    case ( "UD1" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud1 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud1wrap

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud1wrap

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud1wrap

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud1wrap

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud1wrap

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud1wrap

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud1wrap



      if ( IHALO < 1 ) then
         write(*,*) 'xxx IHALO must be >= ', 1
         call PRC_MPIstop
      end if
      if ( JHALO < 1 ) then
         write(*,*) 'xxx JHALO must be >= ', 1
         call PRC_MPIstop
      end if


    case ( "CD2" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the cd2 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd2

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd2

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd2

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd2

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd2

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd2

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd2

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd2

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd2

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd2

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd2

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd2

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd2

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd2

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd2

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd2

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd2

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd2



      if ( IHALO < 1 ) then
         write(*,*) 'xxx IHALO must be >= ', 1
         call PRC_MPIstop
      end if
      if ( JHALO < 1 ) then
         write(*,*) 'xxx JHALO must be >= ', 1
         call PRC_MPIstop
      end if


    case ( "UD3" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud3 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud3

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud3

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud3

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud3

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud3

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud3

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud3

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud3

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud3

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud3

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud3

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud3

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud3

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud3

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud3

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud3

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud3

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud3



      if ( IHALO < 2 ) then
         write(*,*) 'xxx IHALO must be >= ', 2
         call PRC_MPIstop
      end if
      if ( JHALO < 2 ) then
         write(*,*) 'xxx JHALO must be >= ', 2
         call PRC_MPIstop
      end if


    case ( "UD3KOREN1993" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud3Koren1993 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud3Koren1993

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud3Koren1993



      if ( IHALO < 2 ) then
         write(*,*) 'xxx IHALO must be >= ', 2
         call PRC_MPIstop
      end if
      if ( JHALO < 2 ) then
         write(*,*) 'xxx JHALO must be >= ', 2
         call PRC_MPIstop
      end if


    case ( "CD4" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the cd4 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd4

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd4

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd4

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd4

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd4

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd4

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd4

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd4

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd4

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd4

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd4

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd4

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd4

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd4

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd4

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd4

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd4

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd4



      if ( IHALO < 2 ) then
         write(*,*) 'xxx IHALO must be >= ', 2
         call PRC_MPIstop
      end if
      if ( JHALO < 2 ) then
         write(*,*) 'xxx JHALO must be >= ', 2
         call PRC_MPIstop
      end if


    case ( "UD5" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the ud5 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_ud5

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_ud5

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_ud5

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_ud5

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_ud5

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_ud5

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_ud5

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_ud5

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_ud5

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_ud5

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_ud5

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_ud5

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_ud5

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_ud5

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_ud5

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_ud5

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_ud5

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_ud5



      if ( IHALO < 3 ) then
         write(*,*) 'xxx IHALO must be >= ', 3
         call PRC_MPIstop
      end if
      if ( JHALO < 3 ) then
         write(*,*) 'xxx JHALO must be >= ', 3
         call PRC_MPIstop
      end if


    case ( "CD6" )
      if( IO_L ) write(IO_FID_LOG,*) '*** the cd6 scheme is used for flux calculation of tracer'

      ATMOS_DYN_FVM_fluxZ_XYZ_tracer => ATMOS_DYN_FVM_fluxZ_XYZ_cd6

      ATMOS_DYN_FVM_fluxX_XYZ_tracer => ATMOS_DYN_FVM_fluxX_XYZ_cd6

      ATMOS_DYN_FVM_fluxY_XYZ_tracer => ATMOS_DYN_FVM_fluxY_XYZ_cd6

      ATMOS_DYN_FVM_fluxZ_XYW_tracer => ATMOS_DYN_FVM_fluxZ_XYW_cd6

      ATMOS_DYN_FVM_fluxX_XYW_tracer => ATMOS_DYN_FVM_fluxX_XYW_cd6

      ATMOS_DYN_FVM_fluxY_XYW_tracer => ATMOS_DYN_FVM_fluxY_XYW_cd6

      ATMOS_DYN_FVM_fluxJ13_XYW_tracer => ATMOS_DYN_FVM_fluxJ13_XYW_cd6

      ATMOS_DYN_FVM_fluxJ23_XYW_tracer => ATMOS_DYN_FVM_fluxJ23_XYW_cd6

      ATMOS_DYN_FVM_fluxZ_UYZ_tracer => ATMOS_DYN_FVM_fluxZ_UYZ_cd6

      ATMOS_DYN_FVM_fluxX_UYZ_tracer => ATMOS_DYN_FVM_fluxX_UYZ_cd6

      ATMOS_DYN_FVM_fluxY_UYZ_tracer => ATMOS_DYN_FVM_fluxY_UYZ_cd6

      ATMOS_DYN_FVM_fluxJ13_UYZ_tracer => ATMOS_DYN_FVM_fluxJ13_UYZ_cd6

      ATMOS_DYN_FVM_fluxJ23_UYZ_tracer => ATMOS_DYN_FVM_fluxJ23_UYZ_cd6

      ATMOS_DYN_FVM_fluxZ_XVZ_tracer => ATMOS_DYN_FVM_fluxZ_XVZ_cd6

      ATMOS_DYN_FVM_fluxX_XVZ_tracer => ATMOS_DYN_FVM_fluxX_XVZ_cd6

      ATMOS_DYN_FVM_fluxY_XVZ_tracer => ATMOS_DYN_FVM_fluxY_XVZ_cd6

      ATMOS_DYN_FVM_fluxJ13_XVZ_tracer => ATMOS_DYN_FVM_fluxJ13_XVZ_cd6

      ATMOS_DYN_FVM_fluxJ23_XVZ_tracer => ATMOS_DYN_FVM_fluxJ23_XVZ_cd6



      if ( IHALO < 3 ) then
         write(*,*) 'xxx IHALO must be >= ', 3
         call PRC_MPIstop
      end if
      if ( JHALO < 3 ) then
         write(*,*) 'xxx JHALO must be >= ', 3
         call PRC_MPIstop
      end if


    case default
       write(*,*) 'xxx scheme is invalid: ', scheme_tracer
       call PRC_MPIstop
    end select

  end subroutine ATMOS_DYN_FVM_flux_setup

end module scale_atmos_dyn_fvm_flux

!--
! vi:set readonly sw=4 ts=8
!
!Local Variables:
!mode: f90
!buffer-read-only: t
!End:
!
!++
