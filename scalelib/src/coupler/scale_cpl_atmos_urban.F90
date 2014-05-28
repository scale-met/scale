!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Urban Surface fluxes
!!
!! @par Description
!!          Surface flux from the atmosphere-urban coupler
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module scale_cpl_atmos_urban
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  use scale_urban_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: CPL_AtmUrb_setup

  abstract interface
     subroutine cau( &
        TR_URB,  & ! (inout)
        TB_URB,  & ! (inout)
        TG_URB,  & ! (inout)
        TC_URB,  & ! (inout)
        QC_URB,  & ! (inout)
        UC_URB,  & ! (inout)
        TRL_URB, & ! (inout)
        TBL_URB, & ! (inout)
        TGL_URB, & ! (inout)
        UST,     & ! (out)
        SHFLX,   & ! (out)
        LHFLX,   & ! (out)
        GHFLX,   & ! (out)
        LSOLAR,  & ! (in)
        TMPA,    & ! (in)
        QA,      & ! (in)
        UA,      & ! (in)
        U1,      & ! (in)
        V1,      & ! (in)
        Z1,      & ! (in)
        SWD,     & ! (in)
        LWD,     & ! (in)
        PREC,    & ! (in)
        DENS,    & ! (in)
        LON,     & ! (in)
        LAT      ) ! (in)
       use scale_precision
       use scale_grid_index
       use scale_urban_grid_index
       implicit none

       real(RP), intent(inout) :: TR_URB
       real(RP), intent(inout) :: TB_URB
       real(RP), intent(inout) :: TG_URB
       real(RP), intent(inout) :: TC_URB
       real(RP), intent(inout) :: QC_URB
       real(RP), intent(inout) :: UC_URB
       real(RP), intent(inout) :: TRL_URB(UKS:UKE)
       real(RP), intent(inout) :: TBL_URB(UKS:UKE)
       real(RP), intent(inout) :: TGL_URB(UKS:UKE)

       real(RP), intent(out) :: UST
       real(RP), intent(out) :: SHFLX
       real(RP), intent(out) :: LHFLX
       real(RP), intent(out) :: GHFLX

       logical,  intent(in) :: LSOLAR
       real(RP), intent(in) :: TMPA
       real(RP), intent(in) :: QA
       real(RP), intent(in) :: UA
       real(RP), intent(in) :: U1
       real(RP), intent(in) :: V1
       real(RP), intent(in) :: Z1
       real(RP), intent(in) :: SWD
       real(RP), intent(in) :: LWD
       real(RP), intent(in) :: PREC
       real(RP), intent(in) :: DENS
       real(RP), intent(in) :: LON
       real(RP), intent(in) :: LAT
     end subroutine cau
  end interface
  procedure(cau), pointer :: CPL_AtmUrb => NULL()
  public :: CPL_AtmUrb

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

  subroutine CPL_AtmUrb_setup( CPL_TYPE_AtmUrb )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef CAL
    use NAME(scale_cpl_atmos_urban_, CAU,), only: &
       NAME(CPL_AtmUrb_, CAU, _setup), &
       NAME(CPL_AtmUrb_, CAU,)
#else
    use scale_cpl_atmos_urban_bulk, only: &
       CPL_AtmUrb_bulk_setup, &
       CPL_AtmUrb_bulk
#endif
    implicit none

    character(len=*), intent(in) :: CPL_TYPE_AtmUrb
    !---------------------------------------------------------------------------

    select case( CPL_TYPE_AtmUrb )
    case ( 'BULK' )
       call CPL_AtmUrb_bulk_setup( CPL_TYPE_AtmUrb )
       CPL_AtmUrb => CPL_AtmUrb_bulk
    end select

  end subroutine CPL_AtmUrb_setup

end module scale_cpl_atmos_urban
