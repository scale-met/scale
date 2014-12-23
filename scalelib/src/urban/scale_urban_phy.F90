!-------------------------------------------------------------------------------
!> module URBAN / Physics
!!
!! @par Description
!!          urban physics module
!!
!! @author Team SCALE
!!
!! @par History
!<
!-------------------------------------------------------------------------------
module scale_urban_phy
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
  public :: URBAN_PHY_setup

  abstract interface
     subroutine urb( &
           TR_URB_t,    &
           TB_URB_t,    &
           TG_URB_t,    &
           TC_URB_t,    &
           QC_URB_t,    &
           UC_URB_t,    &
           TRL_URB_t,   &
           TBL_URB_t,   &
           TGL_URB_t,   &
           RAINR_URB_t, &
           RAINB_URB_t, &
           RAING_URB_t, &
           ROFF_URB_t,  &
           UST,         &
           ALBD_LW,     &
           ALBD_SW,     &
           MWFLX,       &
           MUFLX,       &
           MVFLX,       &
           SHFLX,       &
           LHFLX,       &
           GHFLX,       &
           Z0M,         &
           Z0H,         &
           Z0E,         &
           U10,         &
           V10,         &
           T2,          &
           Q2,          &
           TMPA,        &
           PRES,        &
           W1,          &
           U1,          &
           V1,          &
           DENS,        &
           QA,          &
           Z1,          &
           PBL,         &
           PRSS,        &
           LWD,         &
           SWD,         &
           PREC,        &
           TR_URB,      &
           TB_URB,      &
           TG_URB,      &
           TC_URB,      &
           QC_URB,      &
           UC_URB,      &
           TRL_URB,     &
           TBL_URB,     &
           TGL_URB,     &
           RAINR_URB,   &
           RAINB_URB,   &
           RAING_URB,   &
           ROFF_URB,    &
           LON,         &
           LAT,         &
           NOWDATE,     &
           dt           )
       use scale_precision
       use scale_grid_index
       use scale_urban_grid_index
       implicit none

       real(RP), intent(out) :: TR_URB_t   (IA,JA)
       real(RP), intent(out) :: TB_URB_t   (IA,JA)
       real(RP), intent(out) :: TG_URB_t   (IA,JA)
       real(RP), intent(out) :: TC_URB_t   (IA,JA)
       real(RP), intent(out) :: QC_URB_t   (IA,JA)
       real(RP), intent(out) :: UC_URB_t   (IA,JA)
       real(RP), intent(out) :: TRL_URB_t  (UKS:UKE,IA,JA)
       real(RP), intent(out) :: TBL_URB_t  (UKS:UKE,IA,JA)
       real(RP), intent(out) :: TGL_URB_t  (UKS:UKE,IA,JA)
       real(RP), intent(out) :: RAINR_URB_t(IA,JA)
       real(RP), intent(out) :: RAINB_URB_t(IA,JA)
       real(RP), intent(out) :: RAING_URB_t(IA,JA)
       real(RP), intent(out) :: ROFF_URB_t (IA,JA)

       real(RP), intent(out) :: UST    (IA,JA)
       real(RP), intent(out) :: ALBD_LW(IA,JA)
       real(RP), intent(out) :: ALBD_SW(IA,JA)
       real(RP), intent(out) :: MWFLX  (IA,JA)
       real(RP), intent(out) :: MUFLX  (IA,JA)
       real(RP), intent(out) :: MVFLX  (IA,JA)
       real(RP), intent(out) :: SHFLX  (IA,JA)
       real(RP), intent(out) :: LHFLX  (IA,JA)
       real(RP), intent(out) :: GHFLX  (IA,JA)
       real(RP), intent(out) :: Z0M    (IA,JA)
       real(RP), intent(out) :: Z0H    (IA,JA)
       real(RP), intent(out) :: Z0E    (IA,JA)
       real(RP), intent(out) :: U10    (IA,JA)
       real(RP), intent(out) :: V10    (IA,JA)
       real(RP), intent(out) :: T2     (IA,JA)
       real(RP), intent(out) :: Q2     (IA,JA)

       real(RP), intent(in) :: TMPA  (IA,JA)
       real(RP), intent(in) :: PRES  (IA,JA)
       real(RP), intent(in) :: W1    (IA,JA)
       real(RP), intent(in) :: U1    (IA,JA)
       real(RP), intent(in) :: V1    (IA,JA)
       real(RP), intent(in) :: DENS  (IA,JA)
       real(RP), intent(in) :: QA    (IA,JA)
       real(RP), intent(in) :: Z1    (IA,JA)
       real(RP), intent(in) :: PBL   (IA,JA)
       real(RP), intent(in) :: PRSS  (IA,JA)
       real(RP), intent(in) :: LWD   (IA,JA)
       real(RP), intent(in) :: SWD   (IA,JA)
       real(RP), intent(in) :: PREC  (IA,JA)

       real(RP), intent(in) :: TR_URB   (IA,JA)
       real(RP), intent(in) :: TB_URB   (IA,JA)
       real(RP), intent(in) :: TG_URB   (IA,JA)
       real(RP), intent(in) :: TC_URB   (IA,JA)
       real(RP), intent(in) :: QC_URB   (IA,JA)
       real(RP), intent(in) :: UC_URB   (IA,JA)
       real(RP), intent(in) :: TRL_URB  (UKS:UKE,IA,JA)
       real(RP), intent(in) :: TBL_URB  (UKS:UKE,IA,JA)
       real(RP), intent(in) :: TGL_URB  (UKS:UKE,IA,JA)
       real(RP), intent(in) :: RAINR_URB(IA,JA)
       real(RP), intent(in) :: RAINB_URB(IA,JA)
       real(RP), intent(in) :: RAING_URB(IA,JA)
       real(RP), intent(in) :: ROFF_URB (IA,JA)

       real(RP), intent(in) :: LON
       real(RP), intent(in) :: LAT
       integer,  intent(in) :: NOWDATE(6)
       real(RP), intent(in) :: dt
     end subroutine urb
  end interface
  procedure(urb), pointer :: URBAN_PHY => NULL()
  public :: URBAN_PHY

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

  subroutine URBAN_PHY_setup( URBAN_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    use scale_urban_phy_slc, only: &
       URBAN_PHY_SLC_setup, &
       URBAN_PHY_SLC
    implicit none

    character(len=*), intent(in) :: URBAN_TYPE
    !---------------------------------------------------------------------------

    select case( URBAN_TYPE )
    case ( 'SLC' )
       call URBAN_PHY_SLC_setup( URBAN_TYPE )
       URBAN_PHY => URBAN_PHY_SLC
    case default
       write(*,*) 'xxx invalid Urban type(', trim(URBAN_TYPE), '). CHECK!'
       call PRC_MPIstop
    end select

  end subroutine URBAN_PHY_setup

end module scale_urban_phy
