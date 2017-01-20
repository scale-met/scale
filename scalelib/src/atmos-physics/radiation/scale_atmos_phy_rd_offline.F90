!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          for offline usage (input radiation flux from the file)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_rd_offline
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_offline_setup
  public :: ATMOS_PHY_RD_offline

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
  !> Setup
  subroutine ATMOS_PHY_RD_offline_setup( RD_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: RD_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[RADIATION] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Offline radiation process'

    if ( RD_TYPE /= 'OFFLINE' ) then
       write(*,*) 'xxx RD_TYPE is not OFFLINE. Check!'
       call PRC_MPIstop
    endif

    return
  end subroutine ATMOS_PHY_RD_offline_setup

  !-----------------------------------------------------------------------------
  !> Radiation main
  subroutine ATMOS_PHY_RD_offline( &
       DENS, RHOT, QTRC,      &
       CZ, FZ,                &
       fact_ocean,            &
       fact_land,             &
       fact_urban,            &
       temp_sfc, albedo_land, &
       solins, cosSZA,        &
       flux_rad,              &
       flux_rad_top,          &
       flux_rad_sfc_dn        )
!       Jval                   )
    use scale_grid_index
    use scale_tracer
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    real(RP), intent(in)  :: DENS           (KA,IA,JA)
    real(RP), intent(in)  :: RHOT           (KA,IA,JA)
    real(RP), intent(in)  :: QTRC           (KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ             (  KA,IA,JA)    ! UNUSED
    real(RP), intent(in)  :: FZ             (0:KA,IA,JA)
    real(RP), intent(in)  :: fact_ocean     (IA,JA)
    real(RP), intent(in)  :: fact_land      (IA,JA)
    real(RP), intent(in)  :: fact_urban     (IA,JA)
    real(RP), intent(in)  :: temp_sfc       (IA,JA)
    real(RP), intent(in)  :: albedo_land    (IA,JA,2)
    real(RP), intent(in)  :: solins         (IA,JA)
    real(RP), intent(in)  :: cosSZA         (IA,JA)
    real(RP), intent(out) :: flux_rad       (KA,IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_top   (IA,JA,2,2,2)
    real(RP), intent(out) :: flux_rad_sfc_dn(IA,JA,2,2)
!    real(RP), intent(out) :: Jval        (KA,IA,JA,CH_QA_photo)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Radiation(offline)'
    write(*,*) 'xxx This subroutine is never called. STOP.'
    call PRC_MPIstop

    return
  end subroutine ATMOS_PHY_RD_offline

end module scale_atmos_phy_rd_offline
