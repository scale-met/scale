!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)  [new]
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf_wrap
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_index
  use mod_tracer
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_wrap_setup
  public :: ATMOS_PHY_SF_wrap

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

  ! surface flux
  real(RP), private, allocatable :: SFLX_MOMZ(:,:)
  real(RP), private, allocatable :: SFLX_MOMX(:,:)
  real(RP), private, allocatable :: SFLX_MOMY(:,:)
  real(RP), private, allocatable :: SFLX_POTT(:,:)
  real(RP), private, allocatable :: SFLX_QV(:,:)


  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_wrap_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACEFLUX]/Categ[ATMOS]'

    allocate( SFLX_MOMZ(IA,JA) )
    allocate( SFLX_MOMX(IA,JA) )
    allocate( SFLX_MOMY(IA,JA) )
    allocate( SFLX_POTT(IA,JA) )
    allocate( SFLX_QV(IA,JA) )

    call ATMOS_PHY_sf_init()
    call ATMOS_PHY_sf_wrap( .true., .false. )

    return
  end subroutine ATMOS_PHY_SF_wrap_setup

  !-----------------------------------------------------------------------------
  ! calculation flux
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_wrap( &
       update_flag, &
       history_flag &
       )
    use mod_time, only: &
       dtsf => TIME_DTSEC_ATMOS_PHY_SF, &
       NOWSEC => TIME_NOWDAYSEC
    use mod_const, only: &
       CPdry  => CONST_CPdry,  &
       LH0    => CONST_LH0
    use mod_history, only: &
       HIST_in
    use mod_grid, only: &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       DENS_tp, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag

    ! monitor
    real(RP) :: SHFLX(IA,JA) ! sensible heat flux [W/m2]
    real(RP) :: LHFLX(IA,JA) ! latent   heat flux [W/m2]

    integer :: i, j

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux'

    if ( update_flag ) then
       call ATMOS_PHY_SF_main( &
            SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
            DENS, MOMZ, MOMX, MOMY,                              & ! (in)
            NOWSEC                                               ) ! (out)

       !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
       do j = JS, JE
       do i = IS, IE
          SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
          LHFLX(i,j) = SFLX_QV  (i,j) * LH0
       end do
       end do

       if ( present(history_flag) ) then
       if ( history_flag ) then
          call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2', dtsf )
          call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2', dtsf )
       end if
       end if

    end if

    !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SFLX_POTT(i,j) &
              + SFLX_QV(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
            + SFLX_QV(i,j) * RCDZ(KS)
       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
            + SFLX_MOMZ(i,j) * RFDZ(KS)
       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
            + SFLX_MOMX(i,j) * RCDZ(KS)
       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
            + SFLX_MOMY(i,j) * RCDZ(KS)
       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
            + SFLX_QV(i,j) * RCDZ(KS)
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sf
