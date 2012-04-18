!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          DYCOMS-II Parametarized Radiative heating
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-26 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
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
  public :: ATMOS_PHY_RD_setup
  public :: ATMOS_PHY_RD

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  subroutine ATMOS_PHY_RD_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ DYCOMS-II Parametarized Radiative heating'

    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  ! Parametarized Radiative heating
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD
    use mod_const, only: &
       CPdry => CONST_CPdry
    use mod_time, only: &
       dtrd => TIME_DTSEC_ATMOS_PHY_RD
    use mod_grid, only : &
       CZ   => GRID_CZ,  &
       CDZ  => GRID_CDZ, &
       RCDZ => GRID_RCDZ
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_total,   &
       DENS, &
       RHOT, &
       QTRC
    implicit none

    real(8) :: RHOT_t  (KA,IA,JA) ! tendency rho*theta     [K*kg/m3/s]
    real(8) :: EFLX_rad(KA,IA,JA) ! Radiative heating flux [J/m2/s]
    real(8) :: Zi      (1 ,IA,JA) ! Cloud top height [m]

    real(8) :: QTOT ! Qv + Qc + Qr [kg/kg]

    real(8) :: Qbelow, Qabove ! scaled LWP (above/below layer)
    real(8) :: H              ! Heaviside Step Function

    real(8), parameter :: kappa = 85.0D0  ! scaling factor for LWP [m2/kg]
    real(8), parameter :: a     =  1.0D0  ! [K/m**-1/3]
    real(8), parameter :: F0    = 70.0D0  ! Upward [J/m2/s]
    real(8), parameter :: F1    = 22.0D0  ! [K/m**-1/3]
    real(8), parameter :: D     = 3.75D-6 ! divergence of large scale horizontal winds [1/s]

    real(8) :: dQ, dZ
    integer :: k_cldtop

    integer :: k, k2, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation'

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          EFLX_rad(k,i,j) = 0.D0
          RHOT_t  (k,i,j) = 0.D0
          Zi      (1,i,j) = 0.D0
       enddo

       ! diagnose cloud top
       k_cldtop = -1
       do k = KS, KE
          QTOT = QTRC(k,i,j,I_QV) &
               + QTRC(k,i,j,I_QC) &
               + QTRC(k,i,j,I_QR)

          if( QTOT < 8.D-3 ) exit ! above cloud
          k_cldtop = k
       enddo

       if( k_cldtop == -1 ) cycle ! no cloud

       Zi(1,i,j) = CZ(k_cldtop)

       do k = KS, KE

          Qbelow = 0.D0
          Qabove = 0.D0
          do k2 = KS, KE
             dQ = kappa * DENS(k2,i,j) * CDZ(k2) * ( QTRC(k2,i,j,I_QC) + QTRC(k2,i,j,I_QR) )

             if ( k2 <= k ) then ! below layer
                Qbelow = Qbelow + dQ
             else                ! above layer
                Qabove = Qabove + dQ
             endif
          enddo

          dZ = CZ(k)-CZ(k_cldtop)

          ! Heaviside Step Function
          if ( dZ > 0.D0 ) then
             H = 1.0D0
          elseif( dZ < 0.D0 ) then
             H = 0.0D0
          else
             H = 0.5D0
          endif

          EFLX_rad(k,i,j) = F0 * exp( -Qabove ) &
                          + F1 * exp( -Qbelow ) &
                          + a * DENS(k_cldtop,i,j) * CPdry * D * H &
                          * ( 0.25D0 * dZ  * dZ**(1.D0/3.D0) &
                            + CZ(k_cldtop) * dZ**(1.D0/3.D0) )

          RHOT_t(k,i,j) = EFLX_rad(k,i,j) / CPdry * RCDZ(k)

          RHOT(k,i,j) = RHOT(k,i,j) + RHOT_t(k,i,j) * dtrd
       enddo

    enddo
    enddo

    ! fill KHALO
    do j  = JS, JE
    do i  = IS, IE
       RHOT(   1:KS-1,i,j) = RHOT(KS,i,j)
       RHOT(KE+1:KA,  i,j) = RHOT(KE,i,j)
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars8( RHOT(:,:,:), 5 )
    call COMM_wait ( RHOT(:,:,:), 5 )

    call ATMOS_vars_total

    call HIST_in( EFLX_rad(:,:,:), 'EFLX_rd',   'Radiative heating flux', 'J/m2/s',    '3D', dtrd )
    call HIST_in( RHOT_t  (:,:,:), 'RHOT_t_rd', 'tendency of RHOT in rd', 'K*kg/m3/s', '3D', dtrd )
    call HIST_in( Zi      (:,:,:), 'Zi',        'Cloud top height',       'm',         '2D', dtrd )

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
