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
!! @li      2012-03-26 (H.Yashiro)  [new]
!! @li      2012-06-10 (Y.Miyamoto) [mod] bug-fix and some modification
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
  include 'inc_precision.h'
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
  real(RP), private, save :: F0    = 70.00_RP  ! Upward [J/m2/s]
  real(RP), private, save :: F1    = 22.00_RP  ! [K/m**-1/3]
  real(RP), private, save :: Dval  = 3.75E-6_RP ! divergence of large scale horizontal winds [1/s]
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_RD
    implicit none

    real(RP) :: ATMOS_RD_F0
    real(RP) :: ATMOS_RD_F1
    real(RP) :: ATMOS_RD_Dval
  
    integer :: ierr

    NAMELIST / PARAM_ATMOS_PHY_RD / &
       ATMOS_RD_F0, &
       ATMOS_RD_F1, &
       ATMOS_RD_Dval
    !---------------------------------------------------------------------------

    ATMOS_RD_F0 = F0
    ATMOS_RD_F1 = F1
    ATMOS_RD_Dval = Dval

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ DYCOMS-II Parametarized Radiative heating'

    if ( trim(ATMOS_TYPE_PHY_RD) .ne. 'DYCOMSII' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_RD is not DYCOMSII. Check!'
       call PRC_MPIstop
    end if


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD,iostat=ierr)

    if( ierr < 0 ) then !--- missing
      if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD. Check!'
      call PRC_MPIstop
    endif 
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_RD)

    F0 = ATMOS_RD_F0
    F1 = ATMOS_RD_F1
    Dval = ATMOS_RD_Dval

    return
  end subroutine ATMOS_PHY_RD_setup

  !-----------------------------------------------------------------------------
  ! Parametarized Radiative heating
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD
    use mod_const, only: &
       Rdry    => CONST_Rdry,   &
       CPdry   => CONST_CPdry,  &
       RovCP   => CONST_RovCP,  &
       RovCV   => CONST_RovCV,  &
       CPovCV  => CONST_CPovCV, &
       EPSTvap => CONST_CPdry,  &
       P00     => CONST_PRE00
    use mod_time, only: &
       dtrd => TIME_DTSEC_ATMOS_PHY_RD
    use mod_grid, only : &
       CZ   => GRID_CZ,  &
       FZ   => GRID_FZ,  &
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

    real(RP) :: TEMP_t  (KA,IA,JA) ! tendency rho*theta     [K*kg/m3/s]
    real(RP) :: EFLX_rad(KA,IA,JA) ! Radiative heating flux [J/m2/s]
    real(RP) :: Zi      (1 ,IA,JA) ! Cloud top height [m]

    real(RP) :: QTOT ! Qv + Qc + Qr [kg/kg]

    real(RP) :: Qbelow, Qabove ! scaled LWP (above/below layer)

    real(RP), parameter :: kappa = 85.00_RP  ! scaling factor for LWP [m2/kg]
    real(RP), parameter :: a     =  1.00_RP  ! [K/m**-1/3]

    real(RP) :: pres

    real(RP) :: dQ, dZ, QWSUM
    integer :: k_cldtop

    integer :: k, k2, i, j, iq

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation'

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          EFLX_rad(k,i,j) = 0.0_RP
          TEMP_t  (k,i,j) = 0.0_RP
          Zi      (1,i,j) = 0.0_RP
       enddo

       ! diagnose cloud top
       k_cldtop = -1
       do k = KS, KE
           QTOT = 0.0_RP
           do iq = QQS, QWE
              QTOT = QTOT + QTRC(k,i,j,iq)
           enddo 
          if( QTOT < 8.E-3_RP ) exit ! above cloud
          k_cldtop = k
       enddo

       if( k_cldtop == -1 ) cycle ! no cloud

       Zi(1,i,j) = CZ(k_cldtop)

       do k = KS-1, KE

          Qbelow = 0.0_RP
          Qabove = 0.0_RP
          do k2 = KS, KE
!             dQ = kappa * CDZ(k2) * DENS(k2,i,j) * ( QTRC(k2,i,j,I_QC) + QTRC(k2,i,j,I_QR) )
             QWSUM = 0.0_RP
             do iq = QWS, QWE
                QWSUM = QWSUM + QTRC(k2,i,j,iq)
             enddo
             dQ = kappa * CDZ(k2) * DENS(k2,i,j) * QWSUM
             if ( k2 <= k ) then ! below layer
                Qbelow = Qbelow + dQ
             else                ! above layer
                Qabove = Qabove + dQ
             endif
          enddo

          EFLX_rad(k,i,j) = F0 * exp( -Qabove ) &
                          + F1 * exp( -Qbelow )
       enddo

       do k = k_cldtop, KE
          dZ = FZ(k)-CZ(k_cldtop)
          QTOT = 0.0_RP
          do iq = QQS, QWE
             QTOT = QTOT + QTRC(k,i,j,iq)
          enddo 
          EFLX_rad(k,i,j) = EFLX_rad(k,i,j) &
                          + a * DENS(k_cldtop,i,j)*( 1.0_RP-QTOT ) * CPdry * Dval &
                          * ( 0.250_RP * dZ  * dZ**(1.0_RP/3.0_RP) &
                            + CZ(k_cldtop) * dZ**(1.0_RP/3.0_RP) )
       enddo

       do k = KS, KE
          TEMP_t(k,i,j) = &
             - ( EFLX_rad(k,i,j) - EFLX_rad(k-1,i,j) ) / CPdry * RCDZ(k) !/ DENS(k,i,j)

          RHOT(k,i,j) = RHOT(k,i,j) + dtrd * ( 1.0_RP - RovCP ) &
                      * ( P00/(RHOT(k,i,j)*Rdry) )**RovCV * TEMP_t(k,i,j)!*DENS(k,i,j)
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

    call HIST_in( EFLX_rad(:,:,:), 'EFLX_rd',   'Radiative heating flux', 'J/m2/s',    dtrd )
    call HIST_in( TEMP_t  (:,:,:), 'TEMP_t_rd', 'tendency of temp in rd', 'K*kg/m3/s', dtrd )
    call HIST_in( Zi      (1,:,:), 'Zi',        'Cloud top height',       'm',         dtrd )

    return
  end subroutine ATMOS_PHY_RD

end module mod_atmos_phy_rd
