!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process
!!          DYCOMS-II Parametarized Radiative heating
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-03-26 (H.Yashiro)  [new]
!! @li      2012-06-10 (Y.Miyamoto) [mod] bug-fix and some modification
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_rd_dycoms2
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_dycoms2_setup
  public :: ATMOS_PHY_RD_dycoms2

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

  logical, private :: first = .true.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_dycoms2_setup( RD_TYPE )
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG, &
       IO_L, &
       IO_SYSCHR
    use mod_process, only: &
       PRC_MPIstop
    implicit none
    character(len=IO_SYSCHR), intent(in) :: RD_TYPE


    real(RP) :: ATMOS_RD_F0
    real(RP) :: ATMOS_RD_F1
    real(RP) :: ATMOS_RD_Dval

    NAMELIST / PARAM_ATMOS_PHY_RD_DYCOMS2 / &
       ATMOS_RD_F0, &
       ATMOS_RD_F1, &
       ATMOS_RD_Dval

    integer :: ierr
    !---------------------------------------------------------------------------

    ATMOS_RD_F0 = F0
    ATMOS_RD_F1 = F1
    ATMOS_RD_Dval = Dval

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-RD]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ DYCOMS-II Parametarized Radiative heating'

    if ( RD_TYPE /= 'DYCOMSII' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_RD_TYPE is not DYCOMSII. Check!'
       call PRC_MPIstop
    endif


    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD_DYCOMS2,iostat=ierr)

    if( ierr < 0 ) then !--- missing
      if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_RD. Check!'
      call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_RD_DYCOMS2)

    F0 = ATMOS_RD_F0
    F1 = ATMOS_RD_F1
    Dval = ATMOS_RD_Dval

    return
  end subroutine ATMOS_PHY_RD_dycoms2_setup

  !-----------------------------------------------------------------------------
  ! Parametarized Radiative heating
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_RD_dycoms2( &
       flux_rad, flux_top,  & ! [out]
       solins, cosSZA,      & ! [out]
       DENS, RHOT, QTRC,    & ! [in]
       temp_sfc, param_sfc, & ! [in]
       CZ, FZ, CDZ, RCDZ,   & ! [in]
       REAL_lon, REAL_lat,  & ! [in]
       TIME_NOWDATE         ) ! [in]
    use mod_const, only: &
       Rdry    => CONST_Rdry,   &
       CPdry   => CONST_CPdry,  &
       RovCP   => CONST_RovCP,  &
       RovCV   => CONST_RovCV,  &
       CPovCV  => CONST_CPovCV, &
       EPSTvap => CONST_CPdry,  &
       P00     => CONST_PRE00,  &
       UNDEF   => CONST_UNDEF
    use mod_time, only: &
       dtrd => TIME_DTSEC_ATMOS_PHY_RD
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_history, only: &
       HIST_in
    use mod_atmos_phy_rd_common, only: &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up
    implicit none
    real(RP), intent(out) :: flux_rad(KA,IA,JA,2,2)
    real(RP), intent(out) :: flux_top(IA,JA,2)
    real(RP), intent(out) :: solins(IA,JA)
    real(RP), intent(out) :: cosSZA(IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: temp_sfc(IA,JA)
    real(RP), intent(in)  :: param_sfc(5)
    real(RP), intent(in)  :: CZ(KA)
    real(RP), intent(in)  :: FZ(KA-1)
    real(RP), intent(in)  :: CDZ(KA)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: REAL_lon(IA,JA)
    real(RP), intent(in)  :: REAL_lat(IA,JA)
    integer , intent(in)  :: TIME_NOWDATE(6)

    real(RP) :: EFLX_rad(KA,IA,JA) ! Radiative heating flux [J/m2/s]
    real(RP) :: Zi      (   IA,JA) ! Cloud top height [m]

    real(RP) :: QTOT ! Qv + Qc + Qr [kg/kg]

    real(RP) :: Qbelow, Qabove ! scaled LWP (above/below layer)

    real(RP), parameter :: kappa = 85.00_RP  ! scaling factor for LWP [m2/kg]
    real(RP), parameter :: a     =  1.00_RP  ! [K/m**-1/3]

    real(RP) :: pres

    real(RP) :: dQ, dZ, QWSUM
    integer :: k_cldtop

    integer :: k, k2, i, j, iq
    real(RP) :: TEMP_t_out  (KA,IA,JA) ! tendency rho*theta     [K/day]

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation (DYCOMS-II'

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          Zi      (  i,j) = 0.0_RP
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

       Zi(i,j) = CZ(k_cldtop)

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

          flux_rad(k,i,j,I_LW,I_up) = F0 * exp( -Qabove ) &
                                    + F1 * exp( -Qbelow )
       enddo

       do k = k_cldtop, KE
          dZ = FZ(k)-CZ(k_cldtop)
          QTOT = 0.0_RP
          do iq = QQS, QWE
             QTOT = QTOT + QTRC(k,i,j,iq)
          enddo
          flux_rad(k,i,j,I_LW,I_up) = flux_rad(k,i,j,I_LW,I_up) &
                          + a * DENS(k_cldtop,i,j)*( 1.0_RP-QTOT ) * CPdry * Dval &
                          * ( 0.250_RP * dZ  * dZ**(1.0_RP/3.0_RP) &
                              + CZ(k_cldtop) * dZ**(1.0_RP/3.0_RP) )
       enddo

     enddo
     enddo

     flux_rad(:,:,:,I_SW,I_dn) = 0.0_RP
     flux_rad(:,:,:,I_SW,I_up) = 0.0_RP
     flux_rad(:,:,:,I_LW,I_dn) = 0.0_RP
     flux_top(:,:,:) = UNDEF
     solins(:,:) = UNDEF
     cosSZA(:,:) = UNDEF

     if ( .not. first ) then
        call HIST_in( Zi(:,:), 'Zi', 'Cloud top height', 'm', dtrd )
        first = .false.
     endif


    return
  end subroutine ATMOS_PHY_RD_dycoms2

end module mod_atmos_phy_rd_dycoms2
