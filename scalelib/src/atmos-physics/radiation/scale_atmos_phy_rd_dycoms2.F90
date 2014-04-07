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
module scale_atmos_phy_rd_dycoms2
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
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: RD_TYPE

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
       DENS, RHOT, QTRC,      &
       CZ, FZ,                &
       oceanfrc,              &
       temp_sfc, albedo_land, &
       solins, cosSZA,        &
       flux_rad,              &
       flux_rad_top           )
    use scale_const, only: &
       Rdry    => CONST_Rdry,   &
       CPdry   => CONST_CPdry,  &
       EPSTvap => CONST_CPdry,  &
       P00     => CONST_PRE00,  &
       UNDEF   => CONST_UNDEF
    use scale_time, only: &
       dtrd => TIME_DTSEC_ATMOS_PHY_RD
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    use scale_history, only: &
       HIST_in
    use scale_atmos_phy_rd_common, only: &
       I_LW, &
       I_up
    implicit none

    real(RP), intent(in)  :: DENS        (KA,IA,JA)
    real(RP), intent(in)  :: RHOT        (KA,IA,JA)
    real(RP), intent(in)  :: QTRC        (KA,IA,JA,QA)
    real(RP), intent(in)  :: CZ          (  KA,IA,JA)
    real(RP), intent(in)  :: FZ          (0:KA,IA,JA)
    real(RP), intent(in)  :: oceanfrc    (IA,JA)      ! UNUSED
    real(RP), intent(in)  :: temp_sfc    (IA,JA)      ! UNUSED
    real(RP), intent(in)  :: albedo_land (IA,JA,2)    ! UNUSED
    real(RP), intent(in)  :: solins      (IA,JA)      ! UNUSED
    real(RP), intent(in)  :: cosSZA      (IA,JA)      ! UNUSED
    real(RP), intent(out) :: flux_rad    (KA,IA,JA,2,2)
    real(RP), intent(out) :: flux_rad_top(IA,JA,2)    ! UNDEFINED

    real(RP), parameter :: kappa = 85.00_RP  ! scaling factor for LWP [m2/kg]
    real(RP), parameter :: a     =  1.00_RP  ! [K/m**-1/3]

    real(RP) :: CDZ(KA)
    real(RP) :: Zi(IA,JA) ! Cloud top height [m]
    real(RP) :: dZ, dZ_CBRT
    integer  :: k_cldtop

    real(RP) :: QTOT ! Qv + Qc + Qr [kg/kg]
    real(RP) :: Qbelow, Qabove ! scaled LWP (above/below layer)
    real(RP) :: dQ, QWSUM

    integer :: k, k2, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation (DYCOMS-II'

    flux_rad    (:,:,:,:,:) = 0.0_RP
    flux_rad_top(:,:,:)     = UNDEF

    Zi(:,:) = 0.0_RP ! for history

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          CDZ(k) = FZ(k,i,j) - FZ(k-1,i,j)
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

       Zi(i,j) = CZ(k_cldtop,i,j)

       do k = KS-1, KE
          Qbelow = 0.0_RP
          Qabove = 0.0_RP

          do k2 = KS, KE
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
           dZ      = FZ(k,i,j)-CZ(k_cldtop,i,j)
           dZ_CBRT = dZ**(1.0_RP/3.0_RP)

           QTOT = 0.0_RP
           do iq = QQS, QWE
              QTOT = QTOT + QTRC(k,i,j,iq)
           enddo

           flux_rad(k,i,j,I_LW,I_up) = flux_rad(k,i,j,I_LW,I_up) &
                                     + a * DENS(k_cldtop,i,j)*( 1.0_RP-QTOT ) * CPdry * Dval &
                                     * ( 0.250_RP * dZ * dZ_CBRT + CZ(k_cldtop,i,j) * dZ_CBRT )
        enddo

     enddo
     enddo

     if ( .not. first ) then
        call HIST_in( Zi(:,:), 'Zi', 'Cloud top height', 'm', dtrd )
        first = .false.
     endif


    return
  end subroutine ATMOS_PHY_RD_dycoms2

end module scale_atmos_phy_rd_dycoms2
