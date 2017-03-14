!-------------------------------------------------------------------------------
!> module LAND / Surface fluxes with thin-slab land model
!!
!! @par Description
!!          Surface flux with thin-slab land model
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
module scale_land_sfc_thin_slab
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_SFC_THIN_SLAB_setup
  public :: LAND_SFC_THIN_SLAB

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
  integer,  private :: LAND_SFC_THIN_SLAB_itr_max = 100 ! maximum iteration number

  real(RP), private :: LAND_SFC_THIN_SLAB_dTS_max = 5.0E-2_RP ! maximum delta surface temperature [K/s]
  real(RP), private :: LAND_SFC_THIN_SLAB_res_min = 1.0E+0_RP ! minimum value of residual
  real(RP), private :: LAND_SFC_THIN_SLAB_err_min = 1.0E-2_RP ! minimum value of error
  real(RP), private :: LAND_SFC_THIN_SLAB_dreslim = 1.0E+2_RP ! limiter of d(residual)

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_SFC_THIN_SLAB_setup( LAND_TYPE )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: LAND_TYPE

    NAMELIST / PARAM_LAND_SFC_THIN_SLAB / &
       LAND_SFC_THIN_SLAB_itr_max, &
       LAND_SFC_THIN_SLAB_dTS_max, &
       LAND_SFC_THIN_SLAB_res_min, &
       LAND_SFC_THIN_SLAB_err_min, &
       LAND_SFC_THIN_SLAB_dreslim

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[THIN-SLAB] / Categ[LAND SFC] / Origin[SCALElib]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_SFC_THIN_SLAB,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_SFC_THIN_SLAB. Check!'
       call PRC_MPIstop
    endif
    if( IO_NML ) write(IO_FID_NML,nml=PARAM_LAND_SFC_THIN_SLAB)

    return
  end subroutine LAND_SFC_THIN_SLAB_setup

  !-----------------------------------------------------------------------------
  subroutine LAND_SFC_THIN_SLAB( &
        LST_t,      &
        ZMFLX,      &
        XMFLX,      &
        YMFLX,      &
        SHFLX,      &
        LHFLX,      &
        GHFLX,      &
        U10,        &
        V10,        &
        T2,         &
        Q2,         &
        TMPA,       &
        PRSA,       &
        WA,         &
        UA,         &
        VA,         &
        RHOA,       &
        QVA,        &
        Z1,         &
        PBL,        &
        PRSS,       &
        LWD,        &
        SWD,        &
        TG,         &
        LST,        &
        QVEF,       &
        ALB_LW,     &
        ALB_SW,     &
        DZG,        &
        TCS,        &
        Z0M,        &
        Z0H,        &
        Z0E,        &
        dt          )
    use scale_process, only: &
      PRC_myrank,  &
      PRC_MPIstop
    use scale_const, only: &
      Rdry  => CONST_Rdry,  &
      CPdry => CONST_CPdry, &
      STB   => CONST_STB
    use scale_landuse, only: &
      LANDUSE_fact_land
    use scale_atmos_hydrometeor, only: &
      HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
    use scale_atmos_saturation, only: &
      qsat => ATMOS_SATURATION_pres2qsat_all
    use scale_bulkflux, only: &
      BULKFLUX
    implicit none

    ! parameters
    real(RP), parameter :: dTS0     = 1.0E-4_RP ! delta surface temp.

    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0E+0_RP ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5E+0_RP ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1E+0_RP ! factor b in Tomita (2009)

    ! arguments
    real(RP), intent(out) :: LST_t(IA,JA) ! tendency of LST
    real(RP), intent(out) :: ZMFLX(IA,JA) ! z-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: XMFLX(IA,JA) ! x-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: YMFLX(IA,JA) ! y-momentum flux at the surface [kg/m2/s]
    real(RP), intent(out) :: SHFLX(IA,JA) ! sensible heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: LHFLX(IA,JA) ! latent heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: GHFLX(IA,JA) ! ground heat flux at the surface [J/m2/s]
    real(RP), intent(out) :: U10  (IA,JA) ! velocity u at 10m [m/s]
    real(RP), intent(out) :: V10  (IA,JA) ! velocity v at 10m [m/s]
    real(RP), intent(out) :: T2   (IA,JA) ! temperature at 2m [K]
    real(RP), intent(out) :: Q2   (IA,JA) ! water vapor at 2m [kg/kg]

    real(RP), intent(in) :: TMPA(IA,JA) ! temperature at the lowest atmospheric layer [K]
    real(RP), intent(in) :: PRSA(IA,JA) ! pressure at the lowest atmospheric layer [Pa]
    real(RP), intent(in) :: WA  (IA,JA) ! velocity w at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: UA  (IA,JA) ! velocity u at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: VA  (IA,JA) ! velocity v at the lowest atmospheric layer [m/s]
    real(RP), intent(in) :: RHOA(IA,JA) ! density at the lowest atmospheric layer [kg/m3]
    real(RP), intent(in) :: QVA (IA,JA) ! ratio of water vapor mass to total mass at the lowest atmospheric layer [kg/kg]
    real(RP), intent(in) :: Z1  (IA,JA) ! cell center height at the lowest atmospheric layer [m]
    real(RP), intent(in) :: PBL (IA,JA) ! the top of atmospheric mixing layer [m]
    real(RP), intent(in) :: PRSS(IA,JA) ! pressure at the surface [Pa]
    real(RP), intent(in) :: LWD (IA,JA) ! downward long-wave radiation flux at the surface [J/m2/s]
    real(RP), intent(in) :: SWD (IA,JA) ! downward short-wave radiation flux at the surface [J/m2/s]

    real(RP), intent(in) :: TG    (IA,JA) ! soil temperature [K]
    real(RP), intent(in) :: LST   (IA,JA) ! land surface temperature [K]
    real(RP), intent(in) :: QVEF  (IA,JA) ! efficiency of evaporation [0-1]
    real(RP), intent(in) :: ALB_LW(IA,JA) ! surface albedo for LW [0-1]
    real(RP), intent(in) :: ALB_SW(IA,JA) ! surface albedo for SW [0-1]
    real(RP), intent(in) :: DZG   (IA,JA) ! soil depth [m]
    real(RP), intent(in) :: TCS   (IA,JA) ! thermal conductivity for soil [J/m/K/s]
    real(RP), intent(in) :: Z0M   (IA,JA) ! roughness length for momemtum [m]
    real(RP), intent(in) :: Z0H   (IA,JA) ! roughness length for heat [m]
    real(RP), intent(in) :: Z0E   (IA,JA) ! roughness length for vapor [m]
    real(DP), intent(in) :: dt            ! delta time

    ! works
    real(RP) :: LST1(IA,JA)

    real(RP) :: res    ! residual
    real(RP) :: dres   ! d(residual)/dLST
    real(RP) :: oldres ! residual in previous step
    real(RP) :: redf   ! reduced factor

    real(RP) :: Ustar, Ustar10, Ustar2, dUstar ! friction velocity [m]
    real(RP) :: Tstar, Tstar10, Tstar2, dTstar ! friction potential temperature [K]
    real(RP) :: Qstar, Qstar10, Qstar2, dQstar ! friction water vapor mass ratio [kg/kg]
    real(RP) :: Uabs,  Uabs10,  Uabs2,  dUabs  ! modified absolute velocity [m/s]
    real(RP) :: QVsat, dQVsat ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: QVS, dQVS     ! water vapor mixing ratio at surface [kg/kg]

    real(RP) :: LHV(IA,JA)    ! latent heat of vaporization [J/kg]

    integer  :: i, j, n
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land surface step: Thin-Slab'

    call HYDROMETEOR_LHV( LHV(:,:), TMPA(:,:) )

    ! copy land surfce temperature for iteration
    do j = JS, JE
    do i = IS, IE
      LST1(i,j) = LST(i,j)
    end do
    end do

    ! update surface temperature
    do j = JS, JE
    do i = IS, IE

      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then

        redf   = 1.0_RP
        oldres = huge(0.0_RP)

        ! modified Newton-Raphson method (Tomita 2009)
        do n = 1, LAND_SFC_THIN_SLAB_itr_max

          call qsat( QVsat,     & ! [OUT]
                     LST1(i,j), & ! [IN]
                     PRSS(i,j)  ) ! [IN]
          call qsat( dQVsat,         & ! [OUT]
                     LST1(i,j)+dTS0, & ! [IN]
                     PRSS(i,j)       ) ! [IN]

          QVS  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat
          dQVS = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * dQVsat

          call BULKFLUX( &
              Ustar,     & ! [OUT]
              Tstar,     & ! [OUT]
              Qstar,     & ! [OUT]
              Uabs,      & ! [OUT]
              TMPA(i,j), & ! [IN]
              LST1(i,j), & ! [IN]
              PRSA(i,j), & ! [IN]
              PRSS(i,j), & ! [IN]
              QVA (i,j), & ! [IN]
              QVS,       & ! [IN]
              UA  (i,j), & ! [IN]
              VA  (i,j), & ! [IN]
              Z1  (i,j), & ! [IN]
              PBL (i,j), & ! [IN]
              Z0M (i,j), & ! [IN]
              Z0H (i,j), & ! [IN]
              Z0E (i,j)  ) ! [IN]

          call BULKFLUX( &
              dUstar,         & ! [OUT]
              dTstar,         & ! [OUT]
              dQstar,         & ! [OUT]
              dUabs,          & ! [OUT]
              TMPA(i,j),      & ! [IN]
              LST1(i,j)+dTS0, & ! [IN]
              PRSA(i,j),      & ! [IN]
              PRSS(i,j),      & ! [IN]
              QVA (i,j),      & ! [IN]
              dQVS,           & ! [IN]
              UA  (i,j),      & ! [IN]
              VA  (i,j),      & ! [IN]
              Z1  (i,j),      & ! [IN]
              PBL (i,j),      & ! [IN]
              Z0M (i,j),      & ! [IN]
              Z0H (i,j),      & ! [IN]
              Z0E (i,j)       ) ! [IN]

          ! calculation for residual
          res = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
              + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST1(i,j)**4 ) &
              + CPdry    * RHOA(i,j) * Ustar * Tstar &
              + LHV(i,j) * RHOA(i,j) * Ustar * Qstar &
              - 2.0_RP * TCS(i,j) * ( LST1(i,j) - TG(i,j) ) / DZG(i,j)

          ! calculation for d(residual)/dLST
          dres = -4.0_RP * ( 1.0_RP - ALB_LW(i,j) ) * STB * LST1(i,j)**3 &
               + CPdry    * RHOA(i,j) * ( (dUstar-Ustar)/dTS0 * Tstar + Ustar * (dTstar-Tstar)/dTS0 ) &
               + LHV(i,j) * RHOA(i,j) * ( (dUstar-Ustar)/dTS0 * Qstar + Ustar * (dQstar-Qstar)/dTS0 ) &
               - 2.0_RP * TCS(i,j) / DZG(i,j)

          ! convergence test with residual and error levels
          if( abs( res      ) < LAND_SFC_THIN_SLAB_res_min .or. &
              abs( res/dres ) < LAND_SFC_THIN_SLAB_err_min      ) then
            exit
          end if

          ! stop iteration to prevent numerical error
          if( abs(dres) * LAND_SFC_THIN_SLAB_dreslim < abs(res) ) then
            exit
          end if

          ! calculate reduced factor
          if( dres < 0.0_RP ) then
            if( abs(res) > abs(oldres) ) then
              redf = max( TFa*abs(redf), redf_min )
            else
              redf = min( TFb*abs(redf), redf_max )
            end if
          else
            redf = -1.0_RP
          end if

          ! estimate next surface temperature
          LST1(i,j) = LST1(i,j) - redf * res / dres

          ! save residual in this step
          oldres = res

        end do

        ! update land surface temperature with limitation
        LST1(i,j) = min( max( LST1(i,j), &
                              LST (i,j) - LAND_SFC_THIN_SLAB_dTS_max * dt ), &
                              LST (i,j) + LAND_SFC_THIN_SLAB_dTS_max * dt )

        if( n > LAND_SFC_THIN_SLAB_itr_max ) then
          ! land surface temperature was not converged
          if( IO_L ) write(IO_FID_LOG,'(A)'       ) 'Warning: land surface tempearture was not converged.'
          if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
          if( IO_L ) write(IO_FID_LOG,'(A,I32)'   ) 'DEBUG --- PRC_myrank                         [no unit] :', PRC_myrank
          if( IO_L ) write(IO_FID_LOG,'(A,I32)'   ) 'DEBUG --- number of i                        [no unit] :', i
          if( IO_L ) write(IO_FID_LOG,'(A,I32)'   ) 'DEBUG --- number of j                        [no unit] :', j
          if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
          if( IO_L ) write(IO_FID_LOG,'(A,I32)'   ) 'DEBUG --- loop number                        [no unit] :', n
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- Residual                           [J/m2/s]  :', res
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- delta Residual                     [J/m2/s]  :', dres
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') ''
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- temperature                        [K]       :', TMPA  (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- pressure                           [Pa]      :', PRSA  (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- velocity w                         [m/s]     :', WA    (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- velocity u                         [m/s]     :', UA    (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- velocity v                         [m/s]     :', VA    (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- density                            [kg/m3]   :', RHOA  (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- water vapor mass ratio             [kg/kg]   :', QVA   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- cell center height                 [m]       :', Z1    (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- atmospheric mixing layer height    [m]       :', PBL   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- pressure at the surface            [Pa]      :', PRSS  (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- downward long-wave radiation       [J/m2/s]  :', LWD   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- downward short-wave radiation      [J/m2/s]  :', SWD   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- soil temperature                   [K]       :', TG    (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- land surface temperature           [K]       :', LST   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- efficiency of evaporation          [0-1]     :', QVEF  (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- surface albedo for LW              [0-1]     :', ALB_LW(i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- surface albedo for SW              [0-1]     :', ALB_SW(i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- soil depth                         [m]       :', DZG   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- thermal conductivity for soil      [J/m/K/s] :', TCS   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- roughness length for momemtum      [m]       :', Z0M   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- roughness length for heat          [m]       :', Z0H   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- roughness length for vapor         [m]       :', Z0E   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A)'       ) ''
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- latent heat                        [J/kg]    :', LHV   (i,j)
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- friction velocity                  [m]       :', Ustar
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- friction potential temperature     [K]       :', Tstar
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- friction water vapor mass ratio    [kg/kg]   :', Qstar
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- d(friction velocity)               [m]       :', dUstar
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- d(friction potential temperature)  [K]       :', dTstar
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- d(friction water vapor mass ratio) [kg/kg]   :', dQstar
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- modified absolute velocity         [m/s]     :', Uabs
          if( IO_L ) write(IO_FID_LOG,'(A,F32.16)') 'DEBUG --- next land surface temperature      [K]       :', LST1  (i,j)

          ! check NaN
          if ( .NOT. ( res > -1.0_RP .OR. res < 1.0_RP ) ) then ! must be NaN
             write(*,*) 'xxx NaN is detected for land surface temperature.'
             write(*,*) ''
             write(*,*) 'DEBUG --- PRC_myrank                                   :', PRC_myrank
             write(*,*) 'DEBUG --- number of i                                  :', i
             write(*,*) 'DEBUG --- number of j                                  :', j
             write(*,*) ''
             write(*,*) 'DEBUG --- Residual                           [J/m2/s]  :', res
             write(*,*) 'DEBUG --- delta Residual                     [J/m2/s]  :', dres
             write(*,*) ''
             write(*,*) 'DEBUG --- temperature                        [K]       :', TMPA  (i,j)
             write(*,*) 'DEBUG --- pressure                           [Pa]      :', PRSA  (i,j)
             write(*,*) 'DEBUG --- velocity w                         [m/s]     :', WA    (i,j)
             write(*,*) 'DEBUG --- velocity u                         [m/s]     :', UA    (i,j)
             write(*,*) 'DEBUG --- velocity v                         [m/s]     :', VA    (i,j)
             write(*,*) 'DEBUG --- density                            [kg/m3]   :', RHOA  (i,j)
             write(*,*) 'DEBUG --- water vapor mass ratio             [kg/kg]   :', QVA   (i,j)
             write(*,*) 'DEBUG --- cell center height                 [m]       :', Z1    (i,j)
             write(*,*) 'DEBUG --- atmospheric mixing layer height    [m]       :', PBL   (i,j)
             write(*,*) 'DEBUG --- pressure at the surface            [Pa]      :', PRSS  (i,j)
             write(*,*) 'DEBUG --- downward long-wave radiation       [J/m2/s]  :', LWD   (i,j)
             write(*,*) 'DEBUG --- downward short-wave radiation      [J/m2/s]  :', SWD   (i,j)
             write(*,*) ''
             write(*,*) 'DEBUG --- soil temperature                   [K]       :', TG    (i,j)
             write(*,*) 'DEBUG --- land surface temperature           [K]       :', LST   (i,j)
             write(*,*) 'DEBUG --- efficiency of evaporation          [0-1]     :', QVEF  (i,j)
             write(*,*) 'DEBUG --- surface albedo for LW              [0-1]     :', ALB_LW(i,j)
             write(*,*) 'DEBUG --- surface albedo for SW              [0-1]     :', ALB_SW(i,j)
             write(*,*) 'DEBUG --- soil depth                         [m]       :', DZG   (i,j)
             write(*,*) 'DEBUG --- thermal conductivity for soil      [J/m/K/s] :', TCS   (i,j)
             write(*,*) 'DEBUG --- roughness length for momemtum      [m]       :', Z0M   (i,j)
             write(*,*) 'DEBUG --- roughness length for heat          [m]       :', Z0H   (i,j)
             write(*,*) 'DEBUG --- roughness length for vapor         [m]       :', Z0E   (i,j)
             write(*,*) ''
             write(*,*) 'DEBUG --- latent heat                        [J/kg]    :', LHV   (i,j)
             write(*,*) 'DEBUG --- friction velocity                  [m]       :', Ustar
             write(*,*) 'DEBUG --- friction potential temperature     [K]       :', Tstar
             write(*,*) 'DEBUG --- friction water vapor mass ratio    [kg/kg]   :', Qstar
             write(*,*) 'DEBUG --- d(friction velocity)               [m]       :', dUstar
             write(*,*) 'DEBUG --- d(friction potential temperature)  [K]       :', dTstar
             write(*,*) 'DEBUG --- d(friction water vapor mass ratio) [kg/kg]   :', dQstar
             write(*,*) 'DEBUG --- modified absolute velocity         [m/s]     :', Uabs
             write(*,*) 'DEBUG --- next land surface temperature      [K]       :', LST1  (i,j)

             call PRC_MPIstop
          endif

        end if

      end if

      ! calculate tendency
      LST_t(i,j) = ( LST1(i,j) - LST(i,j) ) / dt

    end do
    end do

    ! calculate surface flux
    do j = JS, JE
    do i = IS, IE

      if( LANDUSE_fact_land(i,j) > 0.0_RP ) then

        call qsat( QVsat,     & ! [OUT]
                   LST1(i,j), & ! [IN]
                   PRSS(i,j)  ) ! [IN]

        QVS  = ( 1.0_RP - QVEF(i,j) ) * QVA(i,j) + QVEF(i,j) * QVsat

        call BULKFLUX( &
            Ustar,     & ! [OUT]
            Tstar,     & ! [OUT]
            Qstar,     & ! [OUT]
            Uabs,      & ! [OUT]
            TMPA(i,j), & ! [IN]
            LST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVS,       & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            Z1  (i,j), & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ! for 10m wind
        call BULKFLUX( &
            Ustar10,   & ! [OUT]
            Tstar10,   & ! [OUT]
            Qstar10,   & ! [OUT]
            Uabs10,    & ! [OUT]
            TMPA(i,j), & ! [IN]
            LST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVS,       & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            10.0_RP,   & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ! for 2m temperature / mixing ratio
        call BULKFLUX( &
            Ustar2,    & ! [OUT]
            Tstar2,    & ! [OUT]
            Qstar2,    & ! [OUT]
            Uabs2,     & ! [OUT]
            TMPA(i,j), & ! [IN]
            LST1(i,j), & ! [IN]
            PRSA(i,j), & ! [IN]
            PRSS(i,j), & ! [IN]
            QVA (i,j), & ! [IN]
            QVS,       & ! [IN]
            UA  (i,j), & ! [IN]
            VA  (i,j), & ! [IN]
            2.0_RP,    & ! [IN]
            PBL (i,j), & ! [IN]
            Z0M (i,j), & ! [IN]
            Z0H (i,j), & ! [IN]
            Z0E (i,j)  ) ! [IN]

        ZMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * WA(i,j)
        XMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * UA(i,j)
        YMFLX(i,j) = -RHOA(i,j) * Ustar**2 / Uabs * VA(i,j)
        SHFLX(i,j) = -CPdry    * RHOA(i,j) * Ustar * Tstar
        LHFLX(i,j) = -LHV(i,j) * RHOA(i,j) * Ustar * Qstar
        GHFLX(i,j) = -2.0_RP * TCS(i,j) * ( LST1(i,j) - TG(i,j) ) / DZG(i,j)

        ! calculation for residual
        res = ( 1.0_RP - ALB_SW(i,j) ) * SWD(i,j) &
            + ( 1.0_RP - ALB_LW(i,j) ) * ( LWD(i,j) - STB * LST1(i,j)**4 ) &
            - SHFLX(i,j) - LHFLX(i,j) + GHFLX(i,j)

        ! put residual in ground heat flux
        GHFLX(i,j) = GHFLX(i,j) - res

        ! diagnostic variables
        U10(i,j) = Ustar / Ustar10 * UA(i,j)
        V10(i,j) = Ustar / Ustar10 * VA(i,j)
        T2 (i,j) = ( 1.0_RP - Tstar / Tstar2 ) * LST1(i,j) + Tstar / Tstar2 * TMPA(i,j)
        Q2 (i,j) = ( 1.0_RP - Qstar / Qstar2 ) * QVS       + Qstar / Qstar2 * QVA (i,j)

      else

        ! not calculate surface flux
        ZMFLX(i,j) = 0.0_RP
        XMFLX(i,j) = 0.0_RP
        YMFLX(i,j) = 0.0_RP
        SHFLX(i,j) = 0.0_RP
        LHFLX(i,j) = 0.0_RP
        GHFLX(i,j) = 0.0_RP
        U10  (i,j) = 0.0_RP
        V10  (i,j) = 0.0_RP
        T2   (i,j) = 0.0_RP
        Q2   (i,j) = 0.0_RP

      end if

    end do
    end do

    return
  end subroutine LAND_SFC_THIN_SLAB

end module scale_land_sfc_thin_slab
