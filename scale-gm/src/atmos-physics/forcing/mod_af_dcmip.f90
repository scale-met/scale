!-------------------------------------------------------------------------------
!> Module DCMIP2016 Physics forcing
!!
!! @par Description
!!          This module contains subroutines for forcing term of DCMIP2016 physics
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_af_dcmip
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index
  use scale_tracer

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: af_dcmip_init
  public :: af_dcmip

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private :: USE_Kessler         = .false.
  logical, private :: USE_SimpleMicrophys = .false.
  logical, private :: SM_Latdepend_SST    = .false. ! more option for SimpleMicrophysics
  logical, private :: SM_LargeScaleCond   = .false. ! more option for SimpleMicrophysics
  logical, private :: SM_PBL_Bryan        = .false. ! more option for SimpleMicrophysics
  logical, private :: USE_ToyChemistry    = .false.
  logical, public  :: USE_HeldSuarez      = .false.

  integer :: vlayer
  integer :: kdim

  integer :: I_QV
  integer :: I_QC
  integer :: I_QR

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_dcmip_init
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       HYDROMETEOR_regist => ATMOS_HYDROMETEOR_regist
    use mod_atmos_admin, only: &
       ATMOS_PHY_CH_TYPE
    use mod_chemvar, only: &
       NCHEM_MAX
    use mod_af_heldsuarez, only: &
       AF_heldsuarez_init
    use mod_bndcnd, only: &
       Tsfc => tem_sfc
    use mod_grd, only: &
       GRD_LAT,  &
       GRD_s
    use mod_gtl, only: &
       GTL_clip_region_1layer
    use scale_const, only:    &
       omega => CONST_OHM,    &
       pi    => CONST_PI,     &
       rair  => CONST_Rdry,   &
       a     => CONST_RADIUS, &
       rh2o  => CONST_Rvap
    implicit none

    logical :: SET_RJ2012          = .false.
    logical :: SET_DCMIP2016_11    = .false.
    logical :: SET_DCMIP2016_12    = .false.
    logical :: SET_DCMIP2016_13    = .false.
    logical :: SET_DCMIP2016_DRY   = .false.
    logical :: SET_DCMIP2016_LSC   = .false. ! large scale condensation
    logical :: SET_DCMIP2016_NOSST = .false.
    logical :: SET_DCMIP2016_21    = .false.

    namelist /FORCING_DCMIP_PARAM/ &
       SET_RJ2012,          &
       SET_DCMIP2016_11,    &
       SET_DCMIP2016_12,    &
       SET_DCMIP2016_13,    &
       SET_DCMIP2016_DRY,   &
       SET_DCMIP2016_LSC,   &
       SET_DCMIP2016_NOSST, &
       SET_DCMIP2016_21,    &
       USE_Kessler,         &
       USE_SimpleMicrophys, &
       SM_Latdepend_SST,    &
       SM_LargeScaleCond,   &
       SM_PBL_Bryan,        &
       USE_ToyChemistry,    &
       USE_HeldSuarez

    real(RP) :: dphi2
    real(RP) :: q0 = 0.021_RP           ! Maximum specific humidity for baro test
    real(RP) :: latw                    ! Halfwidth for  for baro test
    real(RP) :: etav                    ! Auxiliary variable for baro test
    real(RP) :: SST_TC   = 302.15_RP      ! Constant Value for SST
    real(RP) :: u0       = 35.0_RP          ! Zonal wind constant for moist baro test
    real(RP) :: eta0     = 0.252_RP         ! Center of jets (hybrid) for baro test
    real(RP) :: T00      = 288.0_RP         ! Horizontal mean T at surface for moist baro test
    real(RP) :: zvir  ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
    real(RP) :: lat    (ADM_gall,ADM_lall)

    integer :: ierr, ij, l
    !---------------------------------------------------------------------------

    !--- read parameters
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[af_dcmip]/Category[nhm forcing]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=FORCING_DCMIP_PARAM,iostat=ierr)
    if ( ierr < 0 ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** FORCING_DCMIP_PARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*,*) 'xxx Not appropriate names in namelist FORCING_DCMIP_PARAM. STOP.'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_NML,nml=FORCING_DCMIP_PARAM)

    ! overwrite setting
    if ( SET_RJ2012 ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of Reed and Jablonowski (2012)'
       USE_Kessler         = .false.
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .true.
       USE_ToyChemistry    = .false.

    elseif( SET_DCMIP2016_11 ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016 Case 1-1 (Moist baroclinic wave with terminator chemistry)'
       USE_Kessler         = .true.
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .false.
       USE_ToyChemistry    = .true.

       if ( SET_DCMIP2016_DRY ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: DRY condition'
          USE_Kessler         = .false.
       endif

       if ( SET_DCMIP2016_LSC ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: USE SM_LargeScaleCond'
          USE_Kessler         = .false.
          SM_LargeScaleCond   = .true.
       endif

       if ( SET_DCMIP2016_NOSST ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: Only Precipitation'
          USE_SimpleMicrophys = .false.
          SM_Latdepend_SST    = .false.
       endif

    elseif( SET_DCMIP2016_12 ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016 Case 1-2 (Idealized tropical cyclone)'
       USE_Kessler         = .true.
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .true.
       SM_LargeScaleCond   = .false.
       USE_ToyChemistry    = .false.

       if ( SET_DCMIP2016_DRY ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: DRY condition'
          USE_Kessler         = .false.
       endif

       if ( SET_DCMIP2016_LSC ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: USE SM_LargeScaleCond'
          USE_Kessler         = .false.
          SM_LargeScaleCond   = .true.
       endif

       if ( SET_DCMIP2016_NOSST ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: Only Precipitation'
          USE_SimpleMicrophys = .false.
          SM_Latdepend_SST    = .false.
       endif

    elseif( SET_DCMIP2016_13 ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016 Case 1-3 (Mesoscale storm)'
       USE_Kessler         = .true.
       USE_SimpleMicrophys = .false.
       SM_Latdepend_SST    = .false.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       USE_ToyChemistry    = .false.

       if ( SET_DCMIP2016_DRY ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: DRY condition'
          USE_Kessler         = .false.
       endif

       if ( SET_DCMIP2016_LSC ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: USE SM_LargeScaleCond'
          USE_Kessler         = .false.
          USE_SimpleMicrophys = .true.
          SM_LargeScaleCond   = .true.
       endif

    elseif( SET_DCMIP2016_21 ) then

       if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016 Case 2-1 (Moist Held-Suarez test)'
       USE_Kessler         = .true.
       USE_SimpleMicrophys = .true.
       SM_Latdepend_SST    = .false.
       SM_LargeScaleCond   = .false.
       SM_PBL_Bryan        = .false.
       USE_ToyChemistry    = .false.
       USE_HeldSuarez      = .true.

       if ( SET_DCMIP2016_LSC ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** Force setting of DCMIP2016: USE SM_LargeScaleCond'
          USE_Kessler         = .false.
          SM_LargeScaleCond   = .true.
       endif

       call AF_heldsuarez_init( moist_case = .true. )

    endif

    if ( USE_Kessler ) then
       call HYDROMETEOR_regist( 2, 0,                            & ! [IN]
                                (/ "QV", "QC", "QR" /),          & ! [IN]
                                (/ "water vapor mass ratio", "cloud water mass ratio", "rain water mass ratio " /), & ! [IN]
                                (/ "kg/kg", "kg/kg", "kg/kg" /), & ! [IN]
                                I_QV                             ) ! [OUT]
       I_QC = I_QV + 1
       I_QR = I_QV + 2
    else
       call HYDROMETEOR_regist( 0, 0,                           & ! [IN]
                                (/ "QV" /),                     & ! [IN]
                                (/ "water vapor mass ratio" /), & ! [IN]
                                (/ "kg/kg" /),                  & ! [IN]
                                I_QV                            ) ! [OUT]
    end if

    if( IO_L ) write(IO_FID_LOG,*) '*** Final Settings of FORCING_DCMIP_PARAM'
    if( IO_L ) write(IO_FID_LOG,*) '+ USE_Kessler         : ', USE_Kessler
    if( IO_L ) write(IO_FID_LOG,*) '+ USE_SimpleMicrophys : ', USE_SimpleMicrophys
    if( IO_L ) write(IO_FID_LOG,*) '+ SM_LargeScaleCond   : ', SM_LargeScaleCond
    if( IO_L ) write(IO_FID_LOG,*) '+ SM_Latdepend_SST    : ', SM_Latdepend_SST
    if( IO_L ) write(IO_FID_LOG,*) '+ SM_PBL_Bryan        : ', SM_PBL_Bryan
    if( IO_L ) write(IO_FID_LOG,*) '+ USE_ToyChemistry    : ', USE_ToyChemistry
    if( IO_L ) write(IO_FID_LOG,*) '+ USE_HeldSuarez      : ', USE_HeldSuarez

    ! initial value of the tracer is set in mod_prgvar - mod_ideal_init
    if ( USE_ToyChemistry ) then
       if ( ATMOS_PHY_CH_TYPE == 'PASSIVE' ) then
          if ( NCHEM_MAX /= 2 ) then
             write(*,*) 'xxx Not appropriate number of passive tracer. STOP.', NCHEM_MAX
             call PRC_abort
          endif
       else
          write(*,*) 'xxx ATMOS_PHY_CH_TYPE must be set to PASSIVE. STOP.', trim(ATMOS_PHY_CH_TYPE)
          call PRC_abort
       endif
    endif
!
! Set Sea Surface Temperature (constant for tropical cyclone)
! Tsurf needs to be dependent on latitude for moist baroclinic wave test
! Tsurf needs to be constant for tropical cyclone test
!
    lat(:,:) = GRD_s(:,ADM_KNONE,:,GRD_LAT)
    latw = 2.0_RP * pi / 9.0_RP
    etav = (1._RP-eta0) * 0.5_RP * pi
    zvir   = (rh2o/rair) - 1._RP

    do l=1, ADM_lall
       if ( USE_HeldSuarez ) then ! Moist H-S Test
          dphi2 = ( 26.D0 / 180.D0 * pi )**2
          do ij = 1, ADM_gall
             Tsfc(ij,l) = 29.D0 * exp( -0.5D0 * lat(ij,l) * lat(ij,l) / dphi2 ) + 271.D0
          enddo
       else
          if ( SM_Latdepend_SST ) then ! Moist Baroclinic Wave Test
             do ij=1, ADM_gall
                Tsfc(ij,l) = (T00 + pi*u0/rair * 1.5_RP * sin(etav) * (cos(etav))**0.5_RP *  &
        ((-2._RP*(sin(lat(ij,l)))**6 * ((cos(lat(ij,l)))**2 + 1._RP/3._RP) + 10._RP/63._RP) * &
                     u0 * (cos(etav))**1.5_RP  +  &
  (8._RP/5._RP*(cos(lat(ij,l)))**3 * ((sin(lat(ij,l)))**2 + 2._RP/3._RP) - pi/4._RP)*a*omega*0.5_RP ))/ &
                    (1._RP+zvir*q0*exp(-(lat(ij,l)/latw)**4))
             end do
          else ! Tropical Cyclone Test
             do ij=1, ADM_gall
                Tsfc(ij,l) = SST_TC
             end do
          end if
       endif
    enddo

    vlayer = ADM_vlayer
    kdim   = ADM_kall

    return
  end subroutine af_dcmip_init

  !-----------------------------------------------------------------------------
  subroutine af_dcmip( &
       ijdim,   &
       lat,     &
       lon,     &
       alt,     &
       alth,    &
       rho,     &
       pre,     &
       tem,     &
       vx,      &
       vy,      &
       vz,      &
       q,       &
       ein,     &
       pre_sfc, &
       tem_sfc, &
       fvx,     &
       fvy,     &
       fvz,     &
       fe,      &
       fq,      &
       precip,  &
       ix,      &
       iy,      &
       iz,      &
       jx,      &
       jy,      &
       jz,      &
       dt       )
    use scale_const, only: &
       d2r   => CONST_D2R,   &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
    use mod_chemvar, only: &
       NCHEM_STR, &
       NCHEM_END
    use mod_simple_physics, only: &
       simple_physics
    use mod_af_heldsuarez, only: &
       AF_heldsuarez
    use Terminator, only: &
       tendency_Terminator
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: lat    (ijdim)
    real(RP), intent(in)  :: lon    (ijdim)
    real(RP), intent(in)  :: alt    (ijdim,kdim)
    real(RP), intent(in)  :: alth   (ijdim,kdim)
    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: vx     (ijdim,kdim)
    real(RP), intent(in)  :: vy     (ijdim,kdim)
    real(RP), intent(in)  :: vz     (ijdim,kdim)
    real(RP), intent(in)  :: q      (ijdim,kdim,QA)
    real(RP), intent(in)  :: ein    (ijdim,kdim)
    real(RP), intent(in)  :: pre_sfc(ijdim)
    real(RP), intent(in)  :: tem_sfc(ijdim)
    real(RP), intent(out) :: fvx    (ijdim,kdim)
    real(RP), intent(out) :: fvy    (ijdim,kdim)
    real(RP), intent(out) :: fvz    (ijdim,kdim)
    real(RP), intent(out) :: fe     (ijdim,kdim)
    real(RP), intent(out) :: fq     (ijdim,kdim,QA)
    real(RP), intent(out) :: precip (ijdim)
    real(RP), intent(in)  :: ix     (ijdim)
    real(RP), intent(in)  :: iy     (ijdim)
    real(RP), intent(in)  :: iz     (ijdim)
    real(RP), intent(in)  :: jx     (ijdim)
    real(RP), intent(in)  :: jy     (ijdim)
    real(RP), intent(in)  :: jz     (ijdim)
    real(DP), intent(in)  :: dt

    ! for kessler
    real(RP) :: theta(vlayer) ! potential temperature (K)
    real(RP) :: qv   (vlayer) ! water vapor mixing ratio (gm/gm)
    real(RP) :: qc   (vlayer) ! cloud water mixing ratio (gm/gm)
    real(RP) :: qr   (vlayer) ! rain  water mixing ratio (gm/gm)
    real(RP) :: rhod (vlayer) ! dry air density (not mean state as in KW) (kg/m^3)
    real(RP) :: pk   (vlayer) ! Exner function (p/p0)**(R/cp)
    real(RP) :: z    (vlayer) ! heights of thermo. levels in the grid column (m)
    real(RP) :: qd   (vlayer)
    real(RP) :: cv   (vlayer)

    ! for simple physics
    real(DP) :: t      (ijdim,vlayer)   ! Temperature at full-model level (K)
    real(DP) :: qvv    (ijdim,vlayer)   ! Specific Humidity at full-model level (kg/kg)
    real(DP) :: u      (ijdim,vlayer)   ! Zonal wind at full-model level (m/s)
    real(DP) :: v      (ijdim,vlayer)   ! Meridional wind at full-model level (m/s)
    real(DP) :: pmid   (ijdim,vlayer)   ! Pressure is full-model level (Pa)
    real(DP) :: pint   (ijdim,vlayer+1) ! Pressure at model interfaces (Pa)
    real(DP) :: pdel   (ijdim,vlayer)   ! Layer thickness (Pa)
    real(DP) :: rpdel  (ijdim,vlayer)   ! Reciprocal of layer thickness (1/Pa)
    real(DP) :: ps     (ijdim)          ! surface pressure output [dummy]
    real(DP) :: ts     (ijdim)
    real(DP) :: lat_rad(ijdim)
    real(DP) :: precip2(ijdim)

    ! for toy-chemistory
    real(DP) :: lat_deg, lon_deg
    real(DP) :: cl, cl2
    real(DP) :: cl_f, cl2_f
    real(DP) :: qvd

    ! for Held-Suarez
    real(RP) :: gvx(ijdim,kdim)
    real(RP) :: gvy(ijdim,kdim)
    real(RP) :: gvz(ijdim,kdim)
    real(RP) :: ge (ijdim,kdim)

    integer :: kmin, kmax

    integer  :: ij, k, kk
    !---------------------------------------------------------------------------

    kmin = ADM_kmin
    kmax = ADM_kmax

    fvx(:,:)   = 0.0_RP
    fvy(:,:)   = 0.0_RP
    fvz(:,:)   = 0.0_RP
    fe (:,:)   = 0.0_RP
    fq (:,:,:) = 0.0_RP

    precip (:) = 0.0_RP
    precip2(:) = 0.0_RP

    if ( USE_Kessler ) then
       do ij = 1, ijdim

          qd   (:) = 1.0_RP               &
                   - q(ij,kmin:kmax,I_QV) &
                   - q(ij,kmin:kmax,I_QC) &
                   - q(ij,kmin:kmax,I_QR)

          qv   (:) = q  (ij,kmin:kmax,I_QV) / qd(:)
          qc   (:) = q  (ij,kmin:kmax,I_QC) / qd(:)
          qr   (:) = q  (ij,kmin:kmax,I_QR) / qd(:)
          rhod (:) = rho(ij,kmin:kmax)      * qd(:)

          pk   (:) = ( pre(ij,kmin:kmax) / PRE00 )**( Rdry / CPdry )
          theta(:) = tem(ij,kmin:kmax) / pk(:)
          z    (:) = alt(ij,kmin:kmax)

          call kessler( theta(:),  & ! [INOUT]
                        qv   (:),  & ! [INOUT]
                        qc   (:),  & ! [INOUT]
                        qr   (:),  & ! [INOUT]
                        rhod (:),  & ! [IN]
                        pk   (:),  & ! [IN]
                        dt,        & ! [IN]
                        z    (:),  & ! [IN]
                        vlayer,    & ! [IN]
                        precip(ij) ) ! [OUT]

          qd(:) = 1.0_RP &
                / ( 1.0_RP + qv(:) + qc(:) + qr(:) )
          qv(:) = qv(:) * qd(:)
          qc(:) = qc(:) * qd(:)
          qr(:) = qr(:) * qd(:)

          cv(:) = qd(:) * CVdry     &
                + qv(:) * TRACER_CV(I_QV) &
                + qc(:) * TRACER_CV(I_QC) &
                + qr(:) * TRACER_CV(I_QR)

          fq(ij,kmin:kmax,I_QV) = fq(ij,kmin:kmax,I_QV) + ( qv(:) - q(ij,kmin:kmax,I_QV) ) / dt
          fq(ij,kmin:kmax,I_QC) = fq(ij,kmin:kmax,I_QC) + ( qc(:) - q(ij,kmin:kmax,I_QC) ) / dt
          fq(ij,kmin:kmax,I_QR) = fq(ij,kmin:kmax,I_QR) + ( qr(:) - q(ij,kmin:kmax,I_QR) ) / dt
          fe(ij,kmin:kmax)      = fe(ij,kmin:kmax)      + ( cv(:)*theta(:)*pk(:) - ein(ij,kmin:kmax) ) / dt
       enddo
    endif

    if ( USE_SimpleMicrophys ) then

       do k = 1, vlayer
          kk = kmax - k + 1 ! reverse order

          t   (:,k) = tem(:,kk)
          qvv (:,k) = q  (:,kk,I_QV)
          u   (:,k) = vx (:,kk) * ix(:) &
                    + vy (:,kk) * iy(:) &
                    + vz (:,kk) * iz(:)
          v   (:,k) = vx (:,kk) * jx(:) &
                    + vy (:,kk) * jy(:) &
                    + vz (:,kk) * jz(:)
          pmid(:,k) = pre(:,kk)
       enddo

       pint(:,1) = 0.0_RP
       do k = 2, vlayer
          kk = kmax - k + 1 ! reverse order
          pint(:,k) = pre(:,kk) * exp( log( pre(:,kk+1) / pre(:,kk) )                     &
                                     * (alth(:,kk+1)-alt(:,kk)) / (alt(:,kk+1)-alt(:,kk)) )
       enddo
       pint(:,vlayer+1) = pre_sfc(:)

       do k = 1, vlayer
          pdel (:,k) = pint(:,k+1) - pint(:,k)
          rpdel(:,k) = 1.0_RP / pdel(:,k)
       enddo

       ps     (:) = pre_sfc(:)
       ts     (:) = real( tem_sfc(:), kind=DP )
       lat_rad(:) = real( lat    (:), kind=DP )

       call simple_physics( ijdim,             & ! [IN]
                            vlayer,            & ! [IN]
                            dt,                & ! [IN]
                            lat_rad(:),        & ! [IN]
                            ts    (:),        & ! [IN]
                            t     (:,:),       & ! [INOUT]
                            qvv   (:,:),       & ! [INOUT]
                            u     (:,:),       & ! [INOUT]
                            v     (:,:),       & ! [INOUT]
                            pmid  (:,:),       & ! [INOUT] but not changed
                            pint  (:,:),       & ! [INOUT] but not changed
                            pdel  (:,:),       & ! [INOUT] but not changed
                            rpdel (:,:),       & ! [INOUT] but not changed
                            ps    (:),         & ! [INOUT] but not changed
                            precip2(:),        & ! [OUT]
                            SM_LargeScaleCond, & ! [IN]
                            SM_PBL_Bryan,      & ! [IN]
                            USE_HeldSuarez     ) ! [IN]

       do k = 1, vlayer
          kk = kmax - k + 1 ! reverse order

          fvx(:,kk) = fvx(:,kk) + ( u(:,k) * ix(:) + v(:,k) * jx(:) - vx(:,kk) ) / dt
          fvy(:,kk) = fvy(:,kk) + ( u(:,k) * iy(:) + v(:,k) * jy(:) - vy(:,kk) ) / dt
          fvz(:,kk) = fvz(:,kk) + ( u(:,k) * iz(:) + v(:,k) * jz(:) - vz(:,kk) ) / dt
       enddo

       do k = 1, vlayer
          kk = kmax - k + 1 ! reverse order

          do ij = 1, ijdim
             if ( USE_Kessler ) then
                qv(k) = qvv(ij,k)
                qc(k) = q  (ij,kk,I_QC)
                qr(k) = q  (ij,kk,I_QR)

                qd(k) = 1.0_RP &
                      - qv(k)  &
                      - qc(k)  &
                      - qr(k)

                cv(k) = qd(k) * CVdry     &
                      + qv(k) * TRACER_CV(I_QV) &
                      + qc(k) * TRACER_CV(I_QC) &
                      + qr(k) * TRACER_CV(I_QR)
             else
                qv(k) = qvv(ij,k)

                qd(k) = 1.0_RP &
                      - qv(k)

                cv(k) = qd(k) * CVdry     &
                      + qv(k) * TRACER_CV(I_QV)
             endif

             fq(ij,kk,I_QV) = fq(ij,kk,I_QV) + ( qv(k) - q(ij,kk,I_QV) ) / dt
             fe(ij,kk)      = fe(ij,kk)      + ( cv(k) * t(ij,k) - ein(ij,kk) ) / dt
          enddo
       enddo

       precip(:) = precip(:) + precip2(:)

    endif

    if ( USE_ToyChemistry ) then
       do k  = kmin, kmax
       do ij = 1,    ijdim
          lat_deg = real( lat(ij) / d2r          , kind=DP )
          lon_deg = real( lon(ij) / d2r          , kind=DP )
          qvd     = real( 1.0_RP - q(ij,k,I_QV)  , kind=DP )
          cl      = real( q(ij,k,NCHEM_STR) / qvd, kind=DP )
          cl2     = real( q(ij,k,NCHEM_END) / qvd, kind=DP )

          call tendency_Terminator( lat_deg, & ! [IN]
                                    lon_deg, & ! [IN]
                                    cl,      & ! [IN]
                                    cl2,     & ! [IN]
                                    dt,      & ! [IN]
                                    cl_f,    & ! [OUT]
                                    cl2_f    ) ! [OUT]

          fq(ij,k,NCHEM_STR) = fq(ij,k,NCHEM_STR) + real( cl_f  * qvd, kind=RP )
          fq(ij,k,NCHEM_END) = fq(ij,k,NCHEM_END) + real( cl2_f * qvd, kind=RP )
       enddo
       enddo
    endif

    if ( USE_HeldSuarez ) then
       call AF_heldsuarez( ijdim,    & ! [IN]
                           lat(:),   & ! [IN]
                           pre(:,:), & ! [IN]
                           tem(:,:), & ! [IN]
                           vx (:,:), & ! [IN]
                           vy (:,:), & ! [IN]
                           vz (:,:), & ! [IN]
                           gvx(:,:), & ! [OUT]
                           gvy(:,:), & ! [OUT]
                           gvz(:,:), & ! [OUT]
                           ge (:,:)  ) ! [OUT]

       ! overwrite tendency from simple_physics
       fvx(:,:) = gvx(:,:)
       fvy(:,:) = gvy(:,:)
       fvz(:,:) = gvz(:,:)

       fe (:,:) = fe (:,:) + ge (:,:)
    endif

    return
  end subroutine af_dcmip

end module mod_af_dcmip
