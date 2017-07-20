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
  use scale_stdio
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
  logical, private :: USE_HeldSuarez      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_dcmip_init
    use scale_process, only: &
       PRC_MPIstop
    use mod_runconf, only: &
       CHEM_TYPE, &
       NCHEM_MAX
    use mod_af_heldsuarez, only: &
       AF_heldsuarez_init
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

    integer :: ierr
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
       call PRC_MPIstop
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
       if ( CHEM_TYPE == 'PASSIVE' ) then
          if ( NCHEM_MAX /= 2 ) then
             write(*,*) 'xxx Not appropriate number of passive tracer. STOP.', NCHEM_MAX
             call PRC_MPIstop
          endif
       else
          write(*,*) 'xxx CHEM_TYPE must be set to PASSIVE. STOP.', trim(CHEM_TYPE)
          call PRC_MPIstop
       endif
    endif

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
    use mod_adm, only: &
       vlayer => ADM_vlayer, &
       kdim   => ADM_kall,   &
       kmin   => ADM_kmin,   &
       kmax   => ADM_kmax
    use scale_const, only: &
       d2r   => CONST_D2R,   &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
    use mod_runconf, only: &
       TRC_VMAX,  &
       RAIN_TYPE, &
       I_QV,      &
       I_QC,      &
       I_QR,      &
       NCHEM_STR, &
       NCHEM_END, &
       CVW
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
    real(RP), intent(in)  :: q      (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: ein    (ijdim,kdim)
    real(RP), intent(in)  :: pre_sfc(ijdim)
    real(RP), intent(out) :: fvx    (ijdim,kdim)
    real(RP), intent(out) :: fvy    (ijdim,kdim)
    real(RP), intent(out) :: fvz    (ijdim,kdim)
    real(RP), intent(out) :: fe     (ijdim,kdim)
    real(RP), intent(out) :: fq     (ijdim,kdim,TRC_VMAX)
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
    real(RP) :: t    (ijdim,vlayer)   ! Temperature at full-model level (K)
    real(RP) :: qvv  (ijdim,vlayer)   ! Specific Humidity at full-model level (kg/kg)
    real(RP) :: u    (ijdim,vlayer)   ! Zonal wind at full-model level (m/s)
    real(RP) :: v    (ijdim,vlayer)   ! Meridional wind at full-model level (m/s)
    real(RP) :: pmid (ijdim,vlayer)   ! Pressure is full-model level (Pa)
    real(RP) :: pint (ijdim,vlayer+1) ! Pressure at model interfaces (Pa)
    real(RP) :: pdel (ijdim,vlayer)   ! Layer thickness (Pa)
    real(RP) :: rpdel(ijdim,vlayer)   ! Reciprocal of layer thickness (1/Pa)
    real(RP) :: ps   (ijdim)          ! surface pressure output [dummy]
    real(RP) :: precip2 (ijdim)
    integer  :: test

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

    integer  :: ij, k, kk
    !---------------------------------------------------------------------------

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
                + qv(:) * CVW(I_QV) &
                + qc(:) * CVW(I_QC) &
                + qr(:) * CVW(I_QR)

          fq(ij,kmin:kmax,I_QV) = fq(ij,kmin:kmax,I_QV) + ( qv(:) - q(ij,kmin:kmax,I_QV) ) / dt
          fq(ij,kmin:kmax,I_QC) = fq(ij,kmin:kmax,I_QC) + ( qc(:) - q(ij,kmin:kmax,I_QC) ) / dt
          fq(ij,kmin:kmax,I_QR) = fq(ij,kmin:kmax,I_QR) + ( qr(:) - q(ij,kmin:kmax,I_QR) ) / dt
          fe(ij,kmin:kmax)      = fe(ij,kmin:kmax)      + ( cv(:)*theta(:)*pk(:) - ein(ij,kmin:kmax) ) / dt
       enddo
    endif

    if ( USE_SimpleMicrophys ) then
       if ( SM_Latdepend_SST ) then
          test = 1
       else
          test = 0
       endif

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

       ps(:) = pre_sfc(:)

       call simple_physics( ijdim,             & ! [IN]
                            vlayer,            & ! [IN]
                            dt,                & ! [IN]
                            lat   (:),         & ! [IN]
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
                            test,              & ! [IN]
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
             if    ( RAIN_TYPE == 'DRY' ) then
                qv(k) = qvv(ij,k)

                qd(k) = 1.0_RP &
                      - qv(k)

                cv(k) = qd(k) * CVdry     &
                      + qv(k) * CVW(I_QV)
             elseif( RAIN_TYPE == 'WARM' ) then
                qv(k) = qvv(ij,k)
                qc(k) = q  (ij,kk,I_QC)
                qr(k) = q  (ij,kk,I_QR)

                qd(k) = 1.0_RP &
                      - qv(k)  &
                      - qc(k)  &
                      - qr(k)

                cv(k) = qd(k) * CVdry     &
                      + qv(k) * CVW(I_QV) &
                      + qc(k) * CVW(I_QC) &
                      + qr(k) * CVW(I_QR)
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
