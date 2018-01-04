!-------------------------------------------------------------------------------
!> Module Atmospheric Physics forcing
!!
!! @par Description
!!          This module contains subroutines for physical forcing terms
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_af_atmos_phy
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
  public :: af_atmos_phy

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine af_atmos_phy( &
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
       w,       &
       q,       &
       ein,     &
       pre_sfc, &
       fvx,     &
       fvy,     &
       fvz,     &
       fw,      &
       fe,      &
       fq,      &
       frho,    &
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
       kdim   => ADM_kall
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
    use mod_runconf, only: &
       TRC_VMAX,  &
       RAIN_TYPE, &
       NQW_MAX,   &
       I_QV,      &
       I_QC,      &
       I_QR,      &
       I_QI,      &
       I_QS,      &
       I_QG,      &
       NCHEM_STR, &
       NCHEM_END, &
       CVW,       &
       ATMOS_PHY_TYPE
    use mod_af_heldsuarez, only: &
       AF_heldsuarez
    use scale_tracer, only: &
       QA,          &
       TRACER_MASS, &
       TRACER_R,    &
       TRACER_CV,   &
       TRACER_CP
    use scale_atmos_phy_mp, only: &
       QS_MP,                        &
       QE_MP
    use mod_atmos_phy_mp_driver, only: &
       atmos_phy_mp_step
    use mod_mp_vars, only: &
       CCN_rmgrid,  &
       sflx_rain_rmgrid, &
       sflx_snow_rmgrid, &
       Evaporate_rmgrid
    use scale_grid_index, only: &
       IA,          &
       JA,          &
       KS,          &
       KE,          &
       KA
    use mod_atmos_vars, only: &
       RHO_rmgrid  => DENS,       &
       MOMZ_rmgrid => MOMZ,      &
       MOMX_rmgrid => MOMX,      &
       MOMY_rmgrid => MOMY,      &
       RHOT_rmgrid => RHOT,      &
       Q_rmgrid    => QTRC,         &
       ATMOS_vars_calc_diagnostics, &
       QV0 => QV, &
       QC0 => QC, &
       QI0 => QI
    use mod_atmos_phy_bl_driver, only: &
       ATMOS_PHY_BL_step
    use mod_atmos_phy_sf_driver, only: &
       ATMOS_PHY_SF_step
    use mod_atmos_phy_rd_driver, only: &
       ATMOS_PHY_RD_step
    use mod_atmos_admin, only: &
       ATMOS_sw_phy_mp,   &
       ATMOS_sw_phy_rd,   &
       ATMOS_sw_phy_sf,   &
       ATMOS_sw_phy_bl
    use mod_grd_conversion, only: &
       grd_gm2rm,  &
       grd_rm2gm
    use mod_time, only: &
       do_phy_rd => TIME_DOATMOS_PHY_RD
    use mod_af_dcmip, only: &
       USE_HeldSuarez
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
    real(RP), intent(in)  :: w      (ijdim,kdim)
    real(RP), intent(in)  :: q      (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: ein    (ijdim,kdim)
    real(RP), intent(in)  :: pre_sfc(ijdim)
    real(RP), intent(out) :: fvx    (ijdim,kdim)
    real(RP), intent(out) :: fvy    (ijdim,kdim)
    real(RP), intent(out) :: fvz    (ijdim,kdim)
    real(RP), intent(out) :: fw     (ijdim,kdim)
    real(RP), intent(out) :: fe     (ijdim,kdim)
    real(RP), intent(out) :: fq     (ijdim,kdim,TRC_VMAX)
    real(RP), intent(out) :: frho   (ijdim,kdim)
    real(RP), intent(out) :: precip (ijdim)
    real(RP), intent(in)  :: ix     (ijdim)
    real(RP), intent(in)  :: iy     (ijdim)
    real(RP), intent(in)  :: iz     (ijdim)
    real(RP), intent(in)  :: jx     (ijdim)
    real(RP), intent(in)  :: jy     (ijdim)
    real(RP), intent(in)  :: jz     (ijdim)
    real(DP), intent(in)  :: dt

    ! for kessler
    real(RP) :: theta (kdim) ! potential temperature (K)
    real(RP) :: qv    (kdim) ! water vapor mixing ratio (gm/gm)
    real(RP) :: qc    (kdim) ! cloud water mixing ratio (gm/gm)
    real(RP) :: qr    (kdim) ! rain  water mixing ratio (gm/gm)
    real(RP) :: rhod  (kdim) ! dry air density (not mean state as in KW) (kg/m^3)
    real(RP) :: pk    (kdim) ! Exner function (p/p0)**(R/cp)
    real(RP) :: z     (kdim) ! heights of thermo. levels in the grid column (m)
    real(RP) :: qd    (kdim)
    real(RP) :: cv    (kdim)

    ! for Tomita08
    real(RP) :: CCN                 (ijdim,kdim)
    real(RP) :: Evaporate           (ijdim,kdim)
    real(RP) :: sflx_rain           (ijdim)
    real(RP) :: sflx_snow           (ijdim)

    ! for simple physics
    real(RP) :: t    (ijdim,kdim)   ! Temperature at full-model level (K)
    real(RP) :: qvv  (ijdim,kdim)   ! Specific Humidity at full-model level (kg/kg)
    real(RP) :: u    (ijdim,kdim)   ! Zonal wind at full-model level (m/s)
    real(RP) :: v    (ijdim,kdim)   ! Meridional wind at full-model level (m/s)
    real(RP) :: pmid (ijdim,kdim)   ! Pressure is full-model level (Pa)
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

    ! RM-grid
    real(RP) :: frho_rmgrid         (KA,IA,JA)
    real(RP) :: fmomz_rmgrid        (KA,IA,JA)
    real(RP) :: fmomx_rmgrid        (KA,IA,JA)
    real(RP) :: fmomy_rmgrid        (KA,IA,JA)
    real(RP) :: frhot_rmgrid        (KA,IA,JA)
    real(RP) :: frhoq_rmgrid        (KA,IA,JA,QA)
    real(RP) :: frho_rmgrid_tot     (KA,IA,JA)
    real(RP) :: fmomz_rmgrid_tot    (KA,IA,JA)
    real(RP) :: fmomx_rmgrid_tot    (KA,IA,JA)
    real(RP) :: fmomy_rmgrid_tot    (KA,IA,JA)
    real(RP) :: frhot_rmgrid_tot    (KA,IA,JA)
    real(RP) :: frhoq_rmgrid_tot    (KA,IA,JA,QA)

    ! real physics
    real(RP) :: momz                (ijdim,kdim)
    real(RP) :: momx                (ijdim,kdim)
    real(RP) :: momy                (ijdim,kdim)
    real(RP) :: rhot                (ijdim,kdim)

    ! forcing terms
    integer  :: ij, k, kk, iq
    real(RP) :: fq_rmgrid  (ka,ia,ja,qa)
    !---------------------------------------------------------------------------

    fvx(:,:)   = 0.0_RP
    fvy(:,:)   = 0.0_RP
    fvz(:,:)   = 0.0_RP
    fw (:,:)   = 0.0_RP
    fe (:,:)   = 0.0_RP
    fq (:,:,:) = 0.0_RP
    frho(:,:)  = 0.0_RP

    frho_rmgrid_tot (:,:,:)   = 0.0_RP
    fmomz_rmgrid_tot(:,:,:)   = 0.0_RP
    fmomx_rmgrid_tot(:,:,:)   = 0.0_RP
    fmomy_rmgrid_tot(:,:,:)   = 0.0_RP
    frhot_rmgrid_tot(:,:,:)   = 0.0_RP    
    frhoq_rmgrid_tot(:,:,:,:) = 0.0_RP

    precip (:) = 0.0_RP
    precip2(:) = 0.0_RP

    do k = 1, kdim
       u   (:,k) = vx (:,k) * ix(:) + vy (:,k) * iy(:) + vz (:,k) * iz(:)
       v   (:,k) = vx (:,k) * jx(:) + vy (:,k) * jy(:) + vz (:,k) * jz(:)
    end do

    momz(:,:) = rho(:,:) * w(:,:)
    momx(:,:) = rho(:,:) * u(:,:)
    momy(:,:) = rho(:,:) * v(:,:)
    do ij=1, ijdim
       pk   (:) = ( pre(ij,:) / PRE00 )**( Rdry / CPdry )
       theta(:) = tem(ij,:) / pk(:)
       rhot(ij,:) = rho(ij,:) * theta(:)
    enddo

    call grd_gm2rm( rho,  rho_rmgrid,  kdim )
    call grd_gm2rm( momz, momz_rmgrid, kdim )
    call grd_gm2rm( momx, momx_rmgrid, kdim )
    call grd_gm2rm( momy, momy_rmgrid, kdim )
    call grd_gm2rm( rhot, rhot_rmgrid, kdim )
    do iq = 1, QA
       call grd_gm2rm( q(:,:,iq), q_rmgrid(:,:,:,iq), kdim )
    enddo

    call atmos_vars_calc_diagnostics

    if( atmos_sw_phy_mp ) then
       call grd_gm2rm( ccn(:,:),  ccn_rmgrid, kdim )

       call atmos_phy_mp_step(&
            frho_rmgrid,      & ! [OUT]
            fmomz_rmgrid,     & ! [OUT]
            fmomx_rmgrid,     & ! [OUT]
            fmomy_rmgrid,     & ! [OUT]
            frhot_rmgrid,     & ! [OUT]
            frhoq_rmgrid,     & ! [OUT]
            sflx_rain_rmgrid, & ! [OUT]
            sflx_snow_rmgrid )
       
       call grd_rm2gm( sflx_rain_rmgrid(:,:),  sflx_rain, 1 )
       call grd_rm2gm( sflx_snow_rmgrid(:,:),  sflx_snow, 1 )

       frho_rmgrid_tot (:,:,:)   = frho_rmgrid_tot(:,:,:)    + frho_rmgrid (:,:,:)
       fmomz_rmgrid_tot(:,:,:)   = fmomz_rmgrid_tot(:,:,:)   + fmomz_rmgrid(:,:,:)
       fmomx_rmgrid_tot(:,:,:)   = fmomz_rmgrid_tot(:,:,:)   + fmomx_rmgrid(:,:,:)
       fmomy_rmgrid_tot(:,:,:)   = fmomy_rmgrid_tot(:,:,:)   + fmomy_rmgrid(:,:,:)
       frhot_rmgrid_tot(:,:,:)   = frhot_rmgrid_tot(:,:,:)   + frhot_rmgrid(:,:,:)
       frhoq_rmgrid_tot(:,:,:,:) = frhoq_rmgrid_tot(:,:,:,:) + frhoq_rmgrid(:,:,:,:)
       
       precip(:) = sflx_rain(:) + sflx_snow(:)

    end if

    if( atmos_sw_phy_rd ) then
       call atmos_phy_rd_step( &   
            do_phy_rd, &    ![IN]
            frhot_rmgrid )  ![OUT]

       frhot_rmgrid_tot(:,:,:) = frhot_rmgrid_tot(:,:,:) + frhot_rmgrid(:,:,:)            
    endif

    if( atmos_sw_phy_sf ) then
       call atmos_phy_sf_step(&
            frho_rmgrid  (KS,:,:),      & ! [OUT]
            fmomz_rmgrid (KS,:,:),      & ! [OUT]
            fmomx_rmgrid (KS,:,:),      & ! [OUT]
            fmomy_rmgrid (KS,:,:),      & ! [OUT]
            frhot_rmgrid (KS,:,:),      & ! [OUT]
            frhoq_rmgrid (KS,:,:,I_QV)  ) ! [OUT]

       frho_rmgrid_tot(KS,:,:)    = frho_rmgrid_tot(KS,:,:)    + frho_rmgrid (KS,:,:)
       fmomz_rmgrid_tot(KS,:,:)   = fmomz_rmgrid_tot(KS,:,:)   + fmomz_rmgrid(KS,:,:)
       fmomx_rmgrid_tot(KS,:,:)   = fmomz_rmgrid_tot(KS,:,:)   + fmomx_rmgrid(KS,:,:)
       fmomy_rmgrid_tot(KS,:,:)   = fmomy_rmgrid_tot(KS,:,:)   + fmomy_rmgrid(KS,:,:)
       frhot_rmgrid_tot(KS,:,:)   = frhot_rmgrid_tot(KS,:,:)   + frhot_rmgrid(KS,:,:)
       frhoq_rmgrid_tot(KS,:,:,:) = frhoq_rmgrid_tot(KS,:,:,:) + frhoq_rmgrid   (KS,:,:,:)
    endif

    if( atmos_sw_phy_bl ) then
       call atmos_phy_bl_step(&
            fmomx_rmgrid,   & ! [OUT]
            fmomy_rmgrid,   & ! [OUT]
            frhot_rmgrid,   & ! [OUT]
            frhoq_rmgrid  )   ! [OUT]
       fmomx_rmgrid_tot(:,:,:)   = fmomz_rmgrid_tot(:,:,:)   + fmomx_rmgrid(:,:,:)
       fmomy_rmgrid_tot(:,:,:)   = fmomy_rmgrid_tot(:,:,:)   + fmomy_rmgrid(:,:,:)
       frhot_rmgrid_tot(:,:,:)   = frhot_rmgrid_tot(:,:,:)   + frhot_rmgrid(:,:,:)
       frhoq_rmgrid_tot(:,:,:,:) = frhoq_rmgrid_tot(:,:,:,:) + frhoq_rmgrid(:,:,:,:)
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

    call forcing_convert_rm2gm( &
         ijdim,            & ! [IN]
         kdim,             & ! [IN]
         frho_rmgrid_tot,  & ! [IN]
         frhot_rmgrid_tot, & ! [IN]
         fmomz_rmgrid_tot, & ! [IN]
         fmomx_rmgrid_tot, & ! [IN]
         fmomy_rmgrid_tot, & ! [IN]
         frhoq_rmgrid_tot, & ! [IN]
         frho,             & ! [OUT]
         fvx,              & ! [OUT]
         fvy,              & ! [OUT]
         fvz,              & ! [OUT]
         fw,               & ! [OUT]
         fe,               & ! [OUT]
         fq,               & ! [OUT]
         rho,              & ! [IN]
         q,                & ! [IN]
         vx,               & ! [IN]
         vy,               & ! [IN]
         vz,               & ! [IN]
         w,                & ! [IN]
         pre,              & ! [IN]
         tem,              & ! [IN]
         ix,               & ! [IN]
         iy,               & ! [IN]
         iz,               & ! [IN]
         jx,               & ! [IN]
         jy,               & ! [IN]
         jz                & ! [IN]
         )

    return
  end subroutine af_atmos_phy
  !-----------------------------------------------------------------------------
  subroutine forcing_convert_rm2gm(&
       ijdim,        &
       kdim,         &
       frho_rmgrid,  &
       frhot_rmgrid, &
       fmomz_rmgrid, &
       fmomx_rmgrid, &
       fmomy_rmgrid, &
       frhoq_rmgrid, &
       frho,         &
       fvx,          &
       fvy,          &
       fvz,          &
       fw,           &
       fe,           &
       fq,           &
       rho,          &
       q,            &
       vx,           &
       vy,           &
       vz,           &
       w,            &
       pre,          &
       tem,          &
       ix,           &
       iy,           &
       iz,           &
       jx,           &
       jy,           &
       jz            &
       )

    use scale_grid_index, only: &
       IA,          &
       JA,          &
       KA
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       CPdry => CONST_CPdry, &
       CVdry => CONST_CVdry, &
       PRE00 => CONST_PRE00
    use mod_runconf, only: &
       TRC_VMAX,  &
       NQW_MAX,   &
       CVW
    use mod_grd_conversion, only: &
       grd_gm2rm,  &
       grd_rm2gm

    implicit none

    integer,  intent(in) :: ijdim, kdim

    real(RP), intent(in) :: frho_rmgrid         (KA,IA,JA)
    real(RP), intent(in) :: fmomz_rmgrid        (KA,IA,JA)
    real(RP), intent(in) :: fmomx_rmgrid        (KA,IA,JA)
    real(RP), intent(in) :: fmomy_rmgrid        (KA,IA,JA)
    real(RP), intent(in) :: frhot_rmgrid        (KA,IA,JA)
    real(RP), intent(in) :: frhoq_rmgrid        (KA,IA,JA,TRC_VMAX)

    real(RP), intent(out) :: fvx    (ijdim,kdim)
    real(RP), intent(out) :: fvy    (ijdim,kdim)
    real(RP), intent(out) :: fvz    (ijdim,kdim)
    real(RP), intent(out) :: fw     (ijdim,kdim)
    real(RP), intent(out) :: fe     (ijdim,kdim)
    real(RP), intent(out) :: fq     (ijdim,kdim,TRC_VMAX)
    real(RP), intent(out) :: frho   (ijdim,kdim)

    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(in)  :: q      (ijdim,kdim,TRC_VMAX)
    real(RP), intent(in)  :: vx     (ijdim,kdim)
    real(RP), intent(in)  :: vy     (ijdim,kdim)
    real(RP), intent(in)  :: vz     (ijdim,kdim)
    real(RP), intent(in)  :: w      (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: ix     (ijdim)
    real(RP), intent(in)  :: iy     (ijdim)
    real(RP), intent(in)  :: iz     (ijdim)
    real(RP), intent(in)  :: jx     (ijdim)
    real(RP), intent(in)  :: jy     (ijdim)
    real(RP), intent(in)  :: jz     (ijdim)

    real(RP) :: fmomz    (ijdim,kdim)
    real(RP) :: fmomx    (ijdim,kdim)
    real(RP) :: fmomy    (ijdim,kdim)
    real(RP) :: frhot    (ijdim,kdim)
    real(RP) :: fu       (ijdim,kdim)
    real(RP) :: fv       (ijdim,kdim)
    real(RP) :: fcv      (ijdim,kdim)
    real(RP) :: fqd      (ijdim,kdim)
    real(RP) :: frhoq    (ijdim,kdim,TRC_VMAX)
    real(RP) :: u        (ijdim,kdim)
    real(RP) :: v        (ijdim,kdim)
    real(RP) :: qd       (ijdim,kdim)
    real(RP) :: theta    (ijdim,kdim)
    real(RP) :: pk       (ijdim,kdim)
    real(RP) :: cv       (ijdim,kdim)

    integer :: iq, k

    !---------------------------------------------------------------------------

    call grd_rm2gm( frho_rmgrid , frho, kdim )
    call grd_rm2gm( fmomz_rmgrid, fmomz, kdim )
    call grd_rm2gm( fmomx_rmgrid, fmomx, kdim )
    call grd_rm2gm( fmomy_rmgrid, fmomy, kdim )
    call grd_rm2gm( frhot_rmgrid, frhot, kdim )
    do iq = 1, NQW_MAX
       call grd_rm2gm( frhoq_rmgrid(:,:,:,iq), frhoq(:,:,iq), kdim )
    enddo

    do k=1, kdim
       u(:,k) = vx (:,k) * ix(:) + vy (:,k) * iy(:) + vz (:,k) * iz(:)
       v(:,k) = vx (:,k) * jx(:) + vy (:,k) * jy(:) + vz (:,k) * jz(:)
    enddo

    fu(:,:) = fmomx(:,:) / rho(:,:)
    fv(:,:) = fmomy(:,:) / rho(:,:)
    do k=1, kdim
       fvx(:,k) = fu(:,k) * ix(:) + fv(:,k) * jx(:)
       fvy(:,k) = fu(:,k) * iy(:) + fv(:,k) * jy(:)
       fvz(:,k) = fu(:,k) * iz(:) + fv(:,k) * jz(:)
    enddo
    fw(:,:) = fmomz(:,:) / rho(:,:)

    do iq=1, NQW_MAX
       fq(:,:,iq) =  frhoq(:,:,iq) / rho(:,:)
    enddo

    fqd(:,:) = 0.0_RP
    do iq=1, NQW_MAX
       fqd(:,:) = fqd(:,:) - frhoq(:,:,iq) / rho(:,:)
    enddo
    
    fcv(:,:) = CVdry * fqd(:,:)
    do iq=1, NQW_MAX
       fcv(:,:) = fcv(:,:) + fq(:,:,iq) * CVW(iq)
    enddo

    pk(:,:) = ( pre(:,:) / PRE00 )**( Rdry / CPdry )
    theta(:,:) = tem(:,:) / pk(:,:)

    qd(:,:) = 1.0_RP
    do iq=1, NQW_MAX
       qd(:,:) = qd(:,:) - q(:,:,iq)
    enddo

    cv(:,:) = CVdry * qd(:,:)
    do iq=1, NQW_MAX
       cv(:,:) = cv(:,:) + q(:,:,iq) * CVW(iq)
    enddo

    fe(:,:) = theta(:,:) * pk(:,:) * fcv(:,:) + CV * pk(:,:) * frhot(:,:) / rho(:,:)

  end subroutine forcing_convert_rm2gm

end module mod_af_atmos_phy
