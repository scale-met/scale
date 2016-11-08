!-------------------------------------------------------------------------------
!> Module source
!!
!! @par Description
!!          This module is for the caluculation of source terms
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_src
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: src_advection_convergence_momentum
  public :: src_advection_convergence
  public :: src_flux_convergence

  public :: src_pres_gradient
  public :: src_buoyancy

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: I_SRC_horizontal = 1
  integer, public, parameter :: I_SRC_vertical   = 2
  integer, public, parameter :: I_SRC_default    = 3

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  logical, private, parameter :: first_layer_remedy = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Advection convergence for momentum
  subroutine src_advection_convergence_momentum( &
       vx,      vx_pl,      &
       vy,      vy_pl,      &
       vz,      vz_pl,      &
       w,       w_pl,       &
       rhog,    rhog_pl,    &
       rhogvx,  rhogvx_pl,  &
       rhogvy,  rhogvy_pl,  &
       rhogvz,  rhogvz_pl,  &
       rhogw,   rhogw_pl,   &
       grhogvx, grhogvx_pl, &
       grhogvy, grhogvy_pl, &
       grhogvz, grhogvz_pl, &
       grhogw,  grhogw_pl   )
    use scale_const, only: &
       OHM => CONST_OHM
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax,    &
       ADM_KNONE
    use mod_grd, only: &
       GRD_rscale, &
       GRD_XDIR,   &
       GRD_YDIR,   &
       GRD_ZDIR,   &
       GRD_x,      &
       GRD_x_pl,   &
       GRD_cfac,   &
       GRD_dfac
    use mod_vmtr, only: &
       VMTR_getIJ_C2Wfact
    use mod_runconf, only : &
       NON_HYDRO_ALPHA
    implicit none

    real(RP), intent(in)  :: vx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: vz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: vz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: w        (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: w_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)  :: rhog     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhog_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(out) :: grhogvx   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: grhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: grhogvy   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: grhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: grhogvz   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: grhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: grhogw    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: grhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: vvx       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vvx_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vvy       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vvy_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: vvz       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: vvz_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dvvx      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dvvx_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dvvy      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dvvy_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: dvvz      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: dvvz_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: grhogwc   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: grhogwc_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: VMTR_C2Wfact   (ADM_gall   ,ADM_kall,2,ADM_lall   )
    real(RP) :: VMTR_C2Wfact_pl(ADM_gall_pl,ADM_kall,2,ADM_lall_pl)

    real(RP) :: prd, wc

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('____src_advection_convergence_m',2)

    call VMTR_getIJ_C2Wfact( VMTR_C2Wfact, VMTR_C2Wfact_pl )

    !---< merge horizontal velocity & vertical velocity >

    do l = 1, ADM_lall
       do k = ADM_kmin,ADM_kmax
       do g = 1, ADM_gall
          wc = GRD_cfac(k) * w(g,k+1,l) &
             + GRD_dfac(k) * w(g,k  ,l)

          vvx(g,k,l) = vx(g,k,l) + wc * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
          vvy(g,k,l) = vy(g,k,l) + wc * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
          vvz(g,k,l) = vz(g,k,l) + wc * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
       enddo
       enddo

       do g = 1, ADM_gall
          vvx(g,ADM_kmin-1,l) = 0.0_RP
          vvx(g,ADM_kmax+1,l) = 0.0_RP
          vvy(g,ADM_kmin-1,l) = 0.0_RP
          vvy(g,ADM_kmax+1,l) = 0.0_RP
          vvz(g,ADM_kmin-1,l) = 0.0_RP
          vvz(g,ADM_kmax+1,l) = 0.0_RP
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             wc = GRD_cfac(k) * w_pl(g,k+1,l) &
                + GRD_dfac(k) * w_pl(g,k  ,l)

             vvx_pl(g,k,l) = vx_pl(g,k,l) + wc * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
             vvy_pl(g,k,l) = vy_pl(g,k,l) + wc * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
             vvz_pl(g,k,l) = vz_pl(g,k,l) + wc * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale
          enddo
          enddo

          do g = 1, ADM_gall_pl
             vvx_pl(g,ADM_kmin-1,l) = 0.0_RP
             vvx_pl(g,ADM_kmax+1,l) = 0.0_RP
             vvy_pl(g,ADM_kmin-1,l) = 0.0_RP
             vvy_pl(g,ADM_kmax+1,l) = 0.0_RP
             vvz_pl(g,ADM_kmin-1,l) = 0.0_RP
             vvz_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    endif

    !---< advection term for momentum >

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvx,    vvx_pl,    & ! [IN]
                                    dvvx,   dvvx_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvy,    vvy_pl,    & ! [IN]
                                    dvvy,   dvvy_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]

    call src_advection_convergence( rhogvx, rhogvx_pl, & ! [IN]
                                    rhogvy, rhogvy_pl, & ! [IN]
                                    rhogvz, rhogvz_pl, & ! [IN]
                                    rhogw,  rhogw_pl,  & ! [IN]
                                    vvz,    vvz_pl,    & ! [IN]
                                    dvvz,   dvvz_pl,   & ! [OUT]
                                    I_SRC_default      ) ! [IN]


    do l = 1, ADM_lall
       !---< coriolis force >
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          dvvx(g,k,l) = dvvx(g,k,l) - 2.0_RP * rhog(g,k,l) * ( -OHM * vvy(g,k,l) )
          dvvy(g,k,l) = dvvy(g,k,l) - 2.0_RP * rhog(g,k,l) * (  OHM * vvx(g,k,l) )
       enddo
       enddo

       !---< horizontalize & separate vertical velocity >
       do k = ADM_kmin, ADM_kmax
       do g = 1, ADM_gall
          prd = dvvx(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale &
              + dvvy(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale &
              + dvvz(g,k,l) * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale

          grhogvx(g,k,l) = dvvx(g,k,l) - prd * GRD_x(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
          grhogvy(g,k,l) = dvvy(g,k,l) - prd * GRD_x(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
          grhogvz(g,k,l) = dvvz(g,k,l) - prd * GRD_x(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale

          grhogwc(g,k,l) = prd * real(NON_HYDRO_ALPHA,kind=RP)
       enddo
       enddo

       do k = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          grhogw(g,k,l) = ( VMTR_C2Wfact(g,k,1,l) * grhogwc(g,k  ,l) &
                          + VMTR_C2Wfact(g,k,2,l) * grhogwc(g,k-1,l) )
       enddo
       enddo

       do g = 1, ADM_gall
          grhogvx(g,ADM_kmin-1,l) = 0.0_RP
          grhogvx(g,ADM_kmax+1,l) = 0.0_RP
          grhogvy(g,ADM_kmin-1,l) = 0.0_RP
          grhogvy(g,ADM_kmax+1,l) = 0.0_RP
          grhogvz(g,ADM_kmin-1,l) = 0.0_RP
          grhogvz(g,ADM_kmax+1,l) = 0.0_RP
          grhogw (g,ADM_kmin-1,l) = 0.0_RP
          grhogw (g,ADM_kmin  ,l) = 0.0_RP
          grhogw (g,ADM_kmax+1,l) = 0.0_RP
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          !---< coriolis force >
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             dvvx_pl(g,k,l) = dvvx_pl(g,k,l) - 2.0_RP * rhog_pl(g,k,l) * ( -OHM * vvy_pl(g,k,l) )
             dvvy_pl(g,k,l) = dvvy_pl(g,k,l) - 2.0_RP * rhog_pl(g,k,l) * (  OHM * vvx_pl(g,k,l) )
          enddo
          enddo

          !---< horizontalize & separate vertical velocity >
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             prd = dvvx_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale &
                 + dvvy_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale &
                 + dvvz_pl(g,k,l) * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale

             grhogvx_pl(g,k,l) = dvvx_pl(g,k,l) - prd * GRD_x_pl(g,ADM_KNONE,l,GRD_XDIR) / GRD_rscale
             grhogvy_pl(g,k,l) = dvvy_pl(g,k,l) - prd * GRD_x_pl(g,ADM_KNONE,l,GRD_YDIR) / GRD_rscale
             grhogvz_pl(g,k,l) = dvvz_pl(g,k,l) - prd * GRD_x_pl(g,ADM_KNONE,l,GRD_ZDIR) / GRD_rscale

             grhogwc_pl(g,k,l) = prd * real(NON_HYDRO_ALPHA,kind=RP)
          enddo
          enddo

          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             grhogw_pl(g,k,l) = ( VMTR_C2Wfact_pl(g,k,1,l) * grhogwc_pl(g,k  ,l) &
                                + VMTR_C2Wfact_pl(g,k,2,l) * grhogwc_pl(g,k-1,l) )
          enddo
          enddo

          do g = 1, ADM_gall_pl
             grhogvx_pl(g,ADM_kmin-1,l) = 0.0_RP
             grhogvx_pl(g,ADM_kmax+1,l) = 0.0_RP
             grhogvy_pl(g,ADM_kmin-1,l) = 0.0_RP
             grhogvy_pl(g,ADM_kmax+1,l) = 0.0_RP
             grhogvz_pl(g,ADM_kmin-1,l) = 0.0_RP
             grhogvz_pl(g,ADM_kmax+1,l) = 0.0_RP
             grhogw_pl (g,ADM_kmin-1,l) = 0.0_RP
             grhogw_pl (g,ADM_kmin  ,l) = 0.0_RP
             grhogw_pl (g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    else
       grhogvx_pl(:,:,:) = 0.0_RP
       grhogvy_pl(:,:,:) = 0.0_RP
       grhogvz_pl(:,:,:) = 0.0_RP
       grhogw_pl (:,:,:) = 0.0_RP
    endif

    call PROF_rapend('____src_advection_convergence_m',2)

    return
  end subroutine src_advection_convergence_momentum

  !-----------------------------------------------------------------------------
  !> Advection convergence
  subroutine src_advection_convergence( &
       rhogvx,   rhogvx_pl,   &
       rhogvy,   rhogvy_pl,   &
       rhogvz,   rhogvz_pl,   &
       rhogw,    rhogw_pl,    &
       scl,      scl_pl,      &
       grhogscl, grhogscl_pl, &
       fluxtype               )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_afac, &
       GRD_bfac
    implicit none

    real(RP), intent(in)  :: rhogvx     (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vx ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvx_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy     (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vy ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvy_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz     (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vz ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvz_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw      (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w  ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogw_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl        (ADM_gall   ,ADM_kall,ADM_lall   ) ! scalar
    real(RP), intent(in)  :: scl_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: grhogscl   (ADM_gall   ,ADM_kall,ADM_lall   ) ! scalar tendency
    real(RP), intent(out) :: grhogscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in)  :: fluxtype                                      ! scheme type
                                                                           ! I_SRC_default    : horizontal & vertical convergence
                                                                           ! I_SRC_horizontal : horizontal convergence

    real(RP) :: rhogvxscl   (ADM_gall   ,ADM_kall,ADM_lall   ) ! scalar * rho*Vx ( G^1/2 x gam2 )
    real(RP) :: rhogvxscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvyscl   (ADM_gall   ,ADM_kall,ADM_lall   ) ! scalar * rho*Vy ( G^1/2 x gam2 )
    real(RP) :: rhogvyscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvzscl   (ADM_gall   ,ADM_kall,ADM_lall   ) ! scalar * rho*Vz ( G^1/2 x gam2 )
    real(RP) :: rhogvzscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogwscl    (ADM_gall   ,ADM_kall,ADM_lall   ) ! scalar * rho*w  ( G^1/2 x gam2 )
    real(RP) :: rhogwscl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call PROF_rapstart('____src_advection_convergence',2)

    ! rhogvh * scl
!OCL SERIAL
    do l = 1, ADM_lall
!OCL PARALLEL
    do k = 1, ADM_kall
    do g = 1, ADM_gall
       rhogvxscl(g,k,l) = rhogvx(g,k,l) * scl(g,k,l)
       rhogvyscl(g,k,l) = rhogvy(g,k,l) * scl(g,k,l)
       rhogvzscl(g,k,l) = rhogvz(g,k,l) * scl(g,k,l)
    enddo
    enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          rhogvxscl_pl(g,k,l) = rhogvx_pl(g,k,l) * scl_pl(g,k,l)
          rhogvyscl_pl(g,k,l) = rhogvy_pl(g,k,l) * scl_pl(g,k,l)
          rhogvzscl_pl(g,k,l) = rhogvz_pl(g,k,l) * scl_pl(g,k,l)
       enddo
       enddo
       enddo
    endif

    ! rhogw * scl at half level
    if ( fluxtype == I_SRC_default ) then

       do l = 1, ADM_lall
          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall
             rhogwscl(g,k,l) = rhogw(g,k,l) * ( GRD_afac(k) * scl(g,k,  l) &
                                              + GRD_bfac(k) * scl(g,k-1,l) )
          enddo
          enddo
          do g = 1, ADM_gall
             rhogwscl(g,ADM_kmin-1,l) = 0.0_RP
          enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin, ADM_kmax+1
             do g = 1, ADM_gall_pl
                rhogwscl_pl(g,k,l) = rhogw_pl(g,k,l) * ( GRD_afac(k) * scl_pl(g,k  ,l) &
                                                       + GRD_bfac(k) * scl_pl(g,k-1,l) )
             enddo
             enddo
             do g = 1, ADM_gall_pl
                rhogwscl_pl(g,ADM_kmin-1,l) = 0.0_RP
             enddo
          enddo
       endif

    elseif( fluxtype == I_SRC_horizontal ) then

!OCL SERIAL
       do l = 1, ADM_lall
!OCL PARALLEL
       do k = 1, ADM_kall
       do g = 1, ADM_gall
          rhogwscl(g,k,l) = 0.0_RP
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             rhogwscl_pl(g,k,l) = 0.0_RP
          enddo
          enddo
          enddo
       endif

    endif

    !--- flux convergence step
    call src_flux_convergence( rhogvxscl, rhogvxscl_pl, & ! [IN]
                               rhogvyscl, rhogvyscl_pl, & ! [IN]
                               rhogvzscl, rhogvzscl_pl, & ! [IN]
                               rhogwscl,  rhogwscl_pl,  & ! [IN]
                               grhogscl,  grhogscl_pl,  & ! [OUT]
                               fluxtype                 ) ! [IN]

    call PROF_rapend('____src_advection_convergence',2)

    return
  end subroutine src_advection_convergence

  !-----------------------------------------------------------------------------
  !> Flux convergence calculation
  !! 1. Horizontal flux convergence is calculated by using rhovx, rhovy, and
  !!    rhovz which are defined at cell center (vertical) and A-grid (horizontal).
  !! 2. Vertical flux convergence is calculated by using rhovx, rhovy, rhovz, and rhow.
  !! 3. rhovx, rhovy, and rhovz can be replaced by rhovx*h, rhovy*h, and rhovz*h, respectively.
  subroutine src_flux_convergence( &
       rhogvx, rhogvx_pl, &
       rhogvy, rhogvy_pl, &
       rhogvz, rhogvz_pl, &
       rhogw,  rhogw_pl,  &
       grhog,  grhog_pl,  &
       fluxtype           )
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_iall,    &
       ADM_jall,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       GRD_rdgz
    use mod_oprt, only: &
       OPRT_divergence, &
       OPRT_coef_div,   &
       OPRT_coef_div_pl
    use mod_vmtr, only: &
       VMTR_RGSQRTH,     &
       VMTR_RGSQRTH_pl,  &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_C2WfactGz,   &
       VMTR_C2WfactGz_pl
    implicit none

    real(RP), intent(in)  :: rhogvx   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vx ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvx_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvy   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vy ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvy_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogvz   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*Vz ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogvz_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w  ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: grhog    (ADM_gall   ,ADM_kall,ADM_lall   ) ! source
    real(RP), intent(out) :: grhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in)  :: fluxtype ! scheme type
                                      ! I_SRC_default    : horizontal & vertical convergence
                                      ! I_SRC_horizontal : horizontal convergence

    real(RP) :: div_rhogvh   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! horizontal convergence
    real(RP) :: div_rhogvh_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

    real(RP) :: rhogvx_vm   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*vx / vertical metrics
    real(RP) :: rhogvx_vm_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvy_vm   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*vy / vertical metrics
    real(RP) :: rhogvy_vm_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogvz_vm   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*vz / vertical metrics
    real(RP) :: rhogvz_vm_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: rhogw_vmh   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) ! rho*w  / vertical metrics
    real(RP) :: rhogw_vmh_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl)

    real(RP) :: vertical_flag

    integer  :: iall, jall, kall, kmin, kmax, lall

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    call PROF_rapstart('____src_flux_convergence',2)

    iall = ADM_iall
    jall = ADM_jall
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    if ( fluxtype == I_SRC_default ) then ! Default
       vertical_flag = 1.0_RP
    elseif( fluxtype == I_SRC_horizontal ) then ! Horizontal
       vertical_flag = 0.0_RP
    endif

    !$omp parallel default(none),private(i,j,k,l,g), &
    !$omp shared(iall,jall,kall,kmin,kmax,lall,rhogvx_vm,rhogvy_vm,rhogvz_vm,rhogw_vmh, &
    !$omp rhogvx,rhogvy,rhogvz,rhogw,VMTR_RGAM,VMTR_C2WfactGz,VMTR_RGAMH,VMTR_RGSQRTH,vertical_flag)
    do l = 1, lall
       !--- Horizontal flux
       !$omp do schedule(static)
       do k = 1, kall
          g = 1
          do j = 1, jall
          do i = 1, iall
             rhogvx_vm(i,j,k,l) = rhogvx(g,k,l) * VMTR_RGAM(i,j,k,l)
             rhogvy_vm(i,j,k,l) = rhogvy(g,k,l) * VMTR_RGAM(i,j,k,l)
             rhogvz_vm(i,j,k,l) = rhogvz(g,k,l) * VMTR_RGAM(i,j,k,l)
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do nowait

       !--- Vertical flux
       !$omp do schedule(static)
       do k = kmin+1, kmax
          g = 1
          do j = 1, jall
          do i = 1, iall
             rhogw_vmh(i,j,k,l) = ( VMTR_C2WfactGz(i,j,k,1,l) * rhogvx(g,k  ,l) &
                                  + VMTR_C2WfactGz(i,j,k,2,l) * rhogvx(g,k-1,l) &
                                  + VMTR_C2WfactGz(i,j,k,3,l) * rhogvy(g,k  ,l) &
                                  + VMTR_C2WfactGz(i,j,k,4,l) * rhogvy(g,k-1,l) &
                                  + VMTR_C2WfactGz(i,j,k,5,l) * rhogvz(g,k  ,l) &
                                  + VMTR_C2WfactGz(i,j,k,6,l) * rhogvz(g,k-1,l) &
                                  ) * VMTR_RGAMH(i,j,k,l)                       &      ! horizontal contribution
                                + vertical_flag * rhogw(g,k,l) * VMTR_RGSQRTH(i,j,k,l) ! vertical   contribution
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do

       !$omp workshare
       rhogw_vmh(:,:,kmin  ,l) = 0.0_RP
       rhogw_vmh(:,:,kmax+1,l) = 0.0_RP
       !$omp end workshare
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          !--- Horizontal flux
          do k = 1, ADM_kall
          do g = 1, ADM_gall_pl
             rhogvx_vm_pl(g,k,l) = rhogvx_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)
             rhogvy_vm_pl(g,k,l) = rhogvy_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)
             rhogvz_vm_pl(g,k,l) = rhogvz_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)
          enddo
          enddo

          !--- Vertical flux
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             rhogw_vmh_pl(g,k,l) = ( VMTR_C2WfactGz_pl(g,k,1,l) * rhogvx_pl(g,k  ,l) &
                                   + VMTR_C2WfactGz_pl(g,k,2,l) * rhogvx_pl(g,k-1,l) &
                                   + VMTR_C2WfactGz_pl(g,k,3,l) * rhogvy_pl(g,k  ,l) &
                                   + VMTR_C2WfactGz_pl(g,k,4,l) * rhogvy_pl(g,k-1,l) &
                                   + VMTR_C2WfactGz_pl(g,k,5,l) * rhogvz_pl(g,k  ,l) &
                                   + VMTR_C2WfactGz_pl(g,k,6,l) * rhogvz_pl(g,k-1,l) &
                                   ) * VMTR_RGAMH_pl(g,k,l)                          &      ! horizontal contribution
                                 + vertical_flag * rhogw_pl(g,k,l) * VMTR_RGSQRTH_pl(g,k,l) ! vertical   contribution
          enddo
          enddo
          do g = 1, ADM_gall_pl
             rhogw_vmh_pl(g,ADM_kmin  ,l) = 0.0_RP
             rhogw_vmh_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    endif

    !--- Horizontal flux convergence
    call OPRT_divergence( div_rhogvh   (:,:,:,:),   div_rhogvh_pl   (:,:,:), & ! [OUT]
                          rhogvx_vm    (:,:,:,:),   rhogvx_vm_pl    (:,:,:), & ! [IN]
                          rhogvy_vm    (:,:,:,:),   rhogvy_vm_pl    (:,:,:), & ! [IN]
                          rhogvz_vm    (:,:,:,:),   rhogvz_vm_pl    (:,:,:), & ! [IN]
                          OPRT_coef_div(:,:,:,:,:), OPRT_coef_div_pl(:,:,:)  ) ! [IN]

    !--- Total flux convergence
    !$omp parallel default(none),private(i,j,k,l,g), &
    !$omp shared(iall,jall,kall,kmin,kmax,lall,grhog,div_rhogvh,rhogw_vmh,GRD_rdgz)
    do l = 1, lall
       !$omp do schedule(static)
       do k = kmin, kmax
          g = 1
          do j = 1, jall
          do i = 1, iall
             grhog(g,k,l) = - div_rhogvh(i,j,k,l) &
                            - ( rhogw_vmh(i,j,k+1,l)-rhogw_vmh(i,j,k,l) ) * GRD_rdgz(k)
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do

       !$omp workshare
       grhog(:,kmin-1,l) = 0.0_RP
       grhog(:,kmax+1,l) = 0.0_RP
       !$omp end workshare
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             grhog_pl(g,k,l) = - div_rhogvh_pl(g,k,l) &
                               - ( rhogw_vmh_pl(g,k+1,l)-rhogw_vmh_pl(g,k,l) ) * GRD_rdgz(k)
          enddo
          enddo

          do g = 1, ADM_gall_pl
             grhog_pl(g,ADM_kmin-1,l) = 0.0_RP
             grhog_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    endif

    call PROF_rapend('____src_flux_convergence',2)

    return
  end subroutine src_flux_convergence

  !-----------------------------------------------------------------------------
  !> Gradient operator
  subroutine src_pres_gradient( &
       P,      P_pl,      &
       Pgrad,  Pgrad_pl,  &
       Pgradw, Pgradw_pl, &
       gradtype           )
    use mod_adm, only: &
       ADM_nxyz,    &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_iall,    &
       ADM_jall,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax,    &
       ADM_nxyz
    use mod_grd, only: &
       GRD_XDIR, &
       GRD_YDIR, &
       GRD_ZDIR, &
       GRD_rdgz, &
       GRD_rdgzh
    use mod_oprt, only: &
       OPRT_gradient,          &
       OPRT_horizontalize_vec, &
       OPRT_coef_grad,         &
       OPRT_coef_grad_pl
    use mod_vmtr, only: &
       VMTR_GAM2H,       &
       VMTR_GAM2H_pl,    &
       VMTR_RGAM,        &
       VMTR_RGAM_pl,     &
       VMTR_RGAMH,       &
       VMTR_RGAMH_pl,    &
       VMTR_RGSGAM2,     &
       VMTR_RGSGAM2_pl,  &
       VMTR_C2WfactGz,   &
       VMTR_C2WfactGz_pl
    implicit none

    real(RP), intent(in)  :: P        (ADM_gall   ,ADM_kall,ADM_lall   )          ! phi * G^1/2 * gamma^2
    real(RP), intent(in)  :: P_pl     (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: Pgrad    (ADM_gall   ,ADM_kall,ADM_lall   ,ADM_nxyz) ! horizontal gradient
    real(RP), intent(out) :: Pgrad_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl,ADM_nxyz)
    real(RP), intent(out) :: Pgradw   (ADM_gall   ,ADM_kall,ADM_lall   )          ! vertical gradient
    real(RP), intent(out) :: Pgradw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer,  intent(in)  :: gradtype ! scheme type
                                      ! I_SRC_default    : horizontal & vertical gradient
                                      ! I_SRC_horizontal : horizontal gradient

    real(RP) :: P_vm     (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP) :: P_vm_pl  (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP) :: P_vmh    (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ,ADM_nxyz)
    real(RP) :: P_vmh_pl (ADM_gall_pl      ,ADM_kall,ADM_lall_pl,ADM_nxyz)
    real(RP) :: Pgrad0   (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ,ADM_nxyz)
    real(RP) :: Pgrad0_pl(ADM_gall_pl      ,ADM_kall,ADM_lall_pl,ADM_nxyz)

    integer  :: nxyz, iall, jall, kall, kmin, kmax, lall

    integer :: i, j, k, l, d, g
    !---------------------------------------------------------------------------

    call PROF_rapstart('____src_pres_gradient',2)

    nxyz = ADM_nxyz
    iall = ADM_iall
    jall = ADM_jall
    kall = ADM_kall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    !---< horizontal gradient, horizontal contribution >---

    !$omp parallel default(none),private(i,j,k,l,g), &
    !$omp shared(iall,jall,kall,lall,P_vm,P,VMTR_RGAM)
    do l = 1, lall
       !$omp do schedule(static)
       do k = 1, kall
          g = 1
          do j = 1, jall
          do i = 1, iall
             P_vm(i,j,k,l) = P(g,k,l) * VMTR_RGAM(i,j,k,l)
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    if ( ADM_have_pl) then
       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall
       do g = 1, ADM_gall_pl
          P_vm_pl(g,k,l) = P_pl(g,k,l) * VMTR_RGAM_pl(g,k,l)
       enddo
       enddo
       enddo
    endif

    call OPRT_gradient( Pgrad0        (:,:,:,:,:), Pgrad0_pl        (:,:,:,:), & ! [OUT]
                        P_vm          (:,:,:,:),   P_vm_pl          (:,:,:),   & ! [IN]
                        OPRT_coef_grad(:,:,:,:,:), OPRT_coef_grad_pl(:,:,:)    ) ! [IN]

    !---< horizontal gradient, vertical contribution >---

    !$omp parallel default(none),private(i,j,k,l,g), &
    !$omp shared(iall,jall,kmin,kmax,lall,P_vmh,P,VMTR_C2WfactGz,VMTR_RGAMH)
    do l = 1, lall
       !$omp do schedule(static)
       do k = kmin, kmax+1
          g = 1
          do j = 1, jall
          do i = 1, iall
             P_vmh(i,j,k,l,GRD_XDIR) = ( VMTR_C2WfactGz(i,j,k,1,l) * P(g,k  ,l) &
                                       + VMTR_C2WfactGz(i,j,k,2,l) * P(g,k-1,l) ) * VMTR_RGAMH(i,j,k,l)
             P_vmh(i,j,k,l,GRD_YDIR) = ( VMTR_C2WfactGz(i,j,k,3,l) * P(g,k  ,l) &
                                       + VMTR_C2WfactGz(i,j,k,4,l) * P(g,k-1,l) ) * VMTR_RGAMH(i,j,k,l)
             P_vmh(i,j,k,l,GRD_ZDIR) = ( VMTR_C2WfactGz(i,j,k,5,l) * P(g,k  ,l) &
                                       + VMTR_C2WfactGz(i,j,k,6,l) * P(g,k-1,l) ) * VMTR_RGAMH(i,j,k,l)
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    !$omp parallel default(none),private(i,j,k,l,d,g), &
    !$omp shared(nxyz,iall,jall,kmin,kmax,lall,Pgrad,Pgrad0,P_vmh,GRD_rdgz)
    do d = 1, nxyz
    do l = 1, lall
       !$omp do schedule(static)
       do k = kmin, kmax
          g = 1
          do j = 1, jall
          do i = 1, iall
             Pgrad(g,k,l,d) = Pgrad0(i,j,k,l,d) + ( P_vmh(i,j,k+1,l,d) - P_vmh(i,j,k,l,d) ) * GRD_rdgz(k)
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do

       if ( first_layer_remedy ) then !--- At the lowest layer, do not use the extrapolation value
          !$omp workshare
          Pgrad(:,kmin,l,d) = Pgrad(:,kmin+1,l,d)
          !$omp end workshare
       endif

       !$omp workshare
       Pgrad(:,kmin-1,l,d) = 0.0_RP
       Pgrad(:,kmax+1,l,d) = 0.0_RP
       !$omp end workshare
    enddo
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
       do k = ADM_kmin, ADM_kmax+1
       do g = 1, ADM_gall_pl
          P_vmh_pl(g,k,l,GRD_XDIR) = ( VMTR_C2WfactGz_pl(g,k,1,l) * P_pl(g,k  ,l) &
                                     + VMTR_C2WfactGz_pl(g,k,2,l) * P_pl(g,k-1,l) ) * VMTR_RGAMH_pl(g,k,l)
          P_vmh_pl(g,k,l,GRD_YDIR) = ( VMTR_C2WfactGz_pl(g,k,3,l) * P_pl(g,k  ,l) &
                                     + VMTR_C2WfactGz_pl(g,k,4,l) * P_pl(g,k-1,l) ) * VMTR_RGAMH_pl(g,k,l)
          P_vmh_pl(g,k,l,GRD_ZDIR) = ( VMTR_C2WfactGz_pl(g,k,5,l) * P_pl(g,k  ,l) &
                                     + VMTR_C2WfactGz_pl(g,k,6,l) * P_pl(g,k-1,l) ) * VMTR_RGAMH_pl(g,k,l)
       enddo
       enddo
       enddo

       do d = 1, ADM_nxyz
       do l = 1, ADM_lall_pl
          do k = ADM_kmin, ADM_kmax
          do g = 1, ADM_gall_pl
             Pgrad_pl(g,k,l,d) = Pgrad0_pl(g,k,l,d) + ( P_vmh_pl(g,k+1,l,d) - P_vmh_pl(g,k,l,d) ) * GRD_rdgz(k)
          enddo
          enddo

          if ( first_layer_remedy ) then !--- At the lowest layer, do not use the extrapolation value!
             do g = 1, ADM_gall_pl
                Pgrad_pl(g,ADM_kmin,l,d) = Pgrad_pl(g,ADM_kmin+1,l,d)
             enddo
          endif

          do g = 1, ADM_gall_pl
             Pgrad_pl(g,ADM_kmin-1,l,d) = 0.0_RP
             Pgrad_pl(g,ADM_kmax+1,l,d) = 0.0_RP
          enddo
       enddo
       enddo
    endif

    !--- horizontalize
    call OPRT_horizontalize_vec( Pgrad(:,:,:,GRD_XDIR), Pgrad_pl(:,:,:,GRD_XDIR), & ! [INOUT]
                                 Pgrad(:,:,:,GRD_YDIR), Pgrad_pl(:,:,:,GRD_YDIR), & ! [INOUT]
                                 Pgrad(:,:,:,GRD_ZDIR), Pgrad_pl(:,:,:,GRD_ZDIR)  ) ! [INOUT]



    !---< vertical gradient (half level) >---

    if ( gradtype == I_SRC_default ) then

       !$omp parallel default(none),private(i,j,k,l,d,g), &
       !$omp shared(iall,jall,kmin,kmax,lall,Pgradw,P,VMTR_GAM2H,VMTR_RGSGAM2,GRD_rdgzh)
       do l = 1, lall
          !$omp do schedule(static)
          do k = kmin+1, kmax
             g = 1
             do j = 1, jall
             do i = 1, iall
                Pgradw(g,k,l) = VMTR_GAM2H(i,j,k,l) * ( P(g,k  ,l) * VMTR_RGSGAM2(i,j,k  ,l) &
                                                      - P(g,k-1,l) * VMTR_RGSGAM2(i,j,k-1,l) &
                                                      ) * GRD_rdgzh(k)
                g = g + 1
             enddo
             enddo
          enddo
          !$omp end do

          !$omp workshare
          Pgradw(:,kmin-1,l) = 0.0_RP
          Pgradw(:,kmin  ,l) = 0.0_RP
          Pgradw(:,kmax+1,l) = 0.0_RP
          !$omp end workshare
       enddo
       !$omp end parallel

       if ( ADM_have_pl ) then
          do l = 1, ADM_lall_pl
             do k = ADM_kmin+1, ADM_kmax
             do g = 1, ADM_gall_pl
                Pgradw_pl(g,k,l) = VMTR_GAM2H_pl(g,k,l) * ( P_pl(g,k  ,l) * VMTR_RGSGAM2_pl(g,k  ,l) &
                                                          - P_pl(g,k-1,l) * VMTR_RGSGAM2_pl(g,k-1,l) &
                                                          ) * GRD_rdgzh(k)
             enddo
             enddo

             do g = 1, ADM_gall_pl
                Pgradw_pl(g,ADM_kmin-1,l) = 0.0_RP
                Pgradw_pl(g,ADM_kmin  ,l) = 0.0_RP
                Pgradw_pl(g,ADM_kmax+1,l) = 0.0_RP
             enddo
          enddo
       endif

    elseif( gradtype == I_SRC_horizontal ) then

       !$omp parallel workshare
       Pgradw(:,:,:) = 0.0_RP
       !$omp end parallel workshare
       if ( ADM_have_pl ) then
          Pgradw_pl(:,:,:) = 0.0_RP
       endif

    endif

    call PROF_rapend('____src_pres_gradient',2)

    return
  end subroutine src_pres_gradient

  !-----------------------------------------------------------------------------
  !> Calculation of buoyacy force
  !> NOTICE : Upward direction is positive for buoiw.
  subroutine src_buoyancy( &
       rhog,  rhog_pl, &
       buoiw, buoiw_pl )
    use scale_const, only: &
       GRAV => CONST_GRAV
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_iall,    &
       ADM_jall,    &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_vmtr, only: &
       VMTR_C2Wfact,    &
       VMTR_C2Wfact_pl
    implicit none

    real(RP), intent(in)  :: rhog    (ADM_gall   ,ADM_kall,ADM_lall   ) ! density perturbation ( G^1/2 x gam2 )
    real(RP), intent(in)  :: rhog_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: buoiw   (ADM_gall   ,ADM_kall,ADM_lall   ) ! buoyancy force  at half level
    real(RP), intent(out) :: buoiw_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: iall, jall, kmin, kmax, lall

    integer :: i, j, k, l, g
    !---------------------------------------------------------------------------

    call PROF_rapstart('____src_buoyancy',2)

    iall = ADM_iall
    jall = ADM_jall
    kmin = ADM_kmin
    kmax = ADM_kmax
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l,g), &
    !$omp shared(iall,jall,kmin,kmax,lall,buoiw,rhog,VMTR_C2Wfact,GRAV)
    do l = 1, lall
       !$omp do schedule(static)
       do k = kmin+1, kmax
          g = 1
          do j = 1, jall
          do i = 1, iall
             buoiw(g,k,l) = -GRAV * ( VMTR_C2Wfact(i,j,k,1,l) * rhog(g,k  ,l) &
                                    + VMTR_C2Wfact(i,j,k,2,l) * rhog(g,k-1,l) )
             g = g + 1
          enddo
          enddo
       enddo
       !$omp end do

       !$omp workshare
       buoiw(:,kmin-1,l) = 0.0_RP
       buoiw(:,kmin  ,l) = 0.0_RP
       buoiw(:,kmax+1,l) = 0.0_RP
       !$omp end workshare
    enddo
    !$omp end parallel

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             buoiw_pl(g,k,l) = -GRAV * ( VMTR_C2Wfact_pl(g,k,1,l) * rhog_pl(g,k  ,l) &
                                       + VMTR_C2Wfact_pl(g,k,2,l) * rhog_pl(g,k-1,l) )
          enddo
          enddo

          do g = 1, ADM_gall_pl
             buoiw_pl(g,ADM_kmin-1,l) = 0.0_RP
             buoiw_pl(g,ADM_kmin  ,l) = 0.0_RP
             buoiw_pl(g,ADM_kmax+1,l) = 0.0_RP
          enddo
       enddo
    endif

    call PROF_rapend('____src_buoyancy',2)

    return
  end subroutine src_buoyancy

end module mod_src
