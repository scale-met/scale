!-------------------------------------------------------------------------------
!
!+  precipitation transport module
!
!-------------------------------------------------------------------------------
module mod_precip_transport
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       This module contains subroutines w.r.t. precipitation transport
  !       
  ! 
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: 
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.01     2011/10/24   T.Seiki, Imported from NICAM
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  public :: precip_transport_nwater
  !-----------------------------------------------------------------------------
  private
  character(len=7), parameter ::  PRCIP_TRN_ECORRECT = "KIN2KIN"
  real(8), parameter :: ADM_VMISS = 1.d0
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine precip_transport_nwater ( &
       ijdim,                 & !--- IN :
       kdim, kmin, kmax,      & !--- IN :
       rhog,                  & !--- INOUT :
       rhogvx,                & !--- INOUT :
       rhogvy,                & !--- INOUT :
       rhogw,                 & !--- INOUT :
       rhoge,                 & !--- INOUT :
       rhogq,                 & !--- INOUT :
       rho,                   & !--- INOUT :
       tem,                   & !--- INOUT :
       pre,                   & !--- INOUT :
       q,                     & !--- INOUT :
       qd,                    & !--- OUT :
       z,                     & !--- IN : height at the cell-ceter
       zh,                    & !--- IN : height at the 
       dz,                    & !--- IN : height interval
       dzh,                   & !--- IN : height interval
       WPREC,                 & !--- IN :
       precipitating_flag,    & !--- IN :
       precip,                & !--- OUT :
       precip_rhoe,           & !--- OUT :
       precip_lh_heat,        & !--- OUT :
       precip_rhophi,         & !--- OUT :
       precip_rhokin,         & !--- OUT :
       gsgam2,                & !--- IN :
       gsgam2h,               & !--- IN :
       rgs,                   & !--- IN :
       rgsh,                  & !--- IN :
       dt                     & !--- IN :
       )
    !
    use mod_atmos_cnst, only : &
         nqmax => TRC_VMAX, &
         NQW_STR, NQW_END,  &
         NQL_STR, NQL_END,  &
         NQS_STR, NQS_END,  &
         CVW,               &
         LHF,               &
         CNST_CP,           &
         CNST_CV,           &
         CNST_EGRAV,        &
         CNST_RAIR,         &
         CNST_RVAP
    use mod_thrmdyn, only : &
         thrmdyn_qd,        &
         thrmdyn_tempre
    !
    implicit none
    !
    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: kmin
    integer, intent(in)    :: kmax
    real(8), intent(inout) :: rhog(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhoge(1:ijdim,1:kdim)
    real(8), intent(inout) :: rhogq(1:ijdim,1:kdim,1:nqmax)
    !
    real(8), intent(inout) :: rho(1:ijdim,1:kdim)
    real(8), intent(inout) :: tem(1:ijdim,1:kdim)
    real(8), intent(inout) :: pre(1:ijdim,1:kdim)
    real(8), intent(inout) :: q(1:ijdim,1:kdim,1:nqmax)
    real(8), intent(out)   :: qd(1:ijdim,1:kdim)
    !
    real(8), intent(in) :: WPREC(1:ijdim,1:kdim,1:nqmax)
    logical, intent(in) :: precipitating_flag(1:nqmax)
    real(8), intent(in) :: z (1:kdim) ! height
    real(8), intent(in) :: zh(1:kdim) ! height
    real(8), intent(in) :: dz(1:kdim) ! height
    real(8), intent(in) :: dzh(1:kdim)! height
    !
    real(8), intent(out) :: precip(1:ijdim,2) ! rain or snow
    real(8), intent(out) :: precip_rhoe(1:ijdim)
    real(8), intent(out) :: precip_lh_heat(1:ijdim)
    real(8), intent(out) :: precip_rhophi(1:ijdim)
    real(8), intent(out) :: precip_rhokin(1:ijdim)
    !
    real(8), intent(in) :: gsgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsgam2h(1:ijdim,1:kdim)
    real(8), intent(in) :: rgs(1:ijdim,1:kdim)
    real(8), intent(in) :: rgsh(1:ijdim,1:kdim)
    !
    real(8), intent(in) :: dt
    !
    real(8) :: rdt
    !
    real(8) :: rdz(kdim)
    real(8) :: rdzh(kdim)
    real(8) :: GRD_afac(kdim)
    real(8) :: GRD_bfac(kdim)
    real(8) :: GRD_cfac(kdim)
    real(8) :: GRD_dfac(kdim)
    !
    ! local
    real(8) :: fprec(1:ijdim,1:kdim,1:nqmax)
    real(8) :: drho_prec(1:ijdim,1:kdim,1:nqmax)
    !
    real(8) :: rhoeq_prec(1:ijdim,1:kdim)
    real(8) :: fprec_rhoe(1:ijdim,1:kdim,1:nqmax)
    real(8) :: drhoe_prec(1:ijdim,1:kdim)
    real(8) :: rhophiq_prec(1:ijdim,1:kdim)
    real(8) :: rhokin_h_prec(1:ijdim,1:kdim)
    real(8) :: rhokin_v_prec(1:ijdim,1:kdim)
    !
    real(8) :: fprec_rhophi(1:ijdim,1:kdim,1:nqmax)
    real(8) :: fprec_rhokin_h(1:ijdim,1:kdim,1:nqmax)
    real(8) :: fprec_rhokin_v(1:ijdim,1:kdim,1:nqmax)
    real(8) :: drhophi_prec(1:ijdim,1:kdim)
    real(8) :: drhokin_h_prec(1:ijdim,1:kdim)
    real(8) :: drhokin_v_prec(1:ijdim,1:kdim)
    !
    real(8) :: rhouq(1:ijdim,1:kdim)
    real(8) :: rhovq(1:ijdim,1:kdim)
    real(8) :: rhowq(1:ijdim,1:kdim)
    real(8) :: fprec_rhou(1:ijdim,1:kdim)
    real(8) :: fprec_rhov(1:ijdim,1:kdim)
    real(8) :: fprec_rhow(1:ijdim,1:kdim)
    !
    real(8) :: rhogkin(1:ijdim,1:kdim)
    real(8) :: rhogkin_h(1:ijdim,1:kdim)
    real(8) :: rhogkin_v(1:ijdim,1:kdim)
    real(8) :: rhog_h(1:ijdim,1:kdim)
    real(8) :: kin_h0(1:ijdim,1:kdim),kin_h(1:ijdim,1:kdim)
    real(8) :: kin_v0(1:ijdim,1:kdim),kin_v(1:ijdim,1:kdim)
    real(8) :: vx_t(1:ijdim,1:kdim),vy_t(1:ijdim,1:kdim),w_t(1:ijdim,1:kdim)
    !
    real(8) :: wq(1:ijdim,1:kdim), ein(1:ijdim,1:kdim)
    real(8) :: tmp(1:ijdim,1:kdim)
    real(8) :: tmp2(1:ijdim,1:kdim)
    !
    integer :: k, nq, ij
    !

#ifdef _FPCOLL_
call START_COLLECTION("precip_transport_nwater")
#endif

    rdt      = 1.d0/dt
    rdz(:)   = 1.d0/dz(:)
    rdzh(:)  = 1.d0/dzh(:)
    GRD_afac(:) = 1.d0
    GRD_bfac(:) = 1.d0
    GRD_cfac(:) = 1.d0
    GRD_dfac(:) = 1.d0
    !
    call cnvvar_rhokin(&
         ijdim,        & !--- IN : number of horizontal grid
         kdim, kmin, kmax,   & !--- IN : number of vertical grid
         rhog,         &  !--- IN : rho     ( gam2 X G^{1/2}  )
         rhogvx,       &  !--- IN : rho*Vx  ( gam2 X G^{1/2}  )
         rhogvy,       &  !--- IN : rho*Vy  ( gam2 X G^{1/2}  )
         rhogw,        &  !--- IN : rho*w   ( gam2 X G^{1/2}  )
         gsgam2,       &  !--- IN : G^{1/2} at the cell center
         gsgam2h,      &  !--- IN : G^{1/2} at the cell wall
         rhogkin,      &  !--- OUT : 1/2 rho*v^2  ( gam2 X G^{1/2}  )
         rhogkin_h,    &
         rhogkin_v     &
         )
    rhogkin=0.d0
    !
    do nq = 1, nqmax
       !
       if(.not.precipitating_flag(nq)) then
          fprec(:,:,nq) = 0.0D0
          fprec_rhoe(:,:,nq) = 0.0D0
          fprec_rhophi(:,:,nq) = 0.0D0
          fprec_rhokin_h(:,:,nq) = 0.0D0
          fprec_rhokin_v(:,:,nq) = 0.0D0
          cycle
       end if
       !----- mass
       call vadv1d_getflux(             &
            ijdim,                  & !--- in
            kdim,                   &
            kmin,                   &
            kmax,                   &
            dz(:),                  &
            zh(:),                  &
            rhogq(:,:,nq)*rgs(:,:), & !--- in
            WPREC(:,:,nq),          & !--- in
            fprec(:,:,nq),          & !--- out
            dt )                      !--- in
       ! hereafter nq should be active tracers.
       if( nq >= NQW_STR .and. nq <= NQW_END )then
          !--- internal energy
          rhoeq_prec(:,:) = rhog(:,:) &
               * q(:,:,nq) * CVW(nq) * tem(:,:) * rgs(:,:)
          call vadv1d_getflux( &
               ijdim,          & !--- in
               kdim,           &
               kmin,           &
               kmax,           &
               dz(:),          &
               zh(:),          &
               rhoeq_prec(:,:),   & !--- in
               WPREC(:,:,nq),     & !--- in
               fprec_rhoe(:,:,nq),& !--- out
               dt )
          !--- potential energy
          do k=1, kdim
             rhophiq_prec(:,k) = rhog(:,k) &
                  * q(:,k,nq) * CNST_EGRAV * z(k) * rgs(:,k)
          end do
          call vadv1d_getflux(    &
               ijdim,             & !--- in
               kdim,              &
               kmin,              &
               kmax,              &
               dz(:),             &
               zh(:),             &
               rhophiq_prec(:,:), & !--- in
               WPREC(:,:,nq),     & !--- in
               fprec_rhophi(:,:,nq), & !--- out
               dt )
          !--- horizontal kinetic energy
          rhokin_h_prec(:,:) = q(:,:,nq) * rhogkin_h(:,:) * rgs(:,:)
          call vadv1d_getflux(       &
               ijdim,                & !--- in
               kdim,                 &
               kmin,                 &
               kmax,                 &
               dz(:),                &
               zh(:),                &
               rhokin_h_prec(:,:),   & !--- in
               WPREC(:,:,nq),        & !--- in
               fprec_rhokin_h(:,:,nq), & !--- out
               dt )
          !--- vertical kinetic energy
          rhokin_v_prec(:,:) = 0.0D0
          do k = kmin+1,kmax
             rhokin_v_prec(:,k) = 0.5D0* (q(:,k,nq)+q(:,k-1,nq)) * rhogkin_v(:,k) * rgsh(:,k)
          end do
          wq(:,:) = 0.0D0
          do k = kmin+1,kmax-1
             wq(:,k) = 0.5D0*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
          end do
          !<--- half -> full
          fprec_rhokin_v(:,:,nq) = 0.0D0
          call vadv1d_getflux(               &
               ijdim,                        & !--- in
               kdim-1,                       &
               kmin,                         &
               kmax-1,                       &
               dzh(kmin:kmax+1),             &
               z(kmin-1:kmax),               &
               rhokin_v_prec(:,kmin:kmax+1), & !--- in
               wq(:,kmin:kmax+1),            & !--- in
               fprec_rhokin_v(:,kmin-1:kmax,nq),& !--- out
               dt )
       end if
    end do
    !
    !--- tendency of density and energy
    do nq=1, nqmax
       drho_prec(:,kmin-1,nq) = 0.d0
       drho_prec(:,kmax+1,nq) = 0.d0
    end do
    drhoe_prec(:,kmin-1) = 0.d0
    drhoe_prec(:,kmax+1) = 0.d0
    drhophi_prec(:,kmin-1) = 0.d0
    drhophi_prec(:,kmax+1) = 0.d0
    !
    do nq=1, nqmax
       do k = kmin, kmax
          drho_prec(:,k,nq) &
               = - ( fprec(:,k+1,nq) - fprec(:,k,nq) ) * rdz(k)
       end do
    end do
    do k= kmin, kmax
       drhoe_prec(:,k) = 0.0D0
       drhophi_prec(:,k) = 0.0D0
       drhokin_h_prec(:,k) = 0.0D0
       do nq=NQW_STR,NQW_END
          drhoe_prec(:,k) = drhoe_prec(:,k)&
               - ( fprec_rhoe(:,k+1,nq) - fprec_rhoe(:,k,nq) )*rdz(k)
          drhophi_prec(:,k) = drhophi_prec(:,k)&
               - ( fprec_rhophi(:,k+1,nq) - fprec_rhophi(:,k,nq) )*rdz(k)
          drhokin_h_prec(:,k) = drhokin_h_prec(:,k)&
               - ( fprec_rhokin_h(:,k+1,nq) - fprec_rhokin_h(:,k,nq) )*rdz(k)
       end do
    end do
    !
    drhokin_h_prec(:,kmax+1) = 0.0D0
    drhokin_h_prec(:,kmin-1) = 0.0D0
    !
    drhokin_v_prec(:,kmin)=0.0D0
    drhokin_v_prec(:,kmin-1)=0.0D0
    drhokin_v_prec(:,kmax+1)=0.0D0
    do k = kmin+1, kmax
       drhokin_v_prec(:,k) = 0.0D0
       do nq=NQW_STR,NQW_END
          drhokin_v_prec(:,k) = drhokin_v_prec(:,k)&
               - ( fprec_rhokin_v(:,k,nq) - fprec_rhokin_v(:,k-1,nq) )*rdzh(k)
       end do
    end do
    !
    !--- momentum transport: new values due to rain
    !<------ use old values of rho and qr for consistency
    do nq=NQW_STR,NQW_END
       if(.not.precipitating_flag(nq)) then
          cycle
       end if
       do k = 1,kdim 
          rhouq(:,k)              &
               = ( rhogvx(:,k)    &
               +   rhogvy(:,k)  ) &
               *   q(:,k,nq)*rgs(:,k)
          rhovq(:,k)             &
               = ( rhogvx(:,k)   &
               +   rhogvy(:,k) ) &
               *   q(:,k,nq)*rgs(:,k)
       end do
       !
       call vadv1d_getflux(       &
            ijdim,            & !--- in
            kdim,             &
            kmin,             &
            kmax,             &
            dz(:),            &
            zh(:),            &
            rhouq(:,:),       & !--- in
            WPREC(:,:,nq),    & !--- in
            fprec_rhou(:,:),  & !--- out
            dt                & !--- in
            )
       call vadv1d_getflux(  &
            ijdim,               & !--- in
            kdim,             &
            kmin,             &
            kmax,             &
            dz(:),            &
            zh(:),            &
            rhovq(:,:),       & !--- in
            WPREC(:,:,nq),    & !--- in
            fprec_rhov(:,:),  & !--- out
            dt                   & !--- in
            )
       !
       !--- update rhogvx, rhogvy
       do k = kmin, kmax
          rhogvx(:,k) = rhogvx(:,k) - ( fprec_rhou(:,k+1) - fprec_rhou(:,k) ) * rdz(k)
          rhogvy(:,k) = rhogvy(:,k) - ( fprec_rhov(:,k+1) - fprec_rhov(:,k) ) * rdz(k)
       end do
       !=====================
       rhowq(:,:) = 0.0D0
       do k = kmin+1,kmax
          rhowq(:,k) &
               = rhogw(:,k)* 0.5D0* (q(:,k,nq)+q(:,k-1,nq)) *rgsh(:,k)
       end do
       wq(:,:) = 0.0D0
       do k = kmin+1,kmax-1
          wq(:,k) = 0.5D0*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
       end do
       call vadv1d_getflux(            &
            ijdim,                     &
            kdim-1,                    &
            kmin,                      &
            kmax-1,                    &
            dzh(kmin:kmax+1),          &
            z(kmin-1:kmax),            &
            rhowq(:,kmin:kmax+1),      &
            wq(:,kmin:kmax+1),         &
            fprec_rhow(:,kmin-1:kmax), &
            dt )
       !
       !--- update rhogw
       do k = kmin+1, kmax
          rhogw(:,k) = rhogw(:,k) - ( fprec_rhow(:,k) - fprec_rhow(:,k-1) ) &
               * rdzh(k)
       end do
    end do
    !
    ! new values due to rain:
    ! Change in internal energy comes from precipitation of rain and dissipation
    ! of kinetic energy due to drag force. 
    ! See Ooyama(2001) (3.13)
    !
    !--- update rhogqr
    do nq=1, nqmax
       rhogq(:,:,nq) = rhogq(:,:,nq) + drho_prec(:,:,nq)
    end do
    !
    !--- update rhog, rho
    do nq=NQW_STR,NQW_END
       rhog(:,:) = rhog(:,:) + drho_prec(:,:,nq)
    end do
    rho(:,:) = rhog(:,:)/gsgam2(:,:)
    !--- update rhoge, rhogkin
    rhoge(:,:) = rhoge(:,:) + drhoe_prec(:,:) + drhophi_prec(:,:)
    do nq=NQW_STR,NQW_END
       do k=1, kdim
       rhoge(:,k) = rhoge(:,k) &
            - CNST_EGRAV * z(k) * ( drho_prec(:,k,nq) )
       end do
    end do
    !
    rhogkin_h(:,:) = rhogkin_h(:,:)+  drhokin_h_prec(:,:)
    rhogkin_v(:,:) = rhogkin_v(:,:)+  drhokin_v_prec(:,:)
    !
    !--- update q
    do nq = 1, nqmax
       q(:,:,nq) = rhogq(:,:,nq) / rhog(:,:)
    end do
    !--- update qd
    call thrmdyn_qd( &
         qd,         & !--- out
         q           & !--- in
         )
    if(trim(PRCIP_TRN_ECORRECT)=='KIN2EIN') then
       tmp2(:,kmin-1) = 0.0D0
       tmp2(:,kmax+1) = 0.0D0
       do k = kmin, kmax
          tmp2(:,k) = rhogkin_h(:,k)&
               + (GRD_dfac(k  )*rhogkin_v(:,k+1)&
               +  GRD_cfac(k  )*rhogkin_v(:,k  ) )*0.5D0
       end do
       !
       call cnvvar_rhokin(&
            ijdim,        &  !--- IN : number of horizontal grid
            kdim, kmin, kmax,   & !--- IN : number of vertical grid
            rhog,         &  !--- IN : rho     ( gam2 X G^{1/2}  )
            rhogvx,       &  !--- IN : rho*Vx  ( gam2 X G^{1/2}  )
            rhogvy,       &  !--- IN : rho*Vy  ( gam2 X G^{1/2}  )
            rhogw,        &  !--- IN : rho*w   ( gam2 X G^{1/2}  )
            gsgam2,       &  !--- IN : G^{1/2} at the cell center
            gsgam2h,      &  !--- IN : G^{1/2} at the cell wall
            tmp           &  !--- OUT : 1/2 rho*v^2  ( gam2 X G^{1/2}  )
            )
       rhoge(:,:) = rhoge(:,:) + ( tmp2(:,:) - tmp(:,:) )
       !
    else if(trim(PRCIP_TRN_ECORRECT)=='KIN2KIN') then
       !
       !--- fix kinetic energy
       !
       ! Modify C.Kodama 09.07.14 =>
       kin_h0(:,kmin:kmax) = rhogkin_h(:,kmin:kmax)/rhog(:,kmin:kmax)
       vx_t(:,kmin:kmax) = rhogvx(:,kmin:kmax)/rhog(:,kmin:kmax)
       vy_t(:,kmin:kmax) = rhogvy(:,kmin:kmax)/rhog(:,kmin:kmax)
       kin_h(:,kmin:kmax) = ((vx_t(:,kmin:kmax)**2)+(vy_t(:,kmin:kmax)**2))*0.5D0
       where(kin_h(:,kmin:kmax)>1.0D-20)
          vx_t(:,kmin:kmax) = vx_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
          vy_t(:,kmin:kmax) = vy_t(:,kmin:kmax) *  sqrt(abs(kin_h0(:,kmin:kmax)/kin_h(:,kmin:kmax)))
       end where
       rhogvx(:,kmin:kmax) = rhog(:,kmin:kmax) * vx_t(:,kmin:kmax)
       rhogvy(:,kmin:kmax) = rhog(:,kmin:kmax) * vy_t(:,kmin:kmax)
       !
       rhog_h(:,kmin:kmax+1) = 1.0D0
       do k =kmin,kmax+1
          rhog_h(:,k) = 0.5D0 *            &
               ( GRD_afac(k) * rhog(:,k  )/gsgam2(:,k  ) &
               + GRD_bfac(k) * rhog(:,k-1)/gsgam2(:,k-1) &
               ) * gsgam2h(:,k)
       end do
       kin_v0(:,kmin:kmax+1) = rhogkin_v(:,kmin:kmax+1)/rhog_h(:,kmin:kmax+1)
       w_t(:,kmin:kmax+1) = rhogw(:,kmin:kmax+1)/rhog_h(:,kmin:kmax+1)
       kin_v(:,kmin:kmax+1) = (w_t(:,kmin:kmax+1)**2) * 0.5D0
       where(kin_v(:,kmin:kmax+1)>1.0D-20)
          w_t(:,kmin:kmax+1) = w_t(:,kmin:kmax+1) *  sqrt(abs(kin_v0(:,kmin:kmax+1)/kin_v(:,kmin:kmax+1)))
       end where
       rhogw(:,kmin:kmax+1) = rhog_h(:,kmin:kmax+1) * w_t(:,kmin:kmax+1)
       !<= Modify C.Kodama 09.07.14
    else
       write(*,*) 'Error in precip_transport'
    end if
    !
    !--- update temerature & pressure
    ein(:,:) =  rhoge(:,:) / rhog(:,:)
    call thrmdyn_tempre( &
         tem,                  &  !--- OUT : temperature              
         pre,                  &  !--- OUT : pressure
         ein,                  &  !--- IN  : internal energy
         rho,                  &  !--- IN  : density
         qd,                   &  !--- IN  : dry concentration 
         q )                      !--- IN  : water concentration 
    !
    precip(:,:) = 0.0D0
    do nq=NQL_STR,NQL_END
       precip(:,1) = precip(:,1) -fprec(:,kmin,nq)*rdt
    end do
    do nq=NQS_STR,NQS_END
       precip(:,2) = precip(:,2) -fprec(:,kmin,nq)*rdt
    end do
    !
    precip_rhoe(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhoe(:) = precip_rhoe(:) - fprec_rhoe(:,kmin,nq)*rdt
    end do
    precip_lh_heat(:) = 0.0D0
    do nq=NQS_STR,NQS_END
       precip_lh_heat(:) = precip_lh_heat(:) + fprec(:,kmin,nq)*LHF*rdt
    end do
    precip_rhophi(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhophi(:) = precip_rhophi(:) - fprec_rhophi(:,kmin,nq)*rdt
    end do
    precip_rhokin(:) = 0.0D0
    do nq=NQW_STR,NQW_END
       precip_rhokin(:) = precip_rhokin(:) &
            - fprec_rhokin_h(:,kmin,nq)*rdt &
            - fprec_rhokin_v(:,kmin,nq)*rdt
    end do
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("precip_transport_nwater")
#endif

    return
    !
  end subroutine precip_transport_nwater
  !
  subroutine cnvvar_rhokin(&
       ijdim,              & !--- IN : number of horizontal grid
       kdim, kmin, kmax,   & !--- IN : number of vertical grid
       rhog,               &  !--- IN : rho     ( gam2 X G^{1/2}  )
       rhogvx,             &  !--- IN : rho*Vx  ( gam2 X G^{1/2}  )
       rhogvy,             &  !--- IN : rho*Vy  ( gam2 X G^{1/2}  )
       rhogw,              &  !--- IN : rho*w   ( gam2 X G^{1/2}  )
       gsqrtgam2,          &  !--- IN : G^{1/2} at the cell center
       gsqrtgam2h,         &  !--- IN : G^{1/2} at the cell wall
       rhogkin,            &  !--- OUT : 1/2 rho*v^2  ( gam2 X G^{1/2}  )
       rhogkin_h,          &
       rhogkin_v           &
       )
    !------ 
    !------ Calculation of kinetic energy.
    !------ 
    !
    !
    implicit none
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    real(8), intent(in) :: rhog(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvx(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogvy(1:ijdim,1:kdim)
    real(8), intent(in) :: rhogw(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2(1:ijdim,1:kdim)
    real(8), intent(in) :: gsqrtgam2h(1:ijdim,1:kdim)
    real(8), intent(out) :: rhogkin(1:ijdim,1:kdim)
    !
    real(8) :: rhog_h(1:ijdim,1:kdim)
    real(8), intent(out), optional :: rhogkin_h(1:ijdim,1:kdim)
    real(8), intent(out), optional :: rhogkin_v(1:ijdim,1:kdim)
    real(8) :: rhogkin_h0(1:ijdim,1:kdim)
    real(8) :: rhogkin_v0(1:ijdim,1:kdim)
    !
    real(8) :: GRD_afac(kdim)
    real(8) :: GRD_bfac(kdim)
    real(8) :: GRD_cfac(kdim)
    real(8) :: GRD_dfac(kdim)
    !
    integer :: k
    !

#ifdef _FPCOLL_
call START_COLLECTION("cnvvar_rhokin")
#endif

    GRD_afac=1.d0
    GRD_bfac=1.d0
    GRD_cfac=1.d0
    GRD_dfac=1.d0
    !
    !--- rhogkin = gamma^2 * g_sqrt * rho * kin
    !
    !--- horizontal kinetic energy
    rhogkin_h0(:,kmin-1) = ADM_VMISS
    rhogkin_h0(:,kmax+1) = ADM_VMISS
    do k = kmin, kmax
       rhogkin_h0(:,k)  &
            =((rhogvx(:,k)/rhog(:,k))**2 &
            + (rhogvy(:,k)/rhog(:,k))**2 )&
            * rhog(:,k)*0.5D0
    end do
    !
    !--- rhog at the half level
    do k =kmin,kmax+1
       rhog_h(:,k) = 0.5D0 *            &
            ( GRD_afac(k) * rhog(:,k  )/gsqrtgam2(:,k  ) &
            + GRD_bfac(k) * rhog(:,k-1)/gsqrtgam2(:,k-1) &
            ) * gsqrtgam2h(:,k)
    end do
    !
    !--- vertical kinetic energy
    rhogkin_v0(:,kmin-1) = ADM_VMISS
    rhogkin_v0(:,kmin) = 0.0D0
    rhogkin_v0(:,kmax+1) = 0.0D0
    do k = kmin+1, kmax
       rhogkin_v0(:,k) &
            =((rhogw(:,k)/rhog_h(:,k))**2) &
            * rhog_h(:,k)*0.5D0
    end do
    !
    rhogkin(:,kmin-1) = ADM_VMISS
    rhogkin(:,kmax+1) = ADM_VMISS
    do k = kmin, kmax
       rhogkin(:,k) = rhogkin_h0(:,k)&
            + (GRD_dfac(k  )*rhogkin_v0(:,k+1)&
            +  GRD_cfac(k  )*rhogkin_v0(:,k  ) )*0.5D0
    end do
    !
    if(present(rhogkin_h)) then
       rhogkin_h = rhogkin_h0
    end if
    if(present(rhogkin_v)) then
       rhogkin_v = rhogkin_v0
    end if

#ifdef _FPCOLL_
call STOP_COLLECTION("cnvvar_rhokin")
#endif

    !
    return
  end subroutine cnvvar_rhokin
  !
  subroutine vadv1d_getflux( &
       ijdim,            & !--- IN : number of horizontal grids
       kdim,             & !--- IN
       kmin,             & !--- IN
       kmax,             & !--- IN
       dz,               & !--- IN
       zh,               & !--- IN
       rhof,             & !--- IN : rho X f X gam2
       wp,               & !--- IN : verical velocity at the full level
       frhof,            & !--- OUT : vertical flux ( X dt )
       dt                & !--- IN : time interval
       )
    !
    implicit none
    !
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    real(8), intent(in) :: dz(1:kdim)
    real(8), intent(in) :: zh(1:kdim)
    real(8), intent(in) :: rhof(1:ijdim,1:kdim)   ! p-grid
    real(8), intent(in) :: wp(1:ijdim,1:kdim)     ! p-grid
    real(8), intent(out) :: frhof(1:ijdim,1:kdim) ! w-grid
    real(8), intent(in) :: dt
    !
    real(8) :: wh(1:ijdim,1:kdim)          ! w-grid
    real(8) :: rhofh_cell(1:ijdim,1:kdim)  ! w-grid
    !
    integer :: ij, k
    integer :: k2
    integer :: kcell(1:ijdim,1:kdim) ! w-grid
    real(8) :: zdis(1:ijdim,1:kdim)  ! w-grid
    !
    integer :: kcell_max(1:kdim)    ! w-grid
    integer :: kcell_min(1:kdim)    ! w-grid
    real(8) :: zzmax, zzmin
    !
    real(8) :: rhofh(1:ijdim,1:kdim)

#ifdef _FPCOLL_
call START_COLLECTION("vadv1d_getflux")
#endif

    !
    !--- vetical velocity at the half level
    do k = kmin+1, kmax
       wh(:,k) = ( wp(:,k-1) + wp(:,k) ) / 2.d0
    end do
    !------ bottom boundary for wh
    wh(:,kmin) = wp(:,kmin)
    wh(:,kmin-1) = wp(:,kmin-1)
    !------ topboundary for wh : same as inner region
    wh(:,kmax+1) = wp(:,kmax)
    !
    !--- calculation of distance of cell wall during dt
    do k = kmin+1, kmax
       zdis(:,k) = dt * wh(:,k) &
            - dt ** 2 * wh(:,k) &
            * ( wh(:,k+1) - wh(:,k-1) ) / ( dz(k-1) + dz(k) ) / 2.0d0 &
            + dt ** 3 * wh(:,k) * &
            ( ( ( wh(:,k+1) - wh(:,k-1) ) / ( dz(k-1) + dz(k) ) ) ** 2 &
            + wh(:,k) * ( ( ( wh(:,k+1) - wh(:,k) ) / dz(k) &
            - ( wh(:,k) - wh(:,k-1) ) / dz(k-1) ) &
            / ( dz(k-1) + dz(k) ) * 2.d0 ) ) / 6.0d0
    end do
    !--- bottom and top boundary for zdis
    k=kmin-1
    zdis(:,k) = 0.d0 
    k = kmin
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k+1) - wh(:,k) ) / dz(k) / 2.0d0
    k = kmax+1
    zdis(:,k) = dt * wh(:,k) &
         - dt ** 2 * wh(:,k) * ( wh(:,k) - wh(:,k-1) ) / dz(k-1) / 2.0d0
    !
    !--- calculation of kcell
    !------  top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    do k = kmin, kmax+1
       kcell(:,k) = k
    end do
    !
    !------ setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval(zdis(:,k))
       zzmin = minval(zdis(:,k))
       kcell_min(k) = k
       kcell_max(k) = k
       if ( zzmax > 0.d0 ) then
          do k2 = k, kmin, -1
             if ( zh(k2) <= zh(k) - zzmax &
                  .and. zh(k) - zzmax < zh(k2+1) ) then
                kcell_min(k) = k2
                exit
             end if
          end do
       end if
       if ( zzmin < 0.d0 ) then
          do k2 = k, kmax
             if ( zh(k2) <= zh(k) - zzmin &
                  .and. zh(k) - zzmin < zh(k2+1) ) then
                kcell_max(k) = k2
                exit
             end if
          end do
       end if
    end do
    !
    !------ determine the kcell at each point.
    do k = kmin, kmax
       if ( kcell_min(k) == k .and. kcell_max(k) == k ) then
          kcell(:,k) = k
       else
          kcell(:,k) = 0
          do k2 = kcell_min(k), kcell_max(k)
             kcell(:,k) = max ( kcell(:,k), &
                  int ( k2 &
                  * sign ( 1.d0, ( zh(k) - zdis(:,k) ) - zh(k2) ) &
                  * sign ( 1.d0, zh(k2+1) - ( zh(k) - zdis(:,k) ) ) ) )
          end do
       end if
    end do
    !
    !--- set zero value for norain, see if block in do loop
    frhof(:,:) = 0.d0 
    !
    !------ integration in the integer cells
    do k = kmin, kmax
       if ( kcell_min(k) == k .and. kcell_max(k) == k ) then
          ! frhof(:,:) = 0.d0
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             frhof(:,k) &
                  = frhof(:,k) &
                  + rhof(:,k2) * dz(k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.d0 ) / 2.d0 &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.d0 ) / 2.d0 &
                  + ( - 1.d0 ) * ( sign ( 1, k2 - k ) + 1.d0 ) / 2.d0 &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.d0 ) / 2.d0 )
             zdis(:,k) = zdis(:,k) - dz(k2) &
                  * ( ( sign ( 1, k2 - ( kcell(:,k) + 1 ) ) &
                  + 1.d0 ) / 2.d0 &
                  * ( sign ( 1, ( k - 1 ) - k2 ) + 1.d0 ) / 2.d0 &
                  + ( - 1.d0 ) * ( sign ( 1, k2 - k ) + 1.d0 ) / 2.d0 &
                  * ( sign ( 1, ( kcell(:,k) - 1 ) - k2 ) &
                  + 1.d0 ) / 2.d0 )
          end do
       end if
    end do
    !
    ! top boundary: rigid
    !  zdis(:,kmax+1) = 0.d0
    !  frhof(:,kmax+1) = 0.d0
    !
    !--- integration in the non-integer cell ( depending on the scheme below )
    call sl0_upwind_new ( ijdim, kdim, kmin, kmax, rhof, kcell, rhofh_cell)
    !
    frhof(:,:) = frhof(:,:) + rhofh_cell(:,:) * zdis(:,:)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("vadv1d_getflux")
#endif

    return
  end subroutine vadv1d_getflux
  !
  subroutine sl0_upwind_new(&
       ijdim,           & !--- IN : number of horizontal grids
       kdim,            & !--- IN
       kmin,            & !--- IN
       kmax,            & !--- IN
       f,               & !--- IN : f
       kcell,           & !--- IN : target cell number
       fh               & !--- OUT : flux at the cell wall
       )
    !------
    !
    !------ constant in a box
    implicit none
    integer, intent(in) :: ijdim
    integer, intent(in) :: kdim
    integer, intent(in) :: kmin
    integer, intent(in) :: kmax
    real(8), intent(in) :: f(1:ijdim,1:kdim)
    integer, intent(in) :: kcell(1:ijdim,1:kdim)
    real(8), intent(out) :: fh(1:ijdim,1:kdim)
    integer :: ij, k

#ifdef _FPCOLL_
call START_COLLECTION("sl0_upwind_new")
#endif

    !
    fh(:,:) = 0.d0
    !<--- including the condition of [fh(:,kmax+1) = 0.d0] 
    !<--- at the top boundary
    do k = kmin, kmax
       do ij = 1, ijdim
          fh(ij,k) = f(ij,kcell(ij,k))
       end do
    end do
    !
    do k = kmin, kmax
       where ( kcell(:,k) == 0 )
          fh(:,k) =  f(:,kmin)
       end where
       where ( kcell(:,k) == kmax+1 )
          fh(:,k) =  f(:,kmax)
       end where
    end do
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("sl0_upwind_new")
#endif

    return
  end subroutine sl0_upwind_new
  !
end module mod_precip_transport
!-------------------------------------------------------------------------------
