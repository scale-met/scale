!-------------------------------------------------------------------------------
!
!+  precipitation transport module
!
!-------------------------------------------------------------------------------
module mod_precip_transport
  !-----------------------------------------------------------------------------
  !
  !++ Description: 
  !       precipitation transport
  !
  !       
  !++ Current Corresponding Author : H.Tomita
  ! 
  !++ History: NDW6
  !      Version   Date       Comment 
  !      -----------------------------------------------------------------------
  !      0.01     2011/10/24   T.Seiki, Imported from NICAM
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
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

  public :: precip_transport_nwater

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

  character(len=7), parameter :: PRCIP_TRN_ECORRECT = "KIN2KIN"

  real(8),          parameter :: ADM_VMISS = 1.D0

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> precipitation transport
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
       dt                     ) !--- IN :
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
    implicit none

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
    real(8) :: rdz(kdim)
    real(8) :: rdzh(kdim)

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

    real(8) :: kin_h0(1:ijdim,1:kdim),kin_h(1:ijdim,1:kdim)
    real(8) :: kin_v0(1:ijdim,1:kdim),kin_v(1:ijdim,1:kdim)
    real(8) :: vx_t(1:ijdim,1:kdim),vy_t(1:ijdim,1:kdim),w_t(1:ijdim,1:kdim)
    !
    real(8) :: wq(1:ijdim,1:kdim), ein(1:ijdim,1:kdim)

    integer :: ij, k, nq
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("precip_transport_nwater")
#endif

    rdt      = 1.D0 / dt
    rdz(:)   = 1.D0 / dz(:)
    rdzh(:)  = 1.D0 / dzh(:)

    call cnvvar_rhokin(&
      ijdim, kdim, kmin, kmax,       &
      rhogkin, rhogkin_h, rhogkin_v, &
      rhog, rhogvx, rhogvy, rhogw    )

    do nq = 1, nqmax

       if( .NOT. precipitating_flag(nq) ) then
          fprec         (:,:,nq) = 0.D0
          fprec_rhoe    (:,:,nq) = 0.D0
          fprec_rhophi  (:,:,nq) = 0.D0
          fprec_rhokin_h(:,:,nq) = 0.D0
          fprec_rhokin_v(:,:,nq) = 0.D0
          cycle
       endif

       !----- mass
       call vadv1d_getflux( &
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
          do k  = kmin, kmax
          do ij = 1,    ijdim
             rhoeq_prec(ij,k) = rhog(ij,k) * q(ij,k,nq) * CVW(nq) * tem(ij,k) * rgs(ij,k)
          enddo
          enddo

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
          do k  = kmin, kmax
          do ij = 1,    ijdim
             rhophiq_prec(ij,k) = rhog(ij,k) * q(ij,k,nq) * CNST_EGRAV * z(k) * rgs(ij,k)
          enddo
          enddo

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
          do k  = kmin, kmax
          do ij = 1,    ijdim
             rhokin_h_prec(ij,k) = q(ij,k,nq) * rhogkin_h(ij,k) * rgs(ij,k)
          enddo
          enddo

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
          do k = kmin, kmax-1
             rhokin_v_prec(:,k) = 0.5D0 * ( q(:,k+1,nq) + q(:,k,nq) ) * rhogkin_v(:,k) * rgsh(:,k)
             wq           (:,k) = 0.5D0 * ( WPREC(:,k+1,nq) + WPREC(:,k,nq) )
          enddo
          rhokin_v_prec(:,kmin-1) = 0.D0
          rhokin_v_prec(:,kmax  ) = 0.D0
          wq           (:,kmin-1) = 0.D0
          wq           (:,kmax  ) = 0.D0

          call vadv1d_getflux(         &
               ijdim,                  & !--- in
               kdim,                   &
               kmin,                   &
               kmax,                   &
               dzh(:),                  &
               z(:),                   &
               rhokin_v_prec(:,:),     & !--- in
               wq(:,:),                & !--- in
               fprec_rhokin_v(:,:,nq), & !--- out
               dt )

          fprec_rhokin_v(:,kmin-1,nq) = 0.D0
          fprec_rhokin_v(:,kmax  ,nq) = 0.D0
       endif
    enddo


    !--- tendency of density and energy
    do nq = 1,    nqmax
       do k  = kmin, kmax
       do ij = 1,    ijdim
          drho_prec(ij,k,nq) = - ( fprec(ij,k+1,nq)-fprec(ij,k,nq) ) * rdz(k)
       enddo
       enddo

       drho_prec(:,kmin-1,nq) = 0.d0
       drho_prec(:,kmax+1,nq) = 0.d0
    enddo

    do k  = kmin, kmax
    do ij = 1,    ijdim
       drhoe_prec    (ij,k) = 0.D0
       drhophi_prec  (ij,k) = 0.D0
       drhokin_h_prec(ij,k) = 0.D0

       ! sum flux
       do nq = NQW_STR, NQW_END
          drhoe_prec    (ij,k) = drhoe_prec    (ij,k) - ( fprec_rhoe    (ij,k+1,nq)-fprec_rhoe    (ij,k,nq) ) * rdz(k)
          drhophi_prec  (ij,k) = drhophi_prec  (ij,k) - ( fprec_rhophi  (ij,k+1,nq)-fprec_rhophi  (ij,k,nq) ) * rdz(k)
          drhokin_h_prec(ij,k) = drhokin_h_prec(ij,k) - ( fprec_rhokin_h(ij,k+1,nq)-fprec_rhokin_h(ij,k,nq) ) * rdz(k)
       enddo

    enddo
    enddo
    drhoe_prec    (:,kmin-1) = 0.D0
    drhoe_prec    (:,kmax+1) = 0.D0
    drhophi_prec  (:,kmin-1) = 0.D0
    drhophi_prec  (:,kmax+1) = 0.D0
    drhokin_h_prec(:,kmax+1) = 0.D0
    drhokin_h_prec(:,kmin-1) = 0.D0

    do k  = kmin, kmax-1
    do ij = 1,    ijdim
       drhokin_v_prec(ij,k) = 0.D0

       do nq = NQW_STR, NQW_END
          drhokin_v_prec(ij,k) = drhokin_v_prec(ij,k) - ( fprec_rhokin_v(ij,k+1,nq)-fprec_rhokin_v(ij,k,nq) ) * rdzh(k)
       enddo
    enddo
    enddo
    drhokin_v_prec(:,kmin-1) = 0.D0
    drhokin_v_prec(:,kmax  ) = 0.D0

    !--- momentum transport: new values due to rain
    !<------ use old values of rho and qr for consistency
    do nq = NQW_STR, NQW_END

       if( .NOT. precipitating_flag(nq) ) cycle

       do k = 1, kdim 
          rhouq(:,k) = ( rhogvx(:,k) + rhogvy(:,k) ) * q(:,k,nq) * rgs(:,k)
          rhovq(:,k) = ( rhogvx(:,k) + rhogvy(:,k) ) * q(:,k,nq) * rgs(:,k)
       enddo

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
            ijdim,            & !--- in
            kdim,             &
            kmin,             &
            kmax,             &
            dz(:),            &
            zh(:),            &
            rhovq(:,:),       & !--- in
            WPREC(:,:,nq),    & !--- in
            fprec_rhov(:,:),  & !--- out
            dt                & !--- in
            )

       !--- update rhogvx, rhogvy
       do k = kmin, kmax
          rhogvx(:,k) = rhogvx(:,k) - ( fprec_rhou(:,k+1)-fprec_rhou(:,k) ) * rdz(k)
          rhogvy(:,k) = rhogvy(:,k) - ( fprec_rhov(:,k+1)-fprec_rhov(:,k) ) * rdz(k)
       enddo


       rhowq(:,:) = 0.0D0
       do k = kmin+1,kmax
          rhowq(:,k) &
               = rhogw(:,k)* 0.5D0* (q(:,k,nq)+q(:,k-1,nq)) *rgsh(:,k)
       enddo

       wq(:,:) = 0.0D0
       do k = kmin+1,kmax-1
          wq(:,k) = 0.5D0*(WPREC(:,k,nq)+WPREC(:,k-1,nq))
       end do

       call vadv1d_getflux(   &
            ijdim,            & !--- in
            kdim,             &
            kmin,             &
            kmax,             &
            dzh(:),           &
            z(:),             &
            rhowq(:,:),       &
            wq(:,:),          &
            fprec_rhow(:,:),  &
            dt )

       !--- update rhogw
       do k = kmin, kmax-1
          rhogw(:,k) = rhogw(:,k) - ( fprec_rhow(:,k+1)-fprec_rhow(:,k) ) * rdzh(k)
       enddo
    enddo

    ! new values due to rain:
    ! Change in internal energy comes from precipitation of rain and dissipation
    ! of kinetic energy due to drag force. 
    ! See Ooyama(2001) (3.13)

    !--- update rhog, rho
    do nq = NQW_STR,NQW_END
       rhog(:,:) = rhog(:,:) + drho_prec(:,:,nq)
    enddo
    rho(:,:) = rhog(:,:)/gsgam2(:,:)

    !--- update rhogqr
    do nq = 1, nqmax
       rhogq(:,:,nq) = rhogq(:,:,nq) + drho_prec(:,:,nq)

       q(:,:,nq) = rhogq(:,:,nq) / rhog(:,:)
    enddo


    !--- update rhoge, rhogkin
    rhoge(:,:) = rhoge(:,:) + drhoe_prec(:,:) + drhophi_prec(:,:)
    do nq = NQW_STR,NQW_END
    do k = 1, kdim
       rhoge(:,k) = rhoge(:,k) - CNST_EGRAV * z(k) * ( drho_prec(:,k,nq) )
    enddo
    enddo

    rhogkin_h(:,:) = rhogkin_h(:,:) + drhokin_h_prec(:,:)
    rhogkin_v(:,:) = rhogkin_v(:,:) + drhokin_v_prec(:,:)

    !--- update qd
    call thrmdyn_qd( qd, q )

    !--- fix kinetic energy
!    if ( trim(PRCIP_TRN_ECORRECT) == 'KIN2KIN') then
!       do k  = kmin, kmax
!       do ij = 1,    ijdim
!
!          vx_t  (ij,k) = rhogvx(ij,k) / rhog(ij,k)
!          vy_t  (ij,k) = rhogvy(ij,k) / rhog(ij,k)
!
!          kin_h (ij,k) = 0.5D0 * ( vx_t(ij,k)**2 + vy_t(ij,k)**2 )
!          kin_h0(ij,k) = rhogkin_h(ij,k) / rhog(ij,k)
!
!          if ( kin_h(ij,k) > 1.D-20 ) then
!             vx_t(ij,k) = vx_t(ij,k) * sqrt( abs(kin_h0(ij,k)/kin_h(ij,k)) )
!             vy_t(ij,k) = vy_t(ij,k) * sqrt( abs(kin_h0(ij,k)/kin_h(ij,k)) )
!          endif
!
!          rhogvx(ij,k) = rhog(ij,k) * vx_t(ij,k)
!          rhogvy(ij,k) = rhog(ij,k) * vy_t(ij,k)
!       enddo
!       enddo
!
!       do k = kmin , kmax+1
!       do ij = 1,    ijdim
!          w_t   (ij,k) = 2.D0 * rhogw(ij,k) / ( rhog(ij,k+1)+rhog(ij,k) )
!
!          kin_v (ij,k) = 0.5D0 * w_t(ij,k)**2
!          kin_v0(ij,k) = rhogkin_v(ij,k) / ( rhog(ij,k+1)+rhog(ij,k) )
!
!          if ( kin_v(ij,k) > 1.D-20 ) then
!             w_t(ij,k) = w_t(ij,k) * sqrt( abs(kin_v0(ij,k)/kin_v(ij,k)) )
!          endif
!
!          rhogw(ij,k+1) = 0.5D0 * ( rhog(ij,k+1)+rhog(ij,k) ) * w_t(ij,k+1)
!
!       enddo
!       enddo
!    endif

    !--- update temerature & pressure
    ein(:,:) =  rhoge(:,:) / rhog(:,:)
    call thrmdyn_tempre( tem, pre, ein, rho, qd, q )

    precip        (:,:) = 0.D0
    precip_rhoe   (:)   = 0.D0
    precip_rhophi (:)   = 0.D0
    precip_rhokin (:)   = 0.D0
    precip_lh_heat(:)   = 0.D0

    do nq=NQL_STR,NQL_END
       precip(:,1) = precip(:,1) -fprec(:,kmin,nq)*rdt
    enddo

    do nq = NQS_STR,NQS_END
       precip(:,2) = precip(:,2) -fprec(:,kmin,nq)*rdt
    enddo

    do nq=NQS_STR,NQS_END
       precip_lh_heat(:) = precip_lh_heat(:) + fprec(:,kmin,nq) * LHF * rdt
    end do

    do nq = NQW_STR, NQW_END
       precip_rhoe  (:) = precip_rhoe  (:) - fprec_rhoe    (:,kmin,nq) * rdt
       precip_rhophi(:) = precip_rhophi(:) - fprec_rhophi  (:,kmin,nq) * rdt
       precip_rhokin(:) = precip_rhokin(:) - fprec_rhokin_h(:,kmin,nq) * rdt &
                                           - fprec_rhokin_v(:,kmin,nq) * rdt
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("precip_transport_nwater")
#endif

    return
  end subroutine precip_transport_nwater

  !-----------------------------------------------------------------------------
  !> calculate kinetic energy
  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhokin( &
      ijdim, kdim, kmin, kmax,       &
      rhogkin, rhogkin_h, rhogkin_v, &
      rhog, rhogvx, rhogvy, rhogw    )
    implicit none

    integer, intent(in)  :: ijdim
    integer, intent(in)  :: kdim
    integer, intent(in)  :: kmin
    integer, intent(in)  :: kmax

    real(8), intent(out) :: rhogkin  (ijdim,kdim)
    real(8), intent(out) :: rhogkin_h(ijdim,kdim)
    real(8), intent(out) :: rhogkin_v(ijdim,kdim)

    real(8), intent(in)  :: rhog     (ijdim,kdim)
    real(8), intent(in)  :: rhogvx   (ijdim,kdim)
    real(8), intent(in)  :: rhogvy   (ijdim,kdim)
    real(8), intent(in)  :: rhogw    (ijdim,kdim)

    real(8) :: rhog_h

    integer :: ij, k
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("cnvvar_rhokin")
#endif

    !--- horizontal kinetic energy
    do k  = kmin, kmax
    do ij = 1,    ijdim
       rhogkin_h(ij,k) = 0.5D0 * rhog(ij,k) * ( ( rhogvx(ij,k)/rhog(ij,k))**2 &
                                              + ( rhogvy(ij,k)/rhog(ij,k))**2 )
    enddo
    enddo

    do ij = 1,    ijdim
       rhogkin_h(ij,kmin-1) = ADM_VMISS
       rhogkin_h(ij,kmax+1) = ADM_VMISS
    enddo

    !--- vertical kinetic energy
    do k  = kmin, kmax-1
    do ij = 1,    ijdim
       rhog_h = 0.5D0 * ( rhog(ij,k+1) + rhog(ij,k) )

       rhogkin_v(ij,k) = 0.5D0 * rhog_h * ( ( rhogw(ij,k)/rhog_h )**2 )
    enddo
    enddo

    do ij = 1,    ijdim
       rhogkin_v(ij,kmin-1) = 0.D0
       rhogkin_v(ij,kmax)   = 0.D0
    enddo

    do k  = kmin, kmax
    do ij = 1,    ijdim
       rhogkin(ij,k) = rhogkin_h(ij,k)                              &
                     + 0.5D0 * (rhogkin_v(ij,k+1) + rhogkin_v(ij,k) )
    enddo
    enddo

    do ij = 1,    ijdim
       rhogkin(ij,kmin-1) = ADM_VMISS
       rhogkin(ij,kmax+1) = ADM_VMISS
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("cnvvar_rhokin")
#endif

    return
  end subroutine cnvvar_rhokin


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
    rhofh_cell(:,:) = 0.d0
    !<--- including the condition of [fh(:,kmax+1) = 0.d0] 
    !<--- at the top boundary
    do k = kmin, kmax
       do ij = 1, ijdim
          rhofh_cell(ij,k) = rhof(ij,kcell(ij,k))
       end do
    end do
    !
    do k = kmin, kmax
       where ( kcell(:,k) == 0 )
          rhofh_cell(:,k) =  rhof(:,kmin)
       end where
       where ( kcell(:,k) == kmax+1 )
          rhofh_cell(:,k) =  rhof(:,kmax)
       end where
    end do
    !
    frhof(:,:) = frhof(:,:) + rhofh_cell(:,:) * zdis(:,:)
    !

#ifdef _FPCOLL_
call STOP_COLLECTION("vadv1d_getflux")
#endif

    return
  end subroutine vadv1d_getflux

end module mod_precip_transport
