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
       wprcp,                 & !--- IN :
       precipitating_flag,    & !--- IN :
       precip,                & !--- OUT :
       gsgam2,                & !--- IN :
       gsgam2h,               & !--- IN :
       rgs,                   & !--- IN :
       rgsh,                  & !--- IN :
       dt                     ) !--- IN :
    use mod_atmos_cnst, only : &
         CVW,               &
         LHF,               &
         CNST_CP,           &
         CNST_CV,           &
         CNST_EGRAV,        &
         CNST_RAIR,         &
         CNST_RVAP
    use mod_atmos_vars, only: &
       nqmax   => A_QA,  &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG, &
       I_NG
    use mod_thrmdyn, only : &
       thrmdyn_qd, &
       thrmdyn_tempre
    implicit none

    integer, intent(in)    :: ijdim
    integer, intent(in)    :: kdim
    integer, intent(in)    :: kmin
    integer, intent(in)    :: kmax
    real(8), intent(inout) :: rhog  (ijdim,kdim)
    real(8), intent(inout) :: rhogvx(ijdim,kdim)
    real(8), intent(inout) :: rhogvy(ijdim,kdim)
    real(8), intent(inout) :: rhogw (ijdim,kdim)
    real(8), intent(inout) :: rhoge (ijdim,kdim)
    real(8), intent(inout) :: rhogq (ijdim,kdim,nqmax)

    real(8), intent(inout) :: rho(ijdim,kdim)
    real(8), intent(inout) :: tem(ijdim,kdim)
    real(8), intent(inout) :: pre(ijdim,kdim)
    real(8), intent(inout) :: q(ijdim,kdim,1:nqmax)
    real(8), intent(out)   :: qd(ijdim,kdim)

    real(8), intent(in) :: wprcp(ijdim,kdim,1:nqmax)
    logical, intent(in) :: precipitating_flag(1:nqmax)
    real(8), intent(in) :: z (1:kdim) ! height
    real(8), intent(in) :: zh(1:kdim) ! height
    real(8), intent(in) :: dz(1:kdim) ! height
    real(8), intent(in) :: dzh(1:kdim)! height

    real(8), intent(out) :: precip(1:ijdim,2) ! rain or snow

    real(8), intent(in) :: gsgam2(ijdim,kdim)
    real(8), intent(in) :: gsgam2h(ijdim,kdim)
    real(8), intent(in) :: rgs(ijdim,kdim)
    real(8), intent(in) :: rgsh(ijdim,kdim)

    real(8), intent(in) :: dt

    real(8) :: rdt
    real(8) :: rdz(kdim)
    real(8) :: rdzh(kdim)

    real(8) :: rhogkin  (ijdim,kdim)
    real(8) :: rhogkin_h(ijdim,kdim)
    real(8) :: rhogkin_v(ijdim,kdim)

    real(8) :: rhoq     (ijdim,kdim)
    real(8) :: rhoeq    (ijdim,kdim)
    real(8) :: rhophiq  (ijdim,kdim)
    real(8) :: rhokinq_h(ijdim,kdim)
    real(8) :: rhouq    (ijdim,kdim)
    real(8) :: rhovq    (ijdim,kdim)

    real(8) :: whprcp   (ijdim,kdim)
    real(8) :: rhokinq_v(ijdim,kdim)
    real(8) :: rhowq    (ijdim,kdim)

    real(8) :: flx_rhoq    (ijdim,kdim)
    real(8) :: flx_rhoe    (ijdim,kdim)
    real(8) :: flx_rhophi  (ijdim,kdim)
    real(8) :: flx_rhokin_h(ijdim,kdim)
    real(8) :: flx_rhokin_v(ijdim,kdim)
    real(8) :: flx_rhou    (ijdim,kdim)
    real(8) :: flx_rhov    (ijdim,kdim)
    real(8) :: flx_rhow    (ijdim,kdim)
    real(8) :: drho_prec

    real(8) :: kin_h0(ijdim,kdim), kin_h(ijdim,kdim)
    real(8) :: kin_v0(ijdim,kdim), kin_v(ijdim,kdim)
    real(8) :: vx_t(ijdim,kdim)
    real(8) :: vy_t(ijdim,kdim)
    real(8) :: w_t (ijdim,kdim)
    real(8) :: ein (ijdim,kdim)

    integer :: ij, k, nq, n
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

    do n  = 1, 2
    do ij = 1, ijdim
       precip(ij,n) = 0.D0
    enddo
    enddo

    ! new values due to rain:
    ! Change in internal energy comes from precipitation of rain and dissipation
    ! of kinetic energy due to drag force. 
    ! See Ooyama(2001) (3.13)
    do nq = I_QC, I_NG

       do k  = kmin, kmax
       do ij = 1,    ijdim
          rhoq     (ij,k) = rhog     (ij,k) * q(ij,k,nq)
          rhoeq    (ij,k) = rhog     (ij,k) * q(ij,k,nq) * CVW(nq) * tem(ij,k)
          rhophiq  (ij,k) = rhog     (ij,k) * q(ij,k,nq) * CNST_EGRAV * z(k)
          rhokinq_h(ij,k) = rhogkin_h(ij,k) * q(ij,k,nq)
          rhouq    (ij,k) = rhogvx   (ij,k) * q(ij,k,nq)
          rhovq    (ij,k) = rhogvy   (ij,k) * q(ij,k,nq)
       enddo
       enddo

       call vadv1d_getflux( flx_rhoq    (:,:),   &
                            flx_rhoe    (:,:),   &
                            flx_rhophi  (:,:),   &
                            flx_rhokin_h(:,:),   &
                            flx_rhou    (:,:),   &
                            flx_rhov    (:,:),   &
                            rhoq        (:,:),   & !--- mass
                            rhoeq       (:,:),   & !--- internal energy
                            rhophiq     (:,:),   & !--- potential energy
                            rhokinq_h   (:,:),   & !--- horizontal kinetic energy
                            rhouq       (:,:),   & !--- rhogvx
                            rhovq       (:,:),   & !--- rhogvy
                            wprcp       (:,:,nq) )

       do k  = kmin, kmax-1
       do ij = 1,    ijdim
          drho_prec = - ( flx_rhoq(ij,k)-flx_rhoq(ij,k-1) ) * rdz(k)

          rhog (ij,k)    = rhog (ij,k)    + drho_prec
          rhogq(ij,k,nq) = rhogq(ij,k,nq) + drho_prec

          rhoge    (ij,k) = rhoge    (ij,k) - CNST_EGRAV * z(k) * drho_prec &
                                            - ( flx_rhoe    (ij,k)-flx_rhoe    (ij,k-1) ) * rdz(k) &
                                            - ( flx_rhophi  (ij,k)-flx_rhophi  (ij,k-1) ) * rdz(k)
          rhogkin_h(ij,k) = rhogkin_h(ij,k) - ( flx_rhokin_h(ij,k)-flx_rhokin_h(ij,k-1) ) * rdz(k)
          rhogvx   (ij,k) = rhogvx   (ij,k) - ( flx_rhou    (ij,k)-flx_rhou    (ij,k-1) ) * rdz(k)
          rhogvy   (ij,k) = rhogvy   (ij,k) - ( flx_rhov    (ij,k)-flx_rhov    (ij,k-1) ) * rdz(k)
       enddo
       enddo

       if (      nq == I_QC &
            .OR. nq == I_QR ) then

          do ij = 1,    ijdim
             precip(ij,1) = precip(ij,1) - flx_rhoq(ij,kmin-1) * rdt
          enddo

       elseif(      nq == I_QI & 
               .OR. nq == I_QS & 
               .OR. nq == I_QG ) then

          do ij = 1,    ijdim
             precip(ij,2) = precip(ij,2) - flx_rhoq(ij,kmin-1) * rdt
          enddo

       endif

       do k  = kmin, kmax-1
       do ij = 1,    ijdim
          whprcp   (ij,k) = 0.5D0 * ( wprcp(ij,k+1,nq)+wprcp(ij,k,nq) )
          rhokinq_v(ij,k) = rhogkin_v(ij,k) * 0.5D0 * ( q(ij,k+1,nq)+q(ij,k,nq) )
          rhowq    (ij,k) = rhogw    (ij,k) * 0.5D0 * ( q(ij,k+1,nq)+q(ij,k,nq) )
       enddo
       enddo
       do ij = 1, ijdim
          rhowq    (ij,kmin-1) = 0.D0
          rhowq    (ij,kmax  ) = 0.D0
          rhokinq_v(ij,kmin-1) = 0.D0
          rhokinq_v(ij,kmax  ) = 0.D0
          whprcp   (ij,kmin-1) = 0.D0
          whprcp   (ij,kmax  ) = 0.D0
       enddo

       call vadv1d_getflux_h( flx_rhokin_v(:,:), &
                              flx_rhow    (:,:), &
                              rhokinq_v   (:,:), & !--- vertical kinetic energy
                              rhowq       (:,:), & !--- rhogw
                              whprcp      (:,:)  )

       do k  = kmin, kmax-1
       do ij = 1,    ijdim
          rhogkin_v(ij,k) = rhogkin_v(ij,k) - ( flx_rhokin_v(ij,k+1)-flx_rhokin_v(ij,k) ) * rdzh(k)
          rhogw    (ij,k) = rhogw    (ij,k) - ( flx_rhow    (ij,k+1)-flx_rhow    (ij,k) ) * rdzh(k)
       enddo
       enddo

    enddo

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
    do k  = kmin, kmax
    do ij = 1,    ijdim
       rho(ij,k) = rhog(ij,k)
       ein(ij,k) = rhoge(ij,k) / rhog(ij,k)
    enddo
    enddo

    do nq = 1,    nqmax
    do k  = kmin, kmax
    do ij = 1,    ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo

    call thrmdyn_qd    ( qd, q )
    call thrmdyn_tempre( tem, pre, ein, rho, qd, q )

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

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux( &
      flx_rhoq,     &
      flx_rhoe,     &
      flx_rhophi,   &
      flx_rhokin_h, &
      flx_rhou,     &
      flx_rhov,     &
      rhoq,         &
      rhoe,         &
      rhophi,       &
      rhokin_h,     &
      rhou,         &
      rhov,         &
      wprcp         )
    use mod_time, only: &
       dt  => TIME_DTSEC_ATMOS_PHY_MP
    use mod_grid, only: &
       KA  => GRID_KA,  &
       KS  => GRID_KS,  &
       KE  => GRID_KE,  &
       IJA => GRID_IJA, &
       IJS => GRID_IJS, &
       IJE => GRID_IJE, &
       zh  => GRID_FZ,  &
       dz  => GRID_CDZ
    implicit none

    real(8), intent(out) :: flx_rhoq    (IJA,KA) ! flux rho for each q
    real(8), intent(out) :: flx_rhoe    (IJA,KA) ! flux rho * ein for each q
    real(8), intent(out) :: flx_rhophi  (IJA,KA) ! flux rho * phi for each q
    real(8), intent(out) :: flx_rhokin_h(IJA,KA) ! flux rho * kin for each q
    real(8), intent(out) :: flx_rhou    (IJA,KA) ! flux rho * u for each q
    real(8), intent(out) :: flx_rhov    (IJA,KA) ! flux rho * v for each q

    real(8), intent(in)  :: rhoq    (IJA,KA) ! rho for each q
    real(8), intent(in)  :: rhoe    (IJA,KA) ! rho * ein for each q
    real(8), intent(in)  :: rhophi  (IJA,KA) ! rho * phi for each q
    real(8), intent(in)  :: rhokin_h(IJA,KA) ! rho * kin for each q
    real(8), intent(in)  :: rhou    (IJA,KA) ! rho * u for each q
    real(8), intent(in)  :: rhov    (IJA,KA) ! rho * v for each q

    real(8), intent(in)  :: wprcp(IJA,KA) ! vertical velocity at the full level

    real(8) :: flx_rhoq_cell    (IJA,KA)
    real(8) :: flx_rhoe_cell    (IJA,KA)
    real(8) :: flx_rhophi_cell  (IJA,KA)
    real(8) :: flx_rhokin_h_cell(IJA,KA)
    real(8) :: flx_rhou_cell    (IJA,KA)
    real(8) :: flx_rhov_cell    (IJA,KA)

    real(8) :: velz (IJA,KA)
    real(8) :: zdis (IJA,KA)
    integer :: kcell(IJA,KA)
    integer :: kcell_max(KA)
    integer :: kcell_min(KA)
    real(8) :: zzmax, zzmin
    real(8) :: dwdzp1, dwdzm1, dwdzpm, fact

    integer :: ij, k, k2, kk
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("vadv1d_getflux")
#endif

    !--- vetical velocity at the half level
    do k  = KS,  KE-1
    do ij = IJS, IJE
       velz(ij,k) = 0.5D0 * ( wprcp(ij,k+1) + wprcp(ij,k) )
    enddo
    enddo
    do ij = IJS, IJE
       velz(ij,KS-1) = wprcp(ij,KS  )
       velz(ij,KE  ) = wprcp(ij,KE-1)
    enddo

    !--- calculation of distance of cell wall during dt
    do k  = KS,  KE-1
    do ij = IJS, IJE
       dwdzp1 = ( velz(ij,k+1)-velz(ij,k  ) ) / ( dz(k)         )
       dwdzm1 = ( velz(ij,k  )-velz(ij,k-1) ) / (       dz(k-1) )
       dwdzpm = ( velz(ij,k+1)-velz(ij,k-1) ) / ( dz(k)+dz(k-1) )

       zdis(ij,k) = velz(ij,k) * dt                       &
                  - velz(ij,k) * dt**2 * ( dwdzpm ) * 0.5D0 &
                  + velz(ij,k) * dt**3 * ( dwdzpm**2 + 2.D0 * (dwdzp1-dwdzm1) * velz(ij,k) / (dz(k)+dz(k-1)) ) / 6.D0
    enddo
    enddo

    !--- bottom and top boundary for zdis
    k = KS-1
    do ij = IJS, IJE
       dwdzp1 = ( velz(ij,k+1)-velz(ij,k  ) ) / ( dz(k  ) )

       zdis(ij,k) = velz(ij,k) * dt                     &
                  - velz(ij,k) * dt**2 * dwdzp1 * 0.5D0
    enddo
    k = KE
    do ij = IJS, IJE
       dwdzm1 = ( velz(ij,k  )-velz(ij,k-1) ) / ( dz(k-1) )

       zdis(ij,k) = velz(ij,k) * dt                     &
                  - velz(ij,k) * dt**2 * dwdzm1 * 0.5D0
    enddo

    do k  = KS-1, KE
    do ij = IJS, IJE
       kcell(ij,k) = k
    enddo
    enddo

    !------ setup limiter of max and min of kcell
    do k = KS-1, KE-1
       zzmax = maxval( zdis(:,k) )
       zzmin = minval( zdis(:,k) )
 
       do ij = IJS, IJE
          kcell_min(k) = k
          kcell_max(k) = k

          if ( zzmax > 0.D0 ) then
             do k2 = k, KS-1, -1
                if (       zh(k2)   <= zh(k) - zzmax &
                     .AND. zh(k2+1) >  zh(k) - zzmax ) then
                   kcell_min(k) = k2
                   exit
                endif
             enddo
          endif

          if ( zzmin < 0.d0 ) then
             do k2 = k, KE-1
                if (       zh(k2)   <= zh(k) - zzmin &
                     .AND. zh(k2+1) >  zh(k) - zzmin ) then
                   kcell_max(k) = k2
                   exit
                endif
             enddo
          endif

       enddo
    enddo

    !------ determine the kcell at each point.
    do k = KS-1, KE-1
       if (       kcell_min(k) == k &
            .AND. kcell_max(k) == k ) then

          do ij = IJS, IJE
             kcell(ij,k) = k
          enddo

       else

          do ij = IJS, IJE
             kcell(ij,k) = 0
             do k2 = kcell_min(k), kcell_max(k)
                if (       zh(k2  ) < zh(k)-zdis(ij,k) &
                     .AND. zh(k2+1) > zh(k)-zdis(ij,k) ) then

                   kcell(ij,k) = max( kcell(ij,k), k2 )

                endif
             enddo
          enddo

       endif
    enddo

    !--- set zero value for norain, see if block in do loop
    do k  = KS-1, KE
    do ij = IJS, IJE
       flx_rhoq    (ij,k) = 0.D0
       flx_rhoe    (ij,k) = 0.D0
       flx_rhophi  (ij,k) = 0.D0
       flx_rhokin_h(ij,k) = 0.D0
       flx_rhou    (ij,k) = 0.D0
       flx_rhov    (ij,k) = 0.D0
    enddo
    enddo

    !------ integration in the integer cells
    do k = KS-1, KE-1
       if (       kcell_min(k) == k &
            .AND. kcell_max(k) == k ) then
          ! do nothing
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             fact = ( sign(1,k2-kcell(ij,k)-1) + 1.D0 ) / 2.D0 &
                  * ( sign(1,k-k2-1)           + 1.D0 ) / 2.D0 &
                  - ( sign(1,k2-k)             + 1.D0 ) / 2.D0 &
                  * ( sign(1,kcell(ij,k)-k2-1) + 1.D0 ) / 2.D0

             do ij = IJS, IJE
                flx_rhoq    (ij,k) = flx_rhoq    (ij,k) + rhoq    (ij,k) * dz(k2) * fact
                flx_rhoe    (ij,k) = flx_rhoe    (ij,k) + rhoe    (ij,k) * dz(k2) * fact
                flx_rhophi  (ij,k) = flx_rhophi  (ij,k) + rhophi  (ij,k) * dz(k2) * fact
                flx_rhokin_h(ij,k) = flx_rhokin_h(ij,k) + rhokin_h(ij,k) * dz(k2) * fact
                flx_rhou    (ij,k) = flx_rhou    (ij,k) + rhou    (ij,k) * dz(k2) * fact
                flx_rhov    (ij,k) = flx_rhov    (ij,k) + rhov    (ij,k) * dz(k2) * fact

                zdis (ij,k) = zdis (ij,k) - dz(k2) * fact
             enddo
          enddo
       endif
    enddo

    do k  = KS-1, KE-1
    do ij = IJS, IJE

       kk = kcell(ij,k)
       if    ( kcell(ij,k) == 0 ) then
          kk = KS
       elseif( kcell(ij,k) == KE ) then
          kk = KE
       endif

       flx_rhoq_cell    (ij,k) = flx_rhoq    (ij,kk)
       flx_rhoe_cell    (ij,k) = flx_rhoe    (ij,kk)
       flx_rhophi_cell  (ij,k) = flx_rhophi  (ij,kk)
       flx_rhokin_h_cell(ij,k) = flx_rhokin_h(ij,kk)
       flx_rhou_cell    (ij,k) = flx_rhou    (ij,kk)
       flx_rhov_cell    (ij,k) = flx_rhov    (ij,kk)
    enddo
    enddo

    do k  = KS-1, KE
    do ij = IJS, IJE
       flx_rhoq    (ij,k) = flx_rhoq    (ij,k) + flx_rhoq_cell    (ij,k) * zdis(ij,k)
       flx_rhoe    (ij,k) = flx_rhoe    (ij,k) + flx_rhoe_cell    (ij,k) * zdis(ij,k)
       flx_rhophi  (ij,k) = flx_rhophi  (ij,k) + flx_rhophi_cell  (ij,k) * zdis(ij,k)
       flx_rhokin_h(ij,k) = flx_rhokin_h(ij,k) + flx_rhokin_h_cell(ij,k) * zdis(ij,k)
       flx_rhou    (ij,k) = flx_rhou    (ij,k) + flx_rhou_cell    (ij,k) * zdis(ij,k)
       flx_rhov    (ij,k) = flx_rhov    (ij,k) + flx_rhov_cell    (ij,k) * zdis(ij,k)
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("vadv1d_getflux")
#endif

    return
  end subroutine vadv1d_getflux

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_h( &
      flx_rhokin_v, &
      flx_rhow,     &
      rhokin_v,     &
      rhow,         &
      wprcp         )
    use mod_time, only: &
       dt  => TIME_DTSEC_ATMOS_PHY_MP
    use mod_grid, only: &
       KA  => GRID_KA,  &
       KS  => GRID_KS,  &
       KE  => GRID_KE,  &
       IJA => GRID_IJA, &
       IJS => GRID_IJS, &
       IJE => GRID_IJE, &
       z   => GRID_CZ,  &
       dzh => GRID_FDZ
    implicit none

    real(8), intent(out) :: flx_rhokin_v(IJA,KA) ! flux rho * kin for each q
    real(8), intent(out) :: flx_rhow    (IJA,KA) ! flux rho * w for each q

    real(8), intent(in)  :: rhokin_v(IJA,KA) ! rho * kin for each q
    real(8), intent(in)  :: rhow    (IJA,KA) ! rho * w for each q

    real(8), intent(in)  :: wprcp(IJA,KA) ! vertical velocity at the half level

    real(8) :: flx_rhokin_v_cell(IJA,KA)
    real(8) :: flx_rhow_cell    (IJA,KA)

    real(8) :: velz (IJA,KA)
    real(8) :: zdis (IJA,KA)
    integer :: kcell(IJA,KA)
    integer :: kcell_max(KA)
    integer :: kcell_min(KA)
    real(8) :: zzmax, zzmin
    real(8) :: dwdzp1, dwdzm1, dwdzpm, fact

    integer :: ij, k, k2, kk
    !---------------------------------------------------------------------------

#ifdef _FPCOLL_
call START_COLLECTION("vadv1d_getflux")
#endif

    !--- vetical velocity at the full level
    do k  = KS,  KE
    do ij = IJS, IJE
       velz(ij,k) = 0.5D0 * ( wprcp(ij,k) + wprcp(ij,k-1) )
    enddo
    enddo

    !--- calculation of distance of cell wall during dt
    do k  = KS+1, KE-1
    do ij = IJS,  IJE
       dwdzp1 = ( velz(ij,k+1)-velz(ij,k  ) ) / ( dzh(k)         )
       dwdzm1 = ( velz(ij,k  )-velz(ij,k-1) ) / (        dzh(k-1) )
       dwdzpm = ( velz(ij,k+1)-velz(ij,k-1) ) / ( dzh(k)+dzh(k-1) )

       zdis(ij,k) = velz(ij,k) * dt                       &
                  - velz(ij,k) * dt**2 * ( dwdzpm ) * 0.5D0 &
                  + velz(ij,k) * dt**3 * ( dwdzpm**2 + 2.D0 * (dwdzp1-dwdzm1) * velz(ij,k) / (dzh(k)+dzh(k-1)) ) / 6.D0
    enddo
    enddo

    !--- bottom and top boundary for zdis
    k = KS
    do ij = IJS, IJE
       dwdzp1 = ( velz(ij,k+1)-velz(ij,k  ) ) / ( dzh(k  ) )

       zdis(ij,k) = velz(ij,k) * dt                     &
                  - velz(ij,k) * dt**2 * dwdzp1 * 0.5D0
    enddo
    k = KE
    do ij = IJS, IJE
       dwdzm1 = ( velz(ij,k  )-velz(ij,k-1) ) / ( dzh(k-1) )

       zdis(ij,k) = velz(ij,k) * dt                     &
                  - velz(ij,k) * dt**2 * dwdzm1 * 0.5D0
    enddo

    do k  = KS-1, KE
    do ij = IJS, IJE
       kcell(ij,k) = k
    enddo
    enddo

    !------ setup limiter of max and min of kcell
    do k = KS, KE-1
       zzmax = maxval( zdis(:,k) )
       zzmin = minval( zdis(:,k) )
 
       do ij = IJS, IJE
          kcell_min(k) = k
          kcell_max(k) = k

          if ( zzmax > 0.D0 ) then
             do k2 = k, KS, -1
                if (       z(k2)   <= z(k) - zzmax &
                     .AND. z(k2+1) >  z(k) - zzmax ) then
                   kcell_min(k) = k2
                   exit
                endif
             enddo
          endif

          if ( zzmin < 0.d0 ) then
             do k2 = k, KE-1
                if (       z(k2)   <= z(k) - zzmin &
                     .AND. z(k2+1) >  z(k) - zzmin ) then
                   kcell_max(k) = k2
                   exit
                endif
             enddo
          endif

       enddo
    enddo

    !------ determine the kcell at each point.
    do k = KS, KE-1
       if (       kcell_min(k) == k &
            .AND. kcell_max(k) == k ) then

          do ij = IJS, IJE
             kcell(ij,k) = k
          enddo

       else

          do ij = IJS, IJE
             kcell(ij,k) = 0
             do k2 = kcell_min(k), kcell_max(k)
                if (       z(k2  ) < z(k)-zdis(ij,k) &
                     .AND. z(k2+1) > z(k)-zdis(ij,k) ) then

                   kcell(ij,k) = max( kcell(ij,k), k2 )

                endif
             enddo
          enddo

       endif
    enddo

    !--- set zero value for norain, see if block in do loop
    do k  = KS-1, KE
    do ij = IJS, IJE
       flx_rhokin_v(ij,k) = 0.D0
       flx_rhow    (ij,k) = 0.D0
    enddo
    enddo

    !------ integration in the integer cells
    do k = KS, KE-1
       if (       kcell_min(k) == k &
            .AND. kcell_max(k) == k ) then
          ! do nothing
       else
          do k2 = kcell_min(k), kcell_max(k)
             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             fact = ( sign(1,k2-kcell(ij,k)-1) + 1.D0 ) / 2.D0 &
                  * ( sign(1,k-k2-1)           + 1.D0 ) / 2.D0 &
                  - ( sign(1,k2-k)             + 1.D0 ) / 2.D0 &
                  * ( sign(1,kcell(ij,k)-k2-1) + 1.D0 ) / 2.D0

             do ij = IJS, IJE
                flx_rhokin_v(ij,k) = flx_rhokin_v(ij,k) + rhokin_v(ij,k2) * dzh(k2) * fact
                flx_rhow    (ij,k) = flx_rhow    (ij,k) + rhow    (ij,k2) * dzh(k2) * fact

                zdis(ij,k) = zdis(ij,k) - dzh(k2) * fact
             enddo
          enddo
       endif
    enddo

    do k  = KS, KE-1
    do ij = IJS, IJE

       kk = kcell(ij,k)
       if    ( kcell(ij,k) == 0 ) then
          kk = KS
       elseif( kcell(ij,k) == KE ) then
          kk = KE-1
       endif

       flx_rhokin_v_cell(ij,k) = flx_rhokin_v(ij,kk)
       flx_rhow_cell    (ij,k) = flx_rhow    (ij,kk)
    enddo
    enddo

    do k  = KS, KE-1
    do ij = IJS, IJE
       flx_rhokin_v(ij,k) = flx_rhokin_v(ij,k) + flx_rhokin_v_cell(ij,k) * zdis(ij,k)
       flx_rhow    (ij,k) = flx_rhow    (ij,k) + flx_rhow_cell    (ij,k) * zdis(ij,k)
    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("vadv1d_getflux")
#endif

    return
  end subroutine vadv1d_getflux_h

end module mod_precip_transport
