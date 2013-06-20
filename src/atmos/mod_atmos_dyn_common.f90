!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics common
!!
!! @par Description
!!          common subroutines for Atmospheric dynamical process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-18 (S.Nishizawa) [new] splited from dynamical core
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_dyn_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
#ifdef DEBUG
  use mod_debug, only: &
       CHECK
#endif
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_fct

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  ! advection settings
  real(RP), public, parameter :: FACT_N =   7.0_RP / 12.0_RP !  7/12: fourth, 1: second
  real(RP), public, parameter :: FACT_F = - 1.0_RP / 12.0_RP ! -1/12: fourth, 0: second

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  subroutine ATMOS_DYN_fct(      &
       qflx_anti,                &
       phi_in, qflx_hi, qflx_lo, &
       rdz, rdx, rdy, dtrk       )
    use mod_const, only: &
       UNDEF => CONST_UNDEF, &
       IUNDEF => CONST_UNDEF2, &
       EPSILON => CONST_EPS
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_atmos_vars, only: &
       ZDIR, &
       XDIR, &
       YDIR
    implicit none

    real(RP), intent(out) :: qflx_anti(KA,IA,JA,3)

    real(RP), intent(in) :: phi_in(KA,IA,JA) ! physical quantity
    real(RP), intent(in) :: qflx_hi(KA,IA,JA,3)
    real(RP), intent(in) :: qflx_lo(KA,IA,JA,3)

    real(RP), intent(in) :: RDZ(:)
    real(RP), intent(in) :: RDX(:)
    real(RP), intent(in) :: RDY(:)
    real(RP), intent(in) :: dtrk

    ! factor for FCT (work space)


    ! work for FCT
    real(RP) :: phi_lo(KA,IA,JA)
    real(RP) :: pjpls(KA,IA,JA)
    real(RP) :: pjmns(KA,IA,JA)
    real(RP) :: qjpls(KA,IA,JA)
    real(RP) :: qjmns(KA,IA,JA)
    real(RP) :: rjpls(KA,IA,JA)
    real(RP) :: rjmns(KA,IA,JA)

    real(RP) :: zerosw, dirsw

    integer :: IIS, IIE, JJS, JJE
    integer :: k, i, j, ijs

#ifdef DEBUG
    qflx_anti(:,:,:,:) = UNDEF

    pjpls(:,:,:) = UNDEF
    pjmns(:,:,:) = UNDEF
    qjpls(:,:,:) = UNDEF
    qjmns(:,:,:) = UNDEF
    rjpls(:,:,:) = UNDEF
    rjmns(:,:,:) = UNDEF
#endif

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,ZDIR) )
#endif
          qflx_anti(k,i,j,ZDIR) = qflx_hi(k,i,j,ZDIR) - qflx_lo(k,i,j,ZDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS  , JJE
       do i = IIS-1, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,XDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,XDIR) )
#endif
          qflx_anti(k,i,j,XDIR) = qflx_hi(k,i,j,XDIR) - qflx_lo(k,i,j,XDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE
       do i = IIS  , IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_hi(k,i,j,YDIR) )
          call CHECK( __LINE__, qflx_lo(k,i,j,YDIR) )
#endif
          qflx_anti(k,i,j,YDIR) = qflx_hi(k,i,j,YDIR) - qflx_lo(k,i,j,YDIR)
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- update monotone scheme
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, phi_in(k,i,j) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_lo(k  ,i  ,j-1,YDIR) )
#endif
          phi_lo(k,i,j) = phi_in(k,i,j) &
               + dtrk * ( - ( ( qflx_lo(k,i,j,ZDIR)-qflx_lo(k-1,i  ,j  ,ZDIR) ) * RDZ(k) &
                            + ( qflx_lo(k,i,j,XDIR)-qflx_lo(k  ,i-1,j  ,XDIR) ) * RDX(i) &
                            + ( qflx_lo(k,i,j,YDIR)-qflx_lo(k  ,i  ,j-1,YDIR) ) * RDY(j) ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net incoming quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
          pjpls(k,i,j) = dtrk * ( ( max(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) - min(0.0_RP,qflx_anti(k,i,j,ZDIR)) ) * RDZ(k) &
                                + ( max(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) - min(0.0_RP,qflx_anti(k,i,j,XDIR)) ) * RDX(i) &
                                + ( max(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) - min(0.0_RP,qflx_anti(k,i,j,YDIR)) ) * RDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc net outgoing quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k-1,i  ,j  ,ZDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i-1,j  ,XDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j  ,YDIR) )
          call CHECK( __LINE__, qflx_anti(k  ,i  ,j-1,YDIR) )
#endif
          pjmns(k,i,j) = dtrk * ( ( max(0.0_RP,qflx_anti(k,i,j,ZDIR)) - min(0.0_RP,qflx_anti(k-1,i  ,j  ,ZDIR)) ) * RDZ(k) &
                                + ( max(0.0_RP,qflx_anti(k,i,j,XDIR)) - min(0.0_RP,qflx_anti(k  ,i-1,j  ,XDIR)) ) * RDX(i) &
                                + ( max(0.0_RP,qflx_anti(k,i,j,YDIR)) - min(0.0_RP,qflx_anti(k  ,i  ,j-1,YDIR)) ) * RDY(j) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- calc allowable range or quantity change by antidiffusive flux
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
#ifdef DEBUG
          call CHECK( __LINE__, phi_in(k  ,i  ,j  ) )
          call CHECK( __LINE__, phi_in(k-1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(k+1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(k  ,i-1,j  ) )
          call CHECK( __LINE__, phi_in(k  ,i+1,j  ) )
          call CHECK( __LINE__, phi_in(k  ,i  ,j+1) )
          call CHECK( __LINE__, phi_in(k  ,i  ,j-1) )
          call CHECK( __LINE__, phi_lo(k  ,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(k-1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(k+1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(k  ,i-1,j  ) )
          call CHECK( __LINE__, phi_lo(k  ,i+1,j  ) )
          call CHECK( __LINE__, phi_lo(k  ,i  ,j+1) )
          call CHECK( __LINE__, phi_lo(k  ,i  ,j-1) )
#endif
          qjpls(k,i,j) = max( phi_in(k  ,i  ,j  ), &
                              phi_in(k+1,i  ,j  ), &
                              phi_in(k-1,i  ,j  ), &
                              phi_in(k  ,i+1,j  ), &
                              phi_in(k  ,i-1,j  ), &
                              phi_in(k  ,i  ,j+1), &
                              phi_in(k  ,i  ,j-1), &
                              phi_lo(k  ,i  ,j  ), &
                              phi_lo(k+1,i  ,j  ), &
                              phi_lo(k-1,i  ,j  ), &
                              phi_lo(k  ,i+1,j  ), &
                              phi_lo(k  ,i-1,j  ), &
                              phi_lo(k  ,i  ,j+1), &
                              phi_lo(k  ,i  ,j-1) ) &
                       - phi_lo(k,i,j)
          qjmns(k,i,j) = phi_lo(k,i,j) &
                       - min( phi_in(k  ,i  ,j  ), &
                              phi_in(k+1,i  ,j  ), &
                              phi_in(k-1,i  ,j  ), &
                              phi_in(k  ,i-1,j  ), &
                              phi_in(k  ,i+1,j  ), &
                              phi_in(k  ,i  ,j+1), &
                              phi_in(k  ,i  ,j-1), &
                              phi_lo(k  ,i  ,j  ), &
                              phi_lo(k+1,i  ,j  ), &
                              phi_lo(k-1,i  ,j  ), &
                              phi_lo(k  ,i-1,j  ), &
                              phi_lo(k  ,i+1,j  ), &
                              phi_lo(k  ,i  ,j+1), &
                              phi_lo(k  ,i  ,j-1) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, phi_in(KS  ,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KS+1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KS  ,i-1,j  ) )
          call CHECK( __LINE__, phi_in(KS  ,i+1,j  ) )
          call CHECK( __LINE__, phi_in(KS  ,i  ,j+1) )
          call CHECK( __LINE__, phi_in(KS  ,i  ,j-1) )
          call CHECK( __LINE__, phi_lo(KS  ,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KS+1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KS  ,i-1,j  ) )
          call CHECK( __LINE__, phi_lo(KS  ,i+1,j  ) )
          call CHECK( __LINE__, phi_lo(KS  ,i  ,j+1) )
          call CHECK( __LINE__, phi_lo(KS  ,i  ,j-1) )
          call CHECK( __LINE__, phi_in(KE  ,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KE-1,i  ,j  ) )
          call CHECK( __LINE__, phi_in(KE  ,i-1,j  ) )
          call CHECK( __LINE__, phi_in(KE  ,i+1,j  ) )
          call CHECK( __LINE__, phi_in(KE  ,i  ,j+1) )
          call CHECK( __LINE__, phi_in(KE  ,i  ,j-1) )
          call CHECK( __LINE__, phi_lo(KE  ,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KE-1,i  ,j  ) )
          call CHECK( __LINE__, phi_lo(KE  ,i-1,j  ) )
          call CHECK( __LINE__, phi_lo(KE  ,i+1,j  ) )
          call CHECK( __LINE__, phi_lo(KE  ,i  ,j+1) )
          call CHECK( __LINE__, phi_lo(KE  ,i  ,j-1) )
#endif
          qjpls(KS,i,j) = max( phi_in(KS  ,i  ,j  ), &
                               phi_in(KS+1,i  ,j  ), &
                               phi_in(KS  ,i+1,j  ), &
                               phi_in(KS  ,i-1,j  ), &
                               phi_in(KS  ,i  ,j+1), &
                               phi_in(KS  ,i  ,j-1), &
                               phi_lo(KS  ,i  ,j  ), &
                               phi_lo(KS+1,i  ,j  ), &
                               phi_lo(KS  ,i+1,j  ), &
                               phi_lo(KS  ,i-1,j  ), &
                               phi_lo(KS  ,i  ,j+1), &
                               phi_lo(KS  ,i  ,j-1) ) &
                        - phi_lo(KS,i,j)
          qjmns(KS,i,j) = phi_lo(KS,i,j) &
                        - min( phi_in(KS  ,i  ,j  ), &
                               phi_in(KS+1,i  ,j  ), &
                               phi_in(KS  ,i-1,j  ), &
                               phi_in(KS  ,i+1,j  ), &
                               phi_in(KS  ,i  ,j+1), &
                               phi_in(KS  ,i  ,j-1), &
                               phi_lo(KS  ,i  ,j  ), &
                               phi_lo(KS+1,i  ,j  ), &
                               phi_lo(KS  ,i-1,j  ), &
                               phi_lo(KS  ,i+1,j  ), &
                               phi_lo(KS  ,i  ,j+1), &
                               phi_lo(KS  ,i  ,j-1) )
          qjpls(KE,i,j) = max( phi_in(KE  ,i  ,j  ), &
                               phi_in(KE-1,i  ,j  ), &
                               phi_in(KE  ,i+1,j  ), &
                               phi_in(KE  ,i-1,j  ), &
                               phi_in(KE  ,i  ,j+1), &
                               phi_in(KE  ,i  ,j-1), &
                               phi_lo(KE  ,i  ,j  ), &
                               phi_lo(KE-1,i  ,j  ), &
                               phi_lo(KE  ,i+1,j  ), &
                               phi_lo(KE  ,i-1,j  ), &
                               phi_lo(KE  ,i  ,j+1), &
                               phi_lo(KE  ,i  ,j-1) ) &
                        - phi_lo(KE,i,j)
          qjmns(KE,i,j) = phi_lo(KE,i,j) &
                        - min( phi_in(KE  ,i  ,j  ), &
                               phi_in(KE-1,i  ,j  ), &
                               phi_in(KE  ,i-1,j  ), &
                               phi_in(KE  ,i+1,j  ), &
                               phi_in(KE  ,i  ,j+1), &
                               phi_in(KE  ,i  ,j-1), &
                               phi_lo(KE  ,i  ,j  ), &
                               phi_lo(KE-1,i  ,j  ), &
                               phi_lo(KE  ,i-1,j  ), &
                               phi_lo(KE  ,i+1,j  ), &
                               phi_lo(KE  ,i  ,j+1), &
                               phi_lo(KE  ,i  ,j-1) )
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- incoming flux limitation factor [0-1]
       !$omp parallel do private(i,j,k,zerosw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, pjpls(k,i,j) )
          call CHECK( __LINE__, qjpls(k,i,j) )
#endif
          ! if pjpls == 0, zerosw = 1 and rjpls = 0
          zerosw = 0.5_RP - sign( 0.5_RP, pjpls(k,i,j)-EPSILON )
          rjpls(k,i,j) = min( 1.0_RP, qjpls(k,i,j) * ( 1.0_RP-zerosw ) / ( pjpls(k,i,j)-zerosw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       !--- outgoing flux limitation factor [0-1]
       !$omp parallel do private(i,j,k,zerosw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, pjmns(k,i,j) )
          call CHECK( __LINE__, qjmns(k,i,j) )
#endif
          ! if pjmns == 0, zerosw = 1 and rjmns = 0
          zerosw = 0.5_RP - sign( 0.5_RP, pjmns(k,i,j)-EPSILON )
          rjmns(k,i,j) = min( 1.0_RP, qjmns(k,i,j) * ( 1.0_RP-zerosw ) / ( pjmns(k,i,j)-zerosw ) )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    call COMM_vars8( rjpls(:,:,:), 1 )
    call COMM_vars8( rjmns(:,:,:), 2 )
    call COMM_wait ( rjpls(:,:,:), 1 )
    call COMM_wait ( rjmns(:,:,:), 2 )

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !--- update high order flux with antidiffusive flux
       !$omp parallel do private(i,j,k,dirsw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS , KE-1
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,ZDIR) )
          call CHECK( __LINE__, rjpls(k  ,i,j) )
          call CHECK( __LINE__, rjpls(k+1,i,j) )
          call CHECK( __LINE__, rjmns(k  ,i,j) )
          call CHECK( __LINE__, rjmns(k+1,i,j) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,ZDIR) )
          qflx_anti(k,i,j,ZDIR) = qflx_anti(k,i,j,ZDIR) &
                 * ( min( rjpls(k+1,i,j),rjmns(k  ,i,j) ) * (          dirsw ) &
                   + min( rjpls(k  ,i,j),rjmns(k+1,i,j) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif
       !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = IIS, IIE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(KE,i,j,ZDIR) )
          call CHECK( __LINE__, rjpls(KE  ,i,j) )
          call CHECK( __LINE__, rjmns(KE  ,i,j) )
#endif
          qflx_anti(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_anti(KE  ,i,j,ZDIR) = 0.0_RP ! top    boundary
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       if ( IIS == IS ) then
          ijs = IIS-1
       else
          ijs = IIS
       end if
       !$omp parallel do private(i,j,k,dirsw) schedule(static,1) collapse(2)
       do j = JJS, JJE
       do i = ijs, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,XDIR) )
          call CHECK( __LINE__, rjpls(k,i  ,j) )
          call CHECK( __LINE__, rjpls(k,i+1,j) )
          call CHECK( __LINE__, rjmns(k,i  ,j) )
          call CHECK( __LINE__, rjmns(k,i+1,j) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,XDIR) )
          qflx_anti(k,i,j,XDIR) = qflx_anti(k,i,j,XDIR) &
                 * ( min( rjpls(k,i+1,j),rjmns(k,i  ,j) ) * (          dirsw ) &
                   + min( rjpls(k,i  ,j),rjmns(k,i+1,j) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

       if ( JJS == JS ) then
          ijs = JJS-1
       else
          ijs = JJS
       end if
       !$omp parallel do private(i,j,k,dirsw) schedule(static,1) collapse(2)
       do j = ijs, JJE
       do i = IIS, IIE
       do k = KS, KE
#ifdef DEBUG
          call CHECK( __LINE__, qflx_anti(k,i,j,YDIR) )
          call CHECK( __LINE__, rjpls(k,i,j+1) )
          call CHECK( __LINE__, rjpls(k,i,j  ) )
          call CHECK( __LINE__, rjmns(k,i,j  ) )
          call CHECK( __LINE__, rjmns(k,i,j+1) )
#endif
          ! if qflx_anti > 0, dirsw = 1
          dirsw = 0.5_RP + sign( 0.5_RP, qflx_anti(k,i,j,YDIR) )
          qflx_anti(k,i,j,YDIR) = qflx_anti(k,i,j,YDIR) &
                 * ( min( rjpls(k,i,j+1),rjmns(k,i,j  ) ) * (          dirsw ) &
                   + min( rjpls(k,i,j  ),rjmns(k,i,j+1) ) * ( 1.0_RP - dirsw ) &
                   - 1.0_RP )
       enddo
       enddo
       enddo
#ifdef DEBUG
       k = IUNDEF; i = IUNDEF; j = IUNDEF
#endif

    enddo
    enddo

    return
  end subroutine ATMOS_DYN_fct

end module mod_atmos_dyn_common
