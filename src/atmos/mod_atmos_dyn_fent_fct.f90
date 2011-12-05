!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics FENT + FCT
!!
!! @par Description
!!          Dynamical core for Atmospheric process
!!          Full explicit, no terrain + FCT limiter
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!! @li      2011-11-11 (H.Yashiro) [mod] Merged with Y.Miyamoto's
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_dyn
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN
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

  ! time settings
  integer, parameter :: RK = 3 ! order of Runge-Kutta scheme

  ! advection settings
  real(8), parameter :: FACT_N =   7.D0 / 6.D0 !  7/6: fourth, 1: second
  real(8), parameter :: FACT_F = - 1.D0 / 6.D0 ! -1/6: fourth, 0: second

  integer, parameter :: udir = 1
  integer, parameter :: vdir = 2
  integer, parameter :: wdir = 3

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Dynamical Process
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN( dens,   momx,   momy,   momz,   lwpt,   &
                        qtrc,                                   &
                        pres,   velx,   vely,   velz,   temp,   &
                        dens_t, momx_t, momy_t, momz_t, lwpt_t, &
                        qtrc_t                                  )
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_const, only : &
       GRAV   => CONST_GRAV,  &
       Rair   => CONST_Rair,   &
       CPair  => CONST_CPair,  &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       LH0    => CONST_LH0,    &
       Pstd   => CONST_Pstd
    use mod_time, only: &
       TIME_DTSEC_ATMOS_DYN
    use mod_comm, only: &
       COMM_vars
    use mod_grid, only : &
       IA  => GRID_IA, &
       JA  => GRID_JA, &
       KA  => GRID_KA, &
       IS  => GRID_IS, &
       IE  => GRID_IE, &
       JS  => GRID_JS, &
       JE  => GRID_JE, &
       KS  => GRID_KS, &
       KE  => GRID_KE, &
       WS  => GRID_WS, &
       WE  => GRID_WE, &
       RDX => GRID_RDX, &
       RDY => GRID_RDY, &
       RDZ => GRID_RDZ
    use mod_atmos_vars, only: &
       QA => A_QA
    implicit none

    ! prognostic value
    real(8), intent(inout) :: dens(IA,JA,KA)      ! density [kg/m**3]
    real(8), intent(inout) :: momx(IA,JA,KA)      ! momentum (x) [kg/m**3 * m/s]
    real(8), intent(inout) :: momy(IA,JA,KA)      ! momentum (y) [kg/m**3 * m/s]
    real(8), intent(inout) :: momz(IA,JA,KA)      ! momentum (z) [kg/m**3 * m/s]
    real(8), intent(inout) :: lwpt(IA,JA,KA)      ! liquid water potential temperature [K]
    real(8), intent(inout) :: qtrc(IA,JA,KA,QA)   ! tracer mixing ratio   [kg/kg],[1/m3]
    ! diagnostic value
    real(8), intent(inout) :: pres(IA,JA,KA)      ! pressure [Pa]
    real(8), intent(inout) :: velx(IA,JA,KA)      ! velocity (x) [m/s]
    real(8), intent(inout) :: vely(IA,JA,KA)      ! velocity (y) [m/s]
    real(8), intent(inout) :: velz(IA,JA,KA)      ! velocity (z) [m/s]
    real(8), intent(inout) :: temp(IA,JA,KA)      ! Temperature [K]
    ! prognostic tendency
    real(8), intent(out)   :: dens_t(IA,JA,KA)
    real(8), intent(out)   :: momx_t(IA,JA,KA)
    real(8), intent(out)   :: momy_t(IA,JA,KA)
    real(8), intent(out)   :: momz_t(IA,JA,KA)
    real(8), intent(out)   :: lwpt_t(IA,JA,KA)
    real(8), intent(out)   :: qtrc_t(IA,JA,KA,QA)


    ! work
    real(8), allocatable :: dens_w(:,:,:,:)    ! work, density
    real(8), allocatable :: vec   (:,:,:,:)    ! vector quantity = momx + momy + momz
    real(8), allocatable :: vec_w (:,:,:,:)    ! work
    real(8), allocatable :: vec_t (:,:,:,:)    ! tendency
    real(8), allocatable :: scl   (:,:,:,:)    ! scalar quantity = lwpt + qtrc
    real(8), allocatable :: scl_w (:,:,:,:)    ! work
    real(8), allocatable :: scl_t (:,:,:,:)    ! tendency
    real(8), allocatable :: diag  (:,:,:,:)    ! work, diagnostics
    ! divergence
    real(8), allocatable :: ddiv  (:,:,:)      ! density divergence [kg/m**3/s]
    real(8), allocatable :: qdiv  (:,:,:)      ! divergence for any quantity
    ! mass flux
    real(8), allocatable :: mflx_lo  (:,:,:,:) ! rho * vel(x,y,z) @ (u,v,w)-face low  order
    real(8), allocatable :: mflx_hi  (:,:,:,:) ! rho * vel(x,y,z) @ (u,v,w)-face high order
    real(8), allocatable :: mflx_anti(:,:,:,:) ! rho * vel(x,y,z) @ (u,v,w)-face antidiffusive
    real(8), allocatable :: qflx_lo  (:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face low  order
    real(8), allocatable :: qflx_hi  (:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face high order
    real(8), allocatable :: qflx_anti(:,:,:,:) ! rho * vel(x,y,z) * phi @ (u,v,w)-face antidiffusive
    ! flux correction
    real(8), allocatable :: rjpls    (:,:,:,:) ! plus  in (x,y,z)-direction
    real(8), allocatable :: rjmns    (:,:,:,:) ! minus in (x,y,z)-direction

    real(8) :: pjpls, pjmns, pjmax, pjmin, var_l
    real(8) :: midvel

    real(8) :: fp, dfdp, dp
    real(8) :: eps = 1.D-10
    integer, parameter :: itmax = 20 ! max itelation cycle for pressure assumption

    real(8) :: dtrk
    integer :: i, j, k, v, d, rko, it
    !---------------------------------------------------------------------------

    allocate( dens_w(IA,JA,KA,1)         )
    allocate( vec   (IA,JA,KA,udir:wdir) )
    allocate( vec_w (IA,JA,KA,udir:wdir) )
    allocate( vec_t (IA,JA,KA,udir:wdir) )
    allocate( scl   (IA,JA,KA,1+QA)      )
    allocate( scl_w (IA,JA,KA,1+QA)      )
    allocate( scl_t (IA,JA,KA,1+QA)      )
    allocate( diag  (IA,JA,KA,4)         )

    allocate( ddiv  (IA,JA,KA) )
    allocate( qdiv  (IA,JA,KA) )

    allocate( mflx_lo  (IA,JA,KA,udir:wdir) )
    allocate( mflx_hi  (IA,JA,KA,udir:wdir) )
    allocate( mflx_anti(IA,JA,KA,udir:wdir) )
    allocate( qflx_lo  (IA,JA,KA,udir:wdir) )
    allocate( qflx_hi  (IA,JA,KA,udir:wdir) )
    allocate( qflx_anti(IA,JA,KA,udir:wdir) )

    allocate( rjpls    (IA,JA,KA,udir:wdir) )
    allocate( rjmns    (IA,JA,KA,udir:wdir) )

    ! copy to work from global vars
    dens_w(:,:,:,1)    = dens(:,:,:)

    vec  (:,:,:,udir) = momx(:,:,:)
    vec  (:,:,:,vdir) = momy(:,:,:)
    vec  (:,:,:,wdir) = momz(:,:,:)
    vec_w(:,:,:,udir) = momx(:,:,:)
    vec_w(:,:,:,vdir) = momy(:,:,:)
    vec_w(:,:,:,wdir) = momz(:,:,:)

    scl  (:,:,:,1)    = lwpt(:,:,:)
    scl_w(:,:,:,1)    = lwpt(:,:,:)
    do v = 1, QA
       scl  (:,:,:,v+1) = qtrc(:,:,:,v)
       scl_w(:,:,:,v+1) = qtrc(:,:,:,v)
    enddo
    diag(:,:,:,4) = pres(:,:,:) ! first guess

    do rko = 1, RK
       dtrk = TIME_DTSEC_ATMOS_DYN / (RK - rko + 1)

       !--- reset tendency
       dens_t(:,:,:)   = 0.D0
       vec_t (:,:,:,:) = 0.D0
       scl_t (:,:,:,:) = 0.D0


       !##### continuity equation #####
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          mflx_lo(i,j,k,udir) = 0.5D0 * ( dens_w(i+1,j,k,1)+dens_w(i,j,k,1) ) * velx(i,j,k)  &
                              - 0.5D0 * ( dens_w(i+1,j,k,1)-dens_w(i,j,k,1) ) * abs(velx(i,j,k))
          mflx_hi(i,j,k,udir) = ( 0.5D0 * ( dens_w(i+1,j,k,1)+dens_w(i  ,j,k,1) ) * FACT_N   &
                                + 0.5D0 * ( dens_w(i+2,j,k,1)+dens_w(i-1,j,k,1) ) * FACT_F ) &
                                * velx(i,j,k)
       enddo
       enddo
       enddo

       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          mflx_lo(i,j,k,vdir) = 0.5D0 * ( dens_w(i,j+1,k,1)+dens_w(i,j,k,1) ) * vely(i,j,k)  &
                              - 0.5D0 * ( dens_w(i,j+1,k,1)-dens_w(i,j,k,1) ) * abs(vely(i,j,k))
          mflx_hi(i,j,k,vdir) = ( 0.5D0 * ( dens_w(i,j+1,k,1)+dens_w(i,j  ,k,1) ) * FACT_N   &
                                + 0.5D0 * ( dens_w(i,j+2,k,1)+dens_w(i,j-1,k,1) ) * FACT_F ) &
                                * vely(i,j,k)
       enddo
       enddo
       enddo

       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          mflx_lo(i,j,k,wdir) = 0.5D0 * ( dens_w(i,j,k+1,1)+dens_w(i,j,k,1) ) * velz(i,j,k)  &
                              - 0.5D0 * ( dens_w(i,j,k+1,1)-dens_w(i,j,k,1) ) * abs(velz(i,j,k))
          mflx_hi(i,j,k,wdir) = ( 0.5D0 * ( dens_w(i,j,k+1,1)+dens_w(i,j,k  ,1) ) * FACT_N   &
                                + 0.5D0 * ( dens_w(i,j,k+2,1)+dens_w(i,j,k-1,1) ) * FACT_F ) &
                              * velz(i,j,k)
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          mflx_lo(i,j,WS  ,wdir) = 0.D0
          mflx_lo(i,j,WS+1,wdir) = 0.5D0 * ( dens_w(i,j,WS+2,1)+dens_w(i,j,WS+1,1) ) * velz(i,j,WS+1)
          mflx_lo(i,j,WE-1,wdir) = 0.5D0 * ( dens_w(i,j,WE  ,1)+dens_w(i,j,WE-1,1) ) * velz(i,j,WE-1)
          mflx_lo(i,j,WE  ,wdir) = 0.D0

          mflx_hi(i,j,WS  ,wdir) = mflx_lo(i,j,WS  ,wdir) ! bottom boundary
          mflx_hi(i,j,WS+1,wdir) = mflx_lo(i,j,WS+1,wdir) ! just above the bottom boundary
          mflx_hi(i,j,WE-1,wdir) = mflx_lo(i,j,WE-1,wdir) ! just below the top boundary
          mflx_hi(i,j,WE  ,wdir) = mflx_lo(i,j,WE  ,wdir) ! top boundary
       enddo
       enddo

       ! -- flux-divergence
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          ddiv(i,j,k) = ( mflx_lo(i,j,k,udir)-mflx_lo(i-1,j,  k  ,udir) ) * RDX &
                      + ( mflx_lo(i,j,k,vdir)-mflx_lo(i,  j-1,k  ,vdir) ) * RDY &
                      + ( mflx_lo(i,j,k,wdir)-mflx_lo(i,  j,  k-1,wdir) ) * RDZ
       enddo
       enddo
       enddo

       mflx_anti(:,:,:,:) = mflx_hi(:,:,:,:) - mflx_lo(:,:,:,:)

       rjpls(:,:,:,:) = 0.D0
       rjmns(:,:,:,:) = 0.D0

       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var_l = dens(i,j,k) - dtrk * ddiv(i,j,k)

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( dens_w(i  ,j  ,k  ,1), &
                       dens_w(i+1,j  ,k  ,1), &
                       dens_w(i-1,j  ,k  ,1), &
                       dens_w(i  ,j+1,k  ,1), &
                       dens_w(i  ,j-1,k  ,1), &
                       dens_w(i  ,j  ,k+1,1), &
                       dens_w(i  ,j  ,k-1,1)  )

          pjmin = min( dens_w(i  ,j  ,k  ,1), &
                       dens_w(i+1,j  ,k  ,1), &
                       dens_w(i-1,j  ,k  ,1), &
                       dens_w(i  ,j+1,k  ,1), &
                       dens_w(i  ,j-1,k  ,1), &
                       dens_w(i  ,j  ,k+1,1), &
                       dens_w(i  ,j  ,k-1,1)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, mflx_anti(i-1,j  ,k  ,udir) ) - min( 0.D0, mflx_anti(i  ,j  ,k  ,udir) ) &
                + max( 0.D0, mflx_anti(i  ,j-1,k  ,vdir) ) - min( 0.D0, mflx_anti(i  ,j  ,k  ,vdir) ) &
                + max( 0.D0, mflx_anti(i  ,j  ,k-1,wdir) ) - min( 0.D0, mflx_anti(i  ,j  ,k  ,wdir) )
          pjmns = max( 0.D0, mflx_anti(i  ,j  ,k  ,udir) ) - min( 0.D0, mflx_anti(i-1,j  ,k  ,udir) ) &
                + max( 0.D0, mflx_anti(i  ,j  ,k  ,vdir) ) - min( 0.D0, mflx_anti(i  ,j-1,k  ,vdir) ) &
                + max( 0.D0, mflx_anti(i  ,j  ,k  ,wdir) ) - min( 0.D0, mflx_anti(i  ,j  ,k-1,wdir) )

          ! --- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             rjpls(i,j,k,udir) = (pjmax-var_l) / pjpls * abs(velx(i,j,k)+velx(i-1,j  ,k  )) * 0.5D0
             rjpls(i,j,k,vdir) = (pjmax-var_l) / pjpls * abs(vely(i,j,k)+vely(i  ,j-1,k  )) * 0.5D0
             rjpls(i,j,k,wdir) = (pjmax-var_l) / pjpls * abs(velz(i,j,k)+velz(i  ,j  ,k-1)) * 0.5D0
          endif

          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             rjmns(i,j,k,udir) = (var_l-pjmin) / pjmns * abs(velx(i,j,k)+velx(i-1,j  ,k  )) * 0.5D0
             rjmns(i,j,k,vdir) = (var_l-pjmin) / pjmns * abs(vely(i,j,k)+vely(i  ,j-1,k  )) * 0.5D0
             rjmns(i,j,k,wdir) = (var_l-pjmin) / pjmns * abs(velz(i,j,k)+velz(i  ,j  ,k-1)) * 0.5D0
          endif

       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = KS-1, KE
       do j = JS-1, JE
       do i = IS-1, IE
          if ( mflx_anti(i,j,k,udir) >= 0 ) then
             mflx_anti(i,j,k,udir) = mflx_anti(i,j,k,udir) * min( rjpls(i+1,j,k,udir), rjmns(i  ,j,k,udir), 1.D0 )
          else ! mflx_ua(i,j,k) < 0
             mflx_anti(i,j,k,udir) = mflx_anti(i,j,k,udir) * min( rjpls(i  ,j,k,udir), rjmns(i+1,j,k,udir), 1.D0 )
          endif

          if ( mflx_anti(i,j,k,vdir) >= 0 ) then
             mflx_anti(i,j,k,vdir) = mflx_anti(i,j,k,vdir) * min( rjpls(i,j+1,k,vdir), rjmns(i,j  ,k,vdir), 1.D0 )
          else ! mflx_va(i,j,k) < 0
             mflx_anti(i,j,k,vdir) = mflx_anti(i,j,k,vdir) * min( rjpls(i,j  ,k,vdir), rjmns(i,j+1,k,vdir), 1.D0 )
          endif

          if ( mflx_anti(i,j,k,wdir) >= 0 ) then
             mflx_anti(i,j,k,wdir) = mflx_anti(i,j,k,wdir) * min( rjpls(i,j,k+1,wdir), rjmns(i,j,k  ,wdir), 1.D0 )
          else ! mflx_wa(i,j,k) < 0
             mflx_anti(i,j,k,wdir) = mflx_anti(i,j,k,wdir) * min( rjpls(i,j,k  ,wdir), rjmns(i,j,k+1,wdir), 1.D0 )
          endif
       enddo
       enddo
       enddo

       ! -- modify flux-divergence with antidiffusive fluxes
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          ddiv(i,j,k) = ddiv(i,j,k) + ( mflx_anti(i,j,k,udir)-mflx_anti(i-1,j  ,k  ,udir) ) * RDX &
                                    + ( mflx_anti(i,j,k,vdir)-mflx_anti(i  ,j-1,k  ,vdir) ) * RDY &
                                    + ( mflx_anti(i,j,k,wdir)-mflx_anti(i  ,j  ,k-1,wdir) ) * RDZ
          dens_t(i,j,k) = - ddiv(i,j,k)
       enddo
       enddo
       enddo


       !##### momentum equation #####

       ! -- make mass fluxes ( momentum(x) * vel )

       ! at (x, y, layer)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE+1
          midvel = 0.5D0 * ( velx(i,j,k)+velx(i-1,j,k) ) ! on face

          qflx_lo(i,j,k,udir) = 0.5D0 * ( vec_w(i,j,k,udir)+vec_w(i-1,j,k,udir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i,j,k,udir)-vec_w(i-1,j,k,udir) ) * abs(midvel)
          qflx_hi(i,j,k,udir) = ( 0.5D0 * ( vec_w(i  ,j,k,udir)+vec_w(i-1,j,k,udir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+1,j,k,udir)+vec_w(i-2,j,k,udir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do k = KS,   KE
       do j = JS-1, JE
       do i = IS,   IE
          midvel = 0.5D0 * ( vely(i+1,j,k)+vely(i,j,k) ) ! on face

          qflx_lo(i,j,k,vdir) = 0.5D0 * ( vec_w(i,j+1,k,udir)+vec_w(i,j,k,udir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i,j+1,k,udir)-vec_w(i,j,k,udir) ) * abs(midvel)
          qflx_hi(i,j,k,vdir) = ( 0.5D0 * ( vec_w(i,j+1,k,udir)+vec_w(i,j  ,k,udir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i,j+2,k,udir)+vec_w(i,j-1,k,udir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          midvel = 0.5D0 * ( velz(i+1,j,k)+velz(i,j,k) ) ! on face

          qflx_lo(i,j,k,wdir) = 0.5D0 * ( vec_w(i,j,k+1,udir)+vec_w(i,j,k,udir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i,j,k+1,udir)-vec_w(i,j,k,udir) ) * abs(midvel)
          qflx_hi(i,j,k,wdir) = ( 0.5D0 * ( vec_w(i,j,k+1,udir)+vec_w(i,j,k  ,udir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i,j,k+2,udir)+vec_w(i,j,k-1,udir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_lo(i,j,WS  ,wdir) = 0.D0
          qflx_lo(i,j,WS+1,wdir) = 0.25D0 * ( vec_w(i,j,WS+2,udir)+vec_w(i,j,WS+1,udir) ) * ( velz(i+1,j,WS+2)+velz(i,j,WS+1) )
          qflx_lo(i,j,WE-1,wdir) = 0.25D0 * ( vec_w(i,j,WE  ,udir)+vec_w(i,j,WE-1,udir) ) * ( velz(i+1,j,WE  )+velz(i,j,WE-1) )
          qflx_lo(i,j,WE  ,wdir) = 0.D0

          qflx_hi(i,j,WS  ,wdir) = qflx_lo(i,j,WS  ,wdir) ! bottom boundary
          qflx_hi(i,j,WS+1,wdir) = qflx_lo(i,j,WS+1,wdir) ! just above the bottom boundary
          qflx_hi(i,j,WE-1,wdir) = qflx_lo(i,j,WE-1,wdir) ! just below the top boundary
          qflx_hi(i,j,WE  ,wdir) = qflx_lo(i,j,WE  ,wdir) ! top boundary
       enddo
       enddo

       ! -- update flux-divergence with the monotone scheme
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          qdiv(i,j,k) = ( qflx_lo(i+1,j,k,udir)-qflx_lo(i,j,  k  ,udir) ) * RDX &
                      + ( qflx_lo(i  ,j,k,vdir)-qflx_lo(i,j-1,k  ,vdir) ) * RDY &
                      + ( qflx_lo(i  ,j,k,wdir)-qflx_lo(i,j,  k-1,wdir) ) * RDZ
       enddo
       enddo
       enddo

       qflx_anti(:,:,:,:) = qflx_hi(:,:,:,:) - qflx_lo(:,:,:,:)

       rjpls(:,:,:,:) = 0.D0
       rjmns(:,:,:,:) = 0.D0

       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var_l = vec(i,j,k,udir) - dtrk * ( (pres(i+1,j,k)-pres(i,j,k)) * RDX + qdiv(i,j,k) )

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( vec_w(i  ,j  ,k  ,udir), &
                       vec_w(i+1,j  ,k  ,udir), &
                       vec_w(i-1,j  ,k  ,udir), &
                       vec_w(i  ,j+1,k  ,udir), &
                       vec_w(i  ,j-1,k  ,udir), &
                       vec_w(i  ,j  ,k+1,udir), &
                       vec_w(i  ,j  ,k-1,udir)  )

          pjmin = min( vec_w(i  ,j  ,k  ,udir), &
                       vec_w(i+1,j  ,k  ,udir), &
                       vec_w(i-1,j  ,k  ,udir), &
                       vec_w(i  ,j+1,k  ,udir), &
                       vec_w(i  ,j-1,k  ,udir), &
                       vec_w(i  ,j  ,k+1,udir), &
                       vec_w(i  ,j  ,k-1,udir)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i+1,j  ,k  ,udir) ) &
                + max( 0.D0, qflx_anti(i  ,j-1,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k-1,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) )
          pjmns = max( 0.D0, qflx_anti(i+1,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,vdir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,wdir) )

          ! --- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             rjpls(i,j,k,udir) = (pjmax-var_l) / pjpls * abs(velx(i,j,k))

             midvel = 0.25D0 * ( vely(i+1,j  ,k  ) &
                               + vely(i  ,j  ,k  ) &
                               + vely(i+1,j-1,k  ) &
                               + vely(i  ,j-1,k  ) )
             rjpls(i,j,k,vdir) = (pjmax-var_l) / pjpls * abs(midvel)

             midvel = 0.25D0 * ( velz(i+1,j  ,k  ) &
                               + velz(i  ,j  ,k  ) &
                               + velz(i+1,j  ,k-1) &
                               + velz(i  ,j  ,k-1) )
             rjpls(i,j,k,wdir) = (pjmax-var_l) / pjpls * abs(midvel)
          endif

          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             rjmns(i,j,k,udir) = (var_l-pjmin) / pjmns * abs(velx(i,j,k))

             midvel = 0.25D0 * ( vely(i+1,j  ,k  ) &
                               + vely(i  ,j  ,k  ) &
                               + vely(i+1,j-1,k  ) &
                               + vely(i  ,j-1,k  ) )
             rjmns(i,j,k,vdir) = (var_l-pjmin) / pjmns * abs(midvel)

             midvel = 0.25D0 * ( velz(i+1,j  ,k  ) &
                               + velz(i  ,j  ,k  ) &
                               + velz(i+1,j  ,k-1) &
                               + velz(i  ,j  ,k-1) )
             rjmns(i,j,k,wdir) = (var_l-pjmin) / pjmns * abs(midvel)
          endif

       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = KS-1, KE
       do j = JS-1, JE
       do i = IS  , IE+1
          if ( qflx_anti(i,j,k,udir) >= 0 ) then
             qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i  ,j,k,udir), rjmns(i-1,j,k,udir), 1.D0 )
          else ! qflx_anti(i,j,k,udir) < 0
             qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i-1,j,k,udir), rjmns(i  ,j,k,udir), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,vdir) >= 0 ) then
             qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j+1,k,vdir), rjmns(i,j  ,k,vdir), 1.D0 )
          else ! qflx_anti(i,j,k,vdir) < 0
             qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j  ,k,vdir), rjmns(i,j+1,k,vdir), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,wdir) >= 0 ) then
             qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k+1,wdir), rjmns(i,j,k  ,wdir), 1.D0 )
          else ! qflx_anti(i,j,k,wdir) < 0
             qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k  ,wdir), rjmns(i,j,k+1,wdir), 1.D0 )
          endif
       enddo
       enddo
       enddo

       ! -- modify flux-divergence with antidiffusive fluxes
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          qdiv(i,j,k) = qdiv(i,j,k) + ( qflx_anti(i+1,j,k,udir)-qflx_anti(i,j  ,k  ,udir) ) * RDX &
                                    + ( qflx_anti(i  ,j,k,vdir)-qflx_anti(i,j-1,k  ,vdir) ) * RDY &
                                    + ( qflx_anti(i  ,j,k,wdir)-qflx_anti(i,j  ,k-1,wdir) ) * RDZ
          vec_t(i,j,k,udir) = - qdiv(i,j,k) - ( pres(i+1,j,k)-pres(i,j,k) ) * RDX 
       enddo
       enddo
       enddo

       ! -- make mass fluxes ( momentum(y) * vel )

       ! at (u, v, layer)
       do k = KS,   KE
       do j = JS,   JE
       do i = IS-1, IE
          midvel = 0.5D0 * ( velx(i,j+1,k)+velx(i,j,k) ) ! on face

          qflx_lo(i,j,k,udir) = 0.5D0 * ( vec_w(i+1,j,k,vdir)+vec_w(i,j,k,vdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i+1,j,k,vdir)-vec_w(i,j,k,vdir) ) * abs(midvel)
          qflx_hi(i,j,k,udir) = ( 0.5D0 * ( vec_w(i+1,j,k,vdir)+vec_w(i  ,j,k,vdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+2,j,k,vdir)+vec_w(i-1,j,k,vdir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo
       ! at (x, y, layer)
       do k = KS, KE
       do j = JS, JE+1
       do i = IS, IE
          midvel = 0.5D0 * ( vely(i,j,k)+vely(i,j-1,k) ) ! on face

          qflx_lo(i,j,k,vdir) = 0.5D0 * ( vec_w(i,j,k,vdir)+vec_w(i,j-1,k,vdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i,j,k,vdir)-vec_w(i,j-1,k,vdir) ) * abs(midvel)
          qflx_hi(i,j,k,vdir) = ( 0.5D0 * ( vec_w(i,j  ,k,vdir)+vec_w(i,j-1,k,vdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i,j+1,k,vdir)+vec_w(i,j-2,k,vdir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do k = WS+2, WE-2
       do j = JS,   JE
       do i = IS,   IE
          midvel = 0.5D0 * ( velz(i,j+1,k)+velz(i,j,k) ) ! on face

          qflx_lo(i,j,k,wdir) = 0.5D0 * ( vec_w(i,j,k+1,vdir)+vec_w(i,j,k,vdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i,j,k+1,vdir)-vec_w(i,j,k,vdir) ) * abs(midvel)
          qflx_hi(i,j,k,wdir) = ( 0.5D0 * ( vec_w(i,j,k+1,vdir)+vec_w(i,j,k  ,vdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i,j,k+2,vdir)+vec_w(i,j,k-1,vdir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          qflx_lo(i,j,WS  ,wdir) = 0.D0
          qflx_lo(i,j,WS+1,wdir) = 0.25D0 * ( vec_w(i,j,WS+2,vdir)+vec_w(i,j,WS+1,vdir) ) * ( velz(i,j+1,WS+2)+velz(i,j,WS+1) )
          qflx_lo(i,j,WE-1,wdir) = 0.25D0 * ( vec_w(i,j,WE  ,vdir)+vec_w(i,j,WE-1,vdir) ) * ( velz(i,j+1,WE  )+velz(i,j,WE-1) )
          qflx_lo(i,j,WE  ,wdir) = 0.D0

          qflx_hi(i,j,WS  ,wdir) = qflx_lo(i,j,WS  ,wdir) ! bottom boundary
          qflx_hi(i,j,WS+1,wdir) = qflx_lo(i,j,WS+1,wdir) ! just above the bottom boundary
          qflx_hi(i,j,WE-1,wdir) = qflx_lo(i,j,WE-1,wdir) ! just below the top boundary
          qflx_hi(i,j,WE  ,wdir) = qflx_lo(i,j,WE  ,wdir) ! top boundary
       enddo
       enddo

       ! -- update flux-divergence with the monotone scheme
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          qdiv(i,j,k) = ( qflx_lo(i,j  ,k,udir)-qflx_lo(i-1,j,k  ,udir) ) * RDX &
                      + ( qflx_lo(i,j+1,k,vdir)-qflx_lo(i  ,j,k  ,vdir) ) * RDY &
                      + ( qflx_lo(i,j  ,k,wdir)-qflx_lo(i  ,j,k-1,wdir) ) * RDZ
       enddo
       enddo
       enddo

       qflx_anti(:,:,:,:) = qflx_hi(:,:,:,:) - qflx_lo(:,:,:,:)

       rjpls(:,:,:,:) = 0.D0
       rjmns(:,:,:,:) = 0.D0

       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          var_l = vec(i,j,k,vdir) - dtrk * ( (pres(i,j+1,k)-pres(i,j,k)) * RDY + qdiv(i,j,k) )

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( vec_w(i  ,j  ,k  ,vdir), &
                       vec_w(i+1,j  ,k  ,vdir), &
                       vec_w(i-1,j  ,k  ,vdir), &
                       vec_w(i  ,j+1,k  ,vdir), &
                       vec_w(i  ,j-1,k  ,vdir), &
                       vec_w(i  ,j  ,k+1,vdir), &
                       vec_w(i  ,j  ,k-1,vdir)  )

          pjmin = min( vec_w(i  ,j  ,k  ,vdir), &
                       vec_w(i+1,j  ,k  ,vdir), &
                       vec_w(i-1,j  ,k  ,vdir), &
                       vec_w(i  ,j+1,k  ,vdir), &
                       vec_w(i  ,j-1,k  ,vdir), &
                       vec_w(i  ,j  ,k+1,vdir), &
                       vec_w(i  ,j  ,k-1,vdir)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j+1,k  ,vdir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k-1,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) )
          pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,udir) ) &
                + max( 0.D0, qflx_anti(i  ,j+1,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,wdir) )

          ! --- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             midvel = 0.25D0 * ( velx(i  ,j+1,k  ) &
                               + velx(i  ,j  ,k  ) &
                               + velx(i-1,j+1,k  ) &
                               + velx(i-1,j  ,k  ) )
             rjpls(i,j,k,udir) = (pjmax-var_l) / pjpls * abs(midvel)

             rjpls(i,j,k,vdir) = (pjmax-var_l) / pjpls * abs(vely(i,j,k))

             midvel = 0.25D0 * ( velz(i  ,j+1,k  ) &
                               + velz(i  ,j  ,k  ) &
                               + velz(i  ,j+1,k-1) &
                               + velz(i  ,j  ,k-1) )
             rjpls(i,j,k,wdir) = (pjmax-var_l) / pjpls * abs(midvel)
          endif

          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             midvel = 0.25D0 * ( velx(i  ,j+1,k  ) &
                               + velx(i  ,j  ,k  ) &
                               + velx(i-1,j+1,k  ) &
                               + velx(i-1,j  ,k  ) )
             rjmns(i,j,k,udir) = (var_l-pjmin) / pjmns * abs(midvel)

             rjmns(i,j,k,vdir) = (var_l-pjmin) / pjmns * abs(vely(i,j,k))

             midvel = 0.25D0 * ( velz(i  ,j+1,k  ) &
                               + velz(i  ,j  ,k  ) &
                               + velz(i  ,j+1,k-1) &
                               + velz(i  ,j  ,k-1) )
             rjmns(i,j,k,wdir) = (var_l-pjmin) / pjmns * abs(midvel)
          endif

       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = KS-1, KE
       do j = JS,   JE+1
       do i = IS-1, IE
          if ( qflx_anti(i,j,k,udir) >= 0 ) then
             qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i+1,j,k,udir), rjmns(i  ,j,k,udir), 1.D0 )
          else ! qflx_anti(i,j,k,udir) < 0
             qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i  ,j,k,udir), rjmns(i+1,j,k,udir), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,vdir) >= 0 ) then
             qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j  ,k,vdir), rjmns(i,j-1,k,vdir), 1.D0 )
          else ! qflx_anti(i,j,k,vdir) < 0
             qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j-1,k,vdir), rjmns(i,j  ,k,vdir), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,wdir) >= 0 ) then
             qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k+1,wdir), rjmns(i,j,k  ,wdir), 1.D0 )
          else ! qflx_anti(i,j,k,wdir) < 0
             qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k  ,wdir), rjmns(i,j,k+1,wdir), 1.D0 )
          endif
       enddo
       enddo
       enddo

       ! -- modify flux-divergence with antidiffusive fluxes
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          qdiv(i,j,k) = qdiv(i,j,k) + ( qflx_anti(i,j  ,k,udir)-qflx_anti(i-1,j,k  ,udir) ) * RDX &
                                    + ( qflx_anti(i,j+1,k,vdir)-qflx_anti(i  ,j,k  ,vdir) ) * RDY &
                                    + ( qflx_anti(i,j  ,k,wdir)-qflx_anti(i  ,j,k-1,wdir) ) * RDZ
          vec_t(i,j,k,vdir) = - qdiv(i,j,k) - ( pres(i,j+1,k)-pres(i,j,k) ) * RDY 
       enddo
       enddo
       enddo

       ! -- make mass fluxes ( momentum(z) * vel )

       k = WS ! bottom boundary
       do j = JS,   JE
       do i = IS-1, IE
          midvel = velx(i,j,KS)

          qflx_lo(i,j,k,udir) = 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i,j,k,wdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i+1,j,k,wdir)-vec_w(i,j,k,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,udir) = ( 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i  ,j,k,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+2,j,k,wdir)+vec_w(i-1,j,k,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       k = WE ! top boundary
       do j = JS,   JE
       do i = IS-1, IE
          midvel = velx(i,j,KE)

          qflx_lo(i,j,k,udir) = 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i,j,k,wdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i+1,j,k,wdir)-vec_w(i,j,k,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,udir) = ( 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i  ,j,k,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+2,j,k,wdir)+vec_w(i-1,j,k,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       ! at (u, y, interface)
       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS-1, IE
          midvel = 0.5D0 * ( velx(i,j,k+1)+velx(i,j,k) ) ! on face

          qflx_lo(i,j,k,udir) = 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i,j,k,wdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i+1,j,k,wdir)-vec_w(i,j,k,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,udir) = ( 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i  ,j,k,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+2,j,k,wdir)+vec_w(i-1,j,k,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo

       k = WS ! bottom boundary
       do j = JS-1, JE
       do i = IS,   IE
          midvel = vely(i,j,KS)

          qflx_lo(i,j,k,vdir) = 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i,j,k,wdir) ) * midvel     &
                              - 0.5D0 * ( vec_w(i+1,j,k,wdir)-vec_w(i,j,k,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,vdir) = ( 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i  ,j,k,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+2,j,k,wdir)+vec_w(i-1,j,k,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       k = WE ! top boundary
       do j = JS-1, JE
       do i = IS,   IE
          midvel = vely(i,j,KE)

          qflx_lo(i,j,k,vdir) = 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i,j,k,wdir) ) * midvel       &
                              - 0.5D0 * ( vec_w(i+1,j,k,wdir)-vec_w(i,j,k,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,vdir) = ( 0.5D0 * ( vec_w(i+1,j,k,wdir)+vec_w(i  ,j,k,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i+2,j,k,wdir)+vec_w(i-1,j,k,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       ! at (x, v, interface)
       do k = WS+1, WE-1
       do j = JS-1, JE
       do i = IS,   IE
          midvel = 0.5D0 * ( vely(i,j,k+1)+vely(i,j,k) ) ! on face

          qflx_lo(i,j,k,vdir) = 0.5D0 * ( vec_w(i,j+1,k,wdir)+vec_w(i,j,k,wdir) ) * midvel &
                              - 0.5D0 * ( vec_w(i,j+1,k,wdir)-vec_w(i,j,k,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,vdir) = ( 0.5D0 * ( vec_w(i,j+1,k,wdir)+vec_w(i,j  ,k,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i,j+2,k,wdir)+vec_w(i,j-1,k,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       do k = KS+1, KE-1
       do j = JS,   JE
       do i = IS,   IE
          midvel = 0.5D0 * ( velz(i,j,k)+velz(i,j,k-1) ) ! on face

          qflx_lo(i,j,k,wdir) = 0.5D0 * ( vec_w(i,j,k,wdir)+vec_w(i,j,k-1,wdir) ) * midvel &
                              - 0.5D0 * ( vec_w(i,j,k,wdir)-vec_w(i,j,k-1,wdir) ) * abs(midvel)
          qflx_hi(i,j,k,wdir) = ( 0.5D0 * ( vec_w(i,j,k  ,wdir)+vec_w(i,j,k-1,wdir) ) * FACT_N &
                                + 0.5D0 * ( vec_w(i,j,k+1,wdir)+vec_w(i,j,k-2,wdir) ) * FACT_F ) * midvel
       enddo
       enddo
       enddo
       qflx_lo(:,:,KS,wdir) = 0.D0 ! bottom cell center
       qflx_lo(:,:,KE,wdir) = 0.D0 ! top    cell center
       qflx_hi(:,:,KS,wdir) = 0.D0 ! bottom cell center
       qflx_hi(:,:,KE,wdir) = 0.D0 ! top    cell center

       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS,   IE
          qdiv(i,j,k) = ( qflx_lo(i,j,k  ,udir)-qflx_lo(i-1,j  ,k,udir) ) * RDX &
                      + ( qflx_lo(i,j,k  ,vdir)-qflx_lo(i  ,j-1,k,vdir) ) * RDY &
                      + ( qflx_lo(i,j,k+1,wdir)-qflx_lo(i  ,j  ,k,wdir) ) * RDZ
       enddo
       enddo
       enddo

       qflx_anti(:,:,:,:) = qflx_hi(:,:,:,:) - qflx_lo(:,:,:,:)

       rjpls(:,:,:,:) = 0.D0
       rjmns(:,:,:,:) = 0.D0

       do k = WS+1, WE-1
       do j = JS, JE
       do i = IS, IE
          var_l = vec(i,j,k,wdir) - dtrk * ( (pres(i,j,k+1)-pres(i,j,k)) * RDZ + qdiv(i,j,k) &
                                           + (dens_w(i,j,k+1,1)+dens_w(i,j,k,1)) * 0.5D0 * GRAV  )

          ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
          pjmax = max( vec_w(i  ,j  ,k  ,wdir), &
                       vec_w(i+1,j  ,k  ,wdir), &
                       vec_w(i-1,j  ,k  ,wdir), &
                       vec_w(i  ,j+1,k  ,wdir), &
                       vec_w(i  ,j-1,k  ,wdir), &
                       vec_w(i  ,j  ,k+1,wdir), &
                       vec_w(i  ,j  ,k-1,wdir)  )

          pjmin = min( vec_w(i  ,j  ,k  ,wdir), &
                       vec_w(i+1,j  ,k  ,wdir), &
                       vec_w(i-1,j  ,k  ,wdir), &
                       vec_w(i  ,j+1,k  ,wdir), &
                       vec_w(i  ,j-1,k  ,wdir), &
                       vec_w(i  ,j  ,k+1,wdir), &
                       vec_w(i  ,j  ,k-1,wdir)  )

          ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
          pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) &
                + max( 0.D0, qflx_anti(i  ,j-1,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k+1,wdir) )
          pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,udir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,vdir) ) &
                + max( 0.D0, qflx_anti(i  ,j  ,k+1,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) )

          ! --- incoming fluxes at scalar grid points ---
          if ( pjpls > 0 ) then
             midvel = 0.25D0 * ( velx(i  ,j  ,k+1) &
                               + velx(i  ,j  ,k  ) &
                               + velx(i-1,j  ,k+1) &
                               + velx(i-1,j  ,k  ) )
             rjpls(i,j,k,udir) = (pjmax-var_l) / pjpls * abs(midvel)

             midvel = 0.25D0 * ( vely(i  ,j  ,k+1) &
                               + vely(i  ,j  ,k  ) &
                               + vely(i  ,j-1,k+1) &
                               + vely(i  ,j-1,k  ) )
             rjpls(i,j,k,vdir) = (pjmax-var_l) / pjpls * abs(midvel)

             rjpls(i,j,k,wdir) = (pjmax-var_l) / pjpls * abs(velz(i,j,k))
          endif

          ! --- outgoing fluxes at scalar grid points ---
          if ( pjmns > 0 ) then
             midvel = 0.25D0 * ( velx(i  ,j  ,k+1) &
                               + velx(i  ,j  ,k  ) &
                               + velx(i-1,j  ,k+1) &
                               + velx(i-1,j  ,k  ) )
             rjmns(i,j,k,udir) = (var_l-pjmin) / pjmns * abs(midvel)

             midvel = 0.25D0 * ( vely(i  ,j  ,k+1) &
                               + vely(i  ,j  ,k  ) &
                               + vely(i  ,j-1,k+1) &
                               + vely(i  ,j-1,k  ) )
             rjmns(i,j,k,vdir) = (var_l-pjmin) / pjmns * abs(midvel)

             rjmns(i,j,k,wdir) = (var_l-pjmin) / pjmns * abs(velz(i,j,k))
          endif

       enddo
       enddo
       enddo

       ! --- [STEP 7S] limit the antidiffusive flux at velocity grid points ---
       do k = WS+1, WE
       do j = JS-1, JE
       do i = IS-1, IE
          if ( qflx_anti(i,j,k,udir) >= 0 ) then
             qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i+1,j,k,udir), rjmns(i  ,j,k,udir), 1.D0 )
          else ! qflx_anti(i,j,k,udir) < 0
             qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i  ,j,k,udir), rjmns(i+1,j,k,udir), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,vdir) >= 0 ) then
             qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j+1,k,vdir), rjmns(i,j  ,k,vdir), 1.D0 )
          else ! qflx_anti(i,j,k,vdir) < 0
             qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j  ,k,vdir), rjmns(i,j+1,k,vdir), 1.D0 )
          endif

          if ( qflx_anti(i,j,k,wdir) >= 0 ) then
             qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k  ,wdir), rjmns(i,j,k-1,wdir), 1.D0 )
          else ! qflx_anti(i,j,k,wdir) < 0
             qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k-1,wdir), rjmns(i,j,k  ,wdir), 1.D0 )
          endif
       enddo
       enddo
       enddo
        
       ! -- modify flux-divergence with antidiffusive fluxes
       do k = WS+1, WE-1
       do j = JS, JE
       do i = IS, IE
          qdiv(i,j,k) = qdiv(i,j,k) + ( qflx_anti(i,j,k  ,udir)-qflx_anti(i-1,j  ,k,udir) ) * RDX &
                                    + ( qflx_anti(i,j,k  ,vdir)-qflx_anti(i  ,j-1,k,vdir) ) * RDY &
                                    + ( qflx_anti(i,j,k+1,wdir)-qflx_anti(i  ,j  ,k,wdir) ) * RDZ
          vec_t(i,j,k,wdir) = - qdiv(i,j,k) - ( pres(i,j,k+1)-pres(i,j,k) ) * RDZ &
                              - ( dens_w(i,j,k+1,1)+dens_w(i,j,k,1) ) * 0.5D0 * GRAV
       enddo
       enddo
       enddo


       !##### advection of scalar quantity #####
       do v = 1, QA+1 ! lwpt + tracers

          ! -- make mass fluxes ( scalar value * rho vel )

          do k = KS,   KE
          do j = JS,   JE
          do i = IS-1, IE
             qflx_lo(i,j,k,udir) = 0.5D0 * ( scl_w(i+1,j,k,v)+scl_w(i,j,k,v) ) * mflx_lo(i,j,k,udir)       &
                                 - 0.5D0 * ( scl_w(i+1,j,k,v)-scl_w(i,j,k,v) ) * abs(mflx_lo(i,j,k,udir))
             qflx_hi(i,j,k,udir) = ( 0.5D0 * ( scl_w(i+1,j,k,v)+scl_w(i  ,j,k,v) ) * FACT_N   &
                                   + 0.5D0 * ( scl_w(i+2,j,k,v)+scl_w(i-1,j,k,v) ) * FACT_F ) &
                                   * mflx_hi(i,j,k,udir)
          enddo
          enddo
          enddo

          do k = KS,   KE
          do j = JS-1, JE
          do i = IS,   IE
             qflx_lo(i,j,k,vdir) = 0.5D0 * ( scl_w(i,j+1,k,v)+scl_w(i,j,k,v) ) * mflx_lo(i,j,k,vdir)       &
                                 - 0.5D0 * ( scl_w(i,j+1,k,v)-scl_w(i,j,k,v) ) * abs(mflx_lo(i,j,k,vdir))
             qflx_hi(i,j,k,vdir) = ( 0.5D0 * ( scl_w(i,j+1,k,v)+scl_w(i,j  ,k,v) ) * FACT_N   &
                                   + 0.5D0 * ( scl_w(i,j+2,k,v)+scl_w(i,j-1,k,v) ) * FACT_F ) &
                                   * mflx_hi(i,j,k,vdir)
          enddo
          enddo
          enddo

          do k = WS+2, WE-2
          do j = JS,   JE
          do i = IS,   IE
             qflx_lo(i,j,k,wdir) = 0.5D0 * ( scl_w(i,j,k+1,v)+scl_w(i,j,k,v) ) * mflx_lo(i,j,k,wdir)       &
                                 - 0.5D0 * ( scl_w(i,j,k+1,v)-scl_w(i,j,k,v) ) * abs(mflx_lo(i,j,k,wdir))
             qflx_hi(i,j,k,wdir) = ( 0.5D0 * ( scl_w(i,j,k+1,v)+scl_w(i,j,k  ,v) ) * FACT_N   &
                                   + 0.5D0 * ( scl_w(i,j,k+2,v)+scl_w(i,j,k-1,v) ) * FACT_F ) &
                                   * mflx_hi(i,j,k,wdir)
          enddo
          enddo
          enddo

          do j = JS, JE
          do i = IS, IE
             qflx_lo(i,j,WS  ,wdir) = 0.D0
             qflx_lo(i,j,WS+1,wdir) = 0.5D0 * ( scl_w(i,j,WS+2,v)+scl_w(i,j,WS+1,v) ) * mflx_lo(i,j,WS+1,wdir)
             qflx_lo(i,j,WE-1,wdir) = 0.5D0 * ( scl_w(i,j,WE  ,v)+scl_w(i,j,WE-1,v) ) * mflx_lo(i,j,WE-1,wdir)
             qflx_lo(i,j,WE  ,wdir) = 0.D0

             qflx_hi(i,j,WS  ,wdir) = qflx_lo(i,j,WS  ,wdir) ! bottom boundary
             qflx_hi(i,j,WS+1,wdir) = qflx_lo(i,j,WS+1,wdir) ! just above the bottom boundary
             qflx_hi(i,j,WE-1,wdir) = qflx_lo(i,j,WE-1,wdir) ! just below the top boundary
             qflx_hi(i,j,WE  ,wdir) = qflx_lo(i,j,WE  ,wdir) ! top boundary
          enddo
          enddo

          ! -- update flux-divergence with the monotone scheme
          do k = KS, KE
          do j = JS, JE
          do i = IS, IE
             qdiv(i,j,k) = ( qflx_lo(i,j,k,udir)-qflx_lo(i-1,j  ,k  ,udir) ) * RDX &
                         + ( qflx_lo(i,j,k,vdir)-qflx_lo(i  ,j-1,k  ,vdir) ) * RDY &
                         + ( qflx_lo(i,j,k,wdir)-qflx_lo(i  ,j  ,k-1,wdir) ) * RDZ
          enddo
          enddo
          enddo

          qflx_anti(:,:,:,:) = qflx_hi(:,:,:,:) - qflx_lo(:,:,:,:)

          rjpls(:,:,:,:) = 0.D0
          rjmns(:,:,:,:) = 0.D0

          do k = KS, KE
          do j = JS, JE
          do i = IS, IE
             var_l = scl(i,j,k,v) + dtrk * ( scl_w(i,j,k,v)*ddiv(i,j,k) - qdiv(i,j,k) ) / dens_w(i,j,k,1)

             ! --- STEP A: establish allowed extreme max/min in each cell using low order fluxes ---
             pjmax = max( scl_w(i  ,j  ,k  ,v), &
                          scl_w(i+1,j  ,k  ,v), &
                          scl_w(i-1,j  ,k  ,v), &
                          scl_w(i  ,j+1,k  ,v), &
                          scl_w(i  ,j-1,k  ,v), &
                          scl_w(i  ,j  ,k+1,v), &
                          scl_w(i  ,j  ,k-1,v)  )

             pjmin = min( scl_w(i  ,j  ,k  ,v), &
                          scl_w(i+1,j  ,k  ,v), &
                          scl_w(i-1,j  ,k  ,v), &
                          scl_w(i  ,j+1,k  ,v), &
                          scl_w(i  ,j-1,k  ,v), &
                          scl_w(i  ,j  ,k+1,v), &
                          scl_w(i  ,j  ,k-1,v)  )

             ! --- STEP C: compute the total incoming and outgoing antidiffusive fluxes in each cell ---
             pjpls = max( 0.D0, qflx_anti(i-1,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) &
                   + max( 0.D0, qflx_anti(i  ,j-1,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) &
                   + max( 0.D0, qflx_anti(i  ,j  ,k-1,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) )
             pjmns = max( 0.D0, qflx_anti(i  ,j  ,k  ,udir) ) - min( 0.D0, qflx_anti(i-1,j  ,k  ,udir) ) &
                   + max( 0.D0, qflx_anti(i  ,j  ,k  ,vdir) ) - min( 0.D0, qflx_anti(i  ,j-1,k  ,vdir) ) &
                   + max( 0.D0, qflx_anti(i  ,j  ,k  ,wdir) ) - min( 0.D0, qflx_anti(i  ,j  ,k-1,wdir) )

             ! --- incoming fluxes ---
             if ( pjpls > 0 ) then
                rjpls(i,j,k,udir) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,udir)+mflx_lo(i-1,j  ,k  ,udir)) * 0.5D0)
                rjpls(i,j,k,vdir) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,vdir)+mflx_lo(i  ,j-1,k  ,vdir)) * 0.5D0)
                rjpls(i,j,k,wdir) = (pjmax-var_l) / pjpls * abs((mflx_lo(i,j,k,wdir)+mflx_lo(i  ,j  ,k-1,wdir)) * 0.5D0)
             endif

             ! --- outgoing fluxes at scalar grid points ---
             if ( pjmns > 0 ) then
                rjmns(i,j,k,udir) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,udir)+mflx_lo(i-1,j  ,k  ,udir)) * 0.5D0)
                rjmns(i,j,k,vdir) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,vdir)+mflx_lo(i  ,j-1,k  ,vdir)) * 0.5D0)
                rjmns(i,j,k,wdir) = (var_l-pjmin) / pjmns * abs((mflx_lo(i,j,k,wdir)+mflx_lo(i  ,j  ,k-1,wdir)) * 0.5D0)
             endif

          enddo
          enddo
          enddo

          ! --- [STEP 7S] limit the antidiffusive flux ---
          do k = KS-1, KE 
          do j = JS-1, JE 
          do i = IS-1, IE
             if ( qflx_anti(i,j,k,udir) >= 0 ) then
                qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i+1,j,k,udir), rjmns(i  ,j,k,udir), 1.D0 )
             else !if ( mflx_ua(i,j,k) < 0 ) then
                qflx_anti(i,j,k,udir) = qflx_anti(i,j,k,udir) * min( rjpls(i  ,j,k,udir), rjmns(i+1,j,k,udir), 1.D0 )
             endif

             if ( qflx_anti(i,j,k,vdir) >= 0 ) then
                qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j+1,k,vdir), rjmns(i,j  ,k,vdir), 1.D0 )
             else !if ( mflx_va(i,j,k) < 0 ) then
                qflx_anti(i,j,k,vdir) = qflx_anti(i,j,k,vdir) * min( rjpls(i,j  ,k,vdir), rjmns(i,j+1,k,vdir), 1.D0 )
             endif

             if ( qflx_anti(i,j,k,wdir) >= 0 ) then
                qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k+1,wdir), rjmns(i,j,k  ,wdir), 1.D0 )
             else !if ( mflx_wa(i,j,k) < 0 ) then
                qflx_anti(i,j,k,wdir) = qflx_anti(i,j,k,wdir) * min( rjpls(i,j,k  ,wdir), rjmns(i,j,k+1,wdir), 1.D0 )
             endif
          enddo
          enddo
          enddo

          ! -- advection
          do k = KS, KE
          do j = JS, JE
          do i = IS, IE
             qdiv(i,j,k) = qdiv(i,j,k) + ( qflx_anti(i,j,k,udir)-qflx_anti(i-1,j  ,k  ,udir) ) * RDX &
                                       + ( qflx_anti(i,j,k,vdir)-qflx_anti(i  ,j-1,k  ,vdir) ) * RDY &
                                       + ( qflx_anti(i,j,k,wdir)-qflx_anti(i  ,j  ,k-1,wdir) ) * RDZ
             scl_t(i,j,k,v) = ( scl_w(i,j,k,v)*ddiv(i,j,k) - qdiv(i,j,k) ) / dens_w(i,j,k,1)
          enddo
          enddo
          enddo

       enddo ! scalar quantities loop

       !##### time integrations #####
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          dens_w(i,j,k,1) = dens(i,j,k) + dtrk * dens_t(i,j,k)
       enddo
       enddo
       enddo

       do d = udir, wdir
       do k = WS-1, WE+1
       do j = JS, JE
       do i = IS, IE
          vec_w(i,j,k,d) = vec(i,j,k,d) + dtrk * vec_t(i,j,k,d)
       enddo
       enddo
       enddo
       enddo

       do v = 1,  1+QA
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          scl_w(i,j,k,v) = scl(i,j,k,v) + dtrk * scl_t(i,j,k,v)
       enddo
       enddo
       enddo
       enddo

       ! fill IHALO & JHALO
       call COMM_vars( dens_w(:,:,:,:) )
       call COMM_vars( vec_w (:,:,:,:) )
       call COMM_vars( scl_w (:,:,:,:) )

       ! momentum -> velocity
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          diag(i,j,k,1) = 2.D0 * vec_w(i,j,k,udir) / ( dens_w(i+1,j  ,k,1)+dens_w(i,j,k,1) )
          diag(i,j,k,2) = 2.D0 * vec_w(i,j,k,vdir) / ( dens_w(i  ,j+1,k,1)+dens_w(i,j,k,1) )
       enddo
       enddo
       enddo

       do k = WS+1, WE-1
       do j = JS,   JE
       do i = IS,   IE
          diag(i,j,k,3) = 2.D0 * vec_w(i,j,k,wdir) / ( dens_w(i,j,k+1,1)+dens_w(i,j,k,1) )
       enddo
       enddo
       enddo
       diag(:,:,WS,3) = 0.D0 ! bottom boundary
       diag(:,:,WE,3) = 0.D0 ! top    boundary

       ! diagnose pressure, temperature
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
!          diag(i,j,k,4) = Pstd * ( dens_w(i,j,k,1) * scl_w(i,j,k,1) * Rair / Pstd )**CPovCV
          do it = 1, itmax
             fp   = scl_w(i,j,k,1) * ( pres(i,j,k)/Pstd )**RovCP      &
                  + LH0 / CPair * ( scl_w(i,j,k,3) + scl_w(i,j,k,4) ) &
                  - pres(i,j,k) / ( Rair *  dens_w(i,j,k,1) )
             dfdp = RovCP / Pstd * scl_w(i,j,k,1) * ( pres(i,j,k)/Pstd )**(RovCP-1) &
                  - 1.D0 / ( Rair * dens_w(i,j,k,1) )
             dp   = fp / dfdp

             diag(i,j,k,4) = diag(i,j,k,4) - dp
             if ( dp < eps ) exit
          enddo
       enddo
       enddo
       enddo

       call COMM_vars( diag(:,:,:,:) )

       velx(:,:,:) = diag(:,:,:,1)
       vely(:,:,:) = diag(:,:,:,2)
       velz(:,:,:) = diag(:,:,:,3)
       pres(:,:,:) = diag(:,:,:,4)

    enddo ! RK loop

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       dens(i,j,k) = dens(i,j,k) + TIME_DTSEC_ATMOS_DYN * dens_t(i,j,k)
       momx(i,j,k) = momx(i,j,k) + TIME_DTSEC_ATMOS_DYN * vec_t(i,j,k,udir)
       momy(i,j,k) = momy(i,j,k) + TIME_DTSEC_ATMOS_DYN * vec_t(i,j,k,vdir)
       momz(i,j,k) = momz(i,j,k) + TIME_DTSEC_ATMOS_DYN * vec_t(i,j,k,wdir)
       lwpt(i,j,k) = lwpt(i,j,k) + TIME_DTSEC_ATMOS_DYN * scl_t(i,j,k,1)

       momx_t(i,j,k) = vec_t(i,j,k,udir)
       momy_t(i,j,k) = vec_t(i,j,k,vdir)
       momz_t(i,j,k) = vec_t(i,j,k,wdir)
       lwpt_t(i,j,k) = scl_t(i,j,k,1)
    enddo
    enddo
    enddo

    do v = 1,  QA
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       qtrc(i,j,k,v) = qtrc(i,j,k,v) + TIME_DTSEC_ATMOS_DYN * scl_t(i,j,k,v+1)

       qtrc_t(i,j,k,v) = scl_t(i,j,k,v+1)
    enddo
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_DYN

end module mod_atmos_dyn
