!-------------------------------------------------------------------------------
!> module LAND / Physics Bucket
!!
!! @par Description
!!          bucket-type land physics module
!!
!! @author Team SCALE
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_land_phy_bucket
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: LAND_PHY_driver_setup
  public :: LAND_PHY_driver_first
  public :: LAND_PHY_driver_final

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
  real(RP), allocatable :: GHFLX  (:,:)
  real(RP), allocatable :: PRECFLX(:,:)
  real(RP), allocatable :: QVFLX  (:,:)

  real(RP), allocatable :: dz(:)

  ! limiter
  real(RP), private, parameter :: BETA_MAX = 1.0_RP

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_PHY_driver_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_land_vars, only: &
       LAND_TYPE_PHY
    implicit none

    logical  :: dummy

    NAMELIST / PARAM_LAND_BUCKET / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[BUCKET]/Categ[LAND]'

    allocate( GHFLX  (IA,JA) )
    allocate( PRECFLX(IA,JA) )
    allocate( QVFLX  (IA,JA) )

    allocate( dz(LKS-1:LKE+1) )

    ! tentative
    dz(0) = 0.00_RP
    dz(1) = 0.02_RP
    dz(2) = 0.03_RP
    dz(3) = 0.05_RP
    dz(4) = 0.10_RP
    dz(5) = 0.30_RP
    dz(6) = 0.50_RP

    if ( LAND_TYPE_PHY /= 'BUCKET' ) then
       if( IO_L ) write(IO_FID_LOG,*) 'xxx LAND_TYPE_PHY is not BUCKET. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_LAND_BUCKET,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_LAND_BUCKET. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_LAND_BUCKET)

    return
  end subroutine LAND_PHY_driver_setup

  !-----------------------------------------------------------------------------
  !> Physical processes for land submodel
  subroutine LAND_PHY_driver_first
    use mod_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use mod_time, only: &
       dt => TIME_DTSEC_LAND
    use mod_land_vars, only: &
       TG,                 &
       STRG,               &
       ROFF,               &
       QVEF,               &
       I_STRGMAX,          &
       I_STRGCRT,          &
       I_TCS,              &
       I_HCS,              &
       P => LAND_PROPERTY, &
       LAND_vars_fillhalo
    use mod_cpl_vars, only: &
       CPL_getCPL2Lnd
    implicit none

    ! work
    integer :: k, i, j

    real(RP) :: ov(LKS:LKE)
    real(RP) :: iv(LKS:LKE)
    real(RP) :: ld(LKS:LKE)
    real(RP) :: md(LKS:LKE)
    real(RP) :: ud(LKS:LKE)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land step: Bucket'

    call CPL_getCPL2Lnd( GHFLX  (:,:), & ! [OUT]
                         PRECFLX(:,:), & ! [OUT]
                         QVFLX  (:,:)  ) ! [OUT]

    do j = JS, JE
    do i = IS, IE
      ld(:) = 0.0_RP
      md(:) = 0.0_RP
      ud(:) = 0.0_RP

      ! update water storage
      STRG(LKS,i,j) = STRG(LKS,i,j) + ( PRECFLX(i,j) + QVFLX(i,j) ) * dt

      do k = LKS, LKE
        if( STRG(k,i,j) > P(i,j,I_STRGMAX) ) then
          ! update run-off water
          ROFF(i,j) = ROFF(i,j) + STRG(k,i,j) - P(i,j,I_STRGMAX)

          ! modify STRG to the upper limit
          STRG(k,i,j) = P(i,j,I_STRGMAX)
        end if
      end do

      ! update moisture efficiency
      QVEF(i,j) = min( STRG(LKS,i,j)/P(i,j,I_STRGCRT), BETA_MAX )

      ! prepare to solve tridiagonal matrix
      do k = LKS+1, LKE
        ld(k) = -2.0_RP * dt * P(i,j,I_TCS) / ( dz(k) * ( dz(k) + dz(k-1) ) ) &
              / ( ( 1.0_RP - P(i,j,I_STRGMAX) * 1.E-3_RP ) * P(i,j,I_HCS) + STRG(k,i,j) * 1.E-3_RP * DWATR * CL )
      end do
      do k = LKS, LKE
        ud(k) = -2.0_RP * dt * P(i,j,I_TCS) / ( dz(k) * ( dz(k) + dz(k+1) ) ) &
              / ( ( 1.0_RP - P(i,j,I_STRGMAX) * 1.E-3_RP ) * P(i,j,I_HCS) + STRG(k,i,j) * 1.E-3_RP * DWATR * CL )
        md(k) = 1.0_RP - ld(k) - ud(k)
      end do

      iv(:)   = TG(LKS:LKE,i,j)
      iv(LKS) = TG(LKS,i,j) - dt * 2.0_RP * GHFLX(i,j) / dz(k) &
              / ( ( 1.0_RP - P(i,j,I_STRGMAX) * 1.E-3_RP ) * P(i,j,I_HCS) + STRG(LKS,i,j) * 1.E-3_RP * DWATR * CL )
      iv(LKE) = TG(LKE,i,j) - ud(LKE) * TG(LKE+1,i,j)

      call solve_tridiagonal_matrix( &
        ov(:),                     & ! (out)
        iv(:), ld(:), md(:), ud(:) ) ! (in)

      ! update ground temperature
      TG(LKS:LKE,i,j) = ov(:)

    end do
    end do

    call LAND_vars_fillhalo

    return
  end subroutine LAND_PHY_driver_first

  subroutine LAND_PHY_driver_final
    use mod_land_vars, only: &
       TG,                 &
       QVEF,               &
       I_EMIT,             &
       I_ALBG,             &
       I_TCS,              &
       I_DZG,              &
       I_Z0M,              &
       I_Z0H,              &
       I_Z0E,              &
       P => LAND_PROPERTY
    use mod_cpl_vars, only: &
       CPL_putLnd
    implicit none
    !---------------------------------------------------------------------------

    call CPL_putLnd( TG  (LKS,:,:),    & ! [IN]
                     QVEF(:,:),        & ! [IN]
                     P   (:,:,I_EMIT), & ! [IN]
                     P   (:,:,I_ALBG), & ! [IN]
                     P   (:,:,I_TCS),  & ! [IN]
                     P   (:,:,I_DZG),  & ! [IN]
                     P   (:,:,I_Z0M),  & ! [IN]
                     P   (:,:,I_Z0H),  & ! [IN]
                     P   (:,:,I_Z0E)   ) ! [IN]

    return
  end subroutine LAND_PHY_driver_final

  ! solve tridiagonal matrix with Thomas's algorithm
  subroutine solve_tridiagonal_matrix( &
      ov,            & ! (out)
      iv, ld, md, ud ) ! (in)
    implicit none

    ! argument
    real(RP), intent(out) :: ov(:) ! output vector

    real(RP), intent(in) :: iv(:) ! input vector
    real(RP), intent(in) :: ld(:) ! lower diagonal
    real(RP), intent(in) :: md(:) ! middle diagonal
    real(RP), intent(in) :: ud(:) ! upper diagonal

    ! work
    integer :: n, k

    real(RP) :: a(size(ld))
    real(RP) :: b(size(md))
    real(RP) :: c(size(ud))
    real(RP) :: d(size(iv))
    real(RP) :: tmp
    !---------------------------------------------------------------------------

    ! maximum array size
    n = size( ov(:) )

    ! initialize
    a(:) = ld(:)
    b(:) = md(:)
    c(:) = ud(:)
    d(:) = iv(:)

    c(1) = c(1) / b(1)
    d(1) = d(1) / b(1)

    ! foward reduction
    do k = 2, n
      tmp  = b(k) - a(k) * c(k-1)
      c(k) = c(k) / tmp
      d(k) = ( d(k) - a(k) * d(k-1) ) / tmp
    end do

    ! backward substitution
    ov(n) = d(n)
    do k = n-1, 1, -1
      ov(k) = d(k) - c(k) * ov(k+1)
    end do

    return
  end subroutine solve_tridiagonal_matrix

end module mod_land_phy_bucket
