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
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_land_grid_index
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
  logical, private, save :: LAND_LKE_STRG_UPDATE = .false. ! is STRG updated in the lowest level?
  logical, private, save :: LAND_LKE_TG_UPDATE   = .false. ! is TG updated in the lowest level?

  real(RP), private, save, allocatable :: GHFLX  (:,:)
  real(RP), private, save, allocatable :: PRECFLX(:,:)
  real(RP), private, save, allocatable :: QVFLX  (:,:)

  ! limiter
  real(RP), private, parameter :: BETA_MAX = 1.0_RP

contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine LAND_PHY_driver_setup
    use scale_process, only: &
       PRC_MPIstop
    use mod_land_vars, only: &
       LAND_TYPE_PHY
    implicit none

    logical  :: dummy

    NAMELIST / PARAM_LAND_BUCKET / &
       LAND_LKE_STRG_UPDATE, &
       LAND_LKE_TG_UPDATE

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[BUCKET]/Categ[LAND]'

    allocate( GHFLX  (IA,JA) )
    allocate( PRECFLX(IA,JA) )
    allocate( QVFLX  (IA,JA) )

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
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       CL    => CONST_CL
    use scale_time, only: &
       dt => TIME_DTSEC_LAND
    use scale_land_grid, only: &
       dz => GRID_LCDZ
    use mod_land_vars, only: &
       TG,                 &
       STRG,               &
       ROFF,               &
       QVEF,               &
       I_STRGMAX,          &
       I_STRGCRT,          &
       I_TCS,              &
       I_HCS,              &
       I_DFW,              &
       P => LAND_PROPERTY, &
       LAND_vars_fillhalo
    use mod_cpl_vars, only: &
       CPL_getCPL2Lnd
    implicit none

    ! work
    integer :: k, i, j

    real(RP) :: ov(LKS:LKE-1)
    real(RP) :: iv(LKS:LKE-1)
    real(RP) :: ld(LKS:LKE-1)
    real(RP) :: md(LKS:LKE-1)
    real(RP) :: ud(LKS:LKE-1)
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Land step: Bucket'

    call CPL_getCPL2Lnd( GHFLX  (:,:), & ! [OUT]
                         PRECFLX(:,:), & ! [OUT]
                         QVFLX  (:,:)  ) ! [OUT]

    do j = JS, JE
    do i = IS, IE
      ! prepare to solve tridiagonal matrix for STRG
      ld(:) = 0.0_RP
      md(:) = 0.0_RP
      ud(:) = 0.0_RP

      do k = LKS+1, LKE-1
        ld(k) = -2.0_RP * dt * P(i,j,I_DFW) / ( dz(k) * ( dz(k) + dz(k-1) ) )
      end do
      do k = LKS, LKE-1
        ud(k) = -2.0_RP * dt * P(i,j,I_DFW) / ( dz(k) * ( dz(k) + dz(k+1) ) )
        md(k) = 1.0_RP - ld(k) - ud(k)
      end do

      iv(:)     = STRG(LKS:LKE-1,i,j)
      iv(LKS)   = STRG(LKS,i,j) + ( PRECFLX(i,j) + QVFLX(i,j) ) * dt
      iv(LKE-1) = STRG(LKE-1,i,j) - ud(LKE-1) * STRG(LKE,i,j)

      call solve_tridiagonal_matrix( &
        ov(:),                     & ! (out)
        iv(:), ld(:), md(:), ud(:) ) ! (in)

      ! update water storage
      STRG(LKS:LKE-1,i,j) = ov(:)

      ! update lowest STRG
      if( LAND_LKE_STRG_UPDATE ) then
        STRG(LKE,i,j) = STRG(LKE-1,i,j)
      end if

      ! modify STRG to the upper limit
      do k = LKS, LKE
        if( STRG(k,i,j) > P(i,j,I_STRGMAX) ) then
          ! update vertically integral run-off water
          ROFF(i,j)   = ROFF(i,j) + STRG(k,i,j) - P(i,j,I_STRGMAX)
          STRG(k,i,j) = P(i,j,I_STRGMAX)
        end if
      end do

      ! update moisture efficiency
      QVEF(i,j) = min( STRG(LKS,i,j)/P(i,j,I_STRGCRT), BETA_MAX )

      ! prepare to solve tridiagonal matrix for TG
      ld(:) = 0.0_RP
      md(:) = 0.0_RP
      ud(:) = 0.0_RP

      do k = LKS+1, LKE-1
        ld(k) = -2.0_RP * dt * P(i,j,I_TCS) / ( dz(k) * ( dz(k) + dz(k-1) ) ) &
              / ( ( 1.0_RP - P(i,j,I_STRGMAX) / dz(k) * 1.0E-3_RP ) * P(i,j,I_HCS) & ! heat capacity for soil
                + STRG(k,i,j) / dz(k) * 1.0E-3_RP * DWATR * CL ) ! heat capacity for water
      end do
      do k = LKS, LKE-1
        ud(k) = -2.0_RP * dt * P(i,j,I_TCS) / ( dz(k) * ( dz(k) + dz(k+1) ) ) &
              / ( ( 1.0_RP - P(i,j,I_STRGMAX) / dz(k) * 1.0E-3_RP ) * P(i,j,I_HCS) & ! heat capacity for soil
                + STRG(k,i,j) / dz(k) * 1.0E-3_RP * DWATR * CL ) ! heat capacity for water
        md(k) = 1.0_RP - ld(k) - ud(k)
      end do

      iv(:)     = TG(LKS:LKE-1,i,j)
      iv(LKS)   = TG(LKS,i,j) - dt * 2.0_RP * GHFLX(i,j) / dz(k) &
                / ( ( 1.0_RP - P(i,j,I_STRGMAX) / dz(LKS) * 1.0E-3_RP ) * P(i,j,I_HCS) & ! heat capacity for soil
                  + STRG(LKS,i,j) / dz(LKS) * 1.0E-3_RP * DWATR * CL ) ! heat capacity for water
      iv(LKE-1) = TG(LKE-1,i,j) - ud(LKE-1) * TG(LKE,i,j)

      call solve_tridiagonal_matrix( &
        ov(:),                     & ! (out)
        iv(:), ld(:), md(:), ud(:) ) ! (in)

      ! update ground temperature
      TG(LKS:LKE-1,i,j) = ov(:)

      ! update lowest TG
      if( LAND_LKE_TG_UPDATE ) then
        TG(LKE,i,j) = TG(LKE-1,i,j)
      end if

    end do
    end do

    call LAND_vars_fillhalo

    return
  end subroutine LAND_PHY_driver_first

  subroutine LAND_PHY_driver_final
    use scale_land_grid, only: &
       dz => GRID_LCDZ
    use mod_land_vars, only: &
       TG,                 &
       QVEF,               &
       I_TCS,              &
       I_Z0M,              &
       I_Z0H,              &
       I_Z0E,              &
       P => LAND_PROPERTY
    use mod_cpl_vars, only: &
       CPL_putLnd
    implicit none

    real(RP) :: dz_h(IA,JA)
    !---------------------------------------------------------------------------

    dz_h(:,:) = dz(LKS)

    call CPL_putLnd( TG  (LKS,:,:),    & ! [IN]
                     QVEF(:,:),        & ! [IN]
                     P   (:,:,I_TCS),  & ! [IN]
                     dz_h(:,:),        & ! [IN]
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
