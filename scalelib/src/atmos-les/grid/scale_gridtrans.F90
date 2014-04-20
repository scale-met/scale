!-------------------------------------------------------------------------------
!> module GRIDTRANS
!!
!! @par Description
!!          Grid transfer module
!!          Map projection and Terrain-following metrics
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-10-24 (H.Yashiro)  [new] reconstruct from scale_topography
!!
!<
!-------------------------------------------------------------------------------
module scale_gridtrans
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GTRANS_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, allocatable :: GTRANS_MAPF (:,:,:)   !< map factor

  real(RP), public, allocatable :: GTRANS_GSQRT(:,:,:,:) !< transformation metrics from Z to Xi, {G}^1/2
  real(RP), public, allocatable :: GTRANS_J13G (:,:,:,:) !< (1,3) element of Jacobian matrix * {G}^1/2
  real(RP), public, allocatable :: GTRANS_J23G (:,:,:,:) !< (2,3) element of Jacobian matrix * {G}^1/2
  real(RP), public              :: GTRANS_J33G           !< (3,3) element of Jacobian matrix * {G}^1/2

  integer,  public :: I_XYZ = 1 ! at (x,y,z)
  integer,  public :: I_XYW = 2 ! at (x,y,w)
  integer,  public :: I_UYW = 3 ! at (u,y,w)
  integer,  public :: I_XVW = 4 ! at (x,v,w)
  integer,  public :: I_UYZ = 5 ! at (u,y,z)
  integer,  public :: I_XVZ = 6 ! at (x,v,z)
  integer,  public :: I_UVZ = 7 ! at (u,v,z)

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: GTRANS_mapfactor
  private :: GTRANS_terrainfollowing

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine GTRANS_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[GRIDTRANS]/Categ[GRID]'

    allocate( GTRANS_MAPF (IA,JA,4) )

    allocate( GTRANS_GSQRT(KA,IA,JA,7) )
    allocate( GTRANS_J13G (KA,IA,JA,7) )
    allocate( GTRANS_J23G (KA,IA,JA,7) )

    ! calc metrics for orthogonal curvelinear coordinate
    call GTRANS_mapfactor

    ! calc metrics for terrain-following coordinate
    call GTRANS_terrainfollowing

    return
  end subroutine GTRANS_setup

  !-----------------------------------------------------------------------------
  !> Calculate map factor
  subroutine GTRANS_mapfactor
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine GTRANS_mapfactor

  !-----------------------------------------------------------------------------
  !> Calculate G^1/2 & Jacobian
  subroutine GTRANS_terrainfollowing
    use scale_grid, only: &
       GRID_RCDZ, &
       GRID_RCDX, &
       GRID_RCDY, &
       GRID_RFDZ, &
       GRID_RFDX, &
       GRID_RFDY
    use scale_grid_real, only: &
       REAL_CZ, &
       REAL_FZ
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: REAL_CZ_U (  KA,IA,JA) !< Z coordinate [m] at (u,y,z)
    real(RP) :: REAL_CZ_V (  KA,IA,JA) !< Z coordinate [m] at (x,v,z)
    real(RP) :: REAL_CZ_UV(  KA,IA,JA) !< Z coordinate [m] at (u,y,z)
    real(RP) :: REAL_FZ_U (0:KA,IA,JA) !< Z coordinate [m] at (u,y,w)
    real(RP) :: REAL_FZ_V (0:KA,IA,JA) !< Z coordinate [m] at (x,v,w)
    real(RP) :: REAL_FZ_UV(0:KA,IA,JA) !< Z coordinate [m] at (u,v,w)

    integer :: k, i, j
    !---------------------------------------------------------------------------

    ! calc Z-coordinate height at staggered position
    do j = 1, JA
    do i = 1, IA-1
    do k = 1, KA
       REAL_CZ_U(k,i,j) = 0.5D0 * ( REAL_CZ(k,i+1,j) + REAL_CZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA
    do i = 1, IA-1
    do k = 0, KA
       REAL_FZ_U(k,i,j) = 0.5D0 * ( REAL_FZ(k,i+1,j) + REAL_FZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = 1, KA
       REAL_CZ_V(k,i,j) = 0.5D0 * ( REAL_CZ(k,i,j+1) + REAL_CZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA
    do k = 0, KA
       REAL_FZ_V(k,i,j) = 0.5D0 * ( REAL_FZ(k,i,j+1) + REAL_FZ(k,i,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
    do k = 1, KA
       REAL_CZ_UV(k,i,j) = 0.25D0 * ( REAL_CZ(k,i+1,j+1) + REAL_CZ(k,i+1,j) &
                                    + REAL_CZ(k,i  ,j+1) + REAL_CZ(k,i  ,j) )
    enddo
    enddo
    enddo

    do j = 1, JA-1
    do i = 1, IA-1
    do k = 0, KA
       REAL_FZ_UV(k,i,j) = 0.25D0 * ( REAL_FZ(k,i+1,j+1) + REAL_FZ(k,i+1,j) &
                                    + REAL_FZ(k,i  ,j+1) + REAL_FZ(k,i  ,j) )
    enddo
    enddo
    enddo

    ! G^1/2
    do j = JS, JE
    do i = IS, IE
       ! at (x,y,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_XYZ) = ( REAL_FZ(k,i,j) - REAL_FZ(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (x,y,w)
       do k = 1, KA-1
          GTRANS_GSQRT(k,i,j,I_XYW) = ( REAL_CZ(k+1,i,j) - REAL_CZ(k,i,j) ) * GRID_RFDZ(k)
       enddo
       GTRANS_GSQRT(KA,i,j,I_XYW) = GTRANS_GSQRT(KA-1,i,j,I_XYW)

       ! at (u,y,w)
       do k = 1, KA-1
          GTRANS_GSQRT(k,i,j,I_UYW) = ( REAL_CZ_U(k+1,i,j) - REAL_CZ_U(k,i,j) ) * GRID_RFDZ(k)
       enddo
       GTRANS_GSQRT(KA,i,j,I_UYW) = GTRANS_GSQRT(KA-1,i,j,I_UYW)

       ! at (x,v,w)
       do k = 1, KA-1
          GTRANS_GSQRT(k,i,j,I_XVW) = ( REAL_CZ_V(k+1,i,j) - REAL_CZ_V(k,i,j) ) * GRID_RFDZ(k)
       enddo
       GTRANS_GSQRT(KA,i,j,I_XVW) = GTRANS_GSQRT(KA-1,i,j,I_XVW)

       ! at (u,y,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_UYZ) = ( REAL_FZ_U(k,i,j) - REAL_FZ_U(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (x,v,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_XVZ) = ( REAL_FZ_V(k,i,j) - REAL_FZ_V(k-1,i,j) ) * GRID_RCDZ(k)
       enddo

       ! at (u,v,z)
       do k = 1, KA
          GTRANS_GSQRT(k,i,j,I_UVZ) = ( REAL_FZ_UV(k,i,j) - REAL_FZ_UV(k-1,i,j) ) * GRID_RCDZ(k)
       enddo
    enddo
    enddo

    call COMM_vars8( GTRANS_GSQRT(:,:,:,1), 1 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,2), 2 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,3), 3 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,4), 4 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,5), 5 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,6), 6 )
    call COMM_vars8( GTRANS_GSQRT(:,:,:,7), 7 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,1), 1 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,2), 2 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,3), 3 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,4), 4 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,5), 5 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,6), 6 )
    call COMM_wait ( GTRANS_GSQRT(:,:,:,7), 7 )

    ! Jacobian * G^1/2
    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       GTRANS_J13G(k,i,j,I_XYZ) = -( REAL_CZ_U (k,i  ,j) - REAL_CZ_U (k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_XYW) = -( REAL_FZ_U (k,i  ,j) - REAL_FZ_U (k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_UYW) = -( REAL_FZ   (k,i+1,j) - REAL_FZ   (k,i  ,j) ) * GRID_RFDX(i)
       GTRANS_J13G(k,i,j,I_XVW) = -( REAL_FZ_UV(k,i  ,j) - REAL_FZ_UV(k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_UYZ) = -( REAL_CZ   (k,i  ,j) - REAL_CZ   (k,i-1,j) ) * GRID_RFDX(i)
       GTRANS_J13G(k,i,j,I_XVZ) = -( REAL_CZ_UV(k,i  ,j) - REAL_CZ_UV(k,i-1,j) ) * GRID_RCDX(i)
       GTRANS_J13G(k,i,j,I_UVZ) = -( REAL_CZ_V (k,i  ,j) - REAL_CZ_V (k,i-1,j) ) * GRID_RFDX(i)
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = 1,  KA
       GTRANS_J23G(k,i,j,I_XYZ) = -( REAL_CZ_V (k,i,j  ) - REAL_CZ_V (k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_XYW) = -( REAL_FZ_V (k,i,j  ) - REAL_FZ_V (k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_UYW) = -( REAL_FZ   (k,i,j+1) - REAL_FZ   (k,i,j  ) ) * GRID_RFDY(j)
       GTRANS_J23G(k,i,j,I_XVW) = -( REAL_FZ_UV(k,i,j  ) - REAL_FZ_UV(k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_UYZ) = -( REAL_CZ_UV(k,i,j  ) - REAL_CZ_UV(k,i,j-1) ) * GRID_RCDY(j)
       GTRANS_J23G(k,i,j,I_XVZ) = -( REAL_CZ   (k,i,j  ) - REAL_CZ   (k,i,j-1) ) * GRID_RFDY(j)
       GTRANS_J23G(k,i,j,I_UVZ) = -( REAL_CZ_U (k,i,j  ) - REAL_CZ_U (k,i,j-1) ) * GRID_RFDY(j)
    enddo
    enddo
    enddo

    GTRANS_J33G = 1.0_RP ! - 1 / G^1/2 * G^1/2

    ! fill IHALO & JHALO
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XYZ),  1 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XYW),  2 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_UYW),  3 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XVW),  4 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_UYZ),  5 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_XVZ),  6 )
    call COMM_vars8( GTRANS_J13G(:,:,:,I_UVZ),  7 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XYZ),  8 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XYW),  9 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_UYW), 10 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XVW), 11 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_UYZ), 12 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_XVZ), 13 )
    call COMM_vars8( GTRANS_J23G(:,:,:,I_UVZ), 14 )

    call COMM_wait ( GTRANS_J13G(:,:,:,I_XYZ),  1 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XYW),  2 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_UYW),  3 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XVW),  4 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_UYZ),  5 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_XVZ),  6 )
    call COMM_wait ( GTRANS_J13G(:,:,:,I_UVZ),  7 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XYZ),  8 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XYW),  9 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_UYW), 10 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XVW), 11 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_UYZ), 12 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_XVZ), 13 )
    call COMM_wait ( GTRANS_J23G(:,:,:,I_UVZ), 14 )

    return
  end subroutine GTRANS_terrainfollowing

end module scale_gridtrans
