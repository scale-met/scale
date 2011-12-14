!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!          Smagolinsky-type
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-xx-xx (S.Iga) [new]
!! @li      2011-12-11 (H.Yashiro) [mod] integrate to SCALE3
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb
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
  public :: ATMOS_PHY_TB
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

  real(8), private, parameter :: Cs        = 0.18D0 ! (Sullivan et al.1994, Nakanishi and Niino)
  real(8), private, parameter :: GAMMA     = 1.D0   ! assume dx=dy=dz
  real(8), private, parameter :: Ustar_min = 1.D-6  ! minimum limit of U*

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Smagorinsky-type turblence
  !> comment:
  !>  1, Pr is given linearly (iga)
  !>  2, dx,dy,dz is assumed to be constant.
  !>  3, gamma is assumed to be 1 now (i.e. dx=dy=dz).
  !>  4, heat flux is not accurate yet. (i.e. energy is not conserved, see *1)
  !>  5, stratification effect is considered.
  !>  6, This routine does not deal with surface flux.
  !>  7, limiter is not implemented.
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB( dens,   pott,   qtrc,                  &
                           velx,   vely,   velz,                  &
                           FLXij_sfc, FLXt_sfc, FLXqv_sfc,        &
                           momx_t, momy_t, momz_t, pott_t, qtrc_t )
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rair   => CONST_Rair,   &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd,   &
       CONST_UNDEF8
    use mod_grid, only : &
       IA  => GRID_IA,  &
       JA  => GRID_JA,  &
       KA  => GRID_KA,  &
       IS  => GRID_IS,  &
       IE  => GRID_IE,  &
       JS  => GRID_JS,  &
       JE  => GRID_JE,  &
       KS  => GRID_KS,  &
       KE  => GRID_KE,  &
       WS  => GRID_WS,  &
       WE  => GRID_WE,  &
       DX  => GRID_DX,  &
       RDX => GRID_RDX, &
       RDY => GRID_RDY, &
       RDZ => GRID_RDZ
    use mod_atmos_vars, only: &
       QA => A_QA
    implicit none

    ! prognostic value
    real(8), intent(in)  :: dens(IA,JA,KA)      ! density [kg/m3]
    real(8), intent(in)  :: velx(IA,JA,KA)      ! velocity(x) [m/s]
    real(8), intent(in)  :: vely(IA,JA,KA)      ! velocity(y) [m/s]
    real(8), intent(in)  :: velz(IA,JA,KA)      ! velocity(z) [m/s]
    real(8), intent(in)  :: pott(IA,JA,KA)      ! potential temperature [K]
    real(8), intent(in)  :: qtrc(IA,JA,KA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    ! surface flux
    real(8), intent(in)  :: FLXij_sfc(IA,JA,3)  ! => FLXij(1:IA,1:JA,WS,1:3,3)
    real(8), intent(in)  :: FLXt_sfc (IA,JA)    ! => FLXt (1:IA,1:JA,WS)
    real(8), intent(in)  :: FLXqv_sfc(IA,JA)    ! => FLXq (1:IA,1:JA,WS,1)

    ! prognostic tendency
    real(8), intent(out) :: momx_t(IA,JA,KA)
    real(8), intent(out) :: momy_t(IA,JA,KA)
    real(8), intent(out) :: momz_t(IA,JA,KA)
    real(8), intent(out) :: pott_t(IA,JA,KA)
    real(8), intent(out) :: qtrc_t(IA,JA,KA,QA)

    real(8), allocatable :: Uij(:,:,:,:,:) ! velocity gradient matrix
    real(8), allocatable :: Sij(:,:,:,:,:) ! deformation rate tensor

    real(8), allocatable :: nu   (:,:,:)   ! eddy-viscosity coefficient
    real(8), allocatable :: nurho(:,:,:)   ! nu * dens

    real(8), allocatable :: FLXij(:,:,:,:,:)
    real(8), allocatable :: FLXt (:,:,:,:)
    real(8), allocatable :: FLXq (:,:,:,:,:)

    real(8)              :: Ri         ! Richardson number
    real(8)              :: Pr         ! Prantle number
    real(8), allocatable :: RPr(:,:,:) ! 1 / Pr

    integer :: i, j, k, ii, jj
    !---------------------------------------------------------------------------

    allocate( Uij  (IA,JA,KA,3,3)  ); Uij  (:,:,:,:,:) = CONST_UNDEF8
    allocate( Sij  (IA,JA,KA,3,3)  ); Sij  (:,:,:,:,:) = CONST_UNDEF8
    allocate( nu   (IA,JA,KA)      ); nu   (:,:,:)     = CONST_UNDEF8
    allocate( nurho(IA,JA,KA)      ); nurho(:,:,:)     = CONST_UNDEF8
    allocate( RPr  (IA,JA,KA)      ); RPr  (:,:,:)     = CONST_UNDEF8
    allocate( FLXij(IA,JA,KA,3,3)  ); FLXij(:,:,:,:,:) = CONST_UNDEF8
    allocate( FLXt (IA,JA,KA,3)    ); FLXt (:,:,:,:)   = CONST_UNDEF8
    allocate( FLXq (IA,JA,KA,3,QA) ); FLXq (:,:,:,:,:) = CONST_UNDEF8

    ! ii == jj : at the center of control volume (cube)
    do k = KS,   KE
    do j = JS-2, JE+2
    do i = IS-1, IE+2
       Uij(i,j,k,1,1) = ( velx(i,j,k)-velx(i-1,j,k) ) * RDX ! du/dx, (x, y, layer)
    enddo
    enddo
    enddo
    do k = KS,   KE
    do j = JS-1, JE+2
    do i = IS-2, IE+2
       Uij(i,j,k,2,2) = ( vely(i,j,k)-vely(i,j-1,k) ) * RDY ! dv/dy, (x, y, layer)
    enddo
    enddo
    enddo
    do k = KS,   KE
    do j = JS-2, JE+2
    do i = IS-2, IE+2
       Uij(i,j,k,3,3) = ( velz(i,j,k)-velz(i,j,k-1) ) * RDZ ! dw/dz, (x, y, layer)
    enddo
    enddo
    enddo
    ! ii /= jj : on the edge   of control volume (cube)
    do k = KS,   KE
    do j = JS-2, JE+2
    do i = IS-2, IE+1
       Uij(i,j,k,1,2) = ( vely(i+1,j,k)-vely(i,j,k) ) * RDX ! dv/dx, (u, v, layer)
    enddo
    enddo
    enddo
    do k = WS+1, WE-1
    do j = JS-2, JE+2
    do i = IS-2, IE+1
       Uij(i,j,k,1,3) = ( velz(i+1,j,k)-velz(i,j,k) ) * RDX ! dw/dx, (u, y, interface)
    enddo
    enddo
    enddo
    Uij(:,:,WS,1,3) = 0.D0
    Uij(:,:,WE,1,3) = 0.D0
    do k = KS,   KE
    do j = JS-2, JE+1
    do i = IS-2, IE+2
       Uij(i,j,k,2,1) = ( velx(i,j+1,k)-velx(i,j,k) ) * RDY ! du/dy, (u, v, layer)
    enddo
    enddo
    enddo
    do k = WS+1, WE-1
    do j = JS-2, JE+1
    do i = IS-2, IE+2
       Uij(i,j,k,2,3) = ( velz(i,j+1,k)-velz(i,j,k) ) * RDY ! dw/dy, (x, v, interface)
    enddo
    enddo
    enddo
    Uij(:,:,WS,2,3) = 0.D0
    Uij(:,:,WE,2,3) = 0.D0
    do k = WS+1, WE-1
    do j = JS-2, JE+2
    do i = IS-2, IE+2
       Uij(i,j,k,3,1) = ( velx(i,j,k+1)-velx(i,j,k) ) * RDZ ! du/dz, (u, y, interface)
       Uij(i,j,k,3,2) = ( vely(i,j,k+1)-vely(i,j,k) ) * RDZ ! dv/dz, (x, v, interface)
    enddo
    enddo
    enddo
    Uij(:,:,WS,3,1) = 0.D0
    Uij(:,:,WS,3,1) = 0.D0
    Uij(:,:,WE,3,2) = 0.D0
    Uij(:,:,WE,3,2) = 0.D0

    do k = KS,   KE
    do j = JS-2, JE+2
    do i = IS-2, IE+2
       Sij(i,j,k,1,1) = Uij(i,j,k,1,1) * dens(i,j,k)
       Sij(i,j,k,2,2) = Uij(i,j,k,2,2) * dens(i,j,k)
       Sij(i,j,k,3,3) = Uij(i,j,k,3,3) * dens(i,j,k)
    enddo
    enddo
    enddo
    do k = KS,   KE
    do j = JS-2, JE+1
    do i = IS-2, IE+1
       Sij(i,j,k,1,2) = 0.125D0 * ( Uij(i,j,k,1,2)+Uij(i,j,k,2,1) )               &
                      * ( dens(i,j,k)+dens(i+1,j,k)+dens(i,j+1,k)+dens(i+1,j+1,k) )
       Sij(i,j,k,2,1) = Sij(i,j,k,1,2)
    enddo
    enddo
    enddo
    do k = WS+1, WE-1
    do j = JS-2, JE+2
    do i = IS-2, IE+1
       Sij(i,j,k,1,3) = 0.125D0 * ( Uij(i,j,k,1,3)+Uij(i,j,k,3,1) )               &
                      * ( dens(i,j,k)+dens(i+1,j,k)+dens(i,j,k+1)+dens(i+1,j,k+1) )
       Sij(i,j,k,3,1) = Sij(i,j,k,1,3)
    enddo
    enddo
    enddo
    do k = WS+1, WE-1
    do j = JS-2, JE+1
    do i = IS-2, IE+2
       Sij(i,j,k,2,3) = 0.125D0 * ( Uij(i,j,k,2,3)+Uij(i,j,k,3,2) )               &
                      * ( dens(i,j,k)+dens(i,j+1,k)+dens(i,j,k+1)+dens(i,j+1,k+1) )
       Sij(i,j,k,3,2) = Sij(i,j,k,2,3)
    enddo
    enddo
    enddo

    ! calculate stratification effect (n2 and Pr are defined at w level)
!    do k = KS, KE
!    do j = JS, JE
!    do i = IS, IE
!       a = 1.D0 / ( 1.D0 + dqsl(i,j,k) * Lv / Cp )
!
!       c = (1+0.61*qw(i,j,k)-1.61*ql(i,j,k))* Lv/Cp /T_ov_pott(i,j,k) - 1.61 * pott(i,j,k)
!
!       betatheta= 1 + 0.61 * qw(i,j,k) - 1.61*ql(i,j,k) - a*c*dqsl(i,j,k)* T_ov_pott(i,j,k)
!
!       betaq= 0.61 * pott(i,j,k) + a * c
!
!       n2(i,j,k) = grav / (sum(pottl(i,j,k:kz+1))/2) * &
!                     (sum(betatheta * ( pottl(i,j,k+1)-pottl(i,j,k) )/dz + &
!                     betaq * ( qw(i,j,k+1)-qw(i,j,k)  )/dz )/4)
!    enddo
!    enddo
!    enddo

    ! calculate  Pr (Pr is defined at w level)
    do k = WS+1, WE-1
    do j = JS-1, JE
    do i = IS-1, IE
       ! Ri = g / rho * drho/dz / dU/dz
       Ri = GRAV * ( dens(i,j,k+1)-dens(i,j,k) ) * RDZ * 2.D0 / ( dens(i,j,k+1)+dens(i,j,k) )                   &
          / max( ( ( velx(i,j,k+1)+velx(i-1,j,k+1)-velx(i,j,k)-velx(i-1,j,k) ) * 0.5D0 * RDZ )**2.D0            &
               + ( ( vely(i,j,k+1)+vely(i,j-1,k+1)-vely(i,j,k)-vely(i,j-1,k) ) * 0.5D0 * RDZ )**2.D0, Ustar_min )

       Ri         = min( 0.25D0, max( 0.D0, Ri ) )
       Pr         = ( 1.D0-4.D0*Ri ) / 3.D0 + 4.D0*Ri
       RPr(i,j,k) = 1.D0 / Pr
    enddo
    enddo
    enddo
    RPr(:,:,WS) = 0.D0
    RPr(:,:,WE) = 0.D0

    do k = WS+1, WE-1
    do j = JS-2, JE+1
    do i = IS-2, IE+1
       nurho(i,j,k) = ( Cs * DX ) ** 2.D0 * ( sqrt( 2.D0 * GAMMA                            &
                    * ( ( Sij(i  ,j  ,k  ,1,1) + Sij(i+1,j+1,k+1,1,1)                       &
                        + Sij(i+1,j  ,k  ,1,1) + Sij(i+1,j+1,k  ,1,1)                       &
                        + Sij(i  ,j+1,k  ,1,1) + Sij(i  ,j+1,k+1,1,1)                       &
                        + Sij(i  ,j  ,k+1,1,1) + Sij(i+1,j  ,k+1,1,1) ) * 0.125D0 ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,1,2) + Sij(i  ,j  ,k+1,1,2) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,1,3) + Sij(i  ,j+1,k  ,1,3) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,2,1) + Sij(i  ,j  ,k+1,2,1) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,2,2) + Sij(i+1,j+1,k+1,2,2)                       &
                        + Sij(i+1,j  ,k  ,2,2) + Sij(i+1,j+1,k  ,2,2)                       &
                        + Sij(i  ,j+1,k  ,2,2) + Sij(i  ,j+1,k+1,2,2)                       &
                        + Sij(i  ,j  ,k+1,2,2) + Sij(i+1,j  ,k+1,2,2) ) * 0.125D0 ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,2,3) + Sij(i+1,j  ,k  ,2,3) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,3,1) + Sij(i  ,j+1,k  ,3,1) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,3,2) + Sij(i+1,j  ,k  ,3,2) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(i  ,j  ,k  ,3,3) + Sij(i+1,j+1,k+1,3,3)                       &
                        + Sij(i+1,j  ,k  ,3,3) + Sij(i+1,j+1,k  ,3,3)                       &
                        + Sij(i  ,j+1,k  ,3,3) + Sij(i  ,j+1,k+1,3,3)                       &
                        + Sij(i  ,j  ,k+1,3,3) + Sij(i+1,j  ,k+1,3,3) ) * 0.125D0 ) ** 2.D0 &
                    ) )
    enddo
    enddo
    enddo
    nurho(:,:,WS) = 0.D0
    nurho(:,:,WE) = 0.D0

    do k = KS,   KE
    do j = JS-1, JE
    do i = IS-1, IE
       FLXij(i,j,k,1,2) = GAMMA * Sij(i,j,k,1,2) * ( nurho(i  ,j  ,k  ) + nurho(i  ,j  ,k-1) )
       FLXij(i,j,k,1,3) = GAMMA * Sij(i,j,k,1,3) * ( nurho(i  ,j  ,k  ) + nurho(i  ,j-1,k  ) )
       FLXij(i,j,k,2,1) = GAMMA * Sij(i,j,k,2,1) * ( nurho(i  ,j  ,k  ) + nurho(i  ,j  ,k-1) )
       FLXij(i,j,k,2,3) = GAMMA * Sij(i,j,k,2,3) * ( nurho(i  ,j  ,k  ) + nurho(i-1,j  ,k  ) )
       FLXij(i,j,k,3,1) = GAMMA * Sij(i,j,k,3,1) * ( nurho(i  ,j  ,k  ) + nurho(i  ,j-1,k  ) )
       FLXij(i,j,k,3,2) = GAMMA * Sij(i,j,k,3,2) * ( nurho(i  ,j  ,k  ) + nurho(i-1,j  ,k  ) )
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE+1
    do i = IS, IE+1
       FLXij(i,j,k,1,1) = GAMMA * Sij(i,j,k,1,1) * ( nurho(i  ,j  ,k  ) + nurho(i-1,j-1,k-1) &
                                                   + nurho(i-1,j  ,k  ) + nurho(i-1,j-1,k  ) &
                                                   + nurho(i  ,j-1,k  ) + nurho(i  ,j-1,k-1) &
                                                   + nurho(i  ,j  ,k-1) + nurho(i-1,j  ,k-1) ) * 0.25D0
       FLXij(i,j,k,2,2) = GAMMA * Sij(i,j,k,2,2) * ( nurho(i  ,j  ,k  ) + nurho(i-1,j-1,k-1) &
                                                   + nurho(i-1,j  ,k  ) + nurho(i-1,j-1,k  ) &
                                                   + nurho(i  ,j-1,k  ) + nurho(i  ,j-1,k-1) &
                                                   + nurho(i  ,j  ,k-1) + nurho(i-1,j  ,k-1) ) * 0.25D0
       FLXij(i,j,k,3,3) = GAMMA * Sij(i,j,k,3,3) * ( nurho(i  ,j  ,k  ) + nurho(i-1,j-1,k-1) &
                                                   + nurho(i-1,j  ,k  ) + nurho(i-1,j-1,k  ) &
                                                   + nurho(i  ,j-1,k  ) + nurho(i  ,j-1,k-1) &
                                                   + nurho(i  ,j  ,k-1) + nurho(i-1,j  ,k-1) ) * 0.25D0
    enddo
    enddo
    enddo
    FLXij(:,:,KE+1,3,3) = FLXij(:,:,KE,3,3)

    do k = WS+1, WE-1
    do j = JS-1, JE
    do i = IS-1, IE
       FLXt(i,j,k,1) = RPr(i,j,k) * ( pott(i,j,k)-pott(i+1,j,k) ) * RDX                             &
                     * 0.25D0 * ( nurho(i,j,k) + nurho(i,j-1,k) + nurho(i,j,k-1) + nurho(i,j-1,k-1) )
       FLXt(i,j,k,2) = RPr(i,j,k) * ( pott(i,j,k)-pott(i,j+1,k) ) * RDY                             &
                     * 0.25D0 * ( nurho(i,j,k) + nurho(i-1,j,k) + nurho(i,j,k-1) + nurho(i-1,j,k-1) )
       FLXt(i,j,k,3) = RPr(i,j,k) * ( pott(i,j,k)-pott(i,j,k+1) ) * RDZ                             &
                     * 0.25D0 * ( nurho(i,j,k) + nurho(i-1,j,k) + nurho(i,j-1,k) + nurho(i-1,j-1,k) )

       FLXq(i,j,k,1,:) = RPr(i,j,k) * ( qtrc(i,j,k,:)-qtrc(i+1,j,k,:) ) * RDX                         &
                       * 0.25D0 * ( nurho(i,j,k) + nurho(i,j-1,k) + nurho(i,j,k-1) + nurho(i,j-1,k-1) )
       FLXq(i,j,k,2,:) = RPr(i,j,k) * ( qtrc(i,j,k,:)-qtrc(i,j+1,k,:) ) * RDY                         &
                       * 0.25D0 * ( nurho(i,j,k) + nurho(i-1,j,k) + nurho(i,j,k-1) + nurho(i-1,j,k-1) )
       FLXq(i,j,k,3,:) = RPr(i,j,k) * ( qtrc(i,j,k,:)-qtrc(i,j,k+1,:) ) * RDZ                         &
                       * 0.25D0 * ( nurho(i,j,k) + nurho(i-1,j,k) + nurho(i,j-1,k) + nurho(i-1,j-1,k) )
    enddo
    enddo
    enddo
    FLXt(:,:,WS,:)   = 0.D0
    FLXt(:,:,WE,:)   = 0.D0
    FLXq(:,:,WS,:,:) = 0.D0
    FLXq(:,:,WE,:,:) = 0.D0

    ! Merge Surface Flux
    FLXij(:,:,WS,1,3) = FLXij_sfc(:,:,1)
    FLXij(:,:,WS,2,3) = FLXij_sfc(:,:,2)
    FLXij(:,:,WS,3,3) = FLXij_sfc(:,:,3)
    FLXt (:,:,WS,3)   = FLXt_sfc (:,:)
    FLXq (:,:,WS,3,1) = FLXqv_sfc(:,:)

    ! tendency
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       momx_t(i,j,k) = - ( FLXij(i+1,j,k,1,1)-FLXij(i,j  ,k  ,1,1) ) * RDX &
                       + ( FLXij(i  ,j,k,1,2)-FLXij(i,j-1,k  ,1,2) ) * RDY &
                       + ( FLXij(i  ,j,k,1,3)-FLXij(i,j  ,k-1,1,3) ) * RDZ
       momy_t(i,j,k) = - ( FLXij(i,j  ,k,2,1)-FLXij(i-1,j,k  ,2,1) ) * RDX &
                       + ( FLXij(i,j+1,k,2,2)-FLXij(i  ,j,k  ,2,2) ) * RDY &
                       + ( FLXij(i,j  ,k,2,3)-FLXij(i  ,j,k-1,2,3) ) * RDZ
       momz_t(i,j,k) = - ( FLXij(i,j,k  ,3,1)-FLXij(i-1,j  ,k,3,1) ) * RDX &
                       + ( FLXij(i,j,k  ,3,2)-FLXij(i  ,j-1,k,3,2) ) * RDY &
                       + ( FLXij(i,j,k+1,3,3)-FLXij(i  ,j  ,k,3,3) ) * RDZ

       pott_t(i,j,k) = - ( ( FLXt(i,j,k,1)-FLXt(i-1,j  ,k  ,1) ) * RDX &
                         + ( FLXt(i,j,k,2)-FLXt(i  ,j-1,k  ,2) ) * RDY &
                         + ( FLXt(i,j,k,3)-FLXt(i  ,j  ,k-1,3) ) * RDZ &
                         ) / dens(i,j,k)

       qtrc_t(i,j,k,:) = - ( ( FLXq(i,j,k,1,:)-FLXq(i-1,j  ,k  ,1,:) ) * RDX &
                           + ( FLXq(i,j,k,2,:)-FLXq(i  ,j-1,k  ,2,:) ) * RDY &
                           + ( FLXq(i,j,k,3,:)-FLXq(i  ,j  ,k-1,3,:) ) * RDZ &
                           ) / dens(i,j,k)
    enddo
    enddo
    enddo

    deallocate( Uij   )
    deallocate( Sij   )
    deallocate( nu    )
    deallocate( nurho )
    deallocate( RPr   )
    deallocate( FLXij )
    deallocate( FLXt  )
    deallocate( FLXq  )

    return
  end subroutine ATMOS_PHY_TB

end module mod_atmos_phy_tb
