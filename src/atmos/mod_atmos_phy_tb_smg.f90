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
!! @li      2011-11-29 (S.Iga)     [new]
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
       Rdry   => CONST_Rdry,   &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd,   &
       CONST_UNDEF8
    use mod_grid, only : &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KA   => GRID_KA,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       WS   => GRID_WS,   &
       WE   => GRID_WE,   &
       DXYZ => GRID_DXYZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RCDZ => GRID_RCDZ
    use mod_atmos_vars, only: &
       QA => A_QA
    implicit none

    ! prognostic value
    real(8), intent(in)  :: dens(KA,IA,JA)      ! density [kg/m3]
    real(8), intent(in)  :: velx(KA,IA,JA)      ! velocity(x) [m/s]
    real(8), intent(in)  :: vely(KA,IA,JA)      ! velocity(y) [m/s]
    real(8), intent(in)  :: velz(KA,IA,JA)      ! velocity(z) [m/s]
    real(8), intent(in)  :: pott(KA,IA,JA)      ! potential temperature [K]
    real(8), intent(in)  :: qtrc(KA,IA,JA,QA)   ! tracer mixing ratio [kg/kg],[1/m3]

    ! surface flux
    real(8), intent(in)  :: FLXij_sfc(IA,JA,3)  ! => FLXij(1:IA,1:JA,WS,1:3,3)
    real(8), intent(in)  :: FLXt_sfc (IA,JA)    ! => FLXt (1:IA,1:JA,WS)
    real(8), intent(in)  :: FLXqv_sfc(IA,JA)    ! => FLXq (1:IA,1:JA,WS,1)

    ! prognostic tendency
    real(8), intent(out) :: momx_t(KA,IA,JA)
    real(8), intent(out) :: momy_t(KA,IA,JA)
    real(8), intent(out) :: momz_t(KA,IA,JA)
    real(8), intent(out) :: pott_t(KA,IA,JA)
    real(8), intent(out) :: qtrc_t(KA,IA,JA,QA)

    real(8) :: Uij(KA,IA,JA,3,3) ! velocity gradient matrix
    real(8) :: Sij(KA,IA,JA,3,3) ! deformation rate tensor

    real(8) :: nurho(KA,IA,JA)   ! nu * dens

    real(8) :: FLXij(KA,IA,JA,3,3)
    real(8) :: FLXt (KA,IA,JA,3)
    real(8) :: FLXq (KA,IA,JA,3,QA)

    real(8) :: Ri            ! Richardson number
    real(8) :: Pr            ! Prantle number
    real(8) :: RPr(KA,IA,JA) ! 1 / Pr

    integer :: i, j, k
    !---------------------------------------------------------------------------

    ! ii == jj : at the center of control volume (cube)
    do j = JS-2, JE+2
    do i = IS-1, IE+2
    do k = KS,   KE
       Uij(k,i,j,1,1) = ( velx(k,i,j)-velx(k,i-1,j) ) * RCDX(i) ! du/dx, (x, y, layer)
    enddo
    enddo
    enddo
    do j = JS-1, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       Uij(k,i,j,2,2) = ( vely(k,i,j)-vely(k,i,j-1) ) * RCDY(j) ! dv/dy, (x, y, layer)
    enddo
    enddo
    enddo
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       Uij(k,i,j,3,3) = ( velz(k,i,j)-velz(k-1,i,j) ) * RCDZ(k) ! dw/dz, (x, y, layer)
    enddo
    enddo
    enddo
    ! ii /= jj : on the edge   of control volume (cube)
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = KS,   KE
       Uij(k,i,j,1,2) = ( vely(k,i+1,j)-vely(k,i,j) ) * RCDX(i) ! dv/dx, (u, v, layer)
    enddo
    enddo
    enddo
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = WS+1, WE-1
       Uij(k,i,j,1,3) = ( velz(k,i+1,j)-velz(k,i,j) ) * RCDX(i) ! dw/dx, (u, y, interface)
    enddo
    enddo
    enddo
    Uij(WS,:,:,1,3) = 0.D0
    Uij(WE,:,:,1,3) = 0.D0
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = KS,   KE
       Uij(k,i,j,2,1) = ( velx(k,i,j+1)-velx(k,i,j) ) * RCDY(j) ! du/dy, (u, v, layer)
    enddo
    enddo
    enddo
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = WS+1, WE-1
       Uij(k,i,j,2,3) = ( velz(k,i,j+1)-velz(k,i,j) ) * RCDY(j) ! dw/dy, (x, v, interface)
    enddo
    enddo
    enddo
    Uij(WS,:,:,2,3) = 0.D0
    Uij(WE,:,:,2,3) = 0.D0
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = WS+1, WE-1
       Uij(k,i,j,3,1) = ( velx(k+1,i,j)-velx(k,i,j) ) * RCDZ(k) ! du/dz, (u, y, interface)
       Uij(k,i,j,3,2) = ( vely(k+1,i,j)-vely(k,i,j) ) * RCDZ(k) ! dv/dz, (x, v, interface)
    enddo
    enddo
    enddo
    Uij(WS,:,:,3,1) = 0.D0
    Uij(WS,:,:,3,1) = 0.D0
    Uij(WE,:,:,3,2) = 0.D0
    Uij(WE,:,:,3,2) = 0.D0

    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       Sij(k,i,j,1,1) = Uij(k,i,j,1,1) * dens(k,i,j)
       Sij(k,i,j,2,2) = Uij(k,i,j,2,2) * dens(k,i,j)
       Sij(k,i,j,3,3) = Uij(k,i,j,3,3) * dens(k,i,j)
    enddo
    enddo
    enddo
    do j = JS-2, JE+1
    do i = IS-2, IE+1
    do k = KS,   KE
       Sij(k,i,j,1,2) = 0.125D0 * ( Uij(k,i,j,1,2)+Uij(k,i,j,2,1) )               &
                      * ( dens(k,i,j)+dens(k,i+1,j)+dens(k,i,j+1)+dens(k,i+1,j+1) )
       Sij(k,i,j,2,1) = Sij(k,i,j,1,2)
    enddo
    enddo
    enddo
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = WS+1, WE-1
       Sij(k,i,j,1,3) = 0.125D0 * ( Uij(k,i,j,1,3)+Uij(k,i,j,3,1) )               &
                      * ( dens(k,i,j)+dens(k,i+1,j)+dens(k+1,i,j)+dens(k+1,i+1,j) )
       Sij(k,i,j,3,1) = Sij(k,i,j,1,3)
    enddo
    enddo
    enddo
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = WS+1, WE-1
       Sij(k,i,j,2,3) = 0.125D0 * ( Uij(k,i,j,2,3)+Uij(k,i,j,3,2) )               &
                      * ( dens(k,i,j)+dens(k,i,j+1)+dens(k+1,i,j)+dens(k+1,i,j+1) )
       Sij(k,i,j,3,2) = Sij(k,i,j,2,3)
    enddo
    enddo
    enddo

    ! calculate stratification effect (n2 and Pr are defined at w level)
!    do j = JS, JE
!    do i = IS, IE
!    do k = KS, KE
!       a = 1.D0 / ( 1.D0 + dqsl(k,i,j) * Lv / Cp )
!
!       c = (1+0.61*qw(k,i,j)-1.61*ql(k,i,j))* Lv/Cp /T_ov_pott(k,i,j) - 1.61 * pott(k,i,j)
!
!       betatheta= 1 + 0.61 * qw(k,i,j) - 1.61*ql(k,i,j) - a*c*dqsl(k,i,j)* T_ov_pott(k,i,j)
!
!       betaq= 0.61 * pott(k,i,j) + a * c
!
!       n2(k,i,j) = grav / (sum(pottl(k,i,j:kz+1))/2) * &
!                     (sum(betatheta * ( pottl(k+1,i,j)-pottl(k,i,j) )/dz + &
!                     betaq * ( qw(k+1,i,j)-qw(k,i,j)  )/dz )/4)
!    enddo
!    enddo
!    enddo

    ! calculate  Pr (Pr is defined at w level)
    do j = JS-1, JE
    do i = IS-1, IE
    do k = WS+1, WE-1
       ! Ri = g / rho * drho/dz / dU/dz
       Ri = GRAV * ( dens(k+1,i,j)-dens(k,i,j) ) * RCDZ(k) * 2.D0 / ( dens(k+1,i,j)+dens(k,i,j) )                   &
          / max( ( ( velx(k+1,i,j)+velx(k+1,i-1,j)-velx(k,i,j)-velx(k,i-1,j) ) * 0.5D0 * RCDZ(k) )**2.D0            &
               + ( ( vely(k+1,i,j)+vely(k+1,i,j-1)-vely(k,i,j)-vely(k,i,j-1) ) * 0.5D0 * RCDZ(k) )**2.D0, Ustar_min )

       Ri         = min( 0.25D0, max( 0.D0, Ri ) )
       Pr         = ( 1.D0-4.D0*Ri ) / 3.D0 + 4.D0*Ri
       RPr(k,i,j) = 1.D0 / Pr
    enddo
    enddo
    enddo
    RPr(WS,:,:) = 0.D0
    RPr(WE,:,:) = 0.D0

    do j = JS-2, JE+1
    do i = IS-2, IE+1
    do k = WS+1, WE-1
       nurho(k,i,j) = ( Cs * DXYZ ) ** 2.D0 * ( sqrt( 2.D0 * GAMMA                          &
                    * ( ( Sij(k  ,i  ,j  ,1,1) + Sij(k+1,i+1,j+1,1,1)                       &
                        + Sij(k  ,i+1,j  ,1,1) + Sij(k  ,i+1,j+1,1,1)                       &
                        + Sij(k  ,i  ,j+1,1,1) + Sij(k+1,i  ,j+1,1,1)                       &
                        + Sij(k+1,i  ,j  ,1,1) + Sij(k+1,i+1,j  ,1,1) ) * 0.125D0 ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,1,2) + Sij(k+1,i  ,j  ,1,2) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,1,3) + Sij(k  ,i  ,j+1,1,3) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,2,1) + Sij(k+1,i  ,j  ,2,1) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,2,2) + Sij(k+1,i+1,j+1,2,2)                       &
                        + Sij(k  ,i+1,j  ,2,2) + Sij(k  ,i+1,j+1,2,2)                       &
                        + Sij(k  ,i  ,j+1,2,2) + Sij(k+1,i  ,j+1,2,2)                       &
                        + Sij(k+1,i  ,j  ,2,2) + Sij(k+1,i+1,j  ,2,2) ) * 0.125D0 ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,2,3) + Sij(k  ,i+1,j  ,2,3) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,3,1) + Sij(k  ,i  ,j+1,3,1) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,3,2) + Sij(k  ,i+1,j  ,3,2) ) * 0.5D0   ) ** 2.D0 &
                    + ( ( Sij(k  ,i  ,j  ,3,3) + Sij(k+1,i+1,j+1,3,3)                       &
                        + Sij(k  ,i+1,j  ,3,3) + Sij(k  ,i+1,j+1,3,3)                       &
                        + Sij(k  ,i  ,j+1,3,3) + Sij(k+1,i  ,j+1,3,3)                       &
                        + Sij(k+1,i  ,j  ,3,3) + Sij(k+1,i+1,j  ,3,3) ) * 0.125D0 ) ** 2.D0 &
                    ) )
    enddo
    enddo
    enddo
    nurho(WS,:,:) = 0.D0
    nurho(WE,:,:) = 0.D0

    do k = KS,   KE
    do j = JS-1, JE
    do i = IS-1, IE
       FLXij(k,i,j,1,2) = GAMMA * Sij(k,i,j,1,2) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i  ,j  ) )
       FLXij(k,i,j,1,3) = GAMMA * Sij(k,i,j,1,3) * ( nurho(k  ,i  ,j  ) + nurho(k  ,i  ,j-1) )
       FLXij(k,i,j,2,1) = GAMMA * Sij(k,i,j,2,1) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i  ,j  ) )
       FLXij(k,i,j,2,3) = GAMMA * Sij(k,i,j,2,3) * ( nurho(k  ,i  ,j  ) + nurho(k  ,i-1,j  ) )
       FLXij(k,i,j,3,1) = GAMMA * Sij(k,i,j,3,1) * ( nurho(k  ,i  ,j  ) + nurho(k  ,i  ,j-1) )
       FLXij(k,i,j,3,2) = GAMMA * Sij(k,i,j,3,2) * ( nurho(k  ,i  ,j  ) + nurho(k  ,i-1,j  ) )
    enddo
    enddo
    enddo

    do k = KS, KE
    do j = JS, JE+1
    do i = IS, IE+1
       FLXij(k,i,j,1,1) = GAMMA * Sij(k,i,j,1,1) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i-1,j-1) &
                                                   + nurho(k  ,i-1,j  ) + nurho(k  ,i-1,j-1) &
                                                   + nurho(k  ,i  ,j-1) + nurho(k-1,i  ,j-1) &
                                                   + nurho(k-1,i  ,j  ) + nurho(k-1,i-1,j  ) ) * 0.25D0
       FLXij(k,i,j,2,2) = GAMMA * Sij(k,i,j,2,2) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i-1,j-1) &
                                                   + nurho(k  ,i-1,j  ) + nurho(k  ,i-1,j-1) &
                                                   + nurho(k  ,i  ,j-1) + nurho(k-1,i  ,j-1) &
                                                   + nurho(k-1,i  ,j  ) + nurho(k-1,i-1,j  ) ) * 0.25D0
       FLXij(k,i,j,3,3) = GAMMA * Sij(k,i,j,3,3) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i-1,j-1) &
                                                   + nurho(k  ,i-1,j  ) + nurho(k  ,i-1,j-1) &
                                                   + nurho(k  ,i  ,j-1) + nurho(k-1,i  ,j-1) &
                                                   + nurho(k-1,i  ,j  ) + nurho(k-1,i-1,j  ) ) * 0.25D0
    enddo
    enddo
    enddo
    FLXij(KE+1,:,:,3,3) = FLXij(KE,:,:,3,3)

    do k = WS+1, WE-1
    do j = JS-1, JE
    do i = IS-1, IE
       FLXt(k,i,j,1) = RPr(k,i,j) * ( pott(k,i,j)-pott(k,i+1,j) ) * RCDX(i)                         &
                     * 0.25D0 * ( nurho(k,i,j) + nurho(k,i,j-1) + nurho(k-1,i,j) + nurho(k-1,i,j-1) )
       FLXt(k,i,j,2) = RPr(k,i,j) * ( pott(k,i,j)-pott(k,i,j+1) ) * RCDY(j)                         &
                     * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k-1,i,j) + nurho(k-1,i-1,j) )
       FLXt(k,i,j,3) = RPr(k,i,j) * ( pott(k,i,j)-pott(k+1,i,j) ) * RCDZ(k)                         &
                     * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k,i,j-1) + nurho(k,i-1,j-1) )

       FLXq(k,i,j,1,:) = RPr(k,i,j) * ( qtrc(k,i,j,:)-qtrc(k,i+1,j,:) ) * RCDX(i)                     &
                       * 0.25D0 * ( nurho(k,i,j) + nurho(k,i,j-1) + nurho(k-1,i,j) + nurho(k-1,i,j-1) )
       FLXq(k,i,j,2,:) = RPr(k,i,j) * ( qtrc(k,i,j,:)-qtrc(k,i,j+1,:) ) * RCDY(j)                     &
                       * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k-1,i,j) + nurho(k-1,i-1,j) )
       FLXq(k,i,j,3,:) = RPr(k,i,j) * ( qtrc(k,i,j,:)-qtrc(k+1,i,j,:) ) * RCDZ(k)                     &
                       * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k,i,j-1) + nurho(k,i-1,j-1) )
    enddo
    enddo
    enddo
    FLXt(WS,:,:,:)   = 0.D0
    FLXt(WE,:,:,:)   = 0.D0
    FLXq(WS,:,:,:,:) = 0.D0
    FLXq(WE,:,:,:,:) = 0.D0

    ! Merge Surface Flux
    FLXij(WS,:,:,1,3) = FLXij_sfc(:,:,1)
    FLXij(WS,:,:,2,3) = FLXij_sfc(:,:,2)
    FLXij(WS,:,:,3,3) = FLXij_sfc(:,:,3)
    FLXt (WS,:,:,3)   = FLXt_sfc (:,:)
    FLXq (WS,:,:,3,1) = FLXqv_sfc(:,:)

    ! tendency
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       momx_t(k,i,j) = - ( FLXij(k,i+1,j,1,1)-FLXij(k  ,i,j  ,1,1) ) * RCDX(i) &
                       + ( FLXij(k,i,  j,1,2)-FLXij(k  ,i,j-1,1,2) ) * RCDY(j) &
                       + ( FLXij(k,i,  j,1,3)-FLXij(k-1,i,j  ,1,3) ) * RCDZ(k)
       momy_t(k,i,j) = - ( FLXij(k,i,j,  2,1)-FLXij(k  ,i-1,j,2,1) ) * RCDX(i) &
                       + ( FLXij(k,i,j+1,2,2)-FLXij(k  ,i  ,j,2,2) ) * RCDY(j) &
                       + ( FLXij(k,i,j,  2,3)-FLXij(k-1,i  ,j,2,3) ) * RCDZ(k)
       momz_t(k,i,j) = - ( FLXij(k,  i,j,3,1)-FLXij(k,i-1,j  ,3,1) ) * RCDX(i) &
                       + ( FLXij(k,  i,j,3,2)-FLXij(k,i  ,j-1,3,2) ) * RCDY(j) &
                       + ( FLXij(k+1,i,j,3,3)-FLXij(k,i  ,j  ,3,3) ) * RCDZ(k)

       pott_t(k,i,j) = - ( ( FLXt(k,i,j,1)-FLXt(k  ,i-1,j  ,1) ) * RCDX(i) &
                         + ( FLXt(k,i,j,2)-FLXt(k  ,i  ,j-1,2) ) * RCDY(j) &
                         + ( FLXt(k,i,j,3)-FLXt(k-1,i  ,j  ,3) ) * RCDZ(k) &
                         ) / dens(k,i,j)

       qtrc_t(k,i,j,:) = - ( ( FLXq(k,i,j,1,:)-FLXq(k  ,i-1,j  ,1,:) ) * RCDX(i) &
                           + ( FLXq(k,i,j,2,:)-FLXq(k  ,i  ,j-1,2,:) ) * RCDY(j) &
                           + ( FLXq(k,i,j,3,:)-FLXq(k-1,i  ,j  ,3,:) ) * RCDZ(k) &
                           ) / dens(k,i,j)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB

end module mod_atmos_phy_tb
