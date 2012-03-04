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
  integer, private, parameter :: I_VELZ = 1
  integer, private, parameter :: I_VELX = 2
  integer, private, parameter :: I_VELY = 3
  integer, private, parameter :: I_POTT = 4

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
  subroutine ATMOS_PHY_TB( FLXij_sfc, FLXt_sfc, FLXqv_sfc )
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd,   &
       CONST_UNDEF8
    use mod_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_TB
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait,  &
       COMM_total
    use mod_grid, only : &
       KA   => GRID_KA,   &
       IA   => GRID_IA,   &
       JA   => GRID_JA,   &
       KS   => GRID_KS,   &
       KE   => GRID_KE,   &
       IS   => GRID_IS,   &
       IE   => GRID_IE,   &
       JS   => GRID_JS,   &
       JE   => GRID_JE,   &
       DXYZ => GRID_DXYZ, &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use mod_atmos_vars, only: &
       var => atmos_var, &
       A_NAME,      &
       VA  => A_VA, &
       QA  => A_QA, &
       I_DENS,      &
       I_MOMX,      &
       I_MOMY,      &
       I_MOMZ,      &
       I_RHOT
    implicit none

    ! surface flux
    real(8), intent(in)  :: FLXij_sfc(IA,JA,3)  ! => FLXij(1:IA,1:JA,KS-1,1:3,3)
    real(8), intent(in)  :: FLXt_sfc (IA,JA)    ! => FLXt (1:IA,1:JA,KS-1)
    real(8), intent(in)  :: FLXqv_sfc(IA,JA)    ! => FLXq (1:IA,1:JA,KS-1,1)

    ! work
    real(8) :: diagvar(KA,IA,JA,4)   ! diagnostic variables (work)
    real(8) :: Uij    (KA,IA,JA,3,3) ! velocity gradient matrix
    real(8) :: Sij    (KA,IA,JA,3,3) ! deformation rate tensor

    real(8) :: nurho(KA,IA,JA)   ! nu * dens

    real(8) :: FLXij(KA,IA,JA,3,3)
    real(8) :: FLXt (KA,IA,JA,3)
    real(8) :: FLXq (KA,IA,JA,3,QA)

    real(8) :: Ri            ! Richardson number
    real(8) :: Pr            ! Prantle number
    real(8) :: RPr(KA,IA,JA) ! 1 / Pr

    integer :: k, i, j, iq, iv
    !---------------------------------------------------------------------------

    ! momentum -> velocity
    do j = JS, JE+1
    do i = IS, IE+1
    do k = KS, KE-1
       diagvar(k,i,j,I_VELZ) = 2.D0 * var(k,i,j,I_MOMZ) &
                             / ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) )
    enddo
    enddo
    enddo
    !OCL XFILL
    do j = JS, JE+1
    do i = IS, IE+1
       diagvar(KS-1,i,j,I_VELZ) = 0.D0
       diagvar(KE  ,i,j,I_VELZ) = 0.D0
    enddo
    enddo

    do j = JS,   JE+1
    do i = IS-1, IE+1
    do k = KS,   KE
       diagvar(k,i,j,I_VELX) = 2.D0 * var(k,i,j,I_MOMX) &
                             / ( var(k,i+1,j,I_DENS)+var(k,i,j,I_DENS) )
    enddo
    enddo
    enddo

    do j = JS-1, JE+1
    do i = IS,   IE+1
    do k = KS,   KE
       diagvar(k,i,j,I_VELY) = 2.D0 * var(k,i,j,I_MOMY) &
                             / ( var(k,i,j+1,I_DENS)+var(k,i,j,I_DENS) )
    enddo
    enddo
    enddo

    ! potential temperature
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       diagvar(k,i,j,I_POTT) = var(k,i,j,I_RHOT) / var(k,i,j,I_DENS) 
    enddo
    enddo
    enddo

    ! ii == jj : at the center of control volume (cube)
    do j = JS-2, JE+2
    do i = IS-1, IE+2
    do k = KS,   KE
       Uij(k,i,j,1,1) = ( diagvar(k,i,j,I_VELX)-diagvar(k,i-1,j,I_VELX) ) * RCDX(i) ! du/dx, (x, y, layer)
    enddo
    enddo
    enddo
    do j = JS-1, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       Uij(k,i,j,2,2) = ( diagvar(k,i,j,I_VELY)-diagvar(k,i,j-1,I_VELY) ) * RCDY(j) ! dv/dy, (x, y, layer)
    enddo
    enddo
    enddo
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       Uij(k,i,j,3,3) = ( diagvar(k,i,j,I_VELZ)-diagvar(k-1,i,j,I_VELZ) ) * RCDZ(k) ! dw/dz, (x, y, layer)
    enddo
    enddo
    enddo
    ! ii /= jj : on the edge   of control volume (cube)
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = KS,   KE
       Uij(k,i,j,1,2) = ( diagvar(k,i+1,j,I_VELY)-diagvar(k,i,j,I_VELY) ) * RCDX(i) ! dv/dx, (u, v, layer)
    enddo
    enddo
    enddo
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = KS, KE-1
       Uij(k,i,j,1,3) = ( diagvar(k,i+1,j,I_VELZ)-diagvar(k,i,j,I_VELZ) ) * RCDX(i) ! dw/dx, (u, y, interface)
    enddo
    enddo
    enddo
    Uij(KS-1,:,:,1,3) = 0.D0
    Uij(KE,:,:,1,3) = 0.D0
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = KS,   KE
       Uij(k,i,j,2,1) = ( diagvar(k,i,j+1,I_VELX)-diagvar(k,i,j,I_VELX) ) * RCDY(j) ! du/dy, (u, v, layer)
    enddo
    enddo
    enddo
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = KS, KE-1
       Uij(k,i,j,2,3) = ( diagvar(k,i,j+1,I_VELZ)-diagvar(k,i,j,I_VELZ) ) * RCDY(j) ! dw/dy, (x, v, interface)
    enddo
    enddo
    enddo
    Uij(KS-1,:,:,2,3) = 0.D0
    Uij(KE,:,:,2,3) = 0.D0
    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS, KE-1
       Uij(k,i,j,3,1) = ( diagvar(k+1,i,j,I_VELX)-diagvar(k,i,j,I_VELX) ) * RCDZ(k) ! du/dz, (u, y, interface)
       Uij(k,i,j,3,2) = ( diagvar(k+1,i,j,I_VELY)-diagvar(k,i,j,I_VELY) ) * RCDZ(k) ! dv/dz, (x, v, interface)
    enddo
    enddo
    enddo
    Uij(KS-1,:,:,3,1) = 0.D0
    Uij(KS-1,:,:,3,1) = 0.D0
    Uij(KE,:,:,3,2) = 0.D0
    Uij(KE,:,:,3,2) = 0.D0

    do j = JS-2, JE+2
    do i = IS-2, IE+2
    do k = KS,   KE
       Sij(k,i,j,1,1) = Uij(k,i,j,1,1) * var(k,i,j,I_DENS)
       Sij(k,i,j,2,2) = Uij(k,i,j,2,2) * var(k,i,j,I_DENS)
       Sij(k,i,j,3,3) = Uij(k,i,j,3,3) * var(k,i,j,I_DENS)
    enddo
    enddo
    enddo
    do j = JS-2, JE+1
    do i = IS-2, IE+1
    do k = KS,   KE
       Sij(k,i,j,1,2) = 0.125D0 * ( Uij(k,i,j,1,2)+Uij(k,i,j,2,1) )   &
                      * ( var(k,i,j  ,I_DENS) + var(k,i+1,j  ,I_DENS) &
                        + var(k,i,j+1,I_DENS) + var(k,i+1,j+1,I_DENS) )
       Sij(k,i,j,2,1) = Sij(k,i,j,1,2)
    enddo
    enddo
    enddo
    do j = JS-2, JE+2
    do i = IS-2, IE+1
    do k = KS, KE-1
       Sij(k,i,j,1,3) = 0.125D0 * ( Uij(k,i,j,1,3)+Uij(k,i,j,3,1) )               &
                      * ( var(k  ,i,j,I_DENS) + var(k  ,i+1,j,I_DENS) &
                        + var(k+1,i,j,I_DENS) + var(k+1,i+1,j,I_DENS) )
       Sij(k,i,j,3,1) = Sij(k,i,j,1,3)
    enddo
    enddo
    enddo
    do j = JS-2, JE+1
    do i = IS-2, IE+2
    do k = KS, KE-1
       Sij(k,i,j,2,3) = 0.125D0 * ( Uij(k,i,j,2,3)+Uij(k,i,j,3,2) )               &
                      * ( var(k  ,i,j,I_DENS) + var(k  ,i,j+1,I_DENS) &
                        + var(k+1,i,j,I_DENS) + var(k+1,i,j+1,I_DENS) )
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
    do k = KS, KE-1
       ! Ri = g / rho * drho/dz / dU/dz
       Ri = GRAV * ( var(k+1,i,j,I_DENS)-var(k,i,j,I_DENS) ) * RCDZ(k)                                     &
          * 2.D0 / ( var(k+1,i,j,I_DENS)+var(k,i,j,I_DENS) )                                               &
          / max( ( ( diagvar(k+1,i,j,I_VELX)+diagvar(k+1,i-1,j,I_VELX)                                     &
                   - diagvar(k  ,i,j,I_VELX)-diagvar(k  ,i-1,j,I_VELX) ) * 0.5D0 * RCDZ(k) )**2            &
               + ( ( diagvar(k+1,i,j,I_VELY)+diagvar(k+1,i,j-1,I_VELY)                                     &
                   - diagvar(k  ,i,j,I_VELY)-diagvar(k  ,i,j-1,I_VELY) ) * 0.5D0 * RCDZ(k) )**2, Ustar_min )

       Ri         = min( 0.25D0, max( 0.D0, Ri ) )
       Pr         = ( 1.D0-4.D0*Ri ) / 3.D0 + 4.D0*Ri
       RPr(k,i,j) = 1.D0 / Pr
    enddo
    enddo
    enddo
    RPr(KS-1,:,:) = 0.D0
    RPr(KE,:,:) = 0.D0

    do j = JS-2, JE+1
    do i = IS-2, IE+1
    do k = KS, KE-1
       nurho(k,i,j) = ( Cs * DXYZ ) **2 * ( sqrt( 2.D0 * GAMMA                          &
                    * ( ( Sij(k  ,i  ,j  ,1,1) + Sij(k+1,i+1,j+1,1,1)                   &
                        + Sij(k  ,i+1,j  ,1,1) + Sij(k  ,i+1,j+1,1,1)                   &
                        + Sij(k  ,i  ,j+1,1,1) + Sij(k+1,i  ,j+1,1,1)                   &
                        + Sij(k+1,i  ,j  ,1,1) + Sij(k+1,i+1,j  ,1,1) ) * 0.125D0 ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,1,2) + Sij(k+1,i  ,j  ,1,2) ) * 0.5D0   ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,1,3) + Sij(k  ,i  ,j+1,1,3) ) * 0.5D0   ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,2,1) + Sij(k+1,i  ,j  ,2,1) ) * 0.5D0   ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,2,2) + Sij(k+1,i+1,j+1,2,2)                   &
                        + Sij(k  ,i+1,j  ,2,2) + Sij(k  ,i+1,j+1,2,2)                   &
                        + Sij(k  ,i  ,j+1,2,2) + Sij(k+1,i  ,j+1,2,2)                   &
                        + Sij(k+1,i  ,j  ,2,2) + Sij(k+1,i+1,j  ,2,2) ) * 0.125D0 ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,2,3) + Sij(k  ,i+1,j  ,2,3) ) * 0.5D0   ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,3,1) + Sij(k  ,i  ,j+1,3,1) ) * 0.5D0   ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,3,2) + Sij(k  ,i+1,j  ,3,2) ) * 0.5D0   ) **2 &
                    + ( ( Sij(k  ,i  ,j  ,3,3) + Sij(k+1,i+1,j+1,3,3)                   &
                        + Sij(k  ,i+1,j  ,3,3) + Sij(k  ,i+1,j+1,3,3)                   &
                        + Sij(k  ,i  ,j+1,3,3) + Sij(k+1,i  ,j+1,3,3)                   &
                        + Sij(k+1,i  ,j  ,3,3) + Sij(k+1,i+1,j  ,3,3) ) * 0.125D0 ) **2 &
                    ) )
    enddo
    enddo
    enddo
    nurho(KS-1,:,:) = 0.D0
    nurho(KE,:,:) = 0.D0

    do j = JS-1, JE
    do i = IS-1, IE
    do k = KS,   KE
       FLXij(k,i,j,1,2) = Sij(k,i,j,1,2) * ( nurho(k,i,j) + nurho(k-1,i  ,j  ) )
       FLXij(k,i,j,1,3) = Sij(k,i,j,1,3) * ( nurho(k,i,j) + nurho(k  ,i  ,j-1) )
       FLXij(k,i,j,2,1) = Sij(k,i,j,2,1) * ( nurho(k,i,j) + nurho(k-1,i  ,j  ) )
       FLXij(k,i,j,2,3) = Sij(k,i,j,2,3) * ( nurho(k,i,j) + nurho(k  ,i-1,j  ) )
       FLXij(k,i,j,3,1) = Sij(k,i,j,3,1) * ( nurho(k,i,j) + nurho(k  ,i  ,j-1) )
       FLXij(k,i,j,3,2) = Sij(k,i,j,3,2) * ( nurho(k,i,j) + nurho(k  ,i-1,j  ) )
    enddo
    enddo
    enddo

    do j = JS, JE+1
    do i = IS, IE+1
    do k = KS, KE
       FLXij(k,i,j,1,1) = Sij(k,i,j,1,1) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i-1,j-1) &
                                           + nurho(k  ,i-1,j  ) + nurho(k  ,i-1,j-1) &
                                           + nurho(k  ,i  ,j-1) + nurho(k-1,i  ,j-1) &
                                           + nurho(k-1,i  ,j  ) + nurho(k-1,i-1,j  ) ) * 0.25D0
       FLXij(k,i,j,2,2) = Sij(k,i,j,2,2) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i-1,j-1) &
                                           + nurho(k  ,i-1,j  ) + nurho(k  ,i-1,j-1) &
                                           + nurho(k  ,i  ,j-1) + nurho(k-1,i  ,j-1) &
                                           + nurho(k-1,i  ,j  ) + nurho(k-1,i-1,j  ) ) * 0.25D0
       FLXij(k,i,j,3,3) = Sij(k,i,j,3,3) * ( nurho(k  ,i  ,j  ) + nurho(k-1,i-1,j-1) &
                                           + nurho(k  ,i-1,j  ) + nurho(k  ,i-1,j-1) &
                                           + nurho(k  ,i  ,j-1) + nurho(k-1,i  ,j-1) &
                                           + nurho(k-1,i  ,j  ) + nurho(k-1,i-1,j  ) ) * 0.25D0
    enddo
    enddo
    enddo
    do j = JS, JE+1
    do i = IS, IE+1
       FLXij(KE+1,i,j,3,3) = FLXij(KE,i,j,3,3)
    enddo
    enddo

    do j = JS-1, JE
    do i = IS-1, IE
    do k = KS,   KE-1
       FLXt(k,i,j,1) = RPr(k,i,j) * ( diagvar(k,i,j,I_POTT)-diagvar(k,i+1,j,I_POTT) ) * RCDX(i)     &
                     * 0.25D0 * ( nurho(k,i,j) + nurho(k,i,j-1) + nurho(k-1,i,j) + nurho(k-1,i,j-1) )
       FLXt(k,i,j,2) = RPr(k,i,j) * ( diagvar(k,i,j,I_POTT)-diagvar(k,i,j+1,I_POTT) ) * RCDY(j)     &
                     * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k-1,i,j) + nurho(k-1,i-1,j) )
       FLXt(k,i,j,3) = RPr(k,i,j) * ( diagvar(k,i,j,I_POTT)-diagvar(k+1,i,j,I_POTT) ) * RCDZ(k)     &
                     * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k,i,j-1) + nurho(k,i-1,j-1) )
    enddo
    enddo
    enddo
    do j = JS, JE+1
    do i = IS, IE+1
       FLXt(KS-1,i,j,:)   = 0.D0
       FLXt(KE  ,i,j,:)   = 0.D0
    enddo
    enddo

    do iq = 1,  QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       FLXq(k,i,j,1,iq) = RPr(k,i,j) * ( var(k,i,j,5+iq)-var(k,i+1,j,5+iq) ) * RCDX(i)                 &
                        * 0.25D0 * ( nurho(k,i,j) + nurho(k,i,j-1) + nurho(k-1,i,j) + nurho(k-1,i,j-1) )
       FLXq(k,i,j,2,iq) = RPr(k,i,j) * ( var(k,i,j,5+iq)-var(k,i,j+1,5+iq) ) * RCDY(j)                 &
                        * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k-1,i,j) + nurho(k-1,i-1,j) )
       FLXq(k,i,j,3,iq) = RPr(k,i,j) * ( var(k,i,j,5+iq)-var(k+1,i,j,5+iq) ) * RCDZ(k)                 &
                        * 0.25D0 * ( nurho(k,i,j) + nurho(k,i-1,j) + nurho(k,i,j-1) + nurho(k,i-1,j-1) )
    enddo
    enddo
    enddo
    enddo
    do iq = 1,  QA
    do j  = JS, JE
    do i  = IS, IE
       FLXq(KS-1,i,j,:,iq) = 0.D0
       FLXq(KE  ,i,j,:,iq) = 0.D0
    enddo
    enddo
    enddo

    ! Merge Surface Flux
    FLXij(KS-1,:,:,1,3) = FLXij_sfc(:,:,1)
    FLXij(KS-1,:,:,2,3) = FLXij_sfc(:,:,2)
    FLXij(KS-1,:,:,3,3) = FLXij_sfc(:,:,3)
    FLXt (KS-1,:,:,3)   = FLXt_sfc (:,:)
    FLXq (KS-1,:,:,3,1) = FLXqv_sfc(:,:)

    !--- update
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       var(k,i,j,I_MOMZ) = var(k,i,j,I_MOMZ) &
                         + dt * - ( ( FLXij(k,  i,j,3,1)-FLXij(k,i-1,j  ,3,1) ) * RCDX(i) &
                                  + ( FLXij(k,  i,j,3,2)-FLXij(k,i  ,j-1,3,2) ) * RCDY(j) &
                                  + ( FLXij(k+1,i,j,3,3)-FLXij(k,i  ,j  ,3,3) ) * RFDZ(k) )
       var(k,i,j,I_MOMX) = var(k,i,j,I_MOMX) &
                         + dt * - ( ( FLXij(k,i+1,j,1,1)-FLXij(k  ,i,j  ,1,1) ) * RFDX(i) &
                                  + ( FLXij(k,i,  j,1,2)-FLXij(k  ,i,j-1,1,2) ) * RCDY(j) &
                                  + ( FLXij(k,i,  j,1,3)-FLXij(k-1,i,j  ,1,3) ) * RCDZ(k) )
       var(k,i,j,I_MOMY) = var(k,i,j,I_MOMY) &
                         + dt * - ( ( FLXij(k,i,j,  2,1)-FLXij(k  ,i-1,j,2,1) ) * RCDX(i) &
                                  + ( FLXij(k,i,j+1,2,2)-FLXij(k  ,i  ,j,2,2) ) * RFDY(j) &
                                  + ( FLXij(k,i,j,  2,3)-FLXij(k-1,i  ,j,2,3) ) * RCDZ(k) )

       var(k,i,j,I_RHOT) = var(k,i,j,I_RHOT) &
                         + dt * - ( ( FLXt(k,i,j,1)-FLXt(k  ,i-1,j  ,1) ) * RCDX(i) &
                                  + ( FLXt(k,i,j,2)-FLXt(k  ,i  ,j-1,2) ) * RCDY(j) &
                                  + ( FLXt(k,i,j,3)-FLXt(k-1,i  ,j  ,3) ) * RCDZ(k) )
    enddo
    enddo
    enddo

    do iq = 1,  QA
    do k  = KS, KE
    do j  = JS, JE
    do i  = IS, IE
       var(k,i,j,5+iq) = var(k,i,j,5+iq) &
                       + dt * ( - ( ( FLXq(k,i,j,1,iq)-FLXq(k  ,i-1,j  ,1,iq) ) * RCDX(i)   &
                                  + ( FLXq(k,i,j,2,iq)-FLXq(k  ,i  ,j-1,2,iq) ) * RCDY(j)   &
                                  + ( FLXq(k,i,j,3,iq)-FLXq(k-1,i  ,j  ,3,iq) ) * RCDZ(k) ) &
                              ) / var(k,i,j,I_DENS)
    enddo
    enddo
    enddo
    enddo

    do iv = 1, VA
       call COMM_vars8( var(:,:,:,iv), iv )
       call COMM_wait ( var(:,:,:,iv), iv )
    enddo

    ! check total mass
    call COMM_total( var(:,:,:,:), A_NAME(:) )

    return
  end subroutine ATMOS_PHY_TB

end module mod_atmos_phy_tb
