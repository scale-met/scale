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
!! @li      2011-11-11 (H.Yashiro) [new] integrate
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

  real(8), parameter :: Cs    = 0.18D0 ! (Sullivan et al.1994, Nakanishi and Niino)
  real(8), parameter :: gamma = 1.D0 ! assume dx=dy=dz

  logical, parameter :: stratification_effect=.true.

  real(8), save      :: Pr ! Prantle num

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
  subroutine ATMOS_PHY_TB( dens,   momx,   momy,   momz,   lwpt,   &
                           qtrc,                                   &
                           pres,   velx,   vely,   velz,   temp,   &
                           dens_t, momx_t, momy_t, momz_t, lwpt_t, &
                           qtrc_t                                  )
    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    use mod_const, only : &
       GRAV   => CONST_GRAV,  &
       Rair   => CONST_Rair,  &
       CPovCV => CONST_CPovCV, &
       Pstd   => CONST_Pstd
    use mod_time, only: &
       TIME_DTSEC_ATMOS_PHY_TB
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
    real(8), intent(inout) :: qtrc(IA,JA,KA,QA) ! tracer mixing ratio   [kg/kg],[1/m3]
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

    use mod_stdio, only: &
       IO_FID_LOG,  &
       IO_L
    implicit none

    real(8), allocatable :: Uij(:,:,:,:,:) ! velocity gradient matrix
    real(8), allocatable :: Sij(:,:,:,:,:) ! deformation rate tensor

    real(8), allocatable :: nu   (:,:,:)   ! eddy-viscosity coefficient
    real(8), allocatable :: nurho(:,:,:)   ! nu * dens

    real(8), allocatable :: FLXij (:,:,:,:,:)
    real(8), allocatable :: FLXene(:,:,:,:)
    real(8), allocatable :: FLXq  (:,:,:,:,:)

    real(8) :: L ! mean grid interval
    real(8) :: Ri ! Rechardson num

    q(IS:IE,JS:JE,KS:KE,1:TRC_VMAX) = ! tracer
    qw(IS:IE,JS:JE,KS:KE) = !total water (vapor + liquid), 
    ql(IS:IE,JS:JE,KS:KE) = !liquid water 
    ! solid water is not considered y                                        et
    pottl(IS:IE,JS:JE,KS:KE) = !  liquid water potential temperature
    pott(IS:IE,JS:JE,KS:KE) = !  potential temperature
    dqsl(IS:IE,JS:JE,KS:KE) = ! dq(sat)/dT| at T=Tl

    T_ov_pott(IS:IE,JS:JE,KS:KE) = ! Exner func., i.e. (p0/p)^(R/Cp)
    Lv = ! latent heat

    integer :: i, j
    !---------------------------------------------------------------------------

    allocate( Uij   (IA,JA,KA,3,3)  )
    allocate( Sij   (IA,JA,KA,3,3)  )
    allocate( nu    (IA,JA,KA)      )
    allocate( nurho (IA,JA,KA)      )
    allocate( FLXij (IA,JA,KA,3,3)  )
    allocate( FLXene(IA,JA,KA,3)    )
    allocate( FLXq  (IA,JA,KA,3,QA) )

    L = GRID_DX

    ! no topography is assumed. if there are topography, velz_w(bottom) is not zero.
    ! bottom boundary
    !velx(:,:,WS) = velx_sfc(:,:)
    !vely(:,:,WS) = vely_sfc(:,:)
    !velz(:,:,WS) = velz_sfc(:,:)


    do k = KS, KE
    do j = JS, JE
    do i = IS, IE+1
       Uij(i,j,k,1,1) = ( velx(i  ,j  ,k  )-velx(i-1,j  ,k  ) ) * RDX ! du/dx, (x, y, layer)

    do k = KS, KE
    do j = JS-1, JE
    do i = IS, IE
       Uij(i,j,k,1,2) = ( vely(i+1,j  ,k  )-vely(i  ,j  ,k  ) ) * RDX ! dv/dx, (u, v, layer)

    do k = WS, WE
    do j = JS, JE
    do i = IS, IE
       Uij(i,j,k,1,3) = ( velz(i+1,j  ,k  )-velz(i  ,j  ,k  ) ) * RDX ! dw/dx, (u, y, interface)

    do k = KS, KE
    do j = JS, JE
    do i = IS-1, IE
       Uij(i,j,k,2,1) = ( velx(i  ,j+1,k  )-velx(i  ,j  ,k  ) ) * RDY ! du/dy, (u, v, layer)

    do k = KS, KE
    do j = JS, JE+1
    do i = IS, IE
       Uij(i,j,k,2,2) = ( vely(i  ,j+1,k  )-vely(i  ,j  ,k  ) ) * RDY ! dv/dy, (x, y, layer)

    do k = WS, WE
    do j = JS, JE
    do i = IS, IE
       Uij(i,j,k,2,3) = ( velz(i  ,j+1,k  )-velz(i  ,j  ,k  ) ) * RDY ! dw/dy, (x, v, interface)

    do k = KS, KE
    do j = JS, JE
    do i = IS-1, IE
       Uij(i,j,k,3,1) = ( velx(i  ,j  ,k+1)-velx(i  ,j  ,k  ) ) * RDZ ! du/dz, (u, y, interface)

    do k = KS, KE
    do j = JS-1, JE
    do i = IS, IE
       Uij(i,j,k,3,2) = ( vely(i  ,j  ,k+1)-vely(i  ,j  ,k  ) ) * RDZ ! dv/dz, (x, v, interface)

    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       Uij(i,j,k,3,3) = ( velz(i  ,j  ,k  )-velz(i  ,j  ,k-1) ) * RDZ ! dw/dz, (x, y, layer)

    ! ii == jj : at the center of control volume (cube)
    ! ii /= ii : on the edge   of control volume (cube)
    do ii = 1, 3
    do jj = 1, 3
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          Sij(i,j,k,(jj-1)*3+ii) = ( Uij(i,j,k,ii,jj) + Uij(i,j,k,jj,ii) ) * 0.5D0 ! 
       enddo
       enddo
       enddo
    enddo
    enddo

    ! fill IHALO & JHALO
    call COMM_vars( Sij(:,:,:,:) )


    ! calculate stratification effect (n2 and Pr are defined at w level)
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       a = 1.D0 / ( 1.D0 + dqsl(i,j,k) * Lv / Cp )

       c = (1+0.61*qw(i,j,k)-1.61*ql(i,j,k))* Lv/Cp /T_ov_pott(i,j,k) - 1.61 * pott(i,j,k)

       betatheta= 1 + 0.61 * qw(i,j,k) - 1.61*ql(i,j,k) - a*c*dqsl(i,j,k)* T_ov_pott(i,j,k)

       betaq= 0.61 * pott(i,j,k) + a * c

       n2(i,j,k) = grav / (sum(pottl(i,j,k:kz+1))/2) * &
                     (sum(betatheta * ( pottl(i,j,k+1)-pottl(i,j,k) )/dz + &
                     betaq * ( qw(i,j,k+1)-qw(i,j,k)  )/dz )/4)
    enddo
    enddo
    enddo

    ! calculate  Pr (Pr is defined at w level)
    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       Ri = GRAV / ( sum(pottl(i,j,k:kz+1) ) / 2 )  &
          * ( pottv(i,j,k+1)-pottv(i,j,k) ) &
          / (   ( velx(i,j,k+1)-velx(i,j,k) ) **2          &
              + ( vely(i,j,k+1)-vely(i,j,k) ) **2, 1d-10 ) &
                  * dz ))
       Ri = min( 0.25D0, max( 0.D0, Ri ) )
       Pr(i,j,k) = 0.33333d0 * ( 1.D0 - 4.D0*Ri ) +  4.D0*Ri
    enddo
    enddo
    enddo


    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       nu(i,j,k) = ( Cs * GRID_DX ) ** 2.D0   &
                 * ( sqrt( ( ( Sij(i,j,k,1,1) ) &
                           * ( Sij(i,j,k,1,1) ) &
                           + ( Sij(i,j,k,2,2) ) &
                           * ( Sij(i,j,k,2,2) ) &
                           + ( Sij(i,j,k,3,3) ) &
                           * ( Sij(i,j,k,3,3) ) &
                           ) * 2.D0             &
                         + ( ( Sij(i-1,j-1,k  ,1,2)+Sij(i-1,j  ,k  ,1,2)+Sij(i  ,j-1,k  ,1,2)+Sij(i,j,k,1,2) ) &
                           * ( Sij(i-1,j-1,k  ,1,2)+Sij(i-1,j  ,k  ,1,2)+Sij(i  ,j-1,k  ,1,2)+Sij(i,j,k,1,2) ) &
                           + ( Sij(i  ,j-1,k-1,2,3)+Sij(i  ,j-1,k  ,2,3)+Sij(i  ,j  ,k-1,2,3)+Sij(i,j,k,2,3) ) &
                           * ( Sij(i  ,j-1,k-1,2,3)+Sij(i  ,j-1,k  ,2,3)+Sij(i  ,j  ,k-1,2,3)+Sij(i,j,k,2,3) ) &
                           + ( Sij(i-1,j  ,k-1,3,1)+Sij(i  ,j  ,k-1,3,1)+Sij(i-1,j  ,k  ,3,1)+Sij(i,j,k,3,1) ) &
                           * ( Sij(i-1,j  ,k-1,3,1)+Sij(i  ,j  ,k-1,3,1)+Sij(i-1,j  ,k  ,3,1)+Sij(i,j,k,3,1) ) &
                           ) * 0.25D0                                                                          &
                         )                     &
                   - ( N2(i,j,k) / Pr(i,j,k) ) &
                   )
    enddo
    enddo
    enddo

    do k = WS, WE
    do j = JS, JE
    do i = IS, IE
       nuw(i,j,k) = ( Cs * GRID_DX ) ** 2.D0   &
    enddo
    enddo
    enddo

    nurho(:,:,:) = nu(:,:,:) * dens(:,:,:)


    do k = KS, KE
    do j = JS, JE
    do i = IS, IE
       FLXij(i,j,k,1,1) = 2.D0 * Sij(i,j,k,1,1) * nurho(i,j,k)
       FLXij(i,j,k,1,2) = 2.D0 * Sij(i,j,k,1,2) * (sum(nurho(ix,iy,kz-1:kz))/2)
       FLXij(i,j,k,1,3) = 2.D0 * Sij(i,j,k,1,3) * (sum(nurho(ix,iy-1:iy,kz))/2)

       FLXij(i,j,k,2,1) = FLXij(i,j,k,1,2) 
       FLXij(i,j,k,2,2) = 2.D0 * Sij(i,j,k,1,1) * nurho(i,j,k)
       FLXij(i,j,k,2,3) = FLXij(i,j,k,3,2) 

       FLXij(i,j,k,3,3) = 2.D0 * Sij(i,j,k,1,1) * nurho(i,j,k)
       FLXij(i,j,k,3,2) = gamma * 2.D0 * Sij(i,j,k,3,2) * (sum(nurho(ix-1:ix,iy,kz))/2)
       FLXij(i,j,k,3,1) = FLXij(i,j,k,1,3) 

       FLXpott(i,j,k,1) = ( pottl(i+1,j  ,k  ) - pottl(i,j,k) ) * RDX / Pr * ( sum( nurho(i,iy-1:jy,kz-1:kz))/4) 
       FLXpott(i,j,k,2) = ( pottl(i  ,j+1,k  ) - pottl(i,j,k) ) * RDY / Pr * ( sum( nurho(i-1:i,j,k-1:kz))/4) 
       FLXpott(i,j,k,3) = ( pottl(i  ,j  ,k+1) - pottl(i,j,k) ) * RDZ / Pr * ( sum( nurho(i-1:ix,iy-1:jy,kz))/4) 

       FLXq(i,j,k,1,:) = (q(ix+1,jy,kz,:)-q((i,j,k,:)) / dx / Pr *  ( sum(nurho(ix,iy-1:jy,kz-1:kz))/4)
       FLXq(i,j,k,2,:) = (q(ix,jy+1,kz,:)-q((i,j,k,:)) / dy / Pr *  ( sum(nurho(ix-1:i,j,k-1:kz))/4)
       FLXq(i,j,k,3,:) = (q(i,j,k+1,:)-q((i,j,k,:)) / dz / Pr *  ( sum(nurho(ix-1:ix,iy-1:jy,kz))/4)
!              FLXdens(i,j,k,1) = (dens(ix+1,jy,kz)-dens((i,j,k)) / dx / Pr *  ( sum(nurho(ix,iy-1:jy,kz-1:kz))/4)
!              FLXdens(i,j,k,2) = (dens(ix,jy+1,kz)-dens((i,j,k)) / dy / Pr *  ( sum(nurho(ix-1:i,j,k-1:kz))/4)
!              FLXdens(i,j,k,3) = (dens(i,j,k+1)-dens((i,j,k)) / dz / Pr *  ( sum(nurho(ix-1:ix,iy-1:jy,kz))/4)



          enddo
       enddo
    enddo


    ! tendency
    do k = KS, KE
    do j = JS, JE
    do i = IS-1, IE
       momx_t(i,j,k) = ( FLXij(i+1,j,k,1,1)-FLXij(i,j  ,k  ,1,1) ) * RDX &
                     + ( FLXij(i  ,j,k,2,1)-FLXij(i,j-1,k  ,2,1) ) * RDY &
                     + ( FLXij(i  ,j,k,3,1)-FLXij(i,j  ,k-1,3,1) ) * RDZ
    enddo
    enddo
    enddo


       momx_t(i,j,k) = ( FLXij(i+1,j,k,1,1)-FLXij(i,j   ,k  ,1,1) ) * RDX &
                     + ( FLXij(i  ,j,k,2,1)-FLXij(i,jy-1,k  ,2,1) ) * RDY &
                     + ( FLXij(i  ,j,k,3,1)-FLXij(i,j   ,k-1,3,1) ) * RDZ

             momy_t(i,j,k) =  &
                  (FLXij(i,j,k,1,1)-FLXij(ix-1,jy,kz,1,1)) /dx + &
                  (FLXij(ix,jy+1,kz,2,1)-FLXij(i,j,k,2,1)) /dy + &
                  (FLXij(i,j,k,3,1)-FLXij(i,j,k-1,3,1)) /dz 
             momz_t(i,j,k) =  &
                  (FLXij(i,j,k,1,1)-FLXij(ix-1,jy,kz,1,1)) /dx + &
                  (FLXij(i,j,k,2,1)-FLXij(ix,jy-1,kz,2,1)) /dy + &
                  (FLXij(i,j,k+1,3,1)-FLXij(i,j,k,3,1)) /dz 
!              dens_t(i,j,k) =  &
!                   (FLXdens(i,j,k,1)-FLXdens(ix-1,jy,kz,1)) /dx + &
!                   (FLXdens(i,j,k,2)-FLXdens(ix,jy-1,kz,2)) /dy + &
!                   (FLXdens(i,j,k,3)-FLXdens(i,j,k-1,3)) /dz 
             rhoq_t(i,j,k,:) =  &
                  (FLXq(i,j,k,1,:)-FLXq(ix-1,jy,kz,1,:)) /dx + &
                  (FLXq(i,j,k,2,:)-FLXq(ix,jy-1,kz,2,:)) /dy + &
                  (FLXq(i,j,k,3,:)-FLXq(i,j,k-1,3,:)) /dz 
!              rho_ene_t(i,j,k) =  &
!                   (FLXene(i,j,k,1)-FLXene(ix-1,jy,kz,1)) /dx + &
!                   (FLXene(i,j,k,2)-FLXene(ix,jy-1,kz,2)) /dy + &
!                   (FLXene(i,j,k,3)-FLXene(i,j,k-1,3)) /dz 
             rho_pott_t(i,j,k) =  &
                  (FLXpott(i,j,k,1)-FLXpott(ix-1,jy,kz,1)) /dx + &
                  (FLXpott(i,j,k,2)-FLXpott(ix,jy-1,kz,2)) /dy + &
                  (FLXpott(i,j,k,3)-FLXpott(i,j,k-1,3)) /dz 

    enddo
    enddo
    enddo


    ! add
    do kz = KS, KE
    do jy = JS+2, JE-2
    do ix = IS+2, IE-2
             momx(i,j,k) = momx(i,j,k) + dt * momx_t(i,j,k)
             momy(i,j,k) = momy(i,j,k) + dt * momy_t(i,j,k)
             momz(i,j,k) = momz(i,j,k) + dt * momz_t(i,j,k)
             pottl(i,j,k) = pottl(i,j,k) + dt * rho_pott_t(i,j,k) / dens(i,j,k) ! not accurate
!              pottl(i,j,k) = pottl(i,j,k) + dt * rho_ene_t(i,j,k) / dens(i,j,k) / Cp / T_ov_pott(i,j,k) ! hydrostatic app. is assumed
             rhoq(i,j,k,:) = rhoq(i,j,k,:) + dt * rhoq_t(i,j,k,:)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB

end module mod_atmos_phy_tb
