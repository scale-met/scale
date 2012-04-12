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
!! @li      2011-11-29 (S.Iga)      [new]
!! @li      2011-12-11 (H.Yashiro)  [mod] integrate to SCALE3
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2012-03-27 (H.Yashiro)  [mod] reconstruction
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb
  !-----------------------------------------------------------------------------
  !
  !++ used modules
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
  public :: ATMOS_PHY_TB_setup
  public :: ATMOS_PHY_TB

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  real(8), private,      save :: TB_length(KA,IA,JA) ! mixing length

  real(8), private, parameter :: Cs        = 0.18D0 ! (Sullivan et al.1994, Nakanishi and Niino)
  real(8), private, parameter :: DUDZ2_min = 1.D-6  ! minimum limit of (dU/dz)**2
  real(8), private, parameter :: FACTC     = 2.D0 / 3.D0

  integer, private, parameter :: ZDIR = 1
  integer, private, parameter :: XDIR = 2
  integer, private, parameter :: YDIR = 3

  !-----------------------------------------------------------------------------
contains

  subroutine ATMOS_PHY_TB_setup
    use mod_grid, only : &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY
    implicit none

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-TB]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Smagorinsky-type Eddy Viscocity Model'

    do j = JS, JE+1
    do i = IS, IE+1
    do k = KS, KE
       TB_length(k,i,j) = ( CDZ(k) * CDX(i) * CDY(j) )**(1.D0/3.D0)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_TB_setup

  !-----------------------------------------------------------------------------
  !> Smagorinsky-type turblence
  !>
  !> comment:
  !>  1, Pr is given linearly (iga)
  !>  2, dx,dy,dz is assumed to be constant.
  !>  3, gamma is assumed to be 1 now (i.e. dx=dy=dz).
  !>  4, heat flux is not accurate yet. (i.e. energy is not conserved, see *1)
  !>  5, stratification effect is considered.
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB
    use mod_const, only : &
       GRAV => CONST_GRAV
    use mod_time, only: &
       dttb => TIME_DTSEC_ATMOS_PHY_TB
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_grid, only : &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       ATMOS_vars_fillhalo, &
       ATMOS_vars_total,    &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    implicit none

    ! tendency
    real(8) :: MOMZ_t(KA,IA,JA)
    real(8) :: MOMX_t(KA,IA,JA)
    real(8) :: MOMY_t(KA,IA,JA)
    real(8) :: RHOT_t(KA,IA,JA)
    real(8) :: QTRC_t(KA,IA,JA,QA)

    ! diagnostic variables
    real(8) :: VELZ(KA,IA,JA)
    real(8) :: VELX(KA,IA,JA)
    real(8) :: VELY(KA,IA,JA)
    real(8) :: POTT(KA,IA,JA)

    ! deformation rate tensor
    real(8) :: Sij_zz(KA,IA,JA)
    real(8) :: Sij_yy(KA,IA,JA)
    real(8) :: Sij_xx(KA,IA,JA)
    real(8) :: Sij_zx(KA,IA,JA)
    real(8) :: Sij_zy(KA,IA,JA)
    real(8) :: Sij_xy(KA,IA,JA)

    real(8) :: nu(KA,IA,JA)       ! eddy viscosity
    real(8) :: nuc(KA,IA,JA)      !
    real(8) :: tke(KA,IA,JA)      ! turburent kinetic energy

    real(8) :: Pr   (KA,IA,JA)    ! Prantle number
    real(8) :: Ri   (KA,IA,JA)    ! Richardson number
    real(8) :: buoy (KA,IA,JA)    ! Buoyancy term
    real(8) :: DUDZ2(KA)          ! (dU/dz)**2 U=sqrt(u**2+v**2)
    real(8) :: temp (KA)

    real(8) :: qflx_sgs(KA,IA,JA,3)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

!    VELZ(:,:,:) = -9.999D30
!    VELX(:,:,:) = -9.999D30
!    VELY(:,:,:) = -9.999D30
!    POTT(:,:,:) = -9.999D30
!
!    Sij_zz(:,:,:) = -9.999D30
!    Sij_yy(:,:,:) = -9.999D30
!    Sij_xx(:,:,:) = -9.999D30
!    Sij_zx(:,:,:) = -9.999D30
!    Sij_zy(:,:,:) = -9.999D30
!    Sij_xy(:,:,:) = -9.999D30
!
!    nu  (:,:,:) = -9.999D30
!    nuc (:,:,:) = -9.999D30
!    tke (:,:,:) = -9.999D30
!    Pr  (:,:,:) = -9.999D30
!    Ri  (:,:,:) = -9.999D30
!    buoy(:,:,:) = -9.999D30
!
!    qflx_sgs(:,:,:,:) = -9.999D30

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: SGS Parameterization'

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! momentum -> velocity
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE-1
          VELZ(k,i,j) = 2.D0 * MOMZ(k,i,j) / ( DENS(k+1,i,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
          VELZ(KS-1,i,j) = 0.D0
          VELZ(KE  ,i,j) = 0.D0
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          VELX(k,i,j) = 2.D0 * MOMX(k,i,j) / ( DENS(k,i+1,j)+DENS(k,i,j) )
       enddo
       enddo
       enddo

       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          VELY(k,i,j) = 2.D0 * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
       enddo
       enddo
       enddo

       ! potential temperature
       do j = JJS-1, JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) 
       enddo
       enddo
       enddo

       ! Stress Tensor
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
          Sij_zz(k,i,j) = ( VELZ(k,i,j)-VELZ(k-1,i,j) ) * RCDZ(k) ! dw/dz, (x, y, layer)
       enddo
       enddo
       enddo

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
          Sij_xx(k,i,j) = ( VELX(k,i,j)-VELX(k,i-1,j) ) * RCDX(i) ! du/dx, (x, y, layer)
       enddo
       enddo
       enddo

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
          Sij_yy(k,i,j) = ( VELY(k,i,j)-VELY(k,i,j-1) ) * RCDY(j) ! dv/dy, (x, y, layer)
       enddo
       enddo
       enddo

       do j = JJS,   JJE+1
       do i = IIS-1, IIE
       do k = KS, KE-1
          Sij_zx(k,i,j) = 0.5D0 * ( ( VELX(k+1,i,j)-VELX(k,i,j) ) * RCDZ(k) & ! du/dz, (u, y, interface)
                                  + ( VELZ(k,i+1,j)-VELZ(k,i,j) ) * RCDX(i) ) ! dw/dx, (u, y, interface)
       enddo
       enddo
       enddo

       do j = JJS-1, JJE
       do i = IIS,   IIE+1
       do k = KS, KE-1
          Sij_zy(k,i,j) = 0.5D0 * ( ( VELY(k+1,i,j)-VELY(k,i,j) ) * RCDZ(k) & ! dv/dz, (x, v, interface)
                                  + ( VELZ(k,i,j+1)-VELZ(k,i,j) ) * RCDY(j) ) ! dw/dy, (x, v, interface)
       enddo
       enddo
       enddo

       do j = JJS-1, JJE
       do i = IIS-1, IIE
       do k = KS, KE
          Sij_xy(k,i,j) = 0.5D0 * ( ( VELY(k,i+1,j)-VELY(k,i,j) ) * RCDX(i) & ! dv/dx, (u, v, layer)
                                  + ( VELX(k,i,j+1)-VELX(k,i,j) ) * RCDY(j) ) ! du/dy, (u, v, layer)
       enddo
       enddo
       enddo


       ! calculate Ri, Pr at x, y, interface
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          do k = KS, KE-1
             DUDZ2(k) = ( 0.5D0 * ( VELX(k+1,i,j)+VELX(k+1,i-1,j)-VELX(k,i,j)-VELX(k,i-1,j) ) * RCDZ(k) )**2 &
                      + ( 0.5D0 * ( VELY(k+1,i,j)+VELY(k+1,i,j-1)-VELY(k,i,j)-VELY(k,i,j-1) ) * RCDZ(k) )**2
          enddo

          do k = KS, KE-1
             DUDZ2(k) = max( DUDZ2(k), DUDZ2_min )
          enddo

          ! Ri = g / rho * drho/dz / dU/dz
          do k = KS, KE-1
             Ri(k,i,j) = GRAV * 2.D0 / ( DENS(k+1,i,j)+DENS(k,i,j) ) * ( DENS(k,i,j)-DENS(k+1,i,j) ) * RFDZ(k) &
                       / DUDZ2(k)
          enddo

          do k = KS, KE-1
             Ri(k,i,j) = min( 0.25D0, max( 0.D0, Ri(k,i,j) ) )
          enddo

          Ri(KS-1,i,j) = Ri(KS  ,i,j)
          Ri(KE  ,i,j) = Ri(KE-1,i,j)

          do k = KS, KE-1
             Pr(k,i,j) = ( 1.D0-4.D0*Ri(k,i,j) ) / 3.D0 + 4.D0*Ri(k,i,j)
          enddo

          Pr( 1:KS-1,i,j) = Pr(KS  ,i,j)
          Pr(KE:KA  ,i,j) = Pr(KE-1,i,j)

          ! buoyancy
          do k = KS, KE-1
             buoy(k,i,j) = GRAV * 2.D0 / ( POTT(k+1,i,j)+POTT(k,i,j) ) * ( POTT(k,i,j)-POTT(k+1,i,j) ) * RFDZ(k) &
                         / Pr(k,i,j)
          enddo
       enddo
       enddo


       ! Turbulent Viscosity
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS, KE-1
             temp(k) = 2.D0 * ( ( ( Sij_zz(k  ,i  ,j  ) + Sij_zz(k+1,i+1,j+1) &
                                  + Sij_zz(k+1,i  ,j  ) + Sij_zz(k+1,i+1,j  ) &
                                  + Sij_zz(k  ,i+1,j  ) + Sij_zz(k  ,i+1,j+1) &
                                  + Sij_zz(k  ,i  ,j+1) + Sij_zz(k+1,i  ,j+1) ) * 0.125D0 )**2 &
                              + ( ( Sij_xx(k  ,i  ,j  ) + Sij_xx(k+1,i+1,j+1) &
                                  + Sij_xx(k+1,i  ,j  ) + Sij_xx(k+1,i+1,j  ) &
                                  + Sij_xx(k  ,i+1,j  ) + Sij_xx(k  ,i+1,j+1) &
                                  + Sij_xx(k  ,i  ,j+1) + Sij_xx(k+1,i  ,j+1) ) * 0.125D0 )**2 &
                              + ( ( Sij_yy(k  ,i  ,j  ) + Sij_yy(k+1,i+1,j+1) &
                                  + Sij_yy(k+1,i  ,j  ) + Sij_yy(k+1,i+1,j  ) &
                                  + Sij_yy(k  ,i+1,j  ) + Sij_yy(k  ,i+1,j+1) &
                                  + Sij_yy(k  ,i  ,j+1) + Sij_yy(k+1,i  ,j+1) ) * 0.125D0 )**2 ) &
                     + 4.D0 * ( ( ( Sij_zx(k  ,i  ,j  ) + Sij_zx(k  ,i  ,j+1) ) * 0.5D0   )**2 &
                              + ( ( Sij_zy(k  ,i  ,j  ) + Sij_zy(k  ,i+1,j  ) ) * 0.5D0   )**2 &
                              + ( ( Sij_xy(k  ,i  ,j  ) + Sij_xy(k+1,i  ,j  ) ) * 0.5D0   )**2 ) &
                     - 0.25D0 * ( buoy(k,i,j) + buoy(k,i+1,j) + buoy(k,i,j+1) + buoy(k,i+1,j+1) )
          enddo

          do k = KS, KE-1
             temp(k) = max( 0.D0, temp(k) )
          enddo

          do k = KS, KE-1
             nu(k,i,j) = ( Cs * TB_length(k,i,j) )**2 * sqrt( temp(k) )
          enddo
       enddo
       enddo
       !OCL XFILL
       do j = JJS, JJE
       do i = IIS, IIE
          nu( 1:KS-1,i,j) = 0.D0
          nu(KE:KA  ,i,j) = 0.D0
       enddo
       enddo

    enddo
    enddo

    call COMM_vars8( nu(:,:,:), 1 )
    call COMM_vars8( Pr(:,:,:), 2 )
    call COMM_wait ( nu(:,:,:), 1 )
    call COMM_wait ( Pr(:,:,:), 2 )

    !##### Start Upadate #####

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! nu at (x, y, layer)
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE
          nuc(k,i,j) = 0.125D0 * ( nu(k  ,i  ,j  ) + nu(k-1,i-1,j-1) &
                                 + nu(k-1,i  ,j  ) + nu(k-1,i-1,j  ) &
                                 + nu(k  ,i-1,j  ) + nu(k  ,i-1,j-1) &
                                 + nu(k  ,i  ,j-1) + nu(k-1,i  ,j-1) )
          tke(k,i,j) = nuc(k,i,j)*nuc(k,i,j) / ( Cs * TB_length(k,i,j) )**2
       enddo
       enddo
       enddo

       !##### momentum equation (z) #####
       ! at (x, y, layer)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          qflx_sgs(k,i,j,ZDIR) = DENS(k,i,j) &
                               * ( -2.D0 * nuc(k,i,j) * Sij_zz(k,i,j) &
                                 + FACTC * tke(k,i,j)                 )
       enddo
       enddo
       enddo
       ! at (u, y, interface)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_sgs(k,i,j,XDIR) = 0.25D0 * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) &
                               * ( -( nu(k,i,j)+nu(k,i,j-1) ) * Sij_zx(k,i,j) )
       enddo
       enddo
       enddo
       ! at (x, v, interface)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_sgs(k,i,j,YDIR) = 0.25D0 * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) &
                               * ( -( nu(k,i,j)+nu(k,i-1,j) ) * Sij_zy(k,i,j) )
       enddo
       enddo
       enddo

       !--- tendency momentum(z)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          MOMZ_t(k,i,j) = - ( ( qflx_sgs(k+1,i,j,ZDIR)-qflx_sgs(k,i  ,j  ,ZDIR) ) * RFDZ(k) &
                            + ( qflx_sgs(k  ,i,j,XDIR)-qflx_sgs(k,i-1,j  ,XDIR) ) * RCDX(i) &
                            + ( qflx_sgs(k  ,i,j,YDIR)-qflx_sgs(k,i  ,j-1,YDIR) ) * RCDY(j) )
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          MOMZ_t(KE,i,j) = 0.D0
       enddo
       enddo

       !##### momentum equation (x) #####
       ! at (u, y, interface)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs(k,i,j,ZDIR) = 0.25D0 * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) &
                               * ( -( nu(k,i,j)+nu(k,i,j-1) ) * Sij_zx(k,i,j) )
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs(KS-1,i,j,ZDIR) = SFLX_MOMX(i,j) ! bottom boundary
          qflx_sgs(KE  ,i,j,ZDIR) = 0.D0 ! top boundary
       enddo
       enddo
       ! at (x, y, layer)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
          qflx_sgs(k,i,j,XDIR) = DENS(k,i,j) &
                               * ( -2.D0 * nuc(k,i,j) * Sij_xx(k,i,j) &
                                 + FACTC * tke(k,i,j)                 )
       enddo
       enddo
       enddo
       ! at (u, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_sgs(k,i,j,YDIR) = 0.25D0 * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
                               * ( -( nu(k,i,j)+nu(k-1,i,j) ) * Sij_xy(k,i,j) )
       enddo
       enddo
       enddo

       !--- tendency momentum(x)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMX_t(k,i,j) = - ( ( qflx_sgs(k,i  ,j,ZDIR)-qflx_sgs(k-1,i,j,  ZDIR) ) * RCDZ(k) &
                            + ( qflx_sgs(k,i+1,j,XDIR)-qflx_sgs(k  ,i,j,  XDIR) ) * RFDX(i) &
                            + ( qflx_sgs(k,i  ,j,YDIR)-qflx_sgs(k  ,i,j-1,YDIR) ) * RCDY(j) )
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####
       ! at (x, v, interface)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs(k,i,j,ZDIR) = 0.25D0 * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) &
                               * ( -( nu(k,i,j)+nu(k,i-1,j) ) * Sij_zy(k,i,j) )
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs(KS-1,i,j,ZDIR) = SFLX_MOMY(i,j) ! bottom boundary
          qflx_sgs(KE  ,i,j,ZDIR) = 0.D0 ! top boundary
       enddo
       enddo

       ! at (u, v, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_sgs(k,i,j,XDIR) = 0.25D0 * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
                               * ( -( nu(k,i,j)+nu(k-1,i,j) ) * Sij_xy(k,i,j) )
       enddo
       enddo
       enddo

       ! at (x, y, layer)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
          qflx_sgs(k,i,j,YDIR) = DENS(k,i,j) &
                               * ( -2.D0 * nuc(k,i,j) * Sij_yy(k,i,j) &
                                 + FACTC * tke(k,i,j)                 )
       enddo
       enddo
       enddo

       !--- tendency momentum(y)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMY_t(k,i,j) = - ( ( qflx_sgs(k,i,j  ,ZDIR)-qflx_sgs(k-1,i  ,j,ZDIR) ) * RCDZ(k) &
                            + ( qflx_sgs(k,i,j  ,XDIR)-qflx_sgs(k  ,i-1,j,XDIR) ) * RCDX(i) &
                            + ( qflx_sgs(k,i,j+1,YDIR)-qflx_sgs(k  ,i  ,j,YDIR) ) * RFDY(j) )
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at (x, y, interface)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs(k,i,j,ZDIR) = 0.125D0 * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                               * ( - ( nu(k,i,j) + nu(k,i-1,j) + nu(k,i,j-1) + nu(k,i-1,j-1) ) &
                                   / ( Pr(k,i,j) )                                             &
                                   * ( POTT(k+1,i,j)-POTT(k,i,j) ) * RFDZ(k)                   )
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs(KS-1,i,j,ZDIR) = SFLX_POTT(i,j)
          qflx_sgs(KE  ,i,j,ZDIR) = 0.0D0
       enddo
       enddo

       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_sgs(k,i,j,XDIR) = 0.5D0 * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * ( - ( nu(k,i,j) + nu(k,i,j-1) + nu(k-1,i,j) + nu(k-1,i,j-1) ) &
                                   / ( Pr(k,i,j) + Pr(k,i+1,j) + Pr(k-1,i,j) + Pr(k-1,i+1,j) ) &
                                   * ( POTT(k,i+1,j)-POTT(k,i,j) ) * RFDX(i)                   )
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_sgs(k,i,j,YDIR) = 0.5D0 * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * ( - ( nu(k,i,j) + nu(k,i-1,j) + nu(k-1,i,j) + nu(k-1,i-1,j) ) &
                                   / ( Pr(k,i,j) + Pr(k,i,j+1) + Pr(k-1,i,j) + Pr(k-1,i,j+1) ) &
                                   * ( POTT(k,i,j+1)-POTT(k,i,j) ) * RFDY(j)                 )
       enddo
       enddo
       enddo

       !--- tendency rho*theta
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          RHOT_t(k,i,j) = - ( ( qflx_sgs(k,i,j,ZDIR)-qflx_sgs(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                            + ( qflx_sgs(k,i,j,XDIR)-qflx_sgs(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                            + ( qflx_sgs(k,i,j,YDIR)-qflx_sgs(k  ,i,  j-1,YDIR) ) * RCDY(j) )
       enddo
       enddo
       enddo

    enddo
    enddo

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          MOMZ(k,i,j) = MOMZ(k,i,j) + dttb * MOMZ_t(k,i,j)
          MOMX(k,i,j) = MOMX(k,i,j) + dttb * MOMX_t(k,i,j)
          MOMY(k,i,j) = MOMY(k,i,j) + dttb * MOMY_t(k,i,j)
          RHOT(k,i,j) = RHOT(k,i,j) + dttb * RHOT_t(k,i,j)
       enddo
       enddo
       enddo

    enddo
    enddo

    !##### Tracers #####
    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! at (x, y, interface)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs(k,i,j,ZDIR) = 0.125D0 * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                               * ( - ( nu(k,i,j) + nu(k,i-1,j) + nu(k,i,j-1) + nu(k,i-1,j-1) ) &
                                   / ( Pr(k,i,j) )                                             &
                                   * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) * RFDZ(k)             )
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs(KS-1,i,j,ZDIR) = 0.0D0
          qflx_sgs(KE  ,i,j,ZDIR) = 0.0D0
       enddo
       enddo

       ! Surface QV Flux
       if ( iq == I_QV ) then
         do j = JJS, JJE
         do i = IIS, IIE
             qflx_sgs(KS-1,i,j,ZDIR) = SFLX_QV(i,j)
          enddo
          enddo
       endif

       ! at (u, y, layer)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_sgs(k,i,j,XDIR) = 0.5D0 * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * ( - ( nu(k,i,j) + nu(k,i,j-1) + nu(k-1,i,j) + nu(k-1,i,j-1) ) &
                                   / ( Pr(k,i,j) + Pr(k,i+1,j) + Pr(k-1,i,j) + Pr(k-1,i+1,j) ) &
                                   * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) * RFDX(i)             )
       enddo
       enddo
       enddo
       ! at (x, v, layer)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_sgs(k,i,j,YDIR) = 0.5D0 * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * ( - ( nu(k,i,j) + nu(k,i-1,j) + nu(k-1,i,j) + nu(k-1,i-1,j) ) &
                                   / ( Pr(k,i,j) + Pr(k,i,j+1) + Pr(k-1,i,j) + Pr(k-1,i,j+1) ) &
                                   * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) * RFDY(j)             )
       enddo
       enddo
       enddo

       !--- tendency tracers
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          QTRC_t(k,i,j,iq) = - ( ( qflx_sgs(k,i,j,ZDIR)-qflx_sgs(k-1,i,  j,  ZDIR) ) * RCDZ(k) &
                               + ( qflx_sgs(k,i,j,XDIR)-qflx_sgs(k  ,i-1,j,  XDIR) ) * RCDX(i) &
                               + ( qflx_sgs(k,i,j,YDIR)-qflx_sgs(k  ,i,  j-1,YDIR) ) * RCDY(j) ) &
                             / DENS(k,i,j)
       enddo
       enddo
       enddo

    enddo
    enddo


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + dttb * QTRC_t(k,i,j,iq)
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE
          if ( QTRC(k,i,j,iq) < 1.D-10 ) then
             QTRC(k,i,j,iq) = 0.D0
          endif
       enddo
       enddo
       enddo

    enddo
    enddo

    enddo ! scalar quantities loop

    call ATMOS_vars_fillhalo

    call ATMOS_vars_total

    call HIST_in( MOMZ_t(:,:,:), 'MOMZ_t_tb', 'tendency of MOMZ in tb', 'kg/m2/s2',  '3D', dttb )
    call HIST_in( MOMX_t(:,:,:), 'MOMX_t_tb', 'tendency of MOMX in tb', 'kg/m2/s2',  '3D', dttb )
    call HIST_in( MOMY_t(:,:,:), 'MOMY_t_tb', 'tendency of MOMY in tb', 'kg/m2/s2',  '3D', dttb )
    call HIST_in( RHOT_t(:,:,:), 'RHOT_t_tb', 'tendency of RHOT in tb', 'K*kg/m3/s', '3D', dttb )
    do iq = 1, QA
       call HIST_in( QTRC_t(:,:,:,iq), AQ_NAME(iq)//'_t_tb', AQ_DESC(iq), AQ_UNIT(iq)//'/s', '3D', dttb )
    enddo

    call HIST_in( tke (:,:,:), 'TKE',  'turburent kinetic energy', 'J/m3', '3D', dttb )
    call HIST_in( nuc (:,:,:), 'NU',   'eddy viscosity',           'm2/s', '3D', dttb )
    call HIST_in( Pr  (:,:,:), 'Pr',   'Prantle number',           'NIL',  '3D', dttb )
    call HIST_in( Ri  (:,:,:), 'Ri',   'Richardson number',        'NIL',  '3D', dttb )

    return
  end subroutine ATMOS_PHY_TB

end module mod_atmos_phy_tb
