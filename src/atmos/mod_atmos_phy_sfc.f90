!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Surface fluxes
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
module mod_atmos_phy_sf
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
  public :: ATMOS_PHY_SF
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
  real(8) :: Z0M0 = 0.d0              !! base
  real(8) :: Z0MR = 0.018d0           !! rough factor
  real(8) :: Z0MS = 0.11d0            !! smooth factor
  real(8) :: Z0H0 = 1.4d-5            !! base
  real(8) :: Z0HR = 0.d0              !! rough factor
  real(8) :: Z0HS = 0.4d0             !! smooth factor
  real(8) :: Z0E0 = 1.3d-4            !! base
  real(8) :: Z0ER = 0.d0              !! rough factor
  real(8) :: Z0ES = 0.62d0            !! smooth factor
  real(8) :: VISAIR = 1.5d-5          !! kinematic viscosity 
  real(8) :: CM0 = 1.d-3              !! bulk coef for ustar
  real(8) :: USTRMN = 1.d-3           !! min(u*)
  real(8) :: Z0MMIN = 1.d-5           !! minimum
  real(8) :: Z0HMIN = 1.d-5           !! minimum
  real(8) :: Z0EMIN = 1.d-5           !! minimum
  real(8), public :: CNST_EGRAV    = 9.80616D0
  real(8), public :: CNST_KARMAN = 0.4D0
  real(8), public :: RV = 4.61D2
  real(8), public :: T0 = 2.7315D2
  real(8), private :: CMMIN = 1.d-5     ! min. bulk coef. of u
  real(8), private :: CHMIN = 1.d-5     ! min. bulk coef. of T
  real(8), private :: CEMIN = 1.d-5     ! min. bulk coef. of q
 ! real(8), private :: CMMAX = 1.d0      ! max. bulk coef. of u
  real(8), private :: CMMAX = 2.5d-3      ! max. bulk coef. of u
  real(8), private :: CHMAX = 1.d0      ! max. bulk coef. of T
  real(8), private :: CEMAX = 1.d0      ! max. bulk coef. of q
  real(8), private :: USMINM = 4.0d0    ! min. wind vel. for V
  real(8), private :: USMINH = 4.0d0    ! min. wind vel. for T
  real(8), private :: USMINE = 4.0d0    ! min. wind vel. for q
  real(8), private :: USMAXM = 1000.d0  ! max wind vel. for V
  real(8), private :: USMAXH = 1000.d0  ! max wind vel. for heat
  real(8), private :: USMAXE = 1000.d0  ! max wind vel. for q

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
  subroutine ATMOS_PHY_SF
!    use mod_stdio, only: &
!       IO_FID_LOG,  &
!       IO_L
    implicit none

    real(8), allocatable :: dens_s (:,:) ! density [kg kg-1]
    real(8), allocatable :: velx_us(:,:) ! wind in x-direction [m s-1]
    real(8), allocatable :: vely_vs(:,:) ! wind in y-direction [m s-1]
    real(8), allocatable :: velz_ws(:,:) ! wind in z-direction [m s-1]
    real(8), allocatable :: lwpt_s (:,:) ! liquid water potential temperature [K]
    real(8), allocatable :: qvap_s (:,:) ! water vapoer [kg kg-1]
    real(8), allocatable :: qclw_s (:,:) ! cloud water mixing ratio [kg kg-1]
    real(8), allocatable :: qraw_s (:,:) ! rain water mixing ratio [kg kg-1]
    real(8), allocatable :: pres_s (:,:) ! air pressure [Pa]

    real(8), allocatable :: SFLXm  (:,:,:)
    real(8), allocatable :: SFLXh  (:,:)
    real(8), allocatable :: SFLXe  (:,:,:)

    real(8), allocatable :: sst  (:,:) ! sea surface temperature [K]
    real(8), allocatable :: qvs  (:,:) ! saturation water vapor mixing ratio [kg kg-1]

! --- not used in SCALE experiments 111203 ---
    real(8), allocatable :: z1     (:,:) ! lowest model level [m]
    real(8), allocatable :: z_sfc  (:,:) ! surface height [m]
! --------------------------------------------

    real(8) :: z0m    ! roughness length [m]
    real(8) :: z0h    ! roughness length [m]
    real(8) :: z0e    ! roughness length [m]
    real(8) :: z_atm  ! height of atmosphere at the lowest model level [m]
    real(8) :: vabs   ! absolute velocity at the lowest model level [m s-1]
    real(8) :: vabs10 ! absolute velocity at 10 m height [m s-1]
    real(8) :: ustar  ! friction velocity [m s-1]
    real(8) :: es     ! water vapor pressure [Pa]
    real(8) :: cmx    ! 
    real(8) :: cmm    ! surface exchange coefficient for momentum 
    real(8) :: cmh    ! surface exchange coefficient for temperature
    real(8) :: cme    ! surface exchange coefficient for moisture variables
    real(8) :: R10m   ! 
    real(8) :: R10h   ! 
    real(8) :: R10e   ! 

    integer :: i, j
    !---------------------------------------------------------------------------

    allocate( dens_s  (IS:IE,JS:JE) )
    allocate( velx_us (IS-1:IE,JS:JE) )
    allocate( vely_vs (IS:IE,JS-1:JE) )
    allocate( velz_ws (IS:IE,JS:JE) )
    allocate( lwpt_s  (IS:IE,JS:JE) )
    allocate( qvap_s  (IS:IE,JS:JE) )
    allocate( qclw_s  (IS:IE,JS:JE) )
    allocate( qraw_s  (IS:IE,JS:JE) )
    allocate( pres_s  (IS:IE,JS:JE) )
    allocate( z1      (IS:IE,JS:JE) )
    allocate( z_sfc   (IS:IE,JS:JE) )
    allocate( sst     (IS:IE,JS:JE) )
    allocate( qvs     (IS:IE,JS:JE) )
    allocate( SFLXm (IS:IE,JS:JE,3) )
    allocate( SFLXh (IS:IE,JS:JE)   )
    allocate( SFLXe (IS:IE,JS:JE,3) )

! --- not used in ocean experiments ---
    z_sfc (:,:) = 0.d0
    z1    (:,:) = dz
    R10h = 1
    R10e = 1

    do j = JS, JE
    do i = IS, IE

       z_atm = z1(i,j)-z_srf(i,j) ! z_srf = 0 : surface height in NICAM
! --- absolute velocity ---
       vabs = ((velx_us(i-1,j)+velx_us(i,j))/2)**2 &
            + ((vely_vs(i,j-1)+vely_vs(i,j))/2)**2 &
            + velz_ws(i,j)**2
! --- velocity at the 10-m height ---
       R10m = 10/z_atm
       vabs10 = vabs * R10m
! --- friction velocity ---
       ustar = max ( sqrt ( CM0 * vabs ), USTRMN )

! --- roughness lengths ---
       z0m = max( Z0M0 + Z0MR * ustar ** 2 / CNST_EGRAV + Z0MS * VISAIR / ustar , Z0MMIN )
       z0h = max( Z0H0 + Z0HR * ustar ** 2 / CNST_EGRAV + Z0HS * VISAIR / ustar , Z0HMIN )
       z0e = max( Z0E0 + Z0ER * ustar ** 2 / CNST_EGRAV + Z0ES * VISAIR / ustar , Z0EMIN )

! --- surface exchange coefficients ---
       cmx  = CNST_KARMAN **2 / log( z_atm / z0m )
       cmm  = max( min( cmx / log( z_atm / z0m ) , CMMAX ) , CMMAX )
       cmh  = max( min( cmx / log( z_atm / z0h ) , CHMAX ) , CHMAX )
       cme  = max( min( cmx / log( z_atm / z0e ) , CEMAX ) , CEMIN )

       cmm = cmm * min( max( vabs , USMINM ) , USMAXM )
       cmh = cmh * min( max( vabs , USMINH ) , USMAXH )
       cme = cme * min( max( vabs , USMINE ) , USMAXE )

       es  = 6.11d2 * dexp( (LV/RV) * ( (sst(i,j)-T0)/(T0*sst(i,j)) ) )
       qvs(i,j) = 6.22d-1 * es / (pres_s(i,j)-es)

! --- surface fluxes ---
       SFLXm(i,j,1) = dens_s(i,j) * vabs10 * cmm * ( velx_us(i-1,j) + velx_us(i,j) )/2 * 10/z_atm 
       SFLXm(i,j,2) = dens_s(i,j) * vabs10 * cmm * ( vely_vs(i,j-1) + vely_vs(i,j) )/2 * 10/z_atm
       SFLXm(i,j,3) = dens_s(i,j) * vabs10 * cmm * velz_ws(i,j) * 10/z_atm
       SFLXh(i,j)   = dens_s(i,j) * vabs10 * cmh * ( sst(i,j) - lwpt_s(i,j) * R10h )
       SFLXe(i,j,1) = dens_s(i,j) * vabs10 * cme * ( qvs(i,j) - qvap_s(i,j) * R10e )
       SFLXe(i,j,2) = 0 !dens_s(i,j) * vabs10 * cme * (qclw_s(i,j) * R10e )
       SFLXe(i,j,3) = 0 !dens_s(i,j) * vabs10 * cme * (qraw_s(i,j) * R10e )
       
    end do
    end do

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sfc
