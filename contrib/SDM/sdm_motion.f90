!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Module to control the motion of super-droplets (advection, sedmentation, precipitation)
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-11 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-11 (S.Shima) [rev] Fixed the bug of interpolation in sdm_getvel
!! @li      2015-06-25 (S.Shima) [add] fapp_start/stop added for performance monitoring
!! @li      2015-06-26 (S.Shima) [add] OCL added for Auto parallelization on K/FX10
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_motion
  use scale_precision
  
  implicit none
  private
  public :: sdm_getvz, sdm_getvel, sdm_move

contains
  subroutine sdm_move(sdm_dtadv,                      &
       sd_num,sd_u,sd_v,sd_vz,sd_x,sd_y,sd_rk)
    use scale_grid, only: &
         DZ 
    use m_sdm_common, only: &
         VALID2INVALID
    ! Input variables
    real(RP), intent(in) :: sdm_dtadv ! time step of {motion of super-droplets} process
    integer, intent(in) :: sd_num ! number of super-droplets
    real(RP), intent(in) :: sd_u(1:sd_num) ! x-direction velocity of super-droplets
    real(RP), intent(in) :: sd_v(1:sd_num) ! y-direction velocity of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! z velocity of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Work variables
    integer :: n           ! index
    real(RP) :: dz_inv
    !---------------------------------------------------------------------
    ! The advection process of Super-Droplets.
    ! now only support uniform grid

    ! Section specification for fapp profiler
    call fapp_start("sdm_move",1,1)

    dz_inv=1.0_RP/DZ
    do n=1,sd_num

       !### skip invalid super-droplets ###!
       
       if( sd_rk(n)<VALID2INVALID ) cycle
       
       !### move super-droplets (explicit) ###!
       sd_x(n)  = sd_x(n)  + sd_u(n)  * real(sdm_dtadv,kind=RP)
       sd_y(n)  = sd_y(n)  + sd_v(n)  * real(sdm_dtadv,kind=RP)
       sd_rk(n) = sd_rk(n) + sd_vz(n) * real(sdm_dtadv,kind=RP) * dz_inv
    enddo
    
    ! Section specification for fapp profiler
    call fapp_stop("sdm_move",1,1)

    return
  end subroutine sdm_move
  !----------------------------------------------------------------------------
  subroutine sdm_getvel(u_scale,v_scale,w_scale,              &
       sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vzw)
    use scale_grid, only: &
         FX => GRID_FX, &
         FY => GRID_FY
    use scale_grid_index, only: &
         IA,JA,KA
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID
      
    ! Input variables
    real(RP), intent(in) :: u_scale(KA,IA,JA) ! x-component velocity
    real(RP), intent(in) :: v_scale(KA,IA,JA) ! y-component velocity
    real(RP), intent(in) :: w_scale(KA,IA,JA) ! z-component velocity
    integer, intent(in) :: sd_num  ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! face index[k/real] of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_vzw(1:sd_num)  ! terminal velocity [in] / z velocity [out] of super-droplets
    ! Output variables
    real(RP), intent(out) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
    real(RP), intent(out) :: sd_u(1:sd_num) ! x-direction velocity of super-droplets
    real(RP), intent(out) :: sd_v(1:sd_num) ! y-direction velocity of super-droplets
    ! Work variables
    real(RP) :: sd_vz ! terminal velocity of super-droplets
    real(RP) :: sd_w  ! z velocity of super-droplets
    real(RP) :: ri    ! real index [i] of super-droplets
    real(RP) :: rj    ! real index [j] of super-droplets
    real(RP) :: rk    ! real index [k] of super-droplets
    real(RP) :: sXm   ! variable for inteporation
    real(RP) :: sXp   ! variable for inteporation
    real(RP) :: sYm   ! variable for inteporation
    real(RP) :: sYp   ! variable for inteporation
    real(RP) :: sZm   ! variable for inteporation
    real(RP) :: sZp   ! variable for inteporation
    integer :: iXm    ! index for inteporation
    integer :: iXp    ! index for inteporation
    integer :: iYm    ! index for inteporation
    integer :: iYp    ! index for inteporation
    integer :: iZm    ! index for inteporation
    integer :: iZp    ! index for inteporation
    integer :: n, i, j, k      ! index
    !-------------------------------------------------------------------
    ! The advection process of Super-Droplets.

    ! Section specification for fapp profiler
    call fapp_start("sdm_getvel",1,1)

    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    do n=1,sd_num

       !### skip invalid super-droplets ###!

       if( sd_rk(n)<VALID2INVALID ) cycle

       !### get the real index of super-droplets ###!

       ri = sd_ri(n)
       rj = sd_rj(n)
       rk = sd_rk(n)
       
       !### x-velocity at the position of super-droplets ###!
       !! interpolation in face grid
       iXm = floor(ri)
       iXp = iXm + 1
       sXm = ri - real(iXm,kind=RP)
       sXp = 1.0_RP - sXm

       iYm = floor(rj+0.5_RP)
       iYp = iYm + 1
       sYm = (rj+0.5_RP) - real(iYm,kind=RP)
       sYp = 1.d0 - sYm

       iZm = floor(rk+0.5_RP)
       iZp = iZm + 1
       sZm = (rk+0.5_RP) - real(iZm,kind=RP)
       sZp = 1.0_RP - sZm

       sd_u(n) = u_scale(iZm,iXm,iYm) * ( SXp * SYp * SZp )             &
            + u_scale(iZm,iXp,iYm) * ( SXm * SYp * SZp )             &
            + u_scale(iZm,iXm,iYp) * ( SXp * SYm * SZp )             &
            + u_scale(iZm,iXp,iYp) * ( SXm * SYm * SZp )             &
            + u_scale(iZp,iXm,iYm) * ( SXp * SYp * SZm )             &
            + u_scale(iZp,iXp,iYm) * ( SXm * SYp * SZm )             &
            + u_scale(iZp,iXm,iYp) * ( SXp * SYm * SZm )             &
            + u_scale(iZp,iXp,iYp) * ( SXm * SYm * SZm )

       !### y-velocity at the position of super-droplets ###!
       !! interpolation in face grid
       iXm = floor(ri+0.5_RP)
       iXp = iXm + 1
       sXm = (ri+0.5_RP) - real(iXm,kind=RP)
       sXp = 1.0_RP - sXm

       iYm = floor(rj)
       iYp = iYm + 1
       sYm = rj - real(iYm,kind=RP)
       sYp = 1.0_RP - sYm

       iZm = floor(rk+0.5_RP)
       iZp = iZm + 1
       sZm = (rk+0.5_RP) - real(iZm,kind=RP)
       sZp = 1.0_RP - sZm

       sd_v(n) = v_scale(iZm,iXm,iYm) * ( SXp * SYp * SZp )             &
            + v_scale(iZm,iXp,iYm) * ( SXm * SYp * SZp )             &
            + v_scale(iZm,iXm,iYp) * ( SXp * SYm * SZp )             &
            + v_scale(iZm,iXp,iYp) * ( SXm * SYm * SZp )             &
            + v_scale(iZp,iXm,iYm) * ( SXp * SYp * SZm )             &
            + v_scale(iZp,iXp,iYm) * ( SXm * SYp * SZm )             &
            + v_scale(iZp,iXm,iYp) * ( SXp * SYm * SZm )             &
            + v_scale(iZp,iXp,iYp) * ( SXm * SYm * SZm )

       !### z velocity ###!
       !### at the position of super-droplets         ###!
       !! interpolation in face grid
       iXm = floor(ri+0.5_RP)
       iXp = iXm + 1
       sXm = (ri+0.5_RP) - real(iXm,kind=RP)
       sXp = 1.0_RP - sXm

       iYm = floor(rj+0.5_RP)
       iYp = iYm + 1
       sYm = (rj+0.5_RP) - real(iYm,kind=RP)
       sYp = 1.0_RP - sYm

       iZm = floor(rk)
       iZp = iZm + 1
       sZm = rk - real(iZm,kind=RP)
       sZp = 1.0_RP - sZm

       sd_w = w_scale(iZm,iXm,iYm) * ( SXp * SYp * SZp )              &
            + w_scale(iZm,iXp,iYm) * ( SXm * SYp * SZp )              &
            + w_scale(iZm,iXm,iYp) * ( SXp * SYm * SZp )              &
            + w_scale(iZm,iXp,iYp) * ( SXm * SYm * SZp )              &
            + w_scale(iZp,iXm,iYm) * ( SXp * SYp * SZm )              &
            + w_scale(iZp,iXp,iYm) * ( SXm * SYp * SZm )              &
            + w_scale(iZp,iXm,iYp) * ( SXp * SYm * SZm )              &
            + w_scale(iZp,iXp,iYp) * ( SXm * SYm * SZm )

       ! the impact of terminal veloicty
       
       sd_vz = sd_vzw(n)

       sd_vzw(n) = sd_w - sd_vz
    end do

    ! Section specification for fapp profiler
    call fapp_stop("sdm_getvel",1,1)

    return
  end subroutine sdm_getvel
  !----------------------------------------------------------------------------
  subroutine sdm_getvz(pres,rhod,temp,            &
       sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,   &
       ilist_s,ilist_m,ilist_l,ptype        )
    ! evaluate the terminal velocities
    use scale_grid_index, only: &
         IA,JA,KA
    use scale_grid, only: &
         FX => GRID_FX, &
         FY => GRID_FY
    use scale_const, only: &
         grav_mks => CONST_GRAV, &
         rhow_mks => CONST_DWATR, &
         pstd_mks => CONST_Pstd, &
         p0 => CONST_PRE00, &
         t0 => CONST_TEM00, &
         cp => CONST_CPdry, &
         rd => CONST_Rdry
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID, O_SIX, vz_b, vz_c, doprecipitation
    ! Input variables
    real(RP), intent(in) :: pres(KA,IA,JA)  ! Pressure
    real(RP), intent(in) :: rhod(KA,IA,JA)  ! dry air density
    real(RP), intent(in) :: temp(KA,IA,JA)  ! temperature
    integer, intent(in) :: sd_num           ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    character(len=6), intent(in) :: ptype   ! process type : 'motion process' or 'stochastic coalescence process'
    ! Output variables
    real(RP), intent(out) :: sd_vz(1:sd_num)! terminal velocity of super-droplets in real space
    integer, intent(out) :: ilist_s(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_m(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_l(1:sd_num)  ! buffer for list vectorization

    ! Local parameters
    real(RP), parameter :: grav = 9.80665E+2_RP       ! gravity acceleration [cm/s^2]
    real(RP), parameter :: rhow = 1.0_RP              ! density of water [g/cm^3]
    real(RP), parameter :: lb4l = 6.62E-6_RP          ! base length[cm] for mean free path  of air molecules
    real(RP), parameter :: pb4l = 1013.25_RP          ! base pressure[hPa] for mean free path of air molecules
    real(RP), parameter :: tb4l = 293.15_RP           ! base temperarure[20C] for mean free path of air molecules
    real(RP), parameter :: visb4l = 1.818E-4_RP       ! base dynamic viscosity[g/(cm*s)] for mean free path of air molecules
    real(RP), parameter :: m2cm = 1.E+2_RP            ! converter of unit [m] => [cm]
    real(RP), parameter :: cm2m = 1.E-2_RP            ! converter of unit [cm] => [m]
    real(RP), parameter :: pa2hpa = 1.E-2_RP          ! converter of unit [Pa] => [hPa]
    real(RP), parameter :: kgm2gcm = 1.E-3_RP         ! converter of unit [kg/m^3] => [g/cm^3]

    ! Work variables
    real(RP) :: sd_dia    ! diameter[cm] of super-droplets
    real(RP) :: sd_l      ! mean free path[cm] of air molecules at the location of super-droplets
    real(RP) :: sd_p      ! pressure[hPa] at the location of super-droplets
    real(RP) :: sd_pt     ! potential temperature[K] at the location of super-droplets
    real(RP) :: sd_rd     ! air density[g/cm^3] at the location of super-droplets
    real(RP) :: sd_sigma  ! water tension[N/m*10^3=g/s^2] of super-droplets
    real(RP) :: sd_t      ! temperature[K] at the location of super-droplets
    real(RP) :: sd_visco  ! dynamic viscosity[g/(cm*s)] at the location of super-droplets
    real(RP) :: bond      ! modified Bond number
    real(RP) :: csc       ! slip correction factor
    real(RP) :: gxdrow    ! = grav * (rhow-sd_rd)
    real(RP) :: nda       ! Davies number
    real(RP) :: npp       ! physical property number
    real(RP) :: nre       ! Reynolds number
    real(RP) :: p0iv      ! inverse of p0
    real(RP) :: rddvcp    ! rd / cp

    real(RP) :: c1        ! temporary
    real(RP) :: c2        ! temporary
    real(RP) :: c3        ! temporary
    real(RP) :: x1        ! temporary
    real(RP) :: x2        ! temporary
    real(RP) :: x3        ! temporary
    real(RP) :: y         ! temporary
    real(RP) :: tau       ! temporary
    real(RP) :: tdeg      ! temporary

    real(RP) :: ri        ! real index [i] of super-droplets
    real(RP) :: rj        ! real index [j] of super-droplets
    real(RP) :: rk        ! real index [k] of super-droplets

    real(RP) :: sXm       ! variable for inteporation
    real(RP) :: sXp       ! variable for inteporation
    real(RP) :: sYm       ! variable for inteporation
    real(RP) :: sYp       ! variable for inteporation
    real(RP) :: sZm       ! variable for inteporation
    real(RP) :: sZp       ! variable for inteporation

    integer :: iXm        ! index for inteporation
    integer :: iXp        ! index for inteporation
    integer :: iYm        ! index for inteporation
    integer :: iYp        ! index for inteporation
    integer :: iZm        ! index for inteporation
    integer :: iZp        ! index for inteporation

    integer :: tlist_s         ! total list number for small
    integer :: tlist_m         ! total list number for middle
    integer :: tlist_l         ! total list number for large

    integer :: scnt            ! counter for small
    integer :: mcnt            ! counter for middle
    integer :: lcnt            ! counter for large

    integer :: i, j, k, n, m   ! index
    !---------------------------------------------------------------------

    ! Section specification for fapp profiler
    call fapp_start("sdm_getvz",1,1)

    if( .not. doprecipitation ) then
       sd_vz(:) = 0.0_RP
       return
    endif

    call sdm_x2ri(sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    
    tlist_s = 0
    tlist_m = 0
    tlist_l = 0

    p0iv   = 1.0_RP/real(p0*pa2hpa,kind=RP)
    rddvcp = real(rd,kind=RP)/real(cp,kind=RP)

    ! Get index list for compressing buffer.
    scnt = 0
    mcnt = 0
    lcnt = 0      

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle

       sd_dia = 2.0_RP * sd_r(n) * m2cm   !! diameter [m=>cm]

       if( sd_dia<=1.9E-3_RP ) then

          !== small cloud droplets ==!

          scnt = scnt + 1
          ilist_s(scnt) = n

       else if( sd_dia>1.9E-3_RP .and. sd_dia<=1.07E-1_RP ) then

          !== large cloud droplets and small raindrops ==!

          mcnt = mcnt + 1
          ilist_m(mcnt) = n

       else

          !== raindrops ==!

          lcnt = lcnt + 1
          ilist_l(lcnt) = n

       end if

    end do

    tlist_s = scnt
    tlist_m = mcnt
    tlist_l = lcnt

    ! Calculate terminal velocity of super-droplets

    if( ptype .eq. 'motion' ) then

       !### small cloud droplets ###!

       if( tlist_s>0 ) then

!OCL NORECURRENCE
          do m=1,tlist_s
             n = ilist_s(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             !== diameter[cm] of droplets ==!
             
             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhod(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhod(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhod(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhod(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhod(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhod(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhod(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhod(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_rd = sd_rd * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p  = pres(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + pres(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + pres(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + pres(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + pres(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + pres(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + pres(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + pres(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_p  = sd_p * pa2hpa

             !== temperarure ==!

             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                  + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0 ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * dsqrt(sd_t/tb4l)

             !== terminal vecocity ==!

             c1  = gxdrow / ( 18.0_RP * sd_visco )
             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             sd_vz(n) = c1 * csc * sd_dia * sd_dia

             sd_vz(n) = sd_vz(n) * cm2m

          end do

       end if

       !### large cloud droplets and small raindrops ###!

       if( tlist_m>0 ) then

!OCL NORECURRENCE
          do m=1,tlist_m
             n = ilist_m(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             !== diameter[cm] of droplets ==!
             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhod(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhod(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhod(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhod(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhod(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhod(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhod(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhod(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_rd = sd_rd * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p  = pres(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + pres(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + pres(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + pres(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + pres(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + pres(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + pres(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + pres(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_p  = sd_p * pa2hpa

             !== temperarure ==!

             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                           + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0 ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * dsqrt(sd_t/tb4l)

             !== Davies number ==!

             c2 = sd_rd * (4.0_RP*gxdrow)/(3.0_RP*sd_visco*sd_visco)

             nda = c2 * sd_dia * sd_dia * sd_dia

             !== Reynolds number ==!

             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             x1 = log(nda)
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_b(1)                                       &
                  + vz_b(2) * (x1)                                &
                  + vz_b(3) * (x2)                                &
                  + vz_b(4) * (x3)                                &
                  + vz_b(5) * (x2*x2)                             &
                  + vz_b(6) * (x3*x2)                             &
                  + vz_b(7) * (x3*x3)

             nre = csc * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do

       end if

       !### raindrops ###!

       if( tlist_l>0 ) then

!OCL NORECURRENCE
          do m=1,tlist_l
             n = ilist_l(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm
             
             !== density of air ==!

             sd_rd = rhod(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhod(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhod(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhod(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhod(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhod(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhod(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhod(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_rd = sd_rd * kgm2gcm
             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p  = pres(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + pres(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + pres(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + pres(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + pres(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + pres(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + pres(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + pres(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_p  = sd_p * pa2hpa

             !== temperarure ==!

             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                  + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0 ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== surface tension [N/m*10^3=g/s^2]   ==!
             !==               (Holten et al. 2005) ==!

             tau = 1.0_RP - sd_t/647.0960_RP
             sd_sigma = 0.2358_RP * exp( 1.2560_RP*log(tau) )      &
                  * ( 1.0_RP - 0.6250_RP*tau )

             if( sd_t<267.5_RP ) then

                tau = dtanh( (sd_t-243.9_RP)/35.35_RP )

                sd_sigma = sd_sigma - 2.854E-3_RP * tau + 1.666E-3_RP

             end if

             sd_sigma = sd_sigma * 1.E+3_RP  !! N/m => g/s^2

             !== modified Bond number ==!

             c3 = (4.0_RP*gxdrow)/(3.0_RP*sd_sigma)

             bond = c3 * sd_dia * sd_dia

             !== physical property number ==!

             npp = (sd_rd*sd_rd) * (sd_sigma*sd_sigma*sd_sigma) &
                  / (gxdrow*sd_visco*sd_visco*sd_visco*sd_visco)

             npp = exp( O_SIX * log(npp) )

             !== Reynolds number ==!

             x1 = log( bond * npp )
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_c(1)                                       &
                  + vz_c(2) * (x1)                                &
                  + vz_c(3) * (x2)                                &
                  + vz_c(4) * (x3)                                &
                  + vz_c(5) * (x2*x2)                             &
                  + vz_c(6) * (x3*x2)

             nre = npp * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          
       end if

    else if( ptype .eq. 'coales' ) then

       !### small cloud droplets ###!

       if( tlist_s>0 ) then

!OCL NORECURRENCE
          do m=1,tlist_s
             n = ilist_s(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation 
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhod(k,i,j) * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p = pres(k,i,j) * pa2hpa

             !== temperarure ==!

             sd_t = temp(k,i,j)

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!
             
             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0 ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.718_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * dsqrt(sd_t/tb4l)

             !== terminal vecocity ==!

             c1  = gxdrow / ( 18.0_RP * sd_visco )
             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             sd_vz(n) = c1 * csc * sd_dia * sd_dia

             sd_vz(n) = sd_vz(n) * cm2m

          end do

       end if

       !### large cloud droplets and small raindrops ###!

       if( tlist_m>0 ) then

!OCL NORECURRENCE
          do m=1,tlist_m
             n = ilist_m(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation 
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhod(k,i,j) * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p = pres(k,i,j) * pa2hpa

             !== temperarure ==!

             sd_t = temp(k,i,j)

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!
             
             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0 ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * dsqrt(sd_t/tb4l)

             !== Davies number ==!

             c2 = sd_rd * (4.0_RP*gxdrow)/(3.0_RP*sd_visco*sd_visco)

             nda = c2 * sd_dia * sd_dia * sd_dia

             !== Reynolds number ==!
             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             x1 = log(nda)
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_b(1)                                       &
                  + vz_b(2) * (x1)                                &
                  + vz_b(3) * (x2)                                &
                  + vz_b(4) * (x3)                                &
                  + vz_b(5) * (x2*x2)                             &
                  + vz_b(6) * (x3*x2)                             &
                  + vz_b(7) * (x3*x3)

             nre = csc * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do

       end if

       !### raindrops ###!

       if( tlist_l>0 ) then

!OCL NORECURRENCE
          do m=1,tlist_l
             n = ilist_l(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation 
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!
             
             sd_rd = rhod(k,i,j) * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p = pres(k,i,j) * pa2hpa

             !== temperarure ==!

             sd_t = temp(k,i,j)

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!
             
             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0 ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== surface tension [N/m*10^3=g/s^2]   ==!
             !==               (Holten et al. 2005) ==!

             tau = 1.0_RP - sd_t/647.0960_RP

             sd_sigma = 0.23580_RP * exp( 1.2560_RP*log(tau) )      &
                  * ( 1.00_RP - 0.6250_RP*tau )

             if( sd_t<267.50_RP ) then

                tau = dtanh( (sd_t-243.90_RP)/35.350_RP )

                sd_sigma = sd_sigma - 2.854E-3_RP * tau + 1.666E-3_RP

             end if

             sd_sigma = sd_sigma * 1.E+3_RP  !! N/m => g/s^2

             !== modified Bond number ==!

             c3 = (4.0_RP*gxdrow)/(3.0_RP*sd_sigma)

             bond = c3 * sd_dia * sd_dia

             !== physical property number ==!

             npp = (sd_rd*sd_rd) * (sd_sigma*sd_sigma*sd_sigma) &
                  / (gxdrow*sd_visco*sd_visco*sd_visco*sd_visco)

             npp = exp( O_SIX * log(npp) )

             !== Reynolds number ==!
             
             x1 = log( bond * npp )
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_c(1)                                       &
                  + vz_c(2) * (x1)                                &
                  + vz_c(3) * (x2)                                &
                  + vz_c(4) * (x3)                                &
                  + vz_c(5) * (x2*x2)                             &
                  + vz_c(6) * (x3*x2)

             nre = npp * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do
                  
       end if

    end if

    ! Section specification for fapp profiler
    call fapp_stop("sdm_getvz",1,1)

    return
  end subroutine sdm_getvz
end module m_sdm_motion
