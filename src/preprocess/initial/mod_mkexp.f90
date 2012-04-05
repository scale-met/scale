!-------------------------------------------------------------------------------
!> module initial
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-03-27 (H.Yashiro)  [mod] change subroutines into one module 
!!
!<
!-------------------------------------------------------------------------------
module mod_mkexp
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
  public :: MKEXP_coldbubble
  public :: MKEXP_warmbubble
  public :: MKEXP_turbdyn
  public :: MKEXP_khwave
  public :: MKEXP_squalline

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
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_coldbubble
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI    => CONST_PI,    &
       RovCP => CONST_RovCP, &
       Pstd  => CONST_Pstd,  &
       P00   => CONST_PRE00
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       CX  => GRID_CX, &
       CY  => GRID_CY, &
       FDZ => GRID_FDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho
    implicit none

    real(8) :: ENV_THETA =  300.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_QTRC  = 0.001D0 ! tracer of environment
    real(8) :: EXT_TBBL  =   -5.D0 ! extremum of temperature in bubble [K]
    real(8) :: ZC_BBL    =    3.D3 ! center location [m]: z
    real(8) :: XC_BBL    =   15.D3 ! center location [m]: x
    real(8) :: YC_BBL    =   15.D3 ! center location [m]: y
    real(8) :: ZR_BBL    =    2.D3 ! bubble radius   [m]: z
    real(8) :: XR_BBL    =    1.D3 ! bubble radius   [m]: x
    real(8) :: YR_BBL    =    1.D3 ! bubble radius   [m]: y

    NAMELIST / PARAM_MKEXP_COLDBUBBLE / &
       ENV_THETA, &
       ENV_QTRC,  &
       EXT_TBBL,  &
       ZC_BBL,    &
       XC_BBL,    &
       YC_BBL,    &
       ZR_BBL,    &
       XR_BBL,    &
       YR_BBL

    real(8) :: dens_k(KA) ! density  [kg/m3]
    real(8) :: pres_k(KA) ! pressure [Pa]
    real(8) :: pott_k(KA) ! potential temperature [K]

    real(8) :: pres_sfc
    real(8) :: pott_sfc
    real(8) :: dist

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[COLDBUBBLE]/Categ[MKEXP]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_COLDBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_COLDBUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_COLDBUBBLE)

    pres_sfc = Pstd
    pott_sfc = ENV_THETA
    do k = KS, KE
       pott_k(k) = ENV_THETA
    enddo

    ! make density profile
    call hydro_buildrho( dens_k(:), pres_k(:), pott_k(:), pres_sfc, pott_sfc )

    ! make cold bubble
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dist = ( (CZ(k)-ZC_BBL)/ZR_BBL )**2.D0 &
            + ( (CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (CY(j)-YC_BBL)/YR_BBL )**2.D0

       if ( dist <= 1.D0 ) then
          RHOT(k,i,j) = dens_k(k) &
                      * ( pott_k(k) + EXT_TBBL * cos( 0.5D0*PI*sqrt(dist) )**2 )
       else
          RHOT(k,i,j) = dens_k(k) * pott_k(k)
       endif

       DENS(k,i,j) = dens_k(k)
       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = 0.D0
       MOMY(k,i,j) = 0.D0

       do iq = 1,  QA
          QTRC(k,i,j,iq) = ENV_QTRC
       enddo

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_coldbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_warmbubble
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,     &
       Rdry   => CONST_Rdry,   &
       RovCP  => CONST_RovCP,  &
       EPSvap => CONST_EPSvap, &
       Pstd   => CONST_Pstd,   &
       P00    => CONST_PRE00,  &
       Rvap   => CONST_Rvap,   &
       CPvap  => CONST_CPvap,  &
       CL     => CONST_CL,     &
       LH0    => CONST_LH00,   &
       PSAT0  => CONST_PSAT0,  &
       T00    => CONST_TEM00
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       CX  => GRID_CX, &
       CY  => GRID_CY, &
       FDZ => GRID_FDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho
    implicit none

    real(8) :: ENV_THETA  = 300.0D0  ! Potential Temperature of environment [K]
    real(8) :: ENV_RH     =  80.0D0  ! Relative Humidity of environment [%]
    real(8) :: LAPS_THETA1=   4.0D-3 ! Lapse rate of Potential Temperature in 1st layer [K m-1]
    real(8) :: LAPS_THETA2=   3.0D-2 ! Lapse rate of Potential Temperature in 2nd layer [K m-1]
    real(8) :: CTH_LEVEL1 =   1.0D3  ! depth of 1st constant potential temperature layer [m]
    real(8) :: CTH_LEVEL2 =  12.0D3  ! depth of 2nd constant potential temperature layer [m]
    real(8) :: EXT_TBBL   =   1.0D0  ! extremum of temperature in bubble [K]
    real(8) :: ZC_BBL     =   0.5D3  ! center location [m]: z
    real(8) :: XC_BBL     =  40.0D3  ! center location [m]: x
    real(8) :: YC_BBL     =  40.0D3  ! center location [m]: y
    real(8) :: ZR_BBL     =   3.0D3  ! bubble radius   [m]: z
    real(8) :: XR_BBL     =  10.0D3  ! bubble radius   [m]: x
    real(8) :: YR_BBL     =  10.0D3  ! bubble radius   [m]: y

    NAMELIST / PARAM_MKEXP_WARMBUBBLE / &
       ENV_THETA,   &
       ENV_RH,      &
       LAPS_THETA1, &
       LAPS_THETA2, &
       CTH_LEVEL1,  &
       CTH_LEVEL2,  &
       EXT_TBBL,    &
       ZC_BBL,      &
       XC_BBL,      &
       YC_BBL,      &
       ZR_BBL,      &
       XR_BBL,      &
       YR_BBL

    real(8) :: dens_k(KA) ! density  [kg/m3]
    real(8) :: pres_k(KA) ! pressure [Pa]
    real(8) :: pott_k(KA) ! potential temperature [K]
    real(8) :: temp_k(KA) ! temperature [K]
    real(8) :: qv_k  (KA) ! water vapor mixing ratio [kg/kg]

    real(8) :: pres_sfc
    real(8) :: pott_sfc
    real(8) :: psat, qsat
    real(8) :: dist
    real(8) :: Tmin = 10.D0

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[WARMBUBBLE]/Categ[MKEXP]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_WARMBUBBLE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_WARMBUBBLE)

    pres_sfc = Pstd
    pott_sfc = ENV_THETA
    do k = KS, KE
       if ( CZ(k) <= CTH_LEVEL1 ) then
          pott_k(k) = ENV_THETA
       else if ( CZ(k) <= CTH_LEVEL2 ) then
          pott_k(k) = pott_k(k-1) + LAPS_THETA1 * ( CZ(k)-CZ(k-1) )
       else
          pott_k(k) = pott_k(k-1) + LAPS_THETA2 * ( CZ(k)-CZ(k-1) )
       end if
    enddo

    ! make density profile
    call hydro_buildrho( dens_k(:), pres_k(:), pott_k(:), pres_sfc, pott_sfc )

    do k = KS, KE
       temp_k(k) = pres_k(k) / ( dens_k(k) * Rdry )
       psat = PSAT0 * ( max(temp_k(k),Tmin)/T00 ) ** ( ( CPvap-CL )/Rvap ) &
            * exp ( LH0/Rvap * ( 1.0D0/T00 - 1.0D0/max(temp_k(k),Tmin) ) )
       qsat = EPSvap * psat / ( pres_k(k) - ( 1.D0-EPSvap )*psat )

       qv_k(k) = ENV_RH*1.D-2 * qsat
    enddo

    ! make warm bubble
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dist = ( (CZ(k)-ZC_BBL)/ZR_BBL )**2.D0 &
            + ( (CX(i)-XC_BBL)/XR_BBL )**2.D0 &
            + ( (CY(j)-YC_BBL)/YR_BBL )**2.D0

       if ( dist <= 1.D0 ) then
          RHOT(k,i,j) = dens_k(k) &
                      * ( pott_k(k) + EXT_TBBL * cos( 0.5D0*PI*sqrt(dist) )**2 )
       else
          RHOT(k,i,j) = dens_k(k) * pott_k(k)
       endif

       DENS(k,i,j) = dens_k(k)
       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = 0.D0
       MOMY(k,i,j) = 0.D0

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo
       QTRC(k,i,j,I_QV) = qv_k(k)

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_warmbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for turbulence experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_turbdyn
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,     &
       Rdry   => CONST_Rdry,   &
       RovCP  => CONST_RovCP,  &
       EPSvap => CONST_EPSvap, &
       Pstd   => CONST_Pstd,   &
       P00    => CONST_PRE00,  &
       Rvap   => CONST_Rvap,   &
       CPvap  => CONST_CPvap,  &
       CL     => CONST_CL,     &
       LH0    => CONST_LH00,   &
       PSAT0  => CONST_PSAT0,  &
       T00    => CONST_TEM00
    use mod_random, only: &
       RANDOM_reset, &
       RANDOM_get
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       CX  => GRID_CX, &
       CY  => GRID_CY, &
       FDZ => GRID_FDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho
    implicit none

    real(8) :: ENV_THETA  = 300.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_LTHETA =   4.D-3! Potential Temperature lapse rate [K m-1]
    real(8) :: ENV_RH     = 50.D0  ! Relative Humidity of environment [%]
    real(8) :: ENV_XVEL   =  5.D0  ! environment x-velocity [m s-1]
    real(8) :: ENV_YVEL   =  0.D0  ! environment y-velocity [m s-1]

    NAMELIST / PARAM_MKEXP_TURBDYN / &
       ENV_THETA,  &
       ENV_LTHETA, &
       ENV_RH,     &
       ENV_XVEL,   &
       ENV_YVEL

    real(8) :: dens_k(KA) ! density  [kg/m3]
    real(8) :: pres_k(KA) ! pressure [Pa]
    real(8) :: pott_k(KA) ! potential temperature [K]
    real(8) :: temp_k(KA) ! absolute temperature [K]
    real(8) :: velx_k(KA) ! x-velocity [m s-1]
    real(8) :: vely_k(KA) ! y-velocity [m s-1]
    real(8) :: qv_k  (KA) ! water vapor mixing ratio [kg/kg]

    real(8) :: pres_sfc
    real(8) :: pott_sfc

    real(8) :: rndm1(KA,IA,JA)
    real(8) :: rndm2(KA,IA,JA)
    real(8) :: fact, useed, tseed
    real(8) :: psat, qsat
    real(8) :: dist
    real(8) :: Tmin = 10.D0

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TURBDYN]/Categ[MKEXP]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_TURBDYN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_TURBDYN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_TURBDYN)

    pres_sfc = Pstd
    pott_sfc = ENV_THETA
    do k = KS, KE
       velx_k(k) = ENV_XVEL
       vely_k(k) = ENV_YVEL
       pott_k(k) = ENV_THETA + ENV_LTHETA * CZ(k)
    enddo

    ! make density profile
    call hydro_buildrho( dens_k(:), pres_k(:), pott_k(:), pres_sfc, pott_sfc )

    do k = KS, KE
       temp_k(k) = pres_k(k) / ( dens_k(k) * Rdry )
       psat = PSAT0 * ( max(temp_k(k),Tmin)/T00 ) ** ( ( CPvap-CL )/Rvap ) &
            * exp ( LH0/Rvap * ( 1.0D0/T00 - 1.0D0/max(temp_k(k),Tmin) ) )
       qsat = EPSvap * psat / ( pres_k(k) - ( 1.D0-EPSvap )*psat )

       qv_k(k) = ENV_RH*1.D-2 * qsat
    enddo

    call RANDOM_get(rndm1)
    call RANDOM_reset
    call RANDOM_get(rndm2)
    useed = ( dsqrt( ENV_XVEL**2+ENV_YVEL**2 ) ) * 1.D-2
    tseed = ENV_THETA * 1.D-2

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       DENS(k,i,j) = dens_k(k)
       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = ( velx_k(k) + useed*rndm1(k,i,j) ) * dens_k(k)
       MOMY(k,i,j) = 0.D0
       RHOT(k,i,j) = ( pott_k(k) + tseed*rndm2(k,i,j) ) * dens_k(k)

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo
       QTRC(k,i,j,I_QV) = qv_k(k)

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*) 

    return
  end subroutine MKEXP_turbdyn

  !-----------------------------------------------------------------------------
  !> Make initial state for Kelvin Helmholtz experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_khwave
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       Pstd => CONST_Pstd
    use mod_random, only: &
       RANDOM_reset, &
       RANDOM_get
    use mod_grid, only : &
       CZ  => GRID_CZ, &
       CX  => GRID_CX, &
       CY  => GRID_CY, &
       FDZ => GRID_FDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho
    implicit none

    real(8) :: ENV_THETA  = 300.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_DTHETA =   5.D0 ! Potential Temperature of environment [K]
    real(8) :: ENV_XVEL2  =  20.D0 ! environment x-velocity in layer 2 [m s-1]
    real(8) :: ENV_XVEL1  =   0.D0 ! environment x-velocity in layer 1 [m s-1]
    real(8) :: LEV_XVEL2  = 2.05D3 ! level at which x-velocity changes [m]
    real(8) :: LEV_XVEL1  = 1.95D3 ! level at which x-velocity changes [m]

    NAMELIST / PARAM_MKEXP_TURBDYN / &
       ENV_THETA,  &
       ENV_DTHETA, &
       ENV_XVEL1,  &
       ENV_XVEL2,  &
       LEV_XVEL1,  &
       LEV_XVEL2

    real(8) :: dens_k(KA) ! density  [kg/m3]
    real(8) :: pres_k(KA) ! pressure [Pa]
    real(8) :: pott_k(KA) ! potential temperature [K]
    real(8) :: velx_k(KA) ! velocity u [m/s]
    real(8) :: qv_k  (KA) ! water vapor mixing ratio [kg/kg]

    real(8) :: pres_sfc
    real(8) :: pott_sfc

    real(8) :: rndm1(KA,IA,JA)
    real(8) :: rndm2(KA,IA,JA)
    real(8) :: fact, useed, tseed

    integer :: ierr
    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[TURBDYN]/Categ[MKEXP]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_TURBDYN,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_TURBDYN. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_TURBDYN)

    pres_sfc = Pstd
    pott_sfc = ENV_THETA
    do k = KS, KE
       if    ( CZ(k) < LEV_XVEL1 ) then
          pott_k(k) = ENV_THETA
          velx_k(k) = ENV_XVEL1
          qv_k  (k) = 1.D-3
       elseif( CZ(k) > LEV_XVEL2 ) then
          pott_k(k) = ENV_THETA + ENV_DTHETA
          velx_k(k) = ENV_XVEL2
          qv_k  (k) = 0.D0
       else
          fact = ( CZ(k)-LEV_XVEL1 ) / ( LEV_XVEL2-LEV_XVEL1 ) 

          pott_k(k) = ENV_THETA + ENV_DTHETA * fact
          velx_k(k) = ENV_XVEL1 + ( ENV_XVEL2-ENV_XVEL1 ) * fact
          qv_k  (k) = 1.D-3 * ( 1.D0 - fact )
       endif
    enddo

    ! make density profile
    call hydro_buildrho( dens_k(:), pres_k(:), pott_k(:), pres_sfc, pott_sfc )

    call RANDOM_get(rndm1)
    call RANDOM_reset
    call RANDOM_get(rndm2)
    useed = ( ENV_XVEL2+ENV_XVEL1 ) * 1.D-2
    tseed = ENV_THETA * 1.D-2

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       DENS(k,i,j) = dens_k(k)
       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = ( velx_k(k) + useed*rndm1(k,i,j) ) * dens_k(k)
       MOMY(k,i,j) = 0.D0
       RHOT(k,i,j) = ( pott_k(k) + tseed*rndm2(k,i,j) ) * dens_k(k)

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo
       QTRC(k,i,j,I_QV) = qv_k(k)

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*) 

    return
  end subroutine MKEXP_khwave

  !-----------------------------------------------------------------------------
  !> Make initial state for squalline experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_squalline
    use mod_stdio, only: &
       IO_get_available_fid, &
       IO_FILECHR,  &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only : &
       PI     => CONST_PI,     &
       GRAV   => CONST_GRAV,   &
       Rdry   => CONST_Rdry,   &
       RovCP  => CONST_RovCP,  &
       RovCV  => CONST_RovCV,  &
       CVovCP => CONST_CVovCP, &
       CPovCV => CONST_CPovCV, &
       EPSvap => CONST_EPSvap, &
       Pstd   => CONST_Pstd,   &
       P00    => CONST_PRE00,  &
       Rvap   => CONST_Rvap,   &
       CPvap  => CONST_CPvap,  &
       CL     => CONST_CL,     &
       LH0    => CONST_LH00,   &
       PSAT0  => CONST_PSAT0,  &
       T00    => CONST_TEM00
    use mod_grid, only : &
       CZ  => GRID_CZ,  &
       CX  => GRID_CX,  &
       CY  => GRID_CY,  &
       FDZ => GRID_FDZ
    use mod_atmos_vars, only: &
       DENS, &
       MOMX, &
       MOMY, &
       MOMZ, &
       RHOT, &
       QTRC
    use mod_atmos_hydrostatic, only: &
       hydro_buildrho
    implicit none

    character(len=IO_FILECHR) :: ENV_IN_SOUNDING_file = ''
    real(8) :: REF_XVEL   =  12.D0  ! reference x-velocity [m s-1]
    real(8) :: REF_YVEL   =  -2.D0  ! reference y-velocity [m s-1]
    real(8) :: EXT_TBBL   =   5.D0  ! extremum of temperature in bubble [K]
    real(8) :: EXT_QBBL   =   8.D-4 ! extremum of water vapor mixing ratio in bubble [kg kg-1
    real(8) :: ZC_BBL     =   3.D3  ! center location [m]: z
    real(8) :: XC_BBL     =  15.D3  ! center location [m]: x
    real(8) :: YC_BBL     =  15.D3  ! center location [m]: y
    real(8) :: ZR_BBL     =   2.5D3 ! bubble radius   [m]: z
    real(8) :: XR_BBL     =   7.D3  ! bubble radius   [m]: x
    real(8) :: YR_BBL     =   6.D3  ! bubble radius   [m]: y
    real(8) :: DY_BBL     =  15.D3  ! distance of bubbles [m]
    integer :: NUM_BBL    =   4     ! number of bubbles [-]

    NAMELIST / PARAM_MKEXP_SQUALLINE / &
       ENV_IN_SOUNDING_file, &
       REF_XVEL,   &
       REF_YVEL,   &
       EXT_TBBL,   &
       EXT_QBBL,   &
       ZC_BBL,     &
       XC_BBL,     &
       YC_BBL,     &
       ZR_BBL,     &
       XR_BBL,     &
       YR_BBL,     &
       DY_BBL,     &
       NUM_BBL

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(8) :: EXP_z   (EXP_klim) ! height      [m]
    real(8) :: EXP_dens(EXP_klim) ! density     [kg/m3]
    real(8) :: EXP_pott(EXP_klim) ! potential temperature [K]
    real(8) :: EXP_u   (EXP_klim) ! velocity u  [m/s]
    real(8) :: EXP_v   (EXP_klim) ! velocity v  [m/s]
    real(8) :: EXP_qv  (EXP_klim) ! water vapor [g/kg]
    real(8) :: EXP_pres(EXP_klim) ! pressure    [Pa]

    real(8) :: dens_k(KA) ! density  [kg/m3]
    real(8) :: velx_k(KA) ! x-velocity [m s-1]
    real(8) :: vely_k(KA) ! y-velocity [m s-1]
    real(8) :: pres_k(KA) ! pressure [Pa]
    real(8) :: pott_k(KA) ! potential temperature [K]
    real(8) :: temp_k(KA) ! temperature [K]
    real(8) :: qv_k  (KA) ! water vapor mixing ratio [kg/kg]

    real(8) :: pres_sfc
    real(8) :: pott_sfc

    real(8) :: fact1, fact2, rdz
    real(8) :: RovP, EXP_dens_s, dens_s, dhyd, dgrd
    real(8) :: dist, YP_BBL

    real(8), parameter :: criteria = 1.D-10
    integer, parameter :: itelim = 100
    integer            :: kref, ite

    integer :: ierr, fid
    integer :: k, i, j, m, iq
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ START MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SQUALLINE]/Categ[MKEXP]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKEXP_SQUALLINE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKEXP_SQUALLINE. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_MKEXP_SQUALLINE)

    !--- prepare sounding profile
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding file:', trim(ENV_IN_SOUNDING_file)
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

       if ( ierr /= 0 ) then
          if( IO_L ) write(*,*) 'xxx Input file not found!'
       endif

       !--- read sounding file till end
       read(fid,*) EXP_pres(1), EXP_pott(1), EXP_qv(1)

       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pressure [hPa]',     EXP_pres(1)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface pot. temp  [K]',     EXP_pott(1)
       if( IO_L ) write(IO_FID_LOG,*) '+++ Surface water vapor [g/kg]', EXP_qv(1)

       do k = 2, EXP_klim
          read(fid,*,iostat=ierr) EXP_z(k), EXP_pott(k), EXP_qv(k), EXP_u(k), EXP_v(k)
          if ( ierr /= 0 ) exit
       enddo

       EXP_kmax = k - 1
    close(fid)

    EXP_pres(1) = EXP_pres(1) * 1.D+2
    EXP_dens(1) = P00 / Rdry / EXP_pott(1) * ( EXP_pres(1)/P00 )**CVovCP

    EXP_z(1)    = 0.D0
    EXP_u(1)    = EXP_u(2)
    EXP_v(1)    = EXP_v(2)
    do k = 1, EXP_klim
       EXP_qv(k) = EXP_qv(k) * 1.D-3
    enddo

    ! kref=2: lowermost layer
    do kref = 2, EXP_kmax
       rdz = 1.D0 / ( EXP_z(kref)-EXP_z(kref-1) )

       EXP_dens_s     = 0.D0
       EXP_dens(kref) = EXP_dens(kref-1) ! first guess

       do ite = 1, itelim
          if ( abs(EXP_dens(kref) - EXP_dens_s) <= criteria ) exit

          EXP_dens_s = EXP_dens(kref)

          dhyd = + ( P00 * ( EXP_dens(kref-1) * Rdry * EXP_pott(kref-1) / P00 )**CPovCV &
                   - P00 * ( EXP_dens_s       * Rdry * EXP_pott(kref)   / P00 )**CPovCV ) * rdz & ! dp/dz
                 - GRAV * 0.5D0 * ( EXP_dens(kref-1) + EXP_dens_s )                               ! rho*g

          dgrd = - P00 * ( Rdry * EXP_pott(kref) / P00 )**CPovCV * rdz &
                 * CPovCV * EXP_dens_s**RovCV                          &
                 - 0.5D0 * GRAV

          EXP_dens(kref) = EXP_dens_s - dhyd/dgrd
       enddo

       EXP_pres(kref) = P00 * ( EXP_dens(kref) * Rdry * EXP_pott(kref) / P00 )**CPovCV
    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Input sounding profiles'
    if( IO_L ) write(IO_FID_LOG,*) '  K       Z[m] rho[kg/m3]      P[Pa]   theta[K]     U[m/s]     V[m/s]  Qv[kg/kg]'
    do k = 1, EXP_kmax
       if( IO_L ) write(IO_FID_LOG,'(1x,i3,7(1x,F10.3))') &
       k, EXP_z(k), EXP_dens(k), EXP_pres(k), EXP_pott(k), EXP_u(k), EXP_v(k), EXP_qv(k)
    enddo

    !--- linear interpolate to model grid
    do k    = KS, KE
    do kref = 2, EXP_kmax

       if (       CZ(k) >  EXP_z(kref-1) &
            .AND. CZ(k) <= EXP_z(kref)   ) then

          fact1 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )
          fact2 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
          rdz   = 1.D0 / ( CZ(k)-EXP_z(kref-1) )

          pott_k(k) = EXP_pott(kref-1) * fact1 &
                    + EXP_pott(kref)   * fact2
          velx_k(k) = EXP_u   (kref-1) * fact1 &
                    + EXP_u   (kref)   * fact2
          vely_k(k) = EXP_v   (kref-1) * fact1 &
                    + EXP_v   (kref)   * fact2
          qv_k(k)   = EXP_qv  (kref-1) * fact1 &
                    + EXP_qv  (kref)   * fact2

          velx_k(k) = velx_k(k) - REF_XVEL
          vely_k(k) = vely_k(k) - REF_YVEL

       endif
    enddo
    enddo
    pres_sfc = EXP_pres(1)
    pott_sfc = EXP_pott(1)

    ! make density profile
    call hydro_buildrho( dens_k(:), pres_k(:), pott_k(:), pres_sfc, pott_sfc )

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Interpolated data'
    if( IO_L ) write(IO_FID_LOG,*) '  K       Z[m] rho[kg/m3]      P[Pa]   theta[K]     U[m/s]     V[m/s]  Qv[kg/kg]'
    do k = KS, KE
       if( IO_L ) write(IO_FID_LOG,'(1x,i3,7(1x,F10.3))') k, CZ(k), dens_k(k), pres_k(k), pott_k(k), velx_k(k), vely_k(k), qv_k(k)
    enddo

    ! make warm bubble
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       do iq = 1, QA
          QTRC(k,i,j,iq) = 0.D0
       enddo
       RHOT(k,i,j) = dens_k(k) * pott_k(k)
       QTRC(k,i,j,I_QV) = qv_k(k)

       do m = 1, NUM_BBL

          YP_BBL = YC_BBL + DY_BBL * (m-1)

!          dist = ( (CZ(k)-ZC_BBL)/ZR_BBL )**2.D0 &
          dist = ( (CX(i)-XC_BBL)/XR_BBL )**2.D0 &
               + ( (CY(j)-YP_BBL)/YR_BBL )**2.D0

          if ( dist <= 1.D0 .and. CZ(k) <= ZR_BBL ) then
             RHOT(k,i,j) = dens_k(k) &
                         * ( pott_k(k) - EXT_TBBL * cos( 0.5D0*PI*sqrt(dist) )**2 )
             QTRC(k,i,j,I_QV) = qv_k(k) &
                         - EXT_QBBL * cos( 0.5D0*PI*sqrt(dist) )**2
          endif

       enddo

       DENS(k,i,j) = dens_k(k)
       MOMZ(k,i,j) = 0.D0
       MOMX(k,i,j) = dens_k(k) * velx_k(k)
       MOMY(k,i,j) = dens_k(k) * vely_k(k)

    enddo
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) '++++++ END MAKING INITIAL DATA ++++++'
    if( IO_L ) write(IO_FID_LOG,*)

    return
  end subroutine MKEXP_squalline
                                                                                              

end module mod_mkexp
