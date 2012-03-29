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

    real(8) :: ENV_THETA  = 300.D0  ! Potential Temperature of environment [K]
    real(8) :: ENV_RH     =  80.D0  ! Relative Humidity of environment [%]
    real(8) :: LAPS_THETA =   5.D-3 ! Lapse rate of Potential Temperature [K m-1]
    real(8) :: CTH_LEVEL  =  12.D3  ! depth of the constant potential temperature layer [m]
    real(8) :: EXT_TBBL   =   5.D0  ! extremum of temperature in bubble [K]
    real(8) :: ZC_BBL     =   3.D3  ! center location [m]: z
    real(8) :: XC_BBL     =  15.D3  ! center location [m]: x
    real(8) :: YC_BBL     =  15.D3  ! center location [m]: y
    real(8) :: ZR_BBL     =   2.D3  ! bubble radius   [m]: z
    real(8) :: XR_BBL     =   4.D3  ! bubble radius   [m]: x
    real(8) :: YR_BBL     =   4.D3  ! bubble radius   [m]: y

    NAMELIST / PARAM_MKEXP_WARMBUBBLE / &
       ENV_THETA,  &
       ENV_RH,     &
       LAPS_THETA, &
       CTH_LEVEL,  &
       EXT_TBBL,   &
       ZC_BBL,     &
       XC_BBL,     &
       YC_BBL,     &
       ZR_BBL,     &
       XR_BBL,     &
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
       if ( CZ(k) < CTH_LEVEL ) then
          pott_k(k) = ENV_THETA
       else
          pott_k(k) = ENV_THETA + LAPS_THETA * ( CZ(k)-CTH_LEVEL )
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
  !> Make initial state for Kelvin Helmholtz experiment
  !-----------------------------------------------------------------------------
  subroutine MKEXP_turbdyn
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
          qv_k  (k) = 1.D-3 * fact
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
  end subroutine MKEXP_turbdyn

end module mod_mkexp