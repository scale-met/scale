!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-12-03 (Y.Miyamoto)  [new]
!! @li      2011-12-11 (H.Yashiro)   [mod] integrate to SCALE3
!! @li      2012-03-23 (H.Yashiro)   [mod] Explicit index parameter inclusion
!! @li      2012-04-10 (Y.Miyamoto)  [mod] introduce coefficients for interpolation
!! @li      2012-09-11 (S.Nishizawa) [mod] bugfix based on the scale document
!! @li      2012-09-12 (Y.Sato)    [renew] constant FLUX version
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_sf
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
  public :: ATMOS_PHY_SF_setup
  public :: ATMOS_PHY_SF

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'
  include 'inc_precision.h'

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
  real(RP), private, parameter :: Cm0   = 1.0E-3_RP  ! bulk coef. for U*
  real(RP), private, parameter :: visck = 1.5E-5_RP  ! kinematic viscosity 

  ! parameters
  real(RP), private, save :: Z00 =   0.0_RP      ! base
  real(RP), private, save :: Z0R = 0.018_RP      ! rough factor
  real(RP), private, save :: Z0S =  0.11_RP      ! smooth factor
  real(RP), private, save :: Zt0 =   1.4E-5_RP
  real(RP), private, save :: ZtR =   0.0_RP
  real(RP), private, save :: ZtS =   0.4_RP
  real(RP), private, save :: Ze0 =   1.3E-4_RP
  real(RP), private, save :: ZeR =   0.0_RP
  real(RP), private, save :: ZeS =  0.62_RP
  real(RP), private, save :: ThS = 300.0_RP

  ! limiter
  real(RP), private, parameter :: Ustar_min =  1.0E-3_RP ! minimum limit of U*

  real(RP), private, parameter :: Z0_min =    1.0E-5_RP ! minimum roughness length of u,v,w
  real(RP), private, parameter :: Zt_min =    1.0E-5_RP !                             T
  real(RP), private, parameter :: Ze_min =    1.0E-5_RP !                             q

  real(RP), private, save      :: Cm_min  =    1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, save      :: Ch_min  =    1.0E-5_RP !                       T
  real(RP), private, save      :: Ce_min  =    1.0E-5_RP !                       q
  real(RP), private, parameter :: Cm_max  =    2.5E-3_RP ! maximum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_max  =    1.0_RP  !                       T
  real(RP), private, parameter :: Ce_max  =    1.0_RP  !                       q

  real(RP), private, save      :: U_minM  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, save      :: U_minH  =    0.0_RP  !                   T
  real(RP), private, save      :: U_minE  =    0.0_RP  !                   q
  real(RP), private, parameter :: U_maxM  =  100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH  =  100.0_RP  !                   T
  real(RP), private, parameter :: U_maxE  =  100.0_RP  !                   q

  real(RP), private, save      :: Cm_const =  0.0011_RP ! constant bulk coef. of u,v,w
  real(RP), private, save      :: Const_SH =  15.0_RP   ! constant surface sensible flux [W/m2]
  real(RP), private, save      :: Const_LH =  115.0_RP  ! constant surface latent flux [W/m2]
  real(RP), private, save      :: Const_Ustar = 0.25_RP ! constant friction velocity [m/s]

  integer(4), private, save    :: FLG_MOM_FLUX = 0      ! 0->Bulk coef. is constant
                                                        ! 1->Friction velocity is constant

  real(RP), private, save      :: Const_FREQ = 24.0_RP ! frequency of sensible heat flux [hour]
  !  SHFLX = Const_SH[W/m^2] * sin( 2*pi*(current time)/Const_FREQ )
  logical, private, save       :: FLG_SH_DIURNAL = .false.
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    real(RP) :: ATMOS_PHY_SF_U_minM ! minimum U_abs for u,v,w
    real(RP) :: ATMOS_PHY_SF_U_minH !                   T
    real(RP) :: ATMOS_PHY_SF_U_minE !                   q
    real(RP) :: ATMOS_PHY_SF_CM_min ! minimum bulk coef. of u,v,w
    real(RP) :: ATMOS_PHY_SF_CH_min !                       T
    real(RP) :: ATMOS_PHY_SF_CE_min !                       q
    real(RP) :: ATMOS_PHY_SF_Z00
    real(RP) :: ATMOS_PHY_SF_Z0R
    real(RP) :: ATMOS_PHY_SF_Z0S
    real(RP) :: ATMOS_PHY_SF_Zt0
    real(RP) :: ATMOS_PHY_SF_ZtR
    real(RP) :: ATMOS_PHY_SF_ZtS
    real(RP) :: ATMOS_PHY_SF_Ze0
    real(RP) :: ATMOS_PHY_SF_ZeR
    real(RP) :: ATMOS_PHY_SF_ZeS
    real(RP) :: ATMOS_PHY_SF_ThS
    real(RP) :: ATMOS_PHY_SF_CM_const
    real(RP) :: ATMOS_PHY_SF_Const_SH
    real(RP) :: ATMOS_PHY_SF_Const_LH
    real(RP) :: ATMOS_PHY_SF_Const_Ustar
    real(RP) :: ATMOS_PHY_SF_Const_FREQ
    integer  :: ATMOS_PHY_SF_FLG_MOM_FLUX
    logical  :: ATMOS_PHY_SF_FLG_SH_DIURNAL

    NAMELIST / PARAM_ATMOS_PHY_SF / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_U_minH, &
       ATMOS_PHY_SF_U_minE, &
       ATMOS_PHY_SF_CM_min, &
       ATMOS_PHY_SF_CH_min, &
       ATMOS_PHY_SF_CE_min, &
       ATMOS_PHY_SF_Z00, &
       ATMOS_PHY_SF_Z0R, &
       ATMOS_PHY_SF_Z0S, &
       ATMOS_PHY_SF_Zt0, &
       ATMOS_PHY_SF_ZtR, &
       ATMOS_PHY_SF_ZtS, &
       ATMOS_PHY_SF_Ze0, &
       ATMOS_PHY_SF_ZeR, &
       ATMOS_PHY_SF_ZeS, &
       ATMOS_PHY_SF_ThS, &
       ATMOS_PHY_SF_CM_const, &
       ATMOS_PHY_SF_Const_SH, &
       ATMOS_PHY_SF_Const_LH, &
       ATMOS_PHY_SF_Const_Ustar, &
       ATMOS_PHY_SF_Const_FREQ, &
       ATMOS_PHY_SF_FLG_MOM_FLUX, &
       ATMOS_PHY_SF_FLG_SH_DIURNAL

    integer :: ierr
    !---------------------------------------------------------------------------

    ATMOS_PHY_SF_U_minM = U_minM
    ATMOS_PHY_SF_U_minH = U_minH
    ATMOS_PHY_SF_U_minE = U_minE
    ATMOS_PHY_SF_CM_min = CM_min
    ATMOS_PHY_SF_CH_min = CH_min
    ATMOS_PHY_SF_CE_min = CE_min
    ATMOS_PHY_SF_Z00    = Z00
    ATMOS_PHY_SF_Z0R    = Z0R
    ATMOS_PHY_SF_Z0S    = Z0S
    ATMOS_PHY_SF_Zt0    = Zt0
    ATMOS_PHY_SF_ZtR    = ZtR
    ATMOS_PHY_SF_ZtS    = ZtS
    ATMOS_PHY_SF_Ze0    = Ze0
    ATMOS_PHY_SF_ZeR    = ZeR
    ATMOS_PHY_SF_ZeS    = ZeS
    ATMOS_PHY_SF_ThS    = ThS
    ATMOS_PHY_SF_CM_const = Cm_const
    ATMOS_PHY_SF_Const_SH = Const_SH
    ATMOS_PHY_SF_Const_LH = Const_LH
    ATMOS_PHY_SF_Const_Ustar = Const_Ustar
    ATMOS_PHY_SF_Const_FREQ = Const_FREQ
    ATMOS_PHY_SF_FLG_MOM_FLUX = FLG_MOM_FLUX
    ATMOS_PHY_SF_FLG_SH_DIURNAL = FLG_SH_DIURNAL

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACE]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF)

    U_minM = ATMOS_PHY_SF_U_minM
    U_minH = ATMOS_PHY_SF_U_minH
    U_minE = ATMOS_PHY_SF_U_minE
    CM_min = ATMOS_PHY_SF_CM_min
    CH_min = ATMOS_PHY_SF_CH_min
    CE_min = ATMOS_PHY_SF_CE_min
    Z00    = ATMOS_PHY_SF_Z00
    Z0R    = ATMOS_PHY_SF_Z0R
    Z0S    = ATMOS_PHY_SF_Z0S
    Zt0    = ATMOS_PHY_SF_Zt0
    ZtR    = ATMOS_PHY_SF_ZtR
    ZtS    = ATMOS_PHY_SF_ZtS
    Ze0    = ATMOS_PHY_SF_Ze0
    ZeR    = ATMOS_PHY_SF_ZeR
    ZeS    = ATMOS_PHY_SF_ZeS
    ThS    = ATMOS_PHY_SF_ThS
    Cm_const = ATMOS_PHY_SF_CM_const
    Const_SH = ATMOS_PHY_SF_Const_SH
    Const_LH = ATMOS_PHY_SF_Const_LH
    Const_Ustar = ATMOS_PHY_SF_Const_Ustar
    Const_FREQ = ATMOS_PHY_SF_Const_FREQ
    FLG_MOM_FLUX = ATMOS_PHY_SF_FLG_MOM_FLUX
    FLG_SH_DIURNAL = ATMOS_PHY_SF_FLG_SH_DIURNAL

    return
  end subroutine ATMOS_PHY_SF_setup

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF
    use mod_const, only : &
       GRAV   => CONST_GRAV,   &
       KARMAN => CONST_KARMAN, &
       Rdry   => CONST_Rdry,   &
       CPdry  => CONST_CPdry,  &
       Rvap   => CONST_Rvap,   &
       RovCP  => CONST_RovCP,  &
       CPovCV => CONST_CPovCV, &
       P00    => CONST_PRE00,  &
       T00    => CONST_TEM00,  &
       LH0    => CONST_LH0,    &
       EPSvap => CONST_EPSvap, &
       PSAT0  => CONST_PSAT0
    use mod_time, only: &
       dttb => TIME_DTSEC_ATMOS_PHY_TB, &
       ctime => TIME_NOWSEC
    use mod_grid, only : &
       CZ  => GRID_CZ
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    use mod_ocean_vars, only: &
       SST
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    implicit none

    ! monitor
    real(RP) :: SHFLX(1,IA,JA) ! sensible heat flux [W/m2]
    real(RP) :: LHFLX(1,IA,JA) ! latent   heat flux [W/m2]

    ! work
    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Ustar ! friction velocity [m/s]
    real(RP) :: Cm    !

    integer :: i, j, iw
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface'

    do j = JS-1, JE
    do i = IS-1, IE

       ! at cell center

       !--- absolute velocity
       Uabs = sqrt( &
              ( MOMZ(KS,i,j)                  )**2 &
            + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
            + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 &
            ) / DENS(KS,i,j) * 0.5_RP

       !--- friction velocity at u, v, and w points
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
        Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
        Ustar = Const_Ustar
       endif

       ! flux
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
        SFLX_MOMZ(i,j) = - Cm_const * min(max(Uabs,U_minM),U_maxM) &
            * MOMZ(KS,i,j) * 0.5_RP
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
        Cm = Ustar*Ustar * Uabs
        SFLX_MOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) &
            * MOMZ(KS,i,j) * 0.5_RP
       endif
       if( FLG_SH_DIURNAL ) then
        SFLX_POTT(i,j) =  Const_SH / CPdry &
                       *  sin( ctime / ( Const_FREQ*3600.0_RP )*2.0_RP*3.1415926535_RP )
       else
        SFLX_POTT(i,j) =  Const_SH / CPdry
       endif
       SFLX_QV  (i,j) =  Const_LH / LH0


       ! at (u, y, layer)
       Uabs = sqrt( &
              ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i+1,j)                                     ) )**2 &
            + ( 2.0_RP *   MOMX(KS,i,j)                                                        )**2 &
            + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
        Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
        Ustar = Const_Ustar
       endif

       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
        Cm = Ustar*Ustar * Uabs
        SFLX_MOMX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) &
            * MOMX(KS,i,j) * 0.5_RP
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
        SFLX_MOMX(i,j) = - min(max(Cm,Cm_min),Cm_min) * min(max(Uabs,U_minM),U_maxM) &
            * MOMX(KS,i,j)
       endif


       ! at (x, v, layer)
       Uabs = sqrt( &
              ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i,j+1)                                     ) )**2 &
            + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(KS,i,j)                                                        )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
        Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
        Ustar = Const_Ustar
       endif

       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
        SFLX_MOMY(i,j) = - Cm_const * min(max(Uabs,U_minM),U_maxM) &
            * MOMY(KS,i,j) * 0.5_RP
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
        Cm = Ustar*Ustar * Uabs
        SFLX_MOMY(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) &
            * MOMY(KS,i,j)
       endif

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       SHFLX(1,i,j) = SFLX_POTT(i,j) * CPdry
       LHFLX(1,i,j) = SFLX_QV  (i,j) * LH0
    enddo
    enddo

    call HIST_in( SHFLX(:,:,:), 'SHFLX', 'sensible heat flux', 'W/m2', '2D', dttb )
    call HIST_in( LHFLX(:,:,:), 'LHFLX', 'latent heat flux',   'W/m2', '2D', dttb )

    return
  end subroutine ATMOS_PHY_SF

end module mod_atmos_phy_sf
