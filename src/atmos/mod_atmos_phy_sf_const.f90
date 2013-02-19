!-------------------------------------------------------------------------------
!> module Atmosphere / Physics Surface fluxes
!!
!! @par Description
!!          Flux from/to bottom wall of atmosphere (surface)
!!          Bulk Method
!!
!! @author Team SCALE
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
  include 'inc_precision.h'
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

  ! limiter
  real(RP), private, save      :: Cm_min  =    1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, parameter :: Cm_max  =    2.5E-3_RP ! maximum bulk coef. of u,v,w

  real(RP), private, save      :: U_minM  =    0.0_RP   ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxM  =  100.0_RP   ! maximum U_abs for u,v,w

  real(RP), private, save      :: Const_Cm =  0.0011_RP ! constant bulk coef. of u,v,w
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
    use mod_atmos_vars, only: &
       ATMOS_TYPE_PHY_SF
    use mod_time, only: &
       NOWSEC => TIME_NOWDAYSEC
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP) :: ATMOS_PHY_SF_U_minM ! minimum U_abs for u,v,w
    real(RP) :: ATMOS_PHY_SF_CM_min ! minimum bulk coef. of u,v,w
    real(RP) :: ATMOS_PHY_SF_Const_CM
    real(RP) :: ATMOS_PHY_SF_Const_SH
    real(RP) :: ATMOS_PHY_SF_Const_LH
    real(RP) :: ATMOS_PHY_SF_Const_Ustar
    real(RP) :: ATMOS_PHY_SF_Const_FREQ
    integer  :: ATMOS_PHY_SF_FLG_MOM_FLUX
    logical  :: ATMOS_PHY_SF_FLG_SH_DIURNAL

    NAMELIST / PARAM_ATMOS_PHY_SF_CONST / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_CM_min, &
       ATMOS_PHY_SF_Const_CM, &
       ATMOS_PHY_SF_Const_SH, &
       ATMOS_PHY_SF_Const_LH, &
       ATMOS_PHY_SF_Const_Ustar, &
       ATMOS_PHY_SF_Const_FREQ, &
       ATMOS_PHY_SF_FLG_MOM_FLUX, &
       ATMOS_PHY_SF_FLG_SH_DIURNAL

    integer :: ierr
    !---------------------------------------------------------------------------

    ATMOS_PHY_SF_U_minM = U_minM
    ATMOS_PHY_SF_CM_min = CM_min
    ATMOS_PHY_SF_Const_CM = Const_Cm
    ATMOS_PHY_SF_Const_SH = Const_SH
    ATMOS_PHY_SF_Const_LH = Const_LH
    ATMOS_PHY_SF_Const_Ustar = Const_Ustar
    ATMOS_PHY_SF_Const_FREQ = Const_FREQ
    ATMOS_PHY_SF_FLG_MOM_FLUX = FLG_MOM_FLUX
    ATMOS_PHY_SF_FLG_SH_DIURNAL = FLG_SH_DIURNAL

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACEFLUX]/Categ[ATMOS]'
    if( IO_L ) write(IO_FID_LOG,*) '*** Constant flux parameter'

    if ( ATMOS_TYPE_PHY_SF /= 'CONST' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_TYPE_PHY_SF is not CONST. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_CONST,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_CONST. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_CONST)

    U_minM = ATMOS_PHY_SF_U_minM
    CM_min = ATMOS_PHY_SF_CM_min
    Const_Cm = ATMOS_PHY_SF_Const_Cm
    Const_SH = ATMOS_PHY_SF_Const_SH
    Const_LH = ATMOS_PHY_SF_Const_LH
    Const_Ustar = ATMOS_PHY_SF_Const_Ustar
    Const_FREQ = ATMOS_PHY_SF_Const_FREQ
    FLG_MOM_FLUX = ATMOS_PHY_SF_FLG_MOM_FLUX
    FLG_SH_DIURNAL = ATMOS_PHY_SF_FLG_SH_DIURNAL

    call ATMOS_PHY_SF_main( &
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
         DENS, MOMZ, MOMX, MOMY,                              & ! (in)
         NOWSEC                                               ) ! (in)

    call COMM_vars8( SFLX_MOMZ(:,:), 1 )
    call COMM_vars8( SFLX_MOMX(:,:), 2 )
    call COMM_vars8( SFLX_MOMY(:,:), 3 )
    call COMM_vars8( SFLX_POTT(:,:), 4 )
    call COMM_vars8( SFLX_QV  (:,:), 5 )

    call COMM_wait ( SFLX_MOMZ(:,:), 1 )
    call COMM_wait ( SFLX_MOMX(:,:), 2 )
    call COMM_wait ( SFLX_MOMY(:,:), 3 )
    call COMM_wait ( SFLX_POTT(:,:), 4 )
    call COMM_wait ( SFLX_QV  (:,:), 5 )

    return
  end subroutine ATMOS_PHY_SF_setup

  !-----------------------------------------------------------------------------
  ! calculation flux
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF
    use mod_time, only: &
       dtsf => TIME_DTSEC_ATMOS_PHY_SF, &
       NOWSEC => TIME_NOWDAYSEC
    use mod_const, only: &
       CPdry  => CONST_CPdry,  &
       LH0    => CONST_LH0
    use mod_comm, only: &
       COMM_vars8, &
       COMM_wait
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY
    use mod_atmos_vars_sf, only: &
       SFLX_MOMZ, &
       SFLX_MOMX, &
       SFLX_MOMY, &
       SFLX_POTT, &
       SFLX_QV
    implicit none

    ! monitor
    real(RP) :: SHFLX(IA,JA) ! sensible heat flux [W/m2]
    real(RP) :: LHFLX(IA,JA) ! latent   heat flux [W/m2]

    integer :: i, j

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux'

    call ATMOS_PHY_SF_main( &
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
         DENS, MOMZ, MOMX, MOMY,                              & ! (in)
         NOWSEC                                               ) ! (out)

    call COMM_vars8( SFLX_MOMZ(:,:), 1 )
    call COMM_vars8( SFLX_MOMX(:,:), 2 )
    call COMM_vars8( SFLX_MOMY(:,:), 3 )
    call COMM_vars8( SFLX_POTT(:,:), 4 )
    call COMM_vars8( SFLX_QV  (:,:), 5 )

    do j = JS, JE
    do i = IS, IE
       SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
       LHFLX(i,j) = SFLX_QV  (i,j) * LH0
    enddo
    enddo

    call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2', dtsf )
    call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2', dtsf )

    call COMM_wait ( SFLX_MOMZ(:,:), 1 )
    call COMM_wait ( SFLX_MOMX(:,:), 2 )
    call COMM_wait ( SFLX_MOMY(:,:), 3 )
    call COMM_wait ( SFLX_POTT(:,:), 4 )
    call COMM_wait ( SFLX_QV  (:,:), 5 )

    return
  end subroutine ATMOS_PHY_SF
  !-----------------------------------------------------------------------------
  ! calculation flux
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_main( &
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
         DENS, MOMZ, MOMX, MOMY,                              & ! (in)
         ctime                                                ) ! (in)
    use mod_const, only: &
       CPdry  => CONST_CPdry,  &
       LH0    => CONST_LH0
    use dc_types, only: &
         DP
    implicit none

    real(RP), intent(out) :: SFLX_MOMZ(IA,JA)
    real(RP), intent(out) :: SFLX_MOMX(IA,JA)
    real(RP), intent(out) :: SFLX_MOMY(IA,JA)
    real(RP), intent(out) :: SFLX_POTT(IA,JA)
    real(RP), intent(out) :: SFLX_QV  (IA,JA)

    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)

    ! work
    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm    !

    real(DP) :: ctime

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = JS, JE
    do i = IS, IE

       ! at cell center

       !--- absolute velocity
       Uabs = sqrt( &
              ( MOMZ(KS,i,j)                  )**2 &
            + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
            + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 &
            ) / DENS(KS,i,j) * 0.5_RP

       !--- Bulk coef. at w, theta, and qv points
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
          Cm = Const_Cm
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
          Cm = min( max(Const_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       endif

       ! flux
       SFLX_MOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) &
            * MOMZ(KS,i,j) * 0.5_RP

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
          Cm = Const_Cm
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
          Cm = min( max(Const_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       endif

       SFLX_MOMX(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMX(KS,i,j)


       ! at (x, v, layer)
       Uabs = sqrt( &
              ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i,j+1)                                     ) )**2 &
            + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(KS,i,j)                                                        )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
          Cm = Const_Cm
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
          Cm = min( max(Const_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       endif

       SFLX_MOMY(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMY(KS,i,j)

    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_main

end module mod_atmos_phy_sf
