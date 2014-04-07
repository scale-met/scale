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
module scale_atmos_phy_sf_rico
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_SF_rico_setup
  public :: ATMOS_PHY_SF_rico

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !

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

  real(RP), private, save      :: U_minM  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxM  =  100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, save      :: U_minE  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxE  =  100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, save      :: U_minH  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH  =  100.0_RP  ! maximum U_abs for u,v,w

  real(RP), private, save      :: Cm_const =  0.0011_RP ! constant bulk coef. of u,v,w
  real(RP), private, save      :: Ce_const =  0.0011_RP ! constant bulk coef. of u,v,w
  real(RP), private, save      :: Ch_const =  0.0011_RP ! constant bulk coef. of u,v,w
  real(RP), private, save      :: Const_SH =  15.0_RP   ! constant surface sensible flux [W/m2]
  real(RP), private, save      :: Const_LH =  115.0_RP  ! constant surface latent flux [W/m2]
  real(RP), private, save      :: Const_Ustar = 0.25_RP ! constant friction velocity [m/s]

  integer(4), private, save    :: FLG_MOM_FLUX = 0      ! 0->Bulk coef. is constant
                                                        ! 1->Friction velocity is constant

  real(RP), private, save      :: Const_FREQ = 24.0_RP  ! frequency of sensible heat flux [hour]
  !  SHFLX = Const_SH[W/m^2] * sin( 2*pi*(current time)/Const_FREQ )
  logical, private, save       :: FLG_SH_DIURNAL = .false.
 
  real(RP), private, save      :: FIXED_PTSST = 298.5_RP! fixed PT of surface [K]
  real(RP), private, save      :: FIXED_SST = 299.8_RP  ! fixed Temperature of surface [K]
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_rico_setup( ATMOS_TYPE_PHY_SF )
    use scale_process, only: &
       PRC_MPIstop
    implicit none
 
    character(len=H_SHORT), intent(in) :: ATMOS_TYPE_PHY_SF

    real(RP) :: ATMOS_PHY_SF_U_minM ! minimum U_abs for u,v,w
    real(RP) :: ATMOS_PHY_SF_CM_min ! minimum bulk coef. of u,v,w
    real(RP) :: ATMOS_PHY_SF_CM_const
    real(RP) :: ATMOS_PHY_SF_CH_const
    real(RP) :: ATMOS_PHY_SF_CE_const
    real(RP) :: ATMOS_PHY_SF_Const_SH
    real(RP) :: ATMOS_PHY_SF_Const_LH
    real(RP) :: ATMOS_PHY_SF_Const_Ustar
    real(RP) :: ATMOS_PHY_SF_Const_FREQ
    integer  :: ATMOS_PHY_SF_FLG_MOM_FLUX
    logical  :: ATMOS_PHY_SF_FLG_SH_DIURNAL
    real(RP) :: ATMOS_PHY_SF_FIXED_PTSST

    NAMELIST / PARAM_ATMOS_PHY_SF_RICO / &
       ATMOS_PHY_SF_U_minM, &
       ATMOS_PHY_SF_CM_min, &
       ATMOS_PHY_SF_CM_const, &
       ATMOS_PHY_SF_CE_const, &
       ATMOS_PHY_SF_CH_const, &
       ATMOS_PHY_SF_Const_SH, &
       ATMOS_PHY_SF_Const_LH, &
       ATMOS_PHY_SF_Const_Ustar, &
       ATMOS_PHY_SF_Const_FREQ, &
       ATMOS_PHY_SF_FLG_MOM_FLUX, &
       ATMOS_PHY_SF_FLG_SH_DIURNAL, &
       ATMOS_PHY_SF_FIXED_PTSST

    integer :: ierr
    !---------------------------------------------------------------------------

    ATMOS_PHY_SF_U_minM = U_minM
    ATMOS_PHY_SF_CM_min = CM_min
    ATMOS_PHY_SF_CM_const = Cm_const
    ATMOS_PHY_SF_CH_const = Ch_const
    ATMOS_PHY_SF_CE_const = Ce_const
    ATMOS_PHY_SF_Const_SH = Const_SH
    ATMOS_PHY_SF_Const_LH = Const_LH
    ATMOS_PHY_SF_Const_Ustar = Const_Ustar
    ATMOS_PHY_SF_Const_FREQ = Const_FREQ
    ATMOS_PHY_SF_FLG_MOM_FLUX = FLG_MOM_FLUX
    ATMOS_PHY_SF_FLG_SH_DIURNAL = FLG_SH_DIURNAL
    ATMOS_PHY_SF_FIXED_PTSST = FIXED_PTSST

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[PHY_SURFACE]/Categ[ATMOS]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_SF_RICO,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_SF_RICO. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_SF_RICO)

    U_minM = ATMOS_PHY_SF_U_minM
    CM_min = ATMOS_PHY_SF_CM_min
    Cm_const = ATMOS_PHY_SF_CM_const
    Ce_const = ATMOS_PHY_SF_CE_const
    Ch_const = ATMOS_PHY_SF_CH_const
    Const_SH = ATMOS_PHY_SF_Const_SH
    Const_LH = ATMOS_PHY_SF_Const_LH
    Const_Ustar = ATMOS_PHY_SF_Const_Ustar
    Const_FREQ = ATMOS_PHY_SF_Const_FREQ
    FLG_MOM_FLUX = ATMOS_PHY_SF_FLG_MOM_FLUX
    FLG_SH_DIURNAL = ATMOS_PHY_SF_FLG_SH_DIURNAL
    FIXED_PTSST = ATMOS_PHY_SF_FIXED_PTSST 

    return
  end subroutine ATMOS_PHY_SF_rico_setup
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_SF_rico( &
         SFLX_MOMZ, SFLX_MOMX, SFLX_MOMY, SFLX_POTT, SFLX_QV, & ! (out)
         DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, SST,             & ! (in)
         CZ, ctime                                            ) ! (in)
    use dc_types, only: &
       DP
    use scale_const, only : &
       CPdry  => CONST_CPdry,  &
       RovCP  => CONST_RovCP,  &
       P00    => CONST_PRE00,  &
       T00    => CONST_TEM00,  &
       EPSvap => CONST_EPSvap, &
       Rvap   => CONST_Rvap,  &
       Rdry   => CONST_Rdry,  &
       PSAT0  => CONST_PSAT0, &
       LH0    => CONST_LH0
    use scale_atmos_thermodyn, only: &
       CPw => AQ_CP
    use scale_atmos_saturation, only : &
       moist_pres2qsat_liq  => ATMOS_SATURATION_pres2qsat_liq
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

    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: SST (1,IA,JA)

    real(RP), intent(in)  :: CZ(KA)

    real(DP), intent(in)  :: ctime

    ! work
    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce    !
    real(RP) :: qv_evap, pres_evap, pres
    real(RP) :: Rtot, qdry, qtot, PT_SST, CPtot, CPovCV
    real(RP), parameter :: pres_sfc = 1015.4E2_RP

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

       !--- Bulk coef. at w, theta, and qv points
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
          Cm = Cm_const
          Ch = Ch_const
          Ce = Ce_const
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
          Cm = min( max(Const_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       endif

       !--- saturation at surface
       qtot = 0.0_RP
       qdry = 1.0_RP
       CPtot = 0.0_RP
       do iw = QQS, QQE
          qdry = qdry - QTRC(KS,i,j,iw)
          qtot = qtot + QTRC(KS,i,j,iw)
          CPtot = CPtot + QTRC(KS,i,j,iw) * CPw(iw)
       enddo
       Rtot = Rdry*qdry + Rvap*QTRC(KS,i,j,I_QV)
       CPtot = CPdry*qdry + CPtot
       CPovCV = CPtot / ( CPtot - Rtot )
       pres      = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV
!       pres_evap = PSAT0 * exp( LH0/Rvap * ( 1.0_RP/T00 - 1.0_RP/SST(1,i,j) ) )
!       qv_evap   = EPSvap * pres_evap / pres_sfc 
!       call moist_pres2qsat_liq( qv_evap, SST(1,i,j), pres )
!       call moist_pres2qsat_liq( qv_evap, SST(1,i,j), pres_sfc )
       call moist_pres2qsat_liq( qv_evap, FIXED_SST, pres_sfc )
       ! flux
       SFLX_MOMZ(i,j) = 0.0_RP
       SFLX_POTT(i,j) =   Ch * min(max(Uabs,U_minH),U_maxH) &
            * ( FIXED_PTSST * DENS(KS,i,j) - RHOT(KS,i,j) )
       SFLX_QV  (i,j) =   Ce * min(max(Uabs,U_minE),U_maxE) &
            * DENS(KS,i,j) * ( qv_evap - qtot )

       ! at (u, y, layer)
       Uabs = sqrt( &
            + ( 2.0_RP *   MOMX(KS,i,j)                                                        )**2 &
            + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
          Cm = Cm_const
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
          Cm = min( max(Const_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       endif

       SFLX_MOMX(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMX(KS,i,j)


       ! at (x, v, layer)
       Uabs = sqrt( &
            + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(KS,i,j)                                                        )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
       if( FLG_MOM_FLUX == 0  ) then     ! Bulk coef. is constant
          Cm = Cm_const
       elseif( FLG_MOM_FLUX == 1  ) then ! friction velocity is constant
          Cm = min( max(Const_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       endif

       SFLX_MOMY(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMY(KS,i,j)

    enddo
    enddo

    return
  end subroutine ATMOS_PHY_SF_rico

end module scale_atmos_phy_sf_rico
