!-------------------------------------------------------------------------------
!> module COUPLER / Atmosphere-Land Surface fluxes
!!
!! @par Description
!!          Surface flux between atmosphere and land with Bulk Method
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-08-31 (T.Yamaura)  [new]
!<
!-------------------------------------------------------------------------------
module mod_cpl_atmos_land
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
  public :: CPL_AtmLnd_setup
  public :: CPL_AtmLnd_solve
  public :: CPL_AtmLnd_unsolve
  public :: CPL_AtmLnd_putAtm
  public :: CPL_AtmLnd_putLnd
  public :: CPL_AtmLnd_getDat2Atm
  public :: CPL_AtmLnd_getDat2Lnd
  public :: CPL_AtmLnd_flushDat2Atm
  public :: CPL_AtmLnd_flushDat2Lnd

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
  private :: bulkcoef_uno
  private :: satmixr
  private :: satvapor

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! surface fluxes for atmosphere
  real(RP), private, save :: SFLX_MOMX (IA,JA) ! momentum flux for x [kg/m2/s]
  real(RP), private, save :: SFLX_MOMY (IA,JA) ! momentum flux for y [kg/m2/s]
  real(RP), private, save :: SFLX_MOMZ (IA,JA) ! momentum flux for z [kg/m2/s]
  real(RP), private, save :: SFLX_SWU  (IA,JA) ! upward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, save :: SFLX_LWU  (IA,JA) ! upward long-wave radiation flux (upward positive) [W/m2]
  real(RP), private, save :: SFLX_SH   (IA,JA) ! sensible heat flux (upward positive) [W/m2]
  real(RP), private, save :: SFLX_LH   (IA,JA) ! latent heat flux (upward positive) [W/m2]
  real(RP), private, save :: SFLX_QVAtm(IA,JA) ! moisture flux for atmosphere [kg/m2/s]
  ! surface fluxes for land
  real(RP), private, save :: SFLX_GH   (IA,JA) ! ground heat flux (upward positive) [W/m2]
  real(RP), private, save :: SFLX_PREC (IA,JA) ! precipitation flux [kg/m2/s]
  real(RP), private, save :: SFLX_QVLnd(IA,JA) ! moisture flux for land [kg/m2/s]

  ! Atmospheric values
  real(RP), private, save :: DENS(KA,IA,JA)    ! air density [kg/m3]
  real(RP), private, save :: MOMX(KA,IA,JA)    ! momentum x [kg/m2/s]
  real(RP), private, save :: MOMY(KA,IA,JA)    ! momentum y [kg/m2/s]
  real(RP), private, save :: MOMZ(KA,IA,JA)    ! momentum z [kg/m2/s]
  real(RP), private, save :: RHOT(KA,IA,JA)    ! rho * theta [K*kg/m3]
  real(RP), private, save :: QTRC(KA,IA,JA,QA) ! ratio of mass of tracer to total mass [kg/kg]
  real(RP), private, save :: PREC(IA,JA)       ! surface precipitation rate [kg/m2/s]
  real(RP), private, save :: SWD (IA,JA)       ! downward short-wave radiation flux (upward positive) [W/m2]
  real(RP), private, save :: LWD (IA,JA)       ! downward long-wave radiation flux (upward positive) [W/m2]

  ! Land values
  real(RP), private, save :: TG   (IA,JA) ! soil temperature [K]
  real(RP), private, save :: QvEfc(IA,JA) ! efficiency of evaporation [no unit]
  real(RP), private, save :: EMIT (IA,JA) ! emissivity in long-wave radiation [no unit]
  real(RP), private, save :: ALB  (IA,JA) ! surface albedo in short-wave radiation [no unit]
  real(RP), private, save :: TCS  (IA,JA) ! thermal conductivity for soil [W/m/K]
  real(RP), private, save :: DZg  (IA,JA) ! soil depth [m]
  real(RP), private, save :: Z00  (IA,JA) ! basic factor for momemtum
  real(RP), private, save :: Z0R  (IA,JA) ! rough factor for momemtum
  real(RP), private, save :: Z0S  (IA,JA) ! smooth factor for momemtum
  real(RP), private, save :: Zt0  (IA,JA) ! basic factor for heat
  real(RP), private, save :: ZtR  (IA,JA) ! rough factor for heat
  real(RP), private, save :: ZtS  (IA,JA) ! smooth factor for heat
  real(RP), private, save :: Ze0  (IA,JA) ! basic factor for moisture
  real(RP), private, save :: ZeR  (IA,JA) ! rough factor for moisture
  real(RP), private, save :: ZeS  (IA,JA) ! smooth factor for moisture

  ! counter
  real(RP), private, save :: CNT_putAtm     ! counter for putAtm
  real(RP), private, save :: CNT_putLnd     ! counter for putLnd
  real(RP), private, save :: CNT_getDat2Atm ! counter for getDat2Atm
  real(RP), private, save :: CNT_getDat2Lnd ! counter for getDat2Lnd

  ! limiter
  real(RP), private, parameter :: res_min   =  1.0E-10_RP ! minimum number of residual in the Newton-Raphson Method
  real(RP), private, parameter :: Ustar_min =  1.0E-3_RP ! minimum limit of U*

  real(RP), private, parameter :: Z0_min =   1.0E-5_RP ! minimum roughness length of u,v,w
  real(RP), private, parameter :: Zt_min =   1.0E-5_RP !                             T
  real(RP), private, parameter :: Ze_min =   1.0E-5_RP !                             q

  real(RP), private, parameter :: Cm_min =   1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_min =   1.0E-5_RP !                       T
  real(RP), private, parameter :: Ce_min =   1.0E-5_RP !                       q
  real(RP), private, parameter :: Cm_max =   2.5E-3_RP ! maximum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_max =   1.0_RP    !                       T
  real(RP), private, parameter :: Ce_max =   1.0_RP    !                       q

  real(RP), private, parameter :: U_minM =   0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_minH =   0.0_RP  !                   T
  real(RP), private, parameter :: U_minE =   0.0_RP  !                   q
  real(RP), private, parameter :: U_maxM = 100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH = 100.0_RP  !                   T
  real(RP), private, parameter :: U_maxE = 100.0_RP  !                   q

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmLnd_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_vars, only: &
       CPL_RESTART_IN_BASENAME
    implicit none

    logical  :: dummy

    NAMELIST / PARAM_CPL_AtmLnd / &
       dummy

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[AtmLnd]/Categ[CPL]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CPL_AtmLnd,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_CPL_AtmLnd. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_CPL_AtmLnd)

    CNT_putAtm     = 0.0_RP
    CNT_putLnd     = 0.0_RP

    call CPL_AtmLnd_flushDat2Atm
    call CPL_AtmLnd_flushDat2Lnd

    return
  end subroutine CPL_AtmLnd_setup

  subroutine CPL_AtmLnd_solve
    use mod_const, only: &
       LH0 => CONST_LH0
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_vars, only: &
       LST
    implicit none

    ! parameters
    integer,  parameter :: nmax     = 1000      ! maximum iteration number
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)

    ! works
    integer :: i, j, n

    real(RP) :: RES  (IA,JA)
    real(RP) :: DRES (IA,JA)
    real(RP) :: pMOMX(IA,JA)
    real(RP) :: pMOMY(IA,JA)
    real(RP) :: pMOMZ(IA,JA)
    real(RP) :: pSWU (IA,JA)
    real(RP) :: pLWU (IA,JA)
    real(RP) :: pSH  (IA,JA)
    real(RP) :: pLH  (IA,JA)
    real(RP) :: pGH  (IA,JA)

    real(RP) :: oldRES (IA,JA) ! RES in previous step
    real(RP) :: redf ! reduced factor

    if( IO_L ) write(IO_FID_LOG,*) '*** CPL solve: Atmos-Land'

    ! average putAtm
    if( int( CNT_putAtm ) == 0 ) then
      if( IO_L ) write(IO_FID_LOG,*) 'Error: divided by zero (CNT_putAtm)'
      call PRC_MPIstop
    end if

    DENS(:,:,:)   = DENS(:,:,:)   / CNT_putAtm
    MOMX(:,:,:)   = MOMX(:,:,:)   / CNT_putAtm
    MOMY(:,:,:)   = MOMY(:,:,:)   / CNT_putAtm
    MOMZ(:,:,:)   = MOMZ(:,:,:)   / CNT_putAtm
    RHOT(:,:,:)   = RHOT(:,:,:)   / CNT_putAtm
    QTRC(:,:,:,:) = QTRC(:,:,:,:) / CNT_putAtm
    PREC(:,:)     = PREC(:,:)     / CNT_putAtm
    SWD(:,:)      = SWD(:,:)      / CNT_putAtm
    LWD(:,:)      = LWD(:,:)      / CNT_putAtm

    ! average putLnd
    if( int( CNT_putLnd ) == 0 ) then
      if( IO_L ) write(IO_FID_LOG,*) 'Error: divided by zero (CNT_putLnd)'
      call PRC_MPIstop
    end if

    TG   (:,:)    = TG   (:,:)    / CNT_putLnd
    QvEfc(:,:)    = QvEfc(:,:)    / CNT_putLnd
    EMIT (:,:)    = EMIT (:,:)    / CNT_putLnd
    ALB  (:,:)    = ALB  (:,:)    / CNT_putLnd
    TCS  (:,:)    = TCS  (:,:)    / CNT_putLnd
    DZg  (:,:)    = DZg  (:,:)    / CNT_putLnd
    Z00  (:,:)    = Z00  (:,:)    / CNT_putLnd
    Z0R  (:,:)    = Z0R  (:,:)    / CNT_putLnd
    Z0S  (:,:)    = Z0S  (:,:)    / CNT_putLnd
    Zt0  (:,:)    = Zt0  (:,:)    / CNT_putLnd
    ZtR  (:,:)    = ZtR  (:,:)    / CNT_putLnd
    ZtS  (:,:)    = ZtS  (:,:)    / CNT_putLnd
    Ze0  (:,:)    = Ze0  (:,:)    / CNT_putLnd
    ZeR  (:,:)    = ZeR  (:,:)    / CNT_putLnd
    ZeS  (:,:)    = ZeS  (:,:)    / CNT_putLnd

    redf        = 1.0_RP
    oldRES(:,:) = 1.0E+5_RP

    do n = 1, nmax

      call ts_residual( &
        RES, DRES,           & ! (out)
        pMOMX, pMOMY, pMOMZ, & ! (out)
        pSWU, pLWU,          & ! (out)
        pSH, pLH, pGH        ) ! (out)

      do j = JS, JE+1
      do i = IS, IE+1

        if( redf < 0.0_RP ) then
          redf = 1.0_RP
        end if

        if( abs(RES(i,j)) > abs(oldRES(i,j)) ) then
          redf = max( TFa*redf, redf_min )
        else
          redf = min( TFb*redf, redf_max )
        end if

        if( DRES(i,j) > 0.0_RP ) then
          redf = -1.0_RP
        end if

        ! update surface temperature
        LST(i,j) = LST(i,j) - redf * RES(i,j)/DRES(i,j)

        ! save residual in this step
        oldRES(i,j) = RES(i,j)

      end do
      end do

      if( maxval(abs(RES(IS:IE,JS:JE))) < res_min ) then
        ! iteration converged
        exit
      end if

    end do

    if( n > nmax ) then
      ! not converged and stop program
      if( IO_L ) write(IO_FID_LOG,*) 'Error: surface tempearture is not converged.'
      if( IO_L ) write(IO_FID_LOG,*) 'LST  :',minval(LST  (:,:)),maxval(LST  (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'RES  :',minval(RES  (:,:)),maxval(RES  (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'DRES :',minval(DRES (:,:)),maxval(DRES (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pMOMX:',minval(pMOMX(:,:)),maxval(pMOMX(:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pMOMY:',minval(pMOMY(:,:)),maxval(pMOMY(:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pMOMZ:',minval(pMOMZ(:,:)),maxval(pMOMZ(:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pSWU :',minval(pSWU (:,:)),maxval(pSWU (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pLWU :',minval(pLWU (:,:)),maxval(pLWU (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pSH  :',minval(pSH  (:,:)),maxval(pSH  (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pLH  :',minval(pLH  (:,:)),maxval(pLH  (:,:))
      if( IO_L ) write(IO_FID_LOG,*) 'pGH  :',minval(pGH  (:,:)),maxval(pGH  (:,:))
      call PRC_MPIstop
    end if

    ! put residual in ground heat flux
    pGH(:,:) = pGH(:,:) - RES(:,:)

    ! save flux
    SFLX_MOMX (:,:) = SFLX_MOMX (:,:) + pMOMX(:,:)
    SFLX_MOMY (:,:) = SFLX_MOMY (:,:) + pMOMY(:,:)
    SFLX_MOMZ (:,:) = SFLX_MOMZ (:,:) + pMOMZ(:,:)
    SFLX_SWU  (:,:) = SFLX_SWU  (:,:) + pSWU (:,:)
    SFLX_LWU  (:,:) = SFLX_LWU  (:,:) + pLWU (:,:)
    SFLX_SH   (:,:) = SFLX_SH   (:,:) + pSH  (:,:)
    SFLX_LH   (:,:) = SFLX_LH   (:,:) + pLH  (:,:)
    SFLX_QVAtm(:,:) = SFLX_QVAtm(:,:) + pLH  (:,:)/LH0

    SFLX_GH   (:,:) = SFLX_GH   (:,:) + pGH  (:,:)
    SFLX_PREC (:,:) = SFLX_PREC (:,:) + PREC (:,:)
    SFLX_QVLnd(:,:) = SFLX_QVLnd(:,:) + pLH  (:,:)/LH0

    ! counter
    CNT_putAtm     = 0.0_RP
    CNT_putLnd     = 0.0_RP
    CNT_getDat2Atm = CNT_getDat2Atm + 1.0_RP
    CNT_getDat2Lnd = CNT_getDat2Lnd + 1.0_RP

    return
  end subroutine CPL_AtmLnd_solve

  subroutine CPL_AtmLnd_unsolve
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '*** CPL unsolve: Atmos-Land'

    return
  end subroutine CPL_AtmLnd_unsolve

  subroutine CPL_AtmLnd_putAtm( &
        pDENS, pMOMX, pMOMY, pMOMZ, pRHOT, & ! (in)
        pQTRC, pPREC, pSWD, pLWD     ) ! (in)
    implicit none

    real(RP), intent(in) :: pDENS  (KA,IA,JA)
    real(RP), intent(in) :: pMOMX  (KA,IA,JA)
    real(RP), intent(in) :: pMOMY  (KA,IA,JA)
    real(RP), intent(in) :: pMOMZ  (KA,IA,JA)
    real(RP), intent(in) :: pRHOT  (KA,IA,JA)
    real(RP), intent(in) :: pQTRC  (KA,IA,JA,QA)
    real(RP), intent(in) :: pPREC  (IA,JA)
    real(RP), intent(in) :: pSWD(IA,JA)
    real(RP), intent(in) :: pLWD(IA,JA)

    DENS(:,:,:)   = DENS(:,:,:)   + pDENS  (:,:,:)
    MOMX(:,:,:)   = MOMX(:,:,:)   + pMOMX  (:,:,:)
    MOMY(:,:,:)   = MOMY(:,:,:)   + pMOMY  (:,:,:)
    MOMZ(:,:,:)   = MOMZ(:,:,:)   + pMOMZ  (:,:,:)
    RHOT(:,:,:)   = RHOT(:,:,:)   + pRHOT  (:,:,:)
    QTRC(:,:,:,:) = QTRC(:,:,:,:) + pQTRC  (:,:,:,:)
    PREC(:,:)     = PREC(:,:)     + pPREC  (:,:)
    SWD (:,:)     = SWD (:,:)     + pSWD(:,:)
    LWD (:,:)     = LWD (:,:)     + pLWD(:,:)

    CNT_putAtm = CNT_putAtm + 1.0_RP
    
    return
  end subroutine CPL_AtmLnd_putAtm

  subroutine CPL_AtmLnd_putLnd( &
      pTG, pQvEfc, pEMIT, pALB, pTCS, pDZg,                & ! (in)
      pZ00, pZ0R, pZ0S, pZt0, pZtR, pZtS, pZe0, pZeR, pZeS ) ! (in)
    implicit none

    real(RP), intent(in) :: pTG   (IA,JA)
    real(RP), intent(in) :: pQvEfc(IA,JA)
    real(RP), intent(in) :: pEMIT (IA,JA)
    real(RP), intent(in) :: pALB  (IA,JA)
    real(RP), intent(in) :: pTCS  (IA,JA)
    real(RP), intent(in) :: pDZg  (IA,JA)
    real(RP), intent(in) :: pZ00  (IA,JA)
    real(RP), intent(in) :: pZ0R  (IA,JA)
    real(RP), intent(in) :: pZ0S  (IA,JA)
    real(RP), intent(in) :: pZt0  (IA,JA)
    real(RP), intent(in) :: pZtR  (IA,JA)
    real(RP), intent(in) :: pZtS  (IA,JA)
    real(RP), intent(in) :: pZe0  (IA,JA)
    real(RP), intent(in) :: pZeR  (IA,JA)
    real(RP), intent(in) :: pZeS  (IA,JA)

    TG   (:,:) = TG   (:,:) + pTG   (:,:)
    QvEfc(:,:) = QvEfc(:,:) + pQvEfc(:,:)
    EMIT (:,:) = EMIT (:,:) + pEMIT (:,:)
    ALB  (:,:) = ALB  (:,:) + pALB  (:,:)
    TCS  (:,:) = TCS  (:,:) + pTCS  (:,:)
    DZg  (:,:) = DZg  (:,:) + pDZg  (:,:)
    Z00  (:,:) = Z00  (:,:) + pZ00  (:,:)
    Z0R  (:,:) = Z0R  (:,:) + pZ0R  (:,:)
    Z0S  (:,:) = Z0S  (:,:) + pZ0S  (:,:)
    Zt0  (:,:) = Zt0  (:,:) + pZt0  (:,:)
    ZtR  (:,:) = ZtR  (:,:) + pZtR  (:,:)
    ZtS  (:,:) = ZtS  (:,:) + pZtS  (:,:)
    Ze0  (:,:) = Ze0  (:,:) + pZe0  (:,:)
    ZeR  (:,:) = ZeR  (:,:) + pZeR  (:,:)
    ZeS  (:,:) = ZeS  (:,:) + pZeS  (:,:)

    CNT_putLnd = CNT_putLnd + 1.0_RP

    return
  end subroutine CPL_AtmLnd_putLnd

  subroutine CPL_AtmLnd_getDat2Atm( &
      pSFLX_MOMX, pSFLX_MOMY, pSFLX_MOMZ, pSFLX_SWU, pSFLX_LWU, & ! (out)
      pSFLX_SH, pSFLX_LH, pSFLX_QVAtm                           ) ! (out)
    implicit none

    real(RP), intent(out) :: pSFLX_MOMX (IA,JA)
    real(RP), intent(out) :: pSFLX_MOMY (IA,JA)
    real(RP), intent(out) :: pSFLX_MOMZ (IA,JA)
    real(RP), intent(out) :: pSFLX_SWU  (IA,JA)
    real(RP), intent(out) :: pSFLX_LWU  (IA,JA)
    real(RP), intent(out) :: pSFLX_SH   (IA,JA)
    real(RP), intent(out) :: pSFLX_LH   (IA,JA)
    real(RP), intent(out) :: pSFLX_QVAtm(IA,JA)

    pSFLX_MOMX (:,:) = SFLX_MOMX (:,:) / CNT_getDat2Atm
    pSFLX_MOMY (:,:) = SFLX_MOMY (:,:) / CNT_getDat2Atm
    pSFLX_MOMZ (:,:) = SFLX_MOMZ (:,:) / CNT_getDat2Atm
    pSFLX_SWU  (:,:) = SFLX_SWU  (:,:) / CNT_getDat2Atm
    pSFLX_LWU  (:,:) = SFLX_LWU  (:,:) / CNT_getDat2Atm
    pSFLX_SH   (:,:) = SFLX_SH   (:,:) / CNT_getDat2Atm
    pSFLX_LH   (:,:) = SFLX_LH   (:,:) / CNT_getDat2Atm
    pSFLX_QVAtm(:,:) = SFLX_QVAtm(:,:) / CNT_getDat2Atm

    return
  end subroutine CPL_AtmLnd_getDat2Atm

  subroutine CPL_AtmLnd_getDat2Lnd( &
      pSFLX_GH, pSFLX_PREC, pSFLX_QVLnd ) ! (out)
    implicit none

    real(RP), intent(out) :: pSFLX_GH   (IA,JA)
    real(RP), intent(out) :: pSFLX_PREC (IA,JA)
    real(RP), intent(out) :: pSFLX_QVLnd(IA,JA)

    pSFLX_GH   (:,:) = SFLX_GH   (:,:) / CNT_getDat2Lnd
    pSFLX_PREC (:,:) = SFLX_PREC (:,:) / CNT_getDat2Lnd
    pSFLX_QVLnd(:,:) = SFLX_QVLnd(:,:) / CNT_getDat2Lnd

    return
  end subroutine CPL_AtmLnd_getDat2Lnd

  subroutine CPL_AtmLnd_flushDat2Atm
    implicit none

    DENS(:,:,:)    = 0.0_RP
    MOMX(:,:,:)    = 0.0_RP
    MOMY(:,:,:)    = 0.0_RP
    MOMZ(:,:,:)    = 0.0_RP
    RHOT(:,:,:)    = 0.0_RP
    QTRC(:,:,:,:)  = 0.0_RP
    PREC(:,:)      = 0.0_RP
    SWD (:,:)      = 0.0_RP
    LWD (:,:)      = 0.0_RP

    SFLX_MOMX (:,:) = 0.0_RP
    SFLX_MOMY (:,:) = 0.0_RP
    SFLX_MOMZ (:,:) = 0.0_RP
    SFLX_SWU  (:,:) = 0.0_RP
    SFLX_LWU  (:,:) = 0.0_RP
    SFLX_SH   (:,:) = 0.0_RP
    SFLX_LH   (:,:) = 0.0_RP
    SFLX_QVAtm(:,:) = 0.0_RP

    CNT_getDat2Atm = 0.0_RP

    return
  end subroutine CPL_AtmLnd_flushDat2Atm

  subroutine CPL_AtmLnd_flushDat2Lnd
    implicit none

    TG   (:,:) = 0.0_RP
    QvEfc(:,:) = 0.0_RP
    EMIT (:,:) = 0.0_RP
    ALB  (:,:) = 0.0_RP
    TCS  (:,:) = 0.0_RP
    DZg  (:,:) = 0.0_RP
    Z00  (:,:) = 0.0_RP
    Z0R  (:,:) = 0.0_RP
    Z0S  (:,:) = 0.0_RP
    Zt0  (:,:) = 0.0_RP
    ZtR  (:,:) = 0.0_RP
    ZtS  (:,:) = 0.0_RP
    Ze0  (:,:) = 0.0_RP
    ZeR  (:,:) = 0.0_RP
    ZeS  (:,:) = 0.0_RP

    SFLX_GH   (:,:) = 0.0_RP
    SFLX_PREC (:,:) = 0.0_RP
    SFLX_QVLnd(:,:) = 0.0_RP

    CNT_getDat2Lnd = 0.0_RP

    return
  end subroutine CPL_AtmLnd_flushDat2Lnd

! --- Private procedure

  subroutine ts_residual( &
      RES, DRES,                & ! (out)
      pMOMX, pMOMY, pMOMZ,      & ! (out)
      pSWU, pLWU, pSH, pLH, pGH ) ! (out)
    use mod_const, only: &
      GRAV   => CONST_GRAV,  &
      CPdry  => CONST_CPdry, &
      Rvap   => CONST_Rvap,  &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0
    use mod_grid, only: &
      CZ => GRID_CZ
    use mod_cpl_vars, only: &
      LST
    implicit none

    real(RP), intent(out) :: RES  (IA,JA)
    real(RP), intent(out) :: DRES (IA,JA)
    real(RP), intent(out) :: pMOMZ(IA,JA)
    real(RP), intent(out) :: pMOMX(IA,JA)
    real(RP), intent(out) :: pMOMY(IA,JA)
    real(RP), intent(out) :: pSWU (IA,JA)
    real(RP), intent(out) :: pLWU (IA,JA)
    real(RP), intent(out) :: pSH  (IA,JA)
    real(RP), intent(out) :: pLH  (IA,JA)
    real(RP), intent(out) :: pGH  (IA,JA)
    
    real(RP), parameter :: dTs   = 1.0E-10_RP ! delta skin temperature
    real(RP), parameter :: Cm0   = 1.0E-3_RP  ! bulk coef. for U*
    real(RP), parameter :: visck = 1.5E-5_RP  ! kinematic viscosity 

    ! work
    real(RP) :: pta(IA,JA)

    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Ustar ! friction velocity [m/s]

    real(RP) :: Z0, Zt, Ze ! roughness length (momentum,heat,tracer) [m]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: dCm, dCh, dCe

    real(RP) :: dLST
    real(RP) :: SatQvs, dSatQvs ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dpLWU, dpGH, dpSH1, dpSH2, dpLH1, dpLH2

    integer :: i, j, iw
    !---------------------------------------------------------------------------

    ! rho*theta -> potential temperature at cell centor
    do j = 1, JA
    do i = 1, IA
      pta(i,j) = RHOT(KS,i,j) / DENS(KS,i,j)
    end do
    end do

    do j = JS, JE
    do i = IS, IE
      ! at (u, y, layer)
      Uabs = sqrt( &
             ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i+1,j)                                     ) )**2 &
           + ( 2.0_RP *   MOMX(KS,i,j)                                                        )**2 &
           + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
           ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
      Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

      Z0 = max( Z00(i,j) + Z0R(i,j)/GRAV * Ustar*Ustar + Z0S(i,j)*visck / Ustar, Z0_min )
      Zt = max( Zt0(i,j) + ZtR(i,j)/GRAV * Ustar*Ustar + ZtS(i,j)*visck / Ustar, Zt_min )
      Ze = max( Ze0(i,j) + ZeR(i,j)/GRAV * Ustar*Ustar + ZeS(i,j)*visck / Ustar, Ze_min )

      call bulkcoef_uno( &
          Cm, Ch, Ce,                         & ! (out)
          ( pta(i,j) + pta(i+1,j) ) * 0.5_RP, & ! (in)
          ( LST(i,j) + LST(i+1,j) ) * 0.5_RP, & ! (in)
          Uabs, CZ(KS), Z0, Zt, Ze            ) ! (in)

      pMOMX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMX(KS,i,j)

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
      ! at (x, v, layer)
      Uabs = sqrt( &
             ( 0.5_RP * ( MOMZ(KS,i,j) + MOMZ(KS,i,j+1)                                     ) )**2 &
           + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
           + ( 2.0_RP *   MOMY(KS,i,j)                                                        )**2 &
           ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
      Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

      Z0 = max( Z00(i,j) + Z0R(i,j)/GRAV * Ustar*Ustar + Z0S(i,j)*visck / Ustar, Z0_min )
      Zt = max( Zt0(i,j) + ZtR(i,j)/GRAV * Ustar*Ustar + ZtS(i,j)*visck / Ustar, Zt_min )
      Ze = max( Ze0(i,j) + ZeR(i,j)/GRAV * Ustar*Ustar + ZeS(i,j)*visck / Ustar, Ze_min )

      call bulkcoef_uno( &
          Cm, Ch, Ce,                         & ! (out)
          ( pta(i,j) + pta(i,j+1) ) * 0.5_RP, & ! (in)
          ( LST(i,j) + LST(i,j+1) ) * 0.5_RP, & ! (in)
          Uabs, CZ(KS), Z0, Zt, Ze            ) ! (in)

      pMOMY(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMY(KS,i,j)

    enddo
    enddo

    do j = JS, JE+1
    do i = IS, IE+1
      ! at cell center
      Uabs = sqrt( &
             ( MOMZ(KS,i,j)                  )**2 &
           + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
           + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 &
           ) / DENS(KS,i,j) * 0.5_RP
      Ustar = max ( sqrt ( Cm0 ) * Uabs , Ustar_min )

      Z0 = max( Z00(i,j) + Z0R(i,j)/GRAV * Ustar*Ustar + Z0S(i,j)*visck / Ustar, Z0_min )
      Zt = max( Zt0(i,j) + ZtR(i,j)/GRAV * Ustar*Ustar + ZtS(i,j)*visck / Ustar, Zt_min )
      Ze = max( Ze0(i,j) + ZeR(i,j)/GRAV * Ustar*Ustar + ZeS(i,j)*visck / Ustar, Ze_min )

      call bulkcoef_uno( &
          Cm, Ch, Ce,              & ! (out)
          pta(i,j), LST(i,j),      & ! (in)
          Uabs, CZ(KS), Z0, Zt, Ze ) ! (in)

      ! saturation at surface
      SatQvs = satmixr( LST(i,j) )

      pMOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * MOMZ(KS,i,j) * 0.5_RP

      pSH(i,j) = CPdry * Ch * min(max(Uabs,U_minH),U_maxH) * ( LST(i,j)*DENS(KS,i,j) - RHOT(KS,i,j) )
      pLH(i,j) = LH0   * Ce * min(max(Uabs,U_minE),U_maxE) * DENS(KS,i,j) * ( SatQvs - QTRC(KS,i,j,I_QV) ) * QvEfc(i,j)
      pGH(i,j) = - 2.0_RP * TCS(i,j) * ( LST(i,j) - TG(i,j) ) / DZg(i,j)

      pSWU(i,j) = - ALB(i,j) * SWD(i,j)
      pLWU(i,j) =   EMIT(i,j) * STB * LST(i,j)**4

      ! calculation for residual
      RES(i,j) = - ( SWD(i,j) + pSWU(i,j) + LWD(i,j) + pLWU(i,j) ) &
                 - pSH(i,j) - pLH(i,j) + pGH(i,j)

      dLST  = LST(i,j) + dTs
      dSatQvs = dSatQvs * LH0 / ( Rvap * LST(i,j)**2 )

      call bulkcoef_uno( &
          dCm, dCh, dCe,           & ! (out)
          pta(i,j), dLST,          & ! (in)
          Uabs, CZ(KS), Z0, Zt, Ze ) ! (in)

      dpSH1 = CPdry * (dCh-Ch)/dTs * min(max(Uabs,U_minH),U_maxH) * ( LST(i,j)*DENS(KS,i,j) - RHOT(KS,i,j) )
      dpSH2 = CPdry * Ch * min(max(Uabs,U_minH),U_maxH) * DENS(KS,i,j)
      dpLH1 = LH0 * (dCe-Ce)/dTs * min(max(Uabs,U_minE),U_maxE) * DENS(KS,i,j) * ( SatQvs - QTRC(KS,i,j,I_QV) ) * QvEfc(i,j)
      dpLH2 = LH0 * Ce * min(max(Uabs,U_minE),U_maxE) * DENS(KS,i,j) * dSatQvs * QvEfc(i,j)

      dpLWU = 4.0_RP * EMIT(i,j) * STB * LST(i,j)**3
      dpGH  = - 2.0_RP * TCS(i,j) / DZg(i,j)

      ! calculation for d(residual)/dTs
      DRES(i,j) = - dpLWU + dpGH - dpSH1 - dpSH2 - dpLH1 - dpLH2

    enddo
    enddo

    return
  end subroutine ts_residual

  subroutine bulkcoef_uno( &
      Cm, Ch, Ce,                    & ! (out)
      pta, pts, za, uabs, z0, zt, ze ) ! (in)
    use mod_const, only : &
      GRAV   => CONST_GRAV,  &
      KARMAN => CONST_KARMAN
    implicit none

    ! argument
    real(RP), intent(out) :: Cm   ! momentum bulk coefficient [no unit]
    real(RP), intent(out) :: Ch   ! heat bulk coefficient [no unit]
    real(RP), intent(out) :: Ce   ! moisture bulk coefficient [no unit]
    
    real(RP), intent(in) :: pta  ! potential tempearature at 1st atm. layer [K]
    real(RP), intent(in) :: pts  ! skin potential temperature [K]
    real(RP), intent(in) :: za   ! height at 1st atm. layer [m]
    real(RP), intent(in) :: uabs ! wind speed at 1st atm. layer [m/s]
    real(RP), intent(in) :: z0   ! roughness length of momentum [m]
    real(RP), intent(in) :: zt   ! roughness length of heat [m]
    real(RP), intent(in) :: ze   ! roughness length of moisture [m]

    ! constant
    integer,  parameter :: nmax = 2 ! maximum iteration number (n=2: Uno et al. 1995)

    ! parameters of bulk transfer coefficient
    real(RP), parameter :: dRiB = 1.0E-10_RP ! delta RiB [no unit]
    real(RP), parameter :: tPrn = 0.74E0_RP  ! turbulent Prandtl number (Businger et al. 1971)
    real(RP), parameter :: LFb  = 9.4E0_RP   ! Louis factor b (Louis 1979)
    real(RP), parameter :: LFbp = 4.7E0_RP   ! Louis factor b' (Louis 1979)
    real(RP), parameter :: LFdm = 7.4E0_RP   ! Louis factor d for momemtum (Louis 1979)
    real(RP), parameter :: LFdh = 5.3E0_RP   ! Louis factor d for heat (Louis 1979)

    ! variables
    integer :: n
    real(RP) :: RiBT, RiB0, dRiB0 ! bulk Richardson number [no unit]
    real(RP) :: C0 ! initial drag coefficient [no unit]
    real(RP) :: psi, dpsi
    real(RP) :: fm, fh
    real(RP) :: dfm, dfh
    real(RP) :: LFcm, LFch
    real(RP) :: res, dres

    RiBT = GRAV * za * ( pta - pts ) / ( ( pts + pta ) * 0.5_RP * uabs**2 )
    C0   = ( KARMAN / log( za/z0 ) )**2
    LFcm = LFb * LFdm * C0 * sqrt( za/z0 )
    LFch = LFb * LFdh * C0 * sqrt( za/z0 )

    RiB0 = RiBT
    do n = 1, nmax
      ! calculate psi
      if( RiB0 < 0.0_RP ) then
        ! unstable condition
        fm = 1.0_RP - ( LFb * RiB0 ) / ( 1.0_RP + LFcm * sqrt( abs( RiB0 ) ) )
        fh = 1.0_RP - ( LFb * RiB0 ) / ( 1.0_RP + LFch * sqrt( abs( RiB0 ) ) )
      else
        ! stable condition
        fm = 1.0_RP / ( 1.0_RP + LFbp * RiB0 ) ** 2
        fh = 1.0_RP / ( 1.0_RP + LFbp * RiB0 ) ** 2
      end if
      psi  = tPrn * log( za/z0 ) * sqrt( fm ) / fh

      ! calculate dpsi
      dRiB0 = RiB0 + dRiB
      if( dRiB0 < 0.0_RP ) then
        ! unstable condition
        dfm = 1.0_RP - ( LFb * dRiB0 ) / ( 1.0_RP + LFcm * sqrt( abs( dRiB0 ) ) )
        dfh = 1.0_RP - ( LFb * dRiB0 ) / ( 1.0_RP + LFch * sqrt( abs( dRiB0 ) ) )
      else
        ! stable condition
        dfm = 1.0_RP / ( 1.0_RP + LFbp * dRiB0 ) ** 2
        dfh = 1.0_RP / ( 1.0_RP + LFbp * dRiB0 ) ** 2
      end if
      dpsi = ( tPrn * log( za/z0 ) * sqrt( dfm ) / dfh - psi ) / dRiB

      res  = RiB0 * tPrn * log( z0/zt ) + ( RiB0 - RiBT ) * psi
      dres = tPrn * log( z0/zt ) + psi + ( RiB0 - RiBT ) * dpsi

      ! update the bulk Richardson number
      RiB0 = RiB0 - res / dres

      if( abs( res ) < res_min ) then
        ! finish iteration
        exit
      end if

    end do

    if( n > 3 .and. n > nmax ) then
      if( IO_L ) write(IO_FID_LOG,*) 'Error: reach maximum iteration.'
    end if

    Cm = C0 * fm
    Ch = C0 * fh / tPrn / ( tPrn * log( z0/zt ) / psi + 1.0_RP )
    Ce = C0 * fh / tPrn / ( tPrn * log( z0/ze ) / psi + 1.0_RP )

    Cm = min( max( Cm, Cm_min ), Cm_max )
    Ch = min( max( Ch, Ch_min ), Ch_max )
    Ce = min( max( Ce, Ce_min ), Ce_max )

    return

  end subroutine bulkcoef_uno

  function satmixr( temp )
    use mod_const, only : &
      P00    => CONST_PRE00, &
      EPSvap => CONST_EPSvap
    implicit none

    ! argument
    real(8), intent(in) :: temp ! temperature [K]
!    real(8), intent(in) :: pres ! pressure [Pa]
    ! function
    real(8) :: satmixr ! saturated mixing ratio [kg/kg]

!    satmixr = EPSvap * satvapor( temp ) / ( pres - esat )
    satmixr = EPSvap * satvapor( temp ) / P00

  end function satmixr

  function satvapor( temp )
    use mod_const, only: &
      T00   => CONST_TEM00, &
      LH0   => CONST_LH0,   &
      Rvap  => CONST_Rvap,  &
      PSAT0 => CONST_PSAT0
    implicit none

    ! argument
    real(8), intent(in) :: temp ! temperature [K]
    ! function
    real(8) :: satvapor ! saturated water vapor pressure [Pa]

    ! Tetens (1930)
    satvapor = PSAT0 * exp( LH0/Rvap * ( 1.0d0/T00 - 1.0d0/temp ) )

  end function satvapor

end module mod_cpl_atmos_land
