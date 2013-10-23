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
  ! limiter
  real(RP), private, parameter :: res_min   =  1.0E-10_RP ! minimum number of residual in the Newton-Raphson Method
  real(RP), private, parameter :: Ustar_min =  1.0E-3_RP ! minimum limit of U*

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

  ! work
  real(RP), private, save :: pDENS(KA,IA,JA)
  real(RP), private, save :: pMOMX(KA,IA,JA)
  real(RP), private, save :: pMOMY(KA,IA,JA)
  real(RP), private, save :: pMOMZ(KA,IA,JA)
  real(RP), private, save :: pRHOT(KA,IA,JA)
  real(RP), private, save :: pQTRC(KA,IA,JA,QA)
  real(RP), private, save :: pPREC(IA,JA)
  real(RP), private, save :: pSWD (IA,JA)
  real(RP), private, save :: pLWD (IA,JA)

  real(RP), private, save :: pTG   (IA,JA)
  real(RP), private, save :: pQvEfc(IA,JA)
  real(RP), private, save :: pEMIT (IA,JA)
  real(RP), private, save :: pALB  (IA,JA)
  real(RP), private, save :: pTCS  (IA,JA)
  real(RP), private, save :: pDZg  (IA,JA)
  real(RP), private, save :: pZ0M  (IA,JA)
  real(RP), private, save :: pZ0H  (IA,JA)
  real(RP), private, save :: pZ0E  (IA,JA)

  real(RP), private, save :: pSFLX_MOMX(IA,JA)
  real(RP), private, save :: pSFLX_MOMY(IA,JA)
  real(RP), private, save :: pSFLX_MOMZ(IA,JA)
  real(RP), private, save :: pSFLX_SWU (IA,JA)
  real(RP), private, save :: pSFLX_LWU (IA,JA)
  real(RP), private, save :: pSFLX_SH  (IA,JA)
  real(RP), private, save :: pSFLX_LH  (IA,JA)
  real(RP), private, save :: pSFLX_GH  (IA,JA)

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
       CPL_flushAtm,            &
       CPL_flushLnd,            &
       CPL_AtmLnd_flushCPL,     &
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

    call CPL_flushAtm
    call CPL_flushLnd
    call CPL_AtmLnd_flushCPL

    return
  end subroutine CPL_AtmLnd_setup

  subroutine CPL_AtmLnd_solve
    use mod_const, only: &
       LH0 => CONST_LH0
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_vars, only: &
       CPL_AtmLnd_putCPL,     &
       CPL_AtmLnd_getAtm2CPL, &
       CPL_AtmLnd_getLnd2CPL, &
       CPL_AtmLnd_flushCPL,   &
       LST
    implicit none

    ! parameters
    integer,  parameter :: nmax     = 100       ! maximum iteration number
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)

    ! works
    integer :: i, j, n

    real(RP) :: RES  (IA,JA)
    real(RP) :: DRES (IA,JA)
    real(RP) :: oldRES(IA,JA) ! RES in previous step
    real(RP) :: redf  (IA,JA) ! reduced factor

    if( IO_L ) write(IO_FID_LOG,*) '*** CPL solve: Atmos-Land'

    call CPL_AtmLnd_getAtm2CPL( &
      pDENS, pMOMX, pMOMY, pMOMZ, pRHOT, & !(out)
      pQTRC, pPREC, pSWD, pLWD           ) ! (out)

    call CPL_AtmLnd_getLnd2CPL( &
      pTG, pQvEfc, pEMIT, & ! (out)
      pALB, pTCS, pDZg,   & ! (out)
      pZ0M, pZ0H, pZ0E    ) ! (out)

    redf  (:,:) = 1.0_RP
    oldRES(:,:) = 1.0E+5_RP

    do n = 1, nmax

      ! calc. surface fluxes
      call ts_residual( &
        RES, DRES ) ! (out)
        
      do j = JS-1, JE+1
      do i = IS-1, IE+1

        if( redf(i,j) < 0.0_RP ) then
          redf(i,j) = 1.0_RP
        end if

        if( abs(RES(i,j)) > abs(oldRES(i,j)) ) then
          redf(i,j) = max( TFa*redf(i,j), redf_min )
        else
          redf(i,j) = min( TFb*redf(i,j), redf_max )
        end if

        if( DRES(i,j) > 0.0_RP ) then
          redf(i,j) = -1.0_RP
        end if

        ! update surface temperature
        LST(i,j) = LST(i,j) - redf(i,j) * RES(i,j)/DRES(i,j)

        ! save residual in this step
        oldRES(i,j) = RES(i,j)

      end do
      end do

      if( maxval(abs(RES(IS-1:IE+1,JS-1:JE+1))) < res_min ) then
        ! iteration converged
        exit
      end if

    end do

    if( n > nmax ) then
      ! not converged and stop program
      if( IO_L ) write(IO_FID_LOG,*) 'Error: surface tempearture is not converged.'
      call PRC_MPIstop
    end if

    ! put residual in ground heat flux
    pSFLX_GH(:,:) = pSFLX_GH(:,:) - RES(:,:)

    call CPL_AtmLnd_putCPL( &
      pSFLX_MOMX, &
      pSFLX_MOMY, &
      pSFLX_MOMZ, &
      pSFLX_SWU,  &
      pSFLX_LWU,  &
      pSFLX_SH,   &
      pSFLX_LH,   &
      pSFLX_GH,   &
      pPREC       )

    call CPL_AtmLnd_flushCPL

    return
  end subroutine CPL_AtmLnd_solve

  subroutine CPL_AtmLnd_unsolve
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_vars, only: &
       CPL_AtmLnd_putCPL,     &
       CPL_AtmLnd_getAtm2CPL, &
       CPL_AtmLnd_getLnd2CPL, &
       CPL_AtmLnd_flushCPL
    implicit none

    if( IO_L ) write(IO_FID_LOG,*) '*** CPL unsolve: Atmos-Land'

    call CPL_AtmLnd_getAtm2CPL( &
      pDENS, pMOMX, pMOMY, pMOMZ, pRHOT, & !(out)
      pQTRC, pPREC, pSWD, pLWD           ) ! (out)

    call CPL_AtmLnd_getLnd2CPL( &
      pTG, pQvEfc, pEMIT, & ! (out)
      pALB, pTCS, pDZg,   & ! (out)
      pZ0M, pZ0H, pZ0E    ) ! (out)

    ! calc. surface fluxes
    call ts_known

    call CPL_AtmLnd_putCPL( &
      pSFLX_MOMX, &
      pSFLX_MOMY, &
      pSFLX_MOMZ, &
      pSFLX_SWU,  &
      pSFLX_LWU,  &
      pSFLX_SH,   &
      pSFLX_LH,   &
      pSFLX_GH,   &
      pPREC       )

    call CPL_AtmLnd_flushCPL

    return
  end subroutine CPL_AtmLnd_unsolve

! --- Private procedure

  subroutine ts_residual( &
      RES, DRES ) ! (out)
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
    
    real(RP), parameter :: dTs   = 1.0E-10_RP ! delta skin temperature

    ! work
    real(RP) :: pta(IA,JA)
    real(RP) :: Uabs ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: dCm, dCh, dCe

    real(RP) :: dLST
    real(RP) :: SatQvs, dSatQvs ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dpSFLX_LWU, dpSFLX_GH
    real(RP) :: dpSFLX_SH1, dpSFLX_SH2
    real(RP) :: dpSFLX_LH1, dpSFLX_LH2

    integer :: i, j
    !---------------------------------------------------------------------------

    ! rho*theta -> potential temperature at cell centor
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      pta(i,j) = pRHOT(KS,i,j) / pDENS(KS,i,j)
    end do
    end do

    ! at (u, y, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( pMOMZ(KS,i,j) + pMOMZ(KS,i+1,j)                                       ) )**2 &
           + ( 2.0_RP *   pMOMX(KS,i,j)                                                           )**2 &
           + ( 0.5_RP * ( pMOMY(KS,i,j-1) + pMOMY(KS,i,j) + pMOMY(KS,i+1,j-1) + pMOMY(KS,i+1,j) ) )**2 &
           ) / ( pDENS(KS,i,j) + pDENS(KS,i+1,j) )

      call bulkcoef_uno( &
          Cm, Ch, Ce,                         & ! (out)
          ( pta(i,j) + pta(i+1,j) ) * 0.5_RP, & ! (in)
          ( LST(i,j) + LST(i+1,j) ) * 0.5_RP, & ! (in)
          Uabs, CZ(KS),                       & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j)     ) ! (in)

      pSFLX_MOMX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMX(KS,i,j)
    enddo
    enddo

    ! at (x, v, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( pMOMZ(KS,i,j) + pMOMZ(KS,i,j+1)                                       ) )**2 &
           + ( 0.5_RP * ( pMOMX(KS,i-1,j) + pMOMX(KS,i,j) + pMOMX(KS,i-1,j+1) + pMOMX(KS,i,j+1) ) )**2 &
           + ( 2.0_RP *   pMOMY(KS,i,j)                                                           )**2 &
           ) / ( pDENS(KS,i,j) + pDENS(KS,i,j+1) )

      call bulkcoef_uno( &
          Cm, Ch, Ce,                         & ! (out)
          ( pta(i,j) + pta(i,j+1) ) * 0.5_RP, & ! (in)
          ( LST(i,j) + LST(i,j+1) ) * 0.5_RP, & ! (in)
          Uabs, CZ(KS),                       & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j)     ) ! (in)

      pSFLX_MOMY(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMY(KS,i,j)
    enddo
    enddo

    ! at cell center
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      Uabs = sqrt( &
             ( pMOMZ(KS,i,j)                   )**2 &
           + ( pMOMX(KS,i-1,j) + pMOMX(KS,i,j) )**2 &
           + ( pMOMY(KS,i,j-1) + pMOMY(KS,i,j) )**2 &
           ) / pDENS(KS,i,j) * 0.5_RP

      call bulkcoef_uno( &
          Cm, Ch, Ce,                     & ! (out)
          pta(i,j), LST(i,j),             & ! (in)
          Uabs, CZ(KS),                   & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j) ) ! (in)

      pSFLX_MOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMZ(KS,i,j) * 0.5_RP

      ! saturation at surface
      SatQvs = satmixr( LST(i,j) )

      pSFLX_SH(i,j) = CPdry * Ch * min(max(Uabs,U_minH),U_maxH) * ( LST(i,j)*pDENS(KS,i,j) - pRHOT(KS,i,j) )
      pSFLX_LH(i,j) = LH0   * Ce * min(max(Uabs,U_minE),U_maxE) * pDENS(KS,i,j) * ( SatQvs - pQTRC(KS,i,j,I_QV) ) * pQvEfc(i,j)
      pSFLX_GH(i,j) = 2.0_RP * pTCS(i,j) * ( pTG(i,j) - LST(i,j) ) / pDZg(i,j)

      pSFLX_SWU(i,j) = pALB(i,j) * pSWD(i,j)
      pSFLX_LWU(i,j) = pEMIT(i,j) * STB * LST(i,j)**4

      ! calculation for residual
      RES(i,j) = pSWD(i,j) - pSFLX_SWU(i,j) + pLWD(i,j) - pSFLX_LWU(i,j) - pSFLX_SH(i,j) - pSFLX_LH(i,j) + pSFLX_GH(i,j)

      dLST  = LST(i,j) + dTs
      dSatQvs = SatQvs * LH0 / ( Rvap * LST(i,j)**2 )

      call bulkcoef_uno( &
          dCm, dCh, dCe,                  & ! (out)
          pta(i,j), dLST,                 & ! (in)
          Uabs, CZ(KS),                   & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j) ) ! (in)

      dpSFLX_SH1 = CPdry * (dCh-Ch)/dTs * min(max(Uabs,U_minH),U_maxH) * ( LST(i,j)*pDENS(KS,i,j) - pRHOT(KS,i,j) )
      dpSFLX_SH2 = CPdry * Ch * min(max(Uabs,U_minH),U_maxH) * pDENS(KS,i,j)
      dpSFLX_LH1 = LH0 * (dCe-Ce)/dTs * min(max(Uabs,U_minE),U_maxE) * pDENS(KS,i,j) * ( SatQvs - pQTRC(KS,i,j,I_QV) ) * pQvEfc(i,j)
      dpSFLX_LH2 = LH0 * Ce * min(max(Uabs,U_minE),U_maxE) * pDENS(KS,i,j) * dSatQvs * pQvEfc(i,j)

      dpSFLX_LWU = 4.0_RP * pEMIT(i,j) * STB * LST(i,j)**3
      dpSFLX_GH  = - 2.0_RP * pTCS(i,j) / pDZg(i,j)

      ! calculation for d(residual)/dTs
      DRES(i,j) = - dpSFLX_LWU - dpSFLX_SH1 - dpSFLX_SH2 - dpSFLX_LH1 - dpSFLX_LH2 + dpSFLX_GH
    enddo
    enddo

    return
  end subroutine ts_residual

  subroutine ts_known
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

    ! work
    real(RP) :: pta(IA,JA)
    real(RP) :: Uabs ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: SatQvs ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    ! rho*theta -> potential temperature at cell centor
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      pta(i,j) = pRHOT(KS,i,j) / pDENS(KS,i,j)
    end do
    end do

    ! at (u, y, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( pMOMZ(KS,i,j) + pMOMZ(KS,i+1,j)                                       ) )**2 &
           + ( 2.0_RP *   pMOMX(KS,i,j)                                                           )**2 &
           + ( 0.5_RP * ( pMOMY(KS,i,j-1) + pMOMY(KS,i,j) + pMOMY(KS,i+1,j-1) + pMOMY(KS,i+1,j) ) )**2 &
           ) / ( pDENS(KS,i,j) + pDENS(KS,i+1,j) )

      call bulkcoef_uno( &
          Cm, Ch, Ce,                         & ! (out)
          ( pta(i,j) + pta(i+1,j) ) * 0.5_RP, & ! (in)
          ( LST(i,j) + LST(i+1,j) ) * 0.5_RP, & ! (in)
          Uabs, CZ(KS),                       & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j)     ) ! (in)

      pSFLX_MOMX(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMX(KS,i,j)
    enddo
    enddo

    ! at (x, v, layer)
    do j = JS, JE
    do i = IS, IE
      Uabs = sqrt( &
             ( 0.5_RP * ( pMOMZ(KS,i,j) + pMOMZ(KS,i,j+1)                                       ) )**2 &
           + ( 0.5_RP * ( pMOMX(KS,i-1,j) + pMOMX(KS,i,j) + pMOMX(KS,i-1,j+1) + pMOMX(KS,i,j+1) ) )**2 &
           + ( 2.0_RP *   pMOMY(KS,i,j)                                                           )**2 &
           ) / ( pDENS(KS,i,j) + pDENS(KS,i,j+1) )

      call bulkcoef_uno( &
          Cm, Ch, Ce,                         & ! (out)
          ( pta(i,j) + pta(i,j+1) ) * 0.5_RP, & ! (in)
          ( LST(i,j) + LST(i,j+1) ) * 0.5_RP, & ! (in)
          Uabs, CZ(KS),                       & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j)     ) ! (in)

      pSFLX_MOMY(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMY(KS,i,j)
    enddo
    enddo

    ! at cell center
    do j = JS-1, JE+1
    do i = IS-1, IE+1
      Uabs = sqrt( &
             ( pMOMZ(KS,i,j)                   )**2 &
           + ( pMOMX(KS,i-1,j) + pMOMX(KS,i,j) )**2 &
           + ( pMOMY(KS,i,j-1) + pMOMY(KS,i,j) )**2 &
           ) / pDENS(KS,i,j) * 0.5_RP

      call bulkcoef_uno( &
          Cm, Ch, Ce,                     & ! (out)
          pta(i,j), LST(i,j),             & ! (in)
          Uabs, CZ(KS),                   & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j) ) ! (in)

      pSFLX_MOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMZ(KS,i,j) * 0.5_RP

      ! saturation at surface
      SatQvs = satmixr( LST(i,j) )

      pSFLX_SH(i,j) = CPdry * Ch * min(max(Uabs,U_minH),U_maxH) * ( LST(i,j)*pDENS(KS,i,j) - pRHOT(KS,i,j) )
      pSFLX_LH(i,j) = LH0   * Ce * min(max(Uabs,U_minE),U_maxE) * pDENS(KS,i,j) * ( SatQvs - pQTRC(KS,i,j,I_QV) ) * pQvEfc(i,j)
      pSFLX_GH(i,j) = 2.0_RP * pTCS(i,j) * ( pTG(i,j) - LST(i,j) ) / pDZg(i,j)

      pSFLX_SWU(i,j) = pALB(i,j) * pSWD(i,j)
      pSFLX_LWU(i,j) = pEMIT(i,j) * STB * LST(i,j)**4
    enddo
    enddo

    return
  end subroutine ts_known

  subroutine bulkcoef_uno( &
      Cm, Ch, Ce,                    & ! (out)
      pta, pts, za, uabs, z0, zt, ze ) ! (in)
    use mod_const, only : &
      GRAV   => CONST_GRAV,  &
      KARMAN => CONST_KARMAN
    use mod_process, only: &
       PRC_MPIstop
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
    real(RP), parameter :: tPrn = 0.74_RP    ! turbulent Prandtl number (Businger et al. 1971)
    real(RP), parameter :: LFb  = 9.4_RP     ! Louis factor b (Louis 1979)
    real(RP), parameter :: LFbp = 4.7_RP     ! Louis factor b' (Louis 1979)
    real(RP), parameter :: LFdm = 7.4_RP     ! Louis factor d for momemtum (Louis 1979)
    real(RP), parameter :: LFdh = 5.3_RP     ! Louis factor d for heat (Louis 1979)

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
        fm = 1.0_RP / ( 1.0_RP + LFbp * RiB0 )**2
        fh = 1.0_RP / ( 1.0_RP + LFbp * RiB0 )**2
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
        dfm = 1.0_RP / ( 1.0_RP + LFbp * dRiB0 )**2
        dfh = 1.0_RP / ( 1.0_RP + LFbp * dRiB0 )**2
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
      if( IO_L ) write(IO_FID_LOG,*) 'Error: reach maximum iteration in the function of bulkcoef_uno.'
      call PRC_MPIstop
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
    satvapor = PSAT0 * exp( LH0/Rvap * ( 1.0_RP/T00 - 1.0_RP/temp ) )

  end function satvapor

end module mod_cpl_atmos_land
