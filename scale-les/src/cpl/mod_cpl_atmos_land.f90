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
  use mod_precision
  use mod_index
  use mod_stdio
  use mod_prof
  use mod_tracer
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
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: bulkcoef_uno
  private :: sfcval_estimate
  private :: satmixr
  private :: satvapor

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! limiter
  real(RP), private, parameter :: Cm_min =   1.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_min =   1.0E-5_RP !                       T
  real(RP), private, parameter :: Ce_min =   1.0E-5_RP !                       q
  real(RP), private, parameter :: Cm_max =   1.0_RP    ! maximum bulk coef. of u,v,w
  real(RP), private, parameter :: Ch_max =   1.0_RP    !                       T
  real(RP), private, parameter :: Ce_max =   1.0_RP    !                       q

  real(RP), private, parameter :: U_minM =   0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_minH =   0.0_RP  !                   T
  real(RP), private, parameter :: U_minE =   0.0_RP  !                   q
  real(RP), private, parameter :: U_maxM = 100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH = 100.0_RP  !                   T
  real(RP), private, parameter :: U_maxE = 100.0_RP  !                   q

  ! work
  real(RP), private, allocatable :: pDENS(:,:,:)
  real(RP), private, allocatable :: pMOMX(:,:,:)
  real(RP), private, allocatable :: pMOMY(:,:,:)
  real(RP), private, allocatable :: pMOMZ(:,:,:)
  real(RP), private, allocatable :: pRHOT(:,:,:)
  real(RP), private, allocatable :: pQTRC(:,:,:,:)
  real(RP), private, allocatable :: pPREC(:,:)
  real(RP), private, allocatable :: pSWD (:,:)
  real(RP), private, allocatable :: pLWD (:,:)

  real(RP), private, allocatable :: pTG   (:,:)
  real(RP), private, allocatable :: pQvEfc(:,:)
  real(RP), private, allocatable :: pEMIT (:,:)
  real(RP), private, allocatable :: pALB  (:,:)
  real(RP), private, allocatable :: pTCS  (:,:)
  real(RP), private, allocatable :: pDZg  (:,:)
  real(RP), private, allocatable :: pZ0M  (:,:)
  real(RP), private, allocatable :: pZ0H  (:,:)
  real(RP), private, allocatable :: pZ0E  (:,:)

  real(RP), private, allocatable :: pSFLX_MOMX(:,:)
  real(RP), private, allocatable :: pSFLX_MOMY(:,:)
  real(RP), private, allocatable :: pSFLX_MOMZ(:,:)
  real(RP), private, allocatable :: pSFLX_SWU (:,:)
  real(RP), private, allocatable :: pSFLX_LWU (:,:)
  real(RP), private, allocatable :: pSFLX_SH  (:,:)
  real(RP), private, allocatable :: pSFLX_LH  (:,:)
  real(RP), private, allocatable :: pSFLX_GH  (:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  subroutine CPL_AtmLnd_setup
    use mod_cpl_vars, only: &
       CPL_flushAtm,            &
       CPL_flushLnd,            &
       CPL_AtmLnd_flushCPL
    implicit none
    !---------------------------------------------------------------------------

    allocate( pDENS(KA,IA,JA) )
    allocate( pMOMX(KA,IA,JA) )
    allocate( pMOMY(KA,IA,JA) )
    allocate( pMOMZ(KA,IA,JA) )
    allocate( pRHOT(KA,IA,JA) )
    allocate( pQTRC(KA,IA,JA,QA) )
    allocate( pPREC(IA,JA) )
    allocate( pSWD (IA,JA) )
    allocate( pLWD (IA,JA) )

    allocate( pTG   (IA,JA) )
    allocate( pQvEfc(IA,JA) )
    allocate( pEMIT (IA,JA) )
    allocate( pALB  (IA,JA) )
    allocate( pTCS  (IA,JA) )
    allocate( pDZg  (IA,JA) )
    allocate( pZ0M  (IA,JA) )
    allocate( pZ0H  (IA,JA) )
    allocate( pZ0E  (IA,JA) )

    allocate( pSFLX_MOMX(IA,JA) )
    allocate( pSFLX_MOMY(IA,JA) )
    allocate( pSFLX_MOMZ(IA,JA) )
    allocate( pSFLX_SWU (IA,JA) )
    allocate( pSFLX_LWU (IA,JA) )
    allocate( pSFLX_SH  (IA,JA) )
    allocate( pSFLX_LH  (IA,JA) )
    allocate( pSFLX_GH  (IA,JA) )

    call CPL_flushAtm
    call CPL_flushLnd
    call CPL_AtmLnd_flushCPL

    return
  end subroutine CPL_AtmLnd_setup

  subroutine CPL_AtmLnd_solve
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_vars, only: &
       CPL_AtmLnd_putCPL,     &
       CPL_AtmLnd_getAtm2CPL, &
       CPL_AtmLnd_getLnd2CPL, &
       LST
    implicit none

    ! parameters
    integer,  parameter :: nmax     = 100       ! maximum iteration number
    real(RP), parameter :: redf_min = 1.0E-2_RP ! minimum reduced factor
    real(RP), parameter :: redf_max = 1.0_RP    ! maximum reduced factor
    real(RP), parameter :: TFa      = 0.5_RP    ! factor a in Tomita (2009)
    real(RP), parameter :: TFb      = 1.1_RP    ! factor b in Tomita (2009)
    real(RP), parameter :: res_min  = 1.0_RP    ! minimum number of residual

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
      call ts_residual( RES, DRES ) ! (out)

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
        LST(i,j)  = LST(i,j) - redf(i,j) * RES(i,j)/DRES(i,j)

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

    return
  end subroutine CPL_AtmLnd_solve

  subroutine CPL_AtmLnd_unsolve
    use mod_process, only: &
       PRC_MPIstop
    use mod_cpl_vars, only: &
       CPL_AtmLnd_putCPL,     &
       CPL_AtmLnd_getAtm2CPL, &
       CPL_AtmLnd_getLnd2CPL
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

    return
  end subroutine CPL_AtmLnd_unsolve

! --- Private procedure

  subroutine ts_residual( RES, DRES ) ! (out)
    use mod_const, only: &
      GRAV   => CONST_GRAV,  &
      CPdry  => CONST_CPdry, &
      RovCP  => CONST_RovCP, &
      Rvap   => CONST_Rvap,  &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0,   &
      P00    => CONST_PRE00
    use mod_grid, only: &
      CZ => GRID_CZ
    use mod_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_cpl_vars, only: &
      LST
    implicit none

    real(RP), intent(out) :: RES  (IA,JA)
    real(RP), intent(out) :: DRES (IA,JA)

    real(RP), parameter :: dTs   = 1.0E-8_RP ! delta skin temperature

    ! work
    real(RP) :: tem(KA,IA,JA) ! temperature [K]
    real(RP) :: pre(KA,IA,JA) ! pressure [Pa]

    real(RP) :: rhos(IA,JA) ! surface density [kg/m3]
    real(RP) :: pres(IA,JA) ! surface pressure [Pa]
    real(RP) :: zs  (IA,JA) ! surface height [m]
    real(RP) :: pta (IA,JA) ! surface air temperature [K]

    real(RP) :: Uabs ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: dCm, dCh, dCe ! bulk transfer coeff. [no unit]

    real(RP) :: SatQvs, dSatQvs ! saturation water vapor mixing ratio at surface [kg/kg]
    real(RP) :: dpSFLX_LWU, dpSFLX_GH
    real(RP) :: dpSFLX_SH,  dpSFLX_LH

    integer :: i, j
    !---------------------------------------------------------------------------

    call THERMODYN_temp_pres( tem  (:,:,:),  & ! [OUT]
                              pre  (:,:,:),  & ! [OUT]
                              pDENS(:,:,:),  & ! [IN]
                              pRHOT(:,:,:),  & ! [IN]
                              pQTRC(:,:,:,:) ) ! [IN]

    ! surface height
    zs(:,:) = 0.0_RP

    call sfcval_estimate( rhos (:,:),   & ! [out]
                          pres (:,:),   & ! [out]
                          pDENS(:,:,:), & ! [in]
                          pre  (:,:,:), & ! [in]
                          zs   (:,:)    ) ! [in]

    do j = 1, JA
    do i = 1, IA
      pta(i,j) = tem(KS,i,j) * ( pres(i,j) / pre(KS,i,j) )**RovCP
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
          CZ(KS), Uabs,                       & ! (in)
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
          CZ(KS), Uabs,                       & ! (in)
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
          CZ(KS), Uabs,                   & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j) ) ! (in)

      pSFLX_MOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMZ(KS,i,j) * 0.5_RP

      ! saturation at surface
      SatQvs = satmixr( LST(i,j), pres(i,j) )

      pSFLX_SH(i,j)  = CPdry * min(max(Uabs,U_minH),U_maxH) * rhos(i,j) * Ch * ( LST(i,j) - pta(i,j) )
      pSFLX_LH(i,j)  = LH0   * min(max(Uabs,U_minE),U_maxE) * rhos(i,j) * pQvEfc(i,j) * Ce * ( SatQvs - pQTRC(KS,i,j,I_QV) )
      pSFLX_GH(i,j)  = -2.0_RP * pTCS(i,j) * ( LST(i,j) - pTG(i,j)  ) / pDZg(i,j)
      pSFLX_SWU(i,j) = pALB(i,j) * pSWD(i,j)
      pSFLX_LWU(i,j) = pEMIT(i,j) * STB * LST(i,j)**4

      ! calculation for residual
      RES(i,j) = pSWD(i,j) - pSFLX_SWU(i,j) + pLWD(i,j) - pSFLX_LWU(i,j) - pSFLX_SH(i,j) - pSFLX_LH(i,j) + pSFLX_GH(i,j)

      dSatQvs = SatQvs * LH0 / ( Rvap * LST(i,j)**2 )

      call bulkcoef_uno( &
          dCm, dCh, dCe,                  & ! (out)
          pta(i,j), LST(i,j)+dTs,         & ! (in)
          CZ(KS), Uabs,                   & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j) ) ! (in)

      dpSFLX_SH  = CPdry * min(max(Uabs,U_minH),U_maxH) * rhos(i,j) &
                 * ( (dCh-Ch)/dTs * ( LST(i,j) - pta(i,j) ) + Ch )
      dpSFLX_LH  = LH0   * min(max(Uabs,U_minE),U_maxE) * rhos(i,j) * pQvEfc(i,j) &
                 * ( (dCe-Ce)/dTs * ( SatQvs - pQTRC(KS,i,j,I_QV) ) + Ce * dSatQvs )
      dpSFLX_GH  = -2.0_RP * pTCS(i,j) / pDZg(i,j)
      dpSFLX_LWU = 4.0_RP * pEMIT(i,j) * STB * LST(i,j)**3

      ! calculation for d(residual)/dTs
      DRES(i,j) = - dpSFLX_LWU - dpSFLX_SH - dpSFLX_LH + dpSFLX_GH
    enddo
    enddo

    return
  end subroutine ts_residual

  subroutine ts_known
    use mod_const, only: &
      GRAV   => CONST_GRAV,  &
      CPdry  => CONST_CPdry, &
      RovCP  => CONST_RovCP, &
      Rvap   => CONST_Rvap,  &
      STB    => CONST_STB,   &
      LH0    => CONST_LH0,   &
      P00    => CONST_PRE00
    use mod_grid, only: &
      CZ => GRID_CZ
    use mod_atmos_thermodyn, only: &
       THERMODYN_temp_pres => ATMOS_THERMODYN_temp_pres
    use mod_cpl_vars, only: &
      LST
    implicit none

    ! work
    real(RP) :: tem(KA,IA,JA) ! temperature [K]
    real(RP) :: pre(KA,IA,JA) ! pressure [Pa]

    real(RP) :: rhos(IA,JA) ! surface density [kg/m3]
    real(RP) :: pres(IA,JA) ! surface pressure [Pa]
    real(RP) :: zs  (IA,JA) ! surface height [m]
    real(RP) :: pta (IA,JA) ! surface air tempearature [K]

    real(RP) :: Uabs ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce ! bulk transfer coeff. [no unit]
    real(RP) :: SatQvs ! saturation water vapor mixing ratio at surface [kg/kg]

    integer :: i, j
    !---------------------------------------------------------------------------

    call THERMODYN_temp_pres( tem  (:,:,:),  & ! [OUT]
                              pre  (:,:,:),  & ! [OUT]
                              pDENS(:,:,:),  & ! [IN]
                              pRHOT(:,:,:),  & ! [IN]
                              pQTRC(:,:,:,:) ) ! [IN]

    ! surface height
    zs(:,:) = 0.0_RP

    call sfcval_estimate( rhos (:,:),   & ! [out]
                          pres (:,:),   & ! [out]
                          pDENS(:,:,:), & ! [in]
                          pre  (:,:,:), & ! [in]
                          zs   (:,:)    ) ! [in]

    do j = 1, JA
    do i = 1, IA
      pta(i,j) = tem(KS,i,j) * ( pres(i,j) / pre(KS,i,j) )**RovCP
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
          CZ(KS), Uabs,                       & ! (in)
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
          CZ(KS), Uabs,                       & ! (in)
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
          CZ(KS), Uabs,                   & ! (in)
          pZ0M(i,j), pZ0H(i,j), pZ0E(i,j) ) ! (in)

      pSFLX_MOMZ(i,j) = - Cm * min(max(Uabs,U_minM),U_maxM) * pMOMZ(KS,i,j) * 0.5_RP

      ! saturation at surface
      SatQvs = satmixr( LST(i,j), pres(i,j) )

      pSFLX_SH(i,j)  = CPdry * min(max(Uabs,U_minH),U_maxH) * rhos(i,j) * Ch * ( LST(i,j) - pta(i,j) )
      pSFLX_LH(i,j)  = LH0   * min(max(Uabs,U_minE),U_maxE) * rhos(i,j) * pQvEfc(i,j) * Ce * ( SatQvs - pQTRC(KS,i,j,I_QV) )
      pSFLX_GH(i,j)  = -2.0_RP * pTCS(i,j) * ( LST(i,j) - pTG(i,j)  ) / pDZg(i,j)
      pSFLX_SWU(i,j) = pALB(i,j) * pSWD(i,j)
      pSFLX_LWU(i,j) = pEMIT(i,j) * STB * LST(i,j)**4

    enddo
    enddo

    return
  end subroutine ts_known

  subroutine bulkcoef_uno( &
      Cm, Ch, Ce,                    & ! (out)
      pta, pts, za, uabs, z0, zt, ze ) ! (in)
    use mod_const, only: &
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
    real(RP), parameter :: tPrn = 0.74_RP    ! turbulent Prandtl number (Businger et al. 1971)
    real(RP), parameter :: LFb  = 9.4_RP     ! Louis factor b (Louis 1979)
    real(RP), parameter :: LFbp = 4.7_RP     ! Louis factor b' (Louis 1979)
    real(RP), parameter :: LFdm = 7.4_RP     ! Louis factor d for momemtum (Louis 1979)
    real(RP), parameter :: LFdh = 5.3_RP     ! Louis factor d for heat (Louis 1979)

    ! work
    real(RP) :: RiB, RiBT ! bulk Richardson number [no unit]
    real(RP) :: C0 ! initial drag coefficient [no unit]
    real(RP) :: fm, fh, t0th, q0qe

    C0    = ( KARMAN / log( za/z0 ) )**2
    RiBT  = GRAV * za * ( pta - pts ) / ( ( pta + pts ) * 0.5_RP * uabs**2 )
    RiB   = RiBT

    if( RiBT >= 0.0_RP ) then
      ! stable condition
      fm = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh = fm
    else
      ! unstable condition
      fm = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0 * sqrt( za/z0 ) * sqrt( abs( RiB ) ) )
      fh = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0 * sqrt( za/z0 ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + log( z0/zt ) / log( za/z0 ) / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + log( z0/ze ) / log( za/z0 ) / sqrt( fm ) * fh )
    RiB  = RiB * t0th

    if( RiBT >= 0.0_RP ) then
      ! stable condition
      fm = 1.0_RP / ( 1.0_RP + LFbp * RiB )**2
      fh = fm
    else
      ! unstable condition
      fm = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdm * C0 * sqrt( za/z0 ) * sqrt( abs( RiB ) ) )
      fh = 1.0_RP - LFb * RiB / ( 1.0_RP + LFb * LFdh * C0 * sqrt( za/z0 ) * sqrt( abs( RiB ) ) )
    end if

    t0th = 1.0_RP / ( 1.0_RP + log( z0/zt ) / log( za/z0 ) / sqrt( fm ) * fh )
    q0qe = 1.0_RP / ( 1.0_RP + log( z0/ze ) / log( za/z0 ) / sqrt( fm ) * fh )

    Cm = C0 * fm
    Ch = C0 * fh * t0th / tPrn
    Ce = C0 * fh * q0qe / tPrn

    Cm = min( max( Cm, Cm_min ), Cm_max )
    Ch = min( max( Ch, Ch_min ), Ch_max )
    Ce = min( max( Ce, Ce_min ), Ce_max )

    return
  end subroutine bulkcoef_uno

  subroutine sfcval_estimate( &
      sfc_rho, sfc_pre, & ! (out)
      rho, pre, zs      ) ! (in)
    use mod_const, only: &
      GRAV => CONST_GRAV
    use mod_grid, only: &
      CZ => GRID_CZ
    implicit none

    ! argument
    real(RP), intent(out) :: sfc_rho(IA,JA)    ! density at surface [kg/m3]
    real(RP), intent(out) :: sfc_pre(IA,JA)    ! pressure at surface [Pa]
    real(RP), intent(in)  :: rho    (KA,IA,JA) ! density [kg/m3]
    real(RP), intent(in)  :: pre    (KA,IA,JA) ! pressure [Pa]
    real(RP), intent(in)  :: zs     (IA,JA)    ! surface height [m]

    ! work
    integer :: i, j

    real(RP) :: zz
    real(RP) :: z1, z2, z3
    real(RP) :: p1, p2, p3
    real(RP) :: lag_intpl

    lag_intpl( zz, z1, p1, z2, p2, z3, p3 )                &
      = ( (zz-z2) * (zz-z3) ) / ( (z1-z2) * (z1-z3) ) * p1 &
      + ( (zz-z1) * (zz-z3) ) / ( (z2-z1) * (z2-z3) ) * p2 &
      + ( (zz-z1) * (zz-z2) ) / ( (z3-z1) * (z3-z2) ) * p3

    ! estimate surface density (extrapolation)
    do j = 1, JA
    do i = 1, IA
      sfc_rho(i,j) = lag_intpl( zs(i,j),                 &
                                CZ(KS  ), rho(KS  ,i,j), &
                                CZ(KS+1), rho(KS+1,i,j), &
                                CZ(KS+2), rho(KS+2,i,j)  )
    end do
    end do

    ! estimate surface pressure (hydrostatic balance)
    do j = 1, JA
    do i = 1, IA
      sfc_pre(i,j) = pre(KS,i,j)                             &
                   + 0.5_RP * ( sfc_rho(i,j) + rho(KS,i,j) ) &
                   * GRAV * ( CZ(KS) - zs(i,j) )
    end do
    end do

    return
  end subroutine sfcval_estimate

  function satmixr( temp, pres )
    use mod_const, only : &
      EPSvap => CONST_EPSvap
    implicit none

    ! argument
    real(RP), intent(in) :: temp ! temperature [K]
    real(RP), intent(in) :: pres ! pressure [Pa]
    ! function
    real(RP) :: satmixr ! saturated mixing ratio [kg/kg]

    satmixr = EPSvap * satvapor( temp ) / ( pres - satvapor( temp ) )

    return
  end function satmixr

  function satvapor( temp )
    use mod_const, only: &
      LH0   => CONST_LH0,   &
      Rvap  => CONST_Rvap,  &
      T00   => CONST_TEM00, &
      PSAT0 => CONST_PSAT0
    implicit none

    ! argument
    real(RP), intent(in) :: temp ! temperature [K]
    ! function
    real(RP) :: satvapor ! saturated water vapor pressure [Pa]

    ! Clasius-Clapeyron Equation
    satvapor = PSAT0 * exp( LH0/Rvap * ( 1.0_RP/T00 - 1.0_RP/temp ) )

    return
  end function satvapor

end module mod_cpl_atmos_land
