!-------------------------------------------------------------------------------
!> module Saturation Adjustment
!!
!! @par Description
!!          Saturation adjustment module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-10-24 (T.Seiki)   [new] Import from NICAM
!! @li      2012-02-10 (H.Yashiro) [mod] Reconstruction
!!
!<
!-------------------------------------------------------------------------------
module mod_satadjust
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_const, only : &
     Rdry   => CONST_Rdry,   &
     CPdry  => CONST_CPdry,  &
     CVdry  => CONST_CVdry,  &
     Rvap   => CONST_Rvap,   &
     CPvap  => CONST_CPvap,  &
     CL     => CONST_CL,     &
     CI     => CONST_CI,     &
     LH00   => CONST_LH00,   &
     LHS00  => CONST_LHS00,  &
     LH0    => CONST_LH0,    &
     LHS0   => CONST_LHS0,   &
     EPSvap => CONST_EPSvap, &
     PSAT0  => CONST_PSAT0,  &
     TEM00  => CONST_TEM00
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  use mod_grid, only: &
     KA => GRID_KA, &
     IA => GRID_IA, &
     JA => GRID_JA, &
     KS => GRID_KS, &
     KE => GRID_KE, &
     IS => GRID_IS, &
     IE => GRID_IE, &
     JS => GRID_JS, &
     JE => GRID_JE, &
     IJA => GRID_IJA, &
     IJS => GRID_IJS, &
     IJE => GRID_IJE
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: moist_psat_water
  public :: moist_psat_ice
  public :: moist_qsat_water
  public :: moist_qsat_ice

  public :: moist_dqsw_dtem_rho  
  public :: moist_dqsi_dtem_rho  
  public :: moist_dqsw_dtem_dpre 
  public :: moist_dqsi_dtem_dpre 

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
  real(8), private, parameter :: TEM_MIN   = 10.D0

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  ! psat : Clasius-Clapeyron: based on CPV, CPL constant
  !-----------------------------------------------------------------------------
  subroutine moist_psat_water ( temp, psat )
    implicit none

    real(8), intent(in)  :: temp(IJA,KA)
    real(8), intent(out) :: psat(IJA,KA)

    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_psat_water")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE
       TEM = max( temp(ij,k), TEM_MIN )

       psat(ij,k) = PSAT0                                 &
                  * ( TEM * RTEM00 )**CPovRvap            &
                  * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_psat_water")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_psat_water

  !-----------------------------------------------------------------------------
  subroutine moist_psat_ice ( temp, psat )
    implicit none

    real(8), intent(in)  :: temp(IJA,KA)
    real(8), intent(out) :: psat(IJA,KA)

    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_psat_ice")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE
       TEM = max( temp(ij,k), TEM_MIN )

       psat(ij,k) = PSAT0                                   &
                   * ( TEM * RTEM00 )**CPovRvap              &
                   * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_psat_ice")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_psat_ice

  !-----------------------------------------------------------------------------
  subroutine moist_qsat_water ( temp, pres, qsat )
    implicit none

    real(8), intent(in)  :: temp(IJA,KA)
    real(8), intent(in)  :: pres(IJA,KA)
    real(8), intent(out) :: qsat(IJA,KA)
    
    real(8) :: psat
    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_qsat_water")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE
       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

       qsat(ij,k) = EPSvap * psat / ( pres(ij,k) - ( 1.D0-EPSvap ) * psat )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_qsat_water")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_qsat_water

  !-----------------------------------------------------------------------------
  subroutine moist_qsat_ice ( temp, pres, qsat )
    implicit none

    real(8), intent(in)  :: temp(IJA,KA)
    real(8), intent(in)  :: pres(IJA,KA)
    real(8), intent(out) :: qsat(IJA,KA)
    
    real(8) :: psat
    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_qsat_ice")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE

       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

       qsat(ij,k) = EPSvap * psat / ( pres(ij,k) - ( 1.D0-EPSvap ) * psat )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_qsat_ice")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_qsat_ice

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  !-----------------------------------------------------------------------------
  subroutine moist_dqsw_dtem_rho( temp, dens, dqsdtem )
    implicit none

    real(8), intent(in)  :: temp   (IJA,KA)
    real(8), intent(in)  :: dens   (IJA,KA)
    real(8), intent(out) :: dqsdtem(IJA,KA)

    real(8) :: psat ! saturation vapor pressure
    real(8) :: lhv  ! latent heat for condensation

    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsw_dtem_rho")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE

       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

       lhv  = LH0 + ( CPvap-CL ) * ( temp(ij,k)-TEM00 )

       dqsdtem(ij,k) = psat / ( dens(ij,k) * Rvap * temp(ij,k) * temp(ij,k) ) &
                      * ( lhv / ( Rvap * temp(ij,k) ) - 1.D0 )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsw_dtem_rho")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsw_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  !-----------------------------------------------------------------------------
  subroutine moist_dqsi_dtem_rho( temp, dens, dqsdtem )
    implicit none

    real(8), intent(in)  :: temp   (IJA,KA)
    real(8), intent(in)  :: dens   (IJA,KA)
    real(8), intent(out) :: dqsdtem(IJA,KA)

    real(8) :: psat ! saturation vapor pressure
    real(8) :: lhv  ! latent heat for condensation

    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsi_dtem_rho")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE

       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

       lhv  = LHS0 + ( CPvap-CI ) * ( temp(ij,k)-TEM00 )

       dqsdtem(ij,k) = psat / ( dens(ij,k) * Rvap * temp(ij,k) * temp(ij,k) ) &
                      * ( lhv / ( Rvap * temp(ij,k) ) - 1.D0 )

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsi_dtem_rho")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsi_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine moist_dqsw_dtem_dpre( temp, pres, dqsdtem, dqsdpre )
    implicit none

    real(8), intent(in)  :: temp   (IJA,KA)
    real(8), intent(in)  :: pres   (IJA,KA)
    real(8), intent(out) :: dqsdtem(IJA,KA)
    real(8), intent(out) :: dqsdpre(IJA,KA)

    real(8) :: psat ! saturation vapor pressure
    real(8) :: lhv  ! latent heat for condensation

    real(8) :: den1, den2 ! denominator
    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsw_dtem_dpre")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CL ) / Rvap
    LHovRvap = LH00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE

       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

       lhv  = LH0 + ( CPvap-CL ) * ( temp(ij,k)-TEM00 )

       den1 = ( pres(ij,k) - (1.D0-EPSvap) * psat ) &
            * ( pres(ij,k) - (1.D0-EPSvap) * psat )
       den2 = den1 * Rvap * temp(ij,k) * temp(ij,k)

       dqsdpre(ij,k) = - EPSvap * psat / den1
       dqsdtem(ij,k) =   EPSvap * psat / den2 * lhv * pres(ij,k)

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsw_dtem_dpre")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsw_dtem_dpre

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  !-----------------------------------------------------------------------------
  subroutine moist_dqsi_dtem_dpre( temp, pres, dqsdtem, dqsdpre )
    implicit none

    real(8), intent(in)  :: temp   (IJA,KA)
    real(8), intent(in)  :: pres   (IJA,KA)
    real(8), intent(out) :: dqsdtem(IJA,KA)
    real(8), intent(out) :: dqsdpre(IJA,KA)

    real(8) :: psat ! saturation vapor pressure
    real(8) :: lhv  ! latent heat for condensation

    real(8) :: den1, den2 ! denominator
    real(8) :: RTEM00, CPovRvap, LHovRvap, TEM

    integer :: ij, k
    !---------------------------------------------------------------------------

    call TIME_rapstart('satadjust')
#ifdef _FPCOLL_
call START_COLLECTION("moist_dqsi_dtem_dpre")
#endif

    RTEM00   = 1.D0 / TEM00
    CPovRvap = ( CPvap - CI ) / Rvap
    LHovRvap = LHS00 / Rvap

    do k  = KS,  KE
    do ij = IJS, IJE

       TEM = max( temp(ij,k), TEM_MIN )

       psat = PSAT0                                   &
            * ( TEM * RTEM00 )**CPovRvap              &
            * exp( LHovRvap * ( RTEM00 - 1.D0/TEM ) )

       lhv  = LHS0 + ( CPvap-CI ) * ( temp(ij,k)-TEM00 )

       den1 = ( pres(ij,k) - (1.D0-EPSvap) * psat ) &
            * ( pres(ij,k) - (1.D0-EPSvap) * psat )
       den2 = den1 * Rvap * temp(ij,k) * temp(ij,k)

       dqsdpre(ij,k) = - EPSvap * psat / den1
       dqsdtem(ij,k) =   EPSvap * psat / den2 * lhv * pres(ij,k)

    enddo
    enddo

#ifdef _FPCOLL_
call STOP_COLLECTION("moist_dqsi_dtem_dpre")
#endif
    call TIME_rapend  ('satadjust')

    return
  end subroutine moist_dqsi_dtem_dpre

end module mod_satadjust
!-------------------------------------------------------------------------------
