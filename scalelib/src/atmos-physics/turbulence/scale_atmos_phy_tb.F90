!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa) [new]
!! @li      2014-03-30 (A.Noda)      [mod] add DNS
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb
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
  abstract interface
     subroutine tb( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, &
       qflx_sgs_rhot, qflx_sgs_rhoq,                &
       RHOQ_t,                                      &
       nu_C, Ri, Pr, N2,                            &
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_QV, &
       GSQRT, J13G, J23G, J33G, MAPF, dt            )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none

       real(RP), intent(out)   :: qflx_sgs_momz(KA,IA,JA,3)
       real(RP), intent(out)   :: qflx_sgs_momx(KA,IA,JA,3)
       real(RP), intent(out)   :: qflx_sgs_momy(KA,IA,JA,3)
       real(RP), intent(out)   :: qflx_sgs_rhot(KA,IA,JA,3)
       real(RP), intent(out)   :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

       real(RP), intent(inout) :: RHOQ_t       (KA,IA,JA,QA) ! tendency of rho * QTRC

       real(RP), intent(out)   :: nu_C         (KA,IA,JA)    ! eddy viscosity (center)
       real(RP), intent(out)   :: Ri           (KA,IA,JA)    ! Richardson number
       real(RP), intent(out)   :: Pr           (KA,IA,JA)    ! Prantle number
       real(RP), intent(out)   :: N2           (KA,IA,JA)    ! squared Brunt-Vaisala frequency

       real(RP), intent(in)    :: MOMZ         (KA,IA,JA)
       real(RP), intent(in)    :: MOMX         (KA,IA,JA)
       real(RP), intent(in)    :: MOMY         (KA,IA,JA)
       real(RP), intent(in)    :: RHOT         (KA,IA,JA)
       real(RP), intent(in)    :: DENS         (KA,IA,JA)
       real(RP), intent(in)    :: QTRC         (KA,IA,JA,QA)

       real(RP), intent(in)    :: SFLX_MW      (IA,JA)
       real(RP), intent(in)    :: SFLX_MU      (IA,JA)
       real(RP), intent(in)    :: SFLX_MV      (IA,JA)
       real(RP), intent(in)    :: SFLX_SH      (IA,JA)
       real(RP), intent(in)    :: SFLX_QV      (IA,JA)

       real(RP), intent(in)    :: GSQRT        (KA,IA,JA,7)  !< vertical metrics {G}^1/2
       real(RP), intent(in)    :: J13G         (KA,IA,JA,7)  !< (1,3) element of Jacobian matrix
       real(RP), intent(in)    :: J23G         (KA,IA,JA,7)  !< (1,3) element of Jacobian matrix
       real(RP), intent(in)    :: J33G                       !< (3,3) element of Jacobian matrix
       real(RP), intent(in)    :: MAPF         (IA,JA,2,4)   !< map factor
       real(DP), intent(in)    :: dt
     end subroutine tb

     subroutine su( &
          CDZ, CDX, CDY, CZ )
       use scale_precision
       use scale_grid_index
       use scale_tracer
       implicit none

       real(RP), intent(in) :: CDZ(KA)
       real(RP), intent(in) :: CDX(IA)
       real(RP), intent(in) :: CDY(JA)
       real(RP), intent(in) :: CZ (KA,IA,JA)
     end subroutine su
  end interface

  procedure(tb), pointer :: ATMOS_PHY_TB       => NULL()

  procedure(su), pointer :: ATMOS_PHY_TB_setup => NULL()

  public :: ATMOS_PHY_TB_config
  public :: ATMOS_PHY_TB
  public :: ATMOS_PHY_TB_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: I_TKE = -1

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
  !> register
  subroutine ATMOS_PHY_TB_config( TB_TYPE )
    use scale_process, only: &
       PRC_MPIstop
#define EXTM(pre, name, post) pre ## name ## post
#define NAME(pre, name, post) EXTM(pre, name, post)
#ifdef TB
    use NAME(scale_atmos_phy_tb_, TB,), only: &
       NAME(ATMOS_PHY_TB_, TB, _config), &
       NAME(ATMOS_PHY_TB_, TB, _setup), &
       NAME(ATMOS_PHY_TB_, TB,)
#else
    use scale_atmos_phy_tb_smg, only: &
       ATMOS_PHY_TB_smg_config, &
       ATMOS_PHY_TB_smg_setup, &
       ATMOS_PHY_TB_smg
    use scale_atmos_phy_tb_d1980, only: &
       ATMOS_PHY_TB_d1980_config, &
       ATMOS_PHY_TB_d1980_setup, &
       ATMOS_PHY_TB_d1980
    use scale_atmos_phy_tb_dns, only: &
       ATMOS_PHY_TB_dns_config, &
       ATMOS_PHY_TB_dns_setup, &
       ATMOS_PHY_TB_dns
    use scale_atmos_phy_tb_mynn, only: &
       ATMOS_PHY_TB_mynn_config, &
       ATMOS_PHY_TB_mynn_setup, &
       ATMOS_PHY_TB_mynn
    use scale_atmos_phy_tb_hybrid, only: &
       ATMOS_PHY_TB_hybrid_config, &
       ATMOS_PHY_TB_hybrid_setup, &
       ATMOS_PHY_TB_hybrid
#endif
    implicit none

    character(len=*), intent(in) :: TB_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** => ', trim(TB_TYPE), ' is selected.'

#ifdef TB
    call NAME(ATMOS_PHY_TB_, TB, _config)( TB_TYPE, & ! [IN]
                                           I_TKE    ) ! [OUT]

    ATMOS_PHY_TB        => NAME(ATMOS_PHY_TB_, TB, )
    ATMOS_PHY_TB_setup  => NAME(ATMOS_PHY_TB_, TB, _setup)
#else
    select case( TB_TYPE )
    case( 'SMAGORINSKY' )

       call ATMOS_PHY_TB_smg_config( TB_TYPE, & ! [IN]
                                     I_TKE    ) ! [OUT]

       ATMOS_PHY_TB       => ATMOS_PHY_TB_smg
       ATMOS_PHY_TB_setup => ATMOS_PHY_TB_smg_setup

    case( 'D1980' )

       call ATMOS_PHY_TB_d1980_config( TB_TYPE, & ! [IN]
                                       I_TKE    ) ! [OUT]

       ATMOS_PHY_TB       => ATMOS_PHY_TB_d1980
       ATMOS_PHY_TB_setup => ATMOS_PHY_TB_d1980_setup

    case( 'DNS' )

       call ATMOS_PHY_TB_dns_config( TB_TYPE, & ! [IN]
                                     I_TKE    ) ! [OUT]

       ATMOS_PHY_TB       => ATMOS_PHY_TB_dns
       ATMOS_PHY_TB_setup => ATMOS_PHY_TB_dns_setup

    case( 'MYNN' )

       call ATMOS_PHY_TB_mynn_config( TB_TYPE, & ! [IN]
                                      I_TKE    ) ! [OUT]

       ATMOS_PHY_TB       => ATMOS_PHY_TB_mynn
       ATMOS_PHY_TB_setup => ATMOS_PHY_TB_mynn_setup

    case('HYBRID')

       call ATMOS_PHY_TB_hybrid_config( TB_TYPE, & ! [IN]
                                        I_TKE    ) ! [OUT]

       ATMOS_PHY_TB       => ATMOS_PHY_TB_hybrid
       ATMOS_PHY_TB_setup => ATMOS_PHY_TB_hybrid_setup

    case('OFF')

       ! do nothing

    case default

       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is invalid: ', TB_TYPE
       call PRC_MPIstop

    end select
#endif

    return
  end subroutine ATMOS_PHY_TB_config

end module scale_atmos_phy_tb
