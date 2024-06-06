!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulent process for DNS
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_tb_dns
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer

#ifdef DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_TB_dns_config
  public :: ATMOS_PHY_TB_dns_setup
  public :: ATMOS_PHY_TB_dns

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
  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP

  real(RP), private :: ATMOS_PHY_TB_DNS_NU = 1.512E-5_RP ! [m2/s] kinematic viscosity coefficient for air at 20degC
! real(RP), private :: mu = 1.8E-5_RP   ! [m2/s] molecular diffusive coefficient for air at 20degC
  real(RP), private :: ATMOS_PHY_TB_DNS_MU = 1.512E-5_RP ! same as NU (needed based on hyposes. see Mellado 2010)

  logical,  private              :: twoD
  integer,  private              :: ISn, IEp

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_TB_dns_config( &
       TYPE_TB,  &
       I_TKE_out )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*), intent(in)  :: TYPE_TB
    integer,          intent(out) :: I_TKE_out
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_dns_config",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_dns_config",*) 'Tracers for Deardorff (1980) 1.5th TKE Model'

    if ( TYPE_TB /= 'DNS' ) then
       LOG_ERROR("ATMOS_PHY_TB_dns_config",*) 'ATMOS_PHY_TB_TYPE is not DNS. Check!'
       call PRC_abort
    endif

    I_TKE_out = -1

    return
  end subroutine ATMOS_PHY_TB_dns_config
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_dns_setup( &
       CDZ, CDX, CDY, CZ )
    use scale_prc, only: &
       PRC_abort
    implicit none

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)

    namelist / PARAM_ATMOS_PHY_TB_DNS / &
       ATMOS_PHY_TB_DNS_NU, &
       ATMOS_PHY_TB_DNS_MU

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_dns_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_dns_setup",*) 'Eddy Viscocity Model for DNS'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_DNS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_dns_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_dns_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_DNS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_DNS)

    twoD = ( IA == 1 )
    if ( twoD ) then
       ISn = 1
       IEp = 1
    else
       ISn = IS - 1
       IEp = IE + 1
    end if

    return
  end subroutine ATMOS_PHY_TB_dns_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_dns( &
       qflx_sgs_MOMZ, qflx_sgs_MOMX, qflx_sgs_MOMY, &
       qflx_sgs_rhot, qflx_sgs_rhoq,                &
       RHOQ_t, nu, Ri, Pr,                          &
       MOMZ, MOMX, MOMY, POTT, DENS, QTRC, N2,      &
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH, SFLX_Q,  &
       FZ, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY,     &
       GSQRT, J13G, J23G, J33G, MAPF, dt            )
    use scale_atmos_grid_cartesC_index
    use scale_tracer
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_atmos_phy_tb_common, only: &
       calc_strain_tensor => ATMOS_PHY_TB_calc_strain_tensor, &
       calc_flux_phi      => ATMOS_PHY_TB_calc_flux_phi       
    implicit none

    ! SGS flux
    real(RP), intent(out) :: qflx_sgs_MOMZ(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_MOMX(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_MOMY(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,3,QA)

    real(RP), intent(inout) :: RHOQ_t(KA,IA,JA,QA) ! tendency of rho * QTRC

    real(RP), intent(out) :: nu(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Ri(KA,IA,JA) ! Richardson number
    real(RP), intent(out) :: Pr(KA,IA,JA) ! Prantle number

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: POTT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)  :: N2(KA,IA,JA)

    real(RP), intent(in)  :: SFLX_MW(IA,JA)
    real(RP), intent(in)  :: SFLX_MU(IA,JA)
    real(RP), intent(in)  :: SFLX_MV(IA,JA)
    real(RP), intent(in)  :: SFLX_SH(IA,JA)
    real(RP), intent(in)  :: SFLX_Q (IA,JA,QA)

    real(RP), intent(in)  :: FZ           (0:KA,IA,JA)
    real(RP), intent(in)  :: FDZ          (KA-1)
    real(RP), intent(in)  :: RCDZ         (KA)
    real(RP), intent(in)  :: RFDZ         (KA-1)
    real(RP), intent(in)  :: CDX          (IA)
    real(RP), intent(in)  :: FDX          (IA-1)
    real(RP), intent(in)  :: CDY          (JA)
    real(RP), intent(in)  :: FDY          (JA-1)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(DP), intent(in)  :: dt

    ! deformation rate tensor
    real(RP) :: S33_C(KA,IA,JA) ! (cell center)
    real(RP) :: S11_C(KA,IA,JA)
    real(RP) :: S22_C(KA,IA,JA)
    real(RP) :: S31_C(KA,IA,JA)
    real(RP) :: S12_C(KA,IA,JA)
    real(RP) :: S23_C(KA,IA,JA)
    real(RP) :: S12_Z(KA,IA,JA) ! (z edge or x-y plane)
    real(RP) :: S23_X(KA,IA,JA) ! (x edge or y-z plane)
    real(RP) :: S31_Y(KA,IA,JA) ! (y edge or z-x plane)
    real(RP) :: S2   (KA,IA,JA) ! |S|^2

    real(RP) :: Kh(KA,IA,JA) ! eddy diffusion

    ! implicit scheme (dummy)
    real(RP) :: a   (KA,IA,JA)
    real(RP) :: b   (KA,IA,JA)
    real(RP) :: c   (KA,IA,JA)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------


    !$acc data copyout(qflx_sgs_MOMZ, qflx_sgs_MOMX, qflx_sgs_MOMY, qflx_sgs_rhot, qflx_sgs_rhoq, &
    !$acc              RHOQ_t,     &
    !$acc              nu, Ri, Pr) &
    !$acc      copyin(MOMZ, MOMX, MOMY, POTT, DENS, QTRC, &
    !$acc             FZ, FDZ, RCDZ, RFDZ, CDX, FDX, CDY, FDY, GSQRT, J13G, J23G, MAPF)  &
    !$acc      create(S33_C, S11_C, S22_C, S31_C, S12_C, S23_C, S12_Z, S23_X, S31_Y, S2, &
    !$acc             Kh, a, b, c)

    LOG_PROGRESS(*) 'atmosphere / physics / turbulence / DNS'

#ifdef DEBUG
    !$acc kernels
    qflx_sgs_MOMZ(:,:,:,:)   = UNDEF
    qflx_sgs_MOMX(:,:,:,:)   = UNDEF
    qflx_sgs_MOMY(:,:,:,:)   = UNDEF
    qflx_sgs_rhot(:,:,:,:)   = UNDEF
    qflx_sgs_rhoq(:,:,:,:,:) = UNDEF

    nu  (:,:,:) = UNDEF
    Ri  (:,:,:) = UNDEF
    Pr  (:,:,:) = UNDEF
    Ri  (:,:,:) = UNDEF
    Kh  (:,:,:) = UNDEF
    !$acc end kernels
#endif

    !$omp parallel workshare
    !$acc kernels
    nu (:,:,:) = 0.0_RP
    Ri (:,:,:) = 0.0_RP
    Pr (:,:,:) = 1.0_RP
    Kh (:,:,:) = ATMOS_PHY_TB_DNS_MU
    !$acc end kernels
    !$omp end parallel workshare 

    !##### Start Upadate #####

    call calc_strain_tensor( &
         S33_C, S11_C, S22_C,           & ! (out)
         S31_C, S12_C, S23_C,           & ! (out)
         S12_Z, S23_X, S31_Y,           & ! (out)
         S2                 ,           & ! (out)
         DENS, MOMZ, MOMX, MOMY,        & ! (in)
         GSQRT, J13G, J23G, J33G, MAPF, & ! (in)
         twoD                           ) ! (in)
        
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (z) #####

       !$omp parallel private(i,j,k)
       ! (cell center)

       if ( twoD ) then
         i = IIS
         !$omp do
         !$acc kernels
         do j = JJS, JJE
         do k = KS+1, KE-1
           qflx_sgs_MOMZ(k,i,j,ZDIR) = DENS(k,i,j) * ( &
             - 2.0_RP *  ATMOS_PHY_TB_DNS_NU * S33_C(k,i,j) &
             * ( S33_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) )
         enddo
         enddo
         !$acc end kernels
       else
         !$omp do collapse(2)
         !$acc kernels         
         do j = JJS, JJE
         do i = IIS, IIE
         do k = KS+1, KE-1
           qflx_sgs_MOMZ(k,i,j,ZDIR) = DENS(k,i,j) * ( &
             - 2.0_RP *  ATMOS_PHY_TB_DNS_NU * S33_C(k,i,j) &
             * ( S33_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) )
         enddo
         enddo
         enddo
         !$acc end kernels
       end if 

       !$omp do
       !$acc kernels
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_MOMZ(KS,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_MOMZ(KE,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       ! (y edge)
       if ( .not. twoD ) then
       !$omp do collapse(2)
       !$acc kernels  
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
         qflx_sgs_MOMZ(k,i,j,XDIR) = - 0.5_RP & ! 2.0 / 4
            * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
            * ATMOS_PHY_TB_DNS_NU                                         &
            * S31_Y(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait
       end if

       ! (x edge)
       !$omp do collapse(2)
       !$acc kernels
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
         qflx_sgs_MOMZ(k,i,j,YDIR) = - 0.5_RP & ! 2.0 / 4
            * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
            * ATMOS_PHY_TB_DNS_NU                                         &
            * S23_X(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       !$omp end parallel

       !##### momentum equation (x) #####

       if ( .not. twoD ) then

       !$omp parallel private(i,j,k)
  
       ! (y edge)

       !$omp do collapse(2)
       !$acc kernels
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
         qflx_sgs_MOMX(k,i,j,ZDIR) = - 0.5_RP & ! 2/4
              * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i+1,j)+DENS(k+1,i+1,j) ) &
              * ATMOS_PHY_TB_DNS_NU                                         &
              * S31_Y(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       !$omp do
       !$acc kernels
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_MOMX(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_MOMX(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       ! (cell center)
       !$omp do collapse(2)
       !$acc kernels
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
          qflx_sgs_MOMX(k,i,j,XDIR) = DENS(k,i,j) * ( &
               - 2.0_RP * ATMOS_PHY_TB_DNS_NU  &
               * ( S11_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) )
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       ! (z edge)
       !$omp do collapse(2)
       !$acc kernels
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_sgs_MOMX(k,i,j,YDIR) = - 0.5_RP & ! 2/4
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ATMOS_PHY_TB_DNS_NU                                         &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       !$omp end parallel

       end if ! twoD

       !##### momentum equation (y) #####

       !$omp parallel private(i,j,k)

       ! (x edge)
       !$omp do collapse(2)
       !$acc kernels       
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_MOMY(k,i,j,ZDIR) = - 0.5_RP & ! 2/4
          * ( DENS(k,i,j)+DENS(k+1,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j+1) ) &
          * ATMOS_PHY_TB_DNS_NU                                         &
          * S23_X(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       !$omp do
       !$acc kernels             
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_MOMY(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_MOMY(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait

       ! (z edge)
       if ( .not. twoD ) then
      !$omp do collapse(2)
      !$acc kernels  
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_sgs_MOMY(k,i,j,XDIR) = - 0.5_RP & ! 2/4
               * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
               * ATMOS_PHY_TB_DNS_NU                                         &
               * S12_Z(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
       !$omp end do nowait
       end if

       ! (z-x plane)
       if ( twoD ) then
         i = IIS
         !$omp do
         !$acc kernels
         do j = JJS, JJE+1
         do k = KS, KE
             qflx_sgs_MOMY(k,i,j,YDIR) = DENS(k,i,j) * ( &
                  - 2.0_RP * ATMOS_PHY_TB_DNS_NU &
                  * ( S22_C(k,i,j) - ( S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) )
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       else
          !$omp do collapse(2)
          !$acc kernels
          do j = JJS, JJE+1
          do i = IIS, IIE
          do k = KS, KE
             qflx_sgs_MOMY(k,i,j,YDIR) = DENS(k,i,j) * ( &
                  - 2.0_RP * ATMOS_PHY_TB_DNS_NU &
                  * ( S22_C(k,i,j) - ( S11_C(k,i,j) + S22_C(k,i,j) + S33_C(k,i,j) ) * OneOverThree ) )
          enddo
          enddo
          enddo
          !$acc end kernels
          !$omp end do nowait
       end if

       !$omp end parallel

       !##### Thermodynamic Equation #####

       call calc_flux_phi( &
            qflx_sgs_rhot,                 & ! [inout]
            DENS, POTT, Kh, 1.0_RP,        & ! [in]
            GSQRT, J13G, J23G, J33G, MAPF, & ! [in]
            .false., .false.,              & ! horizontal, implicit
            a, b, c, dt,                   & ! [in, dummy]
            IIS, IIE, JJS, JJE,            &
            twoD )
    enddo
    enddo

    !##### Tracers #####
    do iq = 1, QA

      if ( .not. TRACER_ADVC(iq) ) cycle

      do JJS = JS, JE, JBLOCK
      JJE = JJS+JBLOCK-1
      do IIS = IS, IE, IBLOCK
      IIE = IIS+IBLOCK-1

      call calc_flux_phi( &
         qflx_sgs_rhoq(:,:,:,:,iq),        & ! [inout]
         DENS, QTRC(:,:,:,iq), Kh, 1.0_RP, & ! [in]
         GSQRT, J13G, J23G, J33G, MAPF,    & ! [in]
         .false., .false.,                 & ! horizontal, implicit
         a, b, c, dt,                      & ! [in, dummy]
         IIS, IIE, JJS, JJE,               &
         twoD )
      
      enddo
      enddo

    enddo ! scalar quantities loop

    !$acc end data

    return
  end subroutine ATMOS_PHY_TB_dns

end module scale_atmos_phy_tb_dns
