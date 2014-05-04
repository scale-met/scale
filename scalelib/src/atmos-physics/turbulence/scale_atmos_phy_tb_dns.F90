!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulent process for DNS
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-03-30 (A.Noda)       [new] produced based on scale_atmos_phy_tb_smg.F90
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb_dns
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer

#ifdef DEBUG
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
  real(RP), private, save :: ATMOS_PHY_TB_DNS_NU = 1.512E-5_RP ! [m2/s] kinematic viscosity coefficient for air at 20degC
! real(RP), private, save :: mu = 1.8E-5_RP   ! [m2/s] molecular diffusive coefficient for air at 20degC
  real(RP), private, save :: ATMOS_PHY_TB_DNS_MU = 1.512E-5_RP ! same as NU (needed based on hyposes. see Mellado 2010)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_dns_setup( &
       TYPE_TB,       &
       CDZ, CDX, CDY, &
       CZ             )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA)

    NAMELIST / PARAM_ATMOS_PHY_TB_DNS / &
         ATMOS_PHY_TB_DNS_NU, &
         ATMOS_PHY_TB_DNS_MU

    integer :: k, i, j

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '+++ Eddy Viscocity Model for DNS'

    if ( TYPE_TB /= 'DNS' ) then
       if ( IO_L ) write(IO_FID_LOG,*) 'xxx ATMOS_PHY_TB_TYPE is not DNS. Check!'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_DNS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB_DNS. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB_DNS)

    return
  end subroutine ATMOS_PHY_TB_dns_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_dns( &
       qflx_sgs_MOMZ, qflx_sgs_MOMX, qflx_sgs_MOMY, &
       qflx_sgs_rhot, qflx_sgs_qtrc,                &
       tke, nu, Ri, Pr,                             &
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          &
       GSQRT, J13G, J23G, J33G                      )
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_grid, only: &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ
    implicit none

    ! SGS flux
    real(RP), intent(out) :: qflx_sgs_MOMZ(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_MOMX(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_MOMY(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_qtrc(KA,IA,JA,QA,3)

    real(RP), intent(out) :: tke(KA,IA,JA) ! TKE
    real(RP), intent(out) :: nu (KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Pr (KA,IA,JA) ! Prantle number
    real(RP), intent(out) :: Ri (KA,IA,JA) ! Richardson number

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix

    real(RP) :: POTT(KA,IA,JA)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

#ifdef DEBUG
    qflx_sgs_MOMZ(:,:,:,:)   = UNDEF
    qflx_sgs_MOMX(:,:,:,:)   = UNDEF
    qflx_sgs_MOMY(:,:,:,:)   = UNDEF
    qflx_sgs_rhot(:,:,:,:)   = UNDEF
    qflx_sgs_qtrc(:,:,:,:,:) = UNDEF

    POTT(:,:,:) = UNDEF
#endif

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: SGS Parameterization'

    tke(:,:,:) = 0.0_RP
    nu (:,:,:) = 0.0_RP
    Pr (:,:,:) = 1.0_RP
    Ri (:,:,:) = 0.0_RP

    ! potential temperature
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE
       POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo

    !##### Start Upadate #####
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       !##### momentum equation (z) #####
       ! (cell center)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
          qflx_sgs_MOMZ(k,i,j,ZDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMZ(k,i,j)-MOMZ(k-1,i,j) ) * RCDZ(k)
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_MOMZ(KS,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_MOMZ(KE,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo

       ! (y edge)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
          qflx_sgs_MOMZ(k,i,j,XDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMZ(k,i+1,j)-MOMZ(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo

       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
          qflx_sgs_MOMZ(k,i,j,YDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMZ(k,i,j+1)-MOMZ(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo

       !##### momentum equation (x) #####
       ! (y edge)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_MOMX(k,i,j,ZDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMX(k+1,i,j)-MOMX(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_MOMX(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_MOMX(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo

       ! (cell center)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
          qflx_sgs_MOMX(k,i,j,XDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMX(k,i,j)-MOMX(k,i-1,j) ) * RCDX(i)
       enddo
       enddo
       enddo

       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_sgs_MOMX(k,i,j,YDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMX(k,i,j+1)-MOMX(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo

       !##### momentum equation (y) #####

       ! (x edge)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_MOMY(k,i,j,ZDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMY(k+1,i,j)-MOMY(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_MOMY(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_MOMY(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo

       ! (z edge)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_sgs_MOMY(k,i,j,XDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMY(k,i+1,j)-MOMY(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo

       ! (z-x plane)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
          qflx_sgs_MOMY(k,i,j,YDIR) = -ATMOS_PHY_TB_DNS_NU * ( MOMY(k,i,j)-MOMY(k,i,j-1) ) * RCDY(j)
       enddo
       enddo
       enddo

       !##### Thermodynamic Equation #####

       ! at x, y ,w
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_rhot(k,i,j,ZDIR) = -0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) ) &
                                    * ATMOS_PHY_TB_DNS_MU * ( POTT(k+1,i,j)-POTT(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo

       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_rhot(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_sgs_rhot(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo

       ! at u, y, z
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_sgs_rhot(k,i,j,XDIR) = -0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                                    * ATMOS_PHY_TB_DNS_MU * ( POTT(k,i+1,j)-POTT(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo

       ! at x, v, z
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_sgs_rhot(k,i,j,YDIR) = -0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
                                    * ATMOS_PHY_TB_DNS_MU * ( POTT(k,i,j+1)-POTT(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo

    enddo
    enddo

    !##### Tracers #####
    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! at x, y ,w
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_qtrc(k,i,j,iq,ZDIR) = -0.5_RP * ( DENS(k+1,i,j)+DENS(k,i,j) ) &
                                       * ATMOS_PHY_TB_DNS_MU * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) * RFDZ(k)
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_qtrc(KS-1,i,j,iq,ZDIR) = 0.0_RP
          qflx_sgs_qtrc(KE  ,i,j,iq,ZDIR) = 0.0_RP
       enddo
       enddo

       ! at u, y, z
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_sgs_qtrc(k,i,j,iq,XDIR) = -0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                                       * ATMOS_PHY_TB_DNS_MU * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) * RFDX(i)
       enddo
       enddo
       enddo

       ! at x, v, z
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_sgs_qtrc(k,i,j,iq,YDIR) = -0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) ) &
                                       * ATMOS_PHY_TB_DNS_MU * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) * RFDY(j)
       enddo
       enddo
       enddo

    enddo
    enddo

    enddo ! scalar quantities loop

    return
  end subroutine ATMOS_PHY_TB_dns

end module scale_atmos_phy_tb_dns
