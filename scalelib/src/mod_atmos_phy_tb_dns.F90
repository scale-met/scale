!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulent process for DNS
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-03-30 (A.Noda)       [new] produced based on mod_atmos_phy_tb_smg.F90
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb_dns
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  use mod_tracer

#ifdef DEBUG
  use mod_debug, only: &
     CHECK
  use mod_const, only: &
     UNDEF => CONST_UNDEF, &
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
  real(RP), private, save :: nu=1.512d-5 ! [m2/s] kinematic viscosity coefficient for air at 20degC
! real(RP), private, save :: mu=1.8d-5   ! [m2/s] molecular diffusive coefficient for air at 20degC
  real(RP), private, save :: mu=1.512d-5   ! same as nu (needed based on hyposes. see Mellado 2010)

  !-----------------------------------------------------------------------------
contains

  subroutine ATMOS_PHY_TB_dns_setup( &
       TYPE_TB, &
       CDZ, CDX, CDY,   &
       FDZ, FDX, FDY,   &
       CZ, FZ )
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    character(len=H_SHORT), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: FDZ(KA-1)
    real(RP), intent(in) :: FDX(IA-1)
    real(RP), intent(in) :: FDY(JA-1)
    real(RP), intent(in) :: CZ(KA)
    real(RP), intent(in) :: FZ(0:KA)

    NAMELIST / PARAM_ATMOS_PHY_TB_DNS / &
         NU, &
         MU

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

  end subroutine ATMOS_PHY_TB_dns_setup

  subroutine ATMOS_PHY_TB_dns( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_qtrc,                & ! (out)
       tke, nu_C, Ri, Pr,                           & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC           ) ! (in)
    use mod_const, only: &
       GRAV => CONST_GRAV
    use mod_grid, only: &
       FDZ  => GRID_FDZ,  &
       FDX  => GRID_FDX,  &
       FDY  => GRID_FDY,  &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY
    implicit none

    ! SGS flux
    real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_qtrc(KA,IA,JA,QA,3)

    real(RP), intent(out) :: tke (KA,IA,JA) ! TKE
    real(RP), intent(out) :: nu_C(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Pr  (KA,IA,JA) ! Prantle number
    real(RP), intent(out) :: Ri  (KA,IA,JA) ! Richardson number

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP) :: POTT(KA,IA,JA)

    integer :: IIS, IIE
    integer :: JJS, JJE

    integer :: k, i, j, iq
    !---------------------------------------------------------------------------

#ifdef DEBUG
    POTT(:,:,:) = UNDEF
    qflx_sgs_momz(:,:,:,:) = UNDEF
    qflx_sgs_momx(:,:,:,:) = UNDEF
    qflx_sgs_momy(:,:,:,:) = UNDEF
    qflx_sgs_rhot(:,:,:,:) = UNDEF
    qflx_sgs_qtrc(:,:,:,:,:) = UNDEF
#endif


    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: SGS Parameterization'

    tke (:,:,:)=0.0
    nu_C(:,:,:)=0.0
    Pr  (:,:,:)=1.0
    Ri  (:,:,:)=0.0

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

!write(*,*)'chkdns1',nu,ks,iis,jjs
!write(*,*)'chkdns1-1-1',momz(ks+1,iis,jjs)
!write(*,*)'chkdns1-1-2',momz(ks,iis,jjs)
!write(*,*)'chkdns1-1-3',rfdz(ks+1)
!write(*,*)'chkdns1-1-4',qflx_sgs_momz(ks+1,iis,jjs,zdir)

       !##### momentum equation (z) #####
       ! (cell center)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS+1, KE-1
!write(*,*)'chkdns1-2',i,j,k,momz(k,i,j),momz(k-1,i,j),rfdz(k)
!         qflx_sgs_momz(k,i,j,ZDIR) = DENS(k,i,j) * ( &
!              - nu*( momz(k,i,j)-momz(k-1,i,j) )*rfdz(k)/dens(k,i,j)
          qflx_sgs_momz(k,i,j,ZDIR) = - nu*( momz(k,i,j)-momz(k-1,i,j) )*rfdz(k)
!write(*,*)'chkdns1-3',qflx_sgs_momz(k,i,j,zdir)
       enddo
       enddo
       enddo
!write(*,*)'chkdns2'
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momz(KS,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momz(KE,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
!write(*,*)'chkdns3'
       ! (y edge)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE-1
!         qflx_sgs_momz(k,i,j,XDIR) = - 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) &
!                              * nu*( momz(k,i+1,j)-momz(k,i,j) )*fdx(k) &
!                                  /( 0.25_RP*( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ))
          qflx_sgs_momz(k,i,j,XDIR) = -nu*( momz(k,i+1,j)-momz(k,i,j) )*fdx(i)
       enddo
       enddo
       enddo
!write(*,*)'chkdns4'
       ! (x edge)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE-1
!         qflx_sgs_momz(k,i,j,YDIR) = - 0.25_RP * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) &
!                              * nu * ( momz(k,i,j+1)-momz(k,i,j) )*fdy(j) &
!                              /( 0.25_RP*(DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) )
          qflx_sgs_momz(k,i,j,YDIR) = -nu*( momz(k,i,j+1)-momz(k,i,j) )*fdy(j)
       enddo
       enddo
       enddo

!write(*,*)'chkdns5'
       !##### momentum equation (x) #####
       ! (y edge)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
!         qflx_sgs_momx(k,i,j,ZDIR) = - 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) &
!                              * nu * ( momx(k+1,i,j)-momx(k,i,j) )*fdz(k) &
!                             /( 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k+1,i,j)+DENS(k+1,i+1,j) ) )          qflx_sgs_momx(k,i,j,ZDIR) = -nu*( momx(k+1,i,j)-momx(k,i,j) )*fdz(k)
       enddo
       enddo
       enddo
!write(*,*)'chkdns6'
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momx(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momx(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo
!write(*,*)'chkdns7'
       ! (cell center)
       do j = JJS, JJE
       do i = IIS, IIE+1
       do k = KS, KE
!         qflx_sgs_momx(k,i,j,XDIR) = -DENS(k,i,j) * nu*( momx(k,i,j)-momx(k,i-1,j) )*fdx(i-1) &
!              /dens(k,i,j)
          qflx_sgs_momx(k,i,j,XDIR) = -nu*( momx(k,i,j)-momx(k,i-1,j) )*fdx(i-1)
       enddo
       enddo
       enddo
!write(*,*)'chkdns8'
       ! (z edge)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
!         qflx_sgs_momx(k,i,j,YDIR) = - 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
!                              * nu*(momx(k,i,j+1)-momx(k,i,j))*fdy(j) &
!                              /( 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) )
          qflx_sgs_momx(k,i,j,YDIR) = -nu*(momx(k,i,j+1)-momx(k,i,j))*fdy(j)
       enddo
       enddo
       enddo

!write(*,*)'chkdns9'
       !##### momentum equation (y) #####
       ! (x edge)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
!         qflx_sgs_momy(k,i,j,ZDIR) = - 0.25_RP * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) &
!                              * nu * (momy(k+1,i,j)-momy(k,i,j))*fdz(k) &
!                              /( 0.25_RP * ( DENS(k,i,j)+DENS(k,i,j+1)+DENS(k+1,i,j)+DENS(k+1,i,j+1) ) )
          qflx_sgs_momy(k,i,j,ZDIR) = -nu*(momy(k+1,i,j)-momy(k,i,j))*fdz(k)
       enddo
       enddo
       enddo
!write(*,*)'chkdns10'
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momy(KS-1,i,j,ZDIR) = 0.0_RP ! bottom boundary
          qflx_sgs_momy(KE  ,i,j,ZDIR) = 0.0_RP ! top boundary
       enddo
       enddo

!write(*,*)'chkdns11'
       ! (z edge)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
!         qflx_sgs_momy(k,i,j,XDIR) = - 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) &
!                              * nu *(momy(k,i,j)-momy(k,i,j-1))*fdy(j-1) &
!                              /( 0.25_RP * ( DENS(k,i,j)+DENS(k,i+1,j)+DENS(k,i,j+1)+DENS(k,i+1,j+1) ) )
          qflx_sgs_momy(k,i,j,XDIR) = -nu*(momy(k,i+1,j)-momy(k,i,j))*fdx(i)
       enddo
       enddo
       enddo

!write(*,*)'chkdns12'
       ! (z-x plane)
       do j = JJS, JJE+1
       do i = IIS, IIE
       do k = KS, KE
!         qflx_sgs_momy(k,i,j,YDIR) = -DENS(k,i,j) * nu
!              *( momy(k,i,j)-momy(k,i,j-1))*fdy(j-1)/dens(k,i,j)
          qflx_sgs_momy(k,i,j,YDIR) = -nu*( momy(k,i,j)-momy(k,i,j-1))*fdy(j-1)
       enddo
       enddo
       enddo

!write(*,*)'chkdns13'
       !##### Thermodynamic Equation #####

       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_rhot(k,i,j,ZDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                               * mu &
                               * ( POTT(k+1,i,j)-POTT(k,i,j) ) * RFDZ(k)
       enddo
       enddo
       enddo
!write(*,*)'chkdns14'
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_rhot(KS-1,i,j,ZDIR) = 0.0_RP
          qflx_sgs_rhot(KE  ,i,j,ZDIR) = 0.0_RP
       enddo
       enddo

!write(*,*)'chkdns15'
       ! (y-z plane)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS, KE
          qflx_sgs_rhot(k,i,j,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * mu &
                               * ( POTT(k,i+1,j)-POTT(k,i,j) ) * RFDX(i)
       enddo
       enddo
       enddo
!write(*,*)'chkdns16'
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS, KE
          qflx_sgs_rhot(k,i,j,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * mu &
                               * ( POTT(k,i,j+1)-POTT(k,i,j) ) * RFDY(j)
       enddo
       enddo
       enddo

    enddo
    enddo
!write(*,*)'chkdns17'

    !##### Tracers #####
    do iq = 1, QA

    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       ! (x-y plane)
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE-1
          qflx_sgs_qtrc(k,i,j,iq,ZDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k+1,i,j) ) &
                               * mu * ( QTRC(k+1,i,j,iq)-QTRC(k,i,j,iq) ) * RFDZ(k)
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_qtrc(KS-1,i,j,iq,ZDIR) = 0.0_RP
          qflx_sgs_qtrc(KE  ,i,j,iq,ZDIR) = 0.0_RP
       enddo
       enddo

       ! (y-z plane)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_sgs_qtrc(k,i,j,iq,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * mu * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) * RFDX(i)
       enddo
       enddo
       enddo
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_sgs_qtrc(k,i,j,iq,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * mu * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) * RFDY(j)
       enddo
       enddo
       enddo
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_qtrc(KS-1,i,j,iq,ZDIR) = 0.0_RP
          qflx_sgs_qtrc(KE  ,i,j,iq,ZDIR) = 0.0_RP
       enddo
       enddo

       ! (y-z plane)
       do j = JJS,   JJE
       do i = IIS-1, IIE
       do k = KS,   KE
          qflx_sgs_qtrc(k,i,j,iq,XDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i+1,j) ) &
                               * mu * ( QTRC(k,i+1,j,iq)-QTRC(k,i,j,iq) ) * RFDX(i)
       enddo
       enddo
       enddo
       ! (z-x plane)
       do j = JJS-1, JJE
       do i = IIS,   IIE
       do k = KS,   KE
          qflx_sgs_qtrc(k,i,j,iq,YDIR) = - 0.5_RP * ( DENS(k,i,j)+DENS(k,i,j+1) ) &
                               * mu * ( QTRC(k,i,j+1,iq)-QTRC(k,i,j,iq) ) * RFDY(j)
       enddo
       enddo
       enddo

    enddo
    enddo

    enddo ! scalar quantities loop
!write(*,*)'chkdns18'

    return

  end subroutine ATMOS_PHY_TB_dns

end module mod_atmos_phy_tb_dns
