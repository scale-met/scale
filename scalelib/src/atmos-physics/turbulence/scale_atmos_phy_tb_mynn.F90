!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Boundary layer turbulence model
!!          Mellor-Yamada Nakanishi-Niino model
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-08-27 (S.Nishizawa) [new]
!!
!! - Reference
!!  - Nakanishi and Niino, 2009:
!!    Development of an improved turbulence closure model for the atmospheric boundary layer.
!!    J. Meteorol. Soc. Japan, 87, 895-912
!!
!<
!-------------------------------------------------------------------------------
module scale_atmos_phy_tb_mynn
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
  public :: ATMOS_PHY_TB_mynn_setup
  public :: ATMOS_PHY_TB_mynn

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
  real(RP), parameter :: OneOverThree = 1.0_RP / 3.0_RP
  real(RP), parameter :: LT_min = 1.0E-6_RP
  real(RP), parameter :: Us_min = 1.0E-6_RP

  real(RP)              :: A1
  real(RP)              :: A2
  real(RP), parameter   :: B1 = 24.0_RP
  real(RP), parameter   :: B2 = 15.0_RP
  real(RP)              :: C1
  real(RP), parameter   :: C2 = 0.75_RP
  real(RP), parameter   :: C3 = 0.352_RP
  real(RP), parameter   :: C5 = 0.2_RP
  real(RP), parameter   :: G1 = 0.235_RP
  real(RP)              :: G2
  real(RP)              :: F1
  real(RP)              :: F2
  real(RP)              :: Rf1
  real(RP)              :: Rf2
  real(RP)              :: Rfc
  real(RP)              :: AF12 !< A1 F1 / A2 F2
  real(RP), parameter   :: PrN = 0.74_RP

  integer :: KE_PBL
  logical :: ATMOS_PHY_TB_MYNN_TKE_INIT = .false. !< set tke with that of level 2 at the first time if .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_mynn_setup( &
       TYPE_TB,       &
       CDZ, CDX, CDY, &
       CZ             )
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    character(len=*), intent(in) :: TYPE_TB

    real(RP), intent(in) :: CDZ(KA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: CZ (KA,IA,JA)

    integer ATMOS_PHY_TB_MYNN_KMAX_PBL !< number of layers for planetary boundary layer (must be equal or less than KMAX)

    NAMELIST / PARAM_ATMOS_PHY_TB_MYNN / &
         ATMOS_PHY_TB_MYNN_TKE_INIT, &
         ATMOS_PHY_TB_MYNN_KMAX_PBL

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[TURBULENCE] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Mellor-Yamada Nakanishi-Niino Model'

    if ( TYPE_TB /= 'MYNN' ) then
       write(*,*) 'xxx ATMOS_PHY_TB_TYPE is not MYNN. Check!'
       call PRC_MPIstop
    endif


    ATMOS_PHY_TB_MYNN_KMAX_PBL = KMAX - 1

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_MYNN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_TB_MYNN. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_TB_MYNN)

    if ( ATMOS_PHY_TB_MYNN_KMAX_PBL > KMAX-1 ) then
       write(*,*) 'xxx ATMOS_PHY_TB_MYNN_KMAX_PBL must be equal or less than KMAX-1'
       call PRC_MPIstop
    end if
    KE_PBL = ATMOS_PHY_TB_MYNN_KMAX_PBL + KHALO


    A1 = B1 * (1.0_RP - 3.0_RP * G1) / 6.0_RP
    A2 = 1.0_RP / (3.0_RP * G1 * B1**(1.0_RP/3.0_RP) * PrN )
    C1 = G1 - 1.0_RP / ( 3.0_RP * A1 * B1**(1.0_RP/3.0_RP) )
    G2 = ( 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + B2 * (1.0_RP - C3) ) / B1
    F1 = B1 * (G1 - C1) + 2.0_RP * A1 * (3.0_RP - 2.0_RP * C2) + 3.0_RP * A2 * (1.0_RP - C2) * (1.0_RP - C5)
    F2 = B1 * (G1 + G2) - 3.0_RP * A1 * (1.0_RP - C2)

    Rf1 = B1 * (G1 - C1) / F1
    Rf2 = B1 * G1 / F2
    Rfc = G1 / (G1 + G2)

    AF12 = A1 * F1 / ( A2 * F2 )

    return
  end subroutine ATMOS_PHY_TB_mynn_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_mynn( &
       qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
       qflx_sgs_rhot, qflx_sgs_rhoq,                & ! (out)
       tke,                                         & ! (inout)
       Nu, Ri, Pr,                                  & ! (out) diagnostic variables
       MOMZ, MOMX, MOMY, RHOT, DENS, QTRC,          & ! (in)
       SFLX_MW, SFLX_MU, SFLX_MV, SFLX_SH,          & ! (in)
       GSQRT, J13G, J23G, J33G, MAPF, dt            ) ! (in)
    use scale_precision
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       EPS    => CONST_EPS
    use scale_grid, only: &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ,  &
       RCDZ => GRID_RCDZ, &
       RFDZ => GRID_RFDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY
    use scale_grid_real, only: &
       FZ => REAL_FZ
    use scale_atmos_refstate, only: &
       PT0 => ATMOS_REFSTATE_POTT
    use scale_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_XVW, &
       I_UYW, &
       I_UYZ, &
       I_XVZ, &
       I_XY, &
       I_UY, &
       I_XV
    use scale_comm, only: &
       COMM_vars8, &
       COMM_wait
    implicit none

    real(RP), intent(out) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP), intent(out) :: qflx_sgs_rhoq(KA,IA,JA,QA,3)

    real(RP), intent(inout) :: tke (KA,IA,JA) ! TKE
    real(RP), intent(out) :: Nu(KA,IA,JA) ! eddy viscosity (center)
    real(RP), intent(out) :: Pr(KA,IA,JA) ! Plandtle number
    real(RP), intent(out) :: Ri(KA,IA,JA) ! Richardson number

    real(RP), intent(in)  :: MOMZ(KA,IA,JA)
    real(RP), intent(in)  :: MOMX(KA,IA,JA)
    real(RP), intent(in)  :: MOMY(KA,IA,JA)
    real(RP), intent(in)  :: RHOT(KA,IA,JA)
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA)

    real(RP), intent(in)  :: SFLX_MW(IA,JA)
    real(RP), intent(in)  :: SFLX_MU(IA,JA)
    real(RP), intent(in)  :: SFLX_MV(IA,JA)
    real(RP), intent(in)  :: SFLX_SH(IA,JA)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor
    real(RP), intent(in)  :: dt

    real(RP) :: U(KA,IA,JA)    !< velocity in x-direction (full level)
    real(RP) :: V(KA,IA,JA)    !< velocity in y-direction (full level)
    real(RP) :: POTT(KA,IA,JA) !< potential temperature (full level)
    real(RP) :: phiN(KA,IA,JA)


    real(RP) :: n2(KA,IA,JA) !< square of the Brunt-Vaisala frequency, N^2
    real(RP) :: sm(KA,IA,JA) !< stability function for velocity
    real(RP) :: sh(KA,IA,JA) !< stability function for scalars
    real(RP) :: l(KA,IA,JA) !< length scale L
    real(RP) :: q(KA,IA,JA) !< q
    real(RP) :: dudz2(KA,IA,JA) !< (dudz)^2 + (dvdz)^2
    real(RP) :: q2_2(KA,IA,JA)  !< q^2 for level 2

    real(RP) :: a(KA,IA,JA)
    real(RP) :: b(KA,IA,JA)
    real(RP) :: c(KA,IA,JA)
    real(RP) :: d(KA)
    real(RP) :: ap
    real(RP) :: tke_N(KE_PBL,IA,JA)

    real(RP) :: l2q2         !< L^2/q^2
    real(RP) :: q2           !< q^2
    real(RP) :: ac           !< \alpha_c
    real(RP) :: ac2          !< \alpha_c^2
    real(RP) :: p1           !< \Phi_1
    real(RP) :: p2           !< \Phi_2
    real(RP) :: p3           !< \Phi_3
    real(RP) :: p4           !< \Phi_4
    real(RP) :: p5           !< \Phi_5
    real(RP) :: rd25         !< 1/D_2.5
    real(RP) :: gh           !< G_H
    
    real(RP) :: advc

    real(RP) :: sw

    integer :: k, i, j, iq
    integer :: IIS, IIE, JJS, JJE

!OCL XFILL
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE
       qflx_sgs_momz(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS  , JE
    do i = IS  , IE+1
    do k = KS-1, KE+1
       qflx_sgs_momx(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do

!OCL XFILL
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE+1
       qflx_sgs_momy(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE+1
       qflx_sgs_rhot(k,i,j,XDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do iq = 1, QA
    do j = JS  , JE
    do i = IS-1, IE
    do k = KS-1, KE+1
       qflx_sgs_rhoq(k,i,j,iq,XDIR) = 0.0_RP
    end do
    end do
    end do
    end do

!OCL XFILL
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE
       qflx_sgs_momz(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_momx(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS  , JE+1
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_momy(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_rhot(k,i,j,YDIR) = 0.0_RP
    end do
    end do
    end do
!OCL XFILL
    do iq = 1, QA
    do j = JS-1, JE
    do i = IS  , IE
    do k = KS-1, KE+1
       qflx_sgs_rhoq(k,i,j,iq,YDIR) = 0.0_RP
    end do
    end do
    end do
    end do
!OCL XFILL
    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE+1
       qflx_sgs_momz(k,i,j,ZDIR) = 0.0_RP
    end do
    end do
    end do



    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do j = JJS  , JJE+1
       do i = IIS-1, IIE+1
       do k = KS, KE_PBL+1
          U(k,i,j) = 0.5_RP * ( MOMX(k,i,j) + MOMX(k,i-1,j) ) / DENS(k,i,j)
       end do
       end do
       end do
       do j = JJS  , JJE+1
       do i = IIS-1, IIE+1
          U(KS-1,i,j) = 0.0_RP
       end do
       end do

       do j = JJS-1, JJE+1
       do i = IIS  , IIE+1
       do k = KS, KE_PBL+1
          V(k,i,j) = 0.5_RP * ( MOMY(k,i,j) + MOMY(k,i,j-1) ) / DENS(k,i,j)
       end do
       end do
       end do
       do j = JJS-1, JJE+1
       do i = IIS  , IIE+1
          V(KS-1,i,j) = 0.0_RP
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL+1
          POTT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
       end do
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS+1, KE_PBL
          n2(k,i,j) = GRAV * ( POTT(k+1,i,j) - POTT(k-1,i,j)) * J33G &
               / ( (FDZ(k)+FDZ(k-1)) * GSQRT(k,i,j,I_XYZ) * PT0(k,i,j) )
       end do
       end do
       end do
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          n2(KS,i,j) = GRAV * ( POTT(KS+1,i,j) - POTT(KS,i,j)) * J33G &
               * RFDZ(KS) / ( GSQRT(KS,i,j,I_XYZ) * PT0(KS,i,j) )
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL
          dudz2(k,i,j) = ( ( U(k+1,i,j) - U(k-1,i,j) )**2 + ( V(k+1,i,j) - V(k-1,i,j) )**2 ) &
               / ( FDZ(k-1)*GSQRT(k-1,i,j,I_XYZ) + FDZ(k)*GSQRT(k,i,j,I_XYZ) )**2
          sw = sign(0.5_RP, dudz2(k,i,j)-EPS) + 0.5_RP
          dudz2(k,i,j) = dudz2(k,i,j)*sw + 1.E-10_RP*(1.0_RP-sw)
       end do
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL
          Ri(k,i,j) = n2(k,i,j)/dudz2(k,i,j)
       end do
       end do
       end do


       ! length
       call get_length( &
            l, & ! (out)
            DENS, & ! (in)
            tke, n2, & ! (in)
            SFLX_MU, SFLX_MV, SFLX_SH, & ! (in)
            PT0, & ! (in)
            GSQRT(:,:,:,I_XYZ), & ! (in)
            IIS, IIE, JJS, JJE )

       call get_q2_level2( &
            q2_2, & ! (out)
            dudz2, Ri, l, & ! (in)
            IIS, IIE, JJS, JJE ) ! (in)

       if ( ATMOS_PHY_TB_MYNN_TKE_INIT ) then
          do j = JJS, JJE+1
          do i = IIS, IIE+1
          do k = KS, KE_PBL
             tke(k,i,j) = q2_2(k,i,j) * 0.5_RP
          end do
          end do
          end do
          do j = JJS, JJE+1
          do i = IIS, IIE+1
          do k = KE_PBL+1, KE
             tke(k,i,j) = 0.0_RP
          end do
          end do
          end do
          call COMM_vars8(tke, 1)
          call COMM_wait (tke, 1)
       end if

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL
          sw = sign(0.5_RP, tke(k,i,j)-EPS) + 0.5_RP
          q(k,i,j) = sqrt( tke(k,i,j)*2.0_RP*sw + 1.E-10_RP*(1.0_RP-sw) )
       end do
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL

          ! level 2.5
          ac = min(q(k,i,j)/sqrt(q2_2(k,i,j)), 1.0_RP)
          ac2 = ac**2
          l2q2 = ( l(k,i,j) / q(k,i,j) )**2
          gh = - n2(k,i,j) * l2q2

          p1 = 1.0_RP - 3.0_RP * ac2 * A2 * B2 * (1.0_RP-C3) * gh
          p2 = 1.0_RP - 9.0_RP * ac2 * A1 * A2 * (1.0_RP-C2) * gh
          p3 = p1 + 9.0_RP * ac2 * A2**2 * (1.0_RP-C2) * (1.0_RP-C5) * gh
          p4 = p1 - 12.0_RP * ac2 * A1 * A2 * (1.0_RP-C2) * gh
          p5 = 6.0_RP * ac2 * A1**2 * dudz2(k,i,j) * l2q2

          rd25 = 1.0_RP / max(p2 * p4 + p5 * p3, 1.E-20_RP)
          sm(k,i,j) = ac * A1 * (p3 - 3.0_RP * C1 * p4) * rd25
          sh(k,i,j) = ac * A2 * (p2 + 3.0_RP * C1 * p5) * rd25
       end do
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL
          Nu(k,i,j) = l(k,i,j) * q(k,i,j) * sm(k,i,j)
       end do
       end do
       end do
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KE_PBL+1, KE
          Nu(k,i,j) = 0.0_RP
       end do
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KS, KE_PBL
          sw = 0.5_RP - sign(0.5_RP, sh(k,i,j)-EPS)
          Pr(k,i,j) = sm(k,i,j) / (sh(k,i,j)+sw) * (1.0_RP-sw) &
                    + 1.0_RP * sw
       end do
       end do
       end do
       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KE_PBL+1, KE
          Pr(k,i,j) = 1.0_RP
       end do
       end do
       end do

       do j = JJS, JJE+1
       do i = IIS, IIE+1
       do k = KE_PBL+1, KE
          Ri(k,i,j) = 0.0_RP
       end do
       end do
       end do



       ! time integration

       !  for velocities
       do j = JJS, JJE+1
       do i = IIS, IIE+1
          ap = - dt * 0.5_RP * ( DENS(KS  ,i,j)*Nu(KS  ,i,j) &
                               + DENS(KS+1,i,j)*Nu(KS+1,i,j) ) &
             * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
          a(KS,i,j) = ap * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) ) 
          c(KS,i,j) = 0.0_RP
          b(KS,i,j) = - a(KS,i,j) + 1.0_RP
          do k = KS+1, KE_PBL-1
             c(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) ) 
             ap = - dt * 0.5_RP * ( DENS(k  ,i,j)*Nu(k  ,i,j) &
                                  + DENS(k+1,i,j)*Nu(k+1,i,j) ) &
                * RFDZ(k) / GSQRT(k,i,j,I_XYW)
             a(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) ) 
             b(k,i,j) = - a(k,i,j) - c(k,i,j) + 1.0_RP
          end do
          a(KE_PBL,i,j) = 0.0_RP
          c(KE_PBL,i,j) = ap * RCDZ(KE_PBL) / ( DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) ) 
          b(KE_PBL,i,j) = - c(KE_PBL,i,j) + 1.0_RP
       end do
       end do


       ! integration U
       do j = JJS, JJE
       do i = IIS, IIE+1
          do k = KS, KE_PBL
             d(k) = U(k,i,j)
          end do
          call diffusion_solver( &
               phiN(:,i,j),                    & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d ) ! (in)
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE_PBL-1
          qflx_sgs_momx(k,i,j,ZDIR) = - 0.03125_RP & ! 1/4/4/2
               * ( DENS(k,i,j) + DENS(k+1,i,j) + DENS(k,i+1,j) + DENS(k+1,i+1,j) ) &
               * ( Nu(k,i,j) + Nu(k+1,i,j) + Nu(k,i+1,j) + Nu(k+1,i+1,j) ) &
               * ( (phiN(k+1,i,j)+phiN(k+1,i+1,j)) - (phiN(k,i,j)+phiN(k,i+1,j)) ) &
               * J33G * RFDZ(k) / GSQRT(k,i,j,I_UYW)
       end do
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momx(KS-1,i,j,ZDIR) = 0.0_RP
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KE_PBL, KE
          qflx_sgs_momx(k,i,j,ZDIR) = 0.0_RP
       end do
       end do
       end do


       ! integration V
       do j = JJS, JJE+1
       do i = IIS, IIE
          do k = KS, KE_PBL
             d(k) = V(k,i,j)
          end do
          call diffusion_solver( &
               phiN(:,i,j),                    & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d ) ! (in)
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE_PBL-1
          qflx_sgs_momy(k,i,j,ZDIR) = - 0.03125_RP & ! 1/4/4/2
               * ( DENS(k,i,j) + DENS(k+1,i,j) + DENS(k,i,j+1) + DENS(k+1,i,j+1) ) &
               * ( Nu(k,i,j) + Nu(k+1,i,j) + Nu(k,i,j+1) + Nu(k+1,i,j+1) ) &
               * ( (phiN(k+1,i,j)+phiN(k+1,i,j+1)) - (phiN(k,i,j)+phiN(k,i,j+1)) ) &
               * J33G * RFDZ(k) / GSQRT(k,i,j,I_XVW)
       end do
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_momy(KS-1,i,j,ZDIR) = 0.0_RP
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KE_PBL, KE
          qflx_sgs_momy(k,i,j,ZDIR) = 0.0_RP
       end do
       end do
       end do


       !  for scalars
       do j = JJS, JJE
       do i = IIS, IIE
          ap = - dt * 0.5_RP * ( DENS(KS  ,i,j)*Nu(KS  ,i,j)/Pr(KS  ,i,j) &
                               + DENS(KS+1,i,j)*Nu(KS+1,i,j)/Pr(KS+1,i,j) ) &
             * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
          a(KS,i,j) = ap * RCDZ(KS) / (DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          c(KS,i,j) = 0.0_RP
          b(KS,i,j) = - a(KS,i,j) + 1.0_RP
          do k = KS+1, KE_PBL-1
             c(k,i,j) = ap * RCDZ(k) / (DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             ap = - dt * 0.5_RP * ( DENS(k  ,i,j)*Nu(k  ,i,j)/Pr(k  ,i,j) &
                                  + DENS(k+1,i,j)*Nu(k+1,i,j)/Pr(k+1,i,j) ) &
                * RFDZ(k) / GSQRT(k,i,j,I_XYW)
             a(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             b(k,i,j) = - a(k,i,j) - c(k,i,j) + 1.0_RP
          end do
          a(KE_PBL,i,j) = 0.0_RP
          c(KE_PBL,i,j) = ap * RCDZ(KE_PBL) / (DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) )
          b(KE_PBL,i,j) = - c(KE_PBL,i,j) + 1.0_RP
       end do
       end do

       ! integration POTT
       do j = JJS, JJE
       do i = IIS, IIE
          do k = KS, KE_PBL
             d(k) = POTT(k,i,j)
          end do
          call diffusion_solver( &
               phiN(:,i,j),                    & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d ) ! (in)
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KS, KE_PBL-1
          qflx_sgs_rhot(k,i,j,ZDIR) = - 0.25_RP & ! 1/2/2
               * ( DENS(k,i,j) + DENS(k+1,i,j) ) &
               * ( Nu(k,i,j)/Pr(k,i,j) + Nu(k+1,i,j)/Pr(k+1,i,j) ) &
               * J33G * ( phiN(k+1,i,j) - PhiN(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XYW)
       end do
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
          qflx_sgs_rhot(KS-1,i,j,ZDIR) = 0.0_RP
       end do
       end do
       do j = JJS, JJE
       do i = IIS, IIE
       do k = KE_PBL, KE
          qflx_sgs_rhot(k,i,j,ZDIR) = 0.0_RP
       end do
       end do
       end do


       ! integration QTRC
       do iq = 1, QA
          do j = JJS, JJE
          do i = IIS, IIE
             do k = KS, KE_PBL
                d(k) = QTRC(k,i,j,iq)
             end do
             call diffusion_solver( &
                  phiN(:,i,j),                    & ! (out)
                  a(:,i,j), b(:,i,j), c(:,i,j), d ) ! (in)
          end do
          end do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE_PBL-1
             qflx_sgs_rhoq(k,i,j,iq,ZDIR) = - 0.25_RP & ! 1/2/2
                  * ( Nu(k,i,j)/Pr(k,i,j) + Nu(k+1,i,j)/Pr(k+1,i,j) ) &
                  * ( DENS(k,i,j) + DENS(k+1,i,j) ) &
                  * J33G * ( phiN(k+1,i,j) - phiN(k,i,j) ) * RFDZ(k) / GSQRT(k,i,j,I_XYW)
          end do
          end do
          end do
          do j = JJS, JJE
          do i = IIS, IIE
             qflx_sgs_rhoq(KS-1,i,j,iq,ZDIR) = 0.0_RP
          end do
          end do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KE_PBL, KE
             qflx_sgs_rhoq(k,i,j,iq,ZDIR) = 0.0_RP
          end do
          end do
          end do
       end do

       ! time integration tke
       do j = JJS, JJE
       do i = IIS, IIE
          ap = - dt * 1.5_RP * ( DENS(KS  ,i,j)*Nu(KS  ,i,j) &
                               + DENS(KS+1,i,j)*Nu(KS+1,i,j) ) &
             * RFDZ(KS) / GSQRT(KS,i,j,I_XYW)
          a(KS,i,j) = ap * RCDZ(KS) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          c(KS,i,j) = 0.0_RP
          b(KS,i,j) = - a(KS,i,j) + 1.0_RP + 2.0_RP * dt * q(KS,i,j) / ( B1 * l(KS,i,j) )
          advc =  0.5_RP * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                * ( ( GSQRT(KS,i  ,j,I_UYZ) * MOMX(KS,i  ,j) * (tke(KS,i+1,j)+tke(KS,i  ,j)) / MAPF(i  ,j,2,I_UY) &
                    - GSQRT(KS,i-1,j,I_UYZ) * MOMX(KS,i-1,j) * (tke(KS,i  ,j)+tke(KS,i-1,j)) / MAPF(i-1,j,2,I_UY) ) * RCDX(i) &
                  + ( GSQRT(KS,i,j  ,I_XVZ) * MOMY(KS,i,j  ) * (tke(KS,i,j+1)+tke(i,i,j  )) / MAPF(i,j  ,1,I_XV) &
                    - GSQRT(KS,i,j-1,I_XVZ) * MOMY(KS,i,j-1) * (tke(KS,i,j  )+tke(i,j,j-1)) / MAPF(i,j-1,1,I_XV) ) * RCDY(j) &
                  + ( ( J13G(KS+1,i,j,I_XYZ) * (MOMX(KS+1,i,j)+MOMX(KS+1,i-1,j)) * tke(KS+1,i,j) &
                      - J13G(KS  ,i,j,I_XYZ) * (MOMX(KS  ,i,j)+MOMX(KS  ,i-1,j)) * tke(KS  ,i,j) ) / MAPF(i,j,2,I_XY) &
                    + ( J23G(KS+1,i,j,I_XYZ) * (MOMY(KS+1,i,j)+MOMY(KS+1,i,j-1)) * tke(KS+1,i,j) &
                      - J23G(KS  ,i,j,I_XYZ) * (MOMY(KS  ,i,j)+MOMY(KS  ,i,j-1)) * tke(KS  ,i,j) ) / MAPF(i,j,1,I_XY) &
                    ) * RFDZ(KS) &
                    + ( J33G * MOMZ(KS  ,i,j) * (tke(KS+1,i,j)+tke(KS  ,i,j)) ) &
                      * RCDZ(KS) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) ) &
                  ) / ( DENS(KS,i,j) * GSQRT(KS,i,j,I_XYZ) )
          d(KS) = tke(KS,i,j) + dt * ( Nu(KS,i,j) * (dudz2(KS,i,j) - n2(KS,i,j)/Pr(KS,i,j)) &
                                     - advc )
          do k = KS+1, KE_PBL-1
             c(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             ap = - dt * 1.5_RP * ( DENS(k  ,i,j)*Nu(k  ,i,j) &
                                  + DENS(k+1,i,j)*Nu(k+1,i,j) ) &
                * RFDZ(k) / GSQRT(k,i,j,I_XYW)
             a(k,i,j) = ap * RCDZ(k) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             b(k,i,j) = - a(k,i,j) - c(k,i,j) + 1.0_RP + 2.0_RP * dt * q(k,i,j) / ( B1 * l(k,i,j) )
             advc =  0.5_RP * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
                * ( ( GSQRT(k,i  ,j,I_UYZ) * MOMX(k,i  ,j) * (tke(k,i+1,j)+tke(k,i  ,j)) / MAPF(i  ,j,2,I_UY) &
                    - GSQRT(k,i-1,j,I_UYZ) * MOMX(k,i-1,j) * (tke(k,i  ,j)+tke(k,i-1,j)) / MAPF(i-1,j,2,I_UY) ) * RCDX(i) &
                  + ( GSQRT(k,i,j  ,I_XVZ) * MOMY(k,i,j  ) * (tke(k,i,j+1)+tke(i,i,j  )) / MAPF(i,j  ,1,I_XV) &
                    - GSQRT(k,i,j-1,I_XVZ) * MOMY(k,i,j-1) * (tke(k,i,j  )+tke(i,j,j-1)) / MAPF(i,j-1,1,I_XV) ) * RCDY(j) &
                  + ( ( J13G(k+1,i,j,I_XYZ) * (MOMX(k+1,i,j)+MOMX(k+1,i-1,j)) * tke(k+1,i,j) &
                      - J13G(k-1,i,j,I_XYZ) * (MOMX(k-1,i,j)+MOMX(k-1,i-1,j)) * tke(k-1,i,j) ) / MAPF(i,j,2,I_XY) &
                    + ( J23G(k+1,i,j,I_XYZ) * (MOMY(k+1,i,j)+MOMY(k+1,i,j-1)) * tke(k+1,i,j) &
                      - J23G(k-1,i,j,I_XYZ) * (MOMY(k-1,i,j)+MOMY(k-1,i,j-1)) * tke(k-1,i,j) ) / MAPF(i,j,1,I_XY) &
                    ) / (FDZ(k)+FDZ(k-1)) &
                    + ( J33G * MOMZ(k  ,i,j) * (tke(k+1,i,j)+tke(k  ,i,j)) &
                      - J33G * MOMZ(k-1,i,j) * (tke(k  ,i,j)+tke(k-1,i,j)) ) &
                      * RCDZ(k) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) ) &
                  ) / ( DENS(k,i,j) * GSQRT(k,i,j,I_XYZ) )
             d(k) = tke(k,i,j) &
                  + dt * ( Nu(k,i,j) * (dudz2(k,i,j) - n2(k,i,j)/Pr(k,i,j)) &
                         - advc )
          end do
          a(KE_PBL,i,j) = 0.0_RP
          c(KE_PBL,i,j) = ap * RCDZ(KE_PBL) / ( DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) )
          b(KE_PBL,i,j) = - c(KE_PBL,i,j) + 1.0_RP + 2.0_RP * dt * q(KE_PBL,i,j) / ( B1 * l(KE_PBL,i,j) )
          advc =  0.5_RP * MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) &
               * ( ( GSQRT(KE_PBL,i  ,j,I_UYZ) * MOMX(KE_PBL,i  ,j) * (tke(KE_PBL,i+1,j)+tke(KE_PBL,i  ,j)) / MAPF(i  ,j,2,I_UY) &
                   - GSQRT(KE_PBL,i-1,j,I_UYZ) * MOMX(KE_PBL,i-1,j) * (tke(KE_PBL,i  ,j)+tke(KE_PBL,i-1,j)) / MAPF(i-1,j,2,I_UY) ) &
                   * RCDX(i) &
                 + ( GSQRT(KE_PBL,i,j  ,I_XVZ) * MOMY(KE_PBL,i,j  ) * (tke(KE_PBL,i,j+1)+tke(i,i,j  )) / MAPF(i,j  ,1,I_XV) &
                   - GSQRT(KE_PBL,i,j-1,I_XVZ) * MOMY(KE_PBL,i,j-1) * (tke(KE_PBL,i,j  )+tke(i,j,j-1)) / MAPF(i,j-1,1,I_XV) ) &
                   * RCDY(j) &
                 + ( ( J13G(KE_PBL  ,i,j,I_XYZ) * (MOMX(KE_PBL  ,i,j)+MOMX(KE_PBL  ,i-1,j)) * tke(KE_PBL  ,i,j) &
                     - J13G(KE_PBL-1,i,j,I_XYZ) * (MOMX(KE_PBL-1,i,j)+MOMX(KE_PBL-1,i-1,j)) * tke(KE_PBL-1,i,j) ) &
                     / MAPF(i,j,2,I_XY) &
                   + ( J23G(KE_PBL  ,i,j,I_XYZ) * (MOMY(KE_PBL  ,i,j)+MOMY(KE_PBL  ,i,j-1)) * tke(KE_PBL  ,i,j) &
                     - J23G(KE_PBL-1,i,j,I_XYZ) * (MOMY(KE_PBL-1,i,j)+MOMY(KE_PBL-1,i,j-1)) * tke(KE_PBL-1,i,j) ) &
                     / MAPF(i,j,1,I_XY) &
                   ) * RFDZ(KE_PBL) &
                   + ( - J33G * MOMZ(KE_PBL-1,i,j) * (tke(KE_PBL  ,i,j)+tke(KE_PBL-1,i,j)) ) &
                     * RCDZ(KE_PBL) / ( MAPF(i,j,1,I_XY) * MAPF(i,j,2,I_XY) ) &
                 ) / ( DENS(KE_PBL,i,j) * GSQRT(KE_PBL,i,j,I_XYZ) )
          d(KE_PBL) = tke(KE_PBL,i,j) + dt * ( Nu(KE_PBL,i,j) * (dudz2(KE_PBL,i,j) - n2(KE_PBL,i,j)/Pr(KE_PBL,i,j)) &
                                             - advc )
          call diffusion_solver( &
               tke_N(:,i,j),     & ! (out)
               a(:,i,j), b(:,i,j), c(:,i,j), d ) ! (in)
       end do
       end do
    end do
    end do

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE_PBL
          tke(k,i,j) = max(tke_N(k,i,j), 1.0E-10_RP)
       end do
       do k = KE_PBL+1, KE
          tke(k,i,j) = 0.0_RP
       end do
    end do
    end do

    ATMOS_PHY_TB_MYNN_TKE_INIT = .false.

    return
  end subroutine ATMOS_PHY_TB_mynn

  subroutine diffusion_solver( &
       phi, &
       a, b, c, d )
    implicit none
    real(RP), intent(out) :: phi(KA)
    real(RP), intent(in)  :: a(KA)
    real(RP), intent(in)  :: b(KA)
    real(RP), intent(in)  :: c(KA)
    real(RP), intent(in)  :: d(KA)
    real(RP) :: e(KA)
    real(RP) :: f(KA)
    real(RP) :: denom
    integer :: k

    e(KS) = - a(KS) / b(KS)
    f(KS) =   d(KS) / b(KS)
    do k = KS+1, KE_PBL-1
       denom = b(k) + c(k)*e(k-1)
       e(k) = - a(k) / denom
       f(k) = ( d(k) - c(k)*f(k-1) ) / denom
    end do

    ! flux at the top boundary is zero
    phi(KE_PBL) = ( d(KE_PBL) - c(KE_PBL)*f(KE_PBL-1) ) / ( b(KE_PBL) + c(KE_PBL)*e(KE_PBL-1) ) ! = f(KE_PBL)

    do k = KE_PBL-1, KS, -1
       phi(k) = e(k) * phi(k+1) + f(k)
    end do

    return
  end subroutine diffusion_solver

  subroutine get_length( &
       l, &
       DENS, &
       tke, n2, &
       SFLX_MU, SFLX_MV, SFLX_SH, &
       PT0, &
       GSQRT, &
       IIS, IIE, JJS, JJE )
    use scale_grid, only: &
       CDZ  => GRID_CDZ
    use scale_grid_real, only: &
       FZ => REAL_FZ
    use scale_const, only: &
       GRAV   => CONST_GRAV, &
       KARMAN => CONST_KARMAN, &
       CP     => CONST_CPdry, &
       EPS    => CONST_EPS
    implicit none
    real(RP), intent(out) :: l(KA,IA,JA)
    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: tke(KA,IA,JA)
    real(RP), intent(in) :: n2(KA,IA,JA)
    real(RP), intent(in) :: SFLX_MU(IA,JA)
    real(RP), intent(in) :: SFLX_MV(IA,JA)
    real(RP), intent(in) :: SFLX_SH(IA,JA)
    real(RP), intent(in) :: GSQRT(KA,IA,JA)
    real(RP), intent(in) :: PT0(KA,IA,JA)
    integer,  intent(in) :: IIS
    integer,  intent(in) :: IIE
    integer,  intent(in) :: JJS
    integer,  intent(in) :: JJE

    real(RP) :: ls           !< L_S
    real(RP) :: lt           !< L_T
    real(RP) :: lb           !< L_B
    real(RP) :: rlm          !< 1/L_M
    real(RP) :: rlt          !< 1/L_T

    real(RP) :: q            !< q
    real(RP) :: qc           !< q_c
    real(RP) :: int_q        !< \int q dz
    real(RP) :: int_qz       !< \int qz dz
    real(RP) :: rn2sr         !< 1/N
    real(RP) :: us           !< friction velocity
    real(RP) :: wtg          !< heat flux at the surface
    real(RP) :: zeta         !< height normalized by the Obukhov length

    real(RP) :: z
    real(RP) :: qdz

    real(RP) :: sw, sw2
    integer :: k, i, j

    do j = JJS, JJE+1
    do i = IIS, IIE+1
       int_qz = 0.0_RP
       int_q = 0.0_RP
       do k = KS, KE_PBL
          q = sqrt(tke(k,i,j) * 2.0_RP)
          qdz =  q * CDZ(k) * GSQRT(k,i,j)
          int_qz = int_qz + ((FZ(k,i,j)+FZ(k-1,i,j))*0.5_RP-FZ(KS-1,i,j)) * qdz
          int_q  = int_q + qdz
       end do
       ! LT
       lt = max(0.23_RP * int_qz / (int_q + EPS), LT_min)
       rlt = 1.0_RP / lt

       us = ( (SFLX_MU(i,j)**2 + SFLX_MV(i,j)**2)*0.5_RP )**0.25_RP / DENS(KS,i,j) ! friction velocity
       us = max(us, Us_min)
       wtg = max(SFLX_SH(i,j) / (CP * DENS(KS,i,j)), 0.0_RP) ! surface heat flux
       rlm = - KARMAN * GRAV * wtg / (PT0(KS,i,j) * us**3 )

       do k = KS, KE_PBL
          z = ( FZ(k,i,j)+FZ(k-1,i,j) )*0.5_RP - FZ(KS-1,i,j)
          zeta = z * rlm

          ! LS
          sw = sign(0.5_RP, zeta) + 0.5_RP ! 1 for zeta >= 0, 0 for zeta < 0
          ls = KARMAN * z &
             * ( 1.0_RP / (1.0_RP + 2.7_RP*min(zeta,1.0_RP)*sw ) * sw &
               + ((1.0_RP - 100.0_RP*zeta)*(1.0_RP-sw))**0.2_RP )

          ! LB
          sw = sign(0.5_RP, tke(k,i,j)-EPS) + 0.5_RP
          q = sqrt(tke(k,i,j) * 2.0_RP)*sw + 1.E-10_RP*(1.0_RP-sw)
          qc = (GRAV/PT0(k,i,j)*wtg*lt)**OneOverThree
          sw  = sign(0.5_RP, n2(k,i,j)-EPS) + 0.5_RP ! 1 for dptdz >0, 0 for dptdz < 0
          rn2sr = 1.0_RP / ( sqrt(n2(k,i,j)*sw) + 1.0_RP-sw)
          lb = (1.0_RP + 5.0_RP * sqrt(qc*rn2sr/lt)) * q * rn2sr * sw & ! qc=0 when wtg < 0
             +  999.E10_RP * (1.0_RP-sw)

          ! L
          l(k,i,j) = 1.0_RP / ( 1.0_RP/ls + rlt + 1.0_RP/lb )
       end do
    end do
    end do

    return
  end subroutine get_length

  subroutine get_q2_level2( &
       q2_2, &
       dudz2, Ri, l, &
       IIS, IIE, JJS, JJE)
    use scale_const, only: &
       EPS    => CONST_EPS
    implicit none
    real(RP), intent(out) :: q2_2(KA,IA,JA)
    real(RP), intent(in)  :: dudz2(KA,IA,JA)
    real(RP), intent(in)  :: Ri(KA,IA,JA)
    real(RP), intent(in)  :: l(KA,IA,JA)
    integer,  intent(in)  :: IIS
    integer,  intent(in)  :: IIE
    integer,  intent(in)  :: JJS
    integer,  intent(in)  :: JJE

    real(RP) :: rf           !< Rf
    real(RP) :: sm_2         !< sm for level 2
    real(RP) :: sh_2         !< sh for level 2

    real(RP) :: q2
    real(RP) :: sw
    integer :: k, i, j

    do j = JJS, JJE+1
    do i = IIS, IIE+1
    do k = KS, KE_PBL
       rf = min(0.5_RP / AF12 * ( Ri(k,i,j) &
                              + AF12*Rf1 &
                              - sqrt(Ri(k,i,j)**2 + 2.0_RP*AF12*(Rf1-2.0_RP*Rf2)*Ri(k,i,j) + (AF12*Rf1)**2) ), &
                Rfc)
       sh_2 = 3.0_RP * A2 * (G1+G2) * (Rfc-rf) / (1.0_RP-rf)
       sm_2 = sh_2 * AF12 * (Rf1-rf) / (Rf2-rf)
       q2 = B1 * l(k,i,j)**2 * sm_2 * (1.0_RP-rf) * dudz2(k,i,j)
       sw = sign(0.5_RP,q2-EPS) + 0.5_RP
       q2_2(k,i,j) = q2*sw + 1.E-10_RP*(1.0_RP-sw)
    end do
    end do
    end do

    return
  end subroutine get_q2_level2

end module scale_atmos_phy_tb_mynn
