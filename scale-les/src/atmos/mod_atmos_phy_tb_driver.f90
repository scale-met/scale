!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Turbulence
!!
!! @par Description
!!          Sub-grid scale turbulence process
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-12-05 (S.Nishizawa)       [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_phy_tb_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_grid_index
  use mod_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_TB_driver_setup
  public :: ATMOS_PHY_TB_driver

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
  real(RP), private, allocatable :: MOMZ_t(:,:,:)
  real(RP), private, allocatable :: MOMX_t(:,:,:)
  real(RP), private, allocatable :: MOMY_t(:,:,:)
  real(RP), private, allocatable :: RHOT_t(:,:,:)
  real(RP), private, allocatable :: QTRC_t(:,:,:,:)
  !-----------------------------------------------------------------------------
contains

  subroutine ATMOS_PHY_TB_driver_setup( TB_TYPE )
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only: &
       CDZ => GRID_CDZ, &
       CDX => GRID_CDX, &
       CDY => GRID_CDY, &
       CZ  => GRID_CZ
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB_setup
    implicit none

    character(len=H_SHORT), intent(in) :: TB_TYPE
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Physics-TB]/Categ[ATMOS]'

    allocate( MOMZ_t(KA,IA,JA) )
    allocate( MOMX_t(KA,IA,JA) )
    allocate( MOMY_t(KA,IA,JA) )
    allocate( RHOT_t(KA,IA,JA) )
    allocate( QTRC_t(KA,IA,JA,QA) )

    call ATMOS_PHY_TB_setup( &
         TB_TYPE,       & ! (in)
         CDZ, CDX, CDY, & ! (in)
         CZ             ) ! (in)

    call ATMOS_PHY_TB_driver( .true., .false. )

  end subroutine ATMOS_PHY_TB_driver_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_TB_driver( update_flag, history_flag )
    use mod_time, only: &
       dttb => TIME_DTSEC_ATMOS_PHY_TB
    use mod_history, only: &
       HIST_in
    use mod_grid, only: &
       RCDZ => GRID_RCDZ, &
       RCDX => GRID_RCDX, &
       RCDY => GRID_RCDY, &
       RFDZ => GRID_RFDZ, &
       RFDX => GRID_RFDX, &
       RFDY => GRID_RFDY, &
       CDZ  => GRID_CDZ,  &
       FDZ  => GRID_FDZ
    use mod_gridtrans, only: &
       I_XYZ, &
       I_XYW, &
       I_UYW, &
       I_XVW, &
       I_UYZ, &
       I_XVZ, &
       I_UVZ, &
       GSQRT => GTRANS_GSQRT, &
       J13G  => GTRANS_J13G,  &
       J23G  => GTRANS_J23G,  &
       J33G  => GTRANS_J33G
    use mod_atmos_phy_tb, only: &
       ATMOS_PHY_TB
    use mod_atmos_vars, only: &
       DENS_av, &
       MOMZ_av, &
       MOMX_av, &
       MOMY_av, &
       RHOT_av, &
       QTRC_av, &
       MOMZ_tp, &
       MOMX_tp, &
       MOMY_tp, &
       RHOT_tp, &
       QTRC_tp
    implicit none

    logical, intent(in) :: update_flag
    logical, intent(in), optional :: history_flag

    ! eddy viscosity/diffusion flux
    real(RP) :: qflx_sgs_momz(KA,IA,JA,3)
    real(RP) :: qflx_sgs_momx(KA,IA,JA,3)
    real(RP) :: qflx_sgs_momy(KA,IA,JA,3)
    real(RP) :: qflx_sgs_rhot(KA,IA,JA,3)
    real(RP) :: qflx_sgs_qtrc(KA,IA,JA,QA,3)

    ! diagnostic variables
    real(RP) :: tke(KA,IA,JA) ! TKE
    real(RP) :: nu (KA,IA,JA) ! eddy diffusion
    real(RP) :: Ri (KA,IA,JA) ! Richardoson number
    real(RP) :: Pr (KA,IA,JA) ! Prandtle number

    integer :: k, i, j, iq
    integer :: IIS, IIE, JJS, JJE

    if ( update_flag ) then
       call ATMOS_PHY_TB( &
            qflx_sgs_momz, qflx_sgs_momx, qflx_sgs_momy, & ! (out)
            qflx_sgs_rhot, qflx_sgs_qtrc,                & ! (out)
            tke, nu, Ri, Pr,                             & ! (out) diagnostic variables
            MOMZ_av, MOMX_av, MOMY_av, RHOT_av, DENS_av, QTRC_av, & ! (in)
            GSQRT, J13G, J23G, J33G                      & ! (in)
            )


       do JJS = JS, JE, JBLOCK
       JJE = JJS+JBLOCK-1
       do IIS = IS, IE, IBLOCK
       IIE = IIS+IBLOCK-1
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMZ_t(k,i,j) = - ( &
                  + ( GSQRT(k,i  ,j,I_UYW)*qflx_sgs_momz(k,i  ,j,XDIR) &
                    - GSQRT(k,i-1,j,I_UYW)*qflx_sgs_momz(k,i-1,j,XDIR) ) * RCDX(i) &
                  + ( GSQRT(k,i,j  ,I_XVW)*qflx_sgs_momz(k,i,j  ,YDIR) &
                    - GSQRT(k,i,j-1,I_XVW)*qflx_sgs_momz(k,i,j-1,YDIR) ) * RCDY(j) &
                  + ( J13G(k+1,i,j,I_XYZ)*(qflx_sgs_momz(k+1,i,j,XDIR)+qflx_sgs_momz(k+1,i-1,j,XDIR)) &
                    - J13G(k-1,i,j,I_XYZ)*(qflx_sgs_momz(k-1,i,j,XDIR)-qflx_sgs_momz(k-1,i,j-1,XDIR)) &
                    + J23G(k+1,i,j,I_XYZ)*(qflx_sgs_momz(k+1,i,j,YDIR)+qflx_sgs_momz(k+1,i,j-1,YDIR)) &
                    - J23G(k-1,i,j,I_XYZ)*(qflx_sgs_momz(k-1,i,j,YDIR)+qflx_sgs_momz(k-1,i,j-1,YDIR)) &
                    ) * 0.5_RP / ( CDZ(k+1)+CDZ(k) ) &
                  + J33G * ( qflx_sgs_momz(k+1,i,j,ZDIR) - qflx_sgs_momz(k,i,j,ZDIR) ) * RFDZ(k) &
                  ) / GSQRT(k,i,j,I_XYW)
          end do
          end do
          end do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMX_t(k,i,j) = - ( &
                  + ( GSQRT(k,i+1,j,I_XYZ)*qflx_sgs_momx(k,i+1,j,XDIR) &
                    - GSQRT(k,i  ,j,I_XYZ)*qflx_sgs_momx(k,i  ,j,XDIR) ) * RFDX(i) &
                  + ( GSQRT(k,i,j  ,I_UVZ)*qflx_sgs_momx(k,i,j  ,YDIR) &
                    - GSQRT(k,i,j-1,I_UVZ)*qflx_sgs_momx(k,i,j-1,YDIR) ) * RCDY(j) &
                  + ( J13G(k+1,i,j,I_UYW)*(qflx_sgs_momx(k+1,i+1,j,XDIR)+qflx_sgs_momx(k+1,i,j,XDIR)) &
                    - J13G(k-1,i,j,I_UYW)*(qflx_sgs_momx(k-1,i+1,j,XDIR)+qflx_sgs_momx(k-1,i,j,XDIR)) &
                    + J23G(k+1,i,j,I_UYW)*(qflx_sgs_momx(k+1,i,j,YDIR)+qflx_sgs_momx(k+1,i,j-1,YDIR)) &
                    - J23G(k-1,i,j,I_UYW)*(qflx_sgs_momx(k-1,i,j,YDIR)+qflx_sgs_momx(k-1,i,j-1,YDIR)) &
                    ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
                  + J33G * ( qflx_sgs_momx(k,i,j,ZDIR) - qflx_sgs_momx(k-1,i,j,ZDIR) ) * RCDZ(k) &
                  ) / GSQRT(k,i,j,I_UYZ)
          end do
          end do
          end do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             MOMY_t(k,i,j) = - ( &
                  + ( GSQRT(k,i  ,j,I_UVZ)*qflx_sgs_momy(k,i  ,j,XDIR) &
                    - GSQRT(k,i-1,j,I_UVZ)*qflx_sgs_momy(k,i-1,j,XDIR) ) * RCDX(i) &
                  + ( GSQRT(k,i,j+1,I_XYZ)*qflx_sgs_momy(k,i,j+1,YDIR) &
                    - GSQRT(k,i,j  ,I_XYZ)*qflx_sgs_momy(k,i,j  ,YDIR) ) * RFDY(j) &
                  + ( J13G(k+1,i,j,I_XVW)*(qflx_sgs_momy(k+1,i,j,XDIR)+qflx_sgs_momy(k+1,i-1,j,XDIR)) &
                    - J13G(k-1,i,j,I_XVW)*(qflx_sgs_momy(k-1,i,j,XDIR)+qflx_sgs_momy(k-1,i-1,j,XDIR)) &
                    + J23G(k+1,i,j+1,I_XVW)*(qflx_sgs_momy(k+1,i,j+1,YDIR)+qflx_sgs_momy(k+1,i,j,YDIR)) &
                    - J23G(k-1,i,j+1,I_XVW)*(qflx_sgs_momy(k-1,i,j+1,YDIR)+qflx_sgs_momy(k-1,i,j,YDIR)) &
                    ) * 0.5_RP / ( FDZ(k)+FDZ(k-1) ) &
                  + J33G * ( qflx_sgs_momy(k,i,j  ,ZDIR) - qflx_sgs_momy(k-1,i,j,ZDIR) ) * RCDZ(k) &
                ) / GSQRT(k,i,j,I_XVW)
          end do
          end do
          end do
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             RHOT_t(k,i,j) = - ( &
                  + ( GSQRT(k,i  ,j,I_UYZ)*qflx_sgs_rhot(k,i  ,j,XDIR) &
                    - GSQRT(k,i-1,j,I_UVZ)*qflx_sgs_rhot(k,i-1,j,XDIR) ) * RCDX(i) &
                  + ( GSQRT(k,i,j  ,I_XVZ)*qflx_sgs_rhot(k,i,j  ,YDIR) &
                    - GSQRT(k,i,j-1,I_XVZ)*qflx_sgs_rhot(k,i,j-1,YDIR) ) * RCDY(j) &
                  + ( GSQRT(k  ,i,j,I_XVW)*qflx_sgs_rhot(k  ,i,j,ZDIR) &
                    - GSQRT(k-1,i,j,I_XVW)*qflx_sgs_rhot(k-1,i,j,ZDIR) ) * RCDZ(k) &
                  * ( J13G(k,i,j,I_XYZ) + J23G(k,i,j,I_XYZ) + J33G ) / GSQRT(k,i,j,I_XYZ) &
                ) / GSQRT(k,i,j,I_XYZ)
          end do
          end do
          end do
          do iq = 1, QA
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             QTRC_t(k,i,j,iq) = - ( &
                  + ( GSQRT(k,i  ,j,I_UYZ)*qflx_sgs_qtrc(k,i  ,j,iq,XDIR) &
                    - GSQRT(k,i-1,j,I_UYZ)*qflx_sgs_qtrc(k,i-1,j,iq,XDIR) ) * RCDX(i) &
                  + ( GSQRT(k,i,j  ,I_XVZ)*qflx_sgs_qtrc(k,i,j  ,iq,YDIR) &
                    - GSQRT(k,i,j-1,I_XVZ)*qflx_sgs_qtrc(k,i,j-1,iq,YDIR) ) * RCDY(j) &
                  + ( GSQRT(k  ,i,j,I_XYW)*qflx_sgs_qtrc(k  ,i,j,iq,ZDIR) &
                    - GSQRT(k-1,i,j,I_XYW)*qflx_sgs_qtrc(k-1,i,j,iq,ZDIR) ) * RCDZ(k) &
                  * ( J13G(k,i,j,I_XYZ) + J23G(k,i,j,I_XYZ) + J33G ) / GSQRT(k,i,j,I_XYZ) &
                ) / GSQRT(k,i,j,I_XYZ)
          end do
          end do
          end do
          end do
       end do
       end do

       if ( present(history_flag) ) then
       if ( history_flag ) then
          call HIST_in( tke(:,:,:), 'TKE',  'turburent kinetic energy', 'm2/s2', dttb )
          call HIST_in( nu (:,:,:), 'NU',   'eddy viscosity',           'm2/s',  dttb )
          call HIST_in( Pr (:,:,:), 'Pr',   'Prantle number',           'NIL',   dttb )
          call HIST_in( Ri (:,:,:), 'Ri',   'Richardson number',        'NIL',   dttb )

          call HIST_in( qflx_sgs_momz(:,:,:,ZDIR), 'SGS_ZFLX_MOMZ',   'SGS Z FLUX of MOMZ', 'kg/m/s2', dttb, zdim='half')
          call HIST_in( qflx_sgs_momz(:,:,:,XDIR), 'SGS_XFLX_MOMZ',   'SGS X FLUX of MOMZ', 'kg/m/s2', dttb, xdim='half')
          call HIST_in( qflx_sgs_momz(:,:,:,YDIR), 'SGS_YFLX_MOMZ',   'SGS Y FLUX of MOMZ', 'kg/m/s2', dttb, ydim='half')

          call HIST_in( qflx_sgs_momx(:,:,:,ZDIR), 'SGS_ZFLX_MOMX',   'SGS Z FLUX of MOMX', 'kg/m/s2', dttb, zdim='half')
          call HIST_in( qflx_sgs_momx(:,:,:,XDIR), 'SGS_XFLX_MOMX',   'SGS X FLUX of MOMX', 'kg/m/s2', dttb, xdim='half')
          call HIST_in( qflx_sgs_momx(:,:,:,YDIR), 'SGS_YFLX_MOMX',   'SGS Y FLUX of MOMX', 'kg/m/s2', dttb, ydim='half')

          call HIST_in( qflx_sgs_momy(:,:,:,ZDIR), 'SGS_ZFLX_MOMY',   'SGS Z FLUX of MOMY', 'kg/m/s2', dttb, zdim='half')
          call HIST_in( qflx_sgs_momy(:,:,:,XDIR), 'SGS_XFLX_MOMY',   'SGS X FLUX of MOMY', 'kg/m/s2', dttb, xdim='half')
          call HIST_in( qflx_sgs_momy(:,:,:,YDIR), 'SGS_YFLX_MOMY',   'SGS Y FLUX of MOMY', 'kg/m/s2', dttb, ydim='half')

          call HIST_in( qflx_sgs_rhot(:,:,:,ZDIR), 'SGS_ZFLX_RHOT',   'SGS Z FLUX of RHOT', 'kg K/m2/s', dttb, zdim='half')
          call HIST_in( qflx_sgs_rhot(:,:,:,XDIR), 'SGS_XFLX_RHOT',   'SGS X FLUX of RHOT', 'kg K/m2/s', dttb, xdim='half')
          call HIST_in( qflx_sgs_rhot(:,:,:,YDIR), 'SGS_YFLX_RHOT',   'SGS Y FLUX of RHOT', 'kg K/m2/s', dttb, ydim='half')

          if ( I_QV > 0 ) then
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QV,ZDIR), 'SGS_ZFLX_QV',   'SGS Z FLUX of QV', 'kg/m2 s', dttb, zdim='half')
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QV,XDIR), 'SGS_XFLX_QV',   'SGS X FLUX of QV', 'kg/m2 s', dttb, xdim='half')
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QV,YDIR), 'SGS_YFLX_QV',   'SGS Y FLUX of QV', 'kg/m2 s', dttb, ydim='half')
          endif

          if ( I_QC > 0 ) then
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QC,ZDIR), 'SGS_ZFLX_QC',   'SGS Z FLUX of QC', 'kg/m2 s', dttb, zdim='half')
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QC,XDIR), 'SGS_XFLX_QC',   'SGS X FLUX of QC', 'kg/m2 s', dttb, xdim='half')
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QC,YDIR), 'SGS_YFLX_QC',   'SGS Y FLUX of QC', 'kg/m2 s', dttb, ydim='half')
          endif

          if ( I_QR > 0 ) then
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QR,ZDIR), 'SGS_ZFLX_QR',   'SGS Z FLUX of QR', 'kg/m2 s', dttb, zdim='half')
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QR,XDIR), 'SGS_XFLX_QR',   'SGS X FLUX of QR', 'kg/m2 s', dttb, xdim='half')
             call HIST_in( qflx_sgs_qtrc(:,:,:,I_QR,YDIR), 'SGS_YFLX_QR',   'SGS Y FLUX of QR', 'kg/m2 s', dttb, ydim='half')
          endif

       end if
       end if

    end if

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE-1
       MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) + MOMZ_t(k,i,j)
    end do
    end do
    end do
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       MOMX_tp(k,i,j) = MOMX_tp(k,i,j) + MOMX_t(k,i,j)
       MOMY_tp(k,i,j) = MOMY_tp(k,i,j) + MOMY_t(k,i,j)
       RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + RHOT_t(k,i,j)
    end do
    end do
    end do

    do iq = 1, QA
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QTRC_tp(k,i,j,iq) = QTRC_tp(k,i,j,iq) + QTRC_t(k,i,j,iq)
    end do
    end do
    end do
    end do

    return
  end subroutine ATMOS_PHY_TB_driver


end module mod_atmos_phy_tb_driver
