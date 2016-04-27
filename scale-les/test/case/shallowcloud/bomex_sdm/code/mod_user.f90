!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2013-06-20 (S.Nishizawa)   [new] split from dynamical core
!! @li      2014-12-12 (Y.Sato)        [mod] Modify several bug in forcing term for RICO case
!!
!<
!-------------------------------------------------------------------------------
module mod_user
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
  public :: USER_setup
  public :: USER_step

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
  real(RP), private, allocatable :: MOMZ_LS(:,:)
  real(RP), private, allocatable :: MOMZ_LS_DZ(:,:)
  real(RP), private, allocatable :: QV_LS(:,:)
  logical,  private, save        :: MOMZ_LS_FLG(6)
  real(RP), private, allocatable :: U_GEOS(:)
  real(RP), private, allocatable :: V_GEOS(:)

  real(RP), private, save :: USER_LSsink_D = 0.0_RP         ! large scale sinking parameter [1/s]
  real(RP), private, save :: USER_LSsink_bottom = 0.0_RP    ! large scale sinking parameter [m]
  ! flag
  real(RP), private, allocatable :: Q_rate(:,:,:,:)
  real(RP), private, save :: corioli
  real(RP), private, allocatable :: dQrad(:)

  !--- for surface flux
  real(RP), private, save      :: Cm_min  =    0.0E-5_RP ! minimum bulk coef. of u,v,w
  real(RP), private, parameter :: Cm_max  =    2.5E-3_RP ! maximum bulk coef. of u,v,w
  real(RP), private, save      :: U_minM  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxM  =  100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, save      :: U_minE  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxE  =  100.0_RP  ! maximum U_abs for u,v,w
  real(RP), private, save      :: U_minH  =    0.0_RP  ! minimum U_abs for u,v,w
  real(RP), private, parameter :: U_maxH  =  100.0_RP  ! maximum U_abs for u,v,w
  !-----------------------------------------------------------------------------
  real(RP), private, save      :: USER_Ustar = 0.28_RP ! constant frection velocity [m/s]
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       PI => CONST_PI
    use scale_grid, only: &
       CZ => GRID_CZ, &
       FZ => GRID_FZ
    implicit none

    namelist / PARAM_USER / &
       USER_LSsink_D,          &
       USER_LSsink_bottom,     &
       USER_Ustar

    real(RP) :: zovzb
    real(RP) :: cr, xi, ar, yc, xz

    integer :: k
    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    allocate( MOMZ_LS(KA,2) )
    allocate( MOMZ_LS_DZ(KA,2) )
    allocate( QV_LS(KA,2) )
    allocate( V_GEOS(KA) )
    allocate( U_GEOS(KA) )
    allocate( Q_rate( KA,IA,JA,QA ) )
    allocate( dQrad( KA ) )

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

       do k = KS, KE
          if ( CZ(k) < 0.0_RP ) then
             MOMZ_LS(k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          else if( CZ(k) < 1500.0_RP ) then
             MOMZ_LS(k,1) = - 6.5E-3_RP / 1500.0_RP * CZ(k)
             MOMZ_LS_DZ(k,1) = - 6.5E-3_RP / 1500.0_RP
             !--- fitted by circle
          else if( CZ(k) < 2100.0_RP ) then
             MOMZ_LS(k,1) = - 6.5E-3_RP + 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP ) * ( CZ(k) - 1500.0_RP )
             MOMZ_LS_DZ(k,1) = 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP )
          else
             MOMZ_LS(k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          end if

          if( CZ(k) < 300.0_RP ) then
             QV_LS(k,1) = -1.2E-8_RP
          elseif( CZ(k) < 500.0_RP ) then
             QV_LS(k,1) = - 1.2E-8_RP + 1.2E-8_RP * ( CZ(k) - 300.0_RP ) / 200.0_RP  !--- [kg/kg/s]
          else
             QV_LS(k,1) = 0.0_RP !--- [kg/kg/s]
          endif
       enddo
       do k = KS-1, KE
          if ( FZ(k) < 0.0_RP ) then
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          else if( FZ(k) < 1500.0_RP ) then
             MOMZ_LS(k,2) = - 6.5E-3_RP / 1500.0_RP * FZ(k)
             MOMZ_LS_DZ(k,2) = - 6.5E-3_RP / 1500.0_RP
          else if( FZ(k) < 2100.0_RP ) then
             MOMZ_LS(k,2) = - 6.5E-3_RP + 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP ) * ( FZ(k) - 1500.0_RP )
             MOMZ_LS_DZ(k,1) = 6.5E-3_RP / ( 2100.0_RP - 1500.0_RP )
          else
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          end if
       enddo

       MOMZ_LS_FLG( : ) = .false.
       MOMZ_LS_FLG( I_MOMX ) = .true.
       MOMZ_LS_FLG( I_MOMY ) = .true.
       MOMZ_LS_FLG( I_RHOT ) = .true.
       MOMZ_LS_FLG( I_QTRC ) = .true.
       Q_rate( :,:,:,: ) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = 0.0_RP
          U_GEOS(k) = -10.0_RP + 1.8E-3_RP * CZ(k)
       enddo
       corioli = 3.76E-5_RP

    do k = KS, KE
      if( CZ(k) < 1500.0_RP ) then
        dQrad(k) = - 2.315E-5_RP
      elseif( CZ(k) < 2500.0_RP ) then
        dQrad(k) = - 2.315E-5_RP + 2.315E-5_RP / 1000.0_RP * ( CZ(k) - 1500.0_RP )
      else
        dQrad(k) = 0.0_RP
      endif
    enddo

    call USER_step

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    use mod_atmos_vars, only: &
         DENS,    &
         MOMZ,    &
         MOMX,    &
         MOMY,    &
         RHOT,    &
         QTRC,    &
         DENS_tp, &
         MOMZ_tp, &
         MOMX_tp, &
         MOMY_tp, &
         RHOT_tp, &
         QTRC_tp
    use scale_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use scale_const, only: &
         CPdry => CONST_CPdry, &
         Rdry  => CONST_Rdry, &
         Rvap  => CONST_Rvap, &
         LH0   => CONST_LH0, &
         P00   => CONST_PRE00 
    use scale_atmos_thermodyn, only: &
         CPw   => AQ_CP, &
         ATMOS_THERMODYN_templhv
    use scale_history, only: &
         HIST_in
    use scale_time, only: &
         do_phy_sf => TIME_DOATMOS_PHY_SF, &
         dtsf =>  TIME_DTSEC_ATMOS_PHY_SF
#ifdef _SDM
    use scale_atmos_phy_mp_sdm, only: &
         ATMOS_PHY_MP_sdm_Mixingratio
    use scale_tracer_sdm
#endif
 
    implicit none
    !---------------------------------------------------------------------------
    real(RP) :: ratesum

    real(RP) :: SFLX_MOMZ(IA,JA)
    real(RP) :: SFLX_MOMX(IA,JA)
    real(RP) :: SFLX_MOMY(IA,JA)
    real(RP) :: SFLX_POTT(IA,JA)
    real(RP) :: SFLX_QV  (IA,JA)
    real(RP) :: SHFLX(IA,JA)
    real(RP) :: LHFLX(IA,JA)
    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA), VELY(KA,IA,JA)
    real(RP) :: TEMP_t(KA,IA,JA) ! tendency rho*theta     [K*kg/m3/s]
    real(RP) :: QHYD, RovCP, PRES, CPovCV, CPtot, Rtot, Qdry, Qtot
    real(RP) :: qv_evap, pres_evap, PT_SST, lhv, TEMP
#ifdef _SDM
    real(RP) :: Qe(KA,IA,JA, MP_QA) ! tendency rho*theta     [K*kg/m3/s]
#endif

    integer :: k, i, j, iq, iw
    integer :: IIS, IIE, JJS, JJE
    ! work
    real(RP) :: Uabs  ! absolute velocity at the lowermost atmos. layer [m/s]
    real(RP) :: Cm, Ch, Ce    !


    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       if ( MOMZ_LS_FLG(I_MOMZ) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             WORK(k,i,j) = MOMZ(k,i,j) * 2.0_RP / ( DENS(k+1,i,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-2
             MOMZ_tp(k,i,j) = MOMZ_tp(k,i,j) &
                  - MOMZ_LS(k,2) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RCDZ(k)*0.0_RP
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMZ_tp(KE-1,i,j) = MOMZ_tp(KE-1,i,j) &
                  - MOMZ_LS(KE-1,2) * (           - WORK(KE-1,i,j) ) * RCDZ(KE-1)*0.0_RP
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_MOMX) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS-1, JJE
          do i = IIS,   IIE+1
          do k = KS, KE
             VELY(k,i,j) = 2.0_RP * MOMY(k,i,j) / ( DENS(k,i,j+1)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMX_tp(k,i,j) = MOMX_tp(k,i,j) &
                  + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                  * ( - CORIOLI * V_GEOS(k) &
                      + CORIOLI * 0.25_RP &
                      * ( VELY(k,i,j)+VELY(k,i+1,j)+VELY(k,i,j-1)+VELY(k,i+1,j-1) ) &
                    ) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMX_tp(KE,i,j) = MOMX_tp(KE,i,j) &
                  + 0.5_RP * ( DENS(k,i+1,j)+DENS(k,i,j) ) &
                  *  ( - CORIOLI * V_GEOS(KE) &
                       + CORIOLI * 0.25_RP &
                       * ( VELY(KE,i,j)+VELY(KE,i+1,j)+VELY(KE,i,j-1)+VELY(KE,i+1,j-1) ) &
                     ) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1) 
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_MOMY) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = MOMY(k,i,j) * 2.0_RP / ( DENS(k,i,j+1) + DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE+1
          do i = IIS-1, IIE
          do k = KS, KE
             VELX(k,i,j) = MOMX(k,i,j) * 2.0_RP / ( DENS(k,i+1,j)+DENS(k,i,j) )
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             MOMY_tp(k,i,j) = MOMY_tp(k,i,j) &
                  + 0.5_RP * ( DENS(k,i,j+1)+DENS(k,i,j) )  &
                  * ( + CORIOLI * U_GEOS(k) &
                      - CORIOLI * 0.25_RP &
                      * ( VELX(k,i,j)+VELX(k,i,j+1)+VELX(k,i-1,j)+VELX(k,i-1,j+1) ) &
                    ) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k) 
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS,   JJE+1
          do i = IIS-1, IIE
             MOMY_tp(KE,i,j) = MOMY_tp(KE,i,j) &
                  + 0.5_RP * ( DENS(KE,i,j+1)+DENS(KE,i,j) ) &
                  * ( + CORIOLI * U_GEOS(KE) &
                      - CORIOLI * 0.25_RP  &
                      * ( VELX(KE,i,j)+VELX(KE,i,j+1)+VELX(KE,i-1,j)+VELX(KE,i-1,j+1) ) &
                    ) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_RHOT) ) then
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE
             WORK(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
          do k = KS, KE-1
             RHOT_tp(k,i,j) = RHOT_tp(k,i,j) &
                  - MOMZ_LS(k,1) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RFDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             RHOT_tp(KE,i,j) = RHOT_tp(KE,i,j) &
                  - MOMZ_LS(KE,1) * ( WORK(KE,i,j) - WORK(KE-1,i,j) ) * RFDZ(KE-1)
          enddo
          enddo
       end if

       if ( MOMZ_LS_FLG(I_QTRC) ) then

          !##### advection of scalar quantity #####
#ifdef _SDM
             Q_rate( :,:,:,: ) = 0.0_RP
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                   Q_rate( k,i,j,I_QV ) = 1.0_RP
             enddo
             enddo
             enddo
#else
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE
                ratesum = 0.0_RP
                do iq = 1, QQA
                   ratesum = ratesum + QTRC(k,i,j,iq)
                enddo
                do iq = 1, QQA
                   Q_rate( k,i,j,iq ) = QTRC(k,i,j,iq) / ratesum
                enddo
                do iq = QQA+1, QA
                   Q_rate( k,i,j,iq ) = 0.0_RP
                enddo
             enddo
             enddo
             enddo
#endif

#ifdef _SDM
          do iq = 1, I_QV
#else
          do iq = 1, QA
#endif
             !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
             do k = KS, KE-1
                QTRC_tp(k,i,j,iq) = QTRC_tp(k,i,j,iq) &
                     + QV_LS(k,1) * Q_rate(k,i,j,iq) &
                     - MOMZ_LS(k,1) * ( QTRC(k+1,i,j,iq) - QTRC(k,i,j,iq) ) * RFDZ(k)
             enddo
             enddo
             enddo
             !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
             do j = JJS, JJE
             do i = IIS, IIE
                QTRC_tp(KE,i,j,iq) = QTRC_tp(KE,i,j,iq) &
                     + QV_LS(KE,1) * Q_rate(KE,i,j,iq) &
                     - MOMZ_LS(KE,1) * ( QTRC(KE,i,j,iq) - QTRC(KE-1,i,j,iq) ) * RFDZ(KE-1)
             enddo
             enddo
          end do
       end if

    enddo
    enddo

    ! radiative flux
    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation of BOEMX'

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
            RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + dQrad(k) * DENS(k,i,j) ! [K/s]
       enddo

    enddo
    enddo

    if( do_phy_sf ) then

      if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Surface flux of BOEMX'

#ifdef _SDM
      call ATMOS_PHY_MP_sdm_Mixingratio( Qe, QTRC )
#endif

      do j = JS-1, JE
      do i = IS-1, IE


       ! at cell center

       !--- absolute velocity
       Uabs = sqrt( &
!              ( MOMZ(KS,i,j)                  )**2 &
            + ( MOMX(KS,i-1,j) + MOMX(KS,i,j) )**2 &
            + ( MOMY(KS,i,j-1) + MOMY(KS,i,j) )**2 &
            ) / DENS(KS,i,j) * 0.5_RP

        Cm = min( max(USER_Ustar**2 / Uabs**2, Cm_min), Cm_max )

       !--- saturation at surface
       qtot = 0.0_RP
       qdry = 1.0_RP
       CPtot = 0.0_RP
#ifdef _SDM
       do iw = I_QV, I_QV
          qdry = qdry - QTRC(KS,i,j,iw)
          qtot = qtot + QTRC(KS,i,j,iw)
          CPtot = CPtot + QTRC(KS,i,j,iw) * CPw(iw)
       enddo
       qtot = qtot + Qe(KS,i,j,I_mp_QC)
#else
       do iw = QQS, QQE
          qdry = qdry - QTRC(KS,i,j,iw)
          qtot = qtot + QTRC(KS,i,j,iw)
          CPtot = CPtot + QTRC(KS,i,j,iw) * CPw(iw)
       enddo
#endif
       Rtot   = Rdry*qdry + Rvap*QTRC(KS,i,j,I_QV)
       CPtot  = CPdry*qdry + CPtot
       CPovCV = CPtot / ( CPtot - Rtot )
       PRES   = P00 * ( RHOT(KS,i,j) * Rtot / P00 )**CPovCV
       TEMP   = (RHOT(KS,i,j)/DENS(KS,i,j)) * (PRES/P00)**(Rtot/cptot)
       call ATMOS_THERMODYN_templhv( lhv, TEMP )

       ! flux
       SFLX_MOMZ(i,j) = 0.0_RP
       SFLX_QV  (i,j) =  5.2E-5_RP
       SFLX_POTT(i,j) =  8.0E-3_RP

       ! at (u, y, layer)
       Uabs = sqrt( &
            + ( 2.0_RP *   MOMX(KS,i,j) )**2 &
            + ( 0.5_RP * ( MOMY(KS,i,j-1) + MOMY(KS,i,j) + MOMY(KS,i+1,j-1) + MOMY(KS,i+1,j) ) )**2 &
            ) / ( DENS(KS,i,j) + DENS(KS,i+1,j) )
       Cm = min( max(USER_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       SFLX_MOMX(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMX(KS,i,j)

       ! at (x, v, layer)
       Uabs = sqrt( &
            + ( 0.5_RP * ( MOMX(KS,i-1,j) + MOMX(KS,i,j) + MOMX(KS,i-1,j+1) + MOMX(KS,i,j+1) ) )**2 &
            + ( 2.0_RP *   MOMY(KS,i,j) )**2 &
          ) / ( DENS(KS,i,j) + DENS(KS,i,j+1) )
       Cm = min( max(USER_Ustar**2 / Uabs**2, Cm_min), Cm_max )
       SFLX_MOMY(i,j) = - Cm * min( max(Uabs,U_minM), U_maxM ) * MOMY(KS,i,j)

       SHFLX(i,j) = SFLX_POTT(i,j) * CPdry
       LHFLX(i,j) = SFLX_QV  (i,j) * lhv

      enddo
      enddo

      call HIST_in( SHFLX(:,:), 'SHFLX', 'sensible heat flux', 'W/m2', dtsf )
      call HIST_in( LHFLX(:,:), 'LHFLX', 'latent heat flux',   'W/m2', dtsf )

      !$omp parallel do private(i,j) OMP_SCHEDULE_ collapse(2)
      do j = JS, JE
      do i = IS, IE
       RHOT_tp(KS,i,j) = RHOT_tp(KS,i,j) &
            + ( SFLX_POTT(i,j) &
            + SFLX_QV(i,j) * RHOT(KS,i,j) / DENS(KS,i,j) &
              ) * RCDZ(KS)
       DENS_tp(KS,i,j) = DENS_tp(KS,i,j) &
             + SFLX_QV(i,j) * RCDZ(KS)
       MOMZ_tp(KS,i,j) = MOMZ_tp(KS,i,j) &
             + SFLX_MOMZ(i,j) * RFDZ(KS)
       MOMX_tp(KS,i,j) = MOMX_tp(KS,i,j) &
             + SFLX_MOMX(i,j) * RCDZ(KS)
       MOMY_tp(KS,i,j) = MOMY_tp(KS,i,j) &
             + SFLX_MOMY(i,j) * RCDZ(KS)
       QTRC_tp(KS,i,j,I_QV) = QTRC_tp(KS,i,j,I_QV) &
             + SFLX_QV(i,j) * RCDZ(KS)
      enddo
      enddo

    endif

    return
  end subroutine USER_step

end module mod_user
