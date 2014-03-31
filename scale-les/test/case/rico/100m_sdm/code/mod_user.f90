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
!!
!<
!-------------------------------------------------------------------------------
module mod_user
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
  integer,  private, save :: USER_LS_FLG = 0 !-- 0->no force, 1->dycoms, 2->rico

  real(RP), private, allocatable :: Q_rate(:,:,:,:)
  real(RP), private, save :: corioli
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use mod_process, only: &
       PRC_MPIstop
    use mod_const, only: &
       PI => CONST_PI
    use mod_grid, only: &
       CZ => GRID_CZ, &
       FZ => GRID_FZ
    implicit none

    namelist / PARAM_USER / &
       USER_LSsink_D,          &
       USER_LSsink_bottom,     &
       USER_LS_FLG

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

    if( USER_LS_FLG == 0 ) then  ! no large scale sinking

       MOMZ_LS(:,:) = 0.0_RP
       MOMZ_LS_DZ(:,:) = 0.0_RP
       MOMZ_LS_FLG( : ) = .false.
       QV_LS(:,:) = 0.0_RP
       Q_rate( :,:,:,: ) = 0.0_RP
       V_GEOS(:) = 0.0_RP
       U_GEOS(:) = 0.0_RP
       corioli = 0.0_RP

    elseif( USER_LS_FLG == 1 ) then ! DYCOMS

       do k = KS, KE
          if ( CZ(k) < 0.0_RP ) then
             MOMZ_LS(k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          else if ( CZ(k) < USER_LSsink_bottom ) then
             zovzb = CZ(k) / USER_LSsink_bottom
             MOMZ_LS(k,1) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom &
                  * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,1) = - 0.5_RP * USER_LSsink_d &
                  * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS(k,1) = - USER_LSsink_D * ( CZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,1) = - USER_LSsink_D
          end if
       enddo
       do k = KS-1, KE
          if ( FZ(k) < 0.0_RP ) then
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          else if ( FZ(k) < USER_LSsink_bottom ) then
             zovzb = FZ(k) / USER_LSsink_bottom
             MOMZ_LS(k,2) = - 0.5_RP * USER_LSsink_D * USER_LSsink_bottom &
                  * ( zovzb - sin(PI*zovzb)/PI )
             MOMZ_LS_DZ(k,2) = - 0.5_RP * USER_LSsink_D &
                  * ( 1.0_RP - cos(PI*zovzb) )
          else
             MOMZ_LS(k,2) = - USER_LSsink_D * ( FZ(k) - USER_LSsink_bottom * 0.5_RP )
             MOMZ_LS_DZ(k,2) = - USER_LSsink_D
          end if
       enddo
       MOMZ_LS_FLG( : ) = .true.
       QV_LS(:,:) = 0.0_RP
       Q_rate( :,:,:,: ) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = -5.5_RP
          U_GEOS(k) =  7.0_RP
       enddo
       corioli = 7.6E-5_RP

    elseif( USER_LS_FLG == 2 ) then ! RICO

       do k = KS, KE
          if ( CZ(k) < 0.0_RP ) then
             MOMZ_LS(k,1) = 0.0_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          else if( CZ(k) < 2000.0_RP ) then
             MOMZ_LS(k,1) = - 5.0E-3_RP / 2260.0_RP * CZ(k)
             MOMZ_LS_DZ(k,1) = - 5.0E-3_RP / 2260.0_RP
             !--- fitted by circle
          else if( CZ(k) < 2520.0_RP ) then
             xi   = 2260.0_RP
             ar   = -5.0E-3_RP/xi
             cr   = ar* ( 2520.0_RP-xi )/( 1.0_RP - sqrt( 1.0_RP + ar*ar ) )
             xz   = CZ(k)
             yc   = ar * xi + cr
             MOMZ_LS(k,1) = ar * xi + cr  &
                  - sqrt( cr*cr - ( xz - xi - cr / ar * ( 1.0_RP - sqrt( 1.0_RP + ar*ar ) ) )**2.0_RP )
             MOMZ_LS_DZ(k,1) = - ( xz - 2520._RP )/( MOMZ_LS(k,1) - yc )
          else
             MOMZ_LS(k,1) = - 5.0E-3_RP
             MOMZ_LS_DZ(k,1) = 0.0_RP
          end if
          if( CZ(k) < 2980.0_RP ) then
             QV_LS(k,1) = ( -1.0_RP + 1.3456_RP / 2980.0_RP  &
                  * CZ(k) ) * 1.E-3_RP / 86400.0_RP  !--- [kg/kg/s]
          else
             QV_LS(k,1) = 4.0_RP * 1.E-6_RP * 1.E-3_RP !--- [kg/kg/s]
          endif
       enddo
       do k = KS-1, KE
          if ( FZ(k) < 0.0_RP ) then
             MOMZ_LS(k,2) = 0.0_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          else if( FZ(k) < 2000.0_RP ) then
             MOMZ_LS(k,2) = - 5.0E-3_RP / 2260.0_RP * FZ(k)
             MOMZ_LS_DZ(k,2) = - 5.0E-3_RP / 2260.0_RP
          else if( FZ(k) < 2520.0_RP ) then
             xi   = 2260.0_RP
             ar   = -5.0E-3_RP/xi
             cr   = ar* ( 2520.0_RP-xi )/( 1.0_RP - sqrt( 1.0_RP + ar*ar ) )
             xz   = FZ(k)
             yc   = ar * xi + cr
             MOMZ_LS(k,2) = ar * xi + cr  &
                  - sqrt( cr*cr - ( xz - xi - cr / ar * ( 1.0_RP - sqrt( 1.0_RP + ar*ar ) ) )**2.0_RP )
             MOMZ_LS_DZ(k,2) = - ( xz - 2520._RP )/( MOMZ_LS(k,2) - yc )
          else
             MOMZ_LS(k,2) = - 5.0E-3_RP
             MOMZ_LS_DZ(k,2) = 0.0_RP
          end if
       enddo
       MOMZ_LS_FLG( : ) = .false.
       MOMZ_LS_FLG( I_RHOT ) = .true.
       MOMZ_LS_FLG( I_QTRC ) = .true.
       Q_rate( :,:,:,: ) = 0.0_RP

       do k = KS-1, KE
          V_GEOS(k) = -3.8_RP
          U_GEOS(k) = -9.9_RP + 2.0E-3_RP * CZ(k)
       enddo
       corioli = 4.5E-5_RP

    endif

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
         MOMZ_tp, &
         MOMX_tp, &
         MOMY_tp, &
         RHOT_tp, &
         QTRC_tp
    use mod_grid, only: &
         RCDZ => GRID_RCDZ, &
         RFDZ => GRID_RFDZ
    use mod_const, only: &
         CPdry => CONST_CPdry, &
         Rdry  => CONST_Rdry, &
         Rvap  => CONST_Rvap, &
         LH0   => CONST_LH0, &
         P00   => CONST_PRE00 
    use mod_atmos_thermodyn, only: &
         CPw   => AQ_CP 
 
    implicit none
    !---------------------------------------------------------------------------
    real(RP) :: ratesum

    real(RP) :: WORK(KA,IA,JA)
    real(RP) :: VELX(KA,IA,JA), VELY(KA,IA,JA)
    real(RP) :: TEMP_t(KA,IA,JA) ! tendency rho*theta     [K*kg/m3/s]
    real(RP) :: d_PT, drhot, LWPT, LWPT_a
    real(RP) :: QHYD, RovCP, PRES, CPovCV, CPtot, Rtot, Qdry, Qtot

    integer :: k, i, j, iq
    integer :: IIS, IIE, JJS, JJE


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
                  - MOMZ_LS(k,2) * ( WORK(k+1,i,j) - WORK(k,i,j) ) * RCDZ(k)
          enddo
          enddo
          enddo
          !$omp parallel do private(i,j,k) schedule(static,1) collapse(2)
          do j = JJS, JJE
          do i = IIS, IIE
             MOMZ_tp(KE-1,i,j) = MOMZ_tp(KE-1,i,j) &
                  - MOMZ_LS(KE-1,2) * (           - WORK(KE-1,i,j) ) * RCDZ(KE-1)
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
          if( USER_LS_FLG == 2 ) then ! RICO
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
          endif


          do iq = 1, QA
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
    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Parametarized Radiation of RICO'
    do JJS = JS, JE, JBLOCK
    JJE = JJS+JBLOCK-1
    do IIS = IS, IE, IBLOCK
    IIE = IIS+IBLOCK-1

       do k = KS, KE

            QDRY  = 1.0_RP
            CPTOT = 0.0_RP
            QHYD  = 0.0_RP

            do iq = QQS, QQE
               QDRY  = QDRY - QTRC(k,i,j,iq)
               CPTOT = CPTOT + QTRC(k,i,j,iq) * CPw(iq)
            enddo

            do iq = QWS, QWE
               QHYD = QHYD + QTRC(k,i,j,iq)
            enddo

            RTOT   = Rdry *QDRY + Rvap*QTRC(k,i,j,I_QV)
            CPTOT  = CPdry*QDRY + CPTOT
            CPovCV = CPTOT / ( CPTOT - RTOT )
            PRES   = P00 * ( RHOT(k,i,j) * RTOT / P00 )**CPovCV
            RovCP  = RTOT / CPTOT

            d_PT   = RHOT(k,i,j) / DENS(k,i,j)
            LWPT = RHOT(k,i,j) / DENS(k,i,j)  &
                 - ( LH0 / CPdry * QHYD ) * ( P00 / PRES )**RovCP
            LWPT_a = LWPT - 2.5_RP / 86400.0_RP ! * dt
            drhot = ( LWPT_a &
                 + ( LH0 / CPdry * QHYD ) * ( P00 / PRES )**RovCP &
                 ) * DENS(k,i,j) - RHOT(k,i,j)

            RHOT_tp(k,i,j) = RHOT_tp(k,i,j) + drhot !/ dt
       enddo

    enddo
    enddo

    return
  end subroutine USER_step

end module mod_user
