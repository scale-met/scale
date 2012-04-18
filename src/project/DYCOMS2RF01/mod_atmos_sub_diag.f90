!-------------------------------------------------------------------------------
!> module Atmosphere / Diagnostic Variables
!!
!! @par Description
!!          calc diagnostic variables
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-04-18 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_atmos_sub_diag
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  use mod_stdio, only: &
     IO_FID_LOG,  &
     IO_L
  use mod_time, only: &
     TIME_rapstart, &
     TIME_rapend
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DIAG_setup
  public :: ATMOS_DIAG_history

  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_index.h'
  include 'inc_tracer.h'

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
  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DIAG_setup
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[Diagnostics]/Categ[ATMOS]'

    return
  end subroutine ATMOS_DIAG_setup

  !-----------------------------------------------------------------------------
  ! Parametarized Radiative heating
  !-----------------------------------------------------------------------------
  subroutine ATMOS_DIAG_history
    use mod_process, only : &
       PRC_nmax,   &
       PRC_myrank
    use mod_const, only: &
       Rdry   => CONST_Rdry, &
       CPdry  => CONST_CPdry, &
       RovCP  => CONST_RovCP, &
       CPovCV => CONST_CPovCV, &
       LH0    => CONST_LH0, &
       P00    => CONST_PRE00
    use mod_time, only: &
       TIME_DTSEC
    use mod_grid, only : &
       CZ  => GRID_CZ,  &
       CDZ => GRID_CDZ
    use mod_geometrics, only: &
       area    => GEOMETRICS_area,    &
       totarea => GEOMETRICS_totarea
    use mod_history, only: &
       HIST_in
    use mod_atmos_vars, only: &
       DENS, &
       MOMZ, &
       MOMX, &
       MOMY, &
       RHOT, &
       QTRC
    implicit none

    real(8) :: QT    (KA,IA,JA) ! total  water mixing ratio [g/kg]
    real(8) :: QL    (KA,IA,JA) ! liquid water mixing ratio [g/kg]
    real(8) :: LWPT  (KA,IA,JA) ! liquid water potential temperature [K]
    real(8) :: LWPT_v(KA,IA,JA) ! variance of liquid water potential temperature [K]
    real(8) :: VELZ_v(KA,IA,JA) ! variance of velocity w [m/s]
    real(8) :: VELX_v(KA,IA,JA) ! variance of velocity u [m/s]
    real(8) :: VELY_v(KA,IA,JA) ! variance of velocity v [m/s]
    real(8) :: Zb    (1,IA,JA)  ! cloud-base height [m]
    real(8) :: LWP   (1,IA,JA)  ! liquid water path [g/m2]

    real(8) :: VELZ(KA,IA,JA) ! velocity w [m/s]
    real(8) :: VELX(KA,IA,JA) ! velocity u [m/s]
    real(8) :: VELY(KA,IA,JA) ! velocity v [m/s]
    real(8) :: PRES(KA,IA,JA) ! pressure [Pa]

    real(8) :: statval   (KA,4,0:PRC_nmax-1)
    real(8) :: allstatval(KA,4)

    integer :: ierr
    integer :: k, i, j, p
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** special history: diagnostic variables'

    do j = JS, JE
    do i = IS, IE

       do k = KS, KE
          QL(k,i,j) = ( QTRC(k,i,j,I_QC) &
                      + QTRC(k,i,j,I_QR) ) * 1.D3 ! [kg/kg->g/kg]

          QT(k,i,j) = QTRC(k,i,j,I_QV) * 1.D3 & ! [kg/kg->g/kg]
                    + QL  (k,i,j)

       enddo

       LWP(1,i,j) = 0.D0
       do k = KS, KE
          LWP(1,i,j) = LWP(1,i,j) + QL(k,i,j) * DENS(k,i,j) * CDZ(k)
       enddo

       ! diagnose cloud base
       Zb(1,i,j) = 0.D0
       do k = KS, KE
          if( QL(k,i,j) >= 1.D-2 ) exit ! cloud base >= 0.01[g/kg]
          Zb(1,i,j) = CZ(k)
       enddo

    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       VELZ(k,i,j) = 0.5D0 * ( MOMZ(k-1,i,j)+MOMZ(k,i,j) ) / DENS(k,i,j)
       VELX(k,i,j) = 0.5D0 * ( MOMX(k,i-1,j)+MOMX(k,i,j) ) / DENS(k,i,j)
       VELY(k,i,j) = 0.5D0 * ( MOMY(k,i,j-1)+MOMY(k,i,j) ) / DENS(k,i,j)
       PRES(k,i,j) = P00 * ( RHOT(k,i,j) * Rdry / P00 )**CPovCV
       LWPT(k,i,j) = RHOT(k,i,j) / DENS(k,i,j) - ( LH0 / CPdry * QL(k,i,j)*1.D-3 ) * ( P00 / PRES(k,i,j) )**RovCP
    enddo
    enddo
    enddo

    ! horizontal average
    statval(:,:,PRC_myrank) = 0.D0
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          statval(k,1,PRC_myrank) = statval(k,1,PRC_myrank) + LWPT(k,i,j) * area(1,i,j)
          statval(k,2,PRC_myrank) = statval(k,2,PRC_myrank) + VELZ(k,i,j) * area(1,i,j)
          statval(k,3,PRC_myrank) = statval(k,3,PRC_myrank) + VELX(k,i,j) * area(1,i,j)
          statval(k,4,PRC_myrank) = statval(k,4,PRC_myrank) + VELY(k,i,j) * area(1,i,j)
       enddo
    enddo
    enddo
    do k = KS, KE
       statval(k,1,PRC_myrank) = statval(k,1,PRC_myrank) / totarea
       statval(k,2,PRC_myrank) = statval(k,2,PRC_myrank) / totarea
       statval(k,3,PRC_myrank) = statval(k,3,PRC_myrank) / totarea
       statval(k,4,PRC_myrank) = statval(k,4,PRC_myrank) / totarea
    enddo

    ! MPI broadcast
    call TIME_rapstart('COMM Bcast MPI')
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,1,p),       &
                       KA*4,                 &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
    enddo
    call TIME_rapend  ('COMM Bcast MPI')

    allstatval(:,:) = 0.D0
    do p = 0, PRC_nmax-1
       do k = KS, KE
          allstatval(k,1) = allstatval(k,1) + statval(k,1,p)
          allstatval(k,2) = allstatval(k,2) + statval(k,2,p)
          allstatval(k,3) = allstatval(k,3) + statval(k,3,p)
          allstatval(k,4) = allstatval(k,4) + statval(k,4,p)
       enddo
    enddo
    do k = KS, KE
       allstatval(k,1) = allstatval(k,1) / PRC_nmax
       allstatval(k,2) = allstatval(k,2) / PRC_nmax
       allstatval(k,3) = allstatval(k,3) / PRC_nmax
       allstatval(k,4) = allstatval(k,4) / PRC_nmax
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       LWPT_v(k,i,j) = LWPT(k,i,j) - allstatval(k,1)
       VELZ_v(k,i,j) = VELZ(k,i,j) - allstatval(k,2)
       VELX_v(k,i,j) = VELX(k,i,j) - allstatval(k,3)
       VELY_v(k,i,j) = VELY(k,i,j) - allstatval(k,4)
    enddo
    enddo
    enddo

    call HIST_in( QT    (:,:,:), 'QT',     'Total water',             'g/kg', '3D', TIME_DTSEC )
    call HIST_in( QL    (:,:,:), 'QL',     'Liquid water',            'g/kg', '3D', TIME_DTSEC )
    call HIST_in( LWP   (:,:,:), 'LWP',    'Liquid water path',       'g/m2', '2D', TIME_DTSEC )

    call HIST_in( LWPT  (:,:,:), 'LWPT',   'Liquid water pot. temp.', 'K',    '3D', TIME_DTSEC )

    call HIST_in( LWPT_v(:,:,:), 'LWPT_v', 'Variance of LWPT',        'K',    '3D', TIME_DTSEC )
    call HIST_in( VELZ_v(:,:,:), 'VELZ_v', 'Variance of W',           'm/s',  '3D', TIME_DTSEC )
    call HIST_in( VELX_v(:,:,:), 'VELX_v', 'Variance of U',           'm/s',  '3D', TIME_DTSEC )
    call HIST_in( VELY_v(:,:,:), 'VELY_v', 'Variance of V',           'm/s',  '3D', TIME_DTSEC )

    call HIST_in( Zb    (:,:,:), 'Zb',     'cloud base height',       'm',    '2D', TIME_DTSEC )

    return
  end subroutine ATMOS_DIAG_history

end module mod_atmos_sub_diag
