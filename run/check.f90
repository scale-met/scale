!-------------------------------------------------------------------------------
!> program Check result
!<
!-------------------------------------------------------------------------------
program check
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_history, only: &
       HistoryGet
  implicit none
  include 'mpif.h'

  integer, parameter :: RP = 8

  integer, parameter :: IMAX = 16
  integer, parameter :: JMAX = 16
  integer, parameter :: KMAX = 23

  integer, parameter :: K_zi = 10
  integer, parameter :: Nt = 5

  real(RP) :: DENS(IMAX,JMAX,KMAX)
  real(RP) :: W(IMAX,JMAX,KMAX)
  real(RP) :: U(IMAX,JMAX,KMAX)
  real(RP) :: PT(IMAX,JMAX,KMAX)
  real(RP) :: QV(IMAX,JMAX,KMAX)
  real(RP) :: QC(IMAX,JMAX,KMAX)
  real(RP) :: QR(IMAX,JMAX,KMAX)

  real(RP) :: W2(IMAX,JMAX,KMAX)

  real(RP) :: DENSm(KMAX)
  real(RP) :: MOMZm(KMAX)
  real(RP) :: MOMXm(KMAX)
  real(RP) :: RHOTm(KMAX)
  real(RP) :: RHOQVm(KMAX)
  real(RP) :: RHOQCm(KMAX)
  real(RP) :: RHOQRm(KMAX)

  real(RP) :: W2m(KMAX)


  character(len=128) :: bname

  integer :: PRC_nmax, PRC_myrank
  integer :: ierr

  integer :: i, j, k

!  if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
!     write(*,*) ' xxx Program needs history base name!.'
!     stop
!  end if
!  call get_command_argument(1, bname)
  bname = 'history_23x16x16'


  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, PRC_nmax, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, PRC_myrank, ierr)


  call HistoryGet( DENS(:,:,:), bname, 'DENS', Nt, PRC_myrank )
  call HistoryGet( W   (:,:,:), bname, 'W'   , Nt, PRC_myrank )
  call HistoryGet( U   (:,:,:), bname, 'U'   , Nt, PRC_myrank )
  call HistoryGet( PT  (:,:,:), bname, 'PT'  , Nt, PRC_myrank )
  call HistoryGet( QV  (:,:,:), bname, 'QV'  , Nt, PRC_myrank )
  call HistoryGet( QC  (:,:,:), bname, 'QC'  , Nt, PRC_myrank )
  call HistoryGet( QR  (:,:,:), bname, 'QR'  , Nt, PRC_myrank )

  call mean( DENSm(:), DENS(:,:,:) )
  call mean( MOMZm(:), W(:,:,:)*DENS(:,:,:) )
  call mean( MOMXm(:), U(:,:,:)*DENS(:,:,:) )
  call mean( RHOTm(:), PT(:,:,:)*DENS(:,:,:) )
  call mean( RHOQVm(:), QV(:,:,:)*DENS(:,:,:) )
  call mean( RHOQCm(:), QC(:,:,:)*DENS(:,:,:) )
  call mean( RHOQRm(:), QR(:,:,:)*DENS(:,:,:) )

  do k = 1, K_zi
  do j = 1, JMAX
  do i = 1, IMAX
     W2(i,j,k) = ( W (i,j,k) - MOMZm(k) / DENSm(k) ) ** 2
  end do
  end do
  end do


  call mean( W2m(:), W2(:,:,:)*DENS(:,:,:) )

  if ( PRC_myrank == 0 ) then

     call assert( W2m(:)   /DENSm(:), 0.240_RP,   0.15_RP )
     call assert( MOMXm(:) /DENSm(:), 4.920_RP,   2.5E-3_RP )
     call assert( RHOTm(:) /DENSm(:), 284.996_RP, 7.0E-5_RP )
     call assert( RHOQVm(:)/DENSm(:), 6.43E-3_RP, 4.0E-3_RP )
     call assert( RHOQCm(:)/DENSm(:), 2.47E-5_RP, 0.15_RP )
     call assert( RHOQRm(:)/DENSm(:), 3.1E-8_RP,  1.0_RP )

  end if

  call MPI_Finalize(ierr)
  stop

contains

  subroutine mean( PHIm, PHI )
    real(RP), intent(out) :: PHIm(KMAX)
    real(RP), intent(in)  :: PHI(IMAX,JMAX,KMAX)

    real(RP) :: work(KMAX,0:PRC_nmax-1)

    integer :: ierr
    integer :: i, j, k, p

    work(:,:) = 0.D0
    do k = 1, KMAX
    do j = 1, JMAX
    do i = 1, IMAX
       work(k,PRC_myrank) = work(k,PRC_myrank) + PHI(i,j,k)
    end do
    end do
    end do

    do p = 0, PRC_nmax-1
       call MPI_Bcast( work(:,p),   &
                       KMAX,                 &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
    end do

    PHIm(:) = 0.D0
    do p = 0, PRC_nmax-1
       do k = 1, KMAX
          PHIm(k) = PHIm(k) + work(k,p)
       end do
    end do

    PHIm(:) = PHIm(:) / (IMAX*JMAX*PRC_nmax)

  end subroutine mean

  subroutine assert( PHIm, answer, error )
    real(RP), intent(in) :: PHIm(KMAX)
    real(RP), intent(in) :: answer
    real(RP), intent(in) :: error

    real(RP) :: sum

    integer :: k

    sum = 0.D0
    do k = 1, K_zi
       sum = sum + PHIm(k)
    end do
    sum = sum / K_zi

    if ( abs(answer - sum) / answer .lt. error ) then
       write(*,*) 'OK'
    else
       write(*,*) 'NG: ', sum, abs(answer-sum)/answer
    end if

    return
  end subroutine assert

end program check
