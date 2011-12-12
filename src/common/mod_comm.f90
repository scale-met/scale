!-------------------------------------------------------------------------------
!> module Communication
!!
!! @par Description
!!          MPI module for SCALE3
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new] Imported from SCALE-LES ver.2
!!
!<
!-------------------------------------------------------------------------------
module mod_comm
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mpi
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: COMM_vars
  public :: COMM_stats
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
  subroutine COMM_vars( var )
    use mod_process, only : &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    use mod_grid, only : &
       IA    => GRID_IA,    &
       JA    => GRID_JA,    &
       IHALO => GRID_IHALO, &
       JHALO => GRID_JHALO, &
       IS    => GRID_IS,    &
       IE    => GRID_IE,    &
       JS    => GRID_JS,    &
       JE    => GRID_JE
    implicit none

    real(8), intent(inout) :: var(:,:,:,:)

    real(8), allocatable :: sendpack_P2W(:)
    real(8), allocatable :: sendpack_P2N(:)
    real(8), allocatable :: sendpack_P2E(:)
    real(8), allocatable :: sendpack_P2S(:)
    real(8), allocatable :: recvpack_E2P(:)
    real(8), allocatable :: recvpack_S2P(:)
    real(8), allocatable :: recvpack_W2P(:)
    real(8), allocatable :: recvpack_N2P(:)
    integer              :: packsize_NS
    integer              :: packsize_WE
    integer              :: ksize, vsize

    integer, parameter :: tag = 1
    integer            :: status(MPI_STATUS_SIZE)    ! MPI communication parameters
    integer            :: ierr

    integer :: i, j, k, v, n
    !---------------------------------------------------------------------------

    ksize = size(var(:,:,:,:),3)
    vsize = size(var(:,:,:,:),4)

    packsize_NS = IA * JHALO * ksize * vsize
    packsize_WE = JA * IHALO * ksize * vsize

    allocate( sendpack_P2W(packsize_WE) )
    allocate( sendpack_P2N(packsize_NS) )
    allocate( sendpack_P2E(packsize_WE) )
    allocate( sendpack_P2S(packsize_NS) )
    allocate( recvpack_W2P(packsize_WE) )
    allocate( recvpack_N2P(packsize_NS) )
    allocate( recvpack_E2P(packsize_WE) )
    allocate( recvpack_S2P(packsize_NS) )

    !--- packing packets N<->S
    do v = 1, vsize
    do k = 1, ksize

       do j = JS, JS+JHALO-1
       do i = 1, IA
          n = (v-1)  * IA*JHALO*ksize &
            + (k-1)  * IA*JHALO       &
            + (j-JS) * IA             &
            + i

          sendpack_P2N(n) = var(i,j,k,v)
          recvpack_N2P(n) = var(i,JS,k,v) ! for rigid condition
       enddo
       enddo

       do j = JE-JHALO+1, JE
       do i = 1, IA
          n = (v-1)          * IA*JHALO*ksize &
            + (k-1)          * IA*JHALO       &
            + (j-JE+JHALO-1) * IA             &
            + i

          sendpack_P2S(n) = var(i,j,k,v)
          recvpack_S2P(n) = var(i,JE,k,v) ! for rigid condition
       enddo
       enddo

    enddo
    enddo

    ! From Up To Down HALO communicate
    call mpi_sendrecv( sendpack_P2N(:), packsize_NS,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag, &
                       recvpack_S2P(:), packsize_NS,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag, &
                       MPI_COMM_WORLD, status, ierr                )

    ! From Down To Up HALO communicate
    call mpi_sendrecv( sendpack_P2S(:), packsize_NS,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag, &
                       recvpack_N2P(:), packsize_NS,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag, &
                       MPI_COMM_WORLD, status, ierr                )

    !--- unpacking packets N<->S
    do v = 1, vsize
    do k = 1, ksize

       do j = JE+1, JE+JHALO
       do i = 1, IA
          n = (v-1)    * IA*JHALO*ksize &
            + (k-1)    * IA*JHALO       &
            + (j-JE-1) * IA             &
            + i

          var(i,j,k,v) = recvpack_S2P(n)
       enddo
       enddo

       do j = JS-JHALO, JS-1
       do i = 1, IA
          n = (v-1)        * IA*JHALO*ksize &
            + (k-1)        * IA*JHALO       &
            + (j-JS+JHALO) * IA             &
            + i

          var(i,j,k,v) = recvpack_N2P(n)
       enddo
       enddo

    enddo
    enddo

    !--- packing packets W<->E
    do v = 1, vsize
    do k = 1, ksize

       do i = IS, IS+IHALO-1
       do j = 1, JA
          n = (v-1)  * JA*IHALO*ksize &
            + (k-1)  * JA*IHALO       &
            + (i-IS) * JA             &
            + j

          sendpack_P2W(n) = var(i,j,k,v)
          recvpack_W2P(n) = var(IS,j,k,v) ! for rigid condition
       enddo
       enddo

       do i = IE-IHALO+1, IE
       do j = 1, JA
          n = (v-1)          * JA*IHALO*ksize &
            + (k-1)          * JA*IHALO       &
            + (i-IE+IHALO-1) * JA             &
            + j

          sendpack_P2E(n) = var(i,j,k,v)
          recvpack_E2P(n) = var(IE,j,k,v) ! for rigid condition
       enddo
       enddo

    enddo
    enddo

    ! From Right To Left HALO communicate
    call mpi_sendrecv( sendpack_P2W(:), packsize_WE,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_W), tag, &
                       recvpack_E2P(:), packsize_WE,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_E), tag, &
                       MPI_COMM_WORLD, status, ierr                )

    ! From Left To Right HALO communicate
    call mpi_sendrecv( sendpack_P2E(:), packsize_WE,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_E), tag, &
                       recvpack_W2P(:), packsize_WE,               &
                       MPI_DOUBLE_PRECISION, PRC_next(PRC_W), tag, &
                       MPI_COMM_WORLD, status, ierr                )

    do v = 1, vsize
    do k = 1, ksize

       do i = IE+1, IE+IHALO
       do j = 1, JA
          n = (v-1)    * JA*IHALO*ksize &
            + (k-1)    * JA*IHALO       &
            + (i-IE-1) * JA             &
            + j

          var(i,j,k,v) = recvpack_E2P(n)
       enddo
       enddo

       do i = IS-IHALO, IS-1
       do j = 1, JA
          n = (v-1)        * JA*IHALO*ksize &
            + (k-1)        * JA*IHALO       &
            + (i-IS+IHALO) * JA             &
            + j

          var(i,j,k,v) = recvpack_W2P(n)
       enddo
       enddo

    enddo
    enddo

    deallocate( sendpack_P2W )
    deallocate( sendpack_P2N )
    deallocate( sendpack_P2E )
    deallocate( sendpack_P2S )
    deallocate( recvpack_E2P )
    deallocate( recvpack_S2P )
    deallocate( recvpack_W2P )
    deallocate( recvpack_N2P )

    return
  end subroutine COMM_vars

  !-----------------------------------------------------------------------------
  subroutine COMM_stats( var, varname )
    use mod_stdio, only : &
       IO_FID_LOG, &
       IO_L
    use mod_process, only : &
       PRC_nmax,   &
       PRC_myrank
    use mod_const, only : &
       CONST_UNDEF8, &
       CONST_UNDEF2
    use mod_grid, only : &
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE, &
       KS => GRID_KS, &
       KE => GRID_KE
    implicit none

    real(8),          intent(inout) :: var(:,:,:,:)
    character(len=*), intent(in)    :: varname(:)

    logical, allocatable :: halomask(:,:,:)

    real(8), allocatable :: statval(:,:,:)
    integer, allocatable :: statidx(:,:,:,:)
    real(8), allocatable :: allstatval(:,:)
    integer, allocatable :: allstatidx(:,:,:)
    integer              :: vsize

    integer :: ierr

    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:,:,:,:),4)

    allocate( halomask(IA,JA,KA) )

    halomask(:,:,:) = .false.
    halomask(IS:IE,JS:JE,KS:KE) = .true.

    allocate( statval(  vsize,2,0:PRC_nmax-1) ); statval(:,:,:)   = CONST_UNDEF8
    allocate( statidx(3,vsize,2,0:PRC_nmax-1) ); statidx(:,:,:,:) = CONST_UNDEF2

    allocate( allstatval(  vsize,2) ); allstatval(:,:)   = CONST_UNDEF8
    allocate( allstatidx(1,vsize,2) ); allstatidx(:,:,:) = CONST_UNDEF2

    if( IO_L ) write(IO_FID_LOG,*) '*** Variable Statistics ***'
    do v = 1, vsize
       statval(  v,1,PRC_myrank) = maxval(var(:,:,:,v),mask=halomask)
       statval(  v,2,PRC_myrank) = minval(var(:,:,:,v),mask=halomask)
       statidx(:,v,1,PRC_myrank) = maxloc(var(:,:,:,v),mask=halomask)
       statidx(:,v,2,PRC_myrank) = minloc(var(:,:,:,v),mask=halomask)

! statistics on each node
!       if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,3(I5,A))') '*** MAX = ', &
!                                             statval(  v,1,PRC_myrank),'(', &
!                                             statidx(1,v,1,PRC_myrank),',', &
!                                             statidx(2,v,1,PRC_myrank),',', &
!                                             statidx(3,v,1,PRC_myrank),')'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,3(I5,A))') '*** MIN = ', &
!                                             statval(  v,2,PRC_myrank),'(', &
!                                             statidx(1,v,2,PRC_myrank),',', &
!                                             statidx(2,v,2,PRC_myrank),',', &
!                                             statidx(3,v,2,PRC_myrank),')'
   enddo

    ! MPI broadcast
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,1,p),       &
                       vsize*2,              &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
       call MPI_Bcast( statidx(1,1,1,p),     &
                       3*vsize*2,            &
                       MPI_INTEGER,          &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
    enddo

    do v = 1, vsize
       allstatval(v,1)   = maxval(statval(v,1,:))
       allstatval(v,2)   = minval(statval(v,2,:))
       allstatidx(:,v,1) = maxloc(statval(v,1,:))-1
       allstatidx(:,v,2) = minloc(statval(v,2,:))-1
       if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,4(I5,A))') '  MAX =', &
                                                    allstatval(  v,1), '(', &
                                                    allstatidx(1,v,1), ',', &
                                      statidx(1,v,1,allstatidx(1,v,1)),',', &
                                      statidx(2,v,1,allstatidx(1,v,1)),',', &
                                      statidx(3,v,1,allstatidx(1,v,1)),')'
       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,4(I5,A))') '  MIN =', &
                                                    allstatval(  v,2), '(', &
                                                    allstatidx(1,v,2), ',', &
                                      statidx(1,v,2,allstatidx(1,v,2)),',', &
                                      statidx(2,v,2,allstatidx(1,v,2)),',', &
                                      statidx(3,v,2,allstatidx(1,v,2)),')'
    enddo

    return
  end subroutine COMM_stats

end module mod_comm
