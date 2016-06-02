module mod_comm
  !
  !-----------------------------------------------------------
  ! 2014/11/26: Original (Ryuji Yoshida)
  !
  !-----------------------------------------------------------

  use mpi !external
  use mod_vars

  implicit none

  !--- Public
  !-----------------------------------------------------------
  public :: comm_spawn
  public :: comm_get_parent
  public :: comm_makeyp
  public :: comm_send
  public :: comm_recv
  public :: comm_wait
  public :: comm_barrier
  public :: comm_free

  !--- Private
  !-----------------------------------------------------------
  private :: comm_setup
  private :: comm_abort

  integer, private :: nprocs_parent(2)
  integer, private :: nprocs_child(2)

  integer, private :: ypall_p
  integer, private :: ypall_c

  integer, private :: nreq_p
  integer, private :: nreq_c

  integer, private,  allocatable :: yplist_p(:,:)
  integer, private,  allocatable :: yplist_c(:,:)

  character(len=20) :: argv(2)
 
  !-----------------------------------------------------------
  contains
  !-----------------------------------------------------------


  !-----------------------------------------------------------
  subroutine comm_spawn ( DOMAIN_NUM, LAUNCH_PRC )
    integer, intent(in) :: DOMAIN_NUM
    integer, intent(in) :: LAUNCH_PRC

    integer :: ierr
    integer :: errcodes(1:LAUNCH_PRC)
    character(2) :: dom_num
    ! -----

    write(dom_num,'(I2.2)') DOMAIN_NUM+1
    argv(1) = 'run.d'//dom_num//'.conf'
    argv(2) = ''

    write ( LFID, * ) "COMM/SPAWN: LAUNCH_PRC = ", LAUNCH_PRC
    call MPI_COMM_SPAWN( CMD, argv, LAUNCH_PRC, MPI_INFO_NULL, &
                         MASTER, MPI_COMM_WORLD, INTERCOMM_CHILD, errcodes, ierr )

    write ( LFID, * ) "COMM/SPAWN: INTERCOMM_CHILD = ", INTERCOMM_CHILD
   return
  end subroutine comm_spawn
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_get_parent ()
    integer :: ierr
    ! -----

    call MPI_COMM_GET_PARENT( INTERCOMM_PARENT, ierr )
    call MPI_COMM_REMOTE_SIZE( INTERCOMM_PARENT, PARENT_SIZE, ierr)
    write ( LFID, * ) "COMM/SPAWN: PARENT PROCESS SIZE = ", PARENT_SIZE
    write ( LFID, * ) "COMM/SPAWN: INTERCOMM_PARENT = ", INTERCOMM_PARENT

   return
  end subroutine comm_get_parent
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_setup ( handle )
    integer, intent(in) :: handle
    ! -----

    if ( handle == parent ) then !--- parent
       allocate ( send_buffer( NX, NY, NZ )                )
       allocate ( yplist_p( nprocs_parent(parent), ypmax ) )
       allocate ( ireq_p( ypmax )                            )

    elseif ( handle == child ) then !--- child
       allocate ( recv_buffer ( NX, NY, NZ, ypmax )        )
       allocate ( yplist_c( nprocs_child(child), ypmax )   )
       allocate ( ireq_c( ypmax )                            )

    else
       call comm_abort
    endif

   return
  end subroutine comm_setup
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_makeyp ( nprocs_p, nprocs_c, handle )
    integer, intent(in)    :: nprocs_p
    integer, intent(in)    :: nprocs_c
    integer, intent(in)    :: handle

    integer :: ileng, ierr
    integer :: i, j, k
    logical :: flag
    ! -----

    if ( handle == parent ) then !--- parent
       nprocs_parent(parent) = nprocs_p
       nprocs_child(parent)  = nprocs_c

       call comm_setup ( handle )

       yplist_p(:,:) = -1

       if ( master == myproc ) then
          k = 1
          i = 1
          flag = .false.
          do j = 1, nprocs_child(parent)
             if (flag) then
                k = k + 1
                flag = .false.
             endif

             yplist_p(i,k) = j - 1
             i = i + 1
             if ( i > nprocs_parent(parent) ) then
                i = 1
                flag = .true.
             endif
          enddo

          ypall_p = k
       endif

       ileng = nprocs_parent(parent) * ypmax
       call MPI_BCAST( yplist_p, ileng, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
       call MPI_BCAST( ypall_p,  1,     MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

       write ( LFID, * ) "COMM/MAKEYP[P]: ypall_p", ypall_p
       i = myproc + 1
       do k = 1, ypall_p
          write ( LFID, '(1X,A,I3,A,I3)' ) "COMM/MAKEYP[P]: SEND(me) ", &
                                        myproc, " --> RECV ", yplist_p(i,k)
       enddo

    elseif ( handle == child ) then !--- child
       nprocs_parent(child) = nprocs_p
       nprocs_child(child)  = nprocs_c

       call comm_setup ( handle )

       yplist_c(:,:) = -1

       if ( master == myproc ) then
          k = 1
          i = 1
          do j = 1, nprocs_child(child)
             yplist_c(j,k) = i - 1
             i = i + 1
             if ( i > nprocs_parent(child) ) i = 1
          enddo

          ypall_c = k
       endif

       ileng = nprocs_child(child) * ypmax
       call MPI_BCAST( yplist_c, ileng, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
       call MPI_BCAST( ypall_c,  1,     MPI_INTEGER, master, MPI_COMM_WORLD, ierr )

       write ( LFID, * ) "COMM/MAKEYP[C]: ypall_c", ypall_c
       i = myproc + 1
       do k = 1, ypall_c
          write ( LFID, '(1X,A,I3,A,I3)' ) "COMM/MAKEYP[C]: SEND ", &
                                        yplist_c(i,k), " --> RECV(me) ", myproc
       enddo
    else
       call comm_abort
    endif

   return
  end subroutine comm_makeyp
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_send ( var, handle )
    real(SP), intent(in) :: var(:,:,:)
    integer,  intent(in) :: handle

    integer :: ileng, ierr
    integer :: i
    ! -----

    if ( handle == parent ) then
       nreq_p = 0

       do i = 1, ypall_p
          dist  = yplist_p( myproc+1, i )
          if ( dist >= 0 ) then
             send_buffer(:,:,:) = var(:,:,:)
             ileng  = nx * ny * nz
             tag    = dist + 1
             nreq_p = nreq_p + 1
             call mpi_isend( send_buffer, ileng, MPI_REAL, &
                             dist, tag, INTERCOMM_CHILD, ireq_p(nreq_p), ierr )
             write ( LFID, '(1X,A,I3,A,I3)' ) &
                   "COMM/SEND: to pe =", dist, " with a tag of ", tag
          endif
       enddo
    else
       call comm_abort
    endif

   return
  end subroutine comm_send
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_recv ( handle )
    integer, intent(in) :: handle

    integer :: ileng, ierr
    integer :: i
    ! -----

    if ( handle == child ) then
       nreq_c = 0

       do i = 1, ypall_c
          source = yplist_c( myproc+1, i )
          if ( source >= 0 ) then
             recv_buffer(:,:,:,i) = 0.0D0
             ileng = nx * ny * nz
             tag = myproc + 1
             nreq_c = nreq_c + 1
             call mpi_irecv( recv_buffer(:,:,:,i), ileng, MPI_REAL, &
                             source, tag, INTERCOMM_PARENT, ireq_c(nreq_c), ierr )
             write ( LFID, '(1X,A,I3,A,I3)' ) &
                   "COMM/RECV: from pe =", source, " with a tag of ", tag
          endif
       enddo
    else
       call comm_abort
    endif

   return
  end subroutine comm_recv
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_wait ( var, handle )
    real(SP), intent(out) :: var(:,:,:)
    integer, intent(in) :: handle

    integer :: istatus_p(MPI_STATUS_SIZE,ypall_p)
    integer :: istatus_c(MPI_STATUS_SIZE,ypall_c)
    integer :: ierr

    integer    :: i
    integer(8) :: num
    logical    :: flag
    ! -----

    num  = 0
    flag = .false.

    if ( handle == parent ) then
       do while ( .NOT. flag )
          call mpi_testall( nreq_p, ireq_p, flag, istatus_p, ierr )
          num = num + 1

          if ( num > 100000000 ) then
             write ( LFID, * ) "ERROR: Dead Lock [COMM/WAIT: P]"
             call comm_abort
          endif
       enddo

       var(:,:,:) = UNDEF ! dummy
    elseif ( handle == child ) then
       do while ( .NOT. flag )
          call mpi_testall( nreq_c, ireq_c, flag, istatus_c, ierr )
          num = num + 1

          if ( num > 100000000 ) then
             write ( LFID, * ) "ERROR: Dead Lock [COMM/WAIT: C]"
             call comm_abort
          endif
       enddo

       ! copy from buffer
       do i = 1, ypall_c
          var(:,:,:) = recv_buffer(:,:,:,i)
       enddo
    else
       call comm_abort
    endif

   return
  end subroutine comm_wait
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_barrier ( handle )
    integer, intent(in) :: handle

    integer :: ierr
    ! -----

    if ( handle == parent ) then
       write ( LFID, * ) "COMM/BARRIER: as a PARENT >> PAUSE"
       call MPI_BARRIER(INTERCOMM_CHILD, ierr)
       write ( LFID, * ) " >> RESUME"
    elseif ( handle == child ) then
       write ( LFID, * ) "COMM/BARRIER: as a CHILD >> PAUSE"
       call MPI_BARRIER(INTERCOMM_PARENT, ierr)
       write ( LFID, * ) " >> RESUME"
    else
       call comm_abort
    endif

   return
  end subroutine comm_barrier
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_free ( handle )
    integer, intent(in) :: handle

    integer :: ierr
    ! -----

    if ( handle == parent ) then
       call MPI_COMM_FREE(INTERCOMM_CHILD, ierr)
       write ( LFID, * ) "COMM/FREE: disconnected with child"
    elseif ( handle == child ) then
       call MPI_COMM_FREE(INTERCOMM_PARENT, ierr)
       write ( LFID, * ) "COMM/FREE: disconnected with parent"
    else
       call comm_abort
    endif

   return
  end subroutine comm_free
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine comm_abort ()
    integer :: ierr
    ! -----

    write ( LFID, * ) ""
    write ( LFID, * ) "COMM/ABORT: ERROR"
    close ( LFID )

    if ( INTERCOMM_CHILD /= 0 ) then
       call MPI_ABORT(INTERCOMM_CHILD, 1, ierr)
    endif

    if ( INTERCOMM_PARENT /= 0 ) then
       call MPI_ABORT(INTERCOMM_PARENT, 1, ierr)
    endif

    call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)

   stop
  end subroutine comm_abort
  !-----------------------------------------------------------


  !-----------------------------------------------------------
end module mod_comm
