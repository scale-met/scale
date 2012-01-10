!-------------------------------------------------------------------------------
!> module COMMUNICATION
!!
!! @par Description
!!          MPI module for SCALE3 (Communication Core)
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-10-11 (R.Yoshida) [new]
!! @li      2011-11-11 (H.Yashiro) [mod] Integrate to SCALE3
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
  public :: COMM_setup
  public :: COMM_vars
  public :: COMM_wait
  public :: COMM_stats
  public :: COMM_total
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
  integer, private, save :: MPI_TYPE_WE

  integer, private, save :: COMM_vsize_max = 20

  integer, private, save :: datasize_NS

  integer, private, save :: IREQ_CNT_NS
  integer, private, save :: IREQ_CNT_WE
  integer, private, save :: IREQ_CNT_4C
  integer, private, save :: IREQ_CNT_ALL4

  integer, private, allocatable, save :: ireq_ALL4(:,:)

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  subroutine COMM_setup
    use mod_stdio, only: &
       IO_FID_CONF, &
       IO_FID_LOG,  &
       IO_L
    use mod_process, only: &
       PRC_MPIstop
    use mod_grid, only :    &
       IMAX  => GRID_IMAX,  &
       JMAX  => GRID_JMAX,  &
       IA    => GRID_IA,    &
       KA    => GRID_KA,    &
       IHALO => GRID_IHALO, &
       JHALO => GRID_JHALO, &
       IS    => GRID_IS,    &
       IE    => GRID_IE,    &
       JS    => GRID_JS,    &
       JE    => GRID_JE
    implicit none

    NAMELIST / PARAM_COMM / &
       COMM_vsize_max

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[COMM]/Categ[COMMON]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_COMM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_COMM. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_COMM)

    IREQ_CNT_NS = 2
    IREQ_CNT_WE = 2
    IREQ_CNT_4C = 2 * JHALO
    IREQ_CNT_ALL4 = 2 * IREQ_CNT_NS + 2 * IREQ_CNT_WE

    datasize_NS = IA * KA * JHALO

    call MPI_Type_vector(JMAX, IHALO*KA, IA*KA, MPI_DOUBLE_PRECISION, MPI_TYPE_WE, ierr )
    call MPI_Type_commit(MPI_TYPE_WE, ierr)

    allocate( ireq_ALL4(IREQ_CNT_ALL4,COMM_vsize_max) ); ireq_ALL4(:,:) = 0

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  subroutine COMM_vars( var, vid)
    use mod_stdio, only : &
       IO_FID_LOG, &
       IO_L
    use mod_process, only : &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    use mod_time, only: &
       TIME_rapstart, &
       TIME_rapend
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

    real(8), intent(inout) :: var(:,:,:)
    integer, intent(in)    :: vid

    integer :: ireqc, tag
    integer :: ierr
    !---------------------------------------------------------------------------

    tag = vid * 100
    ireqc = 1

    call TIME_rapstart('COMM_vars')
    !-- From 8-Direction HALO communicate

    ! From S
    call MPI_IRECV( var(:,:,JS-JHALO:JS-1), datasize_NS,          &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag+1, &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From N
    call MPI_IRECV( var(:,:,JE+1:JE+JHALO), datasize_NS,          &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag+2, &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! From W
    call MPI_IRECV( var(1,IE+1,JS), 1,                         &
                    MPI_TYPE_WE, PRC_next(PRC_E), tag+3,       &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! From E
    call MPI_IRECV( var(1,IS-IHALO,JS), 1,                     &
                    MPI_TYPE_WE, PRC_next(PRC_W), tag+4,       &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! To W HALO communicate
    call MPI_ISEND( var(1,IS,JS), 1,                           &
                    MPI_TYPE_WE, PRC_next(PRC_W), tag+3,       &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! To E HALO communicate
    call MPI_ISEND( var(1,IE-IHALO+1,JS), 1,                   &
                    MPI_TYPE_WE, PRC_next(PRC_E), tag+4,       &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr )
    ireqc = ireqc + 1

    ! To S HALO communicate
    call MPI_ISEND( var(:,:,JE-JHALO+1:JE), datasize_NS,          &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_N), tag+1, &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    ! To N HALO communicate
    call MPI_ISEND( var(:,:,JS:JS+JHALO-1), datasize_NS,          &
                    MPI_DOUBLE_PRECISION, PRC_next(PRC_S), tag+2, &
                    MPI_COMM_WORLD, ireq_ALL4(ireqc,vid), ierr    )
    ireqc = ireqc + 1

    call TIME_rapend  ('COMM_vars')

    return
  end subroutine COMM_vars

  !-----------------------------------------------------------------------------
  subroutine COMM_wait( vid )
    use mod_stdio, only : &
       IO_FID_LOG, &
       IO_L
    use mod_time, only: &
       TIME_rapstart, &
       TIME_rapend
    implicit none

    integer, intent(in) :: vid

    integer :: ierr
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM_wait')

    !--- wait packets
    call MPI_WAITALL(IREQ_CNT_ALL4, ireq_ALL4(:,vid), MPI_STATUSES_IGNORE, ierr) ;

    call TIME_rapend  ('COMM_wait')

    return
  end subroutine COMM_wait

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

    logical :: halomask(KA,IA,JA)

    real(8), allocatable :: statval(:,:,:)
    integer, allocatable :: statidx(:,:,:,:)
    real(8), allocatable :: allstatval(:,:)
    integer, allocatable :: allstatidx(:,:,:)
    integer              :: vsize

    integer :: ierr

    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:,:,:,:),4)

    halomask(:,:,:) = .false.
    halomask(KS:KE,IS:IE,JS:JE) = .true.

    allocate( statval(  vsize,2,0:PRC_nmax-1) ); statval(:,:,:)   = CONST_UNDEF8
    allocate( statidx(3,vsize,2,0:PRC_nmax-1) ); statidx(:,:,:,:) = CONST_UNDEF2

    allocate( allstatval(  vsize,2) ); allstatval(:,:)   = CONST_UNDEF8
    allocate( allstatidx(1,vsize,2) ); allstatidx(:,:,:) = CONST_UNDEF2

    if( IO_L ) write(IO_FID_LOG,*)
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

  !-----------------------------------------------------------------------------
  subroutine COMM_total( var, varname )
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
       IA   => GRID_IA,  &
       JA   => GRID_JA,  &
       KA   => GRID_KA,  &
       IS   => GRID_IS,  &
       IE   => GRID_IE,  &
       JS   => GRID_JS,  &
       JE   => GRID_JE,  &
       KS   => GRID_KS,  &
       KE   => GRID_KE,  &
       DXYZ => GRID_DXYZ
    implicit none

    real(8),          intent(inout) :: var(:,:,:,:)
    character(len=*), intent(in)    :: varname(:)

    logical, allocatable :: halomask(:,:,:)

    real(8), allocatable :: statval(:,:)
    real(8), allocatable :: allstatval(:)
    integer              :: vsize

    integer :: ierr

    integer :: v, p
    !---------------------------------------------------------------------------

    vsize = size(var(:,:,:,:),4)

    allocate( halomask(KA,IA,JA) )

    halomask(:,:,:) = .false.
    halomask(KS:KE,IS:IE,JS:JE) = .true.

    allocate( statval   (vsize,0:PRC_nmax-1) ); statval   (:,:) = CONST_UNDEF8
    allocate( allstatval(vsize             ) ); allstatval(:)   = CONST_UNDEF8

    do v = 1, vsize
       statval(v,PRC_myrank) = sum(var(:,:,:,v),mask=halomask)

       ! statistics on each node
!       if( IO_L ) write(IO_FID_LOG,*) '*** [', trim(varname(v)), ']'
!       if( IO_L ) write(IO_FID_LOG,'(1x,A,E17.10,A,I5)') '  SUM = ', statval(v,PRC_myrank), ' at RANK:', PRC_myrank 
    enddo

    ! MPI broadcast
    do p = 0, PRC_nmax-1
       call MPI_Bcast( statval(1,p),         &
                       vsize,                &
                       MPI_DOUBLE_PRECISION, &
                       p,                    &
                       MPI_COMM_WORLD,       &
                       ierr                  )
    enddo

    do v = 1, vsize
       allstatval(v) = sum(statval(v,:))
       if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']',' SUM =', &
                                      allstatval(v) * DXYZ * DXYZ * DXYZ, '[kg * xxx]'
    enddo

    return
  end subroutine COMM_total

end module mod_comm
