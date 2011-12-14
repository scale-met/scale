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
  integer, private, save :: COMM_vsize_max = 20

  integer, private, save :: packsize_NS_max
  integer, private, save :: packsize_WE_max

  real(8), private, allocatable, save :: sendpack_P2W(:)
  real(8), private, allocatable, save :: sendpack_P2N(:)
  real(8), private, allocatable, save :: sendpack_P2E(:)
  real(8), private, allocatable, save :: sendpack_P2S(:)
  real(8), private, allocatable, save :: recvpack_E2P(:)
  real(8), private, allocatable, save :: recvpack_S2P(:)
  real(8), private, allocatable, save :: recvpack_W2P(:)
  real(8), private, allocatable, save :: recvpack_N2P(:)

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
    use mod_grid, only : &
       IA    => GRID_IA,    &
       JA    => GRID_JA,    &
       KA    => GRID_KA,    &
       IHALO => GRID_IHALO, &
       JHALO => GRID_JHALO
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

    packsize_NS_max = IA * JHALO * KA * COMM_vsize_max
    packsize_WE_max = JA * IHALO * KA * COMM_vsize_max

    allocate( sendpack_P2W(packsize_WE_max) )
    allocate( sendpack_P2N(packsize_NS_max) )
    allocate( sendpack_P2E(packsize_WE_max) )
    allocate( sendpack_P2S(packsize_NS_max) )
    allocate( recvpack_W2P(packsize_WE_max) )
    allocate( recvpack_N2P(packsize_NS_max) )
    allocate( recvpack_E2P(packsize_WE_max) )
    allocate( recvpack_S2P(packsize_NS_max) )

    return
  end subroutine COMM_setup

  !-----------------------------------------------------------------------------
  subroutine COMM_vars( var )
    use mod_process, only : &
       PRC_next, &
       PRC_W,    &
       PRC_N,    &
       PRC_E,    &
       PRC_S
    use mod_time, only: &
       TIME_rapstart,                    &
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

    real(8), intent(inout) :: var(:,:,:,:)

    integer              :: packsize_NS
    integer              :: packsize_WE
    integer              :: ksize, vsize

    integer, parameter :: tag = 1
    integer            :: status(MPI_STATUS_SIZE)    ! MPI communication parameters
    integer            :: ierr

    integer :: i, j, k, v, n
    !---------------------------------------------------------------------------

    call TIME_rapstart('COMM_vars')

    ksize = size(var(:,:,:,:),3)
    vsize = size(var(:,:,:,:),4)

    packsize_NS = IA * JHALO * ksize * vsize
    packsize_WE = JA * IHALO * ksize * vsize

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

    call TIME_rapstart('MPI')
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
    call TIME_rapend  ('MPI')

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

    call TIME_rapstart('MPI')
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
    call TIME_rapend  ('MPI')

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

    call TIME_rapend  ('COMM_vars')

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
       IA => GRID_IA, &
       JA => GRID_JA, &
       KA => GRID_KA, &
       IS => GRID_IS, &
       IE => GRID_IE, &
       JS => GRID_JS, &
       JE => GRID_JE, &
       KS => GRID_KS, &
       KE => GRID_KE, &
       DX => GRID_DX, &
       DY => GRID_DY, &
       DZ => GRID_DZ
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

    allocate( halomask(IA,JA,KA) )

    halomask(:,:,:) = .false.
    halomask(IS:IE,JS:JE,KS:KE) = .true.

    allocate( statval(  vsize,0:PRC_nmax-1) ); statval(:,:)   = CONST_UNDEF8

    allocate( allstatval(  vsize) ); allstatval(:)   = CONST_UNDEF8

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
       if( IO_L ) write(IO_FID_LOG,*) '[', trim(varname(v)), ']',' SUM =', allstatval(v)*DX*DY*DZ*1.D-9
    enddo

    return
  end subroutine COMM_total

end module mod_comm
