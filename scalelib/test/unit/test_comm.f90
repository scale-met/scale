module test_comm

  !-----------------------------------------------------------------------------
  use scale_precision
  use scale_grid_index
  use scale_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  public :: test_comm_run

  real(RP) :: data_P
  real(RP) :: data_W
  real(RP) :: data_N
  real(RP) :: data_E
  real(RP) :: data_S
  real(RP) :: data_NW
  real(RP) :: data_NE
  real(RP) :: data_SW
  real(RP) :: data_SE

contains

  subroutine test_comm_run
  use scale_process, only: &
     PRC_myrank
  use scale_rm_process, only: &
     PRC_next, &
     PRC_W, &
     PRC_N, &
     PRC_E, &
     PRC_S, &
     PRC_NW, &
     PRC_NE, &
     PRC_SW, &
     PRC_SE, &
     PRC_HAS_W, &
     PRC_HAS_N, &
     PRC_HAS_E, &
     PRC_HAS_S
    implicit none
    !===========================================================================

    data_P = real(PRC_myrank, RP)
    if ( PRC_HAS_W ) then
       data_W = real(PRC_next(PRC_W), RP)
    else
       data_W = data_P
    end if
    if ( PRC_HAS_N ) then
       data_N = real(PRC_next(PRC_N), RP)
    else
       data_N = data_P
    end if
    if ( PRC_HAS_E ) then
       data_E = real(PRC_next(PRC_E), RP)
    else
       data_E = data_P
    end if
    if ( PRC_HAS_S ) then
       data_S = real(PRC_next(PRC_S), RP)
    else
       data_S = data_P
    end if
    if ( PRC_HAS_N ) then
       if ( PRC_HAS_W ) then
          data_NW = real(PRC_next(PRC_NW), RP)
       else
          data_NW = real(PRC_next(PRC_N), RP)
       end if
       if ( PRC_HAS_E ) then
          data_NE = real(PRC_next(PRC_NE), RP)
       else
          data_NE = real(PRC_next(PRC_N), RP)
       end if
    else
       if ( PRC_HAS_W ) then
          data_NW = real(PRC_next(PRC_W), RP)
       else
          data_NW = data_P
       end if
       if ( PRC_HAS_E ) then
          data_NE = real(PRC_next(PRC_E), RP)
       else
          data_NE = data_P
       end if
    end if
    if ( PRC_HAS_S ) then
       if ( PRC_HAS_W ) then
          data_SW = real(PRC_next(PRC_SW), RP)
       else
          data_SW = real(PRC_next(PRC_S), RP)
       end if
       if ( PRC_HAS_E ) then
          data_SE = real(PRC_next(PRC_SE), RP)
       else
          data_SE = real(PRC_next(PRC_S), RP)
       end if
    else
       if ( PRC_HAS_W ) then
          data_SW = real(PRC_next(PRC_W), RP)
       else
          data_SW = data_P
       end if
       if ( PRC_HAS_E ) then
          data_SE = real(PRC_next(PRC_E), RP)
       else
          data_SE = data_P
       end if
    end if

    call test_vars

    call test_vars8


  end subroutine test_comm_run
!=============================================================================


  subroutine test_vars
    use scale_comm, only: &
         COMM_vars_init, &
         COMM_vars, &
         COMM_wait
    implicit none

    real(RP) :: data(KA,IA,JA)
    integer :: vid = 1

    write(*,*) "Test vars"


    ! MPI
    data(:,:,:) = data_P
    call COMM_vars( data, 1 )
    call COMM_wait( data, 1, .false. )
    call check_vars( data )

    ! MPI PC
    data(:,:,:) = data_P
    call COMM_vars_init( 'testdata', data, vid )
    call COMM_vars( data, vid )
    call COMM_wait( data, vid, .false. )
    call check_vars( data )

  end subroutine test_vars

  subroutine test_vars8
    use scale_process, only: &
         PRC_myrank
    use scale_comm, only: &
         COMM_vars8_init, &
         COMM_vars8, &
         COMM_wait
    implicit none

    integer, parameter :: nmax = 2
!    integer, parameter :: nmax = 5
    real(RP) :: data(KA,IA,JA,nmax)
    integer :: vid(nmax)
    integer :: n

    character(len=9) :: title

    write(*,*) "Test vars8"

    write(title(7:9),'(a1,i2)') "@", PRC_myrank

    ! MPI
    data(:,:,:,:) = data_P
    do n = 1, nmax
       write(title(1:6),'(a4,i02)') "MPI-", n
       call COMM_vars8( data(:,:,:,n), 1 )
       call COMM_wait ( data(:,:,:,n), 1, .false. )
       call check_vars8( data(:,:,:,n), title )
    end do

    ! MPI PC
    data(:,:,:,:) = data_P
    vid(:) = 1
    do n = 1, nmax
       call COMM_vars8_init( 'PC1-testdata', data(:,:,:,n), vid(n) )
    end do
    do n = 1, nmax
       write(title(1:6),'(a4,i02)') "PC1-",n
       call COMM_vars8( data(:,:,:,n), vid(n) )
       call COMM_wait ( data(:,:,:,n), vid(n), .false. )
       call check_vars8( data(:,:,:,n), title )
    end do

    data(:,:,:,:) = data_P
    do n = 1, nmax
       vid(n) = n
       call COMM_vars8_init( 'PC2-testdata', data(:,:,:,n), vid(n) )
    end do
    do n = 1, nmax
       call COMM_vars8( data(:,:,:,n), vid(n) )
    end do
    do n = 1, nmax
       write(title(1:6),'(a4,i02)') "PC2-",n
       call COMM_wait ( data(:,:,:,n), vid(n), .false. )
       call check_vars8( data(:,:,:,n), title )
    end do

  end subroutine test_vars8

  subroutine check_vars( data )
    use dc_test, only: &
         AssertEqual
    implicit none

    real(RP), intent(in) :: data(KA,IA,JA)

    real(RP) :: expect(KA,IA,JA)

    integer :: i, j

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(:,i,j) = data_N
    end do
    end do
    call AssertEqual("North", expect(:,IS:IE,JE+1:JA), data(:,IS:IE,JE+1:JA))

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(:,i,j) = data_S
    end do
    end do
    call AssertEqual("South", expect(:,IS:IE,1:JS-1), data(:,IS:IE,1:JS-1))

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(:,i,j) = data_E
    end do
    end do
    call AssertEqual("East", expect(:,IE+1:IA,JS:JE), data(:,IE+1:IA,JS:JE))

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(:,i,j) = data_W
    end do
    end do
    call AssertEqual("West", expect(:,1:IS-1,JS:JE), data(:,1:IS-1,JS:JE))

  end subroutine check_vars

  subroutine check_vars8( data, title )
    use dc_test, only: &
         AssertEqual
    implicit none

    real(RP), intent(in) :: data(KA,IA,JA)
    character(len=*), intent(in) :: title

    real(RP) :: expect(KA,IA,JA)

    integer :: i, j

    do j = JS, JE
    do i = IS, IE
       expect(:,i,j) = data_P
    end do
    end do

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(:,i,j) = data_N
    end do
    end do

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(:,i,j) = data_S
    end do
    end do

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(:,i,j) = data_E
    end do
    end do

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(:,i,j) = data_W
    end do
    end do

    ! North-East
    do j = JE+1, JA
    do i = IE+1, IA
       expect(:,i,j) = data_NE
    end do
    end do

    ! North-West
    do j = JE+1, JA
    do i = 1, IS-1
       expect(:,i,j) = data_NW
    end do
    end do

    ! South-East
    do j = 1, JS-1
    do i = IE+1, IA
       expect(:,i,j) = data_SE
    end do
    end do

    ! South-West
    do j = 1, JS-1
    do i = 1, IS-1
       expect(:,i,j) = data_SW
    end do
    end do

    call AssertEqual(title, expect, data)

  end subroutine check_vars8

end module test_comm




