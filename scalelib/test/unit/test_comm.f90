module test_comm

  !-----------------------------------------------------------------------------
  use scale_precision
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
    use scale_atmos_grid_cartesC_index
    use scale_prc, only: &
       PRC_myrank
    use scale_prc_cartesC, only: &
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
    use scale_comm_cartesC, only: &
       COMM_regist
    implicit none

    integer, parameter :: KA2 = 10, IA2 = 12, IHALO2 = 4, JA2 = 9, JHALO2 = 3
    integer :: IS2, IE2, JS2, JE2
    integer :: gid
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

    gid = 1
    call test_vars( KA, IA, IS, IE, JA, JS, JE, gid )

    call test_vars8( KA, IA, IS, IE, JA, JS, JE, gid )

    IS2 = IHALO2 + 1
    IE2 = IA2 - IHALO2
    JS2 = JHALO2 + 1
    JE2 = JA2 - JHALO2

    call COMM_regist( KA2, IA2, JA2, IHALO2, JHALO2, gid )

    call test_vars( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

    call test_vars8( KA2, IA2, IS2, IE2, JA2, JS2, JE2, gid )

  end subroutine test_comm_run
!=============================================================================


  subroutine test_vars( KA, IA, IS, IE, JA, JS, JE, gid )
    use scale_prc, only: &
         PRC_myrank
    use scale_comm_cartesC, only: &
         COMM_vars_init, &
         COMM_vars, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in), optional :: gid

    real(RP) :: data3d(KA,IA,JA), data2d(IA,JA)
    integer :: vid

    write(*,*) "Test vars"


    ! MPI
    data2d(:,:) = data_P
    call COMM_vars( data2d, 1, gid=gid )
    call COMM_wait( data2d, 1, .false., gid=gid )
    call check_vars2d( IA, IS, IE, JA, JS, JE, data2d )

    data3d(:,:,:) = data_P
    call COMM_vars( data3d, 1, gid=gid )
    call COMM_wait( data3d, 1, .false., gid=gid )
    call check_vars3d( KA, IA, IS, IE, JA, JS, JE, data3d )

    ! MPI PC
!!$    data2d(:,:) = data_P
!!$    vid = 1
!!$    call COMM_vars_init( 'testdata2d', data2d, vid, gid=gid )
!!$    call COMM_vars( data2d, vid, gid=gid )
!!$    call COMM_wait( data2d, vid, .false., gid=gid )
!!$    call check_vars2d( IA, IS, IE, JA, JS, JE, data2d )

    data3d(:,:,:) = data_P
    vid = 2
    call COMM_vars_init( 'testdata3d', data3d, vid, gid=gid )
    call COMM_vars( data3d, vid, gid=gid )
    call COMM_wait( data3d, vid, .false., gid=gid )
    call check_vars3d( KA, IA, IS, IE, JA, JS, JE, data3d )

  end subroutine test_vars

  subroutine test_vars8( KA, IA, IS, IE, JA, JS, JE, gid )
    use scale_prc, only: &
         PRC_myrank
    use scale_comm_cartesC, only: &
         COMM_vars8_init, &
         COMM_vars8, &
         COMM_wait
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    integer, intent(in) :: gid

    integer, parameter :: nmax = 2
!    integer, parameter :: nmax = 5
    real(RP) :: data3d(KA,IA,JA,nmax), data2d(IA,JA,nmax)
    integer :: vid(nmax)
    integer :: n

    character(len=22) :: title

    write(*,*) "Test vars8"

    write(title(12:22),'(a3,i2.2,a5,i1)') " @ ", PRC_myrank, " gid=", gid

    ! MPI
    data2d(:,:,:) = data_P
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI1(2D)-", n
       call COMM_vars8( data2d(:,:,n), 1, gid=gid )
       call COMM_wait ( data2d(:,:,n), 1, .false., gid=gid )
       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
    end do

    data3d(:,:,:,:) = data_P
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI1(3D)-", n
       call COMM_vars8( data3d(:,:,:,n), 1, gid=gid )
       call COMM_wait ( data3d(:,:,:,n), 1, .false., gid=gid )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,n), title )
    end do

    data2d(:,:,:) = data_P
    do n = 1, nmax
       call COMM_vars8( data2d(:,:,n), n, gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI2(2D)-", n
       call COMM_wait ( data2d(:,:,n), n, .false., gid=gid )
       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
    end do

    data3d(:,:,:,:) = data_P
    do n = 1, nmax
       call COMM_vars8( data3d(:,:,:,n), n, gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a9,i2.2)') "MPI2(3D)-", n
       call COMM_wait ( data3d(:,:,:,n), n, .false., gid=gid )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,n), title )
    end do

    ! MPI PC
!!$    data2d(:,:,:) = data_P
!!$    vid(:) = 1
!!$    do n = 1, nmax
!!$       call COMM_vars8_init( 'PC1-testdata2d', data2d(:,:,n), vid(n), gid=gid )
!!$    end do
!!$    do n = 1, nmax
!!$       write(title(1:11),'(a8,i2.2,a1)') "PC1(2D)-",n," "
!!$       call COMM_vars8( data2d(:,:,n), vid(n), gid=gid )
!!$       call COMM_wait ( data2d(:,:,n), vid(n), .false., gid=gid )
!!$       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
!!$    end do

    data3d(:,:,:,:) = data_P
    vid(:) = 1
    do n = 1, nmax
       call COMM_vars8_init( 'PC1-testdata', data3d(:,:,:,n), vid(n), gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a8,i2.2,a1)') "PC1(3D)-",n," "
       call COMM_vars8( data3d(:,:,:,n), vid(n), gid=gid )
       call COMM_wait ( data3d(:,:,:,n), vid(n), .false., gid=gid )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,n), title )
    end do

!!$    data2d(:,:,:) = data_P
!!$    do n = 1, nmax
!!$       vid(n) = n
!!$       call COMM_vars8_init( 'PC2-testdata2d', data2d(:,:,n), vid(n), gid=gid )
!!$    end do
!!$    do n = 1, nmax
!!$       call COMM_vars8( data2d(:,:,n), vid(n), gid=gid )
!!$    end do
!!$    do n = 1, nmax
!!$       write(title(1:11),'(a8,i02,a1)') "PC2(2D)-",n," "
!!$       call COMM_wait ( data2d(:,:,n), vid(n), .false., gid=gid )
!!$       call check_vars82d( IA, IS, IE, JA, JS, JE, data2d(:,:,n), title )
!!$    end do

    data3d(:,:,:,:) = data_P
    do n = 1, nmax
       vid(n) = n
       call COMM_vars8_init( 'PC2-testdata3d', data3d(:,:,:,n), vid(n), gid=gid )
    end do
    do n = 1, nmax
       call COMM_vars8( data3d(:,:,:,n), vid(n), gid=gid )
    end do
    do n = 1, nmax
       write(title(1:11),'(a8,i02,a1)') "PC2(3D)-",n," "
       call COMM_wait ( data3d(:,:,:,n), vid(n), .false., gid=gid )
       call check_vars83d( KA, IA, IS, IE, JA, JS, JE, data3d(:,:,:,n), title )
    end do

  end subroutine test_vars8

  subroutine check_vars2d( &
       IA, IS, IE, JA, JS, JE, &
       data )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(IA,JA)

    real(RP) :: expect(IA,JA)

    integer :: i, j

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(i,j) = data_N
    end do
    end do
    call AssertEqual("North(2D)", expect(IS:IE,JE+1:JA), data(IS:IE,JE+1:JA))

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(i,j) = data_S
    end do
    end do
    call AssertEqual("South(2D)", expect(IS:IE,1:JS-1), data(IS:IE,1:JS-1))

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(i,j) = data_E
    end do
    end do
    call AssertEqual("East(2D)", expect(IE+1:IA,JS:JE), data(IE+1:IA,JS:JE))

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(i,j) = data_W
    end do
    end do
    call AssertEqual("West(2D)", expect(1:IS-1,JS:JE), data(1:IS-1,JS:JE))

  end subroutine check_vars2d

  subroutine check_vars3d( &
       KA, IA, IS, IE, JA, JS, JE, &
       data )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(KA,IA,JA)

    real(RP) :: expect(KA,IA,JA)

    integer :: i, j

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(:,i,j) = data_N
    end do
    end do
    call AssertEqual("North(3D)", expect(:,IS:IE,JE+1:JA), data(:,IS:IE,JE+1:JA))

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(:,i,j) = data_S
    end do
    end do
    call AssertEqual("South(3D)", expect(:,IS:IE,1:JS-1), data(:,IS:IE,1:JS-1))

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(:,i,j) = data_E
    end do
    end do
    call AssertEqual("East(3D)", expect(:,IE+1:IA,JS:JE), data(:,IE+1:IA,JS:JE))

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(:,i,j) = data_W
    end do
    end do
    call AssertEqual("West(3D)", expect(:,1:IS-1,JS:JE), data(:,1:IS-1,JS:JE))

  end subroutine check_vars3d

  subroutine check_vars82d( &
       IA, IS, IE, JA, JS, JE, &
       data, title )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: IA, IS, IE, JA, JS, JE
    real(RP), intent(in) :: data(IA,JA)
    character(len=*), intent(in) :: title

    real(RP) :: expect(IA,JA)

    integer :: i, j

    do j = JS, JE
    do i = IS, IE
       expect(i,j) = data_P
    end do
    end do

    ! North
    do j = JE+1, JA
    do i = IS, IE
       expect(i,j) = data_N
    end do
    end do

    ! South
    do j = 1, JS-1
    do i = IS, IE
       expect(i,j) = data_S
    end do
    end do

    ! East
    do j = JS, JE
    do i = IE+1, IA
       expect(i,j) = data_E
    end do
    end do

    ! West
    do j = JS, JE
    do i = 1, IS-1
       expect(i,j) = data_W
    end do
    end do

    ! North-East
    do j = JE+1, JA
    do i = IE+1, IA
       expect(i,j) = data_NE
    end do
    end do

    ! North-West
    do j = JE+1, JA
    do i = 1, IS-1
       expect(i,j) = data_NW
    end do
    end do

    ! South-East
    do j = 1, JS-1
    do i = IE+1, IA
       expect(i,j) = data_SE
    end do
    end do

    ! South-West
    do j = 1, JS-1
    do i = 1, IS-1
       expect(i,j) = data_SW
    end do
    end do

    call AssertEqual(title, expect, data)

  end subroutine check_vars82d

  subroutine check_vars83d( &
       KA, IA, IS, IE, JA, JS, JE, &
       data, title )
    use dc_test, only: &
         AssertEqual
    implicit none
    integer, intent(in) :: KA, IA, IS, IE, JA, JS, JE
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

  end subroutine check_vars83d

end module test_comm




