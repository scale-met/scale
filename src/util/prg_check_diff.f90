program check_diff
  implicit none

  ! criterion of difference
  real, parameter :: crt = 1.e-6

  ! check variable
  real :: diff

  ! filename length
  character(len=64) :: fn_a, fn_b

  ! temporal variables
  integer :: ios_a, ios_b
  real :: dat_a, dat_b

  integer :: iargc

  ! initialize
  ios_a=0
  ios_b=0

  ! Read filename
  if(iargc()/=2) then
    write(*,*) 'Error: required two files.'
    call exit(1)
  else
    call getarg(1,fn_a)
    call getarg(2,fn_b)
  end if

  open(unit=10,file=fn_a,access='stream')
  open(unit=11,file=fn_b,access='stream')

  ! loop check difference
  do while(ios_a==0 .and. ios_b==0)
    read(unit=10,iostat=ios_a) dat_a
    read(unit=11,iostat=ios_b) dat_b

    ! check status
    if(ios_a/=0 .or. ios_b/=0) then
      if(is_iostat_end(ios_a) .eqv. .true. .and. &
         is_iostat_end(ios_b) .eqv. .true. ) then
        write(*,*) 'Check complete!'
        call exit(0)
      else
        write(*,*) 'Error: IOSTAT = ', ios_a
        write(*,*) 'Error: IOSTAT = ', ios_b
        call exit(1)
      end if
    end if

    ! find difference
    diff = abs(dat_a - dat_b)
    if(diff > crt) then
      write(*,*), 'Found difference!'
      write(*,fmt='(a8,f12.8)'), 'File A: ', dat_a
      write(*,fmt='(a8,f12.8)'), 'File B: ', dat_b
      call exit(1)
    end if

  end do

  close(unit=10)
  close(unit=20)

end program check_diff
