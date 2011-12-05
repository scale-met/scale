program binary_to_text
  implicit none

  ! filename length
  character(len=64) :: fn

  ! temporal variables
  integer :: ios
  real :: dat

  ! initialize
  ios=0

  ! Read filename
  if(iargc()/=1) then
    write(*,*) 'Error: required 1 files.'
    stop
  else
    call getarg(1,fn)
  end if

  open(unit=10,file=fn,access='stream')
  open(unit=11,file=trim(fn)//'.txt')

  do while(ios==0)
    read(unit=10,iostat=ios) dat

    ! check status
    if(ios/=0) then
      if(is_iostat_end(ios) .eqv. .true.) then
        write(*,*) 'End of File.'
      else
        write(*,*) 'Error: IOSTAT = ', ios
      end if
      exit ! end do-loop
    end if

    ! write data
    write(unit=11,fmt='(f16.8)') dat

  end do

  close(unit=10)
  close(unit=11)

end program binary_to_text
