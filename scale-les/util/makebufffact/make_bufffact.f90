!-------------------------------------------------------------------------------
!> Program SCALE-LES Buffer factor Maker
!!
!! @par Description
!!          This program is make of buffer factor
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaleles_make_bufffact
  !-----------------------------------------------------------------------------
  !
  implicit none

  ! namelists
  integer :: MAXIMUM_ITERATION  = 100
  integer :: NUM_STRETCHED_GRID = 10

  real(8) :: DEFAULT_BUFFFACT = 1.1D0
  real(8) :: GRID_DISTANCE    = 100.0D0
  real(8) :: TOTAL_LENGTH     = 100.0D+3
  real(8) :: LIMIT_ITERATION  = epsilon(0.0D0)

  namelist / MAKE_BUFFFACT / &
    MAXIMUM_ITERATION,  &
    DEFAULT_BUFFFACT,   &
    NUM_STRETCHED_GRID, &
    GRID_DISTANCE,      &
    TOTAL_LENGTH,       &
    LIMIT_ITERATION

  ! works
  integer :: n
  integer :: ierr

  integer :: itr
  integer :: num

  real(8) :: bf
  real(8) :: dx
  real(8) :: L
  real(8) :: EPS

  real(8) :: dbf
  !-----------------------------------------------------------------------------
  
  ! read namelist
  open(unit=10, file='bufffact.conf', &
       status='old', delim='apostrophe')
  read(unit=10, nml=MAKE_BUFFFACT, iostat=ierr)

  if( ierr < 0 ) then !--- missing
    write(*,*) '*** Not found namelist. Default used.'
  elseif( ierr > 0 ) then !--- fatal error
    write(*,*) 'xxx Not appropriate names in namelist MAKE_BUFFFACT. Check!'
    stop
  endif

  close(unit=10)

  ! initialize
  itr = MAXIMUM_ITERATION
  num = NUM_STRETCHED_GRID

  bf  = DEFAULT_BUFFFACT
  dx  = GRID_DISTANCE
  L   = TOTAL_LENGTH
  EPS = LIMIT_ITERATION

  ! output setting parameter
  write (*,*)
  write (*,*) "-------------------------------------------------------------------------------"
  write (*,*) "### SCALE-LES: make setting of buffer factor ###"
  write (*,*) "maximum iteration number                 ", itr
  write (*,*) "number of stretched grid                 ", num
  write (*,*) "default buffer factor                    ", bf
  write (*,*) "non-stretched grid distance [m]          ", dx
  write (*,*) "total length of stretched region [m]     ", L
  write (*,*) "limit of iteration                       ", EPS
  write (*,*) "-------------------------------------------------------------------------------"
  write (*,*)

  write (*,*) "-------------------------------------------------------------------------------"
  write (*,*) "  # of iteration            buffer factor    diff. of buffer factor            "
  write (*,*) "-------------------------------------------------------------------------------"

  ! Newton-Raphson method
  do n = 1, itr
    dbf = ( bf**dble(num+1)*dx - bf*(L+dx) + L ) / ( bf**dble(num)*dble(num+1)*dx - (L+dx) )
    bf  = bf - dbf

    write (*,'(12x,i5,1x,f24.16,2x,e24.16)') n, bf, dbf

    if( dbf < EPS ) then
      exit ! stop iteration
    end if
  end do

  write (*,*) "-------------------------------------------------------------------------------"

  open(unit=20, file='result.txt', &
       form='formatted', status='replace' )
  write(unit=20, fmt='(" BUFFFACT  = ",f20.16,"D0,")')  bf
  close(unit=20)

end program scaleles_make_bufffact
