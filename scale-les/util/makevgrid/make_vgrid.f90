!-------------------------------------------------------------------------------
!> Program SCALE-LES Vertical Grid Maker
!!
!! @par Description
!!          This program is make of vertical grid arrangement
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scaleles_make_vgrid
  !-----------------------------------------------------------------------------
  !
  implicit none
  integer :: num_levs            = 30
  integer :: cnst_nlevs_low      = 10
  integer :: nonlinear_nlevs_mid = 10
  integer :: cnst_nlevs_upp      = 10

  real(4) :: top_height          = 16000.D0
  real(4) :: cnst_thick_low      = 100.D0
  real(4) :: nonlinear_rate_mid  = 8.D0
  real(4) :: cnst_thick_upp      = 1000.D0
  real(4) :: nonlinear_rate_top  = 3.D0

  logical :: nonlinear_top       = .false.
  logical :: force_skip          = .false.

  real(4),       allocatable :: fz   (:)
  real(4),       allocatable :: dz1  (:)
  real(4),       allocatable :: dz2  (:)
  character(1), allocatable :: label(:)

  integer :: i, is, ie
  real(4) :: work
  !-----------------------------------------------------------------------------
  
  namelist /make_vgrid/    &
    num_levs,              &
    top_height,            &
    cnst_thick_low,        &
    cnst_nlevs_low,        &
    nonlinear_rate_mid,    &
    nonlinear_nlevs_mid,   &
    cnst_thick_upp,        &
    cnst_nlevs_upp,        &
    nonlinear_top,         &
    nonlinear_rate_top

  open  ( 20, file='vgrid.conf', status='old', delim='apostrophe' )
  read  ( 20, nml=make_vgrid )
  close ( 20 )

  write (*,*) "### SCALE-LES: make setting of vertical grid arrangement"
  write (*,*) "### Prototype 1:"
  write (*,*) "Number of Total Levels:                  ", num_levs
  write (*,*) "Top Height Request (m):                  ", top_height
  write (*,*) "Constant Thickness [Low Level]:          ", cnst_thick_low
  write (*,*) "Number of Constant Levels [Low Level]:   ", cnst_nlevs_low
  write (*,*) "Increasing Rate (%) [Mid Level]:         ", nonlinear_rate_mid
  write (*,*) "Number of Increasing Levels [Mid Level]: ", nonlinear_nlevs_mid
  write (*,*) "Constant Thickness [Upper Level]:        ", cnst_thick_upp
  write (*,*) "Number of Constant Levels [Upper Level]: ", cnst_nlevs_upp
  write (*,*) "Top Level Non-linear Configuration:      ", nonlinear_top
  write (*,*) "Percent of Increasing (%) [Top Level]:   ", nonlinear_rate_top
  write (*,*) "------------------------------------------------"

  allocate( fz   (0:num_levs) )
  allocate( dz1  (0:num_levs) )
  allocate( dz2  (0:num_levs) )
  allocate( label(0:num_levs) )

  fz(0) = 0.D0

  ! setting check
  work = cnst_thick_low * float(cnst_nlevs_low)
  if( work >= top_height )then
     write(*,*) '"Constant Thickness" or "Number of Constant Levels" is not appropriate'
     stop
  endif
  work = cnst_nlevs_low + nonlinear_nlevs_mid + cnst_nlevs_upp
  if( work > num_levs )then
     write(*,*) "Total of Requested Layers is over the num_levs"
     stop
  endif

  ! constant layer [low]
  is = 1
  ie = is + cnst_nlevs_low - 1
  do i=is, ie
     fz(i) = fz(i-1) + cnst_thick_low
     label(i) = "L"
  enddo

  ! nonlinear layer [mid]
  is = cnst_nlevs_low + 1
  ie = is + nonlinear_nlevs_mid - 1
  do i=is, ie
     fz(i) = fz(i-1) + fz(i-1)*(nonlinear_rate_mid*1.D-2)
     label(i) = "M"
     if( fz(i) > top_height ) then
        write(*,*) "Over Shoot [M]: at",i,fz(i)
        stop
     endif
  enddo

  ! constant layer [upper]
  is = cnst_nlevs_low + nonlinear_nlevs_mid + 1
  ie = is + cnst_nlevs_upp - 1
  if ( ie >= num_levs ) then
     ie = num_levs
     force_skip = .true.
  endif
  do i=is, ie
     fz(i) = fz(i-1) + cnst_thick_upp
     label(i) = "U"
     if( fz(i) > top_height ) then
        write(*,*) "Over Shoot [U]: at",i,fz(i)
        stop
     endif
  enddo

  ! nonlinear layer [top]: optional
  if ( force_skip ) then
     write(*,*) "WARNING: Constant Upper Layer has reached top already!"
  else
     if ( nonlinear_top ) then
        is = cnst_nlevs_low + nonlinear_nlevs_mid + cnst_nlevs_upp + 1
        ie = num_levs
        do i=is, ie
           fz(i) = fz(i-1) + cnst_thick_upp - fz(i-1)*(nonlinear_rate_top*1.D-2)
           label(i) = "T"
           if( fz(i) > top_height*1.2D0 ) then
              write(*,*) "Over Shoot [T]: at",i,fz(i)
              stop
           endif
        enddo
     endif
  endif

  dz1(0) = 0.D0
  do i=1, num_levs
     dz1(i) = fz(i) - fz(i-1)
  enddo
  dz2(0) = 0.D0
  dz2(1) = 0.D0
  do i=2, num_levs
     dz2(i) = dz1(i) - dz1(i-1)
  enddo

  write (*,*) "---------------------------------------------------"
  write (*,*) "  z        fz         dz^1         dz^2  process"
  write (*,*) "---------------------------------------------------"
  do i=0, num_levs
     write (*,'(1X,I3,1X,"|",1X,F10.4,2X,F10.4,2X,F10.4,4X,A)') i, fz(i), dz1(i), dz2(i), label(i)
  enddo
  write (*,*) "---------------------------------------------------"

  open  ( 30, file='vgrid_arrangement.txt', form='formatted', status='replace' )
  write ( 30,'(" FZ(:) = ",5(1X,F10.4,"D0,"))' ) fz(1), fz(2), fz(3), fz(4), fz(5)
  do i=6, num_levs, 5
     write ( 30,'("         ",5(1X,F10.4,"D0,"))' ) fz(i), fz(i+1), fz(i+2), fz(i+3), fz(i+4)
  enddo
  close ( 30 )

  !open  (31, file='test.txt', form='formatted', status='replace')
  !do i=1, num_levs
  !   write (31,'(1X,I3,2X,F10.4)') i, fz(i)
  !enddo
  !close (31)

  stop
  !=============================================================================

end program scaleles_make_vgrid
