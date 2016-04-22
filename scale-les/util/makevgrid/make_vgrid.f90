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
  integer :: num_levs             = 50
  integer :: cnst_nlevs_low       = 10
  integer :: cnst_nlevs_mid       = 10
  integer :: cnst_nlevs_upp       = 10
  integer :: nonlinear_nlevs_mid  = 10
  integer :: nonlinear_nlevs_mid2 = 10
  integer :: nonlinear_nlevs_upp  = 10

  real(4) :: top_height           = 16000.D0
  real(4) :: cnst_thick_low       = 100.D0
  real(4) :: nonlinear_rate_mid   = 8.D0
  real(4) :: nonlinear_rate_mid2  = 8.D0
  real(4) :: cnst_thick_mid       = 500.D0
  real(4) :: nonlinear_rate_upp   = 5.D0
  real(4) :: cnst_thick_upp       = 1000.D0
  real(4) :: nonlinear_rate_top   = 3.D0

  logical :: cnst_middle         = .false.
  logical :: nonlinear_top       = .false.
  logical :: force_skip          = .false.

  real(4),       allocatable :: fz   (:)
  real(4),       allocatable :: dz1  (:)
  real(4),       allocatable :: dz2  (:)
  character(2),  allocatable :: label(:)

  integer :: i, is, ie, integ
  real(4) :: work
  !-----------------------------------------------------------------------------
  
  namelist /make_vgrid/    &
    num_levs,              &
    top_height,            &
    cnst_thick_low,        &
    cnst_nlevs_low,        &
    cnst_thick_mid,        &
    cnst_nlevs_mid,        &
    cnst_thick_upp,        &
    cnst_nlevs_upp,        &
    nonlinear_rate_mid,    &
    nonlinear_nlevs_mid,   &
    nonlinear_rate_mid2,   &
    nonlinear_nlevs_mid2,  &
    nonlinear_nlevs_upp,   &
    nonlinear_rate_upp,    &
    nonlinear_rate_top,    &
    cnst_middle,           &
    nonlinear_top

  open  ( 20, file='vgrid.conf', status='old', delim='apostrophe' )
  read  ( 20, nml=make_vgrid )
  close ( 20 )

  write (*,*) "### SCALE-LES: make setting of vertical grid arrangement"
  write (*,*) "### Prototype 1:"
  write (*,*) "Number of Total Levels:                    ", num_levs
  write (*,*) "Top Height Request (m):                    ", top_height
  write (*,*) "Constant Thickness [Low Level]:            ", cnst_thick_low
  write (*,*) "Number of Constant Levels [Low Level]:     ", cnst_nlevs_low
  write (*,*) "Increasing Rate (%) [Mid Level]:           ", nonlinear_rate_mid
  write (*,*) "Number of Increasing Levels [Mid Level]:   ", nonlinear_nlevs_mid
  write (*,*) "Constant Thickness [Upper Level]:          ", cnst_thick_upp
  write (*,*) "Number of Constant Levels [Upper Level]:   ", cnst_nlevs_upp
  if ( cnst_middle ) then
    write (*,*) "Increasing Rate (%) [Mid Level 2]:         ", nonlinear_rate_mid2
    write (*,*) "Number of Increasing Levels [Mid Level 2]: ", nonlinear_nlevs_mid2
    write (*,*) "Constant Thickness [Mid Level]:            ", cnst_thick_mid
    write (*,*) "Number of Constant Levels [Mid Level]:     ", cnst_nlevs_mid
    write (*,*) "Increasing Rate (%) [Upper Level]:         ", nonlinear_rate_upp
    write (*,*) "Number of Increasing Levels [Upper Level]: ", nonlinear_nlevs_upp
  endif
  if ( nonlinear_top ) then
    write (*,*) "Percent of Increasing (%) [Top Level]:     ", nonlinear_rate_top
  endif
  write (*,*) "Switch for Constant Middle:  ", cnst_middle
  write (*,*) "Switch for Non-linear Top:   ", nonlinear_top
  write (*,*) "------------------------------------------------"

  allocate( fz   (0:num_levs) )
  allocate( dz1  (0:num_levs) )
  allocate( dz2  (0:num_levs) )
  allocate( label(0:num_levs) )

  fz(0) = 0.D0
  integ = 0

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
  integ = integ + cnst_nlevs_low
  do i=is, ie
     fz(i) = fz(i-1) + cnst_thick_low
     label(i) = "LC"
  enddo

  ! nonlinear layer [mid]
  is = integ + 1
  ie = is + nonlinear_nlevs_mid - 1
  integ = integ + nonlinear_nlevs_mid
  do i=is, ie
     fz(i) = fz(i-1) + (fz(i-1)-fz(i-2))*(nonlinear_rate_mid*1.D-2)
     label(i) = "MN"
     if( fz(i) > top_height ) then
        write(*,*) "Over Shoot [MN]: at",i,fz(i)
        stop
     endif
  enddo

  if ( cnst_middle ) then
     ! nonlinear layer [mid:2] for connect cnst layer
     is = integ + 1
     ie = is + nonlinear_nlevs_mid2 - 1
     integ = integ + nonlinear_nlevs_mid2
     do i=is, ie
        fz(i) = fz(i-1) + (fz(i-1)-fz(i-2))*(nonlinear_rate_mid2*1.D-2)
        label(i) = "M2"
        if( fz(i) > top_height ) then
           write(*,*) "Over Shoot [M2]: at",i,fz(i)
           stop
        endif
     enddo

     ! constant layer [mid]
     is = integ + 1
     ie = is + cnst_nlevs_mid - 1
     integ = integ + cnst_nlevs_mid
     do i=is, ie
        fz(i) = fz(i-1) + cnst_thick_mid
        label(i) = "MC"
     enddo

     ! nonlinear layer [upper]
     is = integ + 1
     ie = is + nonlinear_nlevs_upp - 1
     integ = integ + nonlinear_nlevs_upp
     do i=is, ie
        fz(i) = fz(i-1) + (fz(i-1)-fz(i-2))*(nonlinear_rate_upp*1.D-2)
        label(i) = "UN"
        if( fz(i) > top_height ) then
           write(*,*) "Over Shoot [UN]: at",i,fz(i)
           stop
        endif
     enddo
  endif

  ! constant layer [upper]
  is = integ + 1
  if ( nonlinear_top ) then
     ie = is + cnst_nlevs_upp - 1
  else
     ie = num_levs
  endif
  integ = integ + cnst_nlevs_upp
  if ( ie >= num_levs ) then
     ie = num_levs
     force_skip = .true.
  endif
  do i=is, ie
     fz(i) = fz(i-1) + cnst_thick_upp
     label(i) = "UC"
     if( fz(i) > top_height ) then
        write(*,*) "Over Shoot [UC]: at",i,fz(i)
        stop
     endif
  enddo

  ! nonlinear layer [top]: optional
  if ( force_skip ) then
     write(*,*) "WARNING: Constant Upper Layer has reached top already!"
  else
     if ( nonlinear_top ) then
        is = integ + 1
        ie = num_levs
        do i=is, ie
           work = (fz(i-1)-fz(i-2)) - (fz(i-1)-fz(i-2))*(1.0D0 - nonlinear_rate_top*1.D-2)
           fz(i) = fz(i-1) + work
           label(i) = "TN"
           if( fz(i) > top_height*1.2D0 ) then
              write(*,*) "Over Shoot [TN]: at",i,fz(i)
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

  open  (31, file='test.txt', form='formatted', status='replace')
  do i=1, num_levs
     write (31,'(1X,I3,2X,F10.4)') i, fz(i)
  enddo
  close (31)

  stop
  !=============================================================================

end program scaleles_make_vgrid
