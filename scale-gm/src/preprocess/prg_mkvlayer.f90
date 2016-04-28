!-------------------------------------------------------------------------------
!> Program mkrawgrid
!!
!! @par Description
!!          Making vertical grid systems based on lorenz coordinate
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program prg_mkvlayer
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  !=============================================================================

  integer, parameter :: kdum=1
  integer, parameter :: fid=11
  integer :: num_of_layer = 10
  real(RP) :: ztop = 1.E4_RP
  character(256) :: outfname = 'outfile'
  character(256) :: infname = 'infile'
  character(16) :: layer_type = 'POWER'
  real(RP) :: fact = 1.0_RP

  namelist / mkvlayer_cnf / &
       num_of_layer,        & !--- number of layers
       ztop,                & !--- height of model top
       outfname,            & !--- output file name
       layer_type,          & !--- type of layer
       fact,                & !--- factor  if layer_type='POWER'
       infname                !--- input file name if layer_type='GIVEN'

  real(RP),allocatable :: z_c(:)
  real(RP),allocatable :: z_h(:)

  integer :: kmin
  integer :: kmax
  integer :: kall

  open(fid,file='mkvlayer.cnf',status='old',form='formatted')
  read(fid,nml=mkvlayer_cnf)
  close(fid)


  select case(trim(layer_type))
  case('POWER')
     call mk_layer_powerfunc(fact)
  case('GIVEN')
     call mk_layer_given(infname)
  end select
  !
  call output_layer(outfname)

  !=============================================================================
contains
  !-----------------------------------------------------------------------------
  subroutine mk_layer_powerfunc( fact )
    implicit none
    real(RP) :: a
    integer :: k
    real(RP),intent(in) :: fact

    kmin=kdum+1
    kmax=kdum+num_of_layer
    kall=kdum+num_of_layer+kdum
    allocate(z_c(kall))
    allocate(z_h(kall))
    !
    a = ztop/(dble(num_of_layer)**fact)
    !
    do k = kmin, kmax+1
       z_h(k) = a * (dble(k-kmin)**fact)
    enddo
    !
    z_h(kmin-1) = z_h(kmin) - ( z_h(kmin+1) - z_h(kmin) )
    !
    do k= kmin-1, kmax
       z_c(k) = z_h(k) + ( z_h(k+1) - z_h(k) )*0.5_RP
    enddo
    z_c(kmax+1) = z_h(kmax+1) + ( z_h(kmax+1) - z_h(kmax) )*0.5_RP
    !
  end subroutine mk_layer_powerfunc

  subroutine mk_layer_given( infname )
    implicit none
    character(len=*), intent(in) :: infname
    integer,parameter :: fid=10
    integer :: k
    !
    open(fid,file=trim(infname),status='old',form='formatted')
    read(fid,*) num_of_layer
    !
    kmin=kdum+1
    kmax=kdum+num_of_layer
    kall=kdum+num_of_layer+kdum
    allocate(z_c(kall))
    allocate(z_h(kall))
    !
    do k = kmin, kmax+1
       read(fid,*) z_h(k)
    enddo
    !
    z_h(kmin-1) = z_h(kmin) - ( z_h(kmin+1) - z_h(kmin) )
    !
    do k= kmin-1, kmax
       z_c(k) = z_h(k) + ( z_h(k+1) - z_h(k) )*0.5_RP
    enddo
    z_c(kmax+1) = z_h(kmax+1) + ( z_h(kmax+1) - z_h(kmax) )*0.5_RP
    !
    close(fid)
    !
  end subroutine mk_layer_given

  !-----------------------------------------------------------------------------
  subroutine output_layer( outfname )
    implicit none

    character(len=*), intent(in) :: outfname

    character(len=H_LONG) :: fname_all
    character(len=H_LONG) :: fname_def
    integer :: kall, kmin, kmax

    integer :: fid = 10
    integer :: k
    !---------------------------------------------------------------------------

    open( unit = fid,            &
          file = trim(outfname), &
          form = 'unformatted'   )
       write(fid) num_of_layer
       write(fid) z_c
       write(fid) z_h
    close(fid)

    kmin = kdum + 1
    kmax = kdum + num_of_layer
    kall = kdum + num_of_layer + kdum

    if ( kall > 100 ) then
       write(fname_all,'(A5,I3.3,A4)') 'ZSALL', kall,         '.txt'
       write(fname_def,'(A5,I3.3,A4)') 'ZSDEF', num_of_layer, '.txt'
    else
       write(fname_all,'(A5,I2.2,A4)') 'ZSALL', kall,         '.txt'
       write(fname_def,'(A5,I2.2,A4)') 'ZSDEF', num_of_layer, '.txt'
    endif

    open( unit = fid,             &
          file = trim(fname_all), &
          form = 'formatted'      )
       write(fid,'(I4)') kall
       do k = 1, kall
          write(fid,'(f12.3)') z_c(k)
       enddo
    close(fid)

    open( unit = fid,             &
          file = trim(fname_def), &
          form = 'formatted'      )
       write(fid,'(I4)') num_of_layer
       do k = kmin, kmax
          write(fid,'(f12.3)') z_c(k)
       enddo
    close(fid)

    return
  end subroutine output_layer

end program prg_mkvlayer
!-------------------------------------------------------------------------------
