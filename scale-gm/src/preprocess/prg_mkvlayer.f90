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
  use mpi
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  integer, parameter :: kdum = 1
  integer, parameter :: fid  = 11

  integer                :: num_of_layer = 10        ! number of layers
  character(len=H_SHORT) :: layer_type   = 'POWER'   ! type of layer
  real(RP)               :: ztop         = 1.E4_RP   ! height of model top if layer_type='POWER'
  real(RP)               :: fact         = 1.0_RP    ! factor              if layer_type='POWER'
  character(len=H_LONG)  :: infname      = 'infile'  ! input  file name    if layer_type='GIVEN'
  character(len=H_LONG)  :: outfname     = 'outfile' ! output file name

  namelist / mkvlayer_cnf / &
       num_of_layer, &
       layer_type,   &
       ztop,         &
       fact,         &
       infname,      &
       outfname

  real(RP), allocatable :: z_c(:)
  real(RP), allocatable :: z_h(:)

  integer :: kall, kmin, kmax
  integer :: k
  integer :: ierr
  !=============================================================================

  call MPI_Init(ierr)

  open( unit   = fid,            &
        file   = 'mkvlayer.cnf', &
        status = 'old',          &
        form   = 'formatted'     )

     read(fid,nml=mkvlayer_cnf)

  close(fid)

  kmin = kdum + 1
  kmax = kdum + num_of_layer
  kall = kdum + num_of_layer + kdum

  allocate( z_h(kall) )
  allocate( z_c(kall) )

  select case(layer_type)
  case('POWER')

     call mk_layer_powerfunc( ztop, fact )

  case('GIVEN')

     call mk_layer_given( infname )

  end select

  z_h(kmin-1) = z_h(kmin) - ( z_h(kmin+1) - z_h(kmin) )

  do k = kmin-1, kmax
     z_c(k) = z_h(k) + 0.5_RP * ( z_h(k+1) - z_h(k) )
  enddo
  z_c(kmax+1) = z_h(kmax+1) + 0.5_RP * ( z_h(kmax+1) - z_h(kmax) )

  call output_layer(outfname)

  call MPI_Finalize(ierr)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine mk_layer_powerfunc( ztop, fact )
    implicit none

    real(RP), intent(in) :: ztop
    real(RP), intent(in) :: fact

    real(RP) :: a
    integer  :: k
    !---------------------------------------------------------------------------

    a = ztop / real(num_of_layer,kind=RP)**fact
    do k = kmin, kmax+1
       z_h(k) = a * real(k-kmin,kind=RP)**fact
    enddo

    return
  end subroutine mk_layer_powerfunc

  !-----------------------------------------------------------------------------
  subroutine mk_layer_given( infname )
    implicit none

    character(len=*), intent(in) :: infname

    integer, parameter :: fid = 10

    integer :: num_of_layer0
    integer :: k
    !---------------------------------------------------------------------------

    open( unit   = fid,           &
          file   = trim(infname), &
          status = 'old',         &
          form   = 'formatted'    )

       read(fid,*) num_of_layer0

       if ( num_of_layer0 /= num_of_layer ) then
          write(*,*) 'Mismach num_of_layer (input,request)=',num_of_layer0,num_of_layer
       endif

       do k = kmin, kmax+1
          read(fid,*) z_h(k)
       enddo

    close(fid)

    return
  end subroutine mk_layer_given

  !-----------------------------------------------------------------------------
  subroutine output_layer( outfname )
    implicit none

    character(len=*), intent(in) :: outfname

    character(len=H_LONG) :: fname_all
    character(len=H_LONG) :: fname_def

    integer, parameter :: fid = 10

    integer :: k
    !---------------------------------------------------------------------------

    open( unit   = fid,            &
          file   = trim(outfname), &
          status = 'replace',      &
          form   = 'unformatted'   )

       write(fid) num_of_layer
       write(fid) z_c
       write(fid) z_h

    close(fid)

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
