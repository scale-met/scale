!-------------------------------------------------------------------------------
!> Module SNO (RM) mapping
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          SCALE NetCDF Operator (SNO)
!!          module for mapping
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_sno_map
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SNO_map_settable_global
  public :: SNO_map_settable_local

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine SNO_map_settable_global( &
       nprocs_x_in, &
       nprocs_y_in, &
       ngrids_x,    &
       ngrids_y,    &
       nhalos_x,    &
       nhalos_y,    &
       globalmap,   &
       debug        )
    use mod_sno_h, only: &
       I_map_p, &
       I_map_i, &
       I_map_j
    implicit none

    integer,    intent(in)  :: nprocs_x_in                    ! x length of 2D processor topology (input)
    integer,    intent(in)  :: nprocs_y_in                    ! y length of 2D processor topology (input)
    integer,    intent(in)  :: ngrids_x                       ! size of x-axis grids              (global,sometimes including halo)
    integer,    intent(in)  :: ngrids_y                       ! size of y-axis grids              (global,sometimes including halo)
    integer,    intent(in)  :: nhalos_x                       ! size of x-axis halo grids         (global,sometimes have a size)
    integer,    intent(in)  :: nhalos_y                       ! size of y-axis halo grids         (global,sometimes have a size)
    integer(2), intent(out) :: globalmap(ngrids_x,ngrids_y,3) ! mapping table                     (global)
    logical,    intent(in)  :: debug

    integer :: ngrids_x_in ! number of x-axis grids per process (input,sometimes including halo)
    integer :: ngrids_y_in ! number of y-axis grids per process (input,sometimes including halo)

    integer :: px, py
    integer :: ipos, jpos
    integer :: p, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SNO_map_settable_global",*) '[SNO_map_settable_global] Calc relation map for grids and procs (global)'

    jpos = 0
    do py = 1, nprocs_y_in
       ngrids_y_in = ( ngrids_y - 2*nhalos_y ) / nprocs_y_in
       if( py == 1           ) ngrids_y_in = ngrids_y_in + nhalos_y
       if( py == nprocs_y_in ) ngrids_y_in = ngrids_y_in + nhalos_y

       do j = 1, ngrids_y_in
          jpos = jpos + 1

          ipos = 0
          do px = 1, nprocs_x_in
             ngrids_x_in = ( ngrids_x - 2*nhalos_x ) / nprocs_x_in
             if( px == 1           ) ngrids_x_in = ngrids_x_in + nhalos_x
             if( px == nprocs_x_in ) ngrids_x_in = ngrids_x_in + nhalos_x

             do i = 1, ngrids_x_in
                ipos = ipos + 1

                p = (py-1) * nprocs_x_in + px - 1 ! serialize

                globalmap(ipos,jpos,I_map_p) = int(p,kind=2) ! rank
                globalmap(ipos,jpos,I_map_i) = int(i,kind=2) ! i-index
                globalmap(ipos,jpos,I_map_j) = int(j,kind=2) ! j-index
             enddo
          enddo
       enddo
    enddo

    if ( debug ) then
       LOG_INFO("SNO_map_settable_global",*) 'globalmap(rank)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
       do i = 1, ngrids_x
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
       enddo
       LOG_NEWLINE
       do j = 1, ngrids_y
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
          do i = 1, ngrids_x
             if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') globalmap(i,j,I_map_p)
          enddo
          LOG_NEWLINE
       enddo

       LOG_INFO("SNO_map_settable_global",*) 'globalmap(i-index)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
       do i = 1, ngrids_x
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
       enddo
       LOG_NEWLINE
       do j = 1, ngrids_y
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
          do i = 1, ngrids_x
             if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') globalmap(i,j,I_map_i)
          enddo
          LOG_NEWLINE
       enddo

       LOG_INFO("SNO_map_settable_global",*) 'globalmap(j-index)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
       do i = 1, ngrids_x
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
       enddo
       LOG_NEWLINE
       do j = 1, ngrids_y
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
          do i = 1, ngrids_x
             if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') globalmap(i,j,I_map_j)
          enddo
          LOG_NEWLINE
       enddo
    endif

    return
  end subroutine SNO_map_settable_global

  !-----------------------------------------------------------------------------
  subroutine SNO_map_settable_local( &
       nprocs_x_in,  &
       nprocs_y_in,  &
       ngrids_x,     &
       ngrids_y,     &
       ngrids_x_out, &
       ngrids_y_out, &
       ipos,         &
       jpos,         &
       globalmap,    &
       localmap,     &
       readflag,     &
       debug         )
    use mod_sno_h, only: &
       I_map_p, &
       I_map_i, &
       I_map_j
    implicit none

    integer,    intent(in)  :: nprocs_x_in                            ! x length of 2D processor topology (input)
    integer,    intent(in)  :: nprocs_y_in                            ! y length of 2D processor topology (input)
    integer,    intent(in)  :: ngrids_x                               ! size of x-axis grids              (global,sometimes including halo)
    integer,    intent(in)  :: ngrids_y                               ! size of y-axis grids              (global,sometimes including halo)
    integer,    intent(in)  :: ngrids_x_out                           ! size of x-axis grids              (output,sometimes including halo)
    integer,    intent(in)  :: ngrids_y_out                           ! size of y-axis grids              (output,sometimes including halo)
    integer,    intent(in)  :: ipos                                   ! offset of i-index
    integer,    intent(in)  :: jpos                                   ! offset of j-index
    integer(2), intent(in)  :: globalmap(ngrids_x    ,ngrids_y    ,3) ! mapping table                     (global)
    integer(2), intent(out) :: localmap (ngrids_x_out,ngrids_y_out,3) ! mapping table                     (output)
    logical,    intent(out) :: readflag (nprocs_x_in ,nprocs_y_in )   ! flag to read each input file
    logical,    intent(in)  :: debug

    integer :: p, px, py
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SNO_map_settable_local",*) '[SNO_map_settable_local] Calc relation map for grids and procs (local)'

    do j = 1, ngrids_y_out
    do i = 1, ngrids_x_out
       localmap(i,j,I_map_p) = globalmap(i+ipos,j+jpos,I_map_p)
       localmap(i,j,I_map_i) = globalmap(i+ipos,j+jpos,I_map_i)
       localmap(i,j,I_map_j) = globalmap(i+ipos,j+jpos,I_map_j)
    enddo
    enddo

    do py = 1, nprocs_y_in
    do px = 1, nprocs_x_in
       p = (py-1) * nprocs_x_in + px - 1

       readflag(px,py) = .false.

       do j = 1, ngrids_y_out
       do i = 1, ngrids_x_out
          if( localmap(i,j,I_map_p) == p ) readflag(px,py) = .true.
       enddo
       enddo
    enddo
    enddo

    if ( debug ) then
       LOG_INFO("SNO_map_settable_local",*) 'localmap(rank)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
       do i = 1, ngrids_x_out
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
       enddo
       LOG_NEWLINE
       do j = 1, ngrids_y_out
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
          do i = 1, ngrids_x_out
             if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap(i,j,I_map_p)
          enddo
          LOG_NEWLINE
       enddo

       LOG_INFO("SNO_map_settable_local",*) 'localmap(i-index)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
       do i = 1, ngrids_x_out
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
       enddo
       LOG_NEWLINE
       do j = 1, ngrids_y_out
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
          do i = 1, ngrids_x_out
             if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap(i,j,I_map_i)
          enddo
          LOG_NEWLINE
       enddo

       LOG_INFO("SNO_map_settable_local",*) 'localmap(j-index)'
       if( IO_L ) write(IO_FID_LOG,'(1x,A3)',advance='no') "###"
       do i = 1, ngrids_x_out
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') i
       enddo
       LOG_NEWLINE
       do j = 1, ngrids_y_out
          if( IO_L ) write(IO_FID_LOG,'(1x,I3.3)',advance='no') j
          do i = 1, ngrids_x_out
             if( IO_L ) write(IO_FID_LOG,'(1x,I3)',advance='no') localmap(i,j,I_map_j)
          enddo
          LOG_NEWLINE
       enddo

       LOG_INFO("SNO_map_settable_local",*) 'readflag'
       do py = 1, nprocs_y_in
          do px = 1, nprocs_x_in
             if( IO_L ) write(IO_FID_LOG,'(1x,L2)',advance='no') readflag(px,py)
          enddo
          LOG_NEWLINE
       enddo
    endif

    return
  end subroutine SNO_map_settable_local

end module mod_sno_map
