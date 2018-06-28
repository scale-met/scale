!-------------------------------------------------------------------------------
!> Module grid
!!
!! @par Description
!!         This module is for the management of the icosahedral grid system
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_grd_conversion
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_atmos_grid_icoA_index

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: GRD_GM2RM
  public :: GRD_RM2GM

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
  ! > conversion from GM grid to RM grid
  subroutine grd_gm2rm( &
       var_gmgrd, var_rmgrd, kdim )
    implicit none

    integer,  intent(in) :: kdim
    real(RP), intent(in) :: var_gmgrd  (ADM_gall_in,kdim)
    real(RP), intent(out):: var_rmgrd  (kdim,IA,JA)
    integer :: i, j, k, ij

    do j = JS, JE
    do i = IS, IE
       ij = i + (j-1) * ADM_imax
       do k = 1, ADM_kall
          var_rmgrd(k,i,j) = var_gmgrd(ij,k)
       enddo
    enddo
    enddo

    return
  end subroutine grd_gm2rm

  !-----------------------------------------------------------------------------
  ! > conversion from RM grid to GM grid
  subroutine grd_rm2gm( &
       var_rmgrd, var_gmgrd, kdim )
    implicit none

    integer, intent(in) :: kdim
    real(RP), intent(in) :: var_rmgrd  (kdim,IA,JA)
    real(RP), intent(out):: var_gmgrd  (ADM_gall_in,kdim)
    integer :: i, j, k, ij

    do k = 1, kdim
       do j = JS, JE
       do i = IS, IE
          ij = i + (j-1) * ADM_imax
          var_gmgrd(ij,k) = var_rmgrd(k,i,j)
       enddo
       enddo
    enddo
  end subroutine grd_rm2gm

end module mod_grd_conversion

