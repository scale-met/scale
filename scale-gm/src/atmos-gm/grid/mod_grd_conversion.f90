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
  use scale_stdio

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
    use mod_adm, only: &
         ijdim => ADM_gall_in,  &
         kall  => ADM_kall
    use scale_grid_index, only: & 
         KA,    &
         IA,    &
         JA
    implicit none

    integer,  intent(in) :: kdim
    real(RP), intent(in) :: var_gmgrd  (ijdim,kdim)
    real(RP), intent(out):: var_rmgrd  (kdim,IA,JA)
    integer :: i, j, k, ij

    var_rmgrd(:,:,:) = 0.0_RP
    do k=1, kdim
       ij=0
       do j=1, JA
       do i=1, IA
          if( i/=1 .and. j/=1 ) then
             ij=ij+1
             var_rmgrd(k,i,j) = var_gmgrd(ij,k)
          else
             var_rmgrd(k,i,j) = sum(var_gmgrd(:,k)) / ijdim
          end if
       enddo
       enddo
    enddo
  end subroutine grd_gm2rm

  !-----------------------------------------------------------------------------
  ! > conversion from RM grid to GM grid
  subroutine grd_rm2gm( &
       var_rmgrd, var_gmgrd, kdim )
    use mod_adm, only: &
         ijdim => ADM_gall_in
    use scale_grid_index, only: &
         IA,   &
         JA
    implicit none

    integer, intent(in) :: kdim
    real(RP), intent(in) :: var_rmgrd  (kdim,IA,JA)
    real(RP), intent(out):: var_gmgrd  (ijdim,kdim)
    integer :: i, j, k, ij

    do k=1, kdim
       ij=0
       do j=1, JA
       do i=1, IA
          if( i/=1 .and. j/=1 ) then
             ij=ij+1
             var_gmgrd(ij,k) = var_rmgrd(k,i,j)
          end if
       enddo
       enddo
    enddo
  end subroutine grd_rm2gm

end module mod_grd_conversion

