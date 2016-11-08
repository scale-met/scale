!-------------------------------------------------------------------------------
!> Module vertical interpolation
!!
!! @par Description
!!          This module contains vertical interpolation subroutines
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_vintrpl
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: VINTRPL_Z2Xi
  public :: VINTRPL_Xi2Z

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_Z2Xi( &
       var,   &
       var_pl )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       Xi => GRD_gz, &
       GRD_Z,        &
       GRD_vz,       &
       GRD_vz_pl
    implicit none

    real(RP), intent(inout) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: tmp(ADM_kall)
    real(RP) :: Z  (ADM_kall)

    integer  :: g, k, l, kk
    !---------------------------------------------------------------------------

    do l = 1, ADM_lall
    do g = 1, ADM_gall
       Z  (:)          = GRD_vz(g,:,l,GRD_Z)
       tmp(:)          = var   (g,:,l)
       tmp(ADM_kmin-1) = var   (g,ADM_kmin,l)
       tmp(ADM_kmax+1) = var   (g,ADM_kmax,l)

       do k = ADM_kmin, ADM_kmax
          do kk = ADM_kmin, ADM_kmax-1
             if( Z(k) <= Xi(kk) ) exit
          enddo

          kk = max(kk,ADM_kmin+1)

          var(g,k,l) = lagrange3( Z(k), Xi (kk), Xi (kk-1), Xi (kk-2), &
                                        tmp(kk), tmp(kk-1), tmp(kk-2)  )
       enddo

       var(g,ADM_kmin-1,l) = var(g,ADM_kmin,l)
       var(g,ADM_kmax+1,l) = var(g,ADM_kmax,l)
    enddo
    enddo

    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
       do g = 1, ADM_gall_pl
          Z  (:)          = GRD_vz_pl(g,:,l,GRD_Z)
          tmp(:)          = var_pl   (g,:,l)
          tmp(ADM_kmin-1) = var_pl   (g,ADM_kmin,l)
          tmp(ADM_kmax+1) = var_pl   (g,ADM_kmax,l)

          do k = ADM_kmin, ADM_kmax
             do kk = ADM_kmin, ADM_kmax-1
                if( Z(k) <= Xi(kk) ) exit
             enddo

             kk = max(kk,ADM_kmin+1)

             var_pl(g,k,l) = lagrange3( Z(k), Xi (kk), Xi (kk-1), Xi (kk-2), &
                                              tmp(kk), tmp(kk-1), tmp(kk-2)  )
          enddo

          var_pl(g,ADM_kmin-1,l) = var_pl(g,ADM_kmin,l)
          var_pl(g,ADM_kmax+1,l) = var_pl(g,ADM_kmax,l)
       enddo
       enddo
    else
       var_pl(:,:,:) = UNDEF
       var_pl(:,:,:) = UNDEF
    endif

    return
  end subroutine VINTRPL_Z2Xi

  !-----------------------------------------------------------------------------
  subroutine VINTRPL_Xi2Z( &
       var,     &
       var_pl,  &
       use_quad )
    use scale_const, only: &
       UNDEF => CONST_UNDEF, &
       EPS   => CONST_EPS
    use mod_adm, only: &
       ADM_have_pl, &
       ADM_lall,    &
       ADM_lall_pl, &
       ADM_gall,    &
       ADM_gall_pl, &
       ADM_kall,    &
       ADM_kmin,    &
       ADM_kmax
    use mod_grd, only: &
       Xi => GRD_gz, &
       GRD_Z,        &
       GRD_vz,       &
       GRD_vz_pl
    implicit none

    real(RP), intent(inout) :: var   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical,  intent(in)    :: use_quad

    real(RP) :: tmp(ADM_kall)
    real(RP) :: Z  (ADM_kall)

    integer  :: g, k, l, kk
    !---------------------------------------------------------------------------

    if ( use_quad ) then ! 2nd order

       do l = 1, ADM_lall
       do g = 1, ADM_gall
          Z  (:)          = GRD_vz(g,:,l,GRD_Z)
          tmp(:)          = var   (g,:,l)
          tmp(ADM_kmin-1) = var   (g,ADM_kmin,l)
          tmp(ADM_kmax+1) = var   (g,ADM_kmax,l)

          do k = ADM_kmin, ADM_kmax
             if    ( Xi(k) < Z(ADM_kmin-1) ) then
                var(g,k,l) = UNDEF
             elseif( Xi(k) > Z(ADM_kmax+1) ) then
                var(g,k,l) = UNDEF
             else
                do kk = ADM_kmin, ADM_kmax
                   if( Xi(k) <= Z(kk) ) exit
                enddo

                if    ( tmp(kk)   <= UNDEF+EPS ) then
                   var(g,k,l) = tmp(kk-1)
                elseif( tmp(kk-1) <= UNDEF+EPS ) then
                   var(g,k,l) = tmp(kk)
                elseif( kk == ADM_kmin ) then
                   var(g,k,l) = lagrange2( Z(k), Xi (kk), Xi (kk-1), &
                                                 tmp(kk), tmp(kk-1)  )
                else
                   if ( tmp(kk-2) <= UNDEF+EPS ) then
                      var(g,k,l) = lagrange2( Z(k), Xi (kk), Xi (kk-1), &
                                                    tmp(kk), tmp(kk-1)  )
                   else
                      var(g,k,l) = lagrange3( Z(k), Xi (kk), Xi (kk-1), Xi (kk-2), &
                                                    tmp(kk), tmp(kk-1), tmp(kk-2)  )
                   endif
                endif
             endif
          enddo

          var(g,ADM_kmin-1,l) = var(g,ADM_kmin,l)
          var(g,ADM_kmax+1,l) = var(g,ADM_kmax,l)
       enddo
       enddo

       if ( ADM_have_pl ) Then
          do l = 1, ADM_lall_pl
          do g = 1, ADM_gall_pl
             Z  (:)          = GRD_vz_pl(g,:,l,GRD_Z)
             tmp(:)          = var_pl   (g,:,l)
             tmp(ADM_kmin-1) = var_pl   (g,ADM_kmin,l)
             tmp(ADM_kmax+1) = var_pl   (g,ADM_kmax,l)

             do k = ADM_kmin, ADM_kmax
                if    ( Xi(k) < Z(ADM_kmin-1) ) then
                   var_pl(g,k,l) = UNDEF
                elseif( Xi(k) > Z(ADM_kmax+1) ) then
                   var_pl(g,k,l) = UNDEF
                else
                   do kk = ADM_kmin, ADM_kmax
                      if( Xi(k) <= Z(kk) ) exit
                   enddo

                   if    ( tmp(kk)   <= UNDEF+EPS ) then
                      var_pl(g,k,l) = tmp(kk-1)
                   elseif( tmp(kk-1) <= UNDEF+EPS ) then
                      var_pl(g,k,l) = tmp(kk)
                   elseif( kk == ADM_kmin ) then
                      var_pl(g,k,l) = lagrange2( Z(k), Xi (kk), Xi (kk-1), &
                                                       tmp(kk), tmp(kk-1)  )
                   else
                      if ( tmp(kk-2) <= UNDEF+EPS ) then
                         var_pl(g,k,l) = lagrange2( Z(k), Xi (kk), Xi (kk-1), &
                                                          tmp(kk), tmp(kk-1)  )
                      else
                         var_pl(g,k,l) = lagrange3( Z(k), Xi (kk), Xi (kk-1), Xi (kk-2), &
                                                          tmp(kk), tmp(kk-1), tmp(kk-2)  )
                      endif
                   endif
                endif
             enddo

             var_pl(g,ADM_kmin-1,l) = var_pl(g,ADM_kmin,l)
             var_pl(g,ADM_kmax+1,l) = var_pl(g,ADM_kmax,l)
          enddo
          enddo
       else
          var_pl(:,:,:) = UNDEF
          var_pl(:,:,:) = UNDEF
       endif

    else ! 1st order

       do l = 1, ADM_lall
       do g = 1, ADM_gall
          Z  (:)          = GRD_vz(g,:,l,GRD_Z)
          tmp(:)          = var   (g,:,l)
          tmp(ADM_kmin-1) = var   (g,ADM_kmin,l)
          tmp(ADM_kmax+1) = var   (g,ADM_kmax,l)

          do k = ADM_kmin, ADM_kmax
             if    ( Xi(k) < Z(ADM_kmin-1) ) then
                var(g,k,l) = UNDEF
             elseif( Xi(k) > Z(ADM_kmax+1) ) then
                var(g,k,l) = UNDEF
             else
                do kk = ADM_kmin, ADM_kmax
                   if( Xi(k) <= Z(kk) ) exit
                enddo

                if    ( tmp(kk)   <= UNDEF+EPS ) then
                   var(g,k,l) = tmp(kk-1)
                elseif( tmp(kk-1) <= UNDEF+EPS ) then
                   var(g,k,l) = tmp(kk)
                else
                   var(g,k,l) = lagrange2( Z(k), Xi (kk), Xi (kk-1), &
                                                 tmp(kk), tmp(kk-1)  )
                endif
             endif
          enddo

          var(g,ADM_kmin-1,l) = var(g,ADM_kmin,l)
          var(g,ADM_kmax+1,l) = var(g,ADM_kmax,l)
       enddo
       enddo

       if ( ADM_have_pl ) Then
          do l = 1, ADM_lall_pl
          do g = 1, ADM_gall_pl
             Z  (:)          = GRD_vz_pl(g,:,l,GRD_Z)
             tmp(:)          = var_pl   (g,:,l)
             tmp(ADM_kmin-1) = var_pl   (g,ADM_kmin,l)
             tmp(ADM_kmax+1) = var_pl   (g,ADM_kmax,l)

             do k = ADM_kmin, ADM_kmax
                if    ( Xi(k) < Z(ADM_kmin-1) ) then
                   var_pl(g,k,l) = UNDEF
                elseif( Xi(k) > Z(ADM_kmax+1) ) then
                   var_pl(g,k,l) = UNDEF
                else
                   do kk = ADM_kmin, ADM_kmax
                      if( Xi(k) <= Z(kk) ) exit
                   enddo

                   if    ( tmp(kk)   <= UNDEF+EPS ) then
                      var_pl(g,k,l) = tmp(kk-1)
                   elseif( tmp(kk-1) <= UNDEF+EPS ) then
                      var_pl(g,k,l) = tmp(kk)
                   else
                      var_pl(g,k,l) = lagrange2( Z(k), Xi (kk), Xi (kk-1), &
                                                       tmp(kk), tmp(kk-1)  )
                   endif
                endif
             enddo

             var_pl(g,ADM_kmin-1,l) = var_pl(g,ADM_kmin,l)
             var_pl(g,ADM_kmax+1,l) = var_pl(g,ADM_kmax,l)
          enddo
          enddo
       else
          var_pl(:,:,:) = UNDEF
          var_pl(:,:,:) = UNDEF
       endif

    endif

    return
  end subroutine VINTRPL_Xi2Z

  !-----------------------------------------------------------------------------
  function lagrange3( &
       z,          &
       z1, z2, z3, &
       p1, p2, p3  ) &
       result(lagrange)
    implicit none

    real(RP), intent(in) :: z, z1, p1, z2, p2, z3, p3
    real(RP)             :: lagrange
    !---------------------------------------------------------------------------

    lagrange = ( (z-z2)*(z-z3) ) / ( (z1-z2)*(z1-z3) ) * p1 &
             + ( (z-z1)*(z-z3) ) / ( (z2-z1)*(z2-z3) ) * p2 &
             + ( (z-z1)*(z-z2) ) / ( (z3-z1)*(z3-z2) ) * p3

  end function lagrange3

  !-----------------------------------------------------------------------------
  function lagrange2( &
       z,      &
       z1, z2, &
       p1, p2  ) &
       result(lagrange)
    implicit none

    real(RP), intent(in) :: z, z1, p1, z2, p2
    real(RP)             :: lagrange
    !---------------------------------------------------------------------------

    lagrange = (z-z2) / (z1-z2) * p1 &
             + (z-z1) / (z2-z1) * p2

  end function lagrange2

end module mod_vintrpl
