!-------------------------------------------------------------------------------
!
!+  vertical interpolation
!
!-------------------------------------------------------------------------------
module mod_vintrpl
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module contains vertical interpolation subroutines.
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version    Date      Comment
  !      -----------------------------------------------------------------------
  !      0.00       04-02-17  Imported from igdc-4.33
  !                 07-01-24  added VINTRPL_z_level2 (linear interpolation) and
  !                           consider undefined value (K.Suzuki)
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_adm, only: &
     ADM_LOG_FID, &
       ADM_gall,        &
       ADM_kall,        &
       ADM_lall,        &
       ADM_gall_pl,     &
       ADM_lall_pl,     &
       ADM_KNONE,       &
       ADM_kmin,        &
       ADM_kmax,        &
       ADM_have_pl
  use mod_grd, only :   &
       GRD_vz,GRD_vz_pl,&
       GRD_gz,          &
       GRD_ZH,GRD_Z,    &
       GRD_gzh
  use mod_cnst, only :  &
       CNST_UNDEF
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters, variables, and subroutines
  !
  public :: VINTRPL_srfc_val
  public :: VINTRPL_z_level
  public :: VINTRPL_z_level2
  public :: VINTRPL_zstar_level
  public :: VINTRPL_sigma_level
  public :: VINTRPL_mk_sigma
  public :: VINTRPL_half2full
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters, variables, and subroutines
  !
  ! < NONE >
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_srfc_val( &
       var, var_pl,            & !--- IN    : 3d data
       svar, svar_pl,          & !--- OUT : surface data
       z_offset_in             & !--- IN  : z offset
       )
    !------
    !------ Caluculate surface value.( e.g. surface pressure )
    implicit none
    real(RP), intent(in)  :: var(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(in)  :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(out) :: svar(ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(out) :: svar_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    real(RP), intent(in), optional  :: z_offset_in
    !
    real(RP) :: lag_intpl
    real(RP) :: z,z1,p1,z2,p2,z3,p3
    lag_intpl(z,z1,p1,z2,p2,z3,p3)             &
         = ((z-z2)*(z-z3))/((z1-z2)*(z1-z3))*p1&
         + ((z-z1)*(z-z3))/((z2-z1)*(z2-z3))*p2&
         + ((z-z1)*(z-z2))/((z3-z1)*(z3-z2))*p3
    !
    integer :: l,n
    real(RP) :: z_offset
    !
    if(present(z_offset_in)) then
       z_offset = z_offset
    else
       z_offset = 0.0_RP
    endif
    !
    do l=1,ADM_lall
       do n = 1, ADM_gall
          svar(n,ADM_KNONE,l) &
               = lag_intpl(GRD_vz(n,ADM_kmin,l,GRD_ZH)+z_offset,&
               GRD_vz(n,ADM_kmin+2,l,GRD_Z),var(n,ADM_kmin+2,l),&
               GRD_vz(n,ADM_kmin+1,l,GRD_Z),var(n,ADM_kmin+1,l),&
               GRD_vz(n,ADM_kmin  ,l,GRD_Z),var(n,ADM_kmin  ,l))
       enddo
    enddo
    if( ADM_have_pl ) then
       do l=1,ADM_lall_pl
          do n = 1, ADM_gall_pl
             svar_pl(n,ADM_KNONE,l) &
                  = lag_intpl(GRD_vz_pl(n,ADM_kmin,l,GRD_ZH)+z_offset,&
                  GRD_vz_pl(n,ADM_kmin+2,l,GRD_Z),var_pl(n,ADM_kmin+2,l),&
                  GRD_vz_pl(n,ADM_kmin+1,l,GRD_Z),var_pl(n,ADM_kmin+1,l),&
                  GRD_vz_pl(n,ADM_kmin  ,l,GRD_Z),var_pl(n,ADM_kmin  ,l))
          enddo
       enddo
    endif
  end subroutine VINTRPL_srfc_val
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_z_level( v, v_pl, wgrid )
    implicit none
    real(RP), intent(inout) :: v(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(inout) :: v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical, intent(in) :: wgrid

    real(RP) :: tmp(ADM_gall,ADM_kall,ADM_lall)
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer :: l,k,n,kk

    integer :: kp
    real(RP) :: lag_intpl_quadra, lag_intpl_linear
    real(RP) :: z,z1,p1,z2,p2,z3,p3
    lag_intpl_quadra(z,z1,p1,z2,p2,z3,p3)       &
         = ((z-z2)*(z-z3))/((z1-z2)*(z1-z3))*p1 &
         + ((z-z1)*(z-z3))/((z2-z1)*(z2-z3))*p2 &
         + ((z-z1)*(z-z2))/((z3-z1)*(z3-z2))*p3
    lag_intpl_linear(z,z1,p1,z2,p2) &
         = (z-z2)/(z1-z2)*p1 &
         + (z1-z)/(z1-z2)*p2

    if(wgrid) then
       ! 07/01/24 K.Suzuki: consider undefined value
       do l=1,ADM_lall
       do k=ADM_kmin,ADM_kmax
       do n=1,ADM_gall
          if ( v(n,k+1,l) == CNST_UNDEF ) then
            tmp(n,k,l) = v(n,k,l)
          else if ( v(n,k,l) == CNST_UNDEF ) then
            tmp(n,k,l) = v(n,k+1,l)
          else
            tmp(n,k,l) = 0.5_RP*(v(n,k+1,l)+v(n,k,l))
          endif
       enddo
       enddo
       enddo
       tmp(:,ADM_kmin-1,:) = v(:,ADM_kmin,:)
       tmp(:,ADM_kmax+1,:) = v(:,ADM_kmax+1,:)

       if ( ADM_have_pl ) Then
          ! 07/01/24 K.Suzuki: consider undefined value
          do l=1,ADM_lall
          do k=ADM_kmin,ADM_kmax
          do n=1,ADM_gall
            if ( v_pl(n,k+1,l) == CNST_UNDEF ) then
              tmp_pl(n,k,l) = v_pl(n,k,l)
            else if ( v_pl(n,k,l) == CNST_UNDEF ) then
              tmp_pl(n,k,l) = v_pl(n,k+1,l)
            else
              tmp_pl(n,k,l) = 0.5_RP*(v_pl(n,k+1,l)+v_pl(n,k,l))
            endif
          enddo
          enddo
          enddo
          tmp_pl(:,ADM_kmin-1,:) = v_pl(:,ADM_kmin,:)
          tmp_pl(:,ADM_kmax+1,:) = v_pl(:,ADM_kmax+1,:)
       endif
    else
       tmp = v
       tmp_pl = v_pl
    endif

    Do l=1,ADM_lall
       do k=1,ADM_kall
          do n=1, ADM_gall
             if(      (GRD_gz(k)<GRD_vz(n,ADM_kmin,l,GRD_ZH))&
                  .or.(GRD_gz(k)>GRD_vz(n,ADM_kmax+1,l,GRD_ZH)) ) then
                v(n,k,l) = CNST_UNDEF
             else
                do kk = ADM_kmin-1,ADM_kmax+1
                   kp=kk
                   if(GRD_gz(k)<GRD_vz(n,kk,l,GRD_Z)) exit
                enddo
!                if (kp==ADM_kmin-1) then
!                   kp=ADM_kmin
!                elseif(kp==ADM_kmax+1) then
!                   kp=ADM_kmax
!                endif
                if (kp<=ADM_kmin) then
                   kp=ADM_kmin+1
                elseif(kp>=ADM_kmax) then
                   kp=ADM_kmax-1
                endif
                ! 07/01/24 K.Suzuki: consider undefined value
                if ( tmp(n,kp+1,l) == CNST_UNDEF ) then
                  if ( tmp(n,kp,l) == CNST_UNDEF ) then
                    v(n,k,l) = tmp(n,kp-1,l)
                  else if ( tmp(n,kp-1,l) == CNST_UNDEF ) then
                    v(n,k,l) = tmp(n,kp,l)
                  else
                    v(n,k,l) &
                         = lag_intpl_linear(GRD_gz(k),         &
                         GRD_vz(n,kp  ,l,GRD_Z),tmp(n,kp  ,l), &
                         GRD_vz(n,kp-1,l,GRD_Z),tmp(n,kp-1,l)  )
                  endif
                else if ( tmp(n,kp,l) == CNST_UNDEF ) then
                   v(n,k,l) = tmp(n,kp-1,l)
                else if ( tmp(n,kp-1,l) == CNST_UNDEF ) then
                   v(n,k,l) = tmp(n,kp,l)
                else
                   v(n,k,l)                                    &
                        = lag_intpl_quadra(GRD_gz(k),          &
                        GRD_vz(n,kp+1,l,GRD_Z),tmp(n,kp+1,l),  &
                        GRD_vz(n,kp  ,l,GRD_Z),tmp(n,kp  ,l),  &
                        GRD_vz(n,kp-1,l,GRD_Z),tmp(n,kp-1,l) )
                endif
             endif
          enddo
       enddo
    enddo
    if ( ADM_have_pl ) Then
       Do l=1,ADM_lall_pl
          do k=1,ADM_kall
             do n=1, ADM_gall_pl
                if(  (GRD_gz(k)<GRD_vz_pl(n,ADM_kmin,l,GRD_ZH)).or.&
                     (GRD_gz(k)>GRD_vz_pl(n,ADM_kmax+1,l,GRD_ZH)) ) then
                   v_pl(n,k,l) = CNST_UNDEF
                else
                   do kk = ADM_kmin-1,ADM_kmax+1
                     kp=kk
                     if(GRD_gz(k)<GRD_vz_pl(n,kk,l,GRD_Z)) exit
                   enddo
!                   if(kp==ADM_kmin-1) then
!                      kp=ADM_kmin
!                   elseif(kp==ADM_kmax+1) then
!                      kp=ADM_kmax
!                   endif
                   if (kp<=ADM_kmin) then
                      kp=ADM_kmin+1
                   elseif(kp>=ADM_kmax) then
                      kp=ADM_kmax-1
                   endif
                   ! 07/01/24 K.Suzuki: consider undefined value
                   if ( tmp_pl(n,kp+1,l) == CNST_UNDEF ) then
                     if ( tmp_pl(n,kp,l) == CNST_UNDEF ) then
                       v_pl(n,k,l) = tmp_pl(n,kp-1,l)
                     else if ( tmp_pl(n,kp-1,l) == CNST_UNDEF ) then
                       v_pl(n,k,l) = tmp_pl(n,kp,l)
                     else
                       v_pl(n,k,l)                                  &
                              = lag_intpl_linear(GRD_gz(k),         &
                              GRD_vz_pl(n,kp  ,l,GRD_Z),tmp_pl(n,kp  ,l), &
                              GRD_vz_pl(n,kp-1,l,GRD_Z),tmp_pl(n,kp-1,l) )
                     endif
                   else if ( tmp_pl(n,kp,l) == CNST_UNDEF ) then
                      v_pl(n,k,l) = tmp_pl(n,kp-1,l)
                   else if ( tmp_pl(n,kp-1,l) == CNST_UNDEF ) then
                      v_pl(n,k,l) = tmp_pl(n,kp,l)
                   else
                      v_pl(n,k,l)                                       &
                           = lag_intpl_quadra(GRD_gz(k),                &
                           GRD_vz_pl(n,kp+1,l,GRD_Z),tmp_pl(n,kp+1,l),  &
                           GRD_vz_pl(n,kp  ,l,GRD_Z),tmp_pl(n,kp  ,l),  &
                           GRD_vz_pl(n,kp-1,l,GRD_Z),tmp_pl(n,kp-1,l) )
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine VINTRPL_z_level
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_z_level2( v, v_pl, wgrid ) ! 07/01/24 K.Suzuki [add]
    implicit none
    real(RP), intent(inout) :: v(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(inout) :: v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical, intent(in) :: wgrid

    real(RP) :: tmp(ADM_gall,ADM_kall,ADM_lall)
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer :: l,k,n,kk

    integer :: kp
    real(RP) :: lag_intpl

    real(RP) :: z,z1,p1,z2,p2
    lag_intpl(z,z1,p1,z2,p2)  &
         = (z-z2)/(z1-z2)*p1  &
         + (z1-z)/(z1-z2)*p2

    if(wgrid) then
       ! 07/01/24 K.Suzuki: consider undefined value
       do l=1, ADM_lall
       do k=ADM_kmin,ADM_kmax
       do n=1, ADM_gall
         if ( v(n,k+1,l) == CNST_UNDEF ) then
           tmp(n,k,l) = v(n,k,l)
         else if ( v(n,k,l) == CNST_UNDEF ) then
           tmp(n,k,l) = v(n,k+1,l)
         else
           tmp(n,k,l) = 0.5_RP*(v(n,k+1,l)+v(n,k,l))
         endif
       enddo
       enddo
       enddo
       tmp(:,ADM_kmin-1,:) = v(:,ADM_kmin,:)
       tmp(:,ADM_kmax+1,:) = v(:,ADM_kmax+1,:)

       if ( ADM_have_pl ) Then
          do l=1,ADM_lall
          do k=ADM_kmin,ADM_kmax
          do n=1,ADM_gall
            if ( v_pl(n,k+1,l) == CNST_UNDEF ) then
              tmp_pl(n,k,l) = v_pl(n,k,l)
            else if ( v_pl(n,k,l) == CNST_UNDEF ) then
              tmp_pl(n,k,l) = v_pl(n,k+1,l)
            else
              tmp_pl(n,k,l) = 0.5_RP*(v_pl(n,k+1,l)+v_pl(n,k,l))
            endif
          enddo
          enddo
          enddo
          tmp_pl(:,ADM_kmin-1,:) = v_pl(:,ADM_kmin,:)
          tmp_pl(:,ADM_kmax+1,:) = v_pl(:,ADM_kmax+1,:)
       endif
    else
       tmp = v
       tmp_pl = v_pl
    endif

    Do l=1,ADM_lall
       do k=1,ADM_kall
          do n=1, ADM_gall
             if(      (GRD_gz(k)<GRD_vz(n,ADM_kmin,l,GRD_ZH))&
                  .or.(GRD_gz(k)>GRD_vz(n,ADM_kmax+1,l,GRD_ZH)) ) then
                v(n,k,l) = CNST_UNDEF
             else
                do kk = ADM_kmin-1,ADM_kmax+1
                   kp=kk
                   if(GRD_gz(k)<GRD_vz(n,kk,l,GRD_Z)) exit
                enddo
                if (kp<=ADM_kmin) then
                   v(n,k,l)=tmp(n,ADM_kmin,l)
                elseif(kp>=ADM_kmax+1) then
                   v(n,k,l)=tmp(n,ADM_kmax,l)
                else
                  ! 07/01/24 K.Suzuki: consider undefined value
                  if ( tmp(n,kp,l) == CNST_UNDEF ) then
                     v(n,k,l) = tmp(n,kp-1,l)
                  else if ( tmp(n,kp-1,l) == CNST_UNDEF ) then
                     v(n,k,l) = tmp(n,kp,l)
                  else
                     v(n,k,l)                                    &
                          = lag_intpl(GRD_gz(k),                 &
                          GRD_vz(n,kp  ,l,GRD_Z),tmp(n,kp  ,l),  &
                          GRD_vz(n,kp-1,l,GRD_Z),tmp(n,kp-1,l) )
                  endif
                endif
             endif
          enddo
       enddo
    enddo
    if ( ADM_have_pl ) Then
       Do l=1,ADM_lall_pl
          do k=1,ADM_kall
             do n=1, ADM_gall_pl
                if(  (GRD_gz(k)<GRD_vz_pl(n,ADM_kmin,l,GRD_ZH)).or.&
                     (GRD_gz(k)>GRD_vz_pl(n,ADM_kmax+1,l,GRD_ZH)) ) then
                   v_pl(n,k,l) = CNST_UNDEF
                else
                   do kk = ADM_kmin-1,ADM_kmax+1
                     kp=kk
                     if(GRD_gz(k)<GRD_vz_pl(n,kk,l,GRD_Z)) exit
                   enddo
                   if (kp<=ADM_kmin) then
                      v_pl(n,k,l) = tmp_pl(n,ADM_kmin,l)
                   else if(kp>=ADM_kmax+1) then
                      v_pl(n,k,l) = tmp_pl(n,ADM_kmax,l)
                   else
                     ! 07/01/24 K.Suzuki: consider undefined value
                     if ( tmp_pl(n,kp,l) == CNST_UNDEF ) then
                        v_pl(n,k,l) = tmp_pl(n,kp-1,l)
                     else if ( tmp_pl(n,kp-1,l) == CNST_UNDEF ) then
                        v_pl(n,k,l) = tmp_pl(n,kp,l)
                     else
                        v_pl(n,k,l)                                      &
                             = lag_intpl(GRD_gz(k),                      &
                             GRD_vz_pl(n,kp  ,l,GRD_Z),tmp_pl(n,kp  ,l), &
                             GRD_vz_pl(n,kp-1,l,GRD_Z),tmp_pl(n,kp-1,l) )
                     endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif
  end subroutine VINTRPL_z_level2

  !-----------------------------------------------------------------------------
  subroutine VINTRPL_zstar_level( &
       v,    &
       v_pl, &
       wgrid )
    implicit none

    real(RP), intent(inout) :: v   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(inout) :: v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    logical, intent(in)    :: wgrid

    real(RP) :: tmp   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    integer :: n, k, l, kk, kp

    real(RP) :: lag_intpl
    real(RP) :: z, z1, p1, z2, p2, z3, p3
    lag_intpl(z,z1,p1,z2,p2,z3,p3) = ( (z-z2)*(z-z3) ) / ( (z1-z2)*(z1-z3) ) * p1 &
                                   + ( (z-z1)*(z-z3) ) / ( (z2-z1)*(z2-z3) ) * p2 &
                                   + ( (z-z1)*(z-z2) ) / ( (z3-z1)*(z3-z2) ) * p3
    !---------------------------------------------------------------------------

    if ( wgrid ) then
       tmp = v
       v(:,ADM_kmin-1,:) = CNST_UNDEF
       do l=1,ADM_lall
          do n = 1, ADM_gall
             do k = ADM_kmin, ADM_kmax+1
                do kk = 1, ADM_kall-1
                   if ( ( GRD_gzh(kk)   <= GRD_vz(n,k,l,GRD_ZH) ) .and. &
                        ( GRD_gzh(kk+1) >= GRD_vz(n,k,l,GRD_ZH) ) ) then
                      if(kk==1) then
                         kp=2
                      else
                         kp=kk
                      endif
                      v(n,k,l)                         &
                           = lag_intpl(GRD_vz(n,k,l,GRD_ZH),  &
                            GRD_gzh(kp+1),tmp(n,kp+1,l),&
                            GRD_gzh(kp  ),tmp(n,kp  ,l),&
                            GRD_gzh(kp-1),tmp(n,kp-1,l))
                      exit
                   endif
                enddo
             enddo
          enddo
       enddo

       if ( ADM_have_pl ) Then
          tmp_pl = v_pl
          v_pl(:,ADM_kmin-1,:) = CNST_UNDEF
          do l=1,ADM_lall_pl
             do n = 1, ADM_gall_pl
                do k = ADM_kmin, ADM_kmax+1
                   do kk = 1, ADM_kall-1
                      if ( ( GRD_gzh(kk)   <= GRD_vz_pl(n,k,l,GRD_ZH) ) .and. &
                           ( GRD_gzh(kk+1) >=  GRD_vz_pl(n,k,l,GRD_ZH) ) ) then
                         if(kk==1) then
                            kp=2
                         else
                            kp=kk
                         endif
                         v_pl(n,k,l)                     &
                           = lag_intpl(GRD_vz_pl(n,k,l,GRD_ZH), &
                           GRD_gzh(kp+1),tmp_pl(n,kp+1,l),&
                           GRD_gzh(kp  ),tmp_pl(n,kp  ,l),&
                           GRD_gzh(kp-1),tmp_pl(n,kp-1,l))
                         exit
                      endif
                   enddo
                enddo
             enddo
          enddo
       endif

    else !--- full grid

       tmp(:,:,:) = v(:,:,:)
       v(:,ADM_kmin-1,:) = CNST_UNDEF
       v(:,ADM_kmax+1,:) = CNST_UNDEF

       do l = 1, ADM_lall
       do n = 1, ADM_gall
       do k = ADM_kmin, ADM_kmax
          do kk = 1, ADM_kall-1
             if (       GRD_gz(kk)   <= GRD_vz(n,k,l,GRD_Z) &
                  .AND. GRD_gz(kk+1) >= GRD_vz(n,k,l,GRD_Z) ) then

                if ( kk == 1 ) then
                   kp = 2
                else
                   kp = kk
                endif

                v(n,k,l) = lag_intpl( GRD_vz(n,k,l,GRD_Z),         &
                                      GRD_gz(kp+1), tmp(n,kp+1,l), &
                                      GRD_gz(kp  ), tmp(n,kp  ,l), &
                                      GRD_gz(kp-1), tmp(n,kp-1,l)  )
                exit
             endif
          enddo
       enddo
       enddo
       enddo

       if ( ADM_have_pl ) Then
          tmp_pl = v_pl
          v_pl(:,ADM_kmin-1,:) = CNST_UNDEF
          v_pl(:,ADM_kmax+1,:) = CNST_UNDEF
          do l=1,ADM_lall_pl
             do n = 1, ADM_gall_pl
                do k = ADM_kmin, ADM_kmax
                   do kk = 1, ADM_kall-1
                      if ( ( GRD_gz(kk)   <= GRD_vz_pl(n,k,l,GRD_Z) ) .and. &
                           ( GRD_gz(kk+1) >=  GRD_vz_pl(n,k,l,GRD_Z) ) ) then
                         if(kk==1) then
                            kp=2
                         else
                            kp=kk
                         endif
                         v_pl(n,k,l)                           &
                           = lag_intpl(GRD_vz_pl(n,k,l,GRD_Z), &
                           GRD_gz(kp+1),tmp_pl(n,kp+1,l),      &
                           GRD_gz(kp  ),tmp_pl(n,kp  ,l),      &
                           GRD_gz(kp-1),tmp_pl(n,kp-1,l))
                         exit
                      endif
                   enddo
                enddo
             enddo
          enddo
       endif
    endif
    !
  end subroutine VINTRPL_zstar_level
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_half2full( v, v_pl )
    implicit none
    real(RP), intent(inout) :: v(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(inout) :: v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP) :: tmp(ADM_gall,ADM_kall,ADM_lall)
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer :: k

    do k=ADM_kmin,ADM_kmax
       tmp(:,k,:) = 0.5_RP*(v(:,k+1,:)+v(:,k,:))
    enddo
    tmp(:,ADM_kmin-1,:) = v(:,ADM_kmin,:)
    tmp(:,ADM_kmax+1,:) = v(:,ADM_kmax+1,:)
    v = tmp
    !
    if ( ADM_have_pl ) Then
       do k=ADM_kmin,ADM_kmax
          tmp_pl(:,k,:) = 0.5_RP*(v_pl(:,k+1,:)+v_pl(:,k,:))
       enddo
       tmp_pl(:,ADM_kmin-1,:) = v_pl(:,ADM_kmin,:)
       tmp_pl(:,ADM_kmax+1,:) = v_pl(:,ADM_kmax+1,:)
       v_pl = tmp_pl
    endif
    !
  end subroutine VINTRPL_half2full
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_mk_sigma( sigma, sigma_pl, pre, pre_pl, pres, pres_pl )
    !
    implicit none
    !
    real(RP), intent(out) :: sigma(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(out) :: sigma_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: pre(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(in) :: pre_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !
    real(RP), intent(in) :: pres(ADM_gall,ADM_KNONE,ADM_lall)
    real(RP), intent(in) :: pres_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl)
    integer :: l,k,n
    !
    do l = 1, ADM_lall
       do k = 1,ADM_kall
          do n=1, ADM_gall
             sigma(n,k,l) = pre(n,k,l)/pres(n,ADM_KNONE,l)
          enddo
       enddo
    enddo
    if ( ADM_have_pl ) Then
       do l = 1, ADM_lall_pl
          do k = 1,ADM_kall
             do n=1, ADM_gall_pl
                sigma_pl(n,k,l) = pre_pl(n,k,l)/pres_pl(n,ADM_KNONE,l)
             enddo
          enddo
       enddo
    endif
    return
  end subroutine VINTRPL_mk_sigma
  !-----------------------------------------------------------------------------
  subroutine VINTRPL_sigma_level( &
       v, v_pl,                   &
       sigma, sigma_pl,           &
       sigma_lev, MAX_SIGMA,      &
       wgrid )
    implicit none
    real(RP), intent(inout) :: v(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(inout) :: v_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in) :: sigma(ADM_gall,ADM_kall,ADM_lall)
    real(RP), intent(in) :: sigma_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    integer, intent(in) :: MAX_SIGMA
    real(RP), intent(in) :: sigma_lev(MAX_SIGMA)
    logical, intent(in) :: wgrid

    real(RP) :: tmp(ADM_gall,ADM_kall,ADM_lall)
    real(RP) :: tmp_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP) :: a,b
    integer :: l,k,n,kk,kp

    if(wgrid) then
       do k=ADM_kmin,ADM_kmax
          tmp(:,k,:) = 0.5_RP*(v(:,k+1,:)+v(:,k,:))
       enddo
       tmp(:,ADM_kmin-1,:) = v(:,ADM_kmin,:)
       tmp(:,ADM_kmax+1,:) = v(:,ADM_kmax+1,:)

       if ( ADM_have_pl ) Then
          do k=ADM_kmin,ADM_kmax
             tmp_pl(:,k,:) = 0.5_RP*(v_pl(:,k+1,:)+v_pl(:,k,:))
          enddo
          tmp_pl(:,ADM_kmin-1,:) = v_pl(:,ADM_kmin,:)
          tmp_pl(:,ADM_kmax+1,:) = v_pl(:,ADM_kmax+1,:)
       endif
    else
       tmp = v
       tmp_pl = v_pl
    endif

    Do l=1,ADM_lall
       do k=1,MAX_SIGMA
          do n=1, ADM_gall
             do kk = ADM_kmin+1,ADM_kmax
               kp=kk
               if(sigma_lev(k)>sigma(n,kk,l)) exit
             enddo
             a = sigma(n,kp,l)-sigma_lev(k)
             b = sigma_lev(k)-sigma(n,kp-1,l)
             v(n,k,l) = (b*tmp(n,kp,l) + a*tmp(n,kp-1,l))/(a+b)
          enddo
       enddo
    enddo
    if ( ADM_have_pl ) Then
       Do l=1,ADM_lall_pl
          do k=1,MAX_SIGMA
             do n=1, ADM_gall_pl
                do kk = ADM_kmin+1,ADM_kmax
                  kp=kk
                   if(sigma_lev(k)>sigma_pl(n,kk,l)) exit
                enddo
                a = sigma_pl(n,kp,l)-sigma_lev(k)
                b = sigma_lev(k)-sigma_pl(n,kp-1,l)
                v_pl(n,k,l) = (b*tmp_pl(n,kp,l) + a*tmp_pl(n,kp-1,l))/(a+b)
             enddo
          enddo
       enddo
    endif
  end subroutine VINTRPL_sigma_level
  !-----------------------------------------------------------------------------
end module mod_vintrpl
!-------------------------------------------------------------------------------
