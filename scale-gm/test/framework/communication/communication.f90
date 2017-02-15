!-------------------------------------------------------------------------------
!> Program communication
!!
!! @par Description
!!          communication kernel
!!
!! @author NICAM developers, Team SCALE
!<
!-------------------------------------------------------------------------------
program communication
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use dc_log, only: &
     LogInit
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_process, only: &
     PRC_LOCAL_MPIstart, &
     PRC_MPIfinish
  use scale_const, only: &
     CONST_setup
  use mod_adm, only: &
     ADM_setup
  use mod_comm, only: &
     COMM_setup
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scale-gm.h"
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !
  character(len=H_MID), parameter :: MODELNAME = "SCALE-GM ver. "//VERSION

  integer :: myrank
  logical :: ismaster
  !=============================================================================

  !---< MPI start >---
  !---< MPI start >---
  call PRC_LOCAL_MPIstart( myrank,  & ! [OUT]
                           ismaster ) ! [OUT]

  !########## Initial setup ##########

  ! setup standard I/O
  call IO_setup( MODELNAME, .false. )
  call IO_LOG_setup( myrank, ismaster )
  call LogInit( IO_FID_CONF, IO_FID_LOG, IO_L )

  !---< admin module setup >---
  call ADM_setup

  ! setup PROF
  call PROF_setup

  !#############################################################################
  call PROF_setprefx('INIT')
  call PROF_rapstart('Initialize',0)

  !---< cnst module setup >---
  call CONST_setup

  !---< comm module setup >---
  call COMM_setup

  call PROF_rapend('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_Communication',0)

  call communicationtest

  call PROF_rapend('Main_Communication',0)
  !#############################################################################

  call PROF_rapreport

  !--- finalize all process
  call PRC_MPIfinish

  stop
contains
  !-----------------------------------------------------------------------------
  subroutine communicationtest
    use mod_adm, only: &
       ADM_prc_me,   &
       ADM_gall,     &
       ADM_gall_pl,  &
       ADM_kall,     &
       ADM_lall,     &
       ADM_lall_pl,  &
       ADM_lall_pl,  &
       ADM_gall_1d,  &
       ADM_gmin,     &
       ADM_gmax,     &
       ADM_kmin,     &
       ADM_kmax,     &
       ADM_vlink,    &
       ADM_have_sgp, &
       ADM_prc_tab,  &
       ADM_rgn_vnum, &
       ADM_N,        &
       ADM_S
    use mod_comm, only: &
       COMM_data_transfer, &
       COMM_var
    implicit none

    real(RP) :: var   (ADM_gall   ,ADM_kall,ADM_lall   ,4)
    real(RP) :: var_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl,4)

    integer  :: i, j, k, l, ij, rgnid, prc
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ TEST start'

    var   (:,:,:,:) = -999.D0
    var_pl(:,:,:,:) = -999.D0

    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       prc   = ADM_prc_me

       do k = ADM_kmin, ADM_kmax
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = ADM_gall_1d * (j-1) + i

             var(ij,k,l,1) = real(prc,  kind=RP)
             var(ij,k,l,2) = real(rgnid,kind=RP)
             var(ij,k,l,3) = real(i,    kind=RP)
             var(ij,k,l,4) = real(j,    kind=RP)
          enddo
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then
          do k = ADM_kmin, ADM_kmax
             var(1,k,l,:) = -1.D0
          enddo
       endif
    enddo

    do l = 1, ADM_lall_pl
       rgnid = l
       prc   = ADM_prc_me

       do k  = ADM_kmin, ADM_kmax
       do ij = 1, ADM_gall_pl
          var_pl(ij,k,l,1) = real(-prc,  kind=RP)
          var_pl(ij,k,l,2) = real(-rgnid,kind=RP)
          var_pl(ij,k,l,3) = real(-ij,   kind=RP)
          var_pl(ij,k,l,4) = real(-ij,   kind=RP)
       enddo
       enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) "##### (prc,rgnid) #####"
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,1)),',',int(var(ij,k,l,2)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,1)),',',int(var_pl(ij,k,l,2)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) "##### (i,j) #####"
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,3)),',',int(var(ij,k,l,4)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,3)),',',int(var_pl(ij,k,l,4)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo



    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Communication start'

    call COMM_data_transfer( var(:,:,:,:), var_pl(:,:,:,:) )

    if( IO_L ) write(IO_FID_LOG,*) "##### (prc,rgnid) #####"
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,1)),',',int(var(ij,k,l,2)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,1)),',',int(var_pl(ij,k,l,2)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) "##### (i,j) #####"
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,3)),',',int(var(ij,k,l,4)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,3)),',',int(var_pl(ij,k,l,4)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo



    do l = 1, ADM_lall
       rgnid = ADM_prc_tab(l,ADM_prc_me)
       prc   = ADM_prc_me

       if ( ADM_rgn_vnum(ADM_N,rgnid) == ADM_vlink ) then
          do k = ADM_kmin, ADM_kmax
             j  = ADM_gmax+1
             i  = ADM_gmin
             ij = ADM_gall_1d * (j-1) + i

             var(ij,k,l,1) = real(prc,  kind=RP)
             var(ij,k,l,2) = real(rgnid,kind=RP)
             var(ij,k,l,3) = real(i,    kind=RP)
             var(ij,k,l,4) = real(j,    kind=RP)
          enddo
       endif

       if ( ADM_rgn_vnum(ADM_S,rgnid) == ADM_vlink ) then
          do k = ADM_kmin, ADM_kmax
             j  = ADM_gmin
             i  = ADM_gmax+1
             ij = ADM_gall_1d * (j-1) + i

             var(ij,k,l,1) = real(prc,  kind=RP)
             var(ij,k,l,2) = real(rgnid,kind=RP)
             var(ij,k,l,3) = real(i,    kind=RP)
             var(ij,k,l,4) = real(j,    kind=RP)
          enddo
       endif

    enddo

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ pole fill start'

    call COMM_var( var(:,:,:,:), var_pl(:,:,:,:), ADM_kall, 4 )

    if( IO_L ) write(IO_FID_LOG,*) "##### (prc,rgnid) #####"
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,1)),',',int(var(ij,k,l,2)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,1)),',',int(var_pl(ij,k,l,2)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    if( IO_L ) write(IO_FID_LOG,*) "##### (i,j) #####"
    do l  = 1, ADM_lall
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_1d
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       do j = ADM_gall_1d, 1, -1
          if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
          do i = 1, ADM_gall_1d
             ij = ADM_gall_1d * (j-1) + i
             if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                        '(',int(var(ij,k,l,3)),',',int(var(ij,k,l,4)),')'
          enddo
          if( IO_L ) write(IO_FID_LOG,*)
       enddo
    enddo
    enddo

    do l  = 1, ADM_lall_pl
    do k = ADM_kmin, ADM_kmin
       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,'(A9)',advance='no') "        |"
       do i = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(I9)',advance='no') i
       enddo
       if( IO_L ) write(IO_FID_LOG,*)

       if( IO_L ) write(IO_FID_LOG,'(I8,A1)',advance='no') j, "|"
       do ij = 1, ADM_gall_pl
          if( IO_L ) write(IO_FID_LOG,'(A1,I3,A1,I3,A1)',advance='no') &
                     '(',int(var_pl(ij,k,l,3)),',',int(var_pl(ij,k,l,4)),')'
       enddo
       if( IO_L ) write(IO_FID_LOG,*)
    enddo
    enddo

    return
  end subroutine communicationtest

end program communication
