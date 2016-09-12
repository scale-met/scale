!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics - Convert
!!
!! @par Description
!!          Convert module for Cloud Microphysics
!!          Bulk to Bin, and Bin to Bulk
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2015-01-20 (Y.Sato)  [new]
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_mp_convert
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  use scale_atmos_phy_mp_suzuki10, only: &
      nbin, &
      nccn, &
      nspc
  use scale_process, only: &
      PRC_MPIstop, &
      PRC_masterrank, &
      PRC_myrank
  use scale_comm, only: &
      COMM_datatype
  use scale_const, only: &
      PI    => CONST_PI,    &
      EPS   => CONST_EPS,   &
      UNDEF => CONST_UNDEF, &
      DWATR => CONST_DWATR, &
      DICE  => CONST_DICE
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_bulk2bin
  public :: ATMOS_PHY_MP_bin2bulk

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
  !--- Indeces
  integer, parameter :: il = 1
  integer, parameter :: ic = 2
  integer, parameter :: ip = 3
  integer, parameter :: id = 4
  integer, parameter :: iss= 5
  integer, parameter :: ig = 6
  integer, parameter :: ih = 7

contains
  !-----------------------------------------------------------------------------
  !> Bulk to Bin
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_bulk2bin(  &
       xs, xe,              & ! [IN]
       ys, ye,              & ! [IN]
       dims,                & ! [IN]
       it,                  & ! [IN]
       rank,                & ! [IN]
       handle,              & ! [IN]
       basename_org,        & ! [IN]
       dens_org,            & ! [IN]
       qtrc_org             ) ! [INOUT]
    use scale_grid_nest, only: &
       PARENT_IMAX,     &
       PARENT_JMAX
    use scale_specfunc, only: &
       SF_gamma
    use gtool_file, only: &
       FileRead
    use scale_atmos_phy_mp, only: &
       ATMOS_PHY_MP_NAME, &
       QA_MP, &
       QS_MP
    use scale_atmos_hydrometer, only: &
       I_QV
    implicit none

    integer,          intent(in) :: xs, xe
    integer,          intent(in) :: ys, ye
    integer,          intent(in) :: dims(:)
    integer,          intent(in) :: it
    integer,          intent(in) :: rank
    integer,          intent(in) :: handle
    character(LEN=*), intent(in) :: basename_org
    real(RP),         intent(in) :: dens_org(:,:,:)
    real(RP),         intent(inout) :: qtrc_org(:,:,:,:)
!    real(RP),         intent(in) :: dens_org(dims(1)+2,dims(2),dims(3))
!    real(RP),         intent(inout) :: qtrc_org(dims(1)+2,dims(2),dims(3),11)

    real(RP)         :: sigma_sdf(4)
    real(RP)         :: r0_sdf(4)
    real(RP)         :: n0_sdf(4)
    real(RP)         :: rho_sdf(4)
    character(LEN=11),parameter :: fname_micpara="micpara.dat"
    character(LEN=16) :: MP_TYPE_OUTER = "NONE"
    character(LEN=16) :: MP_TYPE_INNER = "NONE"
    integer :: COMM_world

    NAMELIST / PARAM_ATMOS_PHY_MP_BIN2BULK / &
       sigma_sdf, &
       r0_sdf,    &
       n0_sdf,    &
       rho_sdf,   &
       MP_TYPE_OUTER, &
       MP_TYPE_INNER

    real(RP), allocatable :: read3Di(:,:,:)
    real(RP), allocatable :: qc_tmp(:,:,:)
    real(RP), allocatable :: qr_tmp(:,:,:)
    real(RP), allocatable :: qi_tmp(:,:,:)
    real(RP), allocatable :: qs_tmp(:,:,:)
    real(RP), allocatable :: qg_tmp(:,:,:)

    integer  :: fid_micpara
    real(RP) :: coef0, coef1, coef2
    real(RP) :: tmp_hyd, n_hyd, lambda_hyd
    real(RP) :: dummy( nbin ), radc( nbin ) 
    integer  :: k, i, j, iq, iqa, ierr
    integer  :: nnbin, nnspc, nn
    !---------------------------------------------------------------------------

    !--- define coefficients
    coef0 = 4.0_RP/3.0_RP*PI
    coef1 = 4.0_RP/3.0_RP*sqrt(PI/2.0_RP)

    !--- determine the parameters for interpolating SDF from qxx, Nxx of parent domain
    sigma_sdf(1) = 0.2_RP
    sigma_sdf(2) = 0.35_RP
    sigma_sdf(3) = 0.35_RP
    sigma_sdf(4) = 0.35_RP
    r0_sdf(1)    = 5.E-6_RP
    r0_sdf(2)    = 2.61E-6_RP
    r0_sdf(3)    = 5.E-6_RP
    r0_sdf(4)    = 2.61E-6_RP
    n0_sdf(1)    = 8.0E+6_RP
    n0_sdf(2)    = 0.0_RP
    n0_sdf(3)    = 3.0E+6_RP
    n0_sdf(4)    = 4.0E+6_RP
    rho_sdf(1)   = DWATR
    rho_sdf(2)   = DICE
    rho_sdf(3)   = 100.0_RP
    rho_sdf(4)   = 400.0_RP

    !---- initiate
    allocate( read3Di( PARENT_IMAX(handle), PARENT_JMAX(handle), dims(1) ) )
    allocate( qc_tmp(dims(1)+2,dims(2),dims(3)) )
    allocate( qr_tmp(dims(1)+2,dims(2),dims(3)) )
    allocate( qi_tmp(dims(1)+2,dims(2),dims(3)) )
    allocate( qs_tmp(dims(1)+2,dims(2),dims(3)) )
    allocate( qg_tmp(dims(1)+2,dims(2),dims(3)) )

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_BIN2BULK,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) 'xxx Not found namelist. Check!'
       call PRC_MPIstop
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP_BIN2BULK. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP_BIN2BULK)

    if( trim(MP_TYPE_INNER) /= "SUZUKI10" ) then
       write(*,*) 'xxx MP_TYPE_INNER should be SUZUKI10 Check!'
       write(*,*) 'Now MP_TYPE_INNER set as ', MP_TYPE_INNER
       call PRC_MPIstop
    endif

    !--- read micpara.dat (microphysical parameter) and broad cast
    fid_micpara = IO_get_available_fid()
    !--- open parameter of cloud microphysics
    open ( fid_micpara, file = fname_micpara, form = 'formatted', status = 'old', iostat=ierr )

    !--- micpara.dat does not exist
    if( ierr == 0 ) then

      read( fid_micpara,* ) nnspc, nnbin

      if( nnbin /= nbin ) then
         write(*,*) 'xxx nbin in inc_tracer and nbin in micpara.dat is different check!'
         call PRC_MPIstop
      end if

      ! grid parameter
      if( IO_L ) write(IO_FID_LOG,*)  '*** Radius of cloud ****'
      do iq = 1, nbin
        read( fid_micpara,* ) nn, dummy( iq ), radc( iq )
        if( IO_L ) write(IO_FID_LOG,'(a,1x,i3,1x,a,1x,e15.7,1x,a)') &
                  "Radius of ", iq, "th cloud bin (bin center)= ", radc( iq ) , "[m]"
      end do

      close ( fid_micpara )

    else

      write(*,*) 'xxx micpara.dat does not exist. check!'
      call PRC_MPIstop

    endif

    call FileRead( read3Di(:,:,:), BASENAME_ORG, "QV", it, rank )
    do k = 1, dims(1)
       qtrc_org(k+2,xs:xe,ys:ye,I_QV) = read3Di(:,:,k)
    end do

    if( trim(MP_TYPE_OUTER) == "KESSLER" ) then

      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '+++ SDF of Bin model is created from'
      if( IO_L ) write(IO_FID_LOG,*) '+++ Kessler type Bulk microphysical model'

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QC", it, rank )
       do k = 1, dims(1)
          qc_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QR", it, rank )
       do k = 1, dims(1)
          qr_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do

       do k = 1, dims(1)
       do i = xs, xe
       do j = ys, ye

         do iq = 1, nbin
            dummy(iq) = coef1 / sigma_sdf(1) * DWATR * radc( iq )**3 &
                      * exp (  &
                            - ( log( radc(iq) )-log( r0_sdf(1) ) )**2*0.5_RP &
                            / sigma_sdf(1) / sigma_sdf(1) &
                            )
         enddo

         tmp_hyd = 0.0_RP
         do iq = 1, nbin
            tmp_hyd = tmp_hyd + dummy(iq)
         enddo

         coef2 = ( qc_tmp(k+2,i,j)+qr_tmp(k+2,i,j) ) &
               / ( tmp_hyd + ( 0.50_RP - sign(0.50_RP,tmp_hyd-EPS) ) )
         do iq = 1, nbin
           qtrc_org(k+2,i,j,QS_MP+iq) = coef2 * dummy(iq)
         enddo

       enddo
       enddo
       enddo

    elseif( trim(MP_TYPE_OUTER) == "TOMITA08" ) then

       if( IO_L ) write(IO_FID_LOG,*)
       if( IO_L ) write(IO_FID_LOG,*) '+++ SDF of Bin model is created from'
       if( IO_L ) write(IO_FID_LOG,*) '+++ TOMITA08 Bulk microphysical model'

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QC", it, rank )
       do k = 1, dims(1)
          qc_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do
!       do k = 1, dims(1)
!       do i = xs, xe
!       do j = ys, ye
!          qtrc_tmp(k+2,i,j,I_QC) = qtrc_tmp(k+2,i,j,I_QC)*dens_org(k+2,i,j)
!       enddo
!       enddo
!       enddo

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QR", it, rank )
       do k = 1, dims(1)
          qr_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QI", it, rank )
       do k = 1, dims(1)
          qi_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QS", it, rank )
       do k = 1, dims(1)
          qs_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do

       call FileRead( read3Di(:,:,:), BASENAME_ORG, "QG", it, rank )
       do k = 1, dims(1)
          qg_tmp(k+2,xs:xe,ys:ye) = read3Di(:,:,k)
       end do

       if( nspc == 1 ) then !--- put all hydrometeors to liquid (warm bin)

         do k = 1, dims(1)
         do i = xs, xe
         do j = ys, ye

           do iq = 1, nbin
              dummy(iq) = coef1 / sigma_sdf(1) * rho_sdf(1) * radc( iq )**3 &
                        * exp (  &
                              - ( log( radc(iq) )-log( r0_sdf(1) ) )**2*0.5_RP &
                              / sigma_sdf(1) / sigma_sdf(1) &
                              )
           enddo

           tmp_hyd = 0.0_RP
           do iq = 1, nbin
              tmp_hyd = tmp_hyd + dummy(iq)
           enddo

           coef2 = ( &
                     qc_tmp(k+2,i,j)+qr_tmp(k+2,i,j) &
                   + qs_tmp(k+2,i,j)+qi_tmp(k+2,i,j) &
                   + qg_tmp(k+2,i,j) &
                   ) &
                 / ( tmp_hyd + ( 0.50_RP - sign(0.50_RP,tmp_hyd-EPS) ) )
           do iq = 1, nbin
             qtrc_org(k+2,i,j,QS_MP+(il-1)*nbin+iq) = coef2 * dummy(iq)
           enddo

         enddo
         enddo
         enddo

       elseif( nspc == 7 ) then  !--- put each hydrometer to each category (ice bin)

         do k = 1, dims(1)
         do i = xs, xe
         do j = ys, ye

           !--- Rain and Cloud put into liquid bin (log-normal)
           do iq = 1, nbin
              dummy(iq) = coef1 / sigma_sdf(1) * rho_sdf(1) * radc( iq )**3 &
                        * exp (  &
                              - ( log( radc(iq) )-log( r0_sdf(1) ) )**2*0.5_RP &
                              / sigma_sdf(1) / sigma_sdf(1) &
                              )
           enddo

           tmp_hyd = 0.0_RP
           do iq = 1, nbin
              tmp_hyd = tmp_hyd + dummy(iq)
           enddo
 
           coef2 = ( qc_tmp(k+2,i,j)+qr_tmp(k+2,i,j) ) &
                 / ( tmp_hyd + ( 0.50_RP - sign(0.50_RP,tmp_hyd-EPS) ) )
           do iq = 1, nbin
             qtrc_org(k+2,i,j,QS_MP+(il-1)*nbin+iq) = coef2 * dummy(iq)
           enddo

           !--- Ice put into plate bin (log-normal)
           do iq = 1, nbin
              dummy(iq) = coef1 / sigma_sdf(2) * rho_sdf(2) * radc( iq )**3 &
                        * exp (  &
                              - ( log( radc(iq) )-log( r0_sdf(2) ) )**2*0.5_RP &
                              / sigma_sdf(2) / sigma_sdf(2) &
                              )
           enddo

           tmp_hyd = 0.0_RP
           do iq = 1, nbin
              tmp_hyd = tmp_hyd + dummy(iq)
           enddo
 
           coef2 = qi_tmp(k+2,i,j) &
                 / ( tmp_hyd + ( 0.50_RP - sign(0.50_RP,tmp_hyd-EPS) ) )
           do iq = 1, nbin
             qtrc_org(k+2,i,j,QS_MP+(ip-1)*nbin+iq) = coef2 * dummy(iq)
           enddo

           !--- Snow put into snow bin (gamma)
           n_hyd = coef0 * n0_sdf(3) * rho_sdf(3)
           lambda_hyd = ( PI * rho_sdf(3) / 6.0_RP *n0_sdf(3) * SF_gamma(4.0_RP) &
                      / ( qs_tmp(k+2,i,j) &
                        + (0.50_RP-sign(0.50_RP,qs_tmp(k+2,i,j)-EPS)) &
                        ) )**(0.25_RP)
           do iq = 1, nbin
              dummy(iq) = n_hyd * radc( iq )**3 &
                        * exp( -lambda_hyd * 0.5_RP * radc( iq ) )
           enddo
           tmp_hyd = 0.0_RP
           do iq = 1, nbin
              tmp_hyd = tmp_hyd + dummy(iq)
           enddo
 
           coef2 = qs_tmp(k+2,i,j) &
                 / ( tmp_hyd + ( 0.50_RP - sign(0.50_RP,tmp_hyd-EPS) ) )
           do iq = 1, nbin
             qtrc_org(k+2,i,j,QS_MP+(iss-1)*nbin+iq) = coef2 * dummy(iq)
           enddo

           !--- Graupel put into Graupel bin (gamma)
           n_hyd = coef0 * n0_sdf(4) * rho_sdf(4)
           lambda_hyd = ( PI * rho_sdf(4) / 6.0_RP *n0_sdf(4) * SF_gamma(4.0_RP) &
                      / ( qg_tmp(k+2,i,j) &
                        + (0.50_RP-sign(0.50_RP,qg_tmp(k+2,i,j)-EPS)) &
                        ) )**(0.25_RP)
           do iq = 1, nbin
              dummy(iq) = n_hyd * radc( iq )**3 &
                        * exp( -lambda_hyd * 0.5_RP * radc( iq ) )
           enddo
           tmp_hyd = 0.0_RP
           do iq = 1, nbin
              tmp_hyd = tmp_hyd + dummy(iq)
           enddo
 
           coef2 = qg_tmp(k+2,i,j) &
                 / ( tmp_hyd + ( 0.50_RP - sign(0.50_RP,tmp_hyd-EPS) ) )
           do iq = 1, nbin
             qtrc_org(k+2,i,j,QS_MP+(ig-1)*nbin+iq) = coef2 * dummy(iq)
           enddo

         enddo
         enddo
         enddo

       endif

    elseif( trim(MP_TYPE_OUTER) == "SN14" ) then

       write(*,*) 'SN14 is not supported for MP_TYPE_OUTER now'
       write(*,*) 'Please wait'
       call PRC_MPIstop

    elseif( trim(MP_TYPE_OUTER) == "SUZUKI10" ) then

      if( IO_L ) write(IO_FID_LOG,*)
      if( IO_L ) write(IO_FID_LOG,*) '+++ SDF of Bin model is created directory'
      if( IO_L ) write(IO_FID_LOG,*) '+++ from Bin microphysical model'
      do iq = 1, QA_MP
         iqa = QS_MP + iq - 1
         call FileRead( read3Di(:,:,:), BASENAME_ORG, ATMOS_PHY_MP_NAME(iq), it, rank )
         do k = 1, dims(1)
            qtrc_org(k+2,xs:xe,ys:ye,iqa) = read3Di(:,:,k)
         end do
      end do

    else

       write(*,*) 'MP_TYPE_OUTER should be KESSLER, TOMITA08, or SUZUKI10'
       write(*,*) 'Please check! Now MP_TYPE_OUTER set as ', MP_TYPE_OUTER
       call PRC_MPIstop

    endif

    deallocate( read3Di )
    deallocate( qc_tmp )
    deallocate( qr_tmp )
    deallocate( qi_tmp )
    deallocate( qs_tmp )
    deallocate( qg_tmp )

    return
  end subroutine ATMOS_PHY_MP_bulk2bin

  !-----------------------------------------------------------------------------
  !> Bin to Bulk
  subroutine ATMOS_PHY_MP_bin2bulk
    !---------------------------------------------------------------------------

    return
  end subroutine ATMOS_PHY_MP_bin2bulk

end module scale_atmos_phy_mp_convert
