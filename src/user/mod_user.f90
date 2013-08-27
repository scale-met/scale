!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include "inc_precision.h"
  include "inc_index.h"
  include 'inc_tracer.h'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_setup
  public :: USER_step
  
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
  logical, private, save :: USER_do = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine USER_setup
    use mod_stdio, only: &
       IO_FID_CONF
    use mod_process, only: &
       PRC_MPIstop
    implicit none

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[USER]/Categ[MAIN]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_USER)

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    use mod_process, only: &
       PRC_MPIstop
    use mod_atmos_vars, only: &
       DENS,    & 
       MOMX,    & 
       MOMY,    & 
       MOMZ,    & 
       RHOT,    & 
       QTRC
    use mod_atmos_vars_sf, only: &
       PREC,       &
       SWD,        &
       LWD,        &
       SFLX_MOMX,  &
       SFLX_MOMY,  &
       SFLX_MOMZ,  &
       SFLX_SWU,   &
       SFLX_LWU,   &
       SFLX_SH,    &
       SFLX_LH,    &
       SFLX_QVAtm
    use mod_land_vars, only: &
       SFLX_GH,       &
       SFLX_PREC,     &
       SFLX_QVLnd,    &
       TG,            &
       QvEfc,         &
       EMIT,          &
       ALB, TCS, DZg, &
       Z00, Z0R, Z0S, &
       Zt0, ZtR, ZtS, &
       Ze0, ZeR, ZeS
    use mod_cpl_vars, only: &
       LST
    implicit none
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       call PRC_MPIstop
    endif

    if( IO_L ) write(IO_FID_LOG,*) 'DENS       :',minval(DENS(KS,:,:)),              maxval(DENS(KS,:,:))             
    if( IO_L ) write(IO_FID_LOG,*) 'MOMX       :',minval(MOMX(KS,:,:)),              maxval(MOMX(KS,:,:))             
    if( IO_L ) write(IO_FID_LOG,*) 'MOMY       :',minval(MOMY(KS,:,:)),              maxval(MOMY(KS,:,:))             
    if( IO_L ) write(IO_FID_LOG,*) 'MOMZ       :',minval(MOMZ(KS,:,:)),              maxval(MOMZ(KS,:,:))             
    if( IO_L ) write(IO_FID_LOG,*) 'RHOT       :',minval(RHOT(KS,:,:)),              maxval(RHOT(KS,:,:))             
    if( IO_L ) write(IO_FID_LOG,*) 'QTRC       :',minval(QTRC(KS,:,:,I_QV)),         maxval(QTRC(KS,:,:,I_QV))        
    if( IO_L ) write(IO_FID_LOG,*) 'PREC       :',minval(PREC(:,:)),                 maxval(PREC(:,:))                
    if( IO_L ) write(IO_FID_LOG,*) 'SWD        :',minval(SWD (:,:)),                 maxval(SWD (:,:))                
    if( IO_L ) write(IO_FID_LOG,*) 'LWD        :',minval(LWD (:,:)),                 maxval(LWD (:,:))                
    if( IO_L ) write(IO_FID_LOG,*) 'TEMP       :',minval(RHOT(KS,:,:)/DENS(KS,:,:)), maxval(RHOT(KS,:,:)/DENS(KS,:,:))
  
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_MOMX  :',minval(SFLX_MOMX (:,:)),           maxval(SFLX_MOMX (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_MOMY  :',minval(SFLX_MOMY (:,:)),           maxval(SFLX_MOMY (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_MOMZ  :',minval(SFLX_MOMZ (:,:)),           maxval(SFLX_MOMZ (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_SWU   :',minval(SFLX_SWU  (:,:)),           maxval(SFLX_SWU  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_LWU   :',minval(SFLX_LWU  (:,:)),           maxval(SFLX_LWU  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_SH    :',minval(SFLX_SH   (:,:)),           maxval(SFLX_SH   (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_LH    :',minval(SFLX_LH   (:,:)),           maxval(SFLX_LH   (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_QVAtm :',minval(SFLX_QVAtm(:,:)),           maxval(SFLX_QVAtm(:,:))
  
    if( IO_L ) write(IO_FID_LOG,*) 'TG         :',minval(TG   (:,:)),                maxval(TG   (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'QvEfc      :',minval(QvEfc(:,:)),                maxval(QvEfc(:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'EMIT       :',minval(EMIT (:,:)),                maxval(EMIT (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'ALB        :',minval(ALB  (:,:)),                maxval(ALB  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'TCS        :',minval(TCS  (:,:)),                maxval(TCS  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'DZg        :',minval(DZg  (:,:)),                maxval(DZg  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'Z00        :',minval(Z00  (:,:)),                maxval(Z00  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'Z0R        :',minval(Z0R  (:,:)),                maxval(Z0R  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'Z0S        :',minval(Z0S  (:,:)),                maxval(Z0S  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'Zt0        :',minval(Zt0  (:,:)),                maxval(Zt0  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'ZtR        :',minval(ZtR  (:,:)),                maxval(ZtR  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'ZtS        :',minval(ZtS  (:,:)),                maxval(ZtS  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'Ze0        :',minval(Ze0  (:,:)),                maxval(Ze0  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'ZeR        :',minval(ZeR  (:,:)),                maxval(ZeR  (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'ZeS        :',minval(ZeS  (:,:)),                maxval(ZeS  (:,:))
  
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_GH    :',minval(SFLX_GH   (:,:)),           maxval(SFLX_GH   (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_PREC  :',minval(SFLX_PREC (:,:)),           maxval(SFLX_PREC (:,:))
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_QVLnd :',minval(SFLX_QVLnd(:,:)),           maxval(SFLX_QVLnd(:,:))
  
    if( IO_L ) write(IO_FID_LOG,*) 'LST        :',minval(LST(:,:)),                  maxval(LST(:,:))

    return
  end subroutine USER_step

end module mod_user
