module mod_temp
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  implicit none
  include 'inc_precision.h'
  include 'inc_index.h'
  include 'inc_tracer.h'

contains
  subroutine TEMPLOG
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
  
    if( IO_L ) write(IO_FID_LOG,*) 'DENS      :',DENS(KS,IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'MOMX      :',MOMX(KS,IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'MOMY      :',MOMY(KS,IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'MOMZ      :',MOMZ(KS,IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'RHOT      :',RHOT(KS,IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'QTRC      :',QTRC(KS,IS,JS,I_QV)
    if( IO_L ) write(IO_FID_LOG,*) 'PREC      :',PREC(IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SWD       :',SWD (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'LWD       :',LWD (IS,JS)
  
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_MOMX :',SFLX_MOMX (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_MOMY :',SFLX_MOMY (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_MOMZ :',SFLX_MOMZ (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_SWU  :',SFLX_SWU  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_LWU  :',SFLX_LWU  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_SH   :',SFLX_SH   (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_LH   :',SFLX_LH   (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_QVAtm:',SFLX_QVAtm(IS,JS)
  
    if( IO_L ) write(IO_FID_LOG,*) 'TG        :',TG   (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'QvEfc     :',QvEfc(IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'EMIT      :',EMIT (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'ALB       :',ALB  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'TCS       :',TCS  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'DZg       :',DZg  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'Z00       :',Z00  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'Z0R       :',Z0R  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'Z0S       :',Z0S  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'Zt0       :',Zt0  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'ZtR       :',ZtR  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'ZtS       :',ZtS  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'Ze0       :',Ze0  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'ZeR       :',ZeR  (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'ZeS       :',ZeS  (IS,JS)
  
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_GH   :',SFLX_GH   (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_PREC :',SFLX_PREC (IS,JS)
    if( IO_L ) write(IO_FID_LOG,*) 'SFLX_QVLnd:',SFLX_QVLnd(IS,JS)
  
    if( IO_L ) write(IO_FID_LOG,*) 'LST       :',LST(IS,JS)
  
    return
  end subroutine TEMPLOG
end module mod_temp
