program make_include_tracer

 implicit none

  include 'setup.h'

!  integer(4) :: nbin=33, nccn=20, ICEFLG=1
  character(128) :: fname
  character(10) :: fnbin, fnccn
  character(5) :: spec(7) = (/'cloud', 'colmn', 'plate', 'dend.', 'snow ', 'grpl ', 'hail '/)
  character(2) :: tspec(7) = (/'Qc', 'Qi', 'Qp', 'Qd', 'Qs', 'Qg', 'Qh'/)
  integer(4) :: n, nn, ispc, mm, m

  if( nbin < 10 .and. nbin >= 1 ) then
   write(fnbin,'(i1)') nbin
  elseif( nbin < 100 .and. nbin >= 10 ) then
   write(fnbin,'(i2)') nbin
  elseif( nbin < 1000 .and. nbin >= 100 ) then
   write(fnbin,'(i3)') nbin
  endif

  if( nccn < 10 .and. nccn >= 1 ) then
   write(fnccn,'(i1)') nccn
  elseif( nccn < 100 .and. nccn >= 10 ) then
   write(fnccn,'(i2)') nccn
  elseif( nccn < 1000 .and. nccn >= 100 ) then
   write(fnccn,'(i3)') nccn
  endif

  if( ICEFLG == 0 ) then
   ispc = 1
   fname="inc_tracer_hbinw_h"//trim(fnbin)//"_c"//trim(fnccn)//".f90"
  elseif( ICEFLG == 1 ) then
   ispc = 7
   fname="inc_tracer_hbinf_h"//trim(fnbin)//"_c"//trim(fnccn)//".f90"
  endif
  write(*,*) fname

  open(10,file=fname, form="formatted", access="sequential")
  write(10,*)
  write(10,*) "  !------------------------------------------------------------"
  write(10,*) "  !"
  write(10,*) "  !++ scale3 grid parameter for ", trim(fname)
  write(10,*) "  !"
  write(10,*) "  !------------------------------------------------------------"
  write(10,'(a,i4)') "  integer, private, parameter :: nbin =", nbin
  write(10,'(a,i4)') "  integer, private, parameter :: nccn =", nccn
  if( ICEFLG == 1 ) then
   write(10,'(a,i4)') "  integer, private, parameter :: nspc =", nspc
  endif
  write(10,'(a,i4)') "  integer, private, parameter :: QA   =", ispc*nbin+nccn+1
  write(10,*)
  write(10,'(a)')    "  integer, private, parameter :: I_QV =  1" 
  write(10,'(a)')    "  integer, private, parameter :: I_QC =  2" 
  write(10,'(a)')    "  integer, private, parameter :: I_QR =  3" 
  write(10,'(a)')    "  integer, private, parameter :: I_QI =  4" 
  write(10,'(a)')    "  integer, private, parameter :: I_QS =  5" 
  write(10,'(a)')    "  integer, private, parameter :: I_QG =  6" 
  write(10,'(a)')    "  integer, private, parameter :: I_NC =  8" 
  write(10,'(a)')    "  integer, private, parameter :: I_NR =  9" 
  write(10,'(a)')    "  integer, private, parameter :: I_NI =  10" 
  write(10,'(a)')    "  integer, private, parameter :: I_NS =  11" 
  write(10,'(a)')    "  integer, private, parameter :: I_NG =  12" 
  write(10,*)
  write(10,'(a,i4)') "  integer, private, parameter :: QQA  =", ispc*nbin+1
  write(10,'(a)')    "  integer, private, parameter :: QQS  = 1"
  write(10,'(a,i4)') "  integer, private, parameter :: QQE  =", ispc*nbin+1
  write(10,*)
  write(10,'(a)')    "  integer, private, parameter :: QWS  = 2"
  write(10,'(a,i4)') "  integer, private, parameter :: QWE  =", 1*nbin+1
  if( ICEFLG == 1 ) then
   write(10,'(a,i4)') "  integer, private, parameter :: QIS  =", 1*nbin+1+1
   write(10,'(a,i4)') "  integer, private, parameter :: QIE  =", ispc*nbin+1
  else
   write(10,'(a)')    "  integer, private, parameter :: QIS  = 0"
   write(10,'(a)')    "  integer, private, parameter :: QIE  = 0"
  endif
  write(10,*)
  write(10,'(a)')    "  character(len=16), private, save :: AQ_NAME(QA)"
  write(10,'(a)')    "  character(len=64), private, save :: AQ_DESC(QA)"
  write(10,'(a)')    "  character(len=16), private, save :: AQ_UNIT(QA)"
  write(10,*)

  write(10,'(a)')    "  data AQ_NAME(  1)  / 'QV' /"
  m = 1
  mm = 0
  do n = 1, ispc*nbin  
   mm=mm+1
   m=m+1
   nn = (n-1)/nbin+1
   if( mm<9 ) then
    write(10,'(a15,i3,a6,a2,i1,a)')    "  data AQ_NAME(", m, ")  / '", tspec(nn), mm, "' /"
   elseif( mm==9 ) then
    write(10,'(a15,i3,a6,a2,i1,a)')    "  data AQ_NAME(", m, ")  / '", tspec(nn), mm, "' /"
   elseif( mm<99 ) then
    write(10,'(a15,i3,a6,a2,i2,a)')    "  data AQ_NAME(", m, ")  / '", tspec(nn), mm, "' /"
   elseif( mm==99 ) then
    write(10,'(a15,i3,a6,a2,i2,a)')    "  data AQ_NAME(", m, ")  / '", tspec(nn), mm, "' /"
   elseif( mm<=999 ) then
    write(10,'(a15,i3,a6,a2,i3,a)')    "  data AQ_NAME(", m, ")  / '", tspec(nn), mm, "' /"
   endif
   if( mm == nbin ) mm = 0
  enddo

  do n = ispc*nbin+1, ispc*nbin+nccn
   m = m+1
   if( n-ispc*nbin-1<9 ) then
    write(10,'(a15,i3,a,i1,a)')    "  data AQ_NAME(", m, ")  / 'Qa", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1==9 ) then
    write(10,'(a15,i3,a,i2,a)')    "  data AQ_NAME(", m, ")  / 'Qa", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1<99 ) then
    write(10,'(a15,i3,a,i2,a)')    "  data AQ_NAME(", m, ")  / 'Qa", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1==99 ) then
    write(10,'(a15,i3,a,i3,a)')    "  data AQ_NAME(", m, ")  / 'Qa", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1<999 ) then
    write(10,'(a15,i3,a,i3,a)')    "  data AQ_NAME(", m, ")  / 'Qa", n-ispc*nbin, "' /"
   endif
  enddo
  write(10,*)

  write(10,'(a)')    "  data AQ_DESC(  1)  / 'Water Vapor mixing ratio' /"
  m = 1
  mm = 0
  do n = 1, ispc*nbin  
   nn = (n-1)/nbin+1
   m = m+1
   mm = mm+1
    if( mm<9 ) then
     write(10,'(a15,i3,a22,a5,a4,i1,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of ", spec(nn), " bin", mm, "' /"
    elseif( mm==9 ) then
     write(10,'(a15,i3,a22,a5,a4,i1,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of ", spec(nn), " bin", mm, "' /"
    elseif( mm<99 ) then
     write(10,'(a15,i3,a22,a5,a4,i2,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of ", spec(nn), " bin", mm, "' /"
    elseif( mm==99 ) then
     write(10,'(a15,i3,a22,a5,a4,i2,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of ", spec(nn), " bin", mm, "' /"
    elseif( mm<999 ) then
     write(10,'(a15,i3,a22,a5,a4,i3,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of ", spec(nn), " bin", mm, "' /"
    endif
   if( mm == nbin ) mm = 0
  enddo
  do n = ispc*nbin+1, ispc*nbin+nccn
   m = m+1
   if( n-ispc*nbin-1<8 ) then
    write(10,'(a15,i3,a33,i1,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of aerosol bin", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1==8 ) then
    write(10,'(a15,i3,a33,i1,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of aerosol bin", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1<98 ) then
    write(10,'(a15,i3,a33,i2,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of aerosol bin", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1==98 ) then
    write(10,'(a15,i3,a33,i2,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of aerosol bin", n-ispc*nbin, "' /"
   elseif( n-ispc*nbin-1<998 ) then
    write(10,'(a15,i3,a33,i3,a)')    "  data AQ_DESC(", m, ")  / 'Mixing ratio of aerosol bin", n-ispc*nbin, "' /"
   endif
  enddo
  write(10,*)

  write(10,'(a)')    "  data AQ_UNIT(  1)  / 'kg/kg' /"
  m = 1
  do n = 1, ispc*nbin+nccn
   m = m+1
   if( n<9 ) then
    write(10,'(a15,i3,a24,i1,a)')    "  data AQ_UNIT(", m, ")  / 'kg/kg/unit logr' /"
   elseif( n==9 ) then
    write(10,'(a15,i3,a24,i1,a)')    "  data AQ_UNIT(", m, ")  / 'kg/kg/unit logr' /"
   elseif( n<99 ) then
    write(10,'(a15,i3,a24,i2,a)')    "  data AQ_UNIT(", m, ")  / 'kg/kg/unit logr' /"
   elseif( n==99 ) then
    write(10,'(a15,i3,a24,i2,a)')    "  data AQ_UNIT(", m, ")  / 'kg/kg/unit logr' /"
   elseif( n<999 ) then
    write(10,'(a15,i3,a24,i3,a)')    "  data AQ_UNIT(", m, ")  / 'kg/kg/unit logr' /"
   endif
  enddo


  close(10)

endprogram
