#ifndef SCALE_LOG_H
#define SCALE_LOG_H

#ifdef _OPENACC

#define LOG_ERROR(subroutine_name,format) \
  write(*,*) "ERROR  [",subroutine_name,"] universal=",IO_UNIVERSALRANK,", local=",IO_LOCALRANK,", jobID=",IO_JOBID,", domain=",IO_DOMAINID; \
  write(*,*,advance='no') "on rank ",IO_UNIVERSALRANK,": "; \
  write(*,*)

#define LOG_ERROR_CONT(format) \
  write(*,*,advance='no') "on rank ",IO_UNIVERSALRANK,": "; \
  write(*,*)


#define LOG_WARN(subroutine_name,format) \
  if (IO_L) write(*,*,advance='no') "WARN  [",subroutine_name,"]"; \
  if (IO_L) write(*,*)
#define LOG_WARN_CONT(format) \
  if (IO_L) write(*,*,advance='no');		\
  if (IO_L) write(*,*)


#else

#define LOG_ERROR(subroutine_name,format) \
  write(*,'(10A)') "ERROR  [",subroutine_name,"] universal=",IO_UNIVERSALRANK,", local=",IO_LOCALRANK,", jobID=",IO_JOBID,", domain=",IO_DOMAINID; \
  write(*,'(3A)',advance='no') "on rank ",IO_UNIVERSALRANK,": "; \
  write(*,format)

#define LOG_ERROR_CONT(format) \
  write(*,'(3A)',advance='no') "on rank ",IO_UNIVERSALRANK,": "; \
  write(*,format)


#define LOG_WARN(subroutine_name,format) \
  if (IO_L) write(*,'(3A)',advance='no') "WARN  [",subroutine_name,"]"; \
  if (IO_L) write(*,format)
#define LOG_WARN_CONT(format) \
  if (IO_L) write(*,'(5x)',advance='no'); \
  if (IO_L) write(*,format)

#endif


#define LOG_INFO(subroutine_name,format) \
  if (IO_L) write(IO_FID_LOG,'(3A)',advance='no') "INFO  [", subroutine_name, "]"; \
  if (IO_L) write(IO_FID_LOG,format)
#define LOG_INFO_CONT(format) \
  if (IO_L) write(IO_FID_LOG,'(5x)',advance='no'); \
  if (IO_L) write(IO_FID_LOG,format)
#define LOG_INFO_CONTNA(format) \
  if (IO_L) write(IO_FID_LOG,format,advance='no')


#define LOG_NEWLINE \
  if (IO_L) write(IO_FID_LOG,*)


#define LOG_PROGRESS(format) \
  if (IO_L) write(IO_FID_LOG,'(A)',advance='no') "+++++"; \
  if (IO_L) write(IO_FID_LOG,format)


#define LOG_NML(name) \
  if (IO_NML) write(IO_FID_NML,nml=name)


#endif
