#ifndef SCALE_LOG_H
#define SCALE_LOG_H

#define LOG_ERROR(subroutine_name,format) \
  write(*,'(5A)',advance='no') "ERROR on rank ",IO_WHOAMI," [",subroutine_name,"]"; \
  write(*,format)
#define LOG_ERROR_CONT(format) \
  write(*,'(6x,3A)',advance='no') "on rank ",IO_WHOAMI,": "; \
  write(*,format)


#define LOG_WARN(subroutine_name,format) \
  if (IO_L) write(IO_FID_LOG,'(3A)',advance='no') "WARN  [",subroutine_name,"]"; \
  if (IO_L) write(IO_FID_LOG,format)
#define LOG_WARN_CONT(format) \
  if (IO_L) write(IO_FID_LOG,'(5x)',advance='no'); \
  if (IO_L) write(IO_FID_LOG,format)


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
