#ifndef SCALE_LOG_H
#define SCALE_LOG_H

#define LOG_ERROR(subroutine_name,format) \
  write(*,'(3A)',advance='no') "ERROR [",subroutine_name,"] "; \
  write(*,format)
#define LOG_ERROR_CONT(format) \
  write(*,'(3x)',advance='no'); \
  write(*,format)


#define LOG_WARN(subroutine_name,format) \
  if (IO_L) write(IO_FID_LOG,*); \
  if (IO_L) write(IO_FID_LOG,'(3A)',advance='no') "WARN  [",subroutine_name,"] "; \
  if (IO_L) write(IO_FID_LOG,format)


#define LOG_INFO(subroutine_name,format) \
  if (IO_L) write(IO_FID_LOG,'(3A)',advance='no') "INFO  [", subroutine_name, "] "; \
  if (IO_L) write(IO_FID_LOG,format)
#define LOG_INFO_CONT(format) \
  if (IO_L) write(IO_FID_LOG,'(3x)',advance='no'); \
  if (IO_L) write(IO_FID_LOG,format)

#define LOG_NEWLINE \
  if (IO_L) write(IO_FID_LOG,*)


#define LOG_PROGRESS(format) \
  if (IO_L) write(IO_FID_LOG,'(1x,A)',advance='no') "+++ "; \
  if (IO_L) write(IO_FID_LOG,format)


#define LOG_NML(name) \
  if (IO_NML) write(IO_FID_NML,nml=name)


#endif
