!> struct for data infomation
type, BIND(C) :: datainfo
   character(c_char) :: varname(IO_HSHORT)   !< variable name
   character(c_char) :: description(IO_HMID) !< variable description
   character(c_char) :: unit(IO_HSHORT)      !< unit of variable
   character(c_char) :: layername(IO_HSHORT) !< layer name
   character(c_char) :: note(IO_HLONG)       !< longer note of variable
   integer(DP)       :: datasize             !< data size
   integer           :: datatype             !< data type(0,1,2,3)
   integer           :: num_of_layer         !< number of layer
   integer           :: step
   integer(DP)       :: time_start
   integer(DP)       :: time_end
endtype datainfo

!> struct for package infomation
type, BIND(C) :: headerinfo
   character(c_char) :: fname(IO_HLONG)      !< file name
   character(c_char) :: description(IO_HMID) !< variable description
   character(c_char) :: note(IO_HLONG)       !< longer note of variable
   integer           :: num_of_data          !< number of data
   integer           :: fmode                !< file mode(0,1,2)
   integer           :: endiantype           !< endian type(0,1,2)
   integer           :: grid_topology        !< grid topology(0,1,2)
   integer           :: glevel               !< glevel
   integer           :: rlevel               !< rlevel
   integer           :: num_of_rgn           !< number of region
   type(c_ptr)       :: rgnid                !< array of region id
endtype headerinfo

interface

   subroutine fio_mk_fname(fname, base, ext, i, y) BIND(C)
     import c_int, c_char
     character(c_char), intent(in)        :: fname(*)
     character(c_char), intent(in)        :: base(*)
     character(c_char), intent(in)        :: ext(*)
     integer(c_int),    intent(in), VALUE :: i
     integer(c_int),    intent(in), VALUE :: y
   end subroutine fio_mk_fname

   function fio_syscheck() result(ierr) BIND(C)
     import c_int
     integer(c_int) :: ierr
   end function fio_syscheck

   function fio_put_commoninfo(fmode, endiantype, grid_topology, glevel, rlevel, num_of_rgn, rgnid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fmode
     integer(c_int), intent(in), VALUE :: endiantype
     integer(c_int), intent(in), VALUE :: grid_topology
     integer(c_int), intent(in), VALUE :: glevel
     integer(c_int), intent(in), VALUE :: rlevel
     integer(c_int), intent(in), VALUE :: num_of_rgn
     integer(c_int), intent(in)        :: rgnid(*)
     integer(c_int) :: ierr
   end function fio_put_commoninfo

   function fio_put_commoninfo_fromfile(fid, endiantype) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in), VALUE :: endiantype
     integer(c_int) :: ierr
   end function fio_put_commoninfo_fromfile

   function fio_get_pkginfo(hinfo, fid) result(ierr) BIND(C)
     import c_int
     import headerinfo
     type(headerinfo), intent(out)     :: hinfo
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int) :: ierr
   end function fio_get_pkginfo

   function fio_get_datainfo(ditem, fid, did) result(ierr) BIND(C)
     import c_int
     import datainfo
     type(datainfo), intent(out)       :: ditem
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in), VALUE :: did
     integer(c_int) :: ierr
   end function fio_get_datainfo

   function fio_seek_datainfo(fid, varname, step) result(did) BIND(C)
     import c_int, c_char
     integer(c_int),    intent(in), VALUE :: fid
     character(c_char), intent(in)        :: varname(*)
     integer(c_int),    intent(in), VALUE :: step
     integer(c_int) :: did
   end function fio_seek_datainfo

   function fio_fopen(fid, mode) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in), VALUE :: mode
     integer(c_int) :: ierr
   end function fio_fopen

   function fio_fclose(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int) :: ierr
   end function fio_fclose

   function fio_write_pkginfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int) :: ierr
   end function fio_write_pkginfo

   function fio_read_pkginfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int) :: ierr
   end function fio_read_pkginfo

   function fio_read_data(fid, did, data) result(ierr) BIND(C)
     import c_int, c_ptr
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in), VALUE :: did
     type(c_ptr),    intent(in), VALUE :: data
     integer(c_int) :: ierr
   end function fio_read_data

   function fio_register_file(fname) result(fid) BIND(C)
     import c_int, c_char
     character(c_char), intent(in) :: fname(*)
     integer(c_int) :: fid
   end function fio_register_file

   function fio_put_write_pkginfo(fid, description, note) result(ierr) BIND(C)
     import c_int, c_char
     integer(c_int),    intent(in), VALUE :: fid
     character(c_char), intent(in)        :: description(*)
     character(c_char), intent(in)        :: note(*)
     integer(c_int) :: ierr
   end function fio_put_write_pkginfo

   function fio_valid_pkginfo_validrgn(fid, rgnid) result(ierr) BIND(C)
     import c_int
     integer(c_int),    intent(in), VALUE :: fid
     integer(c_int), intent(in)        :: rgnid(*)
     integer(c_int) :: ierr
   end function fio_valid_pkginfo_validrgn

   function fio_put_write_datainfo_data(fid, ditem, data) result(did) BIND(C)
     import c_int, c_ptr
     import datainfo
     integer(c_int), intent(in), VALUE :: fid
     type(datainfo), intent(in)        :: ditem
     type(c_ptr),    intent(in), VALUE :: data
     integer(c_int) :: did
   end function fio_put_write_datainfo_data

   function fio_read_allinfo(fid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int) :: ierr
   end function fio_read_allinfo

   function fio_read_allinfo_validrgn(fid, rgnid) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in)        :: rgnid(*)
     integer(c_int) :: ierr
   end function fio_read_allinfo_validrgn

   function fio_copy_datainfo(fid, fid_org) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in), VALUE :: fid_org
     integer(c_int) :: ierr
   end function fio_copy_datainfo

   function fio_dump_finfo(fid, endiantype, dumptype) result(ierr) BIND(C)
     import c_int
     integer(c_int), intent(in), VALUE :: fid
     integer(c_int), intent(in), VALUE :: endiantype
     integer(c_int), intent(in), VALUE :: dumptype
     integer(c_int) :: ierr
   end function fio_dump_finfo

   subroutine free(ptr) BIND(C)
     import c_ptr
     type(c_ptr), intent(in), VALUE :: ptr
   end subroutine free

end interface
