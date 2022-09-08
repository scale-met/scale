  interface

     function file_open_c(fid, fname, mode, comm) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),           intent(out) :: fid
       character(c_char),        intent(in)  :: fname(*)
       integer(c_int),    VALUE, intent(in)  :: mode
       integer(c_int),    VALUE, intent(in)  :: comm
       integer(c_int)                        :: error
     end function file_open_c

     function file_get_dim_length_c(len, fid, dimname, suppress) result(error) BIND(C)
       import c_int, c_char, c_bool
       integer(c_int),           intent(out) :: len
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: dimname(*)
       logical(c_bool),   VALUE, intent(in)  :: suppress
       integer(c_int)                        :: error
     end function file_get_dim_length_c

     function file_set_option_c(fid, filetype, key, val) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: filetype(*)
       character(c_char),        intent(in)  :: key(*)
       character(c_char),        intent(in)  :: val(*)
       integer(c_int)                        :: error
     end function file_set_option_c

     function file_get_nvars_c(nvars, fid) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),           intent(out) :: nvars
       integer(c_int),    VALUE, intent(in)  :: fid
       integer(c_int)                        :: error
     end function file_get_nvars_c

     function file_get_varname_c(name, fid, vid, len) result(error) BIND(C)
       import c_int, c_char
       character(c_char),        intent(out) :: name(*)
       integer(c_int),    VALUE, intent(in)  :: fid
       integer(c_int),    VALUE, intent(in)  :: vid
       integer(c_int),    VALUE, intent(in)  :: len
       integer(c_int)                        :: error
     end function file_get_varname_c

     function file_get_datainfo_c(dinfo, fid, varname, step, suppress) result(error) BIND(C)
       import c_int, c_char, c_bool
       import datainfo
       type(datainfo),           intent(out) :: dinfo
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: varname(*)
       integer(c_int),    VALUE, intent(in)  :: step
       logical(c_bool),   VALUE, intent(in)  :: suppress
       integer(c_int)                        :: error
     end function file_get_datainfo_c

     function file_get_step_size_c(len, fid, varname) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),           intent(out) :: len
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: varname(*)
       integer(c_int)                        :: error
     end function file_get_step_size_c

     function file_read_data_c(var, dinfo, precision, ntypes, dtype, start, count) result(error) BIND(C)
       import c_int, c_char, c_float, c_ptr
       import datainfo
       type(c_ptr),       VALUE, intent(in)  :: var
       type(datainfo),           intent(in)  :: dinfo
       integer(c_int),    VALUE, intent(in)  :: precision
       integer(c_int),    VALUE, intent(in)  :: ntypes
       integer(c_int),    VALUE, intent(in)  :: dtype
       integer(c_int),           intent(in)  :: start(*)
       integer(c_int),           intent(in)  :: count(*)
       integer(c_int)                        :: error
     end function file_read_data_c

     function file_get_attribute_text_c(value, fid, vname, key, suppress, len) result(error) BIND(C)
       import c_int, c_char, c_bool
       character(c_char),        intent(out) :: value(*)
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: vname(*)
       character(c_char),        intent(in)  :: key(*)
       logical(c_bool),   VALUE, intent(in)  :: suppress
       integer(c_int),    VALUE, intent(in)  :: len
       integer(c_int)                        :: error
     end function file_get_attribute_text_c

     function file_get_attribute_int_c(value, fid, vname, key, suppress, len) result(error) BIND(C)
       import c_int, c_char, c_bool
       integer(c_int),           intent(out) :: value(*)
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: vname(*)
       character(c_char),        intent(in)  :: key(*)
       logical(c_bool),   VALUE, intent(in)  :: suppress
       integer(c_int),    VALUE, intent(in)  :: len
       integer(c_int)                        :: error
     end function file_get_attribute_int_c

     function file_get_attribute_float_c(value, fid, vname, key, suppress, len) result(error) BIND(C)
       import c_int, c_char, c_bool, c_float
       real(c_float),            intent(out) :: value(*)
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: vname(*)
       character(c_char),        intent(in)  :: key(*)
       logical(c_bool),   VALUE, intent(in)  :: suppress
       integer(c_int),    VALUE, intent(in)  :: len
       integer(c_int)                        :: error
     end function file_get_attribute_float_c

     function file_get_attribute_double_c(value, fid, vname, key, suppress, len) result(error) BIND(C)
       import c_int, c_char, c_bool, c_double
       real(c_double),           intent(out) :: value(*)
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: vname(*)
       character(c_char),        intent(in)  :: key(*)
       logical(c_bool),   VALUE, intent(in)  :: suppress
       integer(c_int),    VALUE, intent(in)  :: len
       integer(c_int)                        :: error
     end function file_get_attribute_double_c

     function file_set_attribute_text_c(fid, vname, key, value) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: vname(*)
       character(c_char),        intent(in) :: key(*)
       character(c_char),        intent(in) :: value(*)
       integer(c_int)                       :: error
     end function file_set_attribute_text_c

     function file_set_attribute_int_c(fid, vname, key, value, len) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: vname(*)
       character(c_char),        intent(in) :: key(*)
       integer(c_int),           intent(in) :: value(*)
       integer(c_int),    VALUE, intent(in) :: len
       integer(c_int)                       :: error
     end function file_set_attribute_int_c

     function file_set_attribute_float_c(fid, vname, key, value, len) result(error) BIND(C)
       import c_int, c_char, c_float
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: vname(*)
       character(c_char),        intent(in) :: key(*)
       real(c_float),            intent(in) :: value(*)
       integer(c_int),    VALUE, intent(in) :: len
       integer(c_int)                       :: error
     end function file_set_attribute_float_c

     function file_set_attribute_double_c(fid, vname, key, value, len) result(error) BIND(C)
       import c_int, c_char, c_double
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: vname(*)
       character(c_char),        intent(in) :: key(*)
       real(c_double),           intent(in) :: value(*)
       integer(c_int),    VALUE, intent(in) :: len
       integer(c_int)                       :: error
     end function file_set_attribute_double_c

     function file_add_associatedvariable_c(fid, vname) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: vname(*)
       integer(c_int)                       :: error
     end function file_add_associatedvariable_c

     function file_set_tunits_c(fid, time_units, calendar) result(error) BIND(C)
       import c_int, c_char
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: time_units(*)
       character(c_char),        intent(in) :: calendar(*)
       integer(c_int)                       :: error
     end function file_set_tunits_c

     function file_put_axis_c(fid, name, desc, units, dim_name, dtype, val, size, precision) result(error) BIND(C)
       import c_int, c_char, c_float, c_ptr
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: name(*)
       character(c_char),        intent(in) :: desc(*)
       character(c_char),        intent(in) :: units(*)
       character(c_char),        intent(in) :: dim_name(*)
       integer(c_int),    VALUE, intent(in) :: dtype
       type(c_ptr),       VALUE, intent(in) :: val
       integer(c_int),    VALUE, intent(in) :: size
       integer(c_int),    VALUE, intent(in) :: precision
       integer(c_int)                       :: error
     end function file_put_axis_c

     function file_def_axis_c(fid, name, desc, units, dim_name, dtype, dim_size, bounds) result(error) BIND(C)
       import c_int, c_char, c_double
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: name(*)
       character(c_char),        intent(in) :: desc(*)
       character(c_char),        intent(in) :: units(*)
       character(c_char),        intent(in) :: dim_name(*)
       integer(c_int),    VALUE, intent(in) :: dtype
       integer(c_int),    VALUE, intent(in) :: dim_size
       integer(c_int),    VALUE, intent(in) :: bounds
       integer(c_int)                       :: error
     end function file_def_axis_c

     function file_write_axis_c(fid, name, val, precision, start, count) result(error) BIND(C)
       import c_int, c_char, c_float, c_ptr
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: name(*)
       type(c_ptr),       VALUE, intent(in) :: val
       integer(c_int),    VALUE, intent(in) :: precision
       integer(c_int),           intent(in) :: start(*)
       integer(c_int),           intent(in) :: count(*)
       integer(c_int)                       :: error
     end function file_write_axis_c

     function file_put_associatedcoordinate_c(fid, name, desc, units, dim_names, ndims, dtype, val, precision) result(error) BIND(C)
       import c_int, c_char, c_float, c_ptr
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: name(*)
       character(c_char),        intent(in) :: desc(*)
       character(c_char),        intent(in) :: units(*)
       type(c_ptr),              intent(in) :: dim_names(*)
       integer(c_int),    VALUE, intent(in) :: ndims
       integer(c_int),    VALUE, intent(in) :: dtype
       type(c_ptr),       VALUE, intent(in) :: val
       integer(c_int),    VALUE, intent(in) :: precision
       integer(c_int)                       :: error
     end function file_put_associatedcoordinate_c

     function file_def_associatedcoordinate_c(fid, name, desc, units, dim_names, ndims, dtype) result(error) BIND(C)
       import c_int, c_char, c_ptr
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: name(*)
       character(c_char),        intent(in) :: desc(*)
       character(c_char),        intent(in) :: units(*)
       type(c_ptr),              intent(in) :: dim_names(*)
       integer(c_int),    VALUE, intent(in) :: ndims
       integer(c_int),    VALUE, intent(in) :: dtype
       integer(c_int)                       :: error
     end function file_def_associatedcoordinate_c

     function file_write_associatedcoordinate_c(fid, name, val, ndims, precision, start, count) result(error) BIND(C)
       import c_int, c_char, c_float, c_ptr
       integer(c_int),    VALUE, intent(in) :: fid
       character(c_char),        intent(in) :: name(*)
       type(c_ptr),       VALUE, intent(in) :: val
       integer(c_int),    VALUE, intent(in) :: ndims
       integer(c_int),    VALUE, intent(in) :: precision
       integer(c_int),           intent(in) :: start(*)
       integer(c_int),           intent(in) :: count(*)
       integer(c_int)                       :: error
     end function file_write_associatedcoordinate_c

     function file_add_variable_c(vid, fid, varname, desc, units, stdname, dims, ndims, dtype, tint, tave) result(error) BIND(C)
       import c_int, c_char, c_ptr, c_double, c_bool
       integer(c_int),           intent(out) :: vid
       integer(c_int),    VALUE, intent(in)  :: fid
       character(c_char),        intent(in)  :: varname(*)
       character(c_char),        intent(in)  :: desc(*)
       character(c_char),        intent(in)  :: units(*)
       character(c_char),        intent(in)  :: stdname(*)
       type(c_ptr),              intent(in)  :: dims(*)
       integer(c_int),    VALUE, intent(in)  :: ndims
       integer(c_int),    VALUE, intent(in)  :: dtype
       real(c_double),    VALUE, intent(in)  :: tint
       logical(c_bool),   VALUE, intent(in)  :: tave
       integer(c_int)                        :: error
     end function file_add_variable_c

     function file_write_data_c(fid, vid, var, t_start, t_end, ndims, precision, start, count) result(error) BIND(C)
       import c_int, c_float, c_double, c_ptr
       integer(c_int),    VALUE, intent(in) :: fid
       integer(c_int),    VALUE, intent(in) :: vid
       type(c_ptr),       VALUE, intent(in) :: var
       real(c_double),    VALUE, intent(in) :: t_start
       real(c_double),    VALUE, intent(in) :: t_end
       integer(c_int),    VALUE, intent(in) :: ndims
       integer(c_int),    VALUE, intent(in) :: precision
       integer(c_int),           intent(in) :: start(*)
       integer(c_int),           intent(in) :: count(*)
       integer(c_int)                       :: error
     end function file_write_data_c

     function file_enddef_c(fid) result(error) BIND(C)
       import c_int
       integer(c_int), VALUE, intent(in) :: fid
       integer(c_int)                    :: error
     end function file_enddef_c

     function file_redef_c(fid) result(error) BIND(C)
       import c_int
       integer(c_int), VALUE, intent(in) :: fid
       integer(c_int)                    :: error
     end function file_redef_c

     function file_attach_buffer_c(fid, buf_amount) result(error) BIND(C)
       import c_int, c_int64_t
       integer(c_int),     VALUE, intent(in) :: fid
       integer(c_int64_t), VALUE, intent(in) :: buf_amount
       integer(c_int)                       :: error
     end function file_attach_buffer_c

     function file_detach_buffer_c(fid) result(error) BIND(C)
       import c_int
       integer(c_int), VALUE, intent(in) :: fid
       integer(c_int)                    :: error
     end function file_detach_buffer_c

     function file_flush_c(fid) result(error) BIND(C)
       import c_int
       integer(c_int), VALUE, intent(in) :: fid
       integer(c_int)                    :: error
     end function file_flush_c

     function file_close_c(fid, abort) result(error) BIND(C)
       import c_int, c_bool
       integer(c_int),  VALUE, intent(in) :: fid
       logical(c_bool), VALUE, intent(in) :: abort
       integer(c_int)                     :: error
     end function file_close_c

  end interface
