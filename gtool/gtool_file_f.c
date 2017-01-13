#include "gtool_file.h"

static void fstr2cstr( char   *cstr, // (out)
		       char   *fstr, // (in)
		       int32_t len)
{
  int i;

  if ( cstr != fstr )
    for ( i=0; i<len; i++ )
      cstr[i] = fstr[i];

  for ( i=len-1; i>=0; i-- ) {
    if ( cstr[i] != ' ' ) {
      i += 1;
      break;
    }
  }
  cstr[i] = '\0';
}

static void cstr2fstr( char   *fstr, // (out)
		       char   *cstr, // (in)
		       int32_t len)
{
  int i;

  if ( fstr != cstr )
    for ( i=0; i<len; i++ )
      fstr[i] = cstr[i];

  for ( i=0; i<len; i++ )
    if ( cstr[i] == '\0' ) break;

  for ( ; i < len; i++ )
    fstr[i] = ' ';
}

void file_open_( int32_t *fid,       // (out)
		 char    *fname,     // (in)
		 int32_t *mode,      // (in)
		 int32_t *comm,      // (in)
		 int32_t *error,     // (out)
		 int32_t  fname_len) // (in)
{
  char _fname[File_HLONG+1];
  int32_t len;

  len = fname_len > File_HLONG ? File_HLONG : fname_len;
  fstr2cstr(_fname, fname, len);

  *error = file_open( fid, _fname, *mode, MPI_Comm_f2c(*comm) );
}

void file_set_option_( int32_t *fid,      // (in)
		       char    *filetype, // (in)
		       char    *key,      // (in)
		       char    *val,      // (in)
		       int32_t *error,    // (out)
		       int32_t filetype_len,
		       int32_t key_len,
		       int32_t val_len)
{
  char _filetype[File_HSHORT+1];
  char _key[File_HMID+1];
  char _val[File_HMID+1];
  int32_t len;

  len = filetype_len > File_HSHORT ? File_HSHORT : filetype_len;
  fstr2cstr(_filetype, filetype, len);

  len = key_len > File_HMID ? File_HMID : key_len;
  fstr2cstr(_key, key, len);

  len = val_len > File_HMID ? File_HMID : val_len;
  fstr2cstr(_val, val, len);

  *error = file_set_option(*fid, _filetype, _key, _val);
}

void file_get_datainfo_( datainfo_t *dinfo,       // (out)
			 int32_t    *fid,         // (in)
			 char       *varname,     // (in)
			 int32_t    *step,        // (in)
			 int32_t    *suppress,    // (in)
			 int32_t    *error,       // (out)
			 int32_t     varname_len) // (in)
{
  char _varname[File_HSHORT+1];
  int32_t len;
  int i;

  len = varname_len > File_HSHORT ? File_HSHORT : varname_len;
  fstr2cstr(_varname, varname, len);

  *error = file_get_datainfo( dinfo, *fid, _varname, *step, *suppress );

  cstr2fstr(dinfo->varname, dinfo->varname, File_HSHORT);
  cstr2fstr(dinfo->description, dinfo->description, File_HMID);
  cstr2fstr(dinfo->units, dinfo->units, File_HSHORT);
  cstr2fstr(dinfo->time_units, dinfo->time_units, File_HMID);
  for ( i=0; i<MAX_RANK; i++ )
    cstr2fstr(dinfo->dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT);
}
void file_read_data_( void       *var,       // (out)
		      datainfo_t *dinfo,     // (in)
		      int32_t    *precision, // (in)
		      int32_t    *error)     // (out)
{
  int i;

  fstr2cstr(dinfo->varname, dinfo->varname, File_HSHORT);
  fstr2cstr(dinfo->description, dinfo->description, File_HMID);
  fstr2cstr(dinfo->units, dinfo->units, File_HSHORT);
  fstr2cstr(dinfo->time_units, dinfo->time_units, File_HMID);
  for ( i=0; i<MAX_RANK; i++ )
    fstr2cstr(dinfo->dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT);

  *error = file_read_data( var, dinfo, *precision );
}

void file_read_data_par_( void       *var,       // (out)
		          datainfo_t *dinfo,     // (in)
		          int32_t    *ndims,     // (in)
		          int32_t    *ntypes,    // (in)
		          int32_t    *dtype,     // (in)
		          int32_t    *start,     // (in)
		          int32_t    *count,     // (in)
		          int32_t    *error)     // (out)
{
  int i;
  MPI_Offset ntypes_, start_[4], count_[4];

  fstr2cstr(dinfo->varname, dinfo->varname, File_HSHORT-1);
  fstr2cstr(dinfo->description, dinfo->description, File_HMID-1);
  fstr2cstr(dinfo->units, dinfo->units, File_HSHORT-1);
  fstr2cstr(dinfo->time_units, dinfo->time_units, File_HMID-1);
  for ( i=0; i<MAX_RANK; i++ )
    fstr2cstr(dinfo->dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT-1);

  for (i=0; i<*ndims; i++) {
      start_[i+1] = start[*ndims - i - 1] -  1;
      count_[i+1] = count[*ndims - i - 1];
  }
  ntypes_ = (MPI_Offset)(*ntypes);

  *error = file_read_data_par( var, dinfo, ntypes_, MPI_Type_f2c(*dtype), start_, count_ );
}

void file_get_global_attribute_text_( int32_t *fid,        // (in)
				      char    *key,        // (in)
				      char    *value,      // (out)
				      int32_t *error,      // (out)
				      int32_t  key_len,    // (in)
				      int32_t  value_len ) // (in)
{
  char _key[File_HLONG+1];
  char _value[File_HLONG+1];
  int32_t len;

  len = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, len);

  *error = file_get_global_attribute_text( *fid, _key, _value, value_len );

  len = value_len > File_HLONG ? File_HLONG : value_len;
  cstr2fstr(value, _value, len);
}

void file_get_global_attribute_int_( int32_t *fid,      // (in)
				     char    *key,      // (in)
				     int32_t *len,      // (in)
				     int32_t *value,    // (out)
				     int32_t *error,    // (out)
				     int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_global_attribute_int( *fid, _key, value, (size_t)*len );
}

void file_get_global_attribute_float_( int32_t *fid,      // (in)
				       char    *key,      // (in)
				       int32_t *len,      // (in)
				       float   *value,    // (out)
				       int32_t *error,    // (out)
				       int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_global_attribute_float( *fid, _key, value, (size_t)*len );
}

void file_get_global_attribute_double_( int32_t *fid,      // (in)
					char    *key,      // (in)
					int32_t *len,      // (in)
					double  *value,    // (out)
					int32_t *error,    // (out)
					int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_global_attribute_double( *fid, _key, value, (size_t)*len );
}

void file_set_global_attribute_text_( int32_t *fid,        // (in)
				      char    *key,        // (in)
				      char    *value,      // (in)
				      int32_t *error,      // (out)
				      int32_t  key_len,    // (in)
				      int32_t  value_len ) // (in)
{
  char _key[File_HLONG+1];
  char _value[File_HLONG+1];
  int32_t len;

  len = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, len);

  len = value_len > File_HLONG ? File_HLONG : value_len;
  fstr2cstr(_value, value, len);

  *error = file_set_global_attribute_text( *fid, _key, _value );
}

void file_set_global_attribute_int_( int32_t *fid,      // (in)
				     char    *key,      // (in)
				     int32_t *value,    // (in)
				     int32_t *len,      // (in)
				     int32_t *error,    // (out)
				     int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];

  key_len = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, key_len);

  *error = file_set_global_attribute_int( *fid, _key, value, (size_t)*len );
}

void file_set_global_attribute_float_( int32_t *fid,      // (in)
				       char    *key,      // (in)
				       float   *value,    // (in)
				       int32_t *len,      // (in)
				       int32_t *error,    // (out)
				       int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_global_attribute_float( *fid, _key, value, (size_t)*len );
}

void file_set_global_attribute_double_( int32_t *fid,      // (in)
					char    *key,      // (in)
					double  *value,    // (in)
					int32_t *len,      // (in)
					int32_t *error,    // (out)
					int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_global_attribute_double( *fid, _key, value, (size_t)*len );
}

void file_set_tunits_( int32_t *fid,        // (in)
		       char    *time_units, // (in)
		       int32_t *error,      // (in)
		       int32_t  len)        // (in)
{
  char _time_units[File_HMID+1];

  len = len > File_HMID ? File_HMID : len;
  fstr2cstr(_time_units, time_units, len);

  *error = file_set_tunits( *fid, _time_units );
}

void file_set_tattr_( int32_t *fid,       // (in)
		      char    *vname,     // (in)
		      char    *key,       // (in)
		      char    *val,       // (in)
		      int32_t *error,     // (out)
		      int32_t  vname_len, // (in)
		      int32_t  key_len,   // (in)
		      int32_t  val_len)   // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HSHORT+1];
  char _val[File_HLONG+1];
  int32_t len;

  len = vname_len > File_HLONG ? File_HLONG : vname_len;
  fstr2cstr(_vname, vname, len);

  len = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, len);

  len = val_len > File_HLONG ? File_HLONG : val_len;
  fstr2cstr(_val, val, len);

  *error = file_set_tattr( *fid, _vname, _key, _val );
}

void file_put_axis_( int32_t *fid,          // (in)
		     char    *name,         // (in)
		     char    *desc,         // (in)
		     char    *units,        // (in)
		     char    *dim_name,     // (in)
		     int32_t *dtype,        // (in)
		     void    *val,          // (in)
		     int32_t *size,         // (in)
		     int32_t *precision,    // (in)
		     int32_t *error,        // (out)
		     int32_t  name_len,     // (in)
		     int32_t  desc_len,     // (in)
		     int32_t  units_len,    // (in)
		     int32_t  dim_name_len) // (in)
{
  char _name[File_HSHORT+1];
  char _desc[File_HMID+1];
  char _units[File_HMID+1];
  char _dim_name[File_HSHORT+1];
  int len;

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  len = desc_len > File_HMID ? File_HMID : desc_len;
  fstr2cstr(_desc, desc, len);

  len = units_len > File_HMID ? File_HMID : units_len;
  fstr2cstr(_units, units, len);

  len = dim_name_len > File_HSHORT ? File_HSHORT : dim_name_len;
  fstr2cstr(_dim_name, dim_name, len);

  *error = file_put_axis( *fid, _name, _desc, _units, _dim_name, *dtype, val, *size, *precision );
}

void file_def_axis_( int32_t *fid,          // (in)
		     char    *name,         // (in)
		     char    *desc,         // (in)
		     char    *units,        // (in)
		     char    *dim_name,     // (in)
		     int32_t *dtype,        // (in)
		     int32_t *dim_size,     // (in)
		     int32_t *error,        // (out)
		     int32_t  name_len,     // (in)
		     int32_t  desc_len,     // (in)
		     int32_t  units_len,    // (in)
		     int32_t  dim_name_len) // (in)
{
  char _name[File_HSHORT+1];
  char _desc[File_HMID+1];
  char _units[File_HMID+1];
  char _dim_name[File_HSHORT+1];
  int len;

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  len = desc_len > File_HMID ? File_HMID : desc_len;
  fstr2cstr(_desc, desc, len);

  len = units_len > File_HMID ? File_HMID : units_len;
  fstr2cstr(_units, units, len);

  len = dim_name_len > File_HSHORT ? File_HSHORT : dim_name_len;
  fstr2cstr(_dim_name, dim_name, len);

  *error = file_def_axis( *fid, _name, _desc, _units, _dim_name, *dtype, *dim_size );
}

void file_write_axis_( int32_t *fid,          // (in)
		       char    *name,         // (in)
		       void    *val,          // (in)
		       int32_t *precision,    // (in)
		       int32_t *start,        // (in)
		       int32_t *count,        // (in)
		       int32_t *error,        // (out)
		       int32_t  name_len)     // (in)
{
  char _name[File_HSHORT+1];
  int len;
  MPI_Offset start_[1], count_[1];

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  /* all axes are 1D */
  start_[0] = *start - 1;  /* C index is 0-based */
  count_[0] = *count;

  *error = file_write_axis( *fid, _name, val, *precision, start_, count_ );
}

void file_put_associated_coordinates_( int32_t *fid,          // (in)
				       char    *name,         // (in)
				       char    *desc,         // (in)
				       char    *units,        // (in)
				       char    *dim_names,    // (in)
				       int32_t *ndims,        // (in)
				       int32_t *dtype,        // (in)
				       void    *val,          // (in)
				       int32_t *precision,    // (in)
				       int32_t *error,        // (out)
				       int32_t  name_len,     // (in)
				       int32_t  desc_len,     // (in)
				       int32_t  units_len,    // (in)
				       int32_t  dim_name_len) // (in)
{
  char _name[File_HSHORT+1];
  char _desc[File_HMID+1];
  char _units[File_HMID+1];
  char **_dim_names;
  int len;
  int i;

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  len = desc_len > File_HMID ? File_HMID : desc_len;
  fstr2cstr(_desc, desc, len);

  len = units_len > File_HMID ? File_HMID : units_len;
  fstr2cstr(_units, units, len);

  _dim_names = (char**) malloc(sizeof(char*)*(*ndims));
  len = dim_name_len > File_HSHORT ? File_HSHORT : dim_name_len;
  for ( i=0; i<*ndims; i++ ) {
    _dim_names[i] = (char*) malloc(sizeof(char)*(File_HSHORT+1));
    fstr2cstr(_dim_names[i], dim_names+i*dim_name_len, len);
  }
  *error = file_put_associated_coordinates( *fid, _name, _desc, _units, _dim_names, *ndims, *dtype, val, *precision );
}

void file_def_associated_coordinates_( int32_t *fid,          // (in)
				       char    *name,         // (in)
				       char    *desc,         // (in)
				       char    *units,        // (in)
				       char    *dim_names,    // (in)
				       int32_t *ndims,        // (in)
				       int32_t *dtype,        // (in)
				       int32_t *error,        // (out)
				       int32_t  name_len,     // (in)
				       int32_t  desc_len,     // (in)
				       int32_t  units_len,    // (in)
				       int32_t  dim_name_len) // (in)
{
  char _name[File_HSHORT+1];
  char _desc[File_HMID+1];
  char _units[File_HMID+1];
  char **_dim_names;
  int len;
  int i;

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  len = desc_len > File_HMID ? File_HMID : desc_len;
  fstr2cstr(_desc, desc, len);

  len = units_len > File_HMID ? File_HMID : units_len;
  fstr2cstr(_units, units, len);

  _dim_names = (char**) malloc(sizeof(char*)*(*ndims));
  len = dim_name_len > File_HSHORT ? File_HSHORT : dim_name_len;
  for ( i=0; i<*ndims; i++ ) {
    _dim_names[i] = (char*) malloc(sizeof(char)*(File_HSHORT+1));
    fstr2cstr(_dim_names[i], dim_names+i*dim_name_len, len);
  }
  *error = file_def_associated_coordinates( *fid, _name, _desc, _units, _dim_names, *ndims, *dtype );
}

void file_write_associated_coordinates_( int32_t *fid,          // (in)
				         char    *name,         // (in)
				         void    *val,          // (in)
				         int32_t *precision,    // (in)
				         int32_t *ndims,        // (in)
				         int32_t *start,        // (in)
				         int32_t *count,        // (in)
				         int32_t *error,        // (out)
				         int32_t  name_len)     // (in)
{
  char _name[File_HSHORT+1];
  int i, len;
  MPI_Offset start_[4], count_[4];
  /* all associated coordinates are up to 4D */

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  for (i=0; i<*ndims; i++) {
      start_[i] = start[*ndims - i - 1] -  1;
      count_[i] = count[*ndims - i - 1];
  }
  *error = file_write_associated_coordinates( *fid, _name, val, *precision, start_, count_ );
}

void file_add_variable_( int32_t  *vid,         // (out)
			 int32_t  *fid,         // (in)
			 char     *varname,     // (in)
			 char     *desc,        // (in)
			 char     *units,       // (in)
			 char     *dims,        // (in)
			 int32_t  *ndims,       // (in)
			 int32_t  *dtype,       // (in)
			 real64_t *tint,        // (in)
			 int32_t  *tavg,        // (in)
			 int32_t  *error,       // (out)
			 int32_t   varname_len, // (in)
			 int32_t   desc_len,    // (in)
			 int32_t   units_len,   // (in)
			 int32_t   dims_len)    // (in)
{
  char _varname[File_HSHORT+1];
  char _desc[File_HMID+1];
  char _units[File_HMID+1];
  char **_dims;
  int len;
  int i;

  len = varname_len > File_HSHORT ? File_HSHORT : varname_len;
  fstr2cstr(_varname, varname, len);

  len = desc_len > File_HMID ? File_HMID : desc_len;
  fstr2cstr(_desc, desc, len);

  len = units_len > File_HMID ? File_HMID : units_len;
  fstr2cstr(_units, units, len);

  _dims = (char**) malloc(sizeof(char*)*(*ndims));
  len = dims_len > File_HSHORT ? File_HSHORT : dims_len;
  for ( i=0; i<*ndims; i++ ) {
    _dims[i] = (char*) malloc(sizeof(char)*(File_HSHORT+1));
    fstr2cstr(_dims[i], dims+i*dims_len, len);
  }

  *error = file_add_variable( vid, *fid, _varname, _desc, _units, _dims, *ndims, *dtype, *tint, *tavg );

  for ( i=0; i<*ndims; i++ )
    free( _dims[i] );
  free( _dims );
}

void file_write_data_( int32_t  *fid,       // (in)
                       int32_t  *vid,       // (in)
		       void     *var,       // (in)
		       real64_t *t_start,   // (in)
		       real64_t *t_end,     // (in)
		       int32_t  *precision, // (in)
		       int32_t  *ndims,     // (in)
		       int32_t  *start,     // (in)
		       int32_t  *count,     // (in)
		       int32_t  *error)     // (out)
{
  int i;
  MPI_Offset start_[4], count_[4]; /* assume max ndims is 4 */

  for (i=0; i<*ndims; i++) {
      start_[i] = start[*ndims - i - 1] -  1;
      count_[i] = count[*ndims - i - 1];
  }
  *error = file_write_data( *fid, *vid, var, *t_start, *t_end, *precision, start_, count_ );
}

void file_close_( int32_t *fid ,   // (in)
		  int32_t *error ) // (out)
{
  *error = file_close( *fid );
}

void file_enddef_( int32_t *fid ,   // (in)
		   int32_t *error ) // (out)
{
  *error = file_enddef( *fid );
}

void file_attach_buffer_( int32_t *fid ,       // (in)
		          int32_t *buf_amount, // (out)
		          int32_t *error )     // (out)
{
  *error = file_attach_buffer( *fid, *buf_amount );
}

void file_detach_buffer_( int32_t *fid ,       // (in)
		          int32_t *error )     // (out)
{
  *error = file_detach_buffer( *fid );
}

void file_flush_( int32_t *fid ,   // (in)
		  int32_t *error ) // (out)
{
  *error = file_flush( *fid );
}
