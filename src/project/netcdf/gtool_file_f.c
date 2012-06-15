#include "gtool_file.h"

static void fstr2cstr( char   *cstr, // (out)
		       char   *fstr, // (in)
		       int32_t len)
{
  int i;

  if ( cstr != fstr ) strncpy( cstr, fstr, len );

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

  if ( fstr != cstr ) strncpy( fstr, cstr, len );

  for ( i=0; i<len; i++ )
    if ( cstr[i] == '\0' ) break;

  for ( ; i < len; i++ )
    fstr[i] = ' ';
}

void file_open_( int32_t *fid,       // (out)
		 char    *fname,     // (in)
		 int32_t *mode,      // (in)
		 int32_t *error,     // (out)
		 int32_t  fname_len) // (in)
{
  char _fname[File_HLONG+1];
  int32_t len;

  len = fname_len > File_HLONG ? File_HLONG : fname_len;
  fstr2cstr(_fname, fname, len);

  *error = file_open( fid, _fname, *mode );
}

void file_get_datainfo_( datainfo_t *dinfo,       // (out)
			 int32_t    *fid,         // (in)
			 char       *varname,     // (in)
			 int32_t    *step,        // (in)
			 int32_t    *error,       // (out)
			 int32_t     varname_len) // (in)
{
  char _varname[File_HSHORT+1];
  int32_t len;
  int i;

  len = varname_len > File_HSHORT ? File_HSHORT : varname_len;
  fstr2cstr(_varname, varname, len);

  *error = file_get_datainfo( dinfo, *fid, _varname, *step );

  for ( i=0; i<MAX_RANK; i++ )
    cstr2fstr(dinfo->dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT);
}
void file_read_data_( void       *var,       // (out)
		      datainfo_t *dinfo,     // (in)
		      int32_t    *precision, // (in)
		      int32_t    *error)     // (out)
{
  int i;

  for ( i=0; i<MAX_RANK; i++ )
    fstr2cstr(dinfo->dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT);

  *error = file_read_data( var, dinfo, *precision );
}

void file_set_global_attributes_( int32_t *fid,             // (in)
				  char    *title,           // (in)
				  char    *source,          // (in)
				  char    *institution,     // (in)
				  char    *time_units,      // (in)
				  int32_t *nodeid,          // (in)
				  int32_t *error,           // (out)
				  int32_t  title_len,       // (in)
				  int32_t  source_len,      // (in)
				  int32_t  institution_len, // (in)
				  int32_t  time_units_len)  // (in)
{
  char _title[File_HLONG+1];
  char _source[File_HLONG+1];
  char _institution[File_HLONG+1];
  char _time_units[File_HSHORT+1];
  int32_t len;

  len = title_len > File_HLONG ? File_HLONG : title_len;
  fstr2cstr(_title, title, len);

  len = source_len > File_HLONG ? File_HLONG : source_len;
  fstr2cstr(_source, source, len);

  len = institution_len > File_HLONG ? File_HLONG : institution_len;
  fstr2cstr(_institution, institution, len);

  len = time_units_len > File_HSHORT ? File_HSHORT : time_units_len;
  fstr2cstr(_time_units, time_units, len);

  *error = file_set_global_attributes( *fid, _title, _source, _institution, _time_units, *nodeid );
}

void file_set_dim_info_( int32_t *fid,           // (in)
			 int32_t *ndims,         // (in)
			 char    *dim_name,      // (in)
			 int32_t *dim_size,      // (in)
			 char    *dim_desc,      // (in)
			 char    *dim_units,     // (in)
			 int32_t *dim_type,      // (in)
			 int32_t *error,         // (out)
			 int32_t  dim_name_len,  // (in)
			 int32_t  dim_desc_len,  // (in)
			 int32_t  dim_units_len) // (in)
{
  char **_dim_name, **_dim_desc, **_dim_units;
  int32_t len;
  int i;

  _dim_name = (char**) malloc(sizeof(char*)*(*ndims));
  len = dim_name_len > File_HSHORT ? File_HSHORT : dim_name_len;
  for ( i=0; i<*ndims; i++ ) {
    _dim_name[i] = (char*) malloc(sizeof(char)*(File_HSHORT+1));
    fstr2cstr(_dim_name[i], dim_name+i*dim_name_len, len);
  }

  _dim_desc = (char**) malloc(sizeof(char*)*(*ndims));
  len = dim_desc_len > File_HMID ? File_HMID : dim_desc_len;
  for ( i=0; i<*ndims; i++ ) {
    _dim_desc[i] = (char*) malloc(sizeof(char)*(File_HMID+1));
    fstr2cstr(_dim_desc[i], dim_desc+i*dim_desc_len, len);
  }

  _dim_units = (char**) malloc(sizeof(char*)*(*ndims));
  len = dim_units_len > File_HMID ? File_HMID : dim_units_len;
  for ( i=0; i<*ndims; i++ ) {
    _dim_units[i] = (char*) malloc(sizeof(char)*(File_HMID+1));
    fstr2cstr(_dim_units[i], dim_units+i*dim_units_len, len);
  }

  *error = file_set_dim_info( *fid, *ndims, _dim_name, dim_size, _dim_desc, _dim_units, dim_type );

  for ( i=0; i<*ndims; i++ ){
    free( _dim_name[i] );
    free( _dim_desc[i] );
    free( _dim_units[i] );
  }
  free(_dim_name);
  free(_dim_desc);
  free(_dim_units);
}

void file_put_axis_( int32_t *fid,          // (in)
		     char    *dim_name,     // (in)
		     void    *val,          // (in)
		     int32_t *precision,    // (in)
		     int32_t *error,        // (out)
		     int32_t  dim_name_len) // (in)
{
  char _dim_name[File_HSHORT+1];
  int32_t len;

  len = dim_name_len > File_HSHORT ? File_HSHORT : dim_name_len;
  fstr2cstr(_dim_name, dim_name, len);

  *error = file_put_axis( *fid, _dim_name, val, *precision );
}

void file_put_additional_axis_( int32_t *fid,          // (in)
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

  *error = file_put_additional_axis( *fid, _name, _desc, _units, _dim_name, *dtype, val, *size, *precision );
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
    _dims[i] = (char*) malloc(sizeof(char)*(*ndims));
    fstr2cstr(_dims[i], dims+i*dims_len, len);
  }

  *error = file_add_variable( vid, *fid, _varname, _desc, _units, _dims, *ndims, *dtype, *tint, *tavg );

  for ( i=0; i<*ndims; i++ )
    free( _dims[i] );
  free( _dims );
}

void file_write_data_( int32_t  *vid,       // (in)
		       void     *var,       // (in)
		       real64_t *t_start,   // (in)
		       real64_t *t_end,     // (in)
		       int32_t  *precision, // (in)
		       int32_t  *error)     // (out)
{
  *error = file_write_data( *vid, var, *t_start, *t_end, *precision );
}

void file_close_( int32_t *fid ,   // (in)
		  int32_t *error ) // (out)
{
  *error = file_close( *fid );
}
