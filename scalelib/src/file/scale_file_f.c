#include "scale_file.h"

static void fstr2cstr(       char   *cstr, // (out)
                       const char   *fstr, // (in)
		       const int32_t len)
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

static void cstr2fstr(       char   *fstr, // (out)
                       const char   *cstr, // (in)
                       const int32_t len)
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

void file_open_c_(       int32_t *fid,       // (out)
                   const char    *fname,     // (in)
                   const int32_t *mode,      // (in)
                   const int32_t *comm,      // (in)
                         int32_t *error,     // (out)
                   const int32_t  fname_len) // (in)
{
  char _fname[File_HLONG+1];
  int32_t len;

  len = fname_len > File_HLONG ? File_HLONG : fname_len;
  fstr2cstr(_fname, fname, len);

  *error = file_open_c( fid, _fname, *mode, MPI_Comm_f2c(*comm) );
}

void file_set_option_c_( const int32_t *fid,         // (in)
                         const char    *filetype,    // (in)
                         const char    *key,         // (in)
                         const char    *val,         // (in)
                               int32_t *error,       // (out)
                         const int32_t filetype_len, // (in)
                         const int32_t key_len,      // (in)
                         const int32_t val_len)      // (in)
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

  *error = file_set_option_c(*fid, _filetype, _key, _val);
}

void file_get_nvars_c_( const int32_t *fid,     // (in)
                              int32_t *nvars,   // (out)
                              int32_t *error  ) // (out)
{
  *error = file_get_nvars_c( *fid, nvars );
}

void file_get_varname_c_( const int32_t *fid,          // (in)
			  const int32_t *vid,          // (in)
			        char    *varname,      // (out)
			        int32_t *error,        // (out)
			  const int32_t  varname_len ) // (in)
{
  char _varname[File_HSHORT+1];
  int32_t len;

  *error = file_get_varname_c( *fid, *vid, _varname, varname_len );

  len = varname_len > File_HSHORT ? File_HSHORT : varname_len;
  cstr2fstr(varname, _varname, len);
}

void file_get_datainfo_c_(       datainfo_t *dinfo,       // (out)
                           const int32_t    *fid,         // (in)
                           const char       *varname,     // (in)
			   const int32_t    *step,        // (in)
			   const int32_t    *suppress,    // (in)
			         int32_t    *error,       // (out)
			   const int32_t     varname_len) // (in)
{
  char _varname[File_HSHORT+1];
  int32_t len;
  int i;

  len = varname_len > File_HSHORT ? File_HSHORT : varname_len;
  fstr2cstr(_varname, varname, len);

  *error = file_get_datainfo_c( dinfo, *fid, _varname, *step, *suppress );

  cstr2fstr(dinfo->varname,     dinfo->varname,     File_HSHORT);
  cstr2fstr(dinfo->description, dinfo->description, File_HMID);
  cstr2fstr(dinfo->units,       dinfo->units,       File_HSHORT);
  cstr2fstr(dinfo->time_units,  dinfo->time_units,  File_HMID);
  for ( i=0; i<RANK_MAX; i++ )
    cstr2fstr(dinfo->dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT);
}
void file_read_data_c_(       void       *var,       // (out)
			const datainfo_t *dinfo,     // (in)
			const int32_t    *precision, // (in)
			      int32_t    *error)     // (out)
{
  datainfo_t cdinfo;
  int i;

  cdinfo.datatype = dinfo->datatype;
  cdinfo.rank = dinfo->rank;
  for ( i=0; i<dinfo->rank; i++ ) cdinfo.dim_size[i] = dinfo->dim_size[i];
  cdinfo.step = dinfo->step;
  cdinfo.time_start = dinfo->time_start;
  cdinfo.time_end = dinfo->time_end;
  cdinfo.fid = dinfo->fid;
  fstr2cstr(cdinfo.varname, dinfo->varname, File_HSHORT-1);
  fstr2cstr(cdinfo.description, dinfo->description, File_HMID-1);
  fstr2cstr(cdinfo.units, dinfo->units, File_HSHORT-1);
  fstr2cstr(cdinfo.time_units, dinfo->time_units, File_HMID-1);
  for ( i=0; i<RANK_MAX; i++ )
    fstr2cstr(cdinfo.dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT-1);

  *error = file_read_data_c( var, &cdinfo, *precision );
}

void file_read_data_par_c_(       void       *var,    // (out)
                            const datainfo_t *dinfo,  // (in)
			    const int32_t    *ndims,  // (in)
                            const int32_t    *ntypes, // (in)
                            const int32_t    *dtype,  // (in)
                            const int32_t    *start,  // (in)
                            const int32_t    *count,  // (in)
				  int32_t    *error ) // (out)
{
  int i;
  MPI_Offset ntypes_, start_[4], count_[4];
  datainfo_t cdinfo;

  cdinfo.datatype = dinfo->datatype;
  cdinfo.rank = dinfo->rank;
  for ( i=0; i<dinfo->rank; i++ ) cdinfo.dim_size[i] = dinfo->dim_size[i];
  cdinfo.step = dinfo->step;
  cdinfo.time_start = dinfo->time_start;
  cdinfo.time_end = dinfo->time_end;
  cdinfo.fid = dinfo->fid;
  fstr2cstr(cdinfo.varname, dinfo->varname, File_HSHORT-1);
  fstr2cstr(cdinfo.description, dinfo->description, File_HMID-1);
  fstr2cstr(cdinfo.units, dinfo->units, File_HSHORT-1);
  fstr2cstr(cdinfo.time_units, dinfo->time_units, File_HMID-1);
  for ( i=0; i<RANK_MAX; i++ )
    fstr2cstr(cdinfo.dim_name+i*File_HSHORT, dinfo->dim_name+i*File_HSHORT, File_HSHORT-1);

  for (i=0; i<*ndims; i++) {
      start_[i+1] = start[*ndims - i - 1] -  1;
      count_[i+1] = count[*ndims - i - 1];
  }
  ntypes_ = (MPI_Offset)(*ntypes);

  *error = file_read_data_par_c( var, &cdinfo, ntypes_, MPI_Type_f2c(*dtype), start_, count_ );
}

void file_get_globalattribute_text_c_( const int32_t *fid,        // (in)
				       const char    *key,        // (in)
					     char    *value,      // (out)
				             int32_t *error,      // (out)
				       const int32_t  key_len,    // (in)
				       const int32_t  value_len ) // (in)
{
  char _key[File_HLONG+1];
  char _value[File_HLONG+1];
  int32_t len;

  len = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, len);

  *error = file_get_globalattribute_text_c( *fid, _key, _value, value_len );

  len = value_len > File_HLONG ? File_HLONG : value_len;
  cstr2fstr(value, _value, len);
}

void file_get_globalattribute_int_c_( const int32_t *fid,      // (in)
				      const char    *key,      // (in)
				      const int32_t *len,      // (in)
				            int32_t *value,    // (out)
				            int32_t *error,    // (out)
				      const int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_globalattribute_int_c( *fid, _key, value, (size_t)*len );
}

void file_get_globalattribute_float_c_( const int32_t *fid,      // (in)
					const char    *key,      // (in)
					const int32_t *len,      // (in)
					      float   *value,    // (out)
					      int32_t *error,    // (out)
					const int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_globalattribute_float_c( *fid, _key, value, (size_t)*len );
}

void file_get_globalattribute_double_c_( const int32_t *fid,      // (in)
					 const char    *key,      // (in)
					 const int32_t *len,      // (in)
					       double  *value,    // (out)
					       int32_t *error,    // (out)
					 const int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_globalattribute_double_c( *fid, _key, value, (size_t)*len );
}

void file_set_globalattribute_text_c_( const int32_t *fid,        // (in)
				       const char    *key,        // (in)
				       const char    *value,      // (in)
				             int32_t *error,      // (out)
				       const int32_t  key_len,    // (in)
				       const int32_t  value_len ) // (in)
{
  char _key[File_HLONG+1];
  char _value[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  l = value_len > File_HLONG ? File_HLONG : value_len;
  fstr2cstr(_value, value, l);

  *error = file_set_globalattribute_text_c( *fid, _key, _value );
}

void file_set_globalattribute_int_c_( const int32_t *fid,      // (in)
				      const char    *key,      // (in)
				      const int32_t *value,    // (in)
				      const int32_t *len,      // (in)
				            int32_t *error,    // (out)
				      const int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_globalattribute_int_c( *fid, _key, value, (size_t)*len );
}

void file_set_globalattribute_float_c_( const int32_t *fid,      // (in)
					const char    *key,      // (in)
					const float   *value,    // (in)
					const int32_t *len,      // (in)
					      int32_t *error,    // (out)
					const int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_globalattribute_float_c( *fid, _key, value, (size_t)*len );
}

void file_set_globalattribute_double_c_( const int32_t *fid,      // (in)
					 const char    *key,      // (in)
					 const double  *value,    // (in)
					 const int32_t *len,      // (in)
					       int32_t *error,    // (out)
					 const int32_t  key_len ) // (in)
{
  char _key[File_HLONG+1];
  int32_t l;

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_globalattribute_double_c( *fid, _key, value, (size_t)*len );
}

void file_get_attribute_text_c_( const int32_t *fid,        // (in)
				 const char    *vname,      // (in)
				 const char    *key,        // (in)
     				       char    *value,      // (out)
   				       int32_t *error,      // (out)
				 const int32_t  vname_len,  // (in)
				 const int32_t  key_len,    // (in)
				 const int32_t  value_len ) // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  char _value[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_attribute_text_c( *fid, _vname, _key, _value, value_len );

  l = value_len > File_HLONG ? File_HLONG : value_len;
  cstr2fstr(value, _value, l);
}

void file_get_attribute_int_c_( const int32_t *fid,       // (in)
				const char    *vname,     // (in)
				const char    *key,       // (in)
				const int32_t *len,       // (in)
				      int32_t *value,     // (out)
				      int32_t *error,     // (out)
				const int32_t  vname_len, // (in)
				const int32_t  key_len )  // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_attribute_int_c( *fid, _vname, _key, value, (size_t)*len );
}

void file_get_attribute_float_c_( const int32_t *fid,       // (in)
				  const char    *vname,     // (in)
				  const char    *key,       // (in)
				  const int32_t *len,       // (in)
				        float   *value,     // (out)
				        int32_t *error,     // (out)
				  const int32_t  vname_len, // (in)
				  const int32_t  key_len )  // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_attribute_float_c( *fid, _vname, _key, value, (size_t)*len );
}

void file_get_attribute_double_c_( const int32_t *fid,       // (in)
				   const char    *vname,     // (in)
				   const char    *key,       // (in)
				   const int32_t *len,       // (in)
				         double  *value,     // (out)
				         int32_t *error,     // (out)
				   const int32_t  vname_len, // (in)
				   const int32_t  key_len )  // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_get_attribute_double_c( *fid, _vname, _key, value, (size_t)*len );
}

void file_set_attribute_text_c_( const int32_t *fid,        // (in)
				 const char    *vname,      // (in)
				 const char    *key,        // (in)
				 const char    *value,      // (in)
				       int32_t *error,      // (out)
				 const int32_t  vname_len,  // (in)
				 const int32_t  key_len,    // (in)
				 const int32_t  value_len ) // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  char _value[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  l = value_len > File_HLONG ? File_HLONG : value_len;
  fstr2cstr(_value, value, l);

  *error = file_set_attribute_text_c( *fid, _vname, _key, _value );
}

void file_set_attribute_int_c_( const int32_t *fid,       // (in)
				const char    *vname,     // (in)
				const char    *key,       // (in)
				const int32_t *value,     // (in)
				const int32_t *len,       // (in)
				      int32_t *error,     // (out)
				const int32_t  vname_len, // (in)
				const int32_t  key_len )  // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_attribute_int_c( *fid, _vname, _key, value, (size_t)*len );
}

void file_set_attribute_float_c_( const int32_t *fid,       // (in)
				  const char    *vname,     // (in)
				  const char    *key,       // (in)
				  const float   *value,     // (in)
				  const int32_t *len,       // (in)
				        int32_t *error,     // (out)
				  const int32_t  vname_len, // (in)
				  const int32_t  key_len )  // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_attribute_float_c( *fid, _vname, _key, value, (size_t)*len );
}

void file_set_attribute_double_c_( const int32_t *fid,       // (in)
				   const char    *vname,     // (in)
				   const char    *key,       // (in)
				   const double  *value,     // (in)
				   const int32_t *len,       // (in)
				         int32_t *error,     // (out)
				   const int32_t  vname_len, // (in)
				   const int32_t  key_len )  // (in)
{
  char _vname[File_HSHORT+1];
  char _key[File_HLONG+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  l = key_len > File_HLONG ? File_HLONG : key_len;
  fstr2cstr(_key, key, l);

  *error = file_set_attribute_double_c( *fid, _vname, _key, value, (size_t)*len );
}

void file_add_associatedvariable_c_( const int32_t *fid,        // (in)
				     const char    *vname,      // (in)
				           int32_t *error,      // (out)
				     const int32_t  vname_len ) // (in)
{
  char _vname[File_HSHORT+1];
  int32_t l;

  l = vname_len > File_HSHORT ? File_HSHORT : vname_len;
  fstr2cstr(_vname, vname, l);

  *error = file_add_associatedvariable_c( *fid, _vname );
}

void file_set_tunits_c_( const int32_t *fid,        // (in)
			 const char    *time_units, // (in)
			       int32_t *error,      // (out)
			 const int32_t  len)        // (in)
{
  char _time_units[File_HMID+1];
  int32_t l;

  l = len > File_HMID ? File_HMID : len;
  fstr2cstr(_time_units, time_units, l);

  *error = file_set_tunits_c( *fid, _time_units );
}

void file_put_axis_c_( const int32_t *fid,          // (in)
		       const char    *name,         // (in)
		       const char    *desc,         // (in)
		       const char    *units,        // (in)
		       const char    *dim_name,     // (in)
		       const int32_t *dtype,        // (in)
		       const void    *val,          // (in)
		       const int32_t *size,         // (in)
		       const int32_t *precision,    // (in)
		             int32_t *error,        // (out)
		       const int32_t  name_len,     // (in)
		       const int32_t  desc_len,     // (in)
		       const int32_t  units_len,    // (in)
		       const int32_t  dim_name_len) // (in)
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

  *error = file_put_axis_c( *fid, _name, _desc, _units, _dim_name, *dtype, val, *size, *precision );
}

void file_def_axis_c_( const int32_t *fid,          // (in)
		       const char    *name,         // (in)
		       const char    *desc,         // (in)
		       const char    *units,        // (in)
		       const char    *dim_name,     // (in)
		       const int32_t *dtype,        // (in)
		       const int32_t *dim_size,     // (in)
		             int32_t *error,        // (out)
		       const int32_t  name_len,     // (in)
		       const int32_t  desc_len,     // (in)
		       const int32_t  units_len,    // (in)
		       const int32_t  dim_name_len) // (in)
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

  *error = file_def_axis_c( *fid, _name, _desc, _units, _dim_name, *dtype, *dim_size );
}

void file_write_axis_c_( const int32_t *fid,          // (in)
			 const char    *name,         // (in)
			 const void    *val,          // (in)
			 const int32_t *precision,    // (in)
			 const int32_t *start,        // (in)
			 const int32_t *count,        // (in)
			       int32_t *error,        // (out)
			 const int32_t  name_len)     // (in)
{
  char _name[File_HSHORT+1];
  int len;
  MPI_Offset start_[1], count_[1];

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  /* all axes are 1D */
  start_[0] = *start - 1;  /* C index is 0-based */
  count_[0] = *count;

  *error = file_write_axis_c( *fid, _name, val, *precision, start_, count_ );
}

void file_put_associatedcoordinate_c_( const int32_t *fid,          // (in)
				       const char    *name,         // (in)
				       const char    *desc,         // (in)
				       const char    *units,        // (in)
				       const char    *dim_names,    // (in)
				       const int32_t *ndims,        // (in)
				       const int32_t *dtype,        // (in)
				       const void    *val,          // (in)
				       const int32_t *precision,    // (in)
					      int32_t *error,        // (out)
				       const int32_t  name_len,     // (in)
				       const int32_t  desc_len,     // (in)
				       const int32_t  units_len,    // (in)
				       const int32_t  dim_name_len) // (in)
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
  *error = file_put_associatedcoordinate_c( *fid, _name, _desc, _units, (const char**)_dim_names, *ndims, *dtype, val, *precision );
}

void file_def_associatedcoordinate_c_( const int32_t *fid,          // (in)
				       const char    *name,         // (in)
				       const char    *desc,         // (in)
				       const char    *units,        // (in)
				       const char    *dim_names,    // (in)
				       const int32_t *ndims,        // (in)
				       const int32_t *dtype,        // (in)
				             int32_t *error,        // (out)
				       const int32_t  name_len,     // (in)
				       const int32_t  desc_len,     // (in)
				       const int32_t  units_len,    // (in)
				       const int32_t  dim_name_len) // (in)
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
  *error = file_def_associatedcoordinate_c( *fid, _name, _desc, _units, (const char**)_dim_names, *ndims, *dtype );
}

void file_write_associatedcoordinate_c_( const int32_t *fid,          // (in)
					 const char    *name,         // (in)
					 const void    *val,          // (in)
					 const int32_t *precision,    // (in)
					 const int32_t *ndims,        // (in)
					 const int32_t *start,        // (in)
					 const int32_t *count,        // (in)
					       int32_t *error,        // (out)
					 const int32_t  name_len)     // (in)
{
  char _name[File_HSHORT+1];
  int i, len;
  MPI_Offset start_[4], count_[4];
  /* all associated coordinate are up to 4D */

  len = name_len > File_HSHORT ? File_HSHORT : name_len;
  fstr2cstr(_name, name, len);

  for (i=0; i<*ndims; i++) {
      start_[i] = start[*ndims - i - 1] -  1;
      count_[i] = count[*ndims - i - 1];
  }
  *error = file_write_associatedcoordinate_c( *fid, _name, val, *precision, start_, count_ );
}

void file_add_variable_c_(       int32_t  *vid,         // (out)
			   const int32_t  *fid,         // (in)
			   const char     *varname,     // (in)
			   const char     *desc,        // (in)
			   const char     *units,       // (in)
			   const char     *dims,        // (in)
			   const int32_t  *ndims,       // (in)
			   const int32_t  *dtype,       // (in)
			   const real64_t *tint,        // (in)
			   const int32_t  *tavg,        // (in)
			         int32_t  *error,       // (out)
			   const int32_t   varname_len, // (in)
			   const int32_t   desc_len,    // (in)
			   const int32_t   units_len,   // (in)
			   const int32_t   dims_len)    // (in)
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

  *error = file_add_variable_c( vid, *fid, _varname, _desc, _units, (const char**)_dims, *ndims, *dtype, *tint, *tavg );

  for ( i=0; i<*ndims; i++ )
    free( _dims[i] );
  free( _dims );
}

void file_write_data_c_( const int32_t  *fid,       // (in)
			 const int32_t  *vid,       // (in)
			 const void     *var,       // (in)
			 const real64_t *t_start,   // (in)
			 const real64_t *t_end,     // (in)
			 const int32_t  *precision, // (in)
			 const int32_t  *ndims,     // (in)
			 const int32_t  *start,     // (in)
			 const int32_t  *count,     // (in)
			       int32_t  *error)     // (out)
{
  int i;
  MPI_Offset start_[4], count_[4]; /* assume max ndims is 4 */

  for (i=0; i<*ndims; i++) {
      start_[i] = start[*ndims - i - 1] -  1;
      count_[i] = count[*ndims - i - 1];
  }
  *error = file_write_data_c( *fid, *vid, var, *t_start, *t_end, *precision, start_, count_ );
}

void file_close_c_( const int32_t *fid ,   // (in)
		          int32_t *error ) // (out)
{
  *error = file_close_c( *fid );
}

void file_enddef_c_( const int32_t *fid ,   // (in)
		            int32_t *error ) // (out)
{
  *error = file_enddef_c( *fid );
}

void file_attach_buffer_c_( const int32_t *fid ,       // (in)
			          int32_t *buf_amount, // (out)
			          int32_t *error )     // (out)
{
  *error = file_attach_buffer_c( *fid, *buf_amount );
}

void file_detach_buffer_c_( const int32_t *fid ,       // (in)
			          int32_t *error )     // (out)
{
  *error = file_detach_buffer_c( *fid );
}

void file_flush_c_( const int32_t *fid ,   // (in)
		          int32_t *error ) // (out)
{
  *error = file_flush_c( *fid );
}
