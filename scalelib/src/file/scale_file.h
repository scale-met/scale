#ifndef __SCALE_FILE_H__
#define __SCALE_FILE_H__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


#include "scale_file_const.h"


typedef float (real32_t);
typedef double(real64_t);

/* data item information */
typedef struct{
  char     varname[File_HSHORT];
  char     description[File_HMID];
  char     units[File_HSHORT];
  char     standard_name[File_HMID];
  int32_t  datatype;
  int32_t  rank;
  char     dim_name[File_HSHORT*RANK_MAX];
  int32_t  dim_size[RANK_MAX];
  int32_t  step;
  real64_t time_start;
  real64_t time_end;
  char     time_units[File_HMID];
  char     calendar[File_HSHORT];
  int32_t  natts;
  char     att_name[File_HSHORT*ATT_MAX];
  int32_t  att_type[ATT_MAX];
  int32_t  att_len[ATT_MAX];
  int32_t  fid;
} datainfo_t;


extern int file_open_c(       int  *fid,   // (out)
			const char *fname, // (in)
			const int   mode,  // (in)
			const int   comm); // (in)

extern int file_get_dim_length_c(       int  *len,      // (out)
				  const int   fid,      // (in)
				  const char *dimname); // (in)

extern int file_set_option_c( const int   fid,      // (in)
			      const char *filetype, // (in)
			      const char *key,      // (in)
			      const char *val);     // (in)

extern int file_get_nvars_c(       int *nvars, // (out)
			     const int  fid);  // (in)
			     

extern int file_get_varname_c(       char *name, // (out)
			       const int   fid,  // (in)
			       const int   vid,  // (in)
			       const int   len); // (in)

extern int file_get_datainfo_c(       datainfo_t *dinfo,      // (out)
			        const int         fid,        // (in)
				const char       *varname,    // (in)
				const int         step,       // (in)
				const bool        suppress ); // (in)

extern int file_get_step_size_c(       int  *len,      // (out)
				 const int   fid,      // (in)
				 const char *varname); // (in)
				 

extern int file_read_data_c(       void       *var,       // (out)
			     const datainfo_t *dinfo,     // (in)
			     const int         precision, // (in)
			     const int         ntypes,    // (in)
			     const int         dtype,     // (in)
			     const int        *start,     // (in)
			     const int        *count );   // (in)

extern int file_get_attribute_text_c(       char *value,    // (out)
				      const int   fid,      // (in)
				      const char *vname,    // (in)
				      const char *key,      // (in)
				      const bool  suppress, // (in)
				      const int   len);     // (in)

extern int file_get_attribute_int_c(       int *value,    // (out)
				     const int   fid,      // (in)
				     const char *vname,    // (in)
				     const char *key,      // (in)
				     const bool  suppress, // (in)
				     const int   len);     // (in)

extern int file_get_attribute_float_c(       float *value,    // (out)
				       const int    fid,      // (in)
				       const char  *vname,    // (in)
				       const char  *key,      // (in)
				       const bool   suppress, // (in)
				       const int    len);     // (in)

extern int file_get_attribute_double_c(       double *value,    // (out)
					const int     fid,      // (in)
					const char   *vname,    // (in)
					const char   *key,      // (in)
					const bool    suppress, // (in)
					const int     len);     // (in)

extern int file_set_attribute_text_c( const int   fid,    // (in)
				      const char *vname,  // (in)
				      const char *key,    // (in)
				      const char *value); // (in)

extern int file_set_attribute_int_c( const int   fid,   // (in)
				     const char *vname, // (in)
				     const char *key,   // (in)
				     const int  *value, // (in)
				     const int   len);  // (in)

extern int file_set_attribute_float_c( const int    fid,   // (in)
				       const char  *vname, // (in)
				       const char  *key,   // (in)
				       const float *value, // (in)
				       const int    len);  // (in)

extern int file_set_attribute_double_c( const int     fid,   // (in)
					const char   *vname, // (in)
					const char   *key,   // (in)
					const double *value, // (in)
					const int     len);  // (in)

extern int file_add_associatedvariable_c( const int   fid,    // (in)
					  const char *vname); // (in)

extern int file_set_tunits_c( const int   fid,        // (in)
			      const char *time_units, // (in)
			      const char *calendar);  // (in)

extern int file_put_axis_c( const int   fid,        // (in)
			    const char *name,       // (in)
			    const char *desc,       // (in)
			    const char *units,      // (in)
			    const char *dim_name,   // (in)
			    const int   dtype,      // (in)
			    const void *val,        // (in)
			    const int   size,       // (in)
			    const int   precision); // (in)

extern int file_def_axis_c( const int   fid,      // (in)
			    const char *name,     // (in)
			    const char *desc,     // (in)
			    const char *units,    // (in)
			    const char *dim_name, // (in)
			    const int   dtype,    // (in)
			    const int   dim_size, // (in)
			    const int   bounds);  // (in)

extern int file_write_axis_c( const int   fid,        // (in)
			      const char *name,       // (in)
			      const void *val,        // (in)
			      const int   precision,  // (in)
			      const int  *start,      // (in)
			      const int  *count);     // (in)

extern int file_put_associatedcoordinate_c( const int    fid,        // (in)
					    const char  *name,       // (in)
					    const char  *desc,       // (in)
					    const char  *units,      // (in)
					    const char **dim_names,  // (in)
					    const int    ndims,      // (in)
					    const int    dtype,      // (in)
					    const void  *val,        // (in)
					    const int    precision); // (in)

extern int file_def_associatedcoordinate_c( const int    fid,       // (in)
					    const char  *name,      // (in)
					    const char  *desc,      // (in)
					    const char  *units,     // (in)
					    const char **dim_names, // (in)
					    const int    ndims,     // (in)
					    const int    dtype);    // (in)

extern int file_write_associatedcoordinate_c( const int   fid,       // (in)
					      const char *name,      // (in)
					      const void *val,       // (in)
					      const int   ndims,     // (in)
					      const int   precision, // (in)
					      const int  *start,     // (in)
					      const int  *count);    // (in)

extern int file_add_variable_c(       int    *vid,     // (out)
				const int     fid,     // (in)
				const char   *varname, // (in)
				const char   *desc,    // (in)
				const char   *units,   // (in)
				const char   *stdname, // (in)
				const char  **dims,    // (in)
				const int     ndims,   // (in)
				const int     dtype,   // (in)
				const double  tint,    // (in)
				const bool    tavg);   // (in)

extern int file_write_data_c( const int     fid,       // (in)
			      const int     vid,       // (in)
			      const void   *var,       // (in)
			      const double  t_start,   // (in)
			      const double  t_end,     // (in)
			      const int     ndims,     // (in)
			      const int     precision, // (in)
			      const int    *start,     // (in)
			      const int    *count);    // (in)

extern int file_enddef_c( const int fid ); // (in)
extern int file_redef_c( const int fid ); // (in)

extern int file_attach_buffer_c( const int     fid,          // (in)
				 const int64_t buf_amount ); // (in)

extern int file_detach_buffer_c( const int fid ); // (in)

extern int file_flush_c( const int fid ); // (in)

extern int file_close_c( const int  fid,     // (in)
			 const bool abort ); // (in)

#endif
