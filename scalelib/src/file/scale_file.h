#ifndef __SCALE_FILE_H__
#define __SCALE_FILE_H__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
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


extern int32_t file_open_c(       int32_t  *fid,   // (out)
			    const char     *fname, // (in)
			    const int32_t   mode,  // (in)
			    const MPI_Comm  comm); // (in)

extern int32_t file_get_dim_length_c( const int32_t  fid,       // (in)
				      const char*    dimname,   // (in)
				            int32_t *len     ); // (out)

extern int32_t file_set_option_c( const int32_t  fid,      // (in)
				  const char    *filetype, // (in)
				  const char    *key,      // (in)
				  const char    *val);     // (in)

extern int32_t file_get_nvars_c( const int32_t  fid,     // (in)
				       int32_t *nvars ); // (out)

extern int32_t file_get_varname_c( const int32_t  fid,      // (in)
				   const int32_t  vid,      // (in)
				         char    *name,     // (out)
				   const int32_t  len    ); // (in)

extern int32_t file_get_datainfo_c(       datainfo_t *dinfo,      // (out)
				    const int32_t     fid,        // (in)
				    const char       *varname,    // (in)
				    const int32_t     step,       // (in)
				    const int32_t     suppress ); // (in)

extern int32_t file_read_data_c(       void       *var,       // (out)
				 const datainfo_t *dinfo,     // (in)
				 const int32_t     precision, // (in)
				 const MPI_Offset  ntypes,    // (in)
				 const MPI_Datatype dtype,    // (in)
				 const int32_t     *start,    // (in)
				 const int32_t     *count );   // (in)

extern int32_t file_get_attribute_text_c( const int32_t  fid,      // (in)
					  const char    *vname,    // (in)
					  const char    *key,      // (in)
					  const int32_t  suppress, // (in)
					        char    *value,    // (out)
					  const int32_t  len);     // (in)

extern int32_t file_get_attribute_int_c( const int32_t  fid,      // (in)
					 const char    *vname,    // (in)
					 const char    *key,      // (in)
					 const int32_t  suppress, // (in)
					       int32_t *value,    // (out)
					 const size_t   len);     // (in)

extern int32_t file_get_attribute_float_c( const int32_t  fid,      // (in)
					   const char    *vname,    // (in)
					   const char    *key,      // (in)
					   const int32_t  suppress, // (in)
					         float   *value,    // (out)
					   const size_t   len);     // (in)

extern int32_t file_get_attribute_double_c( const int32_t  fid,      // (in)
					    const char    *vname,    // (in)
					    const char    *key,      // (in)
					    const int32_t  suppress, // (in)
					          double  *value,    // (out)
					    const size_t   len);     // (in)

extern int32_t file_set_attribute_text_c( const int32_t  fid,    // (in)
					  const char    *vname,  // (in)
					  const char    *key,    // (in)
					  const char    *value); // (in)

extern int32_t file_set_attribute_int_c( const int32_t  fid,   // (in)
					 const char    *vname, // (in)
					 const char    *key,   // (in)
					 const int32_t *value, // (in)
					 const size_t   len);  // (in)

extern int32_t file_set_attribute_float_c( const int32_t  fid,   // (in)
					   const char    *vname, // (in)
					   const char    *key,   // (in)
					   const float   *value, // (in)
					   const size_t   len);  // (in)

extern int32_t file_set_attribute_double_c( const int32_t  fid,   // (in)
					    const char    *vname, // (in)
					    const char    *key,   // (in)
					    const double  *value, // (in)
					    const size_t   len);  // (in)

extern int32_t file_add_associatedvariable_c( const int32_t  fid,    // (in)
					      const char    *vname); // (in)

extern int32_t file_set_tunits_c( const int32_t fid,         // (in)
				  const char    *time_units, // (in)
				  const char    *calendar);  // (in)

extern int32_t file_put_axis_c( const int32_t fid,        // (in)
				const char   *name,       // (in)
				const char   *desc,       // (in)
				const char   *units,      // (in)
				const char   *dim_name,   // (in)
				const int32_t dtype,      // (in)
				const void   *val,        // (in)
				const int32_t size,       // (in)
				const int32_t precision); // (in)

extern int32_t file_def_axis_c( const int32_t fid,      // (in)
				const char   *name,     // (in)
				const char   *desc,     // (in)
				const char   *units,    // (in)
				const char   *dim_name, // (in)
				const int32_t dtype,    // (in)
				const int32_t dim_size, // (in)
				const int32_t bounds);  // (in)

extern int32_t file_write_axis_c( const int32_t     fid,        // (in)
				  const char       *name,       // (in)
				  const void       *val,        // (in)
				  const int32_t     precision,  // (in)
				  const MPI_Offset *start,      // (in)
				  const MPI_Offset *count);     // (in)

extern int32_t file_put_associatedcoordinate_c( const int32_t fid,        // (in)
						const char   *name,       // (in)
						const char   *desc,       // (in)
						const char   *units,      // (in)
						const char   **dim_names, // (in)
						const int32_t ndims,      // (in)
						const int32_t dtype,      // (in)
						const void   *val,        // (in)
						const int32_t precision); // (in)

extern int32_t file_def_associatedcoordinate_c( const int32_t fid,        // (in)
						const char   *name,       // (in)
						const char   *desc,       // (in)
						const char   *units,      // (in)
						const char   **dim_names, // (in)
						const int32_t ndims,      // (in)
						const int32_t dtype);     // (in)

extern int32_t file_write_associatedcoordinate_c( const int32_t     fid,        // (in)
						  const char       *name,       // (in)
						  const void       *val,        // (in)
						  const int32_t     precision,  // (in)
						  const MPI_Offset *start,      // (in)
						  const MPI_Offset *count);     // (in)

extern int32_t file_add_variable_c( const int32_t  fid,     // (in)
				    const char    *varname, // (in)
				    const char    *desc,    // (in)
				    const char    *units,   // (in)
				    const char    *stdname, // (in)
				    const char   **dims,    // (in)
				    const int32_t  ndims,   // (in)
				    const int32_t  dtype,   // (in)
				    const real64_t tint,    // (in)
				    const int32_t  tavg,    // (in)
				          int32_t *vid);    // (out)

extern int32_t file_write_data_c( const int32_t   fid,       // (in)
				  const int32_t   vid,       // (in)
				  const void     *var,       // (in)
				  const real64_t  t_start,   // (in)
				  const real64_t  t_end,     // (in)
				  const int32_t   precision, // (in)
				  const int32_t   ndims,     // (in)
				  const int32_t  *start,     // (in)
				  const int32_t  *count);    // (in)

extern int32_t file_enddef_c( const int32_t fid ); // (in)

extern int32_t file_attach_buffer_c( const int32_t fid,          // (in)
				     const int64_t buf_amount ); // (in)

extern int32_t file_detach_buffer_c( const int32_t fid ); // (in)

extern int32_t file_flush_c( const int32_t fid ); // (in)

extern int32_t file_close_c( const int32_t fid ); // (in)

#endif
