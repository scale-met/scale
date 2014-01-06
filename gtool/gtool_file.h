#ifndef __FILE_SCALE_H__
#define __FILE_SCALE_H__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* character length */
#define File_HSHORT  16
#define File_HMID    64
#define File_HLONG  256

/* data type */
#define File_REAL4     0
#define File_REAL8     1
#define File_INTEGER4  2
#define File_INTEGER8  3

/* action type */
#define File_FREAD   0
#define File_FWRITE  1
#define File_FAPPEND 2

/* return type */
#define ERROR_CODE          -1
#define SUCCESS_CODE         0
#define ALREADY_CLOSED_CODE  1
#define ALREADY_EXISTED_CODE 2

/* limit of dimension */
#define MAX_RANK 10


typedef float (real32_t);
typedef double(real64_t);

/* data item information */
typedef struct{
  char     varname[File_HSHORT];
  char     description[File_HMID];
  char     units[File_HSHORT];
  int32_t  datatype;
  int32_t  rank;
  char     dim_name[File_HSHORT*MAX_RANK];
  int32_t  dim_size[MAX_RANK];
  int32_t  step;
  real64_t time_start;
  real64_t time_end;
  int32_t  fid;
} datainfo_t; 


extern int32_t file_open( int32_t *fid,   // (out)
			  char    *fname, // (in)
			  int32_t  mode); // (in)

extern int32_t file_set_option( int32_t  fid,      // (in)
				char    *filetype, // (in)
				char    *key,      // (in)
				char    *val);     // (in)

extern int32_t file_get_datainfo( datainfo_t *dinfo,   // (out)
				  int32_t     fid,     // (in)
				  char       *varname, // (in)
				  int32_t     step);   // (in)

extern int32_t file_read_data( void       *var,        // (out)
			       datainfo_t *dinfo,      // (in)
			       int32_t     precision); // (in)

extern int32_t file_set_global_attributes( int32_t  fid,          // (in)
					   char    *title,        // (in)
					   char    *source,       // (in)
					   char    *institution,  // (in)
					   char    *time_units,   // (in)
					   int32_t  nodeid,       // (in)
					   int32_t *nodeidx,      // (in)
					   int32_t  nodeidx_dim); // (in)

extern int32_t file_put_axis( int32_t fid,        // (in)
			      char   *name,       // (in)
			      char   *desc,       // (in)
			      char   *units,      // (in)
			      char   *dim_name,   // (in)
			      int32_t dtype,      // (in)
			      void   *val,        // (in)
			      int32_t size,       // (in)
			      int32_t precision); // (in)

extern int32_t file_put_associated_coordinates( int32_t fid,        // (in)
						char   *name,       // (in)
						char   *desc,       // (in)
						char   *units,      // (in)
						char   **dim_names, // (in)
						int32_t ndims,      // (in)
						int32_t dtype,      // (in)
						void   *val,        // (in)
						int32_t precision); // (in)

extern int32_t file_add_variable( int32_t *vid,     // (out)
				  int32_t  fid,     // (in)
				  char    *varname, // (in)
				  char    *desc,    // (in)
				  char    *units,   // (in)
				  char   **dims,    // (in)
				  int32_t  ndims,   // (in)
				  int32_t  dtype,   // (in)
				  real64_t tint,    // (in)
				  int32_t  tavg);   // (in)

extern int32_t file_write_data( int32_t  vid,        // (in)
				void    *var,        // (in)
				real64_t t_start,    // (in)
				real64_t t_end,      // (in)
				int32_t  precision); // (in)

extern int32_t file_close( int32_t fid ); // (in)

#endif
