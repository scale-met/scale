#ifndef __FILE_SCALE_H__
#define __FILE_SCALE_H__

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/* character length */
#define File_HSHORT   32
#define File_HMID    128
#define File_HLONG  1024

/* data type */
#define File_REAL4     0
#define File_REAL8     1
#define File_INTEGER2  2
#define File_INTEGER4  3
#define File_INTEGER8  4

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
  char     time_units[File_HMID];
  int32_t  fid;
} datainfo_t;


extern int32_t file_open( int32_t  *fid,   // (out)
			  char     *fname, // (in)
			  int32_t   mode,  // (in)
			  MPI_Comm  comm); // (in)

extern int32_t file_set_option( int32_t  fid,      // (in)
				char    *filetype, // (in)
				char    *key,      // (in)
				char    *val);     // (in)

extern int32_t file_get_datainfo( datainfo_t *dinfo,    // (out)
				  int32_t     fid,      // (in)
				  char       *varname,  // (in)
				  int32_t     step,     // (in)
				  int32_t     suppress);// (in)

extern int32_t file_read_data( void       *var,        // (out)
			       datainfo_t *dinfo,      // (in)
			       int32_t     precision); // (in)

extern int32_t file_read_data_par( void         *var,       // (out)
			           datainfo_t   *dinfo,     // (in)
			           MPI_Offset    ntypes,    // (in)
			           MPI_Datatype  dtype,     // (in)
			           MPI_Offset   *start,     // (in)
			           MPI_Offset   *count);    // (in)


extern int32_t file_get_global_attribute_text( int32_t  fid,   // (in)
					       char    *key,   // (in)
					       char    *value, // (out)
					       int32_t  len);  // (in)

extern int32_t file_get_global_attribute_int( int32_t  fid,   // (in)
					      char    *key,   // (in)
					      int32_t *value, // (out)
					      size_t   len);  // (in)

extern int32_t file_get_global_attribute_float( int32_t  fid,   // (in)
						char    *key,   // (in)
						float   *value, // (out)
						size_t   len);  // (in)

extern int32_t file_get_global_attribute_double( int32_t  fid,   // (in)
						 char    *key,   // (in)
						 double  *value, // (out)
						 size_t   len);  // (in)

extern int32_t file_set_global_attribute_text( int32_t  fid,    // (in)
					       char    *key,    // (in)
					       char    *value); // (in)

extern int32_t file_set_global_attribute_int( int32_t  fid,   // (in)
					      char    *key,   // (in)
					      int32_t *value, // (in)
					      size_t   len);  // (in)

extern int32_t file_set_global_attribute_float( int32_t  fid,   // (in)
						char    *key,   // (in)
						float   *value, // (in)
						size_t   len);  // (in)

extern int32_t file_set_global_attribute_double( int32_t  fid,   // (in)
						 char    *key,   // (in)
						 double  *value, // (in)
						 size_t   len);  // (in)

extern int32_t file_set_tunits( int32_t fid,          // (in)
				char    *time_units); // (in)

extern int32_t file_set_tattr( int32_t  fid,   // (in)
			       char    *vname, // (in)
			       char    *key,   // (in)
			       char    *val);  // (in)

extern int32_t file_put_axis( int32_t fid,        // (in)
			      char   *name,       // (in)
			      char   *desc,       // (in)
			      char   *units,      // (in)
			      char   *dim_name,   // (in)
			      int32_t dtype,      // (in)
			      void   *val,        // (in)
			      int32_t size,       // (in)
			      int32_t precision); // (in)

extern int32_t file_def_axis( int32_t fid,        // (in)
			      char   *name,       // (in)
			      char   *desc,       // (in)
			      char   *units,      // (in)
			      char   *dim_name,   // (in)
			      int32_t dtype,      // (in)
			      int32_t dim_size);  // (in)

extern int32_t file_write_axis( int32_t     fid,        // (in)
			        char       *name,       // (in)
			        void       *val,        // (in)
			        int32_t     precision,  // (in)
			        MPI_Offset *start,      // (in)
			        MPI_Offset *count);     // (in)

extern int32_t file_put_associated_coordinates( int32_t fid,        // (in)
						char   *name,       // (in)
						char   *desc,       // (in)
						char   *units,      // (in)
						char   **dim_names, // (in)
						int32_t ndims,      // (in)
						int32_t dtype,      // (in)
						void   *val,        // (in)
						int32_t precision); // (in)

extern int32_t file_def_associated_coordinates( int32_t fid,        // (in)
						char   *name,       // (in)
						char   *desc,       // (in)
						char   *units,      // (in)
						char   **dim_names, // (in)
						int32_t ndims,      // (in)
						int32_t dtype);     // (in)

extern int32_t file_write_associated_coordinates( int32_t     fid,        // (in)
						  char       *name,       // (in)
						  void       *val,        // (in)
						  int32_t     precision,  // (in)
						  MPI_Offset *start,      // (in)
						  MPI_Offset *count);     // (in)

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

extern int32_t file_write_data( int32_t     fid,       // (in)
                                int32_t     vid,       // (in)
			        void       *var,       // (in)
			        real64_t    t_start,   // (in)
			        real64_t    t_end,     // (in)
			        int32_t     precision, // (in)
			        MPI_Offset *start,     // (in)
			        MPI_Offset *count);    // (in)

extern int32_t file_enddef( int32_t fid ); // (in)

extern int32_t file_attach_buffer( int32_t fid, int32_t buf_amount ); // (in)

extern int32_t file_detach_buffer( int32_t fid ); // (in)

extern int32_t file_flush( int32_t fid ); // (in)

extern int32_t file_close( int32_t fid ); // (in)

#endif
