#include "gtool_file.h"
#ifndef MPI_INCLUDED
#define MPI_INCLUDED
#endif
#include "netcdf.h"

#define RMISS -9.9999e+30
#define TEPS 1e-6
#define NTMAX 102400

#define MIN(a,b) ((a)<(b) ? (a) : (b))

static int32_t ERROR_SUPPRESS = 0;

#define CHECK_ERROR(func)                                       \
  {                                                             \
    int status_ = (func);                                       \
    if (status_ != NC_NOERR) {                                  \
      if ( ! ERROR_SUPPRESS ) {                                 \
        fprintf(stderr, "Error: at l%d in %s\n", __LINE__, __FILE__);   \
        fprintf(stderr, "       %s\n", nc_strerror(status_));   \
      }                                                         \
      return ERROR_CODE;                                        \
    }                                                           \
  }

#ifdef PNETCDF
#include "pnetcdf.h"
#define CHECK_PNC_ERROR(func)                                   \
  {                                                             \
    int status_ = (func);                                       \
    if (status_ != NC_NOERR) {                                  \
      if ( ! ERROR_SUPPRESS ) {                                 \
        fprintf(stderr, "Error: at l%d in %s\n", __LINE__, __FILE__);   \
        fprintf(stderr, "       %s\n", ncmpi_strerror(status_));        \
      }                                                         \
      return ERROR_CODE;                                        \
    }                                                           \
  }
#else
#define CHECK_PNC_ERROR(func)                                   \
  {                                                             \
    fprintf(stderr, "pnetCDF is necessary for shared_mode\n");  \
    fprintf(stderr, "Please re-compile with pnetCDF\n");        \
    return ERROR_CODE;                                          \
  }
#define ncmpi_inq_attid(a,b,c,d) NC2_ERR
#define ncmpi_inq_varid(a,b,c)   NC2_ERR
#define ncmpi_inq_dimid(a,b,c)   NC2_ERR
#endif

#define NCTYPE2TYPE(nctype, type)                               \
  {                                                             \
  switch ( nctype ) {                                           \
  case NC_FLOAT:                                                \
    type = File_REAL4;                                          \
    break;                                                      \
  case NC_DOUBLE:                                               \
    type = File_REAL8;                                          \
    break;                                                      \
  case NC_SHORT:                                                \
    type = File_INTEGER2;                                       \
    break;                                                      \
  default:                                                      \
    fprintf(stderr, "unsupported data type: %d\n", xtype);      \
    return ERROR_CODE;                                          \
  }                                                             \
  }

#define TYPE2NCTYPE(type, nctype)                               \
  {                                                             \
  switch ( type ) {                                             \
  case File_REAL4:                                              \
    nctype = NC_FLOAT;                                          \
    break;                                                      \
  case File_REAL8:                                              \
    nctype = NC_DOUBLE;                                         \
    break;                                                      \
  default:                                                      \
    fprintf(stderr, "unsupported data type: %d\n", xtype);      \
    return ERROR_CODE;                                          \
  }                                                             \
  }


#define DEFAULT_DEFLATE_LEVEL 2

typedef struct {
  int ncid;
  char time_units[File_HMID+1];
  int deflate_level;
#if defined(NETCDF3) || defined(PNETCDF)
  int defmode;
#endif
  int shared_mode;
  char fname[256];  // used for debugging, to be deleted
} fileinfo_t;

typedef struct {
  int ncid;
  int dimid;
  int varid;
  int bndsid;
  int count;
  real64_t t;
  real64_t tint;
  real64_t *tval;
  char name[File_HSHORT+1];
} tdim_t;

typedef struct {
  int ncid;
  int varid;
  tdim_t *t;
  size_t *start;
  size_t *count;
} varinfo_t;

// Keep consistency with "File_nfile_max" in gtool_file.f90
#define FILE_MAX 512
// Keep consistency with "File_nvar_max" in gtool_file.f90
#define VAR_MAX 40960

static fileinfo_t *files[FILE_MAX];
static int nfile = 0;
static varinfo_t *vars[VAR_MAX];
static int nvar = 0;
static tdim_t *tdims[VAR_MAX];
static int nt = 0;


int32_t file_open( int32_t  *fid,     // (out)
                   char     *fname,   // (in)
                   int32_t   mode,    // (in)
                   MPI_Comm  comm )   // (in)
{
  int ncid;
  int len;
  int shared_mode;
  char _fname[File_HLONG+4];

  if ( nfile >= FILE_MAX ) {
    fprintf(stderr, "exceed max number of file limit\n");
    return ERROR_CODE;
  }

  len = strlen(fname);
  strcpy(_fname, fname);
  if (fname[len-3] != '.' || fname[len-2] != 'n' || fname[len-1] != 'c' )
    strcat(_fname, ".nc");

  if ( comm == MPI_COMM_NULL || comm == MPI_COMM_SELF )
    shared_mode = 0;
  else
    shared_mode = 1;

  switch ( mode ) {
  case File_FREAD:
    if ( shared_mode )
      CHECK_PNC_ERROR( ncmpi_open(comm, _fname, NC_NOWRITE, MPI_INFO_NULL, &ncid) )
    else
      CHECK_ERROR( nc_open(_fname, NC_NOWRITE, &ncid) )
    break;
  case File_FWRITE:
    if ( shared_mode )
      CHECK_PNC_ERROR( ncmpi_create(comm, _fname, NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncid) )
    else
#ifdef NETCDF3
      CHECK_ERROR( nc_create(_fname, NC_CLOBBER|NC_64BIT_OFFSET, &ncid) )
#else
      CHECK_ERROR( nc_create(_fname, NC_CLOBBER|NC_NETCDF4, &ncid) )
#endif
    break;
  case File_FAPPEND:
    if ( shared_mode )
      CHECK_PNC_ERROR( ncmpi_open(comm, _fname, NC_WRITE, MPI_INFO_NULL, &ncid) )
    else
      CHECK_ERROR( nc_open(_fname, NC_WRITE, &ncid) )
    break;
  default:
    fprintf(stderr, "invalid mode type\n");
    return ERROR_CODE;
  }

  files[nfile] = (fileinfo_t*) malloc(sizeof(fileinfo_t));
  files[nfile]->ncid = ncid;
  files[nfile]->deflate_level = DEFAULT_DEFLATE_LEVEL;
#if defined(NETCDF3) || defined(PNETCDF)
  if ( mode == File_FWRITE )
    files[nfile]->defmode = 1;
  else
    files[nfile]->defmode = 0;
#endif
  files[nfile]->shared_mode = shared_mode;  /* shared-file I/O mode */
  strcpy(files[nfile]->fname, fname);
  *fid = nfile;
  nfile++;

  return SUCCESS_CODE;
}

int32_t file_set_option( int32_t fid,    // (in)
                         char* filetype, // (in)
                         char* key,      // (in)
                         char* val)      // (in)
{
  if ( strcmp(filetype, "netcdf") != 0 ) return SUCCESS_CODE;

  if ( strcmp(key, "deflate_level") == 0 ) {
    if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
    files[fid]->deflate_level = atoi(val);
    return SUCCESS_CODE;
  } else {
    return ERROR_CODE;
  }
}

int32_t file_get_datainfo( datainfo_t *dinfo,   // (out)
                           int32_t     fid,     // (in)
                           char*       varname, // (in)
                           int32_t     step,    // (in)
                           int32_t     suppress)// (in)
{
  int ncid, varid;
  nc_type xtype;
  int rank;
  int dimids[MAX_RANK], tdim, uldims[NC_MAX_DIMS];
  char name[NC_MAX_NAME+1];
  char *buf;
  size_t size, len;
  int i, n;

  ERROR_SUPPRESS = suppress;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;
  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, varname, &varid) )
  else
    CHECK_ERROR( nc_inq_varid(ncid, varname, &varid) )

  // fid
  dinfo->fid = fid;
  // varname
  strcpy(dinfo->varname, varname);
  if ( files[fid]->shared_mode ) {
    MPI_Offset l;
    // description
    CHECK_PNC_ERROR( ncmpi_inq_attlen  (ncid, varid, "long_name", &l) )
    buf = (char*) malloc(l+1);
    CHECK_PNC_ERROR( ncmpi_get_att_text(ncid, varid, "long_name", buf) )
    for (i=0; i<MIN(File_HMID-1,l); i++)
      dinfo->description[i] = buf[i];
    dinfo->description[i] = '\0';
    free(buf);
    // units
    CHECK_PNC_ERROR( ncmpi_inq_attlen  (ncid, varid, "units", &l) )
    buf = (char*) malloc(l+1);
    CHECK_PNC_ERROR( ncmpi_get_att_text(ncid, varid, "units", buf) )
    for (i=0; i<MIN(File_HSHORT-1,l); i++)
      dinfo->units[i] = buf[i];
    dinfo->units[i] = '\0';
    free(buf);
    // datatype
    CHECK_PNC_ERROR( ncmpi_inq_vartype(ncid, varid, &xtype) )
    NCTYPE2TYPE(xtype, dinfo->datatype);
    // rank
    CHECK_PNC_ERROR( ncmpi_inq_varndims(ncid, varid, &rank) )
    CHECK_PNC_ERROR( ncmpi_inq_vardimid(ncid, varid, dimids) )
#if 1
    CHECK_PNC_ERROR( ncmpi_inq_unlimdim(ncid, uldims) )
    n = 1;
#else
    CHECK_PNC_ERROR( ncmpi_inq_unlimdims(ncid, &n, uldims) )
#endif
  }
  else {
    size_t l;
    // description
    CHECK_ERROR( nc_inq_attlen  (ncid, varid, "long_name", &l) )
    buf = (char*) malloc(l+1);
    CHECK_ERROR( nc_get_att_text(ncid, varid, "long_name", buf) )
    for (i=0; i<MIN(File_HMID-1,l); i++)
      dinfo->description[i] = buf[i];
    dinfo->description[i] = '\0';
    free(buf);
    // units
    CHECK_ERROR( nc_inq_attlen  (ncid, varid, "units", &l) )
    buf = (char*) malloc(l+1);
    CHECK_ERROR( nc_get_att_text(ncid, varid, "units", buf) )
    for (i=0; i<MIN(File_HSHORT-1,l); i++)
      dinfo->units[i] = buf[i];
    dinfo->units[i] = '\0';
    free(buf);
    // datatype
    CHECK_ERROR( nc_inq_vartype(ncid, varid, &xtype) )
    NCTYPE2TYPE(xtype, dinfo->datatype);
    // rank
    CHECK_ERROR( nc_inq_varndims(ncid, varid, &rank) )
    CHECK_ERROR( nc_inq_vardimid(ncid, varid, dimids) )
#ifdef NETCDF3
    CHECK_ERROR( nc_inq_unlimdim(ncid, uldims) )
    n = 1;
#else
    CHECK_ERROR( nc_inq_unlimdims(ncid, &n, uldims) )
#endif
  }

  tdim = -1;
  for ( i=0; i<n; i++ ) {
    if ( uldims[i] == dimids[0] ) {
      tdim = uldims[i];
      break;
    }
  }
  if (rank > MAX_RANK) {
    fprintf(stderr, "rank exceeds limit: %d\n", rank);
    return ERROR_CODE;
  }
  dinfo->rank = tdim >= 0 ? rank -1 : rank; // do not count time dimension
  // dim_name and dim_size
  for (i=0; i<dinfo->rank; i++) {
    // note: C and Fortran orders are opposite
    if ( files[fid]->shared_mode ) {
      MPI_Offset size_;
      CHECK_PNC_ERROR( ncmpi_inq_dim(ncid, dimids[rank-i-1], name, &size_) )
      size = (size_t)size_;
    }
    else
      CHECK_ERROR( nc_inq_dim(ncid, dimids[rank-i-1], name, &size) )
    if ( strlen(name) > File_HSHORT-1 ) {
      fprintf(stderr, "Length of the dimension name (%s) is too long (should be < %d).\n", name, File_HSHORT);
      return ERROR_CODE;
    }
    strncpy(dinfo->dim_name+i*File_HSHORT, name, File_HSHORT);
    dinfo->dim_size[i] = size;
  }

  dinfo->step = step;
  if ( tdim >= 0 ) {
    if ( files[fid]->shared_mode ) {
      MPI_Offset idx[2];
      MPI_Offset l;
      // time_end
      CHECK_PNC_ERROR( ncmpi_inq_dimname(ncid, tdim, name) )
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, name, &varid) )
      idx[0] = step - 1;
      CHECK_PNC_ERROR( ncmpi_get_var1_double_all(ncid, varid, idx, &(dinfo->time_end)) )
      // time_start
      strcat(name, "_bnds");
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, name, &varid) )
      idx[1] = 0;
      CHECK_PNC_ERROR( ncmpi_get_var1_double_all(ncid, varid, idx, &(dinfo->time_start)) )
      // units
      CHECK_PNC_ERROR( ncmpi_inq_attlen  (ncid, varid, "units", &l) )
      buf = (char*) malloc(l+1);
      CHECK_PNC_ERROR( ncmpi_get_att_text(ncid, varid, "units", buf) )
      for (i=0; i<MIN(File_HMID-1,l); i++)
        dinfo->time_units[i] = buf[i];
      dinfo->time_units[i] = '\0';
      free(buf);
    } else {
      size_t idx[2];
      size_t l;
      // time_end
      CHECK_ERROR( nc_inq_dimname(ncid, tdim, name) )
      CHECK_ERROR( nc_inq_varid(ncid, name, &varid) )
      idx[0] = step - 1;
      CHECK_ERROR( nc_get_var1_double(ncid, varid, idx, &(dinfo->time_end)) )
      // time_start
      strcat(name, "_bnds");
      CHECK_ERROR( nc_inq_varid(ncid, name, &varid) )
      idx[1] = 0;
      CHECK_ERROR( nc_get_var1_double(ncid, varid, idx, &(dinfo->time_start)) )
      // units
      CHECK_ERROR( nc_inq_attlen  (ncid, varid, "units", &l) )
      buf = (char*) malloc(l+1);
      CHECK_ERROR( nc_get_att_text(ncid, varid, "units", buf) )
      for (i=0; i<MIN(File_HMID-1,l); i++)
        dinfo->time_units[i] = buf[i];
      dinfo->time_units[i] = '\0';
      free(buf);
    }
  }
  ERROR_SUPPRESS = 0;

  return SUCCESS_CODE;
}

int32_t file_read_data( void       *var,        // (out)
                        datainfo_t *dinfo,      // (in)
                        int32_t     precision)  // (in)
{
  int ncid, varid;
  int rank;
  size_t *start, *count;
  int i;
  int status;

  if ( files[dinfo->fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[dinfo->fid]->ncid;
  CHECK_ERROR( nc_inq_varid(ncid, dinfo->varname, &varid) )

  CHECK_ERROR( nc_inq_varndims(ncid, varid, &rank) )
  start = (size_t*) malloc(sizeof(size_t)*rank);
  count = (size_t*) malloc(sizeof(size_t)*rank);
  for (i=0; i<dinfo->rank; i++) {
    // note: C and Fortran orders are opposite
    start[rank -i-1] = 0;
    count[rank -i-1] = dinfo->dim_size[i];
  }
  if (rank > dinfo->rank) { // have time dimension
    start[0] = dinfo->step - 1;
    count[0] = 1;
  }
  switch ( precision ) {
  case 8:
    status = nc_get_vara_double(ncid, varid, start, count, (double*)var);
    free(start);
    free(count);
    CHECK_ERROR(status);
    break;
  case 4:
    status = nc_get_vara_float(ncid, varid, start, count, (float*)var);
    free(start);
    free(count);
    CHECK_ERROR(status);
    break;
  default:
    free(start);
    free(count);
    fprintf(stderr, "unsupported data precision: %d\n", precision );
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_read_data_par( void         *var,        // (out)
                            datainfo_t   *dinfo,      // (in)
                            MPI_Offset    ntypes,     // (in)
                            MPI_Datatype  dtype,      // (in)
                            MPI_Offset   *start,      // (in)
                            MPI_Offset   *count)      // (in)
{
  int ncid, varid, rank;

  if ( files[dinfo->fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[dinfo->fid]->ncid;
  CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, dinfo->varname, &varid) )

  CHECK_PNC_ERROR( ncmpi_inq_varndims(ncid, varid, &rank) )
  if (rank > dinfo->rank) { // have time dimension
    start[0] = dinfo->step - 1;
    count[0] = 1;
  } else {
    start = start + 1;
    count = count + 1;
  }

  CHECK_PNC_ERROR( ncmpi_iget_vara(ncid, varid, start, count, var, ntypes, dtype, NULL) )

  return SUCCESS_CODE;
}

int32_t file_get_global_attribute_text( int32_t  fid,   // (in)
                                        char    *key,   // (in)
                                        char    *value, // (out)
                                        int32_t  len )  // (in)
{
  int ncid;
  size_t l;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  CHECK_ERROR( nc_inq_attlen(ncid, NC_GLOBAL, key, &l) )
  if ( len < l+1 ) return ERROR_CODE;

  CHECK_ERROR( nc_get_att_text(ncid, NC_GLOBAL, key, value) )
  value[l] = '\0';

  return SUCCESS_CODE;
}

int32_t file_get_global_attribute_int( int32_t  fid,   // (in)
                                       char    *key,   // (in)
                                       int     *value, // (out)
                                       size_t   len )  // (in)
{
  int ncid;
  size_t l;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  CHECK_ERROR( nc_inq_attlen(ncid, NC_GLOBAL, key, &l) )
  if ( len < l ) return ERROR_CODE;
  CHECK_ERROR( nc_get_att_int(ncid, NC_GLOBAL, key, value) )

  return SUCCESS_CODE;
}

int32_t file_get_global_attribute_float( int32_t  fid,   // (in)
                                         char    *key,   // (in)
                                         float   *value, // (out)
                                         size_t   len )  // (in)
{
  int ncid;
  size_t l;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  CHECK_ERROR( nc_inq_attlen(ncid, NC_GLOBAL, key, &l) )
  if ( len < l ) return ERROR_CODE;
  CHECK_ERROR( nc_get_att_float(ncid, NC_GLOBAL, key, value) )

  return SUCCESS_CODE;
}

int32_t file_get_global_attribute_double( int32_t  fid,   // (in)
                                          char    *key,   // (in)
                                          double  *value, // (out)
                                          size_t   len )  // (in)
{
  int ncid;
  size_t l;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  CHECK_ERROR( nc_inq_attlen(ncid, NC_GLOBAL, key, &l) )
  if ( len < l ) return ERROR_CODE;
  CHECK_ERROR( nc_get_att_double(ncid, NC_GLOBAL, key, value) )

  return SUCCESS_CODE;
}

int32_t file_set_global_attribute_text( int32_t  fid,    // (in)
                                        char    *key,    // (in)
                                        char    *value ) // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, NC_GLOBAL, key, strlen(value), value) )
  }
  else {
    CHECK_ERROR( nc_put_att_text(ncid, NC_GLOBAL, key, strlen(value), value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_global_attribute_int( int32_t  fid,   // (in)
                                       char    *key,   // (in)
                                       int     *value, // (in)
                                       size_t   len )  // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_put_att_int(ncid, NC_GLOBAL, key, NC_INT, len, value) )
  }
  else {
    CHECK_ERROR( nc_put_att_int(ncid, NC_GLOBAL, key, NC_INT, len, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_global_attribute_float( int32_t  fid,   // (in)
                                         char    *key,   // (in)
                                         float   *value, // (in)
                                         size_t   len )  // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_put_att_float(ncid, NC_GLOBAL, key, NC_FLOAT, len, value) )
  }
  else {
    CHECK_ERROR( nc_put_att_float(ncid, NC_GLOBAL, key, NC_FLOAT, len, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_global_attribute_double( int32_t  fid,   // (in)
                                          char    *key,   // (in)
                                          double  *value, // (in)
                                          size_t   len )  // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_put_att_double(ncid, NC_GLOBAL, key, NC_DOUBLE, len, value) )
  }
  else {
    CHECK_ERROR( nc_put_att_double(ncid, NC_GLOBAL, key, NC_DOUBLE, len, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_tunits( int32_t fid,         // (in)
                         char    *time_units) // (in)
{
  strcpy(files[fid]->time_units, time_units);

  return SUCCESS_CODE;
}

int32_t file_set_tattr( int32_t  fid,   // (in)
                        char    *vname, // (in)
                        char    *key,   // (in)
                        char    *val)   // (in)
{
  int ncid;
  int varid;
  int attid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    if ( ncmpi_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, key, strlen(val), val) )
  }
  else {
    CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    if ( nc_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
    if (files[fid]->defmode == 0) {
      CHECK_ERROR( nc_redef(ncid) )
      files[fid]->defmode = 1;
    }
#endif

    CHECK_ERROR( nc_put_att_text(ncid, varid, key, strlen(val), val) )
  }

  return SUCCESS_CODE;
}

int32_t file_put_axis( int32_t fid,        // (in)
                       char   *name,       // (in)
                       char   *desc,       // (in)
                       char   *units,      // (in)
                       char   *dim_name,   // (in)
                       int32_t dtype,      // (in)
                       void*   val,        // (in)
                       int32_t size,       // (in)
                       int32_t precision)  // (in)
{
  int ncid, dimid, varid;
  nc_type xtype = -1;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( nc_inq_varid(ncid, name, &varid) == NC_NOERR ) // check if existed
    return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  if ( nc_inq_dimid(ncid, dim_name, &dimid) != NC_NOERR ) // check if existed
    CHECK_ERROR( nc_def_dim(ncid, dim_name, size, &dimid) )

  TYPE2NCTYPE(dtype, xtype);
  CHECK_ERROR( nc_def_var(ncid, name, xtype, 1, &dimid, &varid) )
  CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
  CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) )

#ifdef NETCDF3
  CHECK_ERROR( nc_enddef(ncid) )
  files[fid]->defmode = 0;
#endif

  switch ( precision ) {
  case 8:
    CHECK_ERROR( nc_put_var_double(ncid, varid, (double*)val) )
    break;
  case 4:
    CHECK_ERROR( nc_put_var_float(ncid, varid, (float*)val) )
    break;
  default:
    fprintf(stderr, "unsupported data precision: %d\n", precision);
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_def_axis( int32_t fid,        // (in)
                       char   *name,       // (in)
                       char   *desc,       // (in)
                       char   *units,      // (in)
                       char   *dim_name,   // (in)
                       int32_t dtype,      // (in)
                       int32_t dim_size)   // (in)
{
  int ncid, dimid, varid;
  nc_type xtype = -1;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    if ( ncmpi_inq_varid(ncid, name, &varid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    if ( ncmpi_inq_dimid(ncid, dim_name, &dimid) != NC_NOERR ) // check if existed
      CHECK_PNC_ERROR( ncmpi_def_dim(ncid, dim_name, dim_size, &dimid) )

    TYPE2NCTYPE(dtype, xtype);
    CHECK_PNC_ERROR( ncmpi_def_var(ncid, name, xtype, 1, &dimid, &varid) )
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "units", strlen(units), units) )
  }
  else {
    if ( nc_inq_varid(ncid, name, &varid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
    if (files[fid]->defmode == 0) {
      CHECK_ERROR( nc_redef(ncid) )
      files[fid]->defmode = 1;
    }
#endif

    if ( nc_inq_dimid(ncid, dim_name, &dimid) != NC_NOERR ) // check if existed
      CHECK_ERROR( nc_def_dim(ncid, dim_name, dim_size, &dimid) )

    TYPE2NCTYPE(dtype, xtype);
    CHECK_ERROR( nc_def_var(ncid, name, xtype, 1, &dimid, &varid) )
    CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) )
  }

  return SUCCESS_CODE;
}

int32_t file_write_axis( int32_t     fid,       // (in)
                         char       *name,      // (in)
                         void       *val,       // (in)
                         int32_t     precision, // (in)
                         MPI_Offset *start,     // (in)
                         MPI_Offset *count)     // (in)
{
  int ncid, varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, name, &varid) )
  else
    CHECK_ERROR( nc_inq_varid(ncid, name, &varid) )

#ifdef NETCDF3
  if ( (!files[fid]->shared_mode) && files[fid]->defmode == 1) {
    CHECK_ERROR( nc_enddef(ncid) )
    files[fid]->defmode = 0;
  }
#endif

  switch ( precision ) {
  case 8:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_iput_vara_double(ncid, varid, start, count, val, NULL) )
    else
      CHECK_ERROR( nc_put_var_double(ncid, varid, (double*)val) )
    break;
  case 4:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_iput_vara_float(ncid, varid, start, count, val, NULL) )
    else
      CHECK_ERROR( nc_put_var_float(ncid, varid, (float*)val) )
    break;
  default:
    fprintf(stderr, "unsupported data precision: %d\n", precision);
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_put_associated_coordinates( int32_t fid,        // (in)
                                         char   *name,       // (in)
                                         char   *desc,       // (in)
                                         char   *units,      // (in)
                                         char   **dim_names, // (in)
                                         int32_t ndims,      // (in)
                                         int32_t dtype,      // (in)
                                         void*   val,        // (in)
                                         int32_t precision)  // (in)
{
  int ncid, *dimids, varid;
  nc_type xtype = -1;
  int i;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( nc_inq_varid(ncid, name, &varid) == NC_NOERR ) // check if existed
    return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  dimids = malloc(sizeof(int)*ndims);
  for (i=0; i<ndims; i++)
    CHECK_ERROR( nc_inq_dimid(ncid, dim_names[i], dimids+ndims-i-1) )

  TYPE2NCTYPE(dtype, xtype);

  CHECK_ERROR( nc_def_var(ncid, name, xtype, ndims, dimids, &varid) )
  CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
  CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) )
  free(dimids);

#ifdef NETCDF3
  CHECK_ERROR( nc_enddef(ncid) )
  files[fid]->defmode = 0;
#endif

  switch ( precision ) {
  case 8:
    CHECK_ERROR( nc_put_var_double(ncid, varid, (double*)val) )
    break;
  case 4:
    CHECK_ERROR( nc_put_var_float(ncid, varid, (float*)val) )
    break;
  default:
    fprintf(stderr, "unsupported data precision: %d\n", precision);
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_def_associated_coordinates( int32_t fid,        // (in)
                                         char   *name,       // (in)
                                         char   *desc,       // (in)
                                         char   *units,      // (in)
                                         char   **dim_names, // (in)
                                         int32_t ndims,      // (in)
                                         int32_t dtype)      // (in)
{
  int ncid, *dimids, varid;
  nc_type xtype = -1;
  int i;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    if ( ncmpi_inq_varid(ncid, name, &varid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    dimids = malloc(sizeof(int)*ndims);
    for (i=0; i<ndims; i++)
      CHECK_PNC_ERROR( ncmpi_inq_dimid(ncid, dim_names[i], dimids+ndims-i-1) )

    TYPE2NCTYPE(dtype, xtype);

    CHECK_PNC_ERROR( ncmpi_def_var(ncid, name, xtype, ndims, dimids, &varid) )
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "units", strlen(units), units) )
    free(dimids);
  }
  else {
    if ( nc_inq_varid(ncid, name, &varid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
    if (files[fid]->defmode == 0) {
      CHECK_ERROR( nc_redef(ncid) )
      files[fid]->defmode = 1;
    }
#endif

    dimids = malloc(sizeof(int)*ndims);
    for (i=0; i<ndims; i++)
      CHECK_ERROR( nc_inq_dimid(ncid, dim_names[i], dimids+ndims-i-1) )

    TYPE2NCTYPE(dtype, xtype);

    CHECK_ERROR( nc_def_var(ncid, name, xtype, ndims, dimids, &varid) )
    CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) )
    free(dimids);
  }

  return SUCCESS_CODE;
}

int32_t file_write_associated_coordinates( int32_t     fid,        // (in)
                                           char       *name,       // (in)
                                           void*       val,        // (in)
                                           int32_t     precision,  // (in)
                                           MPI_Offset *start,      // (in)
                                           MPI_Offset *count)      // (in)
{
  int ncid, varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, name, &varid) )
  else
    CHECK_ERROR( nc_inq_varid(ncid, name, &varid) )

#ifdef NETCDF3
  if ( (!files[fid]->shared_mode) && files[fid]->defmode == 1) {
    CHECK_ERROR( nc_enddef(ncid) )
    files[fid]->defmode = 0;
  }
#endif

  switch ( precision ) {
  case 8:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_iput_vara_double(ncid, varid, start, count, (double*)val, NULL) )
    else
      CHECK_ERROR( nc_put_var_double(ncid, varid, (double*)val) )
    break;
  case 4:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_iput_vara_float(ncid, varid, start, count, (float*)val, NULL) )
    else
      CHECK_ERROR( nc_put_var_float(ncid, varid, (float*)val) )
    break;
  default:
    fprintf(stderr, "unsupported data precision: %d\n", precision);
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_add_variable( int32_t *vid,     // (out)
                           int32_t  fid,     // (in)
                           char    *varname, // (in)
                           char    *desc,    // (in)
                           char    *units,   // (in)
                           char   **dims,    // (in)
                           int32_t  ndims,   // (in)
                           int32_t  dtype,   // (in)
                           real64_t tint,    // (in)
                           int32_t  tavg)    // (in)
{
  int ncid, varid, acid, *acdimids;
  int dimids[NC_MAX_DIMS], dimid;
  char tname[File_HSHORT+1];
  int tdimid, tvarid;
  nc_type xtype = -1;
  char buf[File_HMID+1];
  int i, j, k, n, m, err;
  int nndims;
  size_t size;
  double rmiss = RMISS;
  char coord[File_HMID+1];
  int has_assoc;
  int new;

  if ( nvar >= VAR_MAX ) {
    fprintf(stderr, "exceed max number of variable limit\n");
    return ERROR_CODE;
  }

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  vars[nvar] = (varinfo_t*) malloc(sizeof(varinfo_t));
  vars[nvar]->ncid  = ncid;
  vars[nvar]->t = NULL;
  vars[nvar]->start = NULL;
  vars[nvar]->count = NULL;

#if defined(NETCDF3) || defined(PNETCDF)
  if (files[fid]->defmode == 0) {
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_redef(ncid) )
    else
      CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  // get time variable
  if ( tint > 0.0 ) {
    for ( i=0; i<nt; i++ ) {
      if ( tdims[i] != NULL && // still opened
           tdims[i]->ncid == ncid && // same file
           tdims[i]->tint == tint ) { // same time interval
        vars[nvar]->t = tdims[i];
        break;
      }
    }
    if ( vars[nvar]->t == NULL ) {
      tdims[nt] = (tdim_t*) malloc(sizeof(tdim_t));
      tdims[nt]->ncid = ncid;
      tdims[nt]->count = -1;
      tdims[nt]->tint = tint;
      tdims[nt]->tval = (double*) malloc(sizeof(double)*NTMAX);
      // generate name
      if ( nt == 0 )
        strcpy(tname, "time");
      else
        sprintf(tname, "time%d", nt);
      strcpy(tdims[nt]->name, tname);
      // define time dimension and variable
      if ( files[fid]->shared_mode ) {
        CHECK_PNC_ERROR( ncmpi_def_dim(ncid, tname, 0, &tdimid) )
        tdims[nt]->dimid = tdimid;
        CHECK_PNC_ERROR( ncmpi_def_var(ncid, tname, NC_DOUBLE, 1, &tdimid, &tvarid) )
        tdims[nt]->varid = tvarid;
        strcpy(buf, "time");
        CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, tvarid, "long_name", strlen(buf), buf) )
        CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, tvarid, "units", strlen(files[fid]->time_units), files[fid]->time_units) )
        // define boundary variable
        if ( ncmpi_inq_dimid(ncid, "nv", &(dimids[1])) != NC_NOERR ) // first called
          CHECK_PNC_ERROR( ncmpi_def_dim(ncid, "nv", 2, &(dimids[1])) )
        sprintf(buf, "%s_bnds", tname);
        CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, tvarid, "bounds", strlen(buf), buf) )
        dimids[0] = tdimid;
        CHECK_PNC_ERROR( ncmpi_def_var(ncid, buf, NC_DOUBLE, 2, dimids, &tvarid) )
        tdims[nt]->bndsid = tvarid;
        CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, tvarid, "units", strlen(files[fid]->time_units), files[fid]->time_units) )
      }
      else {
        CHECK_ERROR( nc_def_dim(ncid, tname, 0, &tdimid) )
        tdims[nt]->dimid = tdimid;
        CHECK_ERROR( nc_def_var(ncid, tname, NC_DOUBLE, 1, &tdimid, &tvarid) )
        tdims[nt]->varid = tvarid;
        strcpy(buf, "time");
        CHECK_ERROR( nc_put_att_text(ncid, tvarid, "long_name", strlen(buf), buf) )
        CHECK_ERROR( nc_put_att_text(ncid, tvarid, "units", strlen(files[fid]->time_units), files[fid]->time_units) )
        // define boundary variable
        if ( nc_inq_dimid(ncid, "nv", &(dimids[1])) != NC_NOERR ) // first called
          CHECK_ERROR( nc_def_dim(ncid, "nv", 2, &(dimids[1])) )
        sprintf(buf, "%s_bnds", tname);
        CHECK_ERROR( nc_put_att_text(ncid, tvarid, "bounds", strlen(buf), buf) )
        dimids[0] = tdimid;
        CHECK_ERROR( nc_def_var(ncid, buf, NC_DOUBLE, 2, dimids, &tvarid) )
        tdims[nt]->bndsid = tvarid;
        CHECK_ERROR( nc_put_att_text(ncid, tvarid, "units", strlen(files[fid]->time_units), files[fid]->time_units) )
      }

      vars[nvar]->t = tdims[nt];
      nt++;
    }
  }

  // get dimension IDs
  // note: C and Fortran order are opposite
  n = ndims;
  if ( tint > 0.0 ) { // add time dimension
    dimids[0] = vars[nvar]->t->dimid;
    ndims++;
  }
  for (i=ndims-n; i<ndims; i++) dimids[i] = -1;

  has_assoc = 0;
  nndims = 0;
  for (i=0; i<n; i++) {
    //printf("%d %s\n", i, dims[i]);
    if ( files[fid]->shared_mode )
       err = ncmpi_inq_dimid(ncid, dims[i], &dimid);
    else
       err = nc_inq_dimid(ncid, dims[i], &dimid);
    if ( err == NC_NOERR ) {
      //printf("not assoc\n");
      new = 1;
      for (k=0; k<nndims; k++) {
        if (dimid == dimids[k]) {
          new = 0;
          break;
        }
      }
      if (new) {
        dimids[ndims-(++nndims)] = dimid;
      }
    } else {
      //printf("assoc\n");
      if ( files[fid]->shared_mode ) {
        CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, dims[i], &acid) )
        CHECK_PNC_ERROR( ncmpi_inq_varndims(ncid, acid, &m) )
        acdimids = (int*) malloc((sizeof(int)*m));
        CHECK_PNC_ERROR( ncmpi_inq_vardimid(ncid, acid, acdimids) )
      }
      else {
        CHECK_ERROR( nc_inq_varid(ncid, dims[i], &acid) )
        CHECK_ERROR( nc_inq_varndims(ncid, acid, &m) )
        acdimids = (int*) malloc((sizeof(int)*m));
        CHECK_ERROR( nc_inq_vardimid(ncid, acid, acdimids) )
      }
      for (j=m-1; j>=0; j--) {
        new = 1;
        for (k=0; k<ndims; k++) {
          if (acdimids[j] == dimids[k]) {
            new = 0;
            break;
          }
        }
        if (new) {
          if ( nndims >= ndims ) {
            fprintf(stderr, "Error: invalid associated coordinates\n");
            return ERROR_CODE;
          }
          dimids[ndims-(++nndims)] = acdimids[j];
          //nc_inq_dimname(ncid, acdimids[j], tname);
          //printf("add %s\n", tname);
        }
      }
      free(acdimids);
      has_assoc = 1;
    }
  }
  if (nndims != n) {
    fprintf(stderr, "Error: invalid associated coordinates: %d %d\n", ndims, nndims);
    return ERROR_CODE;
  }

  TYPE2NCTYPE(dtype, xtype);
  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_def_var(ncid, varname, xtype, ndims, dimids, &varid) )
    // put variable attribute
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "units", strlen(units), units) )
//    CHECK_PNC_ERROR( ncmpi_put_att_double(ncid, varid, _FillValue, xtype, 1, &rmiss) )
    CHECK_PNC_ERROR( ncmpi_put_att_double(ncid, varid, "missing_value", xtype, 1, &rmiss) )
  }
  else {
    CHECK_ERROR( nc_def_var(ncid, varname, xtype, ndims, dimids, &varid) )
    // put variable attribute
    CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) )
    CHECK_ERROR( nc_put_att_double(ncid, varid, _FillValue, xtype, 1, &rmiss) )
    CHECK_ERROR( nc_put_att_double(ncid, varid, "missing_value", xtype, 1, &rmiss) )
  }
  if ( has_assoc ) {
    strcpy(coord, dims[0]);
    for(i=1; i<n; i++) {
      if (strlen(coord)+strlen(dims[i])+1 < File_HMID) {
        strcat(coord, " ");
        strcat(coord, dims[i]);
      }
    }
    if ( ndims > n && strlen(coord)+6 < File_HMID) {
      strcat(coord, " ");
      strcat(coord, vars[nvar]->t->name);
    }
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "coordinates", strlen(coord), coord) )
    else
      CHECK_ERROR( nc_put_att_text(ncid, varid, "coordinates", strlen(coord), coord) )
  }


  if ( tavg ) {
    sprintf(buf, "%s: mean", vars[nvar]->t->name);
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "cell_methods", strlen(buf), buf) )
    else
      CHECK_ERROR( nc_put_att_text(ncid, varid, "cell_methods", strlen(buf), buf) )
  }

  // set start and count
  vars[nvar]->start = (size_t*) malloc(sizeof(size_t)*ndims);
  vars[nvar]->count = (size_t*) malloc(sizeof(size_t)*ndims);
  for ( i=0; i<ndims; i++ ) {
    if ( files[fid]->shared_mode ) {
      MPI_Offset dimlen;
      CHECK_PNC_ERROR( ncmpi_inq_dimlen(ncid, dimids[i], &dimlen) )
      size = (size_t) dimlen;
    }
    else
      CHECK_ERROR( nc_inq_dimlen(ncid, dimids[i], &size) )
    vars[nvar]->count[i] = size;
    vars[nvar]->start[i] = 0;
  }
  if ( tint > 0.0 ) vars[nvar]->count[0] = 1;

#ifndef NETCDF3
  // set chunk size and deflate level (NetCDF-4 only)
  if ( ! files[fid]->shared_mode && files[fid]->deflate_level > 0 ) {
    CHECK_ERROR( nc_def_var_chunking(ncid, varid, NC_CHUNKED, vars[nvar]->count) )
    CHECK_ERROR( nc_def_var_deflate(ncid, varid, 0, 1, files[fid]->deflate_level) )
  }
#endif

  vars[nvar]->varid = varid;
  *vid = nvar;
  nvar++;

  return SUCCESS_CODE;
}

int32_t file_enddef( int32_t fid ) // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

#if defined(NETCDF3) || defined(PNETCDF)
  if (files[fid]->defmode == 1) {
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_enddef(ncid) )
    else
      CHECK_ERROR( nc_enddef(ncid) )
    files[fid]->defmode = 0;
  }
#endif

  return SUCCESS_CODE;
}

int32_t file_attach_buffer( int32_t fid,
                            int32_t buf_amount ) // (in)
{
  int ncid;
  MPI_Offset buf_amount_ = buf_amount;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_buffer_attach(ncid, buf_amount_) )

  return SUCCESS_CODE;
}

int32_t file_detach_buffer( int32_t fid ) // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_buffer_detach(ncid) )

  return SUCCESS_CODE;
}

int32_t file_flush( int32_t fid ) // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_wait_all(ncid, NC_REQ_ALL, NULL, NULL) )
  else
    CHECK_ERROR( nc_sync(ncid) )

  return SUCCESS_CODE;
}

int32_t file_write_data( int32_t     fid,        // (in)
                         int32_t     vid,        // (in)
                         void       *var,        // (in)
                         real64_t    t_start,    // (in)
                         real64_t    t_end,      // (in)
                         int32_t     precision,  // (in)
                         MPI_Offset *start,      // (in)
                         MPI_Offset *count)      // (in)
{
  int ncid, varid;
  if ( vars[vid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = vars[vid]->ncid;

#ifdef NETCDF3
  if ( files[fid]->defmode == 1 && ! files[fid]->shared_mode ) {
    CHECK_ERROR( nc_enddef(ncid) )
    files[fid]->defmode = 0;
  }
#endif

  varid = vars[vid]->varid;
  if ( vars[vid]->t != NULL ) { // have time dimension
    if ( vars[vid]->t->count < 0 ||  // first time
         t_end > vars[vid]->t->t + TEPS ) { // time goes next step
      vars[vid]->t->count += 1;
      vars[vid]->t->t = t_end;
      if ( vars[vid]->t->count > NTMAX-1 ) {
        fprintf(stderr, "time count exceeds the max limit (%d)\n", NTMAX);
        return ERROR_CODE;
      }
      vars[vid]->t->tval[vars[vid]->t->count] = t_end;
      if ( files[fid]->shared_mode ) { // write a new value to variable time
        MPI_Offset index[2];
        index[0] = (MPI_Offset) vars[vid]->t->count;
        CHECK_PNC_ERROR( ncmpi_put_var1_double_all(ncid, vars[vid]->t->varid, index, &t_end) )
        index[1] = 0;
        CHECK_PNC_ERROR( ncmpi_put_var1_double_all(ncid, vars[vid]->t->bndsid, index, &t_start ) )
        index[1] = 1;
        CHECK_PNC_ERROR( ncmpi_put_var1_double_all(ncid, vars[vid]->t->bndsid, index, &t_end ) )
      } else {
        size_t index[2];
        index[0] = vars[vid]->t->count;
        CHECK_ERROR( nc_put_var1_double(ncid, vars[vid]->t->varid, index, &t_end) )
        index[1] = 0;
        CHECK_ERROR( nc_put_var1_double(ncid, vars[vid]->t->bndsid, index, &t_start) )
        index[1] = 1;
        CHECK_ERROR( nc_put_var1_double(ncid, vars[vid]->t->bndsid, index, &t_end) )
      }
      vars[vid]->start[0] = vars[vid]->t->count;
    } else {
      size_t nt = vars[vid]->t->count + 1;
      int flag, n;
      flag = 1;
      for(n=nt-1;n>=0;n--) {
        if ( fabs(vars[vid]->t->tval[n]-t_end) < TEPS ) {
          vars[vid]->start[0] = n;
          flag = 0;
          break;
        }
      }
      if ( flag ) {
        fprintf(stderr, "cannot find time: %f\n", t_end);
        fprintf(stderr, "  time count is : %d, last time is: %f, diff is: %e\n", vars[vid]->t->count < 0, vars[vid]->t->t, vars[vid]->t->t-t_end);
        fprintf(stderr, "  time is: ");
        for (n=0;n<nt;n++) fprintf(stderr, "%f, ", vars[vid]->t->tval[n]);
        fprintf(stderr, "\n");
        return ERROR_CODE;
      }
    }
    if ( files[fid]->shared_mode ) {
      // add time dimension to start[0] and count[0]
      int i, ndims;
      CHECK_PNC_ERROR( ncmpi_inq_varndims(ncid, varid, &ndims) )
      for (i=ndims-1; i>0; i--) {
        start[i] = start[i-1];
        count[i] = count[i-1];
      }
      start[0] = vars[vid]->start[0];  // start along the time dimension
      count[0] = vars[vid]->count[0];
    }
  }

  switch (precision) {
  case 8:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_bput_vara_double(ncid, varid, start, count, (double*)var, NULL) )
    else
      CHECK_ERROR( nc_put_vara_double(ncid, varid, vars[vid]->start, vars[vid]->count, (double*)var) )
    break;
  case 4:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_bput_vara_float(ncid, varid, start, count, (float*)var, NULL) )
    else
      CHECK_ERROR( nc_put_vara_float(ncid, varid, vars[vid]->start, vars[vid]->count, (float*)var) )
    break;
  default:
    fprintf(stderr, "unsupported data precision: %d\n", precision);
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_close( int32_t fid ) // (in)
{
  int ncid;
  int i;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  for (i=0; i<nvar; i++) {
    if ( vars[i] != NULL && vars[i]->ncid == ncid ) {
      free( vars[i]->start );
      free( vars[i]->count );
      free( vars[i] );
      vars[i] = NULL;
    }
  }

  for (i=0; i<nt; i++) {
    if ( tdims[i] != NULL && tdims[i]->ncid == ncid ) {
      free( tdims[i]->tval );
      free( tdims[i] );
      tdims[i] = NULL;
    }
  }

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_close(ncid) )
  else
    CHECK_ERROR( nc_close(ncid) )

  free( files[fid] );
  files[fid] = NULL;

  return SUCCESS_CODE;
}
