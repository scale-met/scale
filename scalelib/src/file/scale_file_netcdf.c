#include "scale_file.h"
#ifndef MPI_INCLUDED
#define MPI_INCLUDED
#endif
#include "netcdf.h"

#define TEPS 1e-6
#define NTMAX 102400

#define MIN(a,b) ((a)<(b) ? (a) : (b))

static int32_t ERROR_SUPPRESS = 0;

#define CHECK_ERROR(func)                                       \
  {                                                             \
    int status_ = (func);                                       \
    if (status_ != NC_NOERR) {                                  \
      if ( ! ERROR_SUPPRESS ) {                                 \
        fprintf(stderr, "Error: at line %d in %s\n", __LINE__, __FILE__);   \
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
        fprintf(stderr, "Error: at line %d in %s\n", __LINE__, __FILE__);   \
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
  size_t ndims;
  size_t ndims_t;
} varinfo_t;

static fileinfo_t *files[FILE_MAX];
static int nfile = 0;
static varinfo_t *vars[VAR_MAX];
static int nvar = 0;
static tdim_t *tdims[VAR_MAX];
static int nt = 0;


int32_t file_open_c(       int32_t  *fid,     // (out)
		     const char     *fname,   // (in)
		     const int32_t   mode,    // (in)
		     const MPI_Comm  comm )   // (in)
{
  int ncid;
  int len;
  int shared_mode;
  char _fname[File_HLONG+4];
  int add_suffix;

  if ( nfile >= FILE_MAX ) {
    fprintf(stderr, "exceed max number of file limit\n");
    return ERROR_CODE;
  }

  len = strlen(fname);
  strcpy(_fname, fname);

  if ( mode==File_FREAD || mode==File_FAPPEND ) {
    FILE *fp = fopen(_fname, "r");
    if ( fp==NULL ) {
      add_suffix = 1;
    } else {
      fclose(fp);
      add_suffix = 0;
    }
  } else
    add_suffix = 1;

  if ( add_suffix )
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

int32_t file_get_dim_length_c( const int32_t  fid,      // (in)
			       const char*    dimname,  // (in)
			             int32_t *len     ) // (out)
{
  int ncid, dimid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    MPI_Offset l;
    CHECK_PNC_ERROR( ncmpi_inq_dimid(ncid, dimname, &dimid) )
    CHECK_PNC_ERROR( ncmpi_inq_dimlen(ncid, dimid, &l) )
    *len = l;
  } else {
    size_t l;
    CHECK_ERROR( nc_inq_dimid(ncid, dimname, &dimid) )
    CHECK_ERROR( nc_inq_dimlen(ncid, dimid, &l) )
    *len = l;
  }

  return SUCCESS_CODE;
}

int32_t file_set_option_c( const int32_t fid,    // (in)
			   const char* filetype, // (in)
			   const char* key,      // (in)
			   const char* val)      // (in)
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

int32_t file_get_nvars_c( const int32_t  fid,   // (in)
			        int32_t *nvars )// (out)
{
  int ncid;
  int ndims, ngatts, unlimdim;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_inq(ncid, &ndims, nvars, &ngatts, &unlimdim) )
  else
    CHECK_ERROR( nc_inq(ncid, &ndims, nvars, &ngatts, &unlimdim) )

  return SUCCESS_CODE;
}

int32_t file_get_varname_c( const int32_t  fid,  // (in)
			    const int32_t  vid,  // (in)
			          char    *name, // (out)
			    const int32_t  len ) // (in)
{
  int ncid, varid;
  char buf[MAX_NC_NAME+1];
  int i;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid  = files[fid]->ncid;
  varid = vid-1; // index starts from 1 in fortran space

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_inq_varname(ncid, varid, buf) )
  else
    CHECK_ERROR( nc_inq_varname(ncid, varid, buf) )

  for (i=0; i<MIN(len-1,strlen(buf)); i++)
    name[i] = buf[i];
  name[i] = '\0';

  return SUCCESS_CODE;
}

int32_t file_get_datainfo_c(       datainfo_t *dinfo,   // (out)
			     const int32_t     fid,     // (in)
			     const char*       varname, // (in)
			     const int32_t     step,    // (in)
			     const int32_t     suppress)// (in)
{
  int ncid, varid;
  nc_type xtype;
  int rank;
  int dimids[RANK_MAX], tdim, uldims[NC_MAX_DIMS];
  char name[NC_MAX_NAME+1];
  char *buf;
  size_t size;
  int i, n;
  int status;

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
    CHECK_PNC_ERROR( ncmpi_inq_attlen(ncid, varid, "long_name", &l) )
    buf = (char*) malloc(l+1);
    CHECK_PNC_ERROR( ncmpi_get_att_text(ncid, varid, "long_name", buf) )
    for (i=0; i<MIN(File_HMID-1,l); i++)
      dinfo->description[i] = buf[i];
    dinfo->description[i] = '\0';
    free(buf);
    // units
    CHECK_PNC_ERROR( ncmpi_inq_attlen(ncid, varid, "units", &l) )
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
    status = nc_inq_attlen(ncid, varid, "long_name", &l);
    if ( status == NC_NOERR ) {
      buf = (char*) malloc(l+1);
      CHECK_ERROR( nc_get_att_text(ncid, varid, "long_name", buf) )
    } else { // for WRF file
      CHECK_ERROR( nc_inq_attlen(ncid, varid, "description", &l) )
      buf = (char*) malloc(l+1);
      CHECK_ERROR( nc_get_att_text(ncid, varid, "description", buf) )
    }
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
  if (rank > RANK_MAX) {
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
      status = nc_inq_varid(ncid, name, &varid);
      if ( status == NC_NOERR ) {
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
      } else {
	dinfo->time_start = 0.0;
	dinfo->time_end = 0.0;
	dinfo->time_units[0] = '\0';
      }
    }
  } else {
    if ( step > 1 ) { // if variable does not have time dimention, step > 1 should not exist
      fprintf(stderr, "requested step is larger than tdim: step=%d tdim=%d\n", step, tdim);
      return ERROR_CODE;
    }
  }
  ERROR_SUPPRESS = 0;

  return SUCCESS_CODE;
}

int32_t file_read_data_c(       void       *var,       // (out)
			  const datainfo_t *dinfo,     // (in)
			  const int32_t     precision, // (in)
			  const MPI_Offset  ntypes,    // (in)
			  const MPI_Datatype dtype,    // (in)
			  const int32_t     *start,    // (in)
			  const int32_t     *count )   // (in)
{
  int ncid, varid;
  int rank;
  int i;
  int fid;
  size_t *str, *cnt;
  MPI_Offset *strp, *cntp;
  size_t size;
  int l_rescale;

  fid = dinfo->fid;
  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;
  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, dinfo->varname, &varid) )
    CHECK_PNC_ERROR( ncmpi_inq_varndims(ncid, varid, &rank) )
    strp = (MPI_Offset*) malloc(sizeof(MPI_Offset)*rank);
    cntp = (MPI_Offset*) malloc(sizeof(MPI_Offset)*rank);
  } else {
    CHECK_ERROR( nc_inq_varid(ncid, dinfo->varname, &varid) )
    CHECK_ERROR( nc_inq_varndims(ncid, varid, &rank) )
  }
  str = (size_t*) malloc(sizeof(size_t)*rank);
  cnt = (size_t*) malloc(sizeof(size_t)*rank);

  if ( start == NULL || count == NULL ) {
    for (i=0; i<dinfo->rank; i++) {
      // note: C and Fortran orders are opposite
      str[rank -i-1] = 0;
      cnt[rank -i-1] = dinfo->dim_size[i];
    }
  } else {
    for (i=0; i<dinfo->rank; i++) {
      // note: C and Fortran orders are opposite
      str[rank -i-1] = start[i] - 1;
      cnt[rank -i-1] = count[i];
    }
  }
  if (rank > dinfo->rank) { // have time dimension
    str[0] = dinfo->step - 1;
    cnt[0] = 1;
  }

  size = 1;
  for (i=0; i<rank; i++) size *= cnt[i];

  if ( files[fid]->shared_mode ) {
    for (i=0; i<rank; i++) {
      strp[i] = (MPI_Offset) str[i];
      cntp[i] = (MPI_Offset) cnt[i];
    }
    free(str);
    free(cnt);
    CHECK_PNC_ERROR( ncmpi_iget_vara(ncid, varid, strp, cntp, var, ntypes, dtype, NULL) )
    free(strp);
    free(cntp);
    if ( dtype == MPI_FLOAT ) {
      float factor, offset;
      l_rescale = 0;
      if ( ncmpi_get_att_float(ncid, varid, "scale_factor", &factor) != NC_NOERR )
	factor = 1.0f;
      else
	l_rescale = 1;
      if ( ncmpi_get_att_float(ncid, varid, "add_offset", &offset) != NC_NOERR )
	offset = 0.0f;
      else
	l_rescale = 1;
      if ( l_rescale ) for (i=0; i<size; i++) ((float*)var)[i] = ((float*)var)[i] * factor + offset;
    } else if ( dtype == MPI_DOUBLE ) {
      double factor, offset;
      l_rescale = 0;
      if ( ncmpi_get_att_double(ncid, varid, "scale_factor", &factor) != NC_NOERR )
	factor = 1.0;
      else
	l_rescale = 1;
      if ( ncmpi_get_att_double(ncid, varid, "add_offset", &offset) != NC_NOERR )
	offset = 0.0;
      else
	l_rescale = 1;
      if ( l_rescale ) for (i=0; i<size; i++) ((double*)var)[i] = ((double*)var)[i] * factor + offset;
    } else {
      float factor, offset;
      if (    ( ncmpi_get_att_float(ncid, varid, "scale_factor", &factor) == NC_NOERR ) 
           || ( ncmpi_get_att_float(ncid, varid, "add_offset",   &offset) == NC_NOERR ) ) {
	fprintf(stderr, "scale_factor and add_offset is not supported with a MPI derived type\n");
	return ERROR_CODE;
      }
    }
  } else {
    switch ( precision ) {
    case 8:
      CHECK_ERROR( nc_get_vara_double(ncid, varid, str, cnt, (double*)var) )
      {
	double factor, offset;
	l_rescale = 0;
	if ( nc_get_att_double(ncid, varid, "scale_factor", &factor) != NC_NOERR )
	  factor = 1.0;
	else
	  l_rescale = 1;
	if ( nc_get_att_double(ncid, varid, "add_offset", &offset) != NC_NOERR )
	  offset = 0.0;
	else
	  l_rescale = 1;
	if ( l_rescale ) for (i=0; i<size; i++) ((double*)var)[i] = ((double*)var)[i] * factor + offset;
      }
      break;
    case 4:
      CHECK_ERROR( nc_get_vara_float(ncid, varid, str, cnt, (float*)var) )
      {
	float factor, offset;
	l_rescale = 0;
	if ( nc_get_att_float(ncid, varid, "scale_factor", &factor) != NC_NOERR )
	  factor = 1.0f;
	else
	  l_rescale = 1;
	if ( nc_get_att_float(ncid, varid, "add_offset", &offset) != NC_NOERR )
	  offset = 0.0f;
	else
	  l_rescale = 1;
	if ( l_rescale ) for (i=0; i<size; i++) ((float*)var)[i] = ((float*)var)[i] * factor + offset;
      }
      break;
    default:
      free(str);
      free(cnt);
      fprintf(stderr, "unsupported data precision: %d\n", precision );
      return ERROR_CODE;
    }
    free(str);
    free(cnt);
  }


  return SUCCESS_CODE;
}

int32_t file_get_attribute_text_c( const int32_t  fid,   // (in)
				   const char    *vname, // (in)
				   const char    *key,   // (in)
				         char    *value, // (out)
				   const int32_t len)    // (in)
{
  int ncid;
  int varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    MPI_Offset l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    CHECK_PNC_ERROR( ncmpi_inq_attlen(ncid, varid, key, &l) )
    if ( len < l+1 ) return ERROR_CODE;

    CHECK_PNC_ERROR( ncmpi_get_att_text(ncid, varid, key, value) )
    value[l] = '\0';
  }
  else {
    size_t l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    CHECK_ERROR( nc_inq_attlen(ncid, varid, key, &l) )
    if ( len < l+1 ) return ERROR_CODE;

    CHECK_ERROR( nc_get_att_text(ncid, varid, key, value) )
    value[l] = '\0';
  }

  return SUCCESS_CODE;
}

int32_t file_get_attribute_int_c( const int32_t  fid,   // (in)
				  const char    *vname, // (in)
				  const char    *key,   // (in)
				        int     *value, // (out)
				  const size_t   len)   // (in)
{
  int ncid;
  int varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    MPI_Offset l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    CHECK_PNC_ERROR( ncmpi_inq_attlen(ncid, varid, key, &l) )
    if ( len < l ) return ERROR_CODE;
    CHECK_PNC_ERROR( ncmpi_get_att_int(ncid, varid, key, value) )
  }
  else {
    size_t l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    CHECK_ERROR( nc_inq_attlen(ncid, varid, key, &l) )
    if ( len < l ) return ERROR_CODE;
    CHECK_ERROR( nc_get_att_int(ncid, varid, key, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_get_attribute_float_c( const int32_t  fid,   // (in)
				    const char    *vname, // (in)
				    const char    *key,   // (in)
				          float   *value, // (out)
				    const size_t   len)   // (in)
{
  int ncid;
  int varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    MPI_Offset l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    CHECK_PNC_ERROR( ncmpi_inq_attlen(ncid, varid, key, &l) )
    if ( len < l ) return ERROR_CODE;
    CHECK_PNC_ERROR( ncmpi_get_att_float(ncid, varid, key, value) )
  }
  else {
    size_t l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    CHECK_ERROR( nc_inq_attlen(ncid, varid, key, &l) )
    if ( len < l ) return ERROR_CODE;
    CHECK_ERROR( nc_get_att_float(ncid, varid, key, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_get_attribute_double_c( const int32_t  fid,   // (in)
				     const char    *vname, // (in)
				     const char    *key,   // (in)
				           double  *value, // (out)
				     const size_t   len)   // (in)
{
  int ncid;
  int varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    MPI_Offset l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    CHECK_PNC_ERROR( ncmpi_inq_attlen(ncid, varid, key, &l) )
    if ( len < l ) return ERROR_CODE;
    CHECK_PNC_ERROR( ncmpi_get_att_double(ncid, varid, key, value) )
  }
  else {
    size_t l;
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    CHECK_ERROR( nc_inq_attlen(ncid, varid, key, &l) )
    if ( len < l ) return ERROR_CODE;
    CHECK_ERROR( nc_get_att_double(ncid, varid, key, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_attribute_text_c( const int32_t  fid,   // (in)
				   const char    *vname, // (in)
				   const char    *key,   // (in)
				   const char    *val)   // (in)
{
  int ncid;
  int varid;
  int attid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    if ( ncmpi_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, key, strlen(val), val) )
  }
  else {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
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

int32_t file_set_attribute_int_c( const int32_t  fid,   // (in)
				  const char    *vname, // (in)
				  const char    *key,   // (in)
				  const int32_t *value, // (in)
				  const size_t   len )  // (in)
{
  int ncid;
  int varid;
  int attid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    if ( ncmpi_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    CHECK_PNC_ERROR( ncmpi_put_att_int(ncid, varid, key, NC_INT, len, value) )
  }
  else {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    if ( nc_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
    if (files[fid]->defmode == 0) {
      CHECK_ERROR( nc_redef(ncid) )
      files[fid]->defmode = 1;
    }
#endif

    CHECK_ERROR( nc_put_att_int(ncid, varid, key, NC_INT, len, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_attribute_float_c( const int32_t  fid,   // (in)
				    const char    *vname, // (in)
				    const char    *key,   // (in)
				    const float   *value, // (in)
				    const size_t   len )  // (in)
{
  int ncid;
  int varid;
  int attid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    if ( ncmpi_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    CHECK_PNC_ERROR( ncmpi_put_att_float(ncid, varid, key, NC_FLOAT, len, value) )
  }
  else {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    if ( nc_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
    if (files[fid]->defmode == 0) {
      CHECK_ERROR( nc_redef(ncid) )
      files[fid]->defmode = 1;
    }
#endif

    CHECK_ERROR( nc_put_att_float(ncid, varid, key, NC_FLOAT, len, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_set_attribute_double_c( const int32_t  fid,   // (in)
				     const char    *vname, // (in)
				     const char    *key,   // (in)
				     const double  *value, // (in)
				     const size_t   len )  // (in)
{
  int ncid;
  int varid;
  int attid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode ) {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_PNC_ERROR( ncmpi_inq_varid(ncid, vname, &varid) )

    if ( ncmpi_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

    CHECK_PNC_ERROR( ncmpi_put_att_double(ncid, varid, key, NC_DOUBLE, len, value) )
  }
  else {
    if ( strcmp(vname, "global") == 0 ) {
      varid = NC_GLOBAL;
    } else
      CHECK_ERROR( nc_inq_varid(ncid, vname, &varid) )

    if ( nc_inq_attid(ncid, varid, key, &attid) == NC_NOERR ) // check if existed
      return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
    if (files[fid]->defmode == 0) {
      CHECK_ERROR( nc_redef(ncid) )
      files[fid]->defmode = 1;
    }
#endif

    CHECK_ERROR( nc_put_att_double(ncid, varid, key, NC_DOUBLE, len, value) )
  }

  return SUCCESS_CODE;
}

int32_t file_add_associatedvariable_c( const int32_t  fid,   // (in)
				       const char    *vname) // (in)
{
  int ncid, varid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( nc_inq_varid(ncid, vname, &varid) == NC_NOERR ) // check if existed
    return ALREADY_EXISTED_CODE;

#ifdef NETCDF3
  if (files[fid]->defmode == 0) {
    CHECK_ERROR( nc_redef(ncid) )
    files[fid]->defmode = 1;
  }
#endif

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_def_var(ncid, vname, NC_INT, 0, 0, &varid) )
  else
    CHECK_ERROR( nc_def_var(ncid, vname, NC_INT, 0, 0, &varid) )

#ifdef NETCDF3
  CHECK_ERROR( nc_enddef(ncid) )
  files[fid]->defmode = 0;
#endif

  return SUCCESS_CODE;
}

int32_t file_set_tunits_c( const int32_t fid,         // (in)
			   const char    *time_units) // (in)
{
  strcpy(files[fid]->time_units, time_units);

  return SUCCESS_CODE;
}

int32_t file_put_axis_c( const int32_t fid,        // (in)
			 const char   *name,       // (in)
			 const char   *desc,       // (in)
			 const char   *units,      // (in)
			 const char   *dim_name,   // (in)
			 const int32_t dtype,      // (in)
			 const void*   val,        // (in)
			 const int32_t size,       // (in)
			 const int32_t precision)  // (in)
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

int32_t file_def_axis_c( const int32_t fid,        // (in)
			 const char   *name,       // (in)
			 const char   *desc,       // (in)
			 const char   *units,      // (in)
			 const char   *dim_name,   // (in)
			 const int32_t dtype,      // (in)
			 const int32_t dim_size)   // (in)
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

int32_t file_write_axis_c( const int32_t     fid,       // (in)
			   const char       *name,      // (in)
			   const void       *val,       // (in)
			   const int32_t     precision, // (in)
			   const MPI_Offset *start,     // (in)
			   const MPI_Offset *count)     // (in)
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

int32_t file_put_associatedcoordinate_c( const int32_t fid,        // (in)
					 const char   *name,       // (in)
					 const char   *desc,       // (in)
					 const char   *units,      // (in)
					 const char   **dim_names, // (in)
					 const int32_t ndims,      // (in)
					 const int32_t dtype,      // (in)
					 const void*   val,        // (in)
					 const int32_t precision)  // (in)
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

int32_t file_def_associatedcoordinate_c( const int32_t fid,        // (in)
					 const char   *name,       // (in)
					 const char   *desc,       // (in)
					 const char   *units,      // (in)
					 const char   **dim_names, // (in)
					 const int32_t ndims,      // (in)
					 const int32_t dtype)      // (in)
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

int32_t file_write_associatedcoordinate_c( const int32_t     fid,        // (in)
					   const char       *name,       // (in)
					   const void*       val,        // (in)
					   const int32_t     precision,  // (in)
					   const MPI_Offset *start,      // (in)
					   const MPI_Offset *count)      // (in)
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

int32_t file_add_variable_c(       int32_t *vid,     // (out)
			     const int32_t  fid,     // (in)
			     const char    *varname, // (in)
			     const char    *desc,    // (in)
			     const char    *units,   // (in)
			     const char   **dims,    // (in)
			     const int32_t  ndims,   // (in)
			     const int32_t  dtype,   // (in)
			     const real64_t tint,    // (in)
			     const int32_t  tavg)    // (in)
{
  int ncid, varid, acid, *acdimids;
  int dimids[NC_MAX_DIMS], dimid;
  char tname[File_HSHORT+1];
  int tdimid, tvarid;
  nc_type xtype = -1;
  char buf[File_HMID+1];
  int i, j, k, m, err;
  int ndims_t, nndims;
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
  vars[nvar]->ndims = ndims;

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
      if ( tdims[i]       != NULL &&  // still opened
           tdims[i]->ncid == ncid &&  // same file
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
      m=0;
      for (i=0; i<nt; i++) {
        if ( tdims[i] != NULL && tdims[i]->ncid == ncid ) m++;
      }
      if ( m == 0 ) {
        strcpy(tname, "time");
      } else {
        sprintf(tname, "time%d", m);
      }
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
  ndims_t = ndims;
  if ( tint > 0.0 ) { // add time dimension
    dimids[0] = vars[nvar]->t->dimid;
    ndims_t++;
  }
  for (i=ndims_t-ndims; i<ndims_t; i++) dimids[i] = -1;

  has_assoc = 0;
  nndims = 0;
  for (i=0; i<ndims; i++) {
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
        dimids[ndims_t-(++nndims)] = dimid;
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
        for (k=0; k<ndims_t; k++) {
          if (acdimids[j] == dimids[k]) {
            new = 0;
            break;
          }
        }
        if (new) {
          if ( nndims >= ndims_t ) {
            fprintf(stderr, "Error: invalid associated coordinates\n");
            return ERROR_CODE;
          }
          dimids[ndims_t-(++nndims)] = acdimids[j];
          //nc_inq_dimname(ncid, acdimids[j], tname);
          //printf("add %s\n", tname);
        }
      }
      free(acdimids);
      has_assoc = 1;
    }
  }
  if (nndims != ndims) {
    fprintf(stderr, "Error: invalid associated coordinates: %d %d\n", ndims_t, nndims);
    return ERROR_CODE;
  }

  TYPE2NCTYPE(dtype, xtype);
  if ( files[fid]->shared_mode ) {
    CHECK_PNC_ERROR( ncmpi_def_var(ncid, varname, xtype, ndims_t, dimids, &varid) )
    // put variable attribute
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_PNC_ERROR( ncmpi_put_att_text(ncid, varid, "units", strlen(units), units) )
//    CHECK_PNC_ERROR( ncmpi_put_att_double(ncid, varid, _FillValue, xtype, 1, &rmiss) )
    CHECK_PNC_ERROR( ncmpi_put_att_double(ncid, varid, "missing_value", xtype, 1, &rmiss) )
  }
  else {
    CHECK_ERROR( nc_def_var(ncid, varname, xtype, ndims_t, dimids, &varid) )
    // put variable attribute
    CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) )
    CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) )
    CHECK_ERROR( nc_put_att_double(ncid, varid, _FillValue, xtype, 1, &rmiss) )
    CHECK_ERROR( nc_put_att_double(ncid, varid, "missing_value", xtype, 1, &rmiss) )
  }
  if ( has_assoc ) {
    strcpy(coord, dims[0]);
    for(i=1; i<ndims; i++) {
      if (strlen(coord)+strlen(dims[i])+1 < File_HMID) {
        strcat(coord, " ");
        strcat(coord, dims[i]);
      }
    }
    if ( ndims_t > ndims && strlen(coord)+6 < File_HMID) {
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
  vars[nvar]->ndims_t = ndims_t;
  vars[nvar]->start = (size_t*) malloc(sizeof(size_t)*ndims_t);
  vars[nvar]->count = (size_t*) malloc(sizeof(size_t)*ndims_t);
  for ( i=0; i<ndims_t; i++ ) {
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

int32_t file_enddef_c( const int32_t fid ) // (in)
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

int32_t file_attach_buffer_c( const int32_t fid,         // (in)
			      const int32_t buf_amount ) // (in)
{
  int ncid;
  MPI_Offset buf_amount_ = buf_amount;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_buffer_attach(ncid, buf_amount_) )

  return SUCCESS_CODE;
}

int32_t file_detach_buffer_c( const int32_t fid ) // (in)
{
  int ncid;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  if ( files[fid]->shared_mode )
    CHECK_PNC_ERROR( ncmpi_buffer_detach(ncid) )

  return SUCCESS_CODE;
}

int32_t file_flush_c( const int32_t fid ) // (in)
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

int32_t file_write_data_c( const int32_t   fid,       // (in)
			   const int32_t   vid,       // (in)
			   const void     *var,       // (in)
			   const real64_t  t_start,   // (in)
			   const real64_t  t_end,     // (in)
			   const int32_t   precision, // (in)
			   const int32_t   ndims,     // (in)
			   const int32_t  *start,     // (in)
			   const int32_t  *count)     // (in)
{
  int ncid, varid;
  MPI_Offset *str, *cnt;

  if ( vars[vid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = vars[vid]->ncid;

  if ( ndims != vars[vid]->ndims ) {
    fprintf(stderr, "Error: at line %d in %s\n", __LINE__, __FILE__);
    fprintf(stderr, "       dimension size %d is not consistent that was added by file_add_variable %d\n", ndims, vars[vid]->ndims );
    return ERROR_CODE;
  }

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
  }

  if ( files[fid]->shared_mode ) {
    int i;
    int ndims_t = vars[vid]->ndims_t;
    str = (MPI_Offset*) malloc(sizeof(MPI_Offset)*(ndims_t));
    cnt = (MPI_Offset*) malloc(sizeof(MPI_Offset)*(ndims_t));
    if ( vars[vid]->t != NULL ) { // have time dimension
      // add time dimension to start[0] and count[0]
      str[0] = vars[vid]->start[0];  // start along the time dimension
      cnt[0] = vars[vid]->count[0];
      for (i=0; i<ndims; i++) {
        str[ndims_t-i-1] = start[i] - 1;
        cnt[ndims_t-i-1] = count[i];
      }
    } else {
      for (i=0; i<ndims; i++) {
        str[ndims-i-1] = start[i] - 1;
        cnt[ndims-i-1] = count[i];
      }
    }
  }

  switch (precision) {
  case 8:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_bput_vara_double(ncid, varid, str, cnt, (double*)var, NULL) )
    else
      CHECK_ERROR( nc_put_vara_double(ncid, varid, vars[vid]->start, vars[vid]->count, (double*)var) )
    break;
  case 4:
    if ( files[fid]->shared_mode )
      CHECK_PNC_ERROR( ncmpi_bput_vara_float(ncid, varid, str, cnt, (float*)var, NULL) )
    else
      CHECK_ERROR( nc_put_vara_float(ncid, varid, vars[vid]->start, vars[vid]->count, (float*)var) )
    break;
  default:
    fprintf(stderr, "unsupported data precision: %d\n", precision);
    return ERROR_CODE;
  }

  if ( files[fid]->shared_mode) {
    free(str);
    free(cnt);
  }

  return SUCCESS_CODE;
}

int32_t file_close_c( const int32_t fid ) // (in)
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
