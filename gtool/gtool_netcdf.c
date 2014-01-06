#include "netcdf.h"
#include "gtool_file.h"

#define RMISS -9.9999e+30
#define EPS 1e-10

#define CHECK_ERROR(status)					\
  {								\
    if (status != NC_NOERR) {					\
      fprintf(stderr, "Error: %s\n", nc_strerror(status));	\
      return ERROR_CODE;					\
    }								\
  }

#define NCTYPE2TYPE(nctype, type)				\
  {								\
  switch ( nctype ) {						\
  case NC_FLOAT:						\
    type = File_REAL4;						\
    break;							\
  case NC_DOUBLE:						\
    type = File_REAL8;						\
    break;							\
  default:                                                      \
    fprintf(stderr, "unsuppoted data type: %d\n", xtype);	\
    return ERROR_CODE;						\
  }								\
  }

#define TYPE2NCTYPE(type, nctype)				\
  {								\
  switch ( type ) {						\
  case File_REAL4:						\
    nctype = NC_FLOAT;						\
    break;							\
  case File_REAL8:						\
    nctype = NC_DOUBLE;						\
    break;							\
  default:							\
    fprintf(stderr, "unsuppoted data type: %d\n", xtype);	\
    return ERROR_CODE;						\
  }								\
  }


#define DEFAULT_DEFLATE_LEVEL 2

typedef struct {
  int ncid;
  char time_units[File_HSHORT+1];
  int deflate_level;
} fileinfo_t;

typedef struct {
  int ncid;
  int dimid;
  int varid;
  int bndsid;
  int count;
  real64_t t;
  real64_t tint;
  char name[File_HSHORT+1];
} tdim_t;

typedef struct {
  int ncid;
  int varid;
  tdim_t *t;
  size_t *start;
  size_t *count;
} varinfo_t;


#define FILE_MAX 64
#define VAR_MAX 64000

static fileinfo_t *files[FILE_MAX];
static int nfile = 0;
static varinfo_t *vars[VAR_MAX];
static int nvar = 0;
static tdim_t *tdims[VAR_MAX];
static int nt = 0;


int32_t file_open( int32_t *fid,     // (out)
		   char    *fname,   // (in)
		   int32_t  mode )   // (in)
{
  int ncid;
  int len;
  char _fname[File_HLONG+4];

  if ( nfile >= FILE_MAX ) {
    fprintf(stderr, "exceed max number of file limit\n");
    return ERROR_CODE;
  }

  len = strlen(fname);
  strcpy(_fname, fname);
  if (fname[len-3] != '.' || fname[len-2] != 'n' || fname[len-1] != 'c' )
    strcat(_fname, ".nc");

  switch ( mode ) {
  case File_FREAD:
    CHECK_ERROR( nc_open(_fname, NC_NOWRITE, &ncid) );
    break;
  case File_FWRITE:
#ifdef NETCDF3
    CHECK_ERROR( nc_create(_fname, NC_CLOBBER, &ncid) );
#else
    CHECK_ERROR( nc_create(_fname, NC_CLOBBER|NC_NETCDF4, &ncid) );
#endif
    break;
  case File_FAPPEND:
    CHECK_ERROR( nc_open(_fname, NC_WRITE, &ncid) );
    break;
  default:
    fprintf(stderr, "invalid mode type\n");
    return ERROR_CODE;
  }

  files[nfile] = (fileinfo_t*) malloc(sizeof(fileinfo_t));
  files[nfile]->ncid = ncid;
  files[nfile]->deflate_level = DEFAULT_DEFLATE_LEVEL;
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
			   int32_t     step)    // (in)
{
  int ncid, varid;
  nc_type xtype;
  int rank;
  int dimids[MAX_RANK], tdim, uldims[NC_MAX_DIMS];
  char name[File_HSHORT+1];
  size_t size;
  size_t idx[2];
  int i, n;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;
  CHECK_ERROR( nc_inq_varid(ncid, varname, &varid) );

  // fid
  dinfo->fid = fid;
  // varname
  strcpy(dinfo->varname, varname);
  // description
  CHECK_ERROR( nc_get_att_text(ncid, varid, "long_name", dinfo->description) );
  // units
  CHECK_ERROR( nc_get_att_text(ncid, varid, "units", dinfo->units) );
  // datatype
  CHECK_ERROR( nc_inq_vartype(ncid, varid, &xtype) );
  NCTYPE2TYPE(xtype, dinfo->datatype);
  // rank
  CHECK_ERROR( nc_inq_varndims(ncid, varid, &rank) );
  CHECK_ERROR( nc_inq_vardimid(ncid, varid, dimids) );
#ifdef NETCDF3
  CHECK_ERROR( nc_inq_unlimdim(ncid, &n) );
#else
  CHECK_ERROR( nc_inq_unlimdims(ncid, &n, uldims) );
#endif
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
    // note: C and Fortran orders are opposit
    CHECK_ERROR( nc_inq_dim(ncid, dimids[rank-i-1], name, &size) );
    strncpy(dinfo->dim_name+i*File_HSHORT, name, File_HSHORT);
    dinfo->dim_size[i] = size;
  }

  dinfo->step = step;
  if ( tdim >= 0 ) {
    // time_end
    CHECK_ERROR( nc_inq_dimname(ncid, tdim, name) );
    CHECK_ERROR( nc_inq_varid(ncid, name, &varid) );
    idx[0] = step;
    CHECK_ERROR( nc_get_var1_double(ncid, varid, idx, &(dinfo->time_end)) );
    // time_start
    strcat(name, "_bnds");
    CHECK_ERROR( nc_inq_varid(ncid, name, &varid) );
    idx[1] = 0;
    CHECK_ERROR( nc_get_var1_double(ncid, varid, idx, &(dinfo->time_start)) );
  } else {
  }

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
  CHECK_ERROR( nc_inq_varid(ncid, dinfo->varname, &varid) );

  CHECK_ERROR( nc_inq_varndims(ncid, varid, &rank) );
  start = (size_t*) malloc(sizeof(size_t)*rank);
  count = (size_t*) malloc(sizeof(size_t)*rank);
  for (i=0; i<dinfo->rank; i++) {
    // note: C and Fortran orders are opposit
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
    fprintf(stderr, "unsuppoted data precision: %d\n", precision );
    return ERROR_CODE;
  }

  return SUCCESS_CODE;
}

int32_t file_set_global_attributes( int32_t  fid,         // (in)
				    char    *title,       // (in)
				    char    *source,      // (in)
				    char    *institution, // (in)
				    char    *time_units,  // (in)
				    int32_t  nodeid,      // (in)
				    int32_t *nodeidx,     // (in)
				    int32_t  nodeidx_dim) // (in)
{
  int ncid;
  int tmp[1];

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  CHECK_ERROR( nc_put_att_text(ncid, NC_GLOBAL, "title", strlen(title), title) );
  CHECK_ERROR( nc_put_att_text(ncid, NC_GLOBAL, "source", strlen(source), source) );
  CHECK_ERROR( nc_put_att_text(ncid, NC_GLOBAL, "institution", strlen(institution), institution) );
  tmp[0] = nodeid;
  CHECK_ERROR( nc_put_att_int(ncid, NC_GLOBAL, "node_id", NC_INT, 1, tmp) );
  CHECK_ERROR( nc_put_att_int(ncid, NC_GLOBAL, "node_index", NC_INT, nodeidx_dim, nodeidx) );

  strcpy(files[fid]->time_units, time_units);

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
  CHECK_ERROR( nc_redef(ncid) );
#endif

  if ( nc_inq_dimid(ncid, dim_name, &dimid) != NC_NOERR ) // check if existed
    CHECK_ERROR( nc_def_dim(ncid, dim_name, size, &dimid) );

  TYPE2NCTYPE(dtype, xtype);
  CHECK_ERROR( nc_def_var(ncid, name, xtype, 1, &dimid, &varid) );
  CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) );
  CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) );

#ifdef NETCDF3
  CHECK_ERROR( nc_enddef(ncid) );
#endif

  switch ( precision ) {
  case 8:
    CHECK_ERROR( nc_put_var_double(ncid, varid, (double*)val) );
    break;
  case 4:
    CHECK_ERROR( nc_put_var_float(ncid, varid, (float*)val) );
    break;
  default:
    fprintf(stderr, "unsuppoted data precision: %d\n", precision);
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
  CHECK_ERROR( nc_redef(ncid) );
#endif

  dimids = malloc(sizeof(int)*ndims);
  for (i=0; i<ndims; i++)
    CHECK_ERROR( nc_inq_dimid(ncid, dim_names[i], dimids+i) );

  TYPE2NCTYPE(dtype, xtype);

  CHECK_ERROR( nc_def_var(ncid, name, xtype, ndims, dimids, &varid) );
  CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) );
  CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) );
  free(dimids);

#ifdef NETCDF3
  CHECK_ERROR( nc_enddef(ncid) );
#endif

  switch ( precision ) {
  case 8:
    CHECK_ERROR( nc_put_var_double(ncid, varid, (double*)val) );
    break;
  case 4:
    CHECK_ERROR( nc_put_var_float(ncid, varid, (float*)val) );
    break;
  default:
    fprintf(stderr, "unsuppoted data precision: %d\n", precision);
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
  int dimids[NC_MAX_DIMS];
  char tname[File_HSHORT+1];
  int tdimid, tvarid;
  nc_type xtype = -1;
  char buf[File_HMID+1];
  int i, j, n, m;
  size_t size;
  double rmiss = RMISS;
  char coord[File_HMID+1];
  int has_assoc;

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

#ifdef NETCDF3
  CHECK_ERROR( nc_redef(ncid) );
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
      // generate name
      if ( nt == 0 )
	strcpy(tname, "time");
      else
	sprintf(tname, "time%d", nt);
      strcpy(tdims[nt]->name, tname);
      // define time dimension and variable
      CHECK_ERROR( nc_def_dim(ncid, tname, 0, &tdimid) );
      tdims[nt]->dimid = tdimid;
      CHECK_ERROR( nc_def_var(ncid, tname, NC_DOUBLE, 1, &tdimid, &tvarid) );
      tdims[nt]->varid = tvarid;
      strcpy(buf, "time");
      CHECK_ERROR( nc_put_att_text(ncid, tvarid, "long_name", strlen(buf), buf) );
      CHECK_ERROR( nc_put_att_text(ncid, tvarid, "units", strlen(files[fid]->time_units), files[fid]->time_units) );
      // define boundary variable
      if ( nc_inq_dimid(ncid, "nv", &(dimids[1])) != NC_NOERR ) // first called
	CHECK_ERROR( nc_def_dim(ncid, "nv", 2, &(dimids[1])) ); // number of vertices
      sprintf(buf, "%s_bnds", tname);
      CHECK_ERROR( nc_put_att_text(ncid, tvarid, "bounds", strlen(buf), buf) );
      dimids[0] = tdimid;
      CHECK_ERROR( nc_def_var(ncid, buf, NC_DOUBLE, 2, dimids, &tvarid) );
      tdims[nt]->bndsid = tvarid;
      CHECK_ERROR( nc_put_att_text(ncid, tvarid, "units", strlen(files[fid]->time_units), files[fid]->time_units) );

      vars[nvar]->t = tdims[nt];
      nt++;
    }
  }

  // get dimension IDs
  // note: C and Fortran order are opposit
  n = ndims;
  if ( tint > 0.0 ) { // add time dimension
    dimids[0] = vars[nvar]->t->dimid;
    ndims++;
  }

  has_assoc = 0;
  for (i=0; i<n; ) {
    if ( nc_inq_dimid(ncid, dims[i], &(dimids[ndims-i-1])) == NC_NOERR ) {
      i += 1;
    } else {
      CHECK_ERROR( nc_inq_varid(ncid, dims[i], &acid) );
      CHECK_ERROR( nc_inq_varndims(ncid, acid, &m) );
      if ( i+m > ndims ) {
	fprintf(stderr, "Error: invalid associated coordinates\n");
	return ERROR_CODE;
      }
      acdimids = (int*) malloc((sizeof(int)*m));
      CHECK_ERROR( nc_inq_vardimid(ncid, acid, acdimids) );
      for (j=0; j<m; j++) dimids[ndims-i-j-1] = acdimids[j];
      has_assoc = 1;
      i += m;
    }
  }
  TYPE2NCTYPE(dtype, xtype);
  CHECK_ERROR( nc_def_var(ncid, varname, xtype, ndims, dimids, &varid) );

  // put variable attribute
  CHECK_ERROR( nc_put_att_text(ncid, varid, "long_name", strlen(desc), desc) );
  CHECK_ERROR( nc_put_att_text(ncid, varid, "units", strlen(units), units) );
  CHECK_ERROR( nc_put_att_double(ncid, varid, _FillValue, xtype, 1, &rmiss) );
  CHECK_ERROR( nc_put_att_double(ncid, varid, "missing_value", xtype, 1, &rmiss) );
  if ( has_assoc ) {
    strcpy(coord, dims[0]);
    for(i=1; i<n; i++) {
      if (strlen(coord)+strlen(dims[i])+1 < File_HMID) {
	strcat(coord, " ");
	strcat(coord, dims[i]);
      }
    }
    if ( ndims > n && strlen(coord)+6 < File_HMID)
      strcat(coord, " time");
    CHECK_ERROR( nc_put_att_text(ncid, varid, "coordinates", strlen(coord), coord) );
  }


  if ( tavg ) {
    sprintf(buf, "%s: mean", tname);
    CHECK_ERROR( nc_put_att_text(ncid, varid, "cell_methods", strlen(buf), buf) );
  }

  // set start and count
  vars[nvar]->start = (size_t*) malloc(sizeof(size_t)*ndims);
  vars[nvar]->count = (size_t*) malloc(sizeof(size_t)*ndims);
  for ( i=0; i<ndims; i++ ) {
    CHECK_ERROR( nc_inq_dimlen(ncid, dimids[i], &size) );
    vars[nvar]->count[i] = size;
    vars[nvar]->start[i] = 0;
  }
  if ( tint > 0.0 ) vars[nvar]->count[0] = 1;

#ifndef NETCDF3
  // set chunk size and deflate level
  if ( files[fid]->deflate_level > 0 ) {
    CHECK_ERROR( nc_def_var_chunking(ncid, varid, NC_CHUNKED, vars[nvar]->count) );
    CHECK_ERROR( nc_def_var_deflate(ncid, varid, 0, 1, files[fid]->deflate_level) );
  }
#endif

#ifdef NETCDF3
  CHECK_ERROR( nc_enddef(ncid) );
#endif

  vars[nvar]->varid = varid;
  *vid = nvar;
  nvar++;

  return SUCCESS_CODE;
}

int32_t file_write_data( int32_t  vid,        // (in)
			 void    *var,        // (in)
			 real64_t t_start,    // (in)
			 real64_t t_end,      // (in)
			 int32_t  precision)  // (in)
{
  int ncid, varid;

  if ( vars[vid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = vars[vid]->ncid;
  varid = vars[vid]->varid;
  if ( vars[vid]->t != NULL ) { // have time dimension
    if ( vars[vid]->t->count < 0 ||  // first time
	 fabs(t_end - vars[vid]->t->t) > EPS ) { // time goes next step
      vars[vid]->t->count += 1;
      vars[vid]->t->t = t_end;
      size_t index[2];
      index[0] = vars[vid]->t->count;
      CHECK_ERROR( nc_put_var1_double(ncid, vars[vid]->t->varid, index, &t_end) );
      index[1] = 0;
      CHECK_ERROR( nc_put_var1_double(ncid, vars[vid]->t->bndsid, index, &t_start) );
      index[1] = 1;
      CHECK_ERROR( nc_put_var1_double(ncid, vars[vid]->t->bndsid, index, &t_end) );
    }
    vars[vid]->start[0] = vars[vid]->t->count;
  }

  switch (precision) {
  case 8:
    CHECK_ERROR( nc_put_vara_double(ncid, varid, vars[vid]->start, vars[vid]->count, (double*)var) );
    break;
  case 4:
    CHECK_ERROR( nc_put_vara_float(ncid, varid, vars[vid]->start, vars[vid]->count, (float*)var) );
    break;
  default:
    fprintf(stderr, "unsuppoted data precision: %d\n", precision);
    return ERROR_CODE;
  }

  CHECK_ERROR( nc_sync(ncid) );

  return SUCCESS_CODE;
}

int32_t file_close( int32_t fid ) // (in)
{
  int ncid;
  int i;

  if ( files[fid] == NULL ) return ALREADY_CLOSED_CODE;
  ncid = files[fid]->ncid;

  free( files[fid] );
  files[fid] = NULL;

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
      free( tdims[i] );
      tdims[i] = NULL;
    }
  }

  CHECK_ERROR( nc_close(ncid) );

  return SUCCESS_CODE;
}
