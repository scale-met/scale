/*
 * \todo set hlng_tid etc. in poh5_setup().
 */

#include <stdlib.h>

#include "poh5.h"
#include "hdf5.h"

typedef struct{
  int64_t s;
  int64_t e;
} tp_t;

/* #define DBGOUT stderr */
#define DBGOUT stdout

/* strings type */
static hid_t hlng_tid;
static hid_t hmid_tid;
static hid_t hsht_tid;


/*========================================================================*/
/*
 * Declaration of internal routines
 *
 */
int write_gattr(const hid_t file_id, const char *name, const hid_t f_tid, const hid_t sid,
                const hid_t m_tid, const void *vals);
int update_num_of_var( hid_t file_id );
int extend_dataset( hid_t v_gid, char *dname);
int write_dataset( const hid_t   v_gid, const int32_t step,  const char   *dname, 
                   const hid_t   type,  const void   *data);
int write_dataset_1rgn( const hid_t   v_gid, const int32_t step,  const int32_t ll,   
                        const hid_t   type, const void   *data);
int read_dataset( const hid_t   v_gid, const int32_t step,  const char   *dname,
                  const hid_t   type,  void         *data);

void check_h5( int );

/*========================================================================*/
/*
 * Implementation of public routines
 *
 */

/*========================================================================*/
int poh5_setup(){

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_setup:\n");
#endif

  /* strings type */
  hlng_tid = H5Tcopy(H5T_C_S1);  H5Tset_size(hlng_tid, POH5_HLNG);
  hmid_tid = H5Tcopy(H5T_C_S1);  H5Tset_size(hmid_tid, POH5_HMID);
  hsht_tid = H5Tcopy(H5T_C_S1);  H5Tset_size(hsht_tid, POH5_HSHT);

  /* Turn off error handling permanently */
  H5Eset_auto (H5E_DEFAULT, NULL, NULL);

  return 0;
}


/*========================================================================*/
hid_t poh5_open_file( const char* fname, int mode ){
  /* mode is one of [POH5_FREAD, POH5_FWRITE, POH5_FAPPEND] */

  hid_t fid;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_open_file:mode=%d,fname=%s\n",mode,fname);
#endif

  switch (mode){
  case POH5_FWRITE:   /* existing content will be vanished on open.*/
    fid = H5Fcreate(fname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);    break;
  case POH5_FREAD:    /* avoid overwrite action */
    fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);    break;
  case POH5_FAPPEND:  /* overwrite mode */
    fid = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);    break;
  default:
    fprintf(stderr,"Err:poh5_open_file:invalid mode:%d\n",mode);
    exit (1);
  }

  check_h5( fid > 0 );
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_open_file:fid=%d\n",fid);
#endif
 
  return fid;
}


/*========================================================================*/
int poh5_close_file( const hid_t fid ){

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_close_file:fid=%d\n",fid);
#endif
     
  int res = H5Fclose(fid);

  return 0;
}
  

/*========================================================================*/
int poh5_write_global_attr(
                      const hid_t    fid,       /*< [in] poh5 file_id */
                      const int32_t  glevel,        /*< [in] glevel */
                      const int32_t  rlevel,        /*< [in] rlevel */
                      const int32_t  grid_topology, /*< [in] grid_topology */
                      const int      is_complete,   /*< [in] 1 if complete file */
                      const int32_t  num_of_rgn,    /*< [in] num of region this file contains. */
                      const int32_t *rgnid,         /*< [in] array of region id's */
                      const char    *description,   /*< [in] description of this file */
                      const char    *note,          /*< [in] longer note of this file */ 
                      const int32_t  num_of_var)   /*< [in] number of data in this file. */
{
  hsize_t dim1[1];

  int poh5_version=POH5_VERSION;

  /* strings type */
  /* hid_t hlng_tid = H5Tcopy(H5T_C_S1);  H5Tset_size(hlng_tid, POH5_HLNG); */
  /* hid_t hmid_tid = H5Tcopy(H5T_C_S1);  H5Tset_size(hmid_tid, POH5_HMID); */
  /* hid_t hsht_tid = H5Tcopy(H5T_C_S1);  H5Tset_size(hsht_tid, POH5_HSHT); */

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_write_grobal_attr:fid=%d\n",fid);
#endif

  hid_t scalar_sid = H5Screate(H5S_SCALAR);

  check_h5(0== write_gattr(fid, "poh5_version", H5T_STD_I32BE,scalar_sid,H5T_NATIVE_INT, &poh5_version));
  check_h5(0== write_gattr(fid, "glevel",       H5T_STD_I32BE,scalar_sid,H5T_NATIVE_INT, &glevel));
  check_h5(0== write_gattr(fid, "rlevel",       H5T_STD_I32BE,scalar_sid,H5T_NATIVE_INT, &rlevel));
  check_h5(0== write_gattr(fid, "grid_topology",H5T_STD_I32BE,scalar_sid,H5T_NATIVE_INT, &grid_topology));
  check_h5(0== write_gattr(fid, "is_complete",  H5T_STD_I32BE,scalar_sid,H5T_NATIVE_INT, &is_complete));

  check_h5(0== write_gattr(fid, "description", hmid_tid,scalar_sid,hmid_tid, description));
  check_h5(0== write_gattr(fid, "note", hlng_tid,scalar_sid,hlng_tid, note));


  check_h5(0== write_gattr(fid, "num_of_var", H5T_STD_I32BE, scalar_sid, H5T_NATIVE_INT, &num_of_var));
  check_h5(0== write_gattr(fid, "num_of_rgn", H5T_STD_I32BE, scalar_sid, H5T_NATIVE_INT, &num_of_rgn));

  dim1[0] = num_of_rgn;
  check_h5(0== write_gattr(fid, "rgnid", H5T_STD_I32BE, H5Screate_simple(1,dim1,NULL), H5T_NATIVE_INT, rgnid));

  return 0;

};




/*========================================================================*/
int poh5_read_global_attr(
                      const hid_t file_id,    /**< [in] poh5 file_id */
                      int32_t *glevel,        /**< [out] glevel */
                      int32_t *rlevel,        /**< [out] rlevel */
                      int32_t *grid_topology, /**< [out] grid_topology */
                      int     *is_complete,   /**< [out] 1 if complete file */
                      int32_t *num_of_rgn,    /**< [out] num of region this file contains. */
                      int32_t *rgnid[],       /**< [out] array of region id's */
                      char    description[64],   /**< [out] description of this file */
                      char    note[256] ,         /**< [out] longer note of this file */ 
                      int32_t *num_of_var)   /**< [out] number of data in this file. */

{
  herr_t res;
  int poh5_version;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_read_grobal_attr:fid=%d\n",file_id);
#endif

  {
    hid_t aid;
    aid = H5Aopen(file_id, "poh5_version", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, &poh5_version);
    res = H5Aclose(aid);

    /* Is this OK ??? */
    if ( poh5_version != POH5_VERSION ){
      printf("poh5 version mismatch: file:%d, library:%d\n",poh5_version,POH5_VERSION);
      exit(1);
    }

    aid = H5Aopen(file_id, "glevel", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, glevel);
    res = H5Aclose(aid);

    aid = H5Aopen(file_id, "rlevel", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, rlevel);
    res = H5Aclose(aid);

    aid = H5Aopen(file_id, "grid_topology", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, grid_topology);
    res = H5Aclose(aid);
  
    aid = H5Aopen(file_id, "is_complete", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, is_complete);
    res = H5Aclose(aid);

    aid = H5Aopen(file_id, "description",H5P_DEFAULT);
    res = H5Aread(aid, hmid_tid, description);
    res = H5Aclose(aid);

    aid = H5Aopen(file_id, "note", H5P_DEFAULT);
    res = H5Aread(aid, hlng_tid, note);
    res = H5Aclose(aid);

    aid = H5Aopen(file_id, "num_of_rgn", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, num_of_rgn);
    res = H5Aclose(aid);

    aid = H5Aopen(file_id, "num_of_var", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, num_of_var);
    res = H5Aclose(aid);

    size_t nn = *num_of_rgn;
    (*rgnid) = (int32_t *)malloc(nn * sizeof(int32_t));
    aid = H5Aopen(file_id, "rgnid", H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, *rgnid);
    res = H5Aclose(aid);
  }

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_read_global_attr:poh5_version=%d\n", poh5_version);
  fprintf(DBGOUT,"dbg:poh5_read_global_attr:gl=%d,rl=%d\n",*glevel,*rlevel);
#endif
  
  return 0;

};



/*========================================================================*/
hid_t poh5_open_variable(
                    const hid_t file_id, /**<[in] poh5 file id .*/
                    const char *vname)   /**<[in] variable name */
{

  hid_t v_gid;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_open_variable:fid=%d,vname=%s\n",file_id,vname);
#endif

  /* open toplevel group id for variable: "/Var" to open subgroup.*/
  hid_t top_gid = H5Gopen(file_id, "/Var", H5P_DEFAULT);
  check_h5 ( top_gid > 0 );

  /* check if specified variable exist */
  htri_t retval = H5Lexists( top_gid, vname, H5P_DEFAULT );

  if      ( retval>0  ) {    /* exist */
    v_gid = H5Gopen(top_gid, vname, H5P_DEFAULT);
  } else {    /* not exist or something wrong*/
    v_gid = -1;
  }
  
  herr_t res = H5Gclose(top_gid);

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_open_variable:v_gid=%d\n",v_gid);
#endif
  return v_gid;
}



/*========================================================================*/
herr_t poh5_close_variable(
                          const hid_t v_gid) /**<[in] group id .*/
{
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_close_variable:v_gid=%d\n",v_gid);
#endif

  herr_t res = H5Gclose( v_gid );
  return res;
}
  


/*========================================================================*/
hid_t poh5_open_variable_by_idx(
                    const hid_t file_id, /**<[in] poh5 file id.*/
                    const int idx)      /**<[in] variable index */
{

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_open_variable_by_idx:fid=%d,idx=%d\n",file_id,idx);
#endif

  hid_t v_gid = H5Oopen_by_idx(
                                 file_id, /* hid_t loc_id, */
                                 "/Var",      /* const char *group_name, */
                                 H5_INDEX_NAME, /* H5_index_t index_type, */
                                 H5_ITER_INC,   /* H5_iter_order_t order,*/
                                 idx,             /* hsize_t n, */
                                 H5P_DEFAULT );  /* hid_t lapl_id */
  check_h5( v_gid > 0 );

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_open_variable_by_idx:v_gid=%d\n",v_gid);
#endif

  return v_gid;
}

  

/*========================================================================*/
/*
 */

hid_t poh5_create_variable(
                    const hid_t file_id, /**< [in] poh5 file id */
                    const char *name,      /**< [in] variable name */
                    const char *dscr,      /**< [in] variable description */
                    const char *note,         /**< [in] variable long note */
                    const char *unit,      /**< [in] variable unit */
                    const char *lname,    /**< [in] name of vertical layer */
                    const int32_t nlayer, /**< [in] number of layer */
                    const int gall1d,     /**< [in] gall1d */
                    const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                    const int num_of_rgn)   /**< [in] number of regions */
{
  hid_t v_gid;
  herr_t res;

  hid_t scalar_sid = H5Screate(H5S_SCALAR);

  /* group id of this variable */
  {
    hid_t top_gid;
    htri_t retval = H5Lexists( file_id, "/Var", H5P_DEFAULT);
    if      ( retval>0  ) {    /* exist */
      top_gid = H5Gopen(file_id, "/Var", H5P_DEFAULT);
    } else {    /* not exist or something wrong*/
      top_gid = H5Gcreate(file_id, "/Var", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
  
    retval = H5Lexists( top_gid, name, H5P_DEFAULT);
    if ( retval > 0 ) { /* already exist */
      v_gid = H5Gopen(top_gid, name, H5P_DEFAULT);
      check_h5( v_gid > 0 );
      res = H5Gclose(top_gid);
      return v_gid;
    }else{ /* not created yet */
      v_gid = H5Gcreate(top_gid, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      check_h5( v_gid > 0 );
      res = H5Gclose(top_gid);
    }
  }

  /* variable attributes */
  {
    hid_t aid;
    aid = H5Acreate(v_gid, "varname", hsht_tid, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, hsht_tid, name);
    res = H5Aclose(aid);

    aid = H5Acreate(v_gid, "unit", hsht_tid, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, hsht_tid, unit);
    res = H5Aclose(aid);

    aid = H5Acreate(v_gid, "description", hmid_tid, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, hmid_tid, dscr);
    res = H5Aclose(aid);

    aid = H5Acreate(v_gid, "note", hlng_tid, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, hlng_tid, note);
    res = H5Aclose(aid);

    aid = H5Acreate(v_gid, "layername", hsht_tid, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, hsht_tid, lname);
    res = H5Aclose(aid);

    aid = H5Acreate(v_gid, "num_of_layer", H5T_STD_I32BE, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, H5T_NATIVE_INT,&nlayer);
    res = H5Aclose(aid);

    int32_t num_of_steps = 0;
    aid = H5Acreate(v_gid, "num_of_steps", H5T_STD_I32BE, scalar_sid, H5P_DEFAULT, H5P_DEFAULT);
    res = H5Awrite(aid, H5T_NATIVE_INT,&num_of_steps);
    res = H5Aclose(aid);

    /* aid = H5Acreate(v_gid, "datasize", H5T_STD_I64LE,scalar_sid, H5P_DEFAULT, H5P_DEFAULT); */
    /* res = H5Awrite(aid, H5T_STD_I64LE,&datasize); */
  }

  /*
   * 5D dataset, [T,L,K,J,I]
   */
  {
    hsize_t dims[5], max_dims[5], chunk_dims[5];
    dims[0] = 0 ;
    dims[1] = num_of_rgn;
    dims[2] = nlayer;
    dims[3] = gall1d;
    dims[4] = gall1d;

    max_dims[0] = H5S_UNLIMITED;
    max_dims[1] = num_of_rgn;
    max_dims[2] = nlayer;
    max_dims[3] = gall1d;
    max_dims[4] = gall1d;

    /* chunk size should be as large as possible for limited dimension on non-compression.*/
    chunk_dims[0] = 1;
    chunk_dims[1] = num_of_rgn;
    chunk_dims[2] = nlayer;
    chunk_dims[3] = gall1d;
    chunk_dims[4] = gall1d;

#ifdef DEBUG
    {
      fprintf(DBGOUT,"dbg:poh5_create_variable:chunk_dims:\n");
      for(int n=0;n<5;n++){  printf("%d,",(int)chunk_dims[n]); } ; printf("\n");
    }
#endif

    /* Modify dataset creation properties */
    hid_t c_pid = H5Pcreate (H5P_DATASET_CREATE);
    /* First: Chunking */
    res = H5Pset_chunk(c_pid, 5, chunk_dims);  
    /* Second: Bit shuffle. */
#ifdef POH5_DO_SHUF    
    res = H5Pset_shuffle(c_pid);
#endif
    /* Third: Gzip. */
#if POH5_DO_GZIP > 0
    res = H5Pset_deflate (c_pid, POH5_DO_GZIP); 
#endif

    /* datatype for variable data */
    hid_t d_tid;
    switch (dtype){
    case POH5_REAL4:    d_tid = H5T_IEEE_F32BE; break;
    case POH5_REAL8:    d_tid = H5T_IEEE_F64BE; break;
    case POH5_INTEGER4: d_tid = H5T_STD_I32BE; break;
    case POH5_INTEGER8: d_tid = H5T_STD_I64BE; break;
    default: break;
    };

    hid_t d_did = H5Dcreate(v_gid, "data", d_tid, H5Screate_simple(5,dims,max_dims),
                            H5P_DEFAULT, c_pid, H5P_DEFAULT);

    res = H5Pclose(c_pid); 
    res = H5Dclose(d_did);
  }

  /*
   * Dataset for time_period.
   * Note that these must be unlimited size and need to enable chunking.
   */

  {
    hsize_t dims[1], max_dims[1], chunk_dims[1];
    hid_t c_tid = H5Tcreate( H5T_COMPOUND, sizeof(tp_t));
    H5Tinsert(c_tid, "start", HOFFSET(tp_t, s),H5T_STD_I64BE);
    H5Tinsert(c_tid, "end"  , HOFFSET(tp_t, e),H5T_STD_I64BE);
  
    dims[0] = 0;
    max_dims[0] = H5S_UNLIMITED;
    chunk_dims[0] = 1; /* tentative */

#ifdef DEBUG
    {
      fprintf(DBGOUT,"dbg:poh5_create_variable:chunk_dims:\n");
      for(int n=0;n<1;n++){  printf("%d,",(int)chunk_dims[n]); } ; printf("\n");
    }
#endif

    hid_t c_pid = H5Pcreate (H5P_DATASET_CREATE);
    res = H5Pset_chunk(c_pid, 1, chunk_dims);

    hid_t t_did = H5Dcreate(v_gid, "times", c_tid, H5Screate_simple(1,dims,max_dims),
                            H5P_DEFAULT, c_pid, H5P_DEFAULT);
  
    res = H5Pclose(c_pid); 
    res = H5Dclose(t_did); 
  }

  /* update num_of_var */
  res = update_num_of_var( file_id );

  return v_gid;

};



/*========================================================================*/
void poh5_write_variable_data(
                              const hid_t v_gid, /**<[in] poh5 file id .*/
                              const int32_t step,      /**< step count, 1-based */
                              const int64_t ts, /**< start time this data represents */
                              const int64_t te,   /**< end time this data represents */
                              const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                              const void *var_data)   /**< variable data at current step */
{
  int32_t nsteps;
  herr_t res;

  /* nstep is overwritten below, do NOT close aid here */
  hid_t aid = H5Aopen(v_gid, "num_of_steps", H5P_DEFAULT);
  res = H5Aread(aid, H5T_NATIVE_INT, &nsteps);

  {
    /* \todo check step != nstep+1 */
    if (step > nsteps) {
      res = extend_dataset(v_gid, "data");
      res = extend_dataset(v_gid, "times");
      ++nsteps;
    }

    hid_t d_tid; /* datatype ID for data */
    switch (dtype){
    case POH5_REAL4:    d_tid = H5T_NATIVE_FLOAT;  break;
    case POH5_REAL8:    d_tid = H5T_NATIVE_DOUBLE; break;
    case POH5_INTEGER4: d_tid = H5T_NATIVE_INT;    break;
    case POH5_INTEGER8: d_tid = H5T_NATIVE_LLONG;  break;
    default: break;
    };

    res = write_dataset(v_gid, step, "data", d_tid, var_data);
  }
  {
    tp_t tp;
    tp.s = ts;
    tp.e = te;

    /* datatype ID for tp_t */
    hid_t c_tid = H5Tcreate( H5T_COMPOUND, sizeof(tp_t));
    /* Since c_tid is for memory here, use _NATIVE_ */
    H5Tinsert(c_tid, "start", HOFFSET(tp_t, s), H5T_NATIVE_LLONG);
    H5Tinsert(c_tid, "end"  , HOFFSET(tp_t, e), H5T_NATIVE_LLONG);

    res = write_dataset(v_gid, step, "times", c_tid, &tp);
  }  

  /* update numr_of_steps */
  res = H5Awrite(aid, H5T_NATIVE_INT, &nsteps);
  res = H5Aclose(aid);

}



/*========================================================================*/
void poh5_write_variable_data_1rgn(
                              const hid_t v_gid, /**<[in] poh5 file id .*/
                              const int32_t step,      /**< step count, 1-based */
                              const int32_t ll, /**< rgn number, 0-based */
                              const int64_t ts, /**< start time this data represents */
                              const int64_t te,   /**< end time this data represents */
                              const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                              const void *var_data)   /**< variable data at current step */
{
  herr_t res;

  int32_t nsteps;

  /* nstep is overwritten below, do NOT close aid here */
  hid_t aid = H5Aopen(v_gid, "num_of_steps", H5P_DEFAULT);
  res = H5Aread(aid, H5T_NATIVE_INT, &nsteps);

  {
    /* \todo check step != nstep+1 */
    if (step > nsteps) {
      res = extend_dataset(v_gid, "data");
      res = extend_dataset(v_gid, "times");
      ++nsteps;
    }

    hid_t d_tid; /* datatype ID for data */
    switch (dtype){
    case POH5_REAL4:    d_tid = H5T_NATIVE_FLOAT;  break;
    case POH5_REAL8:    d_tid = H5T_NATIVE_DOUBLE; break;
    case POH5_INTEGER4: d_tid = H5T_NATIVE_INT;    break;
    case POH5_INTEGER8: d_tid = H5T_NATIVE_LLONG;  break;
    default: break;
    };

    res = write_dataset_1rgn(v_gid, step, ll, d_tid, var_data);
  }
  {
    tp_t tp;
    tp.s = ts;
    tp.e = te;

    /* datatype ID for tp_t */
    hid_t c_tid = H5Tcreate( H5T_COMPOUND, sizeof(tp_t));
    /* Since c_tid is for memory here, use _NATIVE_ */
    H5Tinsert(c_tid, "start", HOFFSET(tp_t, s), H5T_NATIVE_LLONG);
    H5Tinsert(c_tid, "end"  , HOFFSET(tp_t, e), H5T_NATIVE_LLONG);

    res = write_dataset(v_gid, step, "times", c_tid, &tp);
  }

  /* update num_of_steps */
  res = H5Awrite(aid, H5T_NATIVE_INT, &nsteps);
  res = H5Aclose(aid);
}

  



/*========================================================================*/
int poh5_read_variable_attr(
                   const hid_t v_gid, /**< [in] group id of variable */
                   char    *name,       /**< [out] variable name  */
                   char    *dscr,       /**< [out] variable description */
                   char    *note,       /**< [out] variable long note */
                   char    *unit,       /**< [out] variable unit */
                   char    *lname,      /**< [out] name of vertical layer */
                   int32_t *nlayer,     /**< [out] number of layers */
                   int32_t *nsteps,     /**< [out] number of steps */
                   int     *dtype)      /**< [out] POH5_{REAL4,REAL8,INTEGER4,INTEGER8} */
{                          
  herr_t res;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_read_variable_attr:v_gid=%d\n",v_gid);
#endif

  /* read variable attributes */
  {
    hid_t aid;
    aid = H5Aopen(v_gid, "varname",H5P_DEFAULT);
    res = H5Aread(aid, hmid_tid, name);
    res = H5Aclose(aid);

    aid = H5Aopen(v_gid, "description",H5P_DEFAULT);
    res = H5Aread(aid, hmid_tid, dscr);
    res = H5Aclose(aid);
    
    aid = H5Aopen(v_gid, "note",H5P_DEFAULT);
    res = H5Aread(aid, hlng_tid, note);
    res = H5Aclose(aid);
    
    aid = H5Aopen(v_gid, "unit",H5P_DEFAULT);
    res = H5Aread(aid, hsht_tid, unit);
    res = H5Aclose(aid);
    
    aid = H5Aopen(v_gid, "layername",H5P_DEFAULT);
    res = H5Aread(aid, hsht_tid, lname);
    res = H5Aclose(aid);
    
    aid = H5Aopen(v_gid, "num_of_layer",H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, nlayer);
    res = H5Aclose(aid);

    aid = H5Aopen(v_gid, "num_of_steps",H5P_DEFAULT);
    res = H5Aread(aid, H5T_NATIVE_INT, nsteps);
    res = H5Aclose(aid);
  }
  
  /* get datatype(int or float) and bitsize */
  {
    hid_t did = H5Dopen(v_gid, "data",H5P_DEFAULT);
    hid_t tid = H5Dget_type(did);
    H5T_class_t cls_t = H5Tget_class(tid);
    size_t bsize = H5Tget_size(tid);

    switch (cls_t){
    case H5T_INTEGER:
      switch( bsize ){
      case 4:        *dtype = POH5_INTEGER4;break;
      case 8:        *dtype = POH5_INTEGER8;break;
      default:        break;
      }      break;
    case H5T_FLOAT:
      switch( bsize ){
      case 4:        *dtype = POH5_REAL4;break;
      case 8:        *dtype = POH5_REAL8;break;
      default:        break;
      }      break;
    default:
      exit(1);
    }
    res = H5Tclose(tid);
    res = H5Dclose(did);
  }

  return 0;


}

/*========================================================================*/
/* Caution:
 * ts and te must be allocated BEFORE call 
 */
int poh5_read_variable_time( 
                   const hid_t v_gid, /**< [in] group id of variable */
                   int64_t *ts,       /**< [out] array of start time of each steps */
                   int64_t *te)       /**< [out] array of end time of each steps */
{
  int32_t n,ns;
  herr_t res;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:poh5_read_variable_time:v_gid=%d\n",v_gid);
#endif

  /* read num of steps from variable attributes */
  {
    hid_t aid;
    aid = H5Aopen(v_gid, "num_of_steps",H5P_DEFAULT);
    check_h5( aid > 0 );
    res = H5Aread(aid, H5T_NATIVE_INT, &ns);
    check_h5( res == 0 );
    res = H5Aclose(aid);
  }

  /* read tstart/tend */
  {

    tp_t *tp = (tp_t*)malloc( ns * sizeof(tp_t));
    /* datatype ID for tp_t */
    hid_t c_tid = H5Tcreate( H5T_COMPOUND, sizeof(tp_t));
    /* Since c_tid is for memory here, use _NATIVE_ */
    H5Tinsert(c_tid, "start", HOFFSET(tp_t, s), H5T_NATIVE_LLONG);
    H5Tinsert(c_tid, "end"  , HOFFSET(tp_t, e), H5T_NATIVE_LLONG);

    hid_t  did = H5Dopen(v_gid,"times",H5P_DEFAULT);
    check_h5( did > 0 );
    res = H5Dread(did, c_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, tp);
    check_h5( res == 0 );
    res = H5Dclose(did);
    res = H5Tclose(c_tid);

    for ( n=0; n<ns; n++ ){
      ts[n] = tp[n].s;
      te[n] = tp[n].e;
    }
    free(tp);
  }

#ifdef DEBUG
  {
    fprintf(DBGOUT,"dbg:poh5_read_variable_time:\n");
    for ( int n=0; n<ns; n++ ){
      fprintf(DBGOUT,"ts[%d]=%ld,te[%d]=%ld\n",n,ts[n],n,te[n]);
    }
  }
#endif

  return 0;
}
/*========================================================================*/
int poh5_read_variable_data(
                   const hid_t v_gid, /**< [in] group id of variable */
                   const int32_t step,  /**< [in] step counter, 1-based */
                   int64_t *ts,         /**< [out] start time of this step */
                   int64_t *te,         /**< [out] end time of this step */
                   const int dtype,     /**< [in] POH5_{REAL4,REAL8,INTEGER4,INTEGER8} */
                   void *vdata)         /**< [out] data */

{
  herr_t   res;
  tp_t tp;

#ifdef DEBUG
     fprintf(DBGOUT,"dbg:poh5_read_variable_data:\n");
#endif

  hid_t d_tid;
  switch (dtype){
  case POH5_REAL4:    d_tid = H5T_NATIVE_FLOAT;  break;
  case POH5_REAL8:    d_tid = H5T_NATIVE_DOUBLE; break;
  case POH5_INTEGER4: d_tid = H5T_NATIVE_INT;    break;
  case POH5_INTEGER8: d_tid = H5T_NATIVE_LLONG;  break;
  default: break;
  };

  res = read_dataset( v_gid, step, "data", d_tid, vdata );

  /* datatype ID for tp_t */
  hid_t c_tid = H5Tcreate(H5T_COMPOUND, sizeof(tp_t));
  /* Since c_tid is for memory here, use _NATIVE_ */
  H5Tinsert(c_tid, "start", HOFFSET(tp_t, s),H5T_NATIVE_LLONG);
  H5Tinsert(c_tid, "end"  , HOFFSET(tp_t, e),H5T_NATIVE_LLONG);

  res = read_dataset(v_gid, step, "times", c_tid, &tp);

  *ts = tp.s;
  *te = tp.e;

  return 0;

}

/*========================================================================*/
/*
 *
 * Internal routines.
 *
 */





/*========================================================================*/
/*
 * Update global attribute num_of_var by counting actual number of
 * variable group in the file.
 *
 */
int update_num_of_var( hid_t file_id )
{

  herr_t res;

  int nvar;
  /*
   * Open "/Var" group, and count subgroups
   */
  {
    H5G_info_t ginfo;
    hid_t gid = H5Gopen(file_id, "/Var", H5P_DEFAULT);
    check_h5 ( gid > 0 );
    res = H5Gget_info(gid, &ginfo);
    nvar = ginfo.nlinks; /* discard ginfo hereafter */
    res = H5Gclose(gid);
  }

  /*
   * Overwrite top-level attribute "num_of_var".
   */
  {
    hid_t aid = H5Aopen(file_id, "num_of_var", H5P_DEFAULT);
    res = H5Awrite(aid, H5T_NATIVE_INT, &nvar);
    check_h5( res == 0 );
    res = H5Aclose(aid);
  }
  return 0;

}

/*========================================================================*/
/*
 * extend existing datasets of given variable in T(step) dimension.
 *
 * `dname` should be one of "data", "times", else cause
 * error for non-existing in the group.
 *
 */
int extend_dataset(
                      hid_t v_gid, /**<[in] variable group id */
                      char *dname)   /**<[in] "data" or  "times" */
{
  herr_t res;

  /* get dataset id */
  hid_t did = H5Dopen(v_gid,dname,H5P_DEFAULT);
  check_h5( did > 0 );

  /* get current dataspace */
  hid_t sid = H5Dget_space(did);    /* dataspace handle */
  check_h5( sid > 0 );

  /* get current dimension */
  hsize_t rank = H5Sget_simple_extent_ndims(sid);
  hsize_t *dim = (hsize_t *)malloc(rank * sizeof(hsize_t));
  res  = H5Sget_simple_extent_dims(sid, dim, NULL);
  check_h5( res == rank );

  /* extend dataset */
  dim[0] += 1;
  res = H5Dset_extent(did, dim);
  check_h5 ( res == 0 );

  free(dim);
  res = H5Sclose(sid);
  res = H5Dclose(did);

  return 0;
}


/*========================================================================*/
/*
 * write dataset of given variable into the hyperslab of
 * (step,L,K,J,I) of dataset.
 *
 * `dname` should be either "data" or "times", else HDF5 raises
 * error for non-existing in the group.
 *
 * Note that step is 1-based.
 *
 * \todo Check if `step` is out of range somewhere.
 *
 */
int write_dataset(
                     const hid_t   v_gid, /**<[in] variable group id */
                     const int32_t step,    /**<[in] step count, 1-based */
                     const char   *dname,   /**<[in] "data" or "times" */
                     const hid_t   type,    /**<[in] datatype of given data */
                     const void   *data)    /**<[in] vardata itself. */
{
  herr_t res;
  int n;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:write_dataset:step,dname=%d %s\n",step,dname);
#endif
    
  /* get current dataset */
  hid_t did = H5Dopen(v_gid,dname,H5P_DEFAULT);
  check_h5( did > 0 );

  /* Get current dataspace ID */
  hid_t sid = H5Dget_space(did);
  check_h5 ( sid > 0 );

  /* Get dimensions */
  hsize_t rank = H5Sget_simple_extent_ndims(sid);
  hsize_t *dim = (hsize_t *)malloc(rank * sizeof(hsize_t));
  res  = H5Sget_simple_extent_dims(sid, dim, NULL);

  /* dataspace for hyperslab */
  hsize_t *cnt = (hsize_t *)malloc(rank * sizeof(hsize_t));
  cnt[0]  = 1;
  for ( n=1; n<rank; n++ ){    cnt[n]  = dim[n];  }

  hsize_t *ofs = (hsize_t *)malloc(rank * sizeof(hsize_t));
  ofs[0] = step-1;
  for ( n=1; n<rank; n++ ){    ofs[n] = 0;  }

  res = H5Sselect_hyperslab(sid, H5S_SELECT_SET, ofs, NULL, cnt, NULL);
  hid_t mem_sid = H5Screate_simple(rank, cnt, NULL);
  check_h5 ( mem_sid > 0 );
  
  /* write data */
  res = H5Dwrite(did, type, mem_sid, sid, H5P_DEFAULT, data);
  check_h5 ( res == 0 );

  free(dim);free(ofs);free(cnt);

  res = H5Sclose(sid);
  res = H5Sclose(mem_sid);
  res = H5Dclose(did);

  return 0;
}


/*========================================================================*/
/*
 * write dataset of given variable, 1rgn only, into the hyperslab of
 * (step,ll,K,J,I) of dataset.
 *
 * Note that step is 1-based.
 *
 * See also write_dataset().
 */
int write_dataset_1rgn(
                     const hid_t   v_gid, /**<[in] variable group id */
                     const int32_t step,    /**<[in] step count, 1-based */
                     const int32_t ll,      /**<[in] rgn number L, 0-based. */
                     const hid_t   type,    /**<[in] datatype of given data */
                     const void   *data)    /**<[in] vardata itself. */
{
  herr_t res;
  int n;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:write_dataset_1rgn:step, ll=%d %d\n",step,ll);
#endif
    
  /* get current dataset */
  hid_t did = H5Dopen(v_gid,"data",H5P_DEFAULT);
  check_h5( did > 0 );

  /* Get current dataspace ID */
  hid_t sid = H5Dget_space(did);
  check_h5 ( sid > 0 );

  /* Get dimensions */
  hsize_t rank = H5Sget_simple_extent_ndims(sid);
  hsize_t *dim = (hsize_t *)malloc(rank * sizeof(hsize_t));
  res  = H5Sget_simple_extent_dims(sid, dim, NULL);
  check_h5 ( res != 0 );

  /* dataspace for hyperslab */
  hsize_t *cnt = (hsize_t *)malloc(rank * sizeof(hsize_t));
  cnt[0]  = 1;
  cnt[1]  = 1;
  for ( n=2; n<rank; n++ ){    cnt[n]  = dim[n];  }

  hsize_t *ofs = (hsize_t *)malloc(rank * sizeof(hsize_t));
  ofs[0] = step-1;
  ofs[1] = ll;
  for ( n=2; n<rank; n++ ){    ofs[n] = 0;  }

  res = H5Sselect_hyperslab(sid, H5S_SELECT_SET, ofs, NULL, cnt, NULL);
  hid_t mem_sid = H5Screate_simple(rank, cnt, NULL);
  check_h5 ( mem_sid > 0 );
  
  /* write data */
  res = H5Dwrite(did, type, mem_sid, sid, H5P_DEFAULT, data);
  check_h5 ( res == 0 );
  
  free(dim);free(ofs);free(cnt);

  res = H5Sclose(sid);
  res = H5Sclose(mem_sid);
  res = H5Dclose(did);

  return 0;
}

/*========================================================================*/
/*
 * read dataset of given variable, hyperslab of (step-1,:,:,:,:) of dataset.
 *
 * \todo check unless `step` is out of range somewhere.
 *
 */
int read_dataset(
                 const hid_t   v_gid, /**<[in] variable group id */
                 const int32_t step,    /**<[in] step count, 1-based */
                 const char   *dname,   /**<[in] "data" or "times" */
                 const hid_t   type,    /**<[in] datatype of given data */
                 void         *data)    /**<[out] vardata itself. */
{
  herr_t res;
  int n;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:read_dataset:step=%d\n",step);
#endif

  /* get current dataset */
  hid_t did = H5Dopen(v_gid,dname,H5P_DEFAULT);
  check_h5( did > 0 );

  /* Get current dataspace ID */
  hid_t sid = H5Dget_space(did);
  check_h5 ( sid > 0 );

  /* Get dimensions */
  hsize_t rank = H5Sget_simple_extent_ndims(sid);
  hsize_t *dim = (hsize_t *)malloc(rank * sizeof(hsize_t));
  res  = H5Sget_simple_extent_dims(sid, dim, NULL);
  check_h5 ( res == rank );
  
  /* dataspace for hyperslab */
  hsize_t *cnt = (hsize_t *)malloc(rank * sizeof(hsize_t));
  cnt[0]  = 1;
  for ( n=1; n<rank; n++ ){    cnt[n] = dim[n];  }

  hsize_t *ofs = (hsize_t *)malloc(rank * sizeof(hsize_t));
  ofs[0] = step-1;
  for ( n=1; n<rank; n++ ){    ofs[n] = 0;  }

  res = H5Sselect_hyperslab(sid, H5S_SELECT_SET, ofs, NULL, cnt, NULL);
  hid_t mem_sid = H5Screate_simple(rank, cnt, NULL);
  check_h5 ( mem_sid > 0 );
  
  /* read data */
  res = H5Dread(did, type, mem_sid, sid, H5P_DEFAULT, data);
  check_h5 ( res == 0 );

  free(dim);free(ofs);free(cnt);

  res = H5Sclose(sid);
  res = H5Sclose(mem_sid);
  res = H5Dclose(did);

  return 0;
  
}


/*========================================================================*/
/*
 * Write one global attribute to given file_id.
 * If specified attribute does NOT exist yet, create it, else overwrite it.
 */


int write_gattr(const hid_t fid, const char *name, const hid_t f_tid, const hid_t sid, const hid_t m_tid, const void *vals)
{
  hid_t aid;
  herr_t res;


#ifdef DEBUG
  fprintf(DBGOUT,"dbg:write_gattr_int32:name=%s\n",name);
  fprintf(DBGOUT,"dbg:write_gattr_int32:exist=%d\n",H5Aexists( fid, name ));
#endif
  
  if ( H5Aexists( fid, name ) >0 ){
    aid = H5Aopen(fid, name, H5P_DEFAULT);
  }else{
    aid = H5Acreate(fid, name, f_tid, sid, H5P_DEFAULT, H5P_DEFAULT);
  }
  check_h5( aid > 0 );
  res = H5Awrite(aid, m_tid, vals);
  check_h5 ( res >= 0 );
  res = H5Aclose(aid);
 
  return 0;
}
  



/*========================================================================*/
/*
 * Check if given stat, call H5Eprint() and exit if (!stat)
 *
 * Use this instead of assert()
 */
void check_h5(int stat){

  if ( ! stat ) {
    H5Eprint(H5E_DEFAULT,stderr);
    exit(1);
  }
}
