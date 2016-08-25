/* #define DBGOUT stderr */
#define DBGOUT stdout
/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (with HDF5)                              *
 *                                                                    *
 *  HISTORY                                                           *
 *    1.24      15-11-04  T.Inoue  : [NEW] based on fio.c.             *
 *                                                                    *
 **********************************************************************
 *functions
 * : put/get indicates the communication with the database
 *   read/write indicates the communication with the file
 *<utilities>
 *void           hio_ednchg                   : endian changer
 *void           hio_mk_fname                 : filename generator
 *void           hio_set_str                  : string preprocess
 *<database>
 *int32_t        hio_syscheck                 : check system
 *int32_t        hio_put_commoninfo           : store common informtation
 *static int32_t hio_new_finfo                : add new file structure
 *static int32_t hio_new_datainfo             : add new datainfo structure
 *int32_t        hio_put_pkginfo              : put package information (full)
 *int32_t        hio_get_pkginfo              : get package information (full)
 *int32_t        hio_put_datainfo             : put data information (full)
 *int32_t        hio_get_datainfo             : get data information (full)
 *int32_t        hio_seek_datainfo            : seek data id by varname
 *<file r/w>
 *int32_t        hio_fopen                    : open file IO stream
 *int32_t        hio_fclose                   : close file IO stream
 *int32_t        hio_write_pkginfo            : write package information
 *int32_t        hio_read_pkginfo             : read package information
 *int32_t        hio_write_datainfo           : write data information
 *int32_t        hio_read_datainfo            : read data information
 *int32_t        hio_write_data               : write data array
 *int32_t        hio_write_data_1rgn          : write data array (1 region)
 *int32_t        hio_read_data                : read data array
 *<function suite for fortran program>
 *int32_t        hio_register_file            : register new file
 *int32_t        hio_put_write_pkginfo        : put & write package information (quick put)
 *int32_t        hio_valid_pkginfo            : validate package information with common
 *int32_t        hio_valid_datainfo           : validate data size
 *int32_t        hio_put_write_datainfo_data  : put & write data information and write data
 *int32_t        hio_put_write_datainfo       : put & write data information
 *int32_t        hio_read_allinfo_get_pkginfo : read pkginfo and datainfo and get pkginfo
 *int32_t        hio_copy_datainfo            : allocate and copy datainfo
 *int32_t        hio_dump_finfolist           : dump package summary of all finfo
 *int32_t        hio_dump_finfo               : dump package detail of finfo
 **********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "hio.h"
#include "poh5.h"

/* file ID counter */
int32_t hio_num_of_file = 0;

/* common information for all packages */
hio_commoninfo_t common;

/* package+data+status container */
static hio_fileinfo_t *finfo = NULL;


/* for temporary data (13-04-18 C.Kodama) */ /* \todo check this */
/* int hio_nvar = 0; */
/* char vname[500][HIO_HSHORT]; */

/* declare private routines. Defined in bottom of this file. */
void init_dinfo( hio_datainfo_t dinfo );

/** endian change *****************************************************/
void hio_ednchg( void* avp_pointer,
                 const int32_t ai_size,
                 const int32_t ai_num   )
{
  /* nothing to do for HDF5 */
}

/** filename generator ************************************************/
void hio_mk_fname( char *fname,
                   const char *base,
                   const char *ext,
                   int32_t i,
                   int32_t y )
{
  char _fname[HIO_HLONG];

  switch (y) {
  case 4 :
    sprintf(_fname,"%s.%s%04d.h5",base,ext,i);
    break;
  case 5 :
    sprintf(_fname,"%s.%s%05d.h5",base,ext,i);
    break;
  case 6 :
    sprintf(_fname,"%s.%s%06d.h5",base,ext,i);
    break;
  default :
    break;
  }

  hio_trim_str( fname, _fname, HIO_HLONG-1 );
}

/** string preprocess *************************************************/
void hio_set_str( char *_str,
                  const char *str,
                  const int str_len )
{
  int i;

  strncpy(_str, str, str_len);

  _str[str_len] = '\0'; /* [fix] H.Yashiro 20120621 */
  for( i=str_len-1; i>=0; i-- ) {
    if( _str[i] == ' ' ) {
      _str[i] = '\0';
    } else {
      break;
    }
  }
}

/** string postprocess *************************************************/
void hio_trim_str( char *_str,
                  const char *str,
                  const int str_len )
{
  int i;

  strncpy(_str, str, str_len);

  _str[str_len] = ' ';
  for( i=str_len-1; i>=0; i-- ) {
    if( _str[i] == '\0' ) {
      _str[i] = ' ';
    } else {
      break;
    }
  }
}

/** check system & initialze ******************************************/
int32_t hio_syscheck( void )
{
  int32_t i=1;

  /* Initialize poh5. */
  poh5_setup();

  if ( (sizeof(real32_t)!=4) || (sizeof(real64_t)!=8) ) {
    printf("Data type (real) is inconsistent!\n");
    exit(1);
  }

  /* if ( *(char*)&i ) { */
  /*   system_endiantype = HIO_LITTLE_ENDIAN; */
  /* } else { */
  /*   system_endiantype = HIO_BIG_ENDIAN; */
  /* } */

  /* intitialize */
  common.fmode         = -1;
  common.endiantype    = -1;
  common.grid_topology = -1;
  common.glevel        = -1;
  common.rlevel        = -1;
  common.num_of_rgn    =  0;
  common.rgnid         = NULL;

  return(SUCCESS_CODE);
}

/** put common informtation *******************************************/
int32_t hio_put_commoninfo( int32_t fmode,
                            int32_t endiantype,
                            int32_t grid_topology,
                            int32_t glevel,
                            int32_t rlevel,
                            int32_t num_of_rgn,
                            int32_t rgnid[]        )
{
  int32_t i;

  common.fmode         = fmode;
  common.endiantype    = endiantype;
  common.grid_topology = grid_topology;
  common.glevel        = glevel;
  common.rlevel        = rlevel;
  common.num_of_rgn    = num_of_rgn;
  common.rgnid         = (int32_t *)malloc(num_of_rgn*sizeof(int32_t));
  for( i=0; i<common.num_of_rgn; i++ ) {
    common.rgnid[i] = rgnid[i];
  }


  return(SUCCESS_CODE);
}

/** put common informtation from file *********************************/
int32_t hio_put_commoninfo_fromfile( int32_t fid,
                                     int32_t endiantype )
{
  int32_t i;


  hio_read_pkginfo( fid );

  common.fmode         = finfo[fid].header.fmode;
  common.endiantype    = endiantype;
  common.grid_topology = finfo[fid].header.grid_topology;
  common.glevel        = finfo[fid].header.glevel;
  common.rlevel        = finfo[fid].header.rlevel;
  common.num_of_rgn    = finfo[fid].header.num_of_rgn;
  common.rgnid         = (int32_t *)malloc(common.num_of_rgn*sizeof(int32_t));
  for( i=0; i<common.num_of_rgn; i++ ) {
    common.rgnid[i] = finfo[fid].header.rgnid[i];
  }

  return(SUCCESS_CODE);
}

/** add new file structure ********************************************/
static int32_t hio_new_finfo( void )
{
  int32_t fid;

  /* get file ID */
  fid = hio_num_of_file++;

  /* memory re-allocation (expand by new hio_num_of_file) */
  finfo=(hio_fileinfo_t *)realloc(finfo,sizeof(hio_fileinfo_t)*(hio_num_of_file));

  /* intitialize */
  strcpy(finfo[fid].header.fname,"");
  strcpy(finfo[fid].header.description,"");
  strcpy(finfo[fid].header.note,"");
  finfo[fid].header.fmode         = -1;
  finfo[fid].header.endiantype    = HIO_UNKNOWN_ENDIAN;
  finfo[fid].header.grid_topology = -1;
  finfo[fid].header.glevel        = -1;
  finfo[fid].header.rlevel        = -1;
  finfo[fid].header.num_of_rgn    = 0;
  finfo[fid].header.rgnid         = NULL;
  finfo[fid].header.num_of_var   = 0;

  finfo[fid].dinfo = NULL;

  finfo[fid].status.rwmode    = -1;
  finfo[fid].status.opened    = 0;
  finfo[fid].status.hfid      = -1; /* for negative means error */
  /* finfo[fid].status.eoh.__pos = 0; [add] 20111007 H.Yashiro */

  return(fid);
}

/** add new file structure ********************************************/
static int32_t hio_new_datainfo( int32_t fid )
{
  int32_t did;

  /* get data ID */
  did = finfo[fid].header.num_of_var++;

  /* memory re-allocation (expand by new num_of_var) */
  if ( (finfo[fid].dinfo
       = (hio_datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(hio_datainfo_t)*finfo[fid].header.num_of_var)) == NULL ) {
    printf("Allocation error!\n");
  }

  init_dinfo( finfo[fid].dinfo[did] );
  return(did);
}

/** put package information (full) ************************************/
int32_t hio_put_pkginfo( int32_t fid,
                         hio_headerinfo_t hinfo )
{
  int32_t i;

  hio_set_str( finfo[fid].header.description,hinfo.description,HIO_HMID-1  );
  hio_set_str( finfo[fid].header.note,       hinfo.note,       HIO_HLONG-1 );
  finfo[fid].header.num_of_var   = hinfo.num_of_var;
  finfo[fid].header.fmode         = hinfo.fmode;
  finfo[fid].header.endiantype    = hinfo.endiantype;
  finfo[fid].header.grid_topology = hinfo.grid_topology;
  finfo[fid].header.glevel        = hinfo.glevel;
  finfo[fid].header.rlevel        = hinfo.rlevel;
  finfo[fid].header.num_of_rgn    = hinfo.num_of_rgn;

  finfo[fid].header.rgnid = (int32_t *)realloc(finfo[fid].header.rgnid,
                            sizeof(int32_t)*finfo[fid].header.num_of_rgn);
  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    finfo[fid].header.rgnid[i] = hinfo.rgnid[i];
  }

  return(SUCCESS_CODE);
}

/** get package information (full) ************************************/
hio_headerinfo_t hio_get_pkginfo( int32_t fid )
{
  hio_headerinfo_t hinfo;
  int32_t i;

  hio_trim_str( hinfo.description,finfo[fid].header.description,HIO_HMID-1  );
  hio_trim_str( hinfo.note,       finfo[fid].header.note,HIO_HLONG-1 );
  hinfo.num_of_var   = finfo[fid].header.num_of_var;
  hinfo.fmode         = finfo[fid].header.fmode;
  hinfo.endiantype    = finfo[fid].header.endiantype;
  hinfo.grid_topology = finfo[fid].header.grid_topology;
  hinfo.glevel        = finfo[fid].header.glevel;
  hinfo.rlevel        = finfo[fid].header.rlevel;
  hinfo.num_of_rgn    = finfo[fid].header.num_of_rgn;

  hinfo.rgnid = (int32_t *)malloc(sizeof(int32_t)*hinfo.num_of_rgn);
  for( i=0; i<hinfo.num_of_rgn; i++ ) {
    hinfo.rgnid[i] = finfo[fid].header.rgnid[i];
  }

  return(hinfo);
}

/** put data information (full) ***************************************/
int32_t hio_put_datainfo( int32_t fid,
                          int32_t did,
                          hio_datainfo_t ditem )
{
  hio_set_str( finfo[fid].dinfo[did].varname,    ditem.varname,    HIO_HSHORT-1 );
  hio_set_str( finfo[fid].dinfo[did].description,ditem.description,HIO_HMID-1   );
  hio_set_str( finfo[fid].dinfo[did].unit,       ditem.unit,       HIO_HSHORT-1 );
  hio_set_str( finfo[fid].dinfo[did].layername,  ditem.layername,  HIO_HSHORT-1 );
  hio_set_str( finfo[fid].dinfo[did].note,       ditem.note,       HIO_HLONG-1  );
  finfo[fid].dinfo[did].datasize     = ditem.datasize;
  finfo[fid].dinfo[did].datatype     = ditem.datatype;
  finfo[fid].dinfo[did].num_of_layer = ditem.num_of_layer;

  return(SUCCESS_CODE);
}

/** get data information (full) ***************************************/
hio_datainfo_t hio_get_datainfo( const int32_t fid,
                             const int32_t did  )
{
  hio_datainfo_t ditem;

  hio_trim_str( ditem.varname,    finfo[fid].dinfo[did].varname,    HIO_HSHORT-1 );
  hio_trim_str( ditem.description,finfo[fid].dinfo[did].description,HIO_HMID-1   );
  hio_trim_str( ditem.unit,       finfo[fid].dinfo[did].unit,       HIO_HSHORT-1 );
  hio_trim_str( ditem.layername,  finfo[fid].dinfo[did].layername,  HIO_HSHORT-1 );
  hio_trim_str( ditem.note,       finfo[fid].dinfo[did].note,       HIO_HLONG-1  );
  ditem.datasize     = finfo[fid].dinfo[did].datasize;
  ditem.datatype     = finfo[fid].dinfo[did].datatype;
  ditem.num_of_layer = finfo[fid].dinfo[did].num_of_layer;
  ditem.num_of_step  = finfo[fid].dinfo[did].num_of_step;

  return(ditem);
}

/** get time information(ts,te) ***************************************/
void hio_get_timeinfo( const int32_t fid,
                       const int32_t did,
                       int64_t ts[],
                       int64_t te[])
{
  phid_t vid = poh5_open_variable_by_idx(finfo[fid].status.hfid, did);
  poh5_read_variable_time(vid, ts, te);
}





/** seek data id by varname **********************************/
int32_t hio_seek_datainfo( int32_t fid,
                           char *varname,
                           int32_t step    )
{
  int32_t did;

  char vv[HIO_HSHORT+1];

  hio_set_str( vv, varname, HIO_HSHORT-1);
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_seek_datainfo:fid=%d,varname=%s\n",fid,vv);
  fprintf(DBGOUT,"dbg:hio_seek_datainfo:num_of_var=%d\n",finfo[fid].header.num_of_var);
#endif

  for( did=0; did<finfo[fid].header.num_of_var; did++ ) {
    if ( finfo[fid].dinfo != NULL ) {
      if (    strncmp(finfo[fid].dinfo[did].varname,vv,HIO_HSHORT) == 0 ) {
        return(did);
      }
    }
  }

  return(ERROR_CODE);
}


/** return num_of_var in file *****************************************/
int32_t hio_get_num_of_var( const int32_t fid )
{

  return finfo[fid].header.num_of_var;

}




/** open file IO stream ***********************************************/
int32_t hio_fopen( const int32_t fid, const int32_t mode )
{
  phid_t hfid;
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_fopen:fid=%d,mode=%d\n",fid,mode);
#endif
  if (finfo[fid].status.opened) {
    fprintf(stdout," FIle ( %s ) has been already opened!\n",finfo[fid].header.fname);
    fprintf(stdout," open process will be skipped!\n");
    return(SUCCESS_CODE);
  }

  if ( mode==HIO_FWRITE ) {
    hfid = poh5_open_file( finfo[fid].header.fname, POH5_FWRITE);
 } else if( mode==HIO_FREAD )   { /* avoid overwrite action */
    hfid = poh5_open_file(finfo[fid].header.fname, POH5_FREAD);
  } else if( mode==HIO_FAPPEND ) { /* overwrite mode */
    hfid = poh5_open_file(finfo[fid].header.fname, POH5_FAPPEND );
  }
  if ( hfid < 0 ) {
    fprintf(stdout,"Can not open file : %s!\n",finfo[fid].header.fname);
    exit(1);
  }else{
    finfo[fid].status.hfid = hfid;
  }
  finfo[fid].status.rwmode = mode;
  finfo[fid].status.opened = 1;

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_fopen:hfid=%d\n",hfid);
#endif
  return(SUCCESS_CODE);
}

/** close file IO stream **********************************************/
int32_t hio_fclose( int32_t fid )
{
  int ret = poh5_close_file( finfo[fid].status.hfid );

  finfo[fid].status.opened = 0;

  return(SUCCESS_CODE);
}

/** write package information *****************************************/
int32_t hio_write_pkginfo( int32_t fid )
{

  if(!finfo[fid].status.opened){
    fprintf(stdout,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  poh5_write_global_attr(
                         finfo[fid].status.hfid,
                         finfo[fid].header.glevel,
                         finfo[fid].header.rlevel,
                         finfo[fid].header.grid_topology,
                         finfo[fid].header.fmode,
                         finfo[fid].header.num_of_rgn,
                         finfo[fid].header.rgnid,
                         finfo[fid].header.description,
                         finfo[fid].header.note,
                         finfo[fid].header.num_of_var);


  return(SUCCESS_CODE);
}

/** read package information ******************************************/
int32_t hio_read_pkginfo( int32_t fid )
{
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_read_pkginfo:fid=%d\n",fid);
#endif

  if (!finfo[fid].status.opened) {
    fprintf(stdout,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  poh5_read_global_attr(
                         finfo[fid].status.hfid,
                         &finfo[fid].header.glevel,
                         &finfo[fid].header.rlevel,
                         &finfo[fid].header.grid_topology,
                         &finfo[fid].header.fmode,
                         &finfo[fid].header.num_of_rgn,
                         &finfo[fid].header.rgnid,
                         finfo[fid].header.description,
                         finfo[fid].header.note,
                         &finfo[fid].header.num_of_var);

#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_read_pkginfo:num_of_var=%d\n",finfo[fid].header.num_of_var);
#endif
  return(SUCCESS_CODE);
}

/** write data information ********************************************/
int32_t hio_write_datainfo( int32_t fid,
                            int32_t did  )
{
  if(!finfo[fid].status.opened){
    fprintf(stdout,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  int gall1d = pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2;

  poh5_create_variable( finfo[fid].status.hfid,
                        finfo[fid].dinfo[did].varname,
                        finfo[fid].dinfo[did].description,
                        finfo[fid].dinfo[did].note,
                        finfo[fid].dinfo[did].unit,
                        finfo[fid].dinfo[did].layername,
                        finfo[fid].dinfo[did].num_of_layer,
                        gall1d,
                        finfo[fid].dinfo[did].datatype,
                        finfo[fid].header.num_of_rgn);

  /* \todo Set gall1d in hio_datainfo_t instead of datasize ? */


  return(SUCCESS_CODE);
}

/** read data information *********************************************/
int32_t hio_read_datainfo( int32_t fid )
{
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_read_datainfo:fid=%d\n",fid);
  fprintf(DBGOUT,"dbg:hio_read_datainfo:num_of_var=%d\n",finfo[fid].header.num_of_var);
#endif
  if (!finfo[fid].status.opened) {
    fprintf(stdout,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }


  phid_t vid;
  int ret;
  int nv;

  for( nv=0; nv<finfo[fid].header.num_of_var; nv++ ) {
    vid = poh5_open_variable_by_idx(finfo[fid].status.hfid, nv);
    ret = poh5_read_variable_attr(vid,
                                  finfo[fid].dinfo[nv].varname,
                                  finfo[fid].dinfo[nv].description,
                                  finfo[fid].dinfo[nv].note,
                                  finfo[fid].dinfo[nv].unit,
                                  finfo[fid].dinfo[nv].layername,
                                  &finfo[fid].dinfo[nv].num_of_layer,
                                  &finfo[fid].dinfo[nv].num_of_step,
                                  &finfo[fid].dinfo[nv].datatype);
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_read_datainfo:num_of_var=%d\n",finfo[fid].header.num_of_var);
    fprintf(DBGOUT,"dbg:hio_read_datainfo:nv=%d,varname=%s\n",nv,finfo[fid].dinfo[nv].varname);
#endif
  }

  return(SUCCESS_CODE);
}

/** write data array **************************************************/
int32_t hio_write_data( int32_t fid,
                        int32_t did,
                        int32_t step,
                        int64_t ts,
                        int64_t te,
                        void *data   )
{
  phid_t v_gid = poh5_open_variable( finfo[fid].status.hfid, finfo[fid].dinfo[did].varname);
  poh5_write_variable_data(
                           v_gid,
                           step,
                           ts,
                           te,
                           finfo[fid].dinfo[did].datatype,
                           data);
  poh5_close_variable( v_gid );


  return(SUCCESS_CODE);
}

/** write data array (1 region) ***************************************/
int32_t hio_write_data_1rgn( int32_t fid,
                             int32_t did,
                             int32_t step,
                             int32_t ll,/* rgn number, 0-based */
                             int64_t ts,
                             int64_t te,
                             void *data   )
{
  phid_t v_gid = poh5_open_variable( finfo[fid].status.hfid, finfo[fid].dinfo[did].varname);
  poh5_write_variable_data_1rgn(
                                v_gid,
                                step,
                                ll,
                                ts,
                                te,
                           finfo[fid].dinfo[did].datatype,
                           data);
  poh5_close_variable( v_gid );


  return(SUCCESS_CODE);
}

/** read data array (full size) ***************************************/
int32_t hio_read_data( int32_t fid,
                       int32_t did,
                       const int32_t step,
                       int64_t *ts,
                       int64_t *te,
                       void *data   )
{
  phid_t v_gid
    = poh5_open_variable( finfo[fid].status.hfid, finfo[fid].dinfo[did].varname);
  poh5_read_variable_data(
                   v_gid, /**< [in] group id of variable */
                   step,  /**< [in] step counter */
                   ts,         /**< [out] start time of this step */
                   te,         /**< [out] end time of this step */
                   finfo[fid].dinfo[did].datatype,     /**< [in] HIO_{REAL4,REAL8,INTEGER4,INTEGER8} */
                   data);         /**< [out] data */


  return(SUCCESS_CODE);
}

/** register new file *************************************************/
int32_t hio_register_file( char *fname )
{
  int32_t ierr;
  int32_t fid = -1;
  int i;

  for ( i=0;i<hio_num_of_file; i++){
    if ( strcmp(finfo[i].header.fname,fname) == 0 ){
      fid = i;
    }
  }

  if ( fid < 0 ) {
    /* request new file space */
    fid = hio_new_finfo();
    hio_set_str( finfo[fid].header.fname,fname,HIO_HLONG-1  );
  }

  return(fid);
}

/** put & write package information (quick put) ***********************/
int32_t hio_put_write_pkginfo( int32_t fid,
                               char *description,
                               char *note         )
{
  int32_t i;

  hio_set_str( finfo[fid].header.description,description,HIO_HMID-1  );
  hio_set_str( finfo[fid].header.note,       note,       HIO_HLONG-1 );
  /* use common info */
  finfo[fid].header.fmode         = common.fmode;
  finfo[fid].header.endiantype    = common.endiantype;
  finfo[fid].header.grid_topology = common.grid_topology;
  finfo[fid].header.glevel        = common.glevel;
  finfo[fid].header.rlevel        = common.rlevel;
  finfo[fid].header.num_of_rgn    = common.num_of_rgn;

  finfo[fid].header.rgnid = (int32_t *)realloc(finfo[fid].header.rgnid,
                            sizeof(int32_t)*finfo[fid].header.num_of_rgn);
  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    finfo[fid].header.rgnid[i] = common.rgnid[i];
  }

  hio_write_pkginfo( fid );

  return(SUCCESS_CODE);
}

/** validate package information with common **************************/
int32_t hio_valid_pkginfo( int32_t fid )
{
  int32_t i;

  if(finfo[fid].header.grid_topology!=common.grid_topology) {
    fprintf(stdout,"Warning: grid_topology is not match, %d, %d\n",
                   finfo[fid].header.grid_topology,common.grid_topology);
  }

  if(finfo[fid].header.glevel!=common.glevel) {
    fprintf(stdout,"Warning: glevel is not match, %d, %d\n",
                   finfo[fid].header.glevel,common.glevel              );
  }

  if(finfo[fid].header.rlevel!=common.rlevel) {
    fprintf(stdout,"Warning: rlevel is not match, %d, %d\n",
                   finfo[fid].header.rlevel,common.rlevel              );
  }

  if(finfo[fid].header.num_of_rgn!=common.num_of_rgn) {
    fprintf(stdout,"Warning: num_of_rgn is not match, %d, %d\n",
                   finfo[fid].header.num_of_rgn,common.num_of_rgn      );
  }

  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    if(finfo[fid].header.rgnid[i]!=common.rgnid[i]) {
    fprintf(stdout,"Warning: rgnid[%d] is not match, %d, %d\n",
                   i,finfo[fid].header.rgnid[i],common.rgnid[i]        );
    }
  }
  return(SUCCESS_CODE);
}

/** validate package information with common **************************/
int32_t hio_valid_pkginfo_validrgn( int32_t fid,
                                    int32_t rgnid[] )
{
  int32_t i;

  if(finfo[fid].header.grid_topology!=common.grid_topology) {
    fprintf(stdout,"Warning: grid_topology is not match, %d, %d\n",
                   finfo[fid].header.grid_topology,common.grid_topology);
  }

  if(finfo[fid].header.glevel!=common.glevel) {
    fprintf(stdout,"Warning: glevel is not match, %d, %d\n",
                   finfo[fid].header.glevel,common.glevel              );
  }

  if(finfo[fid].header.rlevel!=common.rlevel) {
    fprintf(stdout,"Warning: rlevel is not match, %d, %d\n",
                   finfo[fid].header.rlevel,common.rlevel              );
  }

  if(finfo[fid].header.num_of_rgn!=common.num_of_rgn) {
    fprintf(stdout,"Warning: num_of_rgn is not match, %d, %d\n",
                   finfo[fid].header.num_of_rgn,common.num_of_rgn      );
  }

  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    if(finfo[fid].header.rgnid[i]!=rgnid[i]) {
    fprintf(stdout,"Warning: rgnid[%d] is not match, %d, %d\n",
                   i,finfo[fid].header.rgnid[i],rgnid[i]               );
    }
  }
  return(SUCCESS_CODE);
}

/** validate data size ************************************************/
int32_t hio_valid_datainfo( int32_t fid )
{
  /* \todo re-implement or is this necessary ?? */

  /* fprintf(stdout,"dbg:hio_valid_datainfo:start\n"); */
 /*  for( int did=0; did<finfo[fid].header.num_of_var; did++ ) { */
 /*    int ns; */
 /*    ns = finfo[fid].dinfo[did].num_of_step; */
 /*    fprintf(stdout,"dbg:varname=%s\n",finfo[fid].dinfo[did].varname); */
 /*    fprintf(stdout,"dbg:ns=%d\n",ns); */
 /*    for( int n=0; n<ns; n++){ */
 /*      fprintf(stdout,"dbg:ts[%d]=%ld\n",ns,finfo[fid].dinfo[did].ts[n]); */
 /*      fprintf(stdout,"dbg:te[%d]=%ld\n",ns,finfo[fid].dinfo[did].te[n]); */
 /*    } */
 /* } */

  /* fprintf(stdout,"dbg:hio_valid_datainfo:end\n"); */
  return(SUCCESS_CODE);
}

/** put & write data information and write data ***********************/
int32_t hio_put_write_datainfo_data( int32_t fid,
                                     hio_datainfo_t ditem,
                                     int32_t step,
                                     int64_t ts,
                                     int64_t te,
                                     void *data        )
{
  int32_t did;

  did = hio_seek_datainfo( fid, ditem.varname,step ); /* step is neglected, only for compatibility. */
  if ( did < 0 ) {
    did = hio_new_datainfo( fid );
  }

  hio_put_datainfo( fid, did, ditem );

  hio_write_pkginfo( fid );
  hio_write_datainfo( fid, did );
  hio_write_data( fid, did, step, ts, te, data );

  return(did);
}

/** put & write data information **************************************/
int32_t hio_put_write_datainfo( int32_t fid,
                                hio_datainfo_t ditem )
{
  int32_t did;
  int32_t dummy;

  did = hio_seek_datainfo( fid, ditem.varname, dummy );
  if ( did < 0 ) {
    did = hio_new_datainfo( fid );
  }
  hio_put_datainfo( fid, did, ditem );
  hio_write_datainfo( fid, did );

  return(did);
}

/** read pkginfo and datainfo *****************************************/
int32_t hio_read_allinfo( int32_t fid )
{
  int32_t i;
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_read_allinfo:fid=%d\n",fid);
#endif
  hio_read_pkginfo( fid );
  hio_valid_pkginfo( fid );
#ifdef DEBUG
  fprintf(DBGOUT,"dbg:hio_read_allinfo:num_of_var=%d\n",finfo[fid].header.num_of_var);
#endif

  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (hio_datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(hio_datainfo_t)*finfo[fid].header.num_of_var)) == NULL ) {
    printf("Allocation error!\n");
  }
  int n;
  for ( n=0; n<finfo[fid].header.num_of_var; n++){
    init_dinfo( finfo[fid].dinfo[n] );
  }

  hio_read_datainfo( fid );
  hio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

/** read pkginfo and datainfo, with validating rgnid ******************/
int32_t hio_read_allinfo_validrgn( int32_t fid,
                                   int32_t rgnid[] )
{
  int32_t i;

  hio_read_pkginfo( fid );
  hio_valid_pkginfo_validrgn( fid, rgnid );

  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (hio_datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(hio_datainfo_t)*finfo[fid].header.num_of_var)) == NULL ) {
    printf("Allocation error!\n");
  }

  hio_read_datainfo( fid );
  hio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

/** allocate and copy datainfo ****************************************/
/* [add] C.Kodama 13-04-18 */
int32_t hio_copy_datainfo( int32_t fid, int32_t fid_org )
{
  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (hio_datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(hio_datainfo_t)*finfo[fid].header.num_of_var)) == NULL ) {
    printf("Allocation error!\n");
  }

  /*hio_read_datainfo( fid );*/
  memcpy( finfo[fid].dinfo, finfo[fid_org].dinfo, sizeof(hio_datainfo_t)*finfo[fid].header.num_of_var );
  hio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

/** dump package summary of all finfo *********************************/
int32_t hio_dump_finfolist( void )
{
  int32_t i;
  int32_t fid;
  char* str_mode[3]     = { "XXXX","SPRT","INTG" };
  char* str_rw[4]       = { "X","R","W","A" };
  char* str_opened[3]   = { "XXX","NO","YES" };
  char* str_endian[4]   = { "XXXXXX","UKNOWN","LITTLE","BIG" };
  char* str_topology[4] = { "XXXX","ICO","LCP","MLCP" };

  printf( "========== common information ==========\n" );
  printf( " MODE ENDIAN GRID GL RL LALL\n" );
  printf( " %4s", str_mode[common.fmode+1] );
  printf( " %6s", str_endian[common.endiantype+1] );
  printf( " %4s", str_topology[common.grid_topology+1] );
  printf( " %2d", common.glevel );
  printf( " %2d", common.rlevel );
  printf( " %4d", common.num_of_rgn );
  printf( "\n" );

  printf( "================== file information ===================\n" );
  printf( " MODE R/W OPEN ENDIAN GRID GL RL LALL ITEMS DESCRIPTION\n" );
  for( fid=0; fid<hio_num_of_file; fid++ ) {
    printf( " %4s", str_mode[finfo[fid].header.fmode+1] );
    printf( " %3s", str_rw[finfo[fid].status.rwmode+1] );
    printf( " %4s", str_opened[finfo[fid].status.opened+1] );
    printf( " %6s", str_endian[finfo[fid].header.endiantype+1] );
    printf( " %4s", str_topology[finfo[fid].header.grid_topology+1] );
    printf( " %2d", finfo[fid].header.glevel );
    printf( " %2d", finfo[fid].header.rlevel );
    printf( " %4d", finfo[fid].header.num_of_rgn );
    printf( " %5d", finfo[fid].header.num_of_var );
    printf( " %s" , finfo[fid].header.description );
    printf( "\n" );
    /*
    for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
      printf( " %5d", finfo[fid].header.rgnid[i] );
    }
    printf( "\n" );
    */
  }
  return(SUCCESS_CODE);
}

/** dump package summary of all finfo *********************************/
int32_t hio_dump_finfo( int32_t fid,
                        int32_t endiantype,
                        int32_t dumptype    )
{

  int32_t i,j,ij,k,l;
  int32_t did;
  int32_t ijall;
  int64_t ijklall,pos;
  void *buf;
  char *c;

  char* str_mode[3]     = { "XXXX","SPRIT","COMPLETE" };
  char* str_topology[4] = { "XXXX","ICOSAHEDRON","LCP","MLCP" };
  char* str_dtype[5]    = { "XXXX","REAL4","REAL8","INTEGER4","INTEGER8" };


  hio_fopen( fid, HIO_FREAD );
  hio_read_pkginfo ( fid );
  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (hio_datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(hio_datainfo_t)*finfo[fid].header.num_of_var)) == NULL ) {
    printf("Allocation error!\n");
  }
  hio_read_datainfo( fid );

  ijall = (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
        * (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2);

  /* package info */
  printf( "============ DATA PACKAGE DEFINITION =============\n" );
  printf( "--- PACKAGE DESCRIPTION : %s\n", finfo[fid].header.description );
  printf( "--- COMPLETE/SPRIT      : %s\n", str_mode[finfo[fid].header.fmode+1] );
  printf( "--- GRID TOPOLOGY       : %s\n", str_topology[finfo[fid].header.grid_topology+1] );
  printf( "--- GLEVEL              : %d\n", finfo[fid].header.glevel );
  printf( "--- RLEVEL              : %d\n", finfo[fid].header.rlevel );
  printf( "--- NUMBER OF GRIDS         \n" );
  printf( "---     FOR EACH REGION : %d\n", ijall );
  printf( "--- NUMBER OF REGION    : %d\n", finfo[fid].header.num_of_rgn );
  printf( "--- INCLUDING REGION ID\n");
  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    if((i%10)==0) { printf("--- "); }
    printf(" %06d",finfo[fid].header.rgnid[i]);
    if((i%10)==9) { printf("\n");   }
  }
  if((i%10)!=0) { printf("\n");   }
  printf( "--- NUMBER OF VARS      : %d\n", finfo[fid].header.num_of_var );
  printf( "--- NOTE                : %s\n", finfo[fid].header.note );
  printf( "\n" );

  printf("============ VARS ATTRIBUTES ======================\n");
  /* fsetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /\* [add] 20111007 H.Yashiro *\/ */
  /* pos = dinfosize;                                         /\* [mod] 20111007 H.Yashiro *\/ */
  for( did=0; did<finfo[fid].header.num_of_var; did++ ) {
    if ( dumptype == HIO_DUMP_ALL || dumptype == HIO_DUMP_ALL_MORE ) {
      /* data info */
      printf( "########## Item ID = %4d ###############\n",did);
      printf( "--- VAR NAME            : %s\n",  finfo[fid].dinfo[did].varname );
      printf( "--- DATA DESCRIPTION    : %s\n",  finfo[fid].dinfo[did].description );
      printf( "--- UNIT                : %s\n",  finfo[fid].dinfo[did].unit );
      printf( "--- VERT. LAYER NAME    : %s\n",  finfo[fid].dinfo[did].layername );
      printf( "--- NUMBER OF LAYERS    : %d\n",  finfo[fid].dinfo[did].num_of_layer );
      printf( "--- NUMBER OF STEPS     : %d\n",  finfo[fid].dinfo[did].num_of_step );
      printf( "--- TIMES                   \n" );
      int n;
      /*d disabled tentatively */
      /* for(n=0;n<finfo[fid].dinfo[did].num_of_step;n++){ */
      /*   if((n%10)==0) { printf("--- "); } */
      /*   printf(" (%ld,%ld) ",finfo[fid].dinfo[did].ts[n],finfo[fid].dinfo[did].te[n]); */
      /*   if((n%10)==9) { printf("\n");   } */
      /* } */
      if((n%10)!=0) { printf("\n");   }
      printf( "--- DATA TYPE           : %s\n",  str_dtype[finfo[fid].dinfo[did].datatype+1] );
      printf( "--- NOTE                : %s\n",  finfo[fid].dinfo[did].note );
    } else {
    /* data info */
      printf( "+%4d:%16s|%16s[%3d]|%dsteps|%s\n",
            did,
            finfo[fid].dinfo[did].varname,
            finfo[fid].dinfo[did].layername,
            finfo[fid].dinfo[did].num_of_layer,
            finfo[fid].dinfo[did].num_of_step,
            str_dtype[finfo[fid].dinfo[did].datatype+1] );
    }

    /* data array */
    if ( dumptype == HIO_DUMP_ALL || dumptype == HIO_DUMP_ALL_MORE) {
      int64_t  *vdatai8 = NULL;
      int32_t  *vdatai4 = NULL;
      real64_t  *vdatar8 = NULL;
      real32_t  *vdatar4 = NULL;
      int32_t step;
      int64_t ts,te;
      int status;

      char *r4form, *r8form;
      char *i4form, *i8form;
      switch( dumptype ){
      case HIO_DUMP_ALL :
        i4form = "%d\n";
        i8form = "%lld\n";
        r4form = "%f\n";
        r8form = "%lf\n";
        break;
      case HIO_DUMP_ALL_MORE :
        i4form = "%d\n";
        i8form = "%lld\n";
        r4form = "%.60e\n";
        r8form = "%.60le\n";
        break;
      }
      /* void dump_variable_data(int32_t,int64_t, int64_t, int, void *); */

      size_t datasize = ijall
        * finfo[fid].dinfo[did].num_of_layer
        * finfo[fid].header.num_of_rgn;


      for ( step=1;step<=finfo[fid].dinfo[did].num_of_step;step++){
        switch(finfo[fid].dinfo[did].datatype){
        case HIO_INTEGER4 :
          printf(" %15s: %s\n","datatype","INTEGER 4byte\n");
          vdatai4 = (int32_t *)malloc(datasize * sizeof(int32_t) );
          status = hio_read_data(fid,did,step,&ts,&te, vdatai4);
          break;
        case HIO_INTEGER8 :
          printf(" %15s: %s\n","datatype","INTEGER 8byte\n");
          vdatai8 = (int64_t *)malloc(datasize * sizeof(int64_t) );
          status = hio_read_data(fid,did,step,&ts,&te, vdatai8);
          break;
        case HIO_REAL4:
          printf(" %15s: %s\n","datatype","REAL 4byte\n");
          vdatar4 = (real32_t *)malloc(datasize * sizeof(real32_t) );
          status = hio_read_data(fid,did,step,&ts,&te, vdatar4);
          break;
        case HIO_REAL8:
          printf(" %15s: %s\n","datatype","REAL 8byte\n");
          vdatar8 = (real64_t *)malloc(datasize * sizeof(real64_t) );
          status = hio_read_data(fid,did,step,&ts,&te, vdatar8);
          break;
        };
        int gall1d = pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2;
        int ijkl = 0;

        for( l=0; l<finfo[fid].header.num_of_rgn; l++ ) {
          for( k=0; k<finfo[fid].dinfo[did].num_of_layer; k++ ) {
            for( j=0; j<gall1d; j++ ) {
              for( i=0; i<gall1d; i++ ) {
                int ij = j*gall1d+i;
                printf("+++ [%8s, %4d, %6d, %3d, %6d] ",
                       finfo[fid].dinfo[did].varname,step,ij,k,l);
                switch(finfo[fid].dinfo[did].datatype){
                case HIO_INTEGER4 :
                  printf(i4form, vdatai4[ijkl]);
                  break;
                case HIO_INTEGER8 :
                  printf(i8form, vdatai8[ijkl]);
                  break;
                case HIO_REAL4 :
                  printf(r4form, vdatar4[ijkl]);
                  break;
                case HIO_REAL8 :
                  printf(r8form, vdatar8[ijkl]);
                  break;

                } /* switch */
                ijkl++;
              } /* i */
            } /* j */
          } /* k */
        } /* l */
      } /* step */
    }
  }
  printf("============ END ITEM ATTRIBUTES ==================\n");
  hio_fclose( fid );

  return(SUCCESS_CODE);
}


void init_dinfo( hio_datainfo_t dinfo ){

  strcpy(dinfo.varname,"");
  strcpy(dinfo.description, "");
  strcpy(dinfo.unit,"");
  strcpy(dinfo.layername,"");
  strcpy(dinfo.note,"");
  dinfo.datasize= -1;
  dinfo.datatype= -1;
  dinfo.num_of_layer = -1;
  dinfo.num_of_step = -1;
  /* dinfo.ts=NULL; */
  /* dinfo.te=NULL; */

}
