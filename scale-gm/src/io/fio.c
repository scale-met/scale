/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (CORE)                                   *
 *                                                                    *
 *  HISTORY                                                           *
 *    0.80      11-07-27  H.Tomita  : [NEW]                           *
 *    0.90      11-08-19  H.Yashiro : Incorporate into NICAM          *
 *    1.00      11-08-25  H.Yashiro : Complete format specification   *
 *    1.10      11-09-07  H.Yashiro : sepalate READ/APPEND mode       *
 *                                    remove \0 in return string      *
 *    1.20      11-10-07  H.Yashiro : sepalate MPI/nonMPI fpos        *
 *                                    thanks to kameyama-san@riken    *
 *    1.21      11-10-28  H.Yashiro : [fix] bad memory allocation in  *
 *                                    fio_write_data_1rgn             *
 *    1.22      12-06-21  H.Yashiro : [fix] bad char tail treatment   *
 *                                    thanks to ohno-san@riken        *
 *    1.23      13-04-18  C.Kodama  : [add] fio_read_datainfo_tmpdata,*
 *                                      fio_register_vname_tmpdata,   *
 *                                      fio_read_data_tmpdata,        *
 *                                      fio_read_allinfo_tmpdata,     *
 *                                      fio_copy_datainfo             *
 *                                                                    *
 **********************************************************************
 *functions
 * : put/get indicates the communication with the database
 *   read/write indicates the communication with the file
 *<utilities>
 *void           fio_ednchg                   : endian changer
 *void           fio_mk_fname                 : filename generator
 *void           fio_set_str                  : string preprocess
 *<database>
 *int32_t        fio_syscheck                 : check system
 *int32_t        fio_put_commoninfo           : store common informtation
 *static int32_t fio_new_finfo                : add new file structure
 *static int32_t fio_new_datainfo             : add new datainfo structure
 *int32_t        fio_put_pkginfo              : put package information (full)
 *int32_t        fio_get_pkginfo              : get package information (full)
 *int32_t        fio_put_datainfo             : put data information (full)
 *int32_t        fio_get_datainfo             : get data information (full)
 *int32_t        fio_seek_datainfo            : seek data id by varname and step
 *<file r/w>
 *int32_t        fio_fopen                    : open file IO stream
 *int32_t        fio_fclose                   : close file IO stream
 *int32_t        fio_write_pkginfo            : write package information
 *int32_t        fio_read_pkginfo             : read package information
 *int32_t        fio_write_datainfo           : write data information
 *int32_t        fio_read_datainfo            : read data information
 *int32_t        fio_read_datainfo_tmpdata    : read data information and store data as tmpdata
 *int32_t        fio_write_data               : write data array
 *int32_t        fio_write_data_1rgn          : write data array (1 region)
 *int32_t        fio_read_data                : read data array
 *int32_t        fio_read_data_tmpdata        : read data array from tmpdata
 *<function suite for fortran program>
 *int32_t        fio_register_file            : register new file
 *int32_t        fio_put_write_pkginfo        : put & write package information (quick put)
 *int32_t        fio_valid_pkginfo            : validate package information with common
 *int32_t        fio_valid_datainfo           : validate data size
 *int32_t        fio_put_write_datainfo_data  : put & write data information and write data
 *int32_t        fio_put_write_datainfo       : put & write data information
 *int32_t        fio_read_allinfo_tmpdata     : read pkginfo and datainfo and store data as tmpdata
 *int32_t        fio_read_allinfo_get_pkginfo : read pkginfo and datainfo and get pkginfo
 *int32_t        fio_copy_datainfo            : allocate and copy datainfo
 *int32_t        fio_dump_finfolist           : dump package summary of all finfo
 *int32_t        fio_dump_finfo               : dump package detail of finfo
 **********************************************************************/
#include "fio.h"

/* file ID counter */
int32_t num_of_file = 0;

/* common information for all packages */
commoninfo_t common;

/* package+data+status container */
fileinfo_t *finfo = NULL;


/* system information */
int32_t system_endiantype    = FIO_UNKNOWN_ENDIAN;
int32_t system_ednchg        = 0;

/* list */
int32_t precision[4] = { 4,8,4,8 };
int32_t dinfosize = sizeof(char)*FIO_HSHORT*3
                  + sizeof(char)*FIO_HMID
                  + sizeof(char)*FIO_HLONG
                  + sizeof(int64_t)*3
                  + sizeof(int32_t)*3;

/* for temporary data (13-04-18 C.Kodama) */
int nvar = 0;
char vname[500][FIO_HSHORT];
struct type_tmpdata{ void **tmpdata; };
struct type_tmpdata *tdata = NULL;

/** endian change *****************************************************/
void fio_ednchg( void* avp_pointer,
                 const int32_t ai_size,
                 const int32_t ai_num   )
{
  int ai_csize, ai_cnum;
  char ac_buf[16];
  char *acp_tmp;
  char *acp_local;

  memset(ac_buf,'\0',sizeof(ac_buf));

  acp_tmp   = avp_pointer;
  acp_local = avp_pointer;
  for( ai_cnum=0; ai_cnum<ai_num; ai_cnum++ ) {
    memcpy(ac_buf, acp_local, ai_size);
    for( ai_csize=0; ai_csize<ai_size; ai_csize++ ) {
      *acp_local = ac_buf[ai_size-ai_csize-1];
      acp_local++;
    }
    acp_tmp += ai_size;
  }
}

/** filename generator ************************************************/
void fio_mk_fname( char *fname,
                   const char *base,
                   const char *ext,
                   int32_t i,
                   int32_t y )
{
  char _fname[FIO_HLONG];

  switch (y) {
  case 4 :
    sprintf(_fname,"%s.%s%04d",base,ext,i);
    break;
  case 5 :
    sprintf(_fname,"%s.%s%05d",base,ext,i);
    break;
  case 6 :
    sprintf(_fname,"%s.%s%06d",base,ext,i);
    break;
  default :
    break;
  }

  fio_trim_str( fname, _fname, FIO_HLONG-1 );
}

/** string preprocess *************************************************/
void fio_set_str( char *_str,
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
void fio_trim_str( char *_str,
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
int32_t fio_syscheck( void )
{
  int32_t i=1;

  if ( (sizeof(real32_t)!=4) || (sizeof(real64_t)!=8) ) {
    printf("Data type (real) is inconsistent!\n");
    exit(1);
  }

  if ( *(char*)&i ) {
    system_endiantype = FIO_LITTLE_ENDIAN;
  } else {
    system_endiantype = FIO_BIG_ENDIAN;
  }

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
int32_t fio_put_commoninfo( int32_t fmode,
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

  /* exchange endian? */
  if ( common.endiantype!=system_endiantype ) {
    system_ednchg = 1;
  }

  return(SUCCESS_CODE);
}

/** put common informtation from file *********************************/
int32_t fio_put_commoninfo_fromfile( int32_t fid,
                                     int32_t endiantype )
{
  int32_t i;

  /* exchange endian? */
  if ( endiantype!=system_endiantype ) {
    system_ednchg = 1;
  }

  fio_read_pkginfo( fid );

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
static int32_t fio_new_finfo( void )
{
  int32_t fid;

  /* get file ID */
  fid = num_of_file++;

  /* memory re-allocation (expand by new num_of_file) */
  finfo=(fileinfo_t *)realloc(finfo,sizeof(fileinfo_t)*(num_of_file));

  /* intitialize */
  strcpy(finfo[fid].header.fname,"");
  strcpy(finfo[fid].header.description,"");
  strcpy(finfo[fid].header.note,"");
  finfo[fid].header.fmode         = -1;
  finfo[fid].header.endiantype    = FIO_UNKNOWN_ENDIAN;
  finfo[fid].header.grid_topology = -1;
  finfo[fid].header.glevel        = -1;
  finfo[fid].header.rlevel        = -1;
  finfo[fid].header.num_of_rgn    = 0;
  finfo[fid].header.rgnid         = NULL;
  finfo[fid].header.num_of_data   = 0;

  finfo[fid].dinfo = NULL;

  finfo[fid].status.rwmode    = -1;
  finfo[fid].status.opened    = 0;
  finfo[fid].status.fp        = NULL;
  /* finfo[fid].status.eoh.__pos = 0; [add] 20111007 H.Yashiro */

  return(fid);
}

/** add new file structure ********************************************/
static int32_t fio_new_datainfo( int32_t fid )
{
  int32_t did;

  /* get data ID */
  did = finfo[fid].header.num_of_data++;

  /* memory re-allocation (expand by new num_of_file) */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }

  return(did);
}

/** put package information (full) ************************************/
int32_t fio_put_pkginfo( int32_t fid,
                         headerinfo_t hinfo )
{
  int32_t i;

  fio_set_str( finfo[fid].header.description,hinfo.description,FIO_HMID-1  );
  fio_set_str( finfo[fid].header.note,       hinfo.note,       FIO_HLONG-1 );
  finfo[fid].header.num_of_data   = hinfo.num_of_data;
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
headerinfo_t fio_get_pkginfo( int32_t fid )
{
  headerinfo_t hinfo;
  int32_t i;

  fio_trim_str( hinfo.description,finfo[fid].header.description,FIO_HMID-1  );
  fio_trim_str( hinfo.note,       finfo[fid].header.note,FIO_HLONG-1 );
  hinfo.num_of_data   = finfo[fid].header.num_of_data;
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
int32_t fio_put_datainfo( int32_t fid,
                          int32_t did,
                          datainfo_t ditem )
{
  fio_set_str( finfo[fid].dinfo[did].varname,    ditem.varname,    FIO_HSHORT-1 );
  fio_set_str( finfo[fid].dinfo[did].description,ditem.description,FIO_HMID-1   );
  fio_set_str( finfo[fid].dinfo[did].unit,       ditem.unit,       FIO_HSHORT-1 );
  fio_set_str( finfo[fid].dinfo[did].layername,  ditem.layername,  FIO_HSHORT-1 );
  fio_set_str( finfo[fid].dinfo[did].note,       ditem.note,       FIO_HLONG-1  );
  finfo[fid].dinfo[did].datasize     = ditem.datasize;
  finfo[fid].dinfo[did].datatype     = ditem.datatype;
  finfo[fid].dinfo[did].num_of_layer = ditem.num_of_layer;
  finfo[fid].dinfo[did].step         = ditem.step;
  finfo[fid].dinfo[did].time_start   = ditem.time_start;
  finfo[fid].dinfo[did].time_end     = ditem.time_end;

  return(SUCCESS_CODE);
}

/** get data information (full) ***************************************/
datainfo_t fio_get_datainfo( int32_t fid,
                             int32_t did  )
{
  datainfo_t ditem;

  fio_trim_str( ditem.varname,    finfo[fid].dinfo[did].varname,    FIO_HSHORT-1 );
  fio_trim_str( ditem.description,finfo[fid].dinfo[did].description,FIO_HMID-1   );
  fio_trim_str( ditem.unit,       finfo[fid].dinfo[did].unit,       FIO_HSHORT-1 );
  fio_trim_str( ditem.layername,  finfo[fid].dinfo[did].layername,  FIO_HSHORT-1 );
  fio_trim_str( ditem.note,       finfo[fid].dinfo[did].note,       FIO_HLONG-1  );
  ditem.datasize     = finfo[fid].dinfo[did].datasize;
  ditem.datatype     = finfo[fid].dinfo[did].datatype;
  ditem.num_of_layer = finfo[fid].dinfo[did].num_of_layer;
  ditem.step         = finfo[fid].dinfo[did].step;
  ditem.time_start   = finfo[fid].dinfo[did].time_start;
  ditem.time_end     = finfo[fid].dinfo[did].time_end;

  return(ditem);
}

/** seek data id by varname and step **********************************/
int32_t fio_seek_datainfo( int32_t fid,
                           char *varname,
                           int32_t step   )
{
  int32_t did;

  for( did=0; did<finfo[fid].header.num_of_data; did++ ) {
    if (    strcmp(finfo[fid].dinfo[did].varname,varname) == 0
         && finfo[fid].dinfo[did].step    == step              ) {
      return(did);
    }
  }

  return(ERROR_CODE);
}

/** open file IO stream ***********************************************/
int32_t fio_fopen( int32_t fid, int32_t mode )
{
  if (finfo[fid].status.opened) {
    fprintf(stderr," FIle ( %s ) has been already opened!\n",finfo[fid].header.fname);
    fprintf(stderr," open process will be skipped!\n");
    return(SUCCESS_CODE);
  }

  if ( mode==FIO_FWRITE ) {
    if ( (finfo[fid].status.fp=fopen(finfo[fid].header.fname,"wb"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
      exit(1);
    }
  } else if( mode==FIO_FREAD ) { /* [mod] H.Yashiro 20110907 avoid overwrite action */
    if ( (finfo[fid].status.fp=fopen(finfo[fid].header.fname,"rb"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
      exit(1);
    }
  } else if( mode==FIO_FAPPEND ) { /* [add] H.Yashiro 20110907 overwrite mode */
    if ( (finfo[fid].status.fp=fopen(finfo[fid].header.fname,"r+b"))==NULL ) {
      fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
      exit(1);
    }
  }
  finfo[fid].status.rwmode = mode;
  finfo[fid].status.opened = 1;

  return(SUCCESS_CODE);
}

/** close file IO stream **********************************************/
int32_t fio_fclose( int32_t fid )
{
  fclose(finfo[fid].status.fp);
  finfo[fid].status.opened = 0;

  return(SUCCESS_CODE);
}

/** write package information *****************************************/
int32_t fio_write_pkginfo( int32_t fid )
{

  int32_t temp32;
  int32_t *array32;

  if(!finfo[fid].status.opened){
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  fseek(finfo[fid].status.fp,0L,SEEK_SET);
  /* description */
  fwrite(finfo[fid].header.description,sizeof(char),FIO_HMID,finfo[fid].status.fp);
  /* note */
  fwrite(finfo[fid].header.note,sizeof(char),FIO_HLONG,finfo[fid].status.fp);
  /* file mode */
  temp32=finfo[fid].header.fmode;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* endian type */
  temp32 = finfo[fid].header.endiantype;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* grid topology */
  temp32 = finfo[fid].header.grid_topology;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* glevel */
  temp32 = finfo[fid].header.glevel;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* rlevel */
  temp32=finfo[fid].header.rlevel;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* number of region & region id */
  temp32=finfo[fid].header.num_of_rgn;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );

  array32=(int32_t *)malloc(sizeof(int32_t)*finfo[fid].header.num_of_rgn);
  memcpy(array32,finfo[fid].header.rgnid,finfo[fid].header.num_of_rgn*sizeof(int32_t));
  if(system_ednchg){ fio_ednchg(array32,sizeof(int32_t),finfo[fid].header.num_of_rgn); }
  fwrite( array32,sizeof(int32_t),finfo[fid].header.num_of_rgn,finfo[fid].status.fp );
  free(array32);
  /* number of data */
  temp32=finfo[fid].header.num_of_data;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );

  /* remember endpoint of pkginfo */
  /* finfo[fid].status.EOH = ftell(finfo[fid].status.fp);     [del] 20111007 H.Yashiro */
  fgetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */

  return(SUCCESS_CODE);
}

/** read package information ******************************************/
int32_t fio_read_pkginfo( int32_t fid )
{
  int32_t i;
  int32_t temp32;
  int32_t *array32;

  if (!finfo[fid].status.opened) {
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  fseek(finfo[fid].status.fp,0L,SEEK_SET);
  /* description */
  fread(finfo[fid].header.description,sizeof(char),FIO_HMID,finfo[fid].status.fp);
  /* note */
  fread(finfo[fid].header.note,sizeof(char),FIO_HLONG,finfo[fid].status.fp);
  /* file mode */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.fmode = temp32;
  /* endian type */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.endiantype = temp32;
  /* grid topology */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.grid_topology = temp32;
  /* glevel */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.glevel = temp32;
  /* rlevel */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.rlevel = temp32;
  /* number of region & region id */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.num_of_rgn = temp32;

  array32=(int32_t *)malloc(sizeof(int32_t)*temp32);
  fread(array32,sizeof(int32_t),temp32,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(array32,sizeof(int32_t),temp32); }
  finfo[fid].header.rgnid = (int32_t *)realloc(finfo[fid].header.rgnid,
                            finfo[fid].header.num_of_rgn*sizeof(int32_t));
  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    finfo[fid].header.rgnid[i] = array32[i];
  }
  free(array32);
  /* number of data */
  fread(&temp32,sizeof(int32_t),1,finfo[fid].status.fp);
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1); }
  finfo[fid].header.num_of_data = temp32;

  /* remember endpoint of pkginfo */
  /* finfo[fid].status.EOH = ftell(finfo[fid].status.fp);     [del] 20111007 H.Yashiro */
  fgetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */

  return(SUCCESS_CODE);
}

/** write data information ********************************************/
int32_t fio_write_datainfo( int32_t fid,
                            int32_t did  )
{
  int32_t temp32;
  int64_t temp64;

  if(!finfo[fid].status.opened){
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  fseek(finfo[fid].status.fp,0L,SEEK_END);
  /* varname */
  fwrite(finfo[fid].dinfo[did].varname,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
  /* description */
  fwrite(finfo[fid].dinfo[did].description,sizeof(char),FIO_HMID,finfo[fid].status.fp);
  /* unit */
  fwrite(finfo[fid].dinfo[did].unit,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
  /* layername */
  fwrite(finfo[fid].dinfo[did].layername,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
  /* note */
  fwrite(finfo[fid].dinfo[did].note,sizeof(char),FIO_HLONG,finfo[fid].status.fp);
  /* datasize */
  temp64 = finfo[fid].dinfo[did].datasize;
  if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
  fwrite( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
  /* datatype */
  temp32 = finfo[fid].dinfo[did].datatype;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* num_of_layer */
  temp32 = finfo[fid].dinfo[did].num_of_layer;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* step */
  temp32 = finfo[fid].dinfo[did].step;
  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
  fwrite( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
  /* time_start */
  temp64 = finfo[fid].dinfo[did].time_start;
  if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
  fwrite( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
  /* time_end */
  temp64 = finfo[fid].dinfo[did].time_end;
  if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
  fwrite( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );

  return(SUCCESS_CODE);
}

/** read data information *********************************************/
int32_t fio_read_datainfo( int32_t fid )
{
  int32_t did;
  int32_t pos;
  int32_t temp32;
  int64_t temp64;

  if (!finfo[fid].status.opened) {
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  /* read all data information from file */
  fsetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */
  pos = 0;                                                 /* [mod] 20111007 H.Yashiro */
  for( did=0; did<finfo[fid].header.num_of_data; did++ ) {
    fseek(finfo[fid].status.fp,pos,SEEK_CUR); /* [mod] 20111007 H.Yashiro */
    /* varname */
    fread(finfo[fid].dinfo[did].varname,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
    /* description */
    fread(finfo[fid].dinfo[did].description,sizeof(char),FIO_HMID,finfo[fid].status.fp);
    /* unit */
    fread(finfo[fid].dinfo[did].unit,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
    /* layername */
    fread(finfo[fid].dinfo[did].layername,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
    /* note */
    fread(finfo[fid].dinfo[did].note,sizeof(char),FIO_HLONG,finfo[fid].status.fp);
    /* datasize */
    fread( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
    finfo[fid].dinfo[did].datasize = temp64;
    /* datatype */
    fread( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
    finfo[fid].dinfo[did].datatype = temp32;
    /* num_of_layer */
    fread( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
    finfo[fid].dinfo[did].num_of_layer = temp32;
    /* step */
    fread( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
    finfo[fid].dinfo[did].step = temp32;
    /* time_start */
    fread( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
    finfo[fid].dinfo[did].time_start = temp64;
    /* time_end */
    fread( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
    finfo[fid].dinfo[did].time_end = temp64;

    /* skip data array */
    pos = finfo[fid].dinfo[did].datasize; /* [mod] 20111007 H.Yashiro */
  }

  return(SUCCESS_CODE);
}

/** read data information and store data as tmpdata *******************/
/* [add] C.Kodama 13-04-18 */
int32_t fio_read_datainfo_tmpdata( int32_t fid )
{
  int32_t did;
  int32_t pos;
  int32_t temp32;
  int64_t temp64;
  int32_t i, ijklall, v;

  if(tdata == NULL){
    tdata = (struct type_tmpdata *)malloc(sizeof(struct type_tmpdata)*(10*(int32_t)(pow(4,common.rlevel)))); /* common.num_of_rgn >= fid */
    for( i=0; i<common.num_of_rgn; i++ ){ tdata[i].tmpdata = NULL; }

  }
  tdata[fid].tmpdata = (void **)malloc(sizeof(void *)*finfo[fid].header.num_of_data);

  if (!finfo[fid].status.opened) {
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  /* read all data information from file */
  fsetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */
  pos = 0;                                                 /* [mod] 20111007 H.Yashiro */
  fseek(finfo[fid].status.fp,pos,SEEK_CUR);
  for( did=0; did<finfo[fid].header.num_of_data; did++ ) {
    /*fseek(finfo[fid].status.fp,pos,SEEK_CUR);*/ /* [mod] 20111007 H.Yashiro */
    /* varname */
    fread(finfo[fid].dinfo[did].varname,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
    /* description */
    fread(finfo[fid].dinfo[did].description,sizeof(char),FIO_HMID,finfo[fid].status.fp);
    /* unit */
    fread(finfo[fid].dinfo[did].unit,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
    /* layername */
    fread(finfo[fid].dinfo[did].layername,sizeof(char),FIO_HSHORT,finfo[fid].status.fp);
    /* note */
    fread(finfo[fid].dinfo[did].note,sizeof(char),FIO_HLONG,finfo[fid].status.fp);
    /* datasize */
    fread( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
    finfo[fid].dinfo[did].datasize = temp64;
    /* datatype */
    fread( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
    finfo[fid].dinfo[did].datatype = temp32;
    /* num_of_layer */
    fread( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
    finfo[fid].dinfo[did].num_of_layer = temp32;
    /* step */
    fread( &temp32,sizeof(int32_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
    finfo[fid].dinfo[did].step = temp32;
    /* time_start */
    fread( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
    finfo[fid].dinfo[did].time_start = temp64;
    /* time_end */
    fread( &temp64,sizeof(int64_t),1,finfo[fid].status.fp );
    if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
    finfo[fid].dinfo[did].time_end = temp64;

    /* skip data array */
    pos = finfo[fid].dinfo[did].datasize; /* [mod] 20111007 H.Yashiro */

    for( v=0; v<nvar; v++ ){
      if( strcmp( vname[v], finfo[fid].dinfo[did].varname ) == 0 ){
        tdata[fid].tmpdata[did] = (void *)malloc(finfo[fid].dinfo[did].datasize);

        fread( tdata[fid].tmpdata[did], finfo[fid].dinfo[did].datasize, 1, finfo[fid].status.fp );

        ijklall = finfo[fid].dinfo[did].datasize
          / precision[finfo[fid].dinfo[did].datatype];
        if(system_ednchg){
          fio_ednchg(tdata[fid].tmpdata[did],precision[finfo[fid].dinfo[did].datatype],ijklall);
        }
        break;
      }
      if( v == nvar-1 ){
        tdata[fid].tmpdata[did] = NULL;
        fseek(finfo[fid].status.fp,pos,SEEK_CUR);
      }
    }
  }

  return(SUCCESS_CODE);
}

/** write data array **************************************************/
int32_t fio_write_data( int32_t fid,
                        int32_t did,
                        void *data   )
{
  int64_t ijklall;
  void *_data;

  ijklall = finfo[fid].dinfo[did].datasize
          / precision[finfo[fid].dinfo[did].datatype];

  fseek(finfo[fid].status.fp,0L,SEEK_END);
  /* data */
  if(system_ednchg){
    _data = malloc(finfo[fid].dinfo[did].datasize);
    memcpy( _data, data, finfo[fid].dinfo[did].datasize);
    fio_ednchg(_data,precision[finfo[fid].dinfo[did].datatype],ijklall);
  }else{
    _data = data;
  }

  fwrite(_data,finfo[fid].dinfo[did].datasize,1,finfo[fid].status.fp);

  if(system_ednchg) { free(_data); }

  return(SUCCESS_CODE);
}

/** write data array (1 region) ***************************************/
int32_t fio_write_data_1rgn( int32_t fid,
                             int32_t did,
                             void *data   )
{
  int64_t ijkall;
  int64_t datasize;
  void *_data;

  ijkall = (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
                   * (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
         * finfo[fid].dinfo[did].num_of_layer;

  datasize = ijkall * precision[finfo[fid].dinfo[did].datatype];

  fseek(finfo[fid].status.fp,0L,SEEK_END);
  /* data */
  if(system_ednchg){
    _data = malloc(datasize);       /* [fix] 20111028 H.Yashiro */
    memcpy( _data, data, datasize); /* [fix] 20111028 H.Yashiro */
    fio_ednchg(_data,precision[finfo[fid].dinfo[did].datatype],ijkall);
  }else{
    _data = data;
  }

  fwrite(_data,datasize,1,finfo[fid].status.fp);

  if(system_ednchg) { free(_data); }

  return(SUCCESS_CODE);
}

/** read data array (full size) ***************************************/
int32_t fio_read_data( int32_t fid,
                       int32_t did,
                       void *data   )
{
  int64_t i;
  int64_t pos;
  int64_t ijklall;

  ijklall = finfo[fid].dinfo[did].datasize
          / precision[finfo[fid].dinfo[did].datatype];

  fsetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */
  pos = 0;                                                 /* [mod] 20111007 H.Yashiro */
  for( i=0; i<did; i++ ) {
    pos += dinfosize + finfo[fid].dinfo[i].datasize;
  }
  pos += dinfosize;

  fseek(finfo[fid].status.fp,pos,SEEK_CUR);

  fread(data,finfo[fid].dinfo[did].datasize,1,finfo[fid].status.fp);
  if(system_ednchg){
    fio_ednchg(data,precision[finfo[fid].dinfo[did].datatype],ijklall);
  }

  return(SUCCESS_CODE);
}

/** read data array from tmpdata **************************************/
/* [add] C.Kodama 13-04-18 */
int32_t fio_read_data_tmpdata( int32_t fid,
                               int32_t did,
                               void *data   )
{
  int64_t i;
  int64_t pos;
  int64_t ijklall;
  int32_t flag_tmpdata;

  flag_tmpdata = 1;
  if( tdata == NULL ){ flag_tmpdata = 0; }
  else {
    if( tdata[fid].tmpdata == NULL ){ flag_tmpdata = 0; }
    else {
      if( tdata[fid].tmpdata[did] == NULL ){ flag_tmpdata = 0; }
    }
  }
  if( flag_tmpdata == 0 ){ return fio_read_data( fid,did,data ); }
  memcpy( data, tdata[fid].tmpdata[did], finfo[fid].dinfo[did].datasize );
  /* to save memory */
  free( tdata[fid].tmpdata[did] );
  tdata[fid].tmpdata[did] = NULL;

  return(SUCCESS_CODE);
}

/** small functions ***************************************************/
int32_t put_header_fname( int32_t fid, char *fname )
{
  fio_set_str( finfo[fid].header.fname,fname,FIO_HLONG-1 );
  return(SUCCESS_CODE);
}

int32_t put_header_description( int32_t fid, char *description )
{
  fio_set_str( finfo[fid].header.description,description,FIO_HMID-1 );
  return(SUCCESS_CODE);
}

int32_t put_header_note( int32_t fid, char *note )
{
  fio_set_str( finfo[fid].header.note,note,FIO_HLONG-1 );
  return(SUCCESS_CODE);
}

int32_t put_header_fmode( int32_t fid, int32_t fmode )
{
  finfo[fid].header.fmode = fmode;
  return(SUCCESS_CODE);
}

int32_t put_header_endiantype( int32_t fid, int32_t endiantype )
{
  finfo[fid].header.endiantype = endiantype;
  return(SUCCESS_CODE);
}

int32_t put_header_grid_topology( int32_t fid, int32_t grid_topology )
{
  finfo[fid].header.grid_topology = grid_topology;
  return(SUCCESS_CODE);
}

int32_t put_header_glevel( int32_t fid, int32_t glevel )
{
  finfo[fid].header.glevel = glevel;
  return(SUCCESS_CODE);
}

int32_t put_header_rlevel( int32_t fid, int32_t rlevel )
{
  finfo[fid].header.rlevel = rlevel;
  return(SUCCESS_CODE);
}

int32_t put_header_rgn( int32_t fid, int32_t num_of_rgn, int32_t rgnid[] )
{
  int32_t i;

  finfo[fid].header.num_of_rgn = num_of_rgn;

  finfo[fid].header.rgnid = (int32_t *)realloc(finfo[fid].header.rgnid,
                            finfo[fid].header.num_of_rgn*sizeof(int32_t));

  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    finfo[fid].header.rgnid[i] = rgnid[i];
  }
  return(SUCCESS_CODE);
}

int32_t put_header_num_of_data( int32_t fid, int32_t num_of_data )
{
  finfo[fid].header.num_of_data = num_of_data;
  return(SUCCESS_CODE);
}

int32_t get_header_fname( int32_t fid, char *fname )
{
  fio_set_str( fname,finfo[fid].header.fname,FIO_HLONG-1 );
  return(SUCCESS_CODE);
}

int32_t get_header_description( int32_t fid, char *description )
{
  fio_set_str( description,finfo[fid].header.description,FIO_HMID-1 );
  return(SUCCESS_CODE);
}

int32_t get_header_note( int32_t fid, char *note )
{
  fio_set_str( note,finfo[fid].header.note,FIO_HLONG-1 );
  return(SUCCESS_CODE);
}

int32_t get_header_fmode( int32_t fid, int32_t fmode )
{
  fmode = finfo[fid].header.fmode;
  return(SUCCESS_CODE);
}

int32_t get_header_endiantype( int32_t fid, int32_t endiantype )
{
  endiantype = finfo[fid].header.fmode;
  return(SUCCESS_CODE);
}

int32_t get_header_grid_topology( int32_t fid, int32_t grid_topology )
{
  grid_topology = finfo[fid].header.grid_topology;
  return(SUCCESS_CODE);
}

int32_t get_header_glevel( int32_t fid, int32_t glevel )
{
  glevel = finfo[fid].header.glevel;
  return(SUCCESS_CODE);
}

int32_t get_header_rlevel( int32_t fid, int32_t rlevel )
{
  rlevel = finfo[fid].header.rlevel;
  return(SUCCESS_CODE);
}

int32_t get_header_rgn( int32_t fid, int32_t num_of_rgn, int32_t *rgnid )
{
  int32_t i;

  num_of_rgn = finfo[fid].header.num_of_rgn;

  rgnid = (int32_t *)malloc(num_of_rgn*sizeof(int32_t));

  for( i=0; i<num_of_rgn; i++ ) {
     rgnid[i]= finfo[fid].header.rgnid[i];
  }
  return(SUCCESS_CODE);
}

int32_t get_header_num_of_data( int32_t fid, int32_t num_of_data )
{
  num_of_data = finfo[fid].header.num_of_data;
  return(SUCCESS_CODE);
}

/** register new file *************************************************/
int32_t fio_register_file( char *fname )
{
  int32_t ierr;
  int32_t fid;

  /* request new file space */
  fid = fio_new_finfo();

  fio_set_str( finfo[fid].header.fname,fname,FIO_HLONG-1  );

  return(fid);
}

/** put & write package information (quick put) ***********************/
int32_t fio_put_write_pkginfo( int32_t fid,
                               char *description,
                               char *note         )
{
  int32_t i;

  fio_set_str( finfo[fid].header.description,description,FIO_HMID-1  );
  fio_set_str( finfo[fid].header.note,       note,       FIO_HLONG-1 );
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

  fio_write_pkginfo( fid );

  return(SUCCESS_CODE);
}

/** validate package information with common **************************/
int32_t fio_valid_pkginfo( int32_t fid )
{
  int32_t i;

  if(finfo[fid].header.grid_topology!=common.grid_topology) {
    fprintf(stderr,"Warning: grid_topology is not match, %d, %d\n",
                   finfo[fid].header.grid_topology,common.grid_topology);
  }

  if(finfo[fid].header.glevel!=common.glevel) {
    fprintf(stderr,"Warning: glevel is not match, %d, %d\n",
                   finfo[fid].header.glevel,common.glevel              );
  }

  if(finfo[fid].header.rlevel!=common.rlevel) {
    fprintf(stderr,"Warning: rlevel is not match, %d, %d\n",
                   finfo[fid].header.rlevel,common.rlevel              );
  }

  if(finfo[fid].header.num_of_rgn!=common.num_of_rgn) {
    fprintf(stderr,"Warning: num_of_rgn is not match, %d, %d\n",
                   finfo[fid].header.num_of_rgn,common.num_of_rgn      );
  }

  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    if(finfo[fid].header.rgnid[i]!=common.rgnid[i]) {
    fprintf(stderr,"Warning: rgnid[%d] is not match, %d, %d\n",
                   i,finfo[fid].header.rgnid[i],common.rgnid[i]        );
    }
  }
  return(SUCCESS_CODE);
}

/** validate package information with common **************************/
int32_t fio_valid_pkginfo_validrgn( int32_t fid,
                                    int32_t rgnid[] )
{
  int32_t i;

  if(finfo[fid].header.grid_topology!=common.grid_topology) {
    fprintf(stderr,"Warning: grid_topology is not match, %d, %d\n",
                   finfo[fid].header.grid_topology,common.grid_topology);
  }

  if(finfo[fid].header.glevel!=common.glevel) {
    fprintf(stderr,"Warning: glevel is not match, %d, %d\n",
                   finfo[fid].header.glevel,common.glevel              );
  }

  if(finfo[fid].header.rlevel!=common.rlevel) {
    fprintf(stderr,"Warning: rlevel is not match, %d, %d\n",
                   finfo[fid].header.rlevel,common.rlevel              );
  }

  if(finfo[fid].header.num_of_rgn!=common.num_of_rgn) {
    fprintf(stderr,"Warning: num_of_rgn is not match, %d, %d\n",
                   finfo[fid].header.num_of_rgn,common.num_of_rgn      );
  }

  for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
    if(finfo[fid].header.rgnid[i]!=rgnid[i]) {
    fprintf(stderr,"Warning: rgnid[%d] is not match, %d, %d\n",
                   i,finfo[fid].header.rgnid[i],rgnid[i]               );
    }
  }
  return(SUCCESS_CODE);
}

/** validate data size ************************************************/
int32_t fio_valid_datainfo( int32_t fid )
{
  int32_t did;
  int32_t ijall;
  int64_t datasize;
  char* str_dtype[5] = { "XXXX","REAL4","REAL8","INTEGER4","INTEGER8" };

  ijall = (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
                  * (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2);

  for( did=0; did<finfo[fid].header.num_of_data; did++ ) {
    datasize = ijall
             * finfo[fid].dinfo[did].num_of_layer
             * finfo[fid].header.num_of_rgn
             * precision[finfo[fid].dinfo[did].datatype];

    if(finfo[fid].dinfo[did].datasize!=datasize) {
      fprintf(stderr,"Warning: datasize is not match, %ld\n",
                     (long)finfo[fid].dinfo[did].datasize     );
      fprintf(stderr,"         datasize must be %d[grid]x%d[layer]x%d[region]x%s=%ld\n",
                     ijall,
                     finfo[fid].dinfo[did].num_of_layer,
                     finfo[fid].header.num_of_rgn,
                     str_dtype[finfo[fid].dinfo[did].datatype+1],
                     (long)datasize                                             );
    }

  }
  return(SUCCESS_CODE);
}

/** put & write data information and write data ***********************/
int32_t fio_put_write_datainfo_data( int32_t fid,
                                     datainfo_t ditem,
                                     void *data        )
{
  int32_t did;

  did = fio_new_datainfo( fid );

  fio_put_datainfo( fid, did, ditem );

  fio_write_pkginfo( fid ); /* update num_of_data */
  fio_write_datainfo( fid, did );
  fio_write_data( fid, did, data );

  return(did);
}

/** put & write data information **************************************/
int32_t fio_put_write_datainfo( int32_t fid,
                                datainfo_t ditem )
{
  int32_t did;

  did = fio_new_datainfo( fid );

  fio_put_datainfo( fid, did, ditem );

  fio_write_pkginfo( fid ); /* update num_of_data */
  fio_write_datainfo( fid, did );

  return(did);
}

/** read pkginfo and datainfo *****************************************/
int32_t fio_read_allinfo( int32_t fid )
{
  int32_t i;

  fio_read_pkginfo( fid );
  fio_valid_pkginfo( fid );

  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }

  fio_read_datainfo( fid );
  fio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

/** read pkginfo and datainfo, with validating rgnid ******************/
int32_t fio_read_allinfo_validrgn( int32_t fid,
                                   int32_t rgnid[] )
{
  int32_t i;

  fio_read_pkginfo( fid );
  fio_valid_pkginfo_validrgn( fid, rgnid );

  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }

  fio_read_datainfo( fid );
  fio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

/** read pkginfo and datainfo and store data as tmpdata ***************/
/* [add] C.Kodama 13-04-18 */
int32_t fio_read_allinfo_tmpdata( int32_t fid )
{
  int32_t i;

  fio_read_pkginfo( fid );
  fio_valid_pkginfo( fid );

  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }

  fio_read_datainfo_tmpdata( fid );
  fio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

void fio_register_vname_tmpdata( const char *vname_in )
{
  strcpy(vname[nvar++],vname_in);
  /* printf("debug:%s\n", vname[nvar-1]);*/
}

/** allocate and copy datainfo ****************************************/
/* [add] C.Kodama 13-04-18 */
int32_t fio_copy_datainfo( int32_t fid, int32_t fid_org )
{
  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }

  /*fio_read_datainfo( fid );*/
  memcpy( finfo[fid].dinfo, finfo[fid_org].dinfo, sizeof(datainfo_t)*finfo[fid].header.num_of_data );
  fio_valid_datainfo( fid );

  return(SUCCESS_CODE);
}

/** dump package summary of all finfo *********************************/
int32_t fio_dump_finfolist( void )
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
  for( fid=0; fid<num_of_file; fid++ ) {
    printf( " %4s", str_mode[finfo[fid].header.fmode+1] );
    printf( " %3s", str_rw[finfo[fid].status.rwmode+1] );
    printf( " %4s", str_opened[finfo[fid].status.opened+1] );
    printf( " %6s", str_endian[finfo[fid].header.endiantype+1] );
    printf( " %4s", str_topology[finfo[fid].header.grid_topology+1] );
    printf( " %2d", finfo[fid].header.glevel );
    printf( " %2d", finfo[fid].header.rlevel );
    printf( " %4d", finfo[fid].header.num_of_rgn );
    printf( " %5d", finfo[fid].header.num_of_data );
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
int32_t fio_dump_finfo( int32_t fid,
                        int32_t endiantype,
                        int32_t dumptype    )
{
  int32_t i,ij,k,l;
  int32_t did;
  int32_t ijall;
  int64_t ijklall,pos;
  void *buf;
  char *c;

  char* str_mode[3]     = { "XXXX","SPRIT","COMPLETE" };
  char* str_endian[4]   = { "XXXXXX","UKNOWN","LITTLE ENDIAN","BIG ENDIAN" };
  char* str_topology[4] = { "XXXX","ICOSAHEDRON","LCP","MLCP" };
  char* str_dtype[5]    = { "XXXX","REAL4","REAL8","INTEGER4","INTEGER8" };

  /* exchange endian? */
  if ( endiantype!=system_endiantype ) {
    system_ednchg = 1;
    printf("%s to %s\n",str_endian[endiantype+1],str_endian[system_endiantype+1]);
  }

  fio_fopen( fid, FIO_FREAD );
  fio_read_pkginfo ( fid );
  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }
  fio_read_datainfo( fid );

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
  printf( "--- NUMBER OF ITEMS     : %d\n", finfo[fid].header.num_of_data );
  printf( "--- DATA ENDIAN TYPE    : %s\n", str_endian[finfo[fid].header.endiantype+1] );
  printf( "--- NOTE                : %s\n", finfo[fid].header.note );
  printf( "\n" );

  printf("============ ITEM ATTRIBUTES ======================\n");
  fsetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */
  pos = dinfosize;                                         /* [mod] 20111007 H.Yashiro */
  for( did=0; did<finfo[fid].header.num_of_data; did++ ) {
    if ( dumptype == FIO_DUMP_ALL || dumptype == FIO_DUMP_ALL_MORE ) {
      /* data info */
      printf( "########## Item ID = %4d ###############\n",did);
      printf( "--- VAR NAME            : %s\n",  finfo[fid].dinfo[did].varname );
      printf( "--- STEP NO.            : %d\n",  finfo[fid].dinfo[did].step );
      printf( "--- DATA DESCRIPTION    : %s\n",  finfo[fid].dinfo[did].description );
      printf( "--- UNIT                : %s\n",  finfo[fid].dinfo[did].unit );
      printf( "--- VERT. LAYER NAME    : %s\n",  finfo[fid].dinfo[did].layername );
      printf( "--- NUMBER OF LAYERS    : %d\n",  finfo[fid].dinfo[did].num_of_layer );
      printf( "--- START TIME [sec]    : %ld\n", (long)finfo[fid].dinfo[did].time_start );
      printf( "---   END TIME [sec]    : %ld\n", (long)finfo[fid].dinfo[did].time_end );
      printf( "--- DATA TYPE           : %s\n",  str_dtype[finfo[fid].dinfo[did].datatype+1] );
      printf( "--- DATA SIZE [byte]    : %ld\n", (long)finfo[fid].dinfo[did].datasize );
      printf( "--- NOTE                : %s\n",  finfo[fid].dinfo[did].note );
    } else {
    /* data info */
    printf( "+%4d:%16s|%6d|%16s[%3d]|%12ld-%12ld|%10ldbyte[%s]\n",
            did,
            finfo[fid].dinfo[did].varname,
            finfo[fid].dinfo[did].step,
            finfo[fid].dinfo[did].layername,
            finfo[fid].dinfo[did].num_of_layer,
            (long)finfo[fid].dinfo[did].time_start,
            (long)finfo[fid].dinfo[did].time_end,
            (long)finfo[fid].dinfo[did].datasize,
            str_dtype[finfo[fid].dinfo[did].datatype+1] );
    }

    /* data array */
    if ( dumptype == FIO_DUMP_ALL ) {
      ijklall = ijall
              * finfo[fid].dinfo[did].num_of_layer
              * finfo[fid].header.num_of_rgn;

      /* skip datainfo */
      fseek(finfo[fid].status.fp,pos,SEEK_CUR); /* [mod] 20111007 H.Yashiro */

      /* read data */
      buf = malloc(finfo[fid].dinfo[did].datasize);
      fread(buf,finfo[fid].dinfo[did].datasize,1,finfo[fid].status.fp);
      if(system_ednchg){
        fio_ednchg(buf,precision[finfo[fid].dinfo[did].datatype],ijklall);
      }

      /* dump data */
      c = (char *)buf;
      for( l=0; l<finfo[fid].header.num_of_rgn; l++ ) {
      for( k=0; k<finfo[fid].dinfo[did].num_of_layer; k++ ) {
      for( ij=0; ij<ijall; ij++ ) {
        printf("+++ [%8s, %4d, %6d, %3d, %6d] ",
        finfo[fid].dinfo[did].varname,finfo[fid].dinfo[did].step,ij,k,l);

        switch(finfo[fid].dinfo[did].datatype){
        case FIO_REAL4:
               printf("%f\n",*(real32_t *)c);
               c += 4;
               break;
        case FIO_REAL8:
               printf("%lf\n",*(real64_t *)c);
               c += 8;
               break;
        case FIO_INTEGER4:
               printf("%d\n",*(int32_t *)c);
               c += 4;
               break;
        case FIO_INTEGER8:
               printf("%lld\n",*(int64_t *)c);
               c += 8;
               break;
        default :
               fprintf(stderr,"xxx Undefined data type! stop\n");
               break;
        }
      }
      }
      }
    }

    /* data array (more large digit)*/
    if ( dumptype == FIO_DUMP_ALL_MORE ) {
      ijklall = ijall
              * finfo[fid].dinfo[did].num_of_layer
              * finfo[fid].header.num_of_rgn;

      /* skip datainfo */
      fseek(finfo[fid].status.fp,pos,SEEK_CUR); /* [mod] 20111007 H.Yashiro */

      /* read data */
      buf = malloc(finfo[fid].dinfo[did].datasize);
      fread(buf,finfo[fid].dinfo[did].datasize,1,finfo[fid].status.fp);
      if(system_ednchg){
        fio_ednchg(buf,precision[finfo[fid].dinfo[did].datatype],ijklall);
      }

      /* dump data */
      c = (char *)buf;
      for( l=0; l<finfo[fid].header.num_of_rgn; l++ ) {
      for( k=0; k<finfo[fid].dinfo[did].num_of_layer; k++ ) {
      for( ij=0; ij<ijall; ij++ ) {
        printf("+++ [%8s, %4d, %6d, %3d, %6d] ",
        finfo[fid].dinfo[did].varname,finfo[fid].dinfo[did].step,ij,k,l);

        switch(finfo[fid].dinfo[did].datatype){
        case FIO_REAL4:
               printf("%.60e\n",*(real32_t *)c);
               c += 4;
               break;
        case FIO_REAL8:
               printf("%.60le\n",*(real64_t *)c);
               c += 8;
               break;
        case FIO_INTEGER4:
               printf("%d\n",*(int32_t *)c);
               c += 4;
               break;
        case FIO_INTEGER8:
               printf("%lld\n",*(int64_t *)c);
               c += 8;
               break;
        default :
               fprintf(stderr,"xxx Undefined data type! stop\n");
               break;
        }
      }
      }
      }
    }

  }
  printf("============ END ITEM ATTRIBUTES ==================\n");
  fio_fclose( fid );

  return(SUCCESS_CODE);
}

