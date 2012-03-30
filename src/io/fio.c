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
 *    1.30      12-03-02  A.Shimada : Apply HDF5 format               *
 *                                                                    *
 **********************************************************************
 *functions
 * : put/get indicates the communication with the database
 *   read/write indicates the communication with the file
 *   MPI-IO is not supported yet.
 *<utilities>
 *void           fio_ednchg         : endian changer
 *void           fio_mk_fname       : filename generator
 *void           fio_set_str        : string preprocess
 *<database>
 *int32_t        fio_syscheck       : check system
 *int32_t        fio_put_commoninfo : store common informtation
 *static int32_t fio_new_finfo      : add new file structure
 *static int32_t fio_new_datainfo   : add new datainfo structure
 *int32_t        fio_put_pkginfo    : put package information (full)
 *int32_t        fio_get_pkginfo    : get package information (full)
 *int32_t        fio_put_datainfo   : put data information (full)
 *int32_t        fio_get_datainfo   : get data information (full)
 *int32_t        fio_seek_datainfo  : seek data id by varname and step
 *<file r/w>
 *int32_t        fio_fopen           : open file IO stream
 *int32_t        fio_fclose          : close file IO stream
 *int32_t        fio_write_pkginfo   : write package information
 *int32_t        fio_read_pkginfo    : read package information
 *int32_t        fio_write_datainfo  : write data information
 *int32_t        fio_read_datainfo   : read data information
 *int32_t        fio_write_data      : write data array
 *int32_t        fio_write_data_1rgn : write data array (1 region)
 *int32_t        fio_read_data       : read data array
 *<function suite for fortran program>
 *int32_t        fio_register_file            : register new file
 *void           fio_new_group                : create new group
 *void           fio_close_group              : close a group
 *int32_t        fio_put_write_pkginfo        : put & write package information (quick put)
 *int32_t        fio_valid_pkginfo            : validate package information with common
 *int32_t        fio_valid_datainfo           : validate data size
 *int32_t        fio_put_write_datainfo_data  : put & write data information and write data
 *int32_t        fio_put_write_datainfo       : put & write data information
 *int32_t        fio_read_allinfo_get_pkginfo : read pkginfo and datainfo and get pkginfo
 *int32_t        fio_dump_finfolist           : dump package summary of all finfo
 *int32_t        fio_dump_finfo               : dump package detail of finfo
 **********************************************************************/
#ifdef CONFIG_HDF5
#include "hdf5.h"
#endif
#include "fio.h"

double ch_elapsed_time = 0;
double cm_elapsed_time = 0;
double io_elapsed_time = 0;
double sdp_elapsed_time = 0;
unsigned long zero = 0;

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

  if ( _str[str_len-1]!=' ' ) {
    _str[str_len] = '\0';
  } else {
    for( i=str_len-1; i>=0; i-- ) {
      if( _str[i] == ' ' ) {
        _str[i] = '\0';
      } else {
        break;
      }
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

  if ( *((char*)&i) ) {
    system_endiantype = FIO_LITTLE_ENDIAN;
/*    printf("This system is Little endian.\n"); */
  } else {
    system_endiantype = FIO_BIG_ENDIAN;
/*    printf("This system is Big endian.\n"); */
  }

  /* intitialize */
  common.use_mpiio     = -1;
  common.fmode         = -1;
  common.endiantype    = -1;
  common.grid_topology = -1;
  common.glevel        = -1;
  common.rlevel        = -1;
#ifdef CONFIG_HDF5
  common.rlevel_i      = -1;
  common.rlevel_j      = -1;
#endif
  common.num_of_rgn    =  0;
  common.rgnid         = NULL;

  return(SUCCESS_CODE);
}

/** put common informtation *******************************************/
int32_t fio_put_commoninfo( int32_t use_mpiio,
                            int32_t fmode,
                            int32_t endiantype,
                            int32_t grid_topology,
                            int32_t glevel,
                            int32_t rlevel,
#ifdef CONFIG_HDF5
                            int32_t rlevel_i,
                            int32_t rlevel_j,
#endif
                            int32_t num_of_rgn,
                            int32_t rgnid[]        )
{
  int32_t i;

  common.use_mpiio     = use_mpiio;
  common.fmode         = fmode;
  common.endiantype    = endiantype;
  common.grid_topology = grid_topology;
  common.glevel        = glevel;
  common.rlevel        = rlevel;
#ifdef CONFIG_HDF5
  common.rlevel_i      = rlevel_i;
  common.rlevel_j      = rlevel_j;
#endif
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
#ifdef CONFIG_HDF5
  finfo[fid].status.file      = -1;
  finfo[fid].tmp_gid          = -1;
#endif
  /* finfo[fid].status.eoh.__pos = 0; /* [add] 20111007 H.Yashiro */
#ifndef NO_MPIIO
  finfo[fid].status.mpi_eoh   = 0;
  finfo[fid].status.mpi_fid   = NULL;
#endif

  return(fid);
}

#ifdef CONFIG_HDF5
/** create new group and return gid ***********************************/
void fio_new_group(int32_t fid, int i, int y)
{
	hid_t gid;
	char gname[FIO_HLONG];

	printf("[DEBUG] %s\n", __func__);

	switch(y) {
	case 4 :
		sprintf(gname,"%s%04d","step",i);
		break;
	case 5 :
		sprintf(gname,"%s%05d","step",i);
		break;
	case 6 :
		sprintf(gname,"%s%06d","step",i);
		break;
	default :
		break;
	}
	
	gid = H5Gcreate(finfo[fid].status.file, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	finfo[fid].tmp_gid = (int32_t)gid;
/*	if(gid > 0) {
		finfo[fid].tmp_gid = (int32_t)gid;
		return (int32_t)gid;
	} else {
		return finfo[fid].tmp_gid;
	}
*/
}

/** get gid **********************************************************/
int32_t fio_get_group(int32_t fid)
{
	printf("[DEBUG] %s\n", __func__);
	return (int32_t)(finfo[fid].tmp_gid);
}

/** close group ******************************************************/
void fio_close_group(int32_t fid)
{
	printf("[DEBUG] %s\n", __func__);
	H5Gclose((hid_t)(finfo[fid].tmp_gid));
	finfo[fid].tmp_gid = -1;
}
#endif

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
#ifdef CONFIG_HDF5
	hid_t file;
#endif
  if (finfo[fid].status.opened) {
    fprintf(stderr," FIle ( %s ) has been already opened!\n",finfo[fid].header.fname);
    fprintf(stderr," open process will be skipped!\n");
    return(SUCCESS_CODE);
  }

  if(common.use_mpiio) {
    /*    if(mode==FIO_FCREATE){
      ierr = MPI_File_delete(finfo[fid].header.fname,mpi_info);
    }
    ierr = MPI_File_open ( MPI_COMM_WORLD,finfo[fid].header.fname, 
			   MPI_MODE_RDWR | MPI_MODE_CREATE, 
			   MPI_INFO_NULL, &(finfo[fid].status.mpi_fid));
    if(ierr!=0){
      fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
      return(ERROR_CODE);
    }
    finfo[fid].status.opened = 1;
    */
  } else {
    if ( mode==FIO_FWRITE ) {
#ifdef CONFIG_HDF5
      file = H5Fopen(finfo[fid].header.fname, H5F_ACC_RDWR, H5P_DEFAULT);
      if(file < 0) {
	      file = H5Fcreate(finfo[fid].header.fname, 
			       H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	      if(file < 0)
		      exit(1);
      }
      finfo[fid].status.file = (int32_t)file;
#else
      if ( (finfo[fid].status.fp=fopen(finfo[fid].header.fname,"wb"))==NULL ) {
        fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
        exit(1);
      }
#endif
    } else if( mode==FIO_FREAD ) { /* [mod] H.Yashiro 20110907 avoid overwrite action */
#ifdef CONFIG_HDF5_INPUT
	    file = H5Fopen(finfo[fid].header.fname, H5F_ACC_RDWR, H5P_DEFAULT);
	    if(file < 0) {
		    fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
		    exit(1);
	    }
	    finfo[fid].status.file = (int32_t)file;
#else
      if ( (finfo[fid].status.fp=fopen(finfo[fid].header.fname,"rb"))==NULL ) {
        fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
        exit(1);
      }
#endif
    } else if( mode==FIO_FAPPEND ) { /* [add] H.Yashiro 20110907 overwrite mode */
      if ( (finfo[fid].status.fp=fopen(finfo[fid].header.fname,"r+b"))==NULL ) {
        fprintf(stderr,"Can not open file : %s!\n",finfo[fid].header.fname);
        exit(1);
      }
    }
    finfo[fid].status.rwmode = mode;
    finfo[fid].status.opened = 1;
  }
  return(SUCCESS_CODE);
}

/** close file IO stream **********************************************/
int32_t fio_fclose( int32_t fid )
{
  if(common.use_mpiio) {
    /*    MPI_File_close(&(finfo[fid].status.mpi_fid));
    finfo[fid].status.opened=0;
    */
  } else {
#ifdef CONFIG_HDF5
    H5Fclose(finfo[fid].status.file);
    if(finfo[fid].status.fp != NULL)
	    fclose(finfo[fid].status.fp);
#else
    fclose(finfo[fid].status.fp);
#endif
    finfo[fid].status.opened = 0;
  }
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

  if(common.use_mpiio) {
    /*
    MPI_File_seek(finfo[fid].status.mpi_fid,0,MPI_SEEK_SET);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    if(myid==0){
      ** file description **
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    finfo[fid].header.description,
			    MAX_CHAR,MPI_CHAR,
			    &status);
      ** endian_type **
      dummy32=finfo[fid].header.endiantype;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      ** grid_topology **
      dummy32=finfo[fid].header.grid_topology;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      ** glevel **
      dummy32=finfo[fid].header.glevel;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      ** rlevel **
      dummy32=finfo[fid].header.rlevel;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      ** number of region **
      num_of_rgn = pow(4,finfo[fid].header.rlevel)*10;
      dummy32=num_of_rgn;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      ** region id **
      rgnid=(int32_t *)malloc(sizeof(int32_t)*num_of_rgn);
      for(i=0;i<num_of_rgn;i++){
	rgnid[i] = i;
      }
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg(rgnid, sizeof(int32_t),num_of_rgn);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    rgnid,
			    num_of_rgn,MPI_INTEGER,
			    &status);
      free(rgnid);
      ** number of data **
      dummy32=finfo[fid].header.num_of_data;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      ** File mode  **
      dummy32=finfo[fid].header.fmode;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);

      MPI_File_get_position(finfo[fid].status.mpi_fid,&(finfo[fid].status.mpi_eoh));
    }
    MPI_Bcast(&(finfo[fid].status.mpi_eoh),sizeof(MPI_Offset), MPI_BYTE, 0, MPI_COMM_WORLD);
    */
  } else {
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
  }
  return(SUCCESS_CODE);
}

/** read package information ******************************************/
int32_t fio_read_pkginfo( int32_t fid )
{
  int32_t i;
  int32_t temp32;
  int32_t *array32;
#ifdef CONFIG_HDF5_INPUT
  hid_t attr;
  hid_t type;
#endif

  if (!finfo[fid].status.opened) {
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }
  if(common.use_mpiio) {
    /*
    MPI_File_seek(finfo[fid].status.mpi_fid,0,MPI_SEEK_SET);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    if(myid==0){
      ** file description **
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   finfo[fid].header.description,
			   MAX_CHAR,MPI_CHAR,
			   &status);
      ** endian type *
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(finfo[fid].header.endiantype),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &(finfo[fid].header.endiantype), sizeof(int32_t),1);	
      }
      * grid topology *
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(finfo[fid].header.grid_topology),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &(finfo[fid].header.grid_topology), sizeof(int32_t),1);	
      }
      * glevel *
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(finfo[fid].header.glevel),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &(finfo[fid].header.glevel), sizeof(int32_t),1);	
      }
      * rlevel *
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(finfo[fid].header.rlevel),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &(finfo[fid].header.rlevel), sizeof(int32_t),1);	
      }
      * number of region *
      num_of_rgn = pow(4,finfo[fid].header.rlevel)*10;
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(num_of_rgn),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &num_of_rgn, sizeof(int32_t),1);	
      }
      * region id *
      rgnid=(int32_t *)malloc(sizeof(int32_t)*num_of_rgn);
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   rgnid,
			   num_of_rgn,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg(rgnid, sizeof(int32_t),num_of_rgn);
      }
      free(rgnid);
      * number of data *
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(finfo[fid].header.num_of_data),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &(finfo[fid].header.num_of_data),sizeof(int32_t),1);
      }
      * File mode  *
      ierr = MPI_File_read(finfo[fid].status.mpi_fid,
			   &(finfo[fid].header.fmode),
			   1,MPI_INTEGER,
			   &status);
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &(finfo[fid].header.fmode),sizeof(int32_t),1);
      }

      MPI_File_get_position(finfo[fid].status.mpi_fid,&(finfo[fid].status.EOH));
    }
    MPI_Bcast(&(finfo[fid].status.EOH),sizeof(MPI_Offset), MPI_BYTE, 0, MPI_COMM_WORLD);
    */
  } else {
#ifdef CONFIG_HDF5_INPUT
    /* description */
	attr = H5Aopen(finfo[fid].status.file, "description", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, finfo[fid].header.description);
	H5Aclose(attr);
    /* note */
	attr = H5Aopen(finfo[fid].status.file, "note", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, finfo[fid].header.note);
	H5Aclose(attr);
    /* file mode */
	attr = H5Aopen(finfo[fid].status.file, "fmode", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.fmode = temp32;
	H5Aclose(attr);
    /* endian type */
	attr = H5Aopen(finfo[fid].status.file, "endiantype", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.endiantype = temp32;
	H5Aclose(attr);
    /* grid topology */
	attr = H5Aopen(finfo[fid].status.file, "grid_topology", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.grid_topology = temp32;
	printf("grid_topology = %d\n", temp32);
	H5Aclose(attr);
    /* glevel */
	attr = H5Aopen(finfo[fid].status.file, "glevel", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.glevel = temp32;
	H5Aclose(attr);
    /* rlevel */
	attr = H5Aopen(finfo[fid].status.file, "rlevel", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.rlevel = temp32;
	H5Aclose(attr);
    /* number of region & region id */
	attr = H5Aopen(finfo[fid].status.file, "num_of_rgn", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.num_of_rgn = temp32;
	H5Aclose(attr);

	array32=(int32_t *)malloc(sizeof(int32_t)*temp32);	
	H5Aopen(finfo[fid].status.file, "rgnid", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, array32);
	if(system_ednchg) {
		fio_ednchg(array32,sizeof(int32_t),temp32);
	}
	finfo[fid].header.rgnid = (int32_t *)realloc(finfo[fid].header.rgnid,
						     finfo[fid].header.num_of_rgn*sizeof(int32_t));
	for( i=0; i<finfo[fid].header.num_of_rgn; i++ ) {
		finfo[fid].header.rgnid[i] = array32[i];
	}
	free(array32);
    /* number of data */
	attr = H5Aopen(finfo[fid].status.file, "num_of_data", H5P_DEFAULT);
	type = H5Aget_type(attr);
	H5Aread(attr, type, &temp32);
	if(system_ednchg) {
		fio_ednchg(&temp32, sizeof(int32_t),1);
	}
	finfo[fid].header.num_of_data = temp32;
	H5Aclose(attr);
#else
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
#endif
  }
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

  if(common.use_mpiio) {

  } else {
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
  }
  return(SUCCESS_CODE);
}

/** read data information *********************************************/
int32_t fio_read_datainfo( int32_t fid )
{
  int32_t did;
  int32_t pos;
  int32_t temp32;
  int64_t temp64;
#ifdef CONFIG_HDF5_INPUT
  hid_t gid, dset, attr, type;
  char dname[FIO_HSHORT];
#endif

  if (!finfo[fid].status.opened) {
    fprintf(stderr,"%s is not open!\n",finfo[fid].header.fname);
    return(ERROR_CODE);
  }

  if(common.use_mpiio) {
    /*
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_File_seek(finfo[fid].status.mpi_fid,0,MPI_SEEK_SET);
    if(myid==0){
      pos = finfo[fid].status.EOH;
      di_headersize=sizeof(char)*MAX_CHAR+sizeof(int32_t)*3+sizeof(int64_t);
      for(i=0;i<finfo[fid].header.num_of_data-1;i++){
	pos += ( di_headersize
		 +precision[finfo[fid].datainfo[i].datatype]
		 *(pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
		 *(pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
		 *finfo[fid].datainfo[i].num_of_layer*pow(4,finfo[fid].header.rlevel)*10 );
	* <- The extension for another topology is needed. Currently, not supported!*
      }
      * set stating point of dataset *
      MPI_File_seek(finfo[fid].status.mpi_fid,pos,MPI_SEEK_SET);

      * varname *
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    finfo[fid].datainfo[did].varname,
			    MAX_CHAR,MPI_CHAR,
			    &status);
      * datatype *
      dummy32=finfo[fid].datainfo[did].datatype;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);	
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      * num of layer  *
      dummy32=finfo[fid].datainfo[did].num_of_layer;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);	
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      * time  *
      dummy32=finfo[fid].datainfo[did].time;
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &dummy32, sizeof(int32_t),1);	
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &dummy32,
			    1,MPI_INTEGER,
			    &status);
      * datasize  *
      di_datasize
	=precision[finfo[fid].datainfo[did].datatype]
	*(pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
	*(pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
	*finfo[fid].datainfo[did].num_of_layer*pow(4,finfo[fid].header.rlevel)*10;
      * <- The extension for another topology is needed. Currently, not supported!*
      if(finfo[fid].header.endiantype!=fio_system_endian){
	fio_ednchg( &di_datasize, sizeof(int64_t),1);	
      }
      ierr = MPI_File_write(finfo[fid].status.mpi_fid,
			    &(di_datasize),
			    1,MPI_INTEGER8,
			    &status);
      MPI_File_get_position(finfo[fid].status.mpi_fid,&(offset0));
    }
    MPI_Bcast(&offset0,sizeof(MPI_Offset), MPI_BYTE, 0, MPI_COMM_WORLD);

    if(finfo[fid].datainfo[did].datatype==FIO_REAL4){
      datatype_for_mpi=MPI_REAL;
    } else if(finfo[fid].datainfo[did].datatype==FIO_REAL8){
      datatype_for_mpi=MPI_DOUBLE_PRECISION;
    } else if(finfo[fid].datainfo[did].datatype==FIO_INTEGER4){
      datatype_for_mpi=MPI_INTEGER;
    } else if(finfo[fid].datainfo[did].datatype==FIO_INTEGER8){
      datatype_for_mpi=MPI_INTEGER8;
    }
    count 
      = (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
      * (pow(2,finfo[fid].header.glevel-finfo[fid].header.rlevel)+2)
      *finfo[fid].datainfo[did].num_of_layer;

    if(finfo[fid].header.endiantype!=fio_system_endian){
      _data=malloc(count*finfo[fid].header.num_of_rgn*precision[finfo[fid].datainfo[did].datatype]);
      fio_ednchg(_data,
		precision[finfo[fid].datainfo[did].datatype],
		count*finfo[fid].header.num_of_rgn);
    } else {
      _data = data;
    }

    for(i=0;i<finfo[fid].header.num_of_rgn;i++){
      offset = offset0 + ( count*precision[finfo[fid].datainfo[did].datatype] )*finfo[fid].header.rgnid[i];
      data_start = count*i;
      dp = (char *)_data;
      dp=dp+data_start*precision[finfo[fid].datainfo[did].datatype];
      MPI_File_write_at(finfo[fid].status.mpi_fid, offset, 
			dp, count, datatype_for_mpi, &status);
    }
    if(finfo[fid].header.endiantype!=fio_system_endian){
      free(_data);
    }
    */
  } else {
#ifdef CONFIG_HDF5_INPUT
	  gid = H5Gopen(finfo[fid].status.file, "restart", H5P_DEFAULT);
	  for(did = 0; did < finfo[fid].header.num_of_data; did++) {
		  sprintf(dname, "DATA%2d", did+1);
		  dset = H5Dopen(gid, dname, H5P_DEFAULT);
		  /* varname */
		  attr = H5Aopen(dset, "varname", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, finfo[fid].dinfo[did].varname);
		  H5Aclose(attr);
		  /* description */
		  attr = H5Aopen(dset, "description", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, finfo[fid].dinfo[did].description);
		  H5Aclose(attr);
		  /* unit */
		  attr = H5Aopen(dset, "unit", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, finfo[fid].dinfo[did].unit);
		  H5Aclose(attr);
		  /* layername */
		  attr = H5Aopen(dset, "layername", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, finfo[fid].dinfo[did].layername);
		  H5Aclose(attr);
		  /* note */
		  attr = H5Aopen(dset, "note", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, finfo[fid].dinfo[did].note);
		  H5Aclose(attr);
		  /* datasize */
		  attr = H5Aopen(dset, "datasize", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, &temp64);
		  if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
		  finfo[fid].dinfo[did].datasize = temp64;
		  H5Aclose(attr);
		  /* datatype */
		  attr = H5Aopen(dset, "datatype", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, &temp32);
		  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
		  finfo[fid].dinfo[did].datatype = temp32;
		  H5Aclose(attr);
		  /* num_of_layer */
		  attr = H5Aopen(dset, "num_of_layer", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, &temp32);
		  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
		  finfo[fid].dinfo[did].num_of_layer = temp32;
		  H5Aclose(attr);
		  /* step */
		  attr = H5Aopen(dset, "step", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, &temp32);
		  if(system_ednchg){ fio_ednchg(&temp32,sizeof(int32_t),1 ); }
		  finfo[fid].dinfo[did].step = temp32;
		  H5Aclose(attr);
		  /* time_start */
		  attr = H5Aopen(dset, "time_start", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, &temp64);
		  if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
		  finfo[fid].dinfo[did].time_start = temp64;
		  H5Aclose(attr);
		  /* time_end */
		  attr = H5Aopen(dset, "time_end", H5P_DEFAULT);
		  type = H5Aget_type(attr);
		  H5Aread(attr, type, &temp64);
		  if(system_ednchg){ fio_ednchg(&temp64,sizeof(int64_t),1 ); }
		  finfo[fid].dinfo[did].time_start = temp64;
		  H5Aclose(attr);

		  H5Dclose(dset);
	  }
	  H5Gclose(gid);
#else
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
#endif
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

  if(common.use_mpiio) {

  } else {
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
  }
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

  ijkall = finfo[fid].header.rlevel * finfo[fid].header.rlevel * finfo[fid].dinfo[did].num_of_layer;

  datasize = ijkall * precision[finfo[fid].dinfo[did].datatype];

  if(common.use_mpiio) {

  } else {
    fseek(finfo[fid].status.fp,0L,SEEK_END);
    /* data */
    if(system_ednchg){ 
      _data = malloc(finfo[fid].dinfo[did].datasize);
      memcpy( _data, data, finfo[fid].dinfo[did].datasize);
      fio_ednchg(_data,precision[finfo[fid].dinfo[did].datatype],ijkall);
    }else{
      _data = data;
    }

    fwrite(_data,datasize,1,finfo[fid].status.fp);

    if(system_ednchg) { free(_data); }    
  }
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
#ifdef CONFIG_HDF5_INPUT
  hid_t gid, dset, type;
  char dname[FIO_HSHORT];
  char *_data;
  int x, y, z;
  int e_size;
#endif

  ijklall = finfo[fid].dinfo[did].datasize
          / precision[finfo[fid].dinfo[did].datatype];

  if(common.use_mpiio) {

  } else {
#ifdef CONFIG_HDF5_INPUT
	  e_size = precision[finfo[fid].dinfo[did].datatype];
	  _data = (char *)malloc(finfo[fid].dinfo[did].datasize);
	  gid = H5Gopen(finfo[fid].status.file, "restart", H5P_DEFAULT);
	  sprintf(dname, "DATA%2d", did+1);
	  dset = H5Dopen(gid, dname, H5P_DEFAULT);
	  type = H5Dget_type(dset);
	  H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, _data);
	  if(system_ednchg) {
		  fio_ednchg(_data, precision[finfo[fid].dinfo[did].datatype], ijklall);
	  }
	  for(y = 0; y <  common.rlevel_j; y++) {
		  for (x = 0; x < common.rlevel_i; x++) {
			  for(z = 0; z < finfo[fid].dinfo[did].num_of_layer; z++) {
				  memcpy((char *)data+(y*common.rlevel_i*finfo[fid].dinfo[did].num_of_layer+x*finfo[fid].dinfo[did].num_of_layer+z)*e_size,
					 _data+(z*common.rlevel_i*common.rlevel_j+x*common.rlevel_j+y)*e_size, e_size);
			  }
		  }
	  }
	  free(_data);
#else
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
#endif
  }
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

#ifndef CONFIG_HDF5
  fio_write_pkginfo( fid );
#endif

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

/** validate data size ************************************************/
int32_t fio_valid_datainfo( int32_t fid )
{
  int32_t did;
  int32_t ijall;
  int64_t datasize;
  char* str_dtype[5] = { "XXXX","REAL4","REAL8","INTEGER4","INTEGER8" };

  ijall = finfo[fid].header.rlevel;

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

inline double get_dtime(void){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return ((double)(tv.tv_sec) + (double)(tv.tv_usec) * 0.001 * 0.001);
}

#ifdef CONFIG_HDF5
static void fio_write_dinfo_attribute(int32_t fid, int32_t did, hid_t dset) {
	hid_t attr, attr_space;
	hid_t itype32, itype64;
	hsize_t dim;
	int64_t datasize;
	int32_t datatype;
	int32_t num_of_layer;
	int32_t step;
	int64_t time_start;
	int64_t time_end;
	int32_t tmp32;
	int64_t tmp64;

	if(common.endiantype == FIO_LITTLE_ENDIAN) {
		itype32 = H5T_STD_I32LE;
		itype64 = H5T_STD_I64LE;
	} else {
		itype32 = H5T_STD_I32BE;
		itype64 = H5T_STD_I64BE;
	}

	/* write varname */
	dim = FIO_HSHORT;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(dset, "varname", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].dinfo[did].varname);
	H5Aclose(attr);

	/* write description */
	dim = FIO_HMID;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(dset, "description", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].dinfo[did].description);
	H5Aclose(attr);

	/* write unit */
	dim = FIO_HSHORT;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(dset, "unit", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].dinfo[did].unit);
	H5Aclose(attr);

	/* write layername */
	dim = FIO_HSHORT;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(dset, "layername", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].dinfo[did].layername);
	H5Aclose(attr);

	/* write note */
	dim = FIO_HLONG;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(dset, "note", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].dinfo[did].note);
	H5Aclose(attr);
	
	/* write datasize */
	datasize = finfo[fid].dinfo[did].datasize;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(dset, "datasize", itype64, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&datasize,sizeof(int64_t),1);}
	H5Awrite(attr, itype64, &datasize);
	H5Aclose(attr);

	/* write datatype */
	datatype = finfo[fid].dinfo[did].datatype;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(dset, "datatype", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&datatype,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &datatype);
	H5Aclose(attr);

	/* write num_of_layer */
	num_of_layer = finfo[fid].dinfo[did].num_of_layer;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(dset, "num_of_layer", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&num_of_layer,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &num_of_layer);
	H5Aclose(attr);

	/* write step */
	step = finfo[fid].dinfo[did].step;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(dset, "step", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&step,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &step);
	H5Aclose(attr);

	/* write time_start */
	time_start = finfo[fid].dinfo[did].time_start;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(dset, "time_start", itype64, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&time_start,sizeof(int64_t),1);}
	H5Awrite(attr, itype64, &time_start);
	H5Aclose(attr);

	/* write time_end */
	time_end = finfo[fid].dinfo[did].time_end;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(dset, "time_end", itype64, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&time_end,sizeof(int64_t),1);}
	H5Awrite(attr, itype64, &time_end);
	H5Aclose(attr);
}

static void fio_write_header_attribute(int32_t fid) {
	hid_t file;
	hid_t attr, attr_space;
	hid_t itype32;
	hsize_t dim;
	int32_t tmp32;

	file = finfo[fid].status.file;

	if(common.endiantype == FIO_LITTLE_ENDIAN) {
		itype32 = H5T_STD_I32LE;
	} else {
		itype32 = H5T_STD_I32BE;
	}

	attr = H5Aopen(file, "num_of_data", H5P_DEFAULT);
	if( attr >= 0) {
		tmp32 = finfo[fid].header.num_of_data;
		attr_space = H5Screate(H5S_SCALAR);
		if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
		H5Awrite(attr, itype32, &tmp32);
		H5Aclose(attr);
		return;
	}

	/* write fname */
	dim = FIO_HLONG;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(file, "fname", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].header.fname);
	H5Aclose(attr);

	/* write description */
	dim = FIO_HMID;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(file, "description", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].header.description);
	H5Aclose(attr);

	/* write note */
	dim = FIO_HLONG;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(file, "note", H5T_C_S1, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr, H5T_C_S1, finfo[fid].header.note);
	H5Aclose(attr);

	/* write num_of_data */
	tmp32 = finfo[fid].header.num_of_data;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "num_of_data", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write fmode */
	tmp32 = finfo[fid].header.fmode;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "fmode", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write endiantype */
	tmp32 = finfo[fid].header.endiantype;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "endiantype", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write grid_topology */
	tmp32 = finfo[fid].header.grid_topology;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "grid_topology", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write glevel */
	tmp32 = finfo[fid].header.glevel;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "glevel", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write rlevel */
	tmp32 = finfo[fid].header.rlevel;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "rlevel", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write num_of_rgn */
	tmp32 = finfo[fid].header.num_of_rgn;
	attr_space = H5Screate(H5S_SCALAR);
	attr = H5Acreate(file, "num_of_rgn", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, &tmp32);
	H5Aclose(attr);

	/* write rgnid */
	dim = finfo[fid].header.num_of_rgn;
	attr_space = H5Screate_simple(1, &dim, NULL);
	attr = H5Acreate(file, "rgnid", itype32, attr_space, H5P_DEFAULT, H5P_DEFAULT);
	if(system_ednchg) {fio_ednchg(&tmp32,sizeof(int32_t),1);}
	H5Awrite(attr, itype32, finfo[fid].header.rgnid);
	H5Aclose(attr);
}
#endif

/** put & write data information and write data ***********************/
int32_t fio_put_write_datainfo_data( int32_t fid,
                                     datainfo_t ditem,
                                     void *data        )
{
  int32_t did;
#ifdef CONFIG_HDF5
  hid_t gid;
  hid_t space, dcpl, dset;
  int rank = 3; /* dimension */
  hid_t type;
  int64_t ijklall;
  void *_data_jik;
  char dname[FIO_HSHORT];
  hsize_t dims[3];
  hsize_t chunk[3];
  char *_data_kij;
  int e_size;
  int x, y, z;
#endif


  did = fio_new_datainfo( fid );

  fio_put_datainfo( fid, did, ditem );

#ifdef CONFIG_HDF5
  dims[0] = finfo[fid].dinfo[did].num_of_layer;
  dims[1] = common.rlevel_i;
  dims[2] = common.rlevel_j;
  chunk[0] = finfo[fid].dinfo[did].num_of_layer;
  chunk[1] = common.rlevel_i;
  chunk[2] = common.rlevel_j;

  _data_kij = malloc(finfo[fid].dinfo[did].datasize);
  e_size = precision[finfo[fid].dinfo[did].datatype];

  if(system_ednchg) {
	  ijklall = finfo[fid].dinfo[did].datasize
		  / precision[finfo[fid].dinfo[did].datatype];
	  _data_jik = malloc(finfo[fid].dinfo[did].datasize);
	  memcpy( _data_jik, data, finfo[fid].dinfo[did].datasize);
	  fio_ednchg(_data_jik,precision[finfo[fid].dinfo[did].datatype],ijklall);
  } else {
	  _data_jik = data;
  }

  if(common.endiantype == FIO_LITTLE_ENDIAN) {
	  if(finfo[fid].dinfo[did].datatype == FIO_REAL4)
		  type = H5T_IEEE_F32LE;
	  else if(finfo[fid].dinfo[did].datatype == FIO_REAL8)
		  type = H5T_IEEE_F64LE;
  } else if(common.endiantype == FIO_BIG_ENDIAN) {
	  if(finfo[fid].dinfo[did].datatype == FIO_REAL4)
		  type = H5T_IEEE_F32BE;
	  else if(finfo[fid].dinfo[did].datatype == FIO_REAL8)
		  type = H5T_IEEE_F64BE;
  }

  /* change index */
  for(z = 0; z <  finfo[fid].dinfo[did].num_of_layer; z++) {
	  for (x = 0; x < common.rlevel_i; x++) {
		  for(y = 0; y < common.rlevel_j; y++) {
			  memcpy((_data_kij+((z*common.rlevel_i*common.rlevel_j+x*common.rlevel_j+y)*e_size)),
				 (((char *)_data_jik)+((y*common.rlevel_i*finfo[fid].dinfo[did].num_of_layer+x*finfo[fid].dinfo[did].num_of_layer+z)*e_size)), e_size);
		  }
	  }
  }

  /* create dataspace */
  space = H5Screate_simple(rank, dims, NULL);

  /* set property and create dataset */
  dcpl = H5Pcreate(H5P_DATASET_CREATE);
#if 0
  H5Pset_deflate(dcpl, 9);
#endif
  H5Pset_szip(dcpl, H5_SZIP_NN_OPTION_MASK, 8);
  H5Pset_chunk(dcpl, rank, chunk);
  if(finfo[fid].tmp_gid == -1) { /* This is restart file. */
	  gid = H5Gopen1(finfo[fid].status.file, "restart");
	  if(gid < 0) {
		  gid = H5Gcreate(finfo[fid].status.file, "restart", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	  }
  } else {
	  gid = finfo[fid].tmp_gid;
  }
  sprintf(dname, "DATA%2d", finfo[fid].header.num_of_data);
  dset = H5Dcreate((hid_t)gid, dname,
		   type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);

  /* write header attribute*/
  fio_write_header_attribute(fid);

  /* create and write datainfo attribute*/
  fio_write_dinfo_attribute(fid, did, dset);

  /* write the data to dataset */
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, _data_kij);

  /* close property, dataset, dataspace */
  H5Pclose(dcpl);
  H5Dclose(dset);
  H5Sclose(space);
  if(finfo[fid].tmp_gid == -1) {
	  H5Gclose(gid);
  }
  if(system_ednchg)
	  free(_data_jik);
  free(_data_kij);
#else

  fio_write_pkginfo( fid ); /* update num_of_data */
  fio_write_datainfo( fid, did );
  fio_write_data( fid, did, data );
#endif

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

/** read pkginfo and datainfo *****************************************/
int32_t fio_read_allinfo_novalid( int32_t fid )
{
  int32_t i;

  fio_read_pkginfo( fid );
  /* fio_valid_pkginfo( fid ); */

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

/** read pkginfo and datainfo and get pkginfo *************************/
headerinfo_t fio_read_allinfo_get_pkginfo( int32_t fid )
{
  headerinfo_t hinfo;
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

  hinfo = fio_get_pkginfo( fid );

  return(hinfo);
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
  printf( " MODE ENDIAN GRID GL RL LALL MPI\n" );
  printf( " %4s", str_mode[common.fmode+1] );
  printf( " %6s", str_endian[common.endiantype+1] );
  printf( " %4s", str_topology[common.grid_topology+1] );
  printf( " %2d", common.glevel );
  printf( " %2d", common.rlevel );
  printf( " %4d", common.num_of_rgn );
  printf( " %3d", common.use_mpiio );
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
  char* str_topology[6] = { "XXXX","ICOSAHEDRON","LCP","MLCP","CARTESIAN","NONE" };
  char* str_dtype[5]    = { "XXXX","REAL4","REAL8","INTEGER4","INTEGER8" };

  /* exchange endian? */
  if ( endiantype!=system_endiantype ) {
    system_ednchg = 1;
    printf("%s to %s\n",str_endian[endiantype+1],str_endian[system_endiantype+1]);
  }
  common.use_mpiio = FIO_MPIIO_NOUSE;

  fio_fopen( fid, FIO_FREAD );
  fio_read_pkginfo ( fid );

  ijall = finfo[fid].header.rlevel;

  /* package info */
  printf( "============ DATA PACKAGE DEFINITION =============\n" );
  printf( "--- PACKAGE DESCRIPTION : %s\n", finfo[fid].header.description );
  printf( "--- COMPLETE/SPRIT      : %s\n", str_mode[finfo[fid].header.fmode+1] );
  printf( "--- GRID TOPOLOGY       : %s\n", str_topology[finfo[fid].header.grid_topology+1] );
  printf( "--- Grid resolution [m] : %d\n", finfo[fid].header.glevel );
  printf( "--- NUMBER OF XxY GRIDS : %d\n", finfo[fid].header.rlevel );
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

  /* memory allocation */
  if ( (finfo[fid].dinfo
       = (datainfo_t *)realloc(finfo[fid].dinfo,
         sizeof(datainfo_t)*finfo[fid].header.num_of_data)) == NULL ) {
    printf("Allocation error!\n");
  }
  fio_read_datainfo( fid );

  printf("============ ITEM ATTRIBUTES ======================\n");
  fsetpos(finfo[fid].status.fp, &(finfo[fid].status.eoh)); /* [add] 20111007 H.Yashiro */
  pos = dinfosize;                                         /* [mod] 20111007 H.Yashiro */
  for( did=0; did<finfo[fid].header.num_of_data; did++ ) {
    if ( dumptype == FIO_DUMP_ALL ) {
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
      printf("+++ [    ITEM, STEP,    Z,    XxY, REGION] VALUE\n");

      c = (char *)buf;
      for( l=0; l<finfo[fid].header.num_of_rgn; l++ ) {
      for( ij=0; ij<ijall; ij++ ) {
      for( k=0; k<finfo[fid].dinfo[did].num_of_layer; k++ ) {
        printf("+++ [%8s, %4d, %4d, %6d, %6d] ",
        finfo[fid].dinfo[did].varname,finfo[fid].dinfo[did].step,k,ij,l);

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

  }
  printf("============ END ITEM ATTRIBUTES ==================\n");
  fio_fclose( fid );

  return(SUCCESS_CODE);
}



