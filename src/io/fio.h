/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (CORE)                                   *
 *                                                                    *
 *  HISTORY                                                           *
 *    0.80      11-07-27  H.Tomita  : [NEW]                           *
 *    0.90      11-08-19  H.Yashiro : Incorporate into NICAM          *
 *    1.00      11-08-25  H.Yashiro : Complete format specification   *
 *    1.20      11-10-07  H.Yashiro : sepalate MPI/nonMPI fpos        *
 *                                    thanks to kameyama-san@riken    *
 *                                                                    *
 **********************************************************************/
#ifndef __FIO_H__
#define __FIO_H__

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef NO_MPIIO
#include <mpi.h>
#endif

#include "fio_def.h"

typedef float (real32_t);
typedef double(real64_t);

/* common information */
typedef struct{
  int32_t use_mpiio;

  int32_t fmode;
  int32_t endiantype;
  int32_t grid_topology;
  int32_t glevel;
  int32_t rlevel;
  int32_t num_of_rgn;
  int32_t *rgnid;
} commoninfo_t;

/* package information */
typedef struct{
  char fname[FIO_HLONG];
  char description[FIO_HMID];
  char note[FIO_HLONG];
  int32_t  num_of_data;

  int32_t fmode;
  int32_t endiantype;
  int32_t grid_topology;
  int32_t glevel;
  int32_t rlevel;
  int32_t num_of_rgn;
  int32_t *rgnid;
} headerinfo_t;

/* data item information */
typedef struct{
  char varname[FIO_HSHORT];
  char description[FIO_HMID];
  char unit[FIO_HSHORT];
  char layername[FIO_HSHORT];
  char note[FIO_HLONG];
  int64_t datasize;
  int32_t datatype;
  int32_t num_of_layer;
  int32_t step;
  int64_t time_start;
  int64_t time_end;
} datainfo_t; 

/* status information */
typedef struct{
  int32_t rwmode;
  int32_t opened;
  FILE *fp;
  fpos_t eoh;
#ifndef NO_MPIIO
  MPI_Offset mpi_eoh;
  MPI_File mpi_fid;
#endif
} statusinfo_t;

/* file structure */
typedef struct{
  headerinfo_t header;
  datainfo_t *dinfo;
  statusinfo_t status;
} fileinfo_t;

extern commoninfo_t common;
extern fileinfo_t *finfo;
extern headerinfo_t *hinfo;
extern datainfo_t *ditem;


/** endian change *****************************************************/
extern void fio_ednchg( void* avp_pointer,
                        const int32_t ai_size,
                        const int32_t ai_num   );

/** filename generator ************************************************/
extern void fio_mk_fname( char *fname,
                          const char *base,
                          const char *ext,
                          int32_t i,
                          int32_t y         );

/** string preprocess *************************************************/
extern void fio_set_str( char *_str,
                         const char *str,
                         const int str_len );

/** string postprocess *************************************************/
extern void fio_trim_str( char *_str,
                         const char *str,
                         const int str_len );

/** check system & initialze ******************************************/
extern int32_t fio_syscheck( void );

/** store common informtation *****************************************/
extern int32_t fio_put_commoninfo( int32_t use_mpiio,
                                   int32_t fmode,
                                   int32_t endiantype,
                                   int32_t grid_topology,
                                   int32_t glevel,
                                   int32_t rlevel,
                                   int32_t num_of_rgn,
                                   int32_t rgnid[]        );

/** add new file structure ********************************************/
static int32_t fio_new_finfo( void );                     /*<internal>*/

/** add new file structure ********************************************/
static int32_t fio_new_datainfo( int32_t fid );           /*<internal>*/

/** put package information (full) ************************************/
extern int32_t fio_put_pkginfo( int32_t fid,
                                headerinfo_t hinfo );

/** get package information (full) ************************************/
extern headerinfo_t fio_get_pkginfo( int32_t fid );

/** put data information (full) ***************************************/
extern int32_t fio_put_datainfo( int32_t fid,
                                 int32_t did,
                                 datainfo_t ditem );

/** get data information (full) ***************************************/
extern datainfo_t fio_get_datainfo( int32_t fid, 
                                    int32_t did  );

/** seek data id by varname and step **********************************/
extern int32_t fio_seek_datainfo( int32_t fid, 
                                  char *varname,
                                  int32_t step   );


/** open file IO stream ***********************************************/
extern int32_t fio_fopen( int32_t fid, int32_t mode );

/** close file IO stream **********************************************/
extern int32_t fio_fclose( int32_t fid );

/** write package information *****************************************/
extern int32_t fio_write_pkginfo( int32_t fid );

/** read package information ******************************************/
extern int32_t fio_read_pkginfo( int32_t fid );

/** write data information ********************************************/
extern int32_t fio_write_datainfo( int32_t fid,
                                   int32_t did  );

/** read data information *********************************************/
extern int32_t fio_read_datainfo( int32_t fid );

/** write data array **************************************************/
extern int32_t fio_write_data( int32_t fid,
                               int32_t did,
                               void *data   );

/** write data array (1 region) ***************************************/
extern int32_t fio_write_data_1rgn( int32_t fid,
                                    int32_t did,
                                    void *data   );

/** read data array (full size) ***************************************/
extern int32_t fio_read_data( int32_t fid,
                              int32_t did,
                              void *data   );


/** register new file *************************************************/
extern int32_t fio_register_file( char *fname );

/** put & write package information (quick put) ***********************/
extern int32_t fio_put_write_pkginfo( int32_t fid,
                                      char *description,
                                      char *note         );

/** validate package information with common **************************/
extern int32_t fio_valid_pkginfo( int32_t fid );

/** put & write data information and write data ***********************/
extern int32_t fio_put_write_datainfo_data( int32_t fid,
                                            datainfo_t ditem,
                                            void *data        );

/** put & write data information **************************************/
extern int32_t fio_put_write_datainfo( int32_t fid,
                                       datainfo_t ditem );

/** read pkginfo and datainfo and get pkginfo *************************/
extern int32_t fio_read_allinfo( int32_t fid );

/** read pkginfo and datainfo and get pkginfo *************************/
extern headerinfo_t fio_read_allinfo_get_pkginfo( int32_t fid );

/** dump package summary of all finfo *********************************/
extern int32_t fio_dump_finfolist( void );

/** dump package summary of all finfo *********************************/
extern int32_t fio_dump_finfo( int32_t fid,
                               int32_t endiantype,
                               int32_t dumptype    );


#endif
