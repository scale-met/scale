/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (with HDF5)                              *
 *                                                                    *
 *  HISTORY                                                           *
 *    1.24      15-11-04  T.Inoue  : [NEW] based on fio.c.             *
 *                                                                    *
 **********************************************************************/
#ifndef __HIO_H__
#define __HIO_H__

#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "poh5.h"

#include "hio_def.h"

typedef float (real32_t);
typedef double(real64_t);

/* common information */
typedef struct{
  int32_t fmode;
  int32_t endiantype;
  int32_t grid_topology;
  int32_t glevel;
  int32_t rlevel;
  int32_t num_of_rgn;
  int32_t *rgnid;
} hio_commoninfo_t;

/* package information */
typedef struct{
  char fname[HIO_HLONG];
  char description[HIO_HMID];
  char note[HIO_HLONG];
  int32_t num_of_var;
  int32_t fmode;
  int32_t endiantype;
  int32_t grid_topology;
  int32_t glevel;
  int32_t rlevel;
  int32_t num_of_rgn;
  int32_t *rgnid;
} hio_headerinfo_t;

/* data item information */
typedef struct{
  char varname[HIO_HSHORT];
  char description[HIO_HMID];
  char unit[HIO_HSHORT];
  char layername[HIO_HSHORT];
  char note[HIO_HLONG];
  int64_t datasize;
  int32_t datatype;
  int32_t num_of_layer;
  int32_t num_of_step;
  /* int64_t *ts; */
  /* int64_t *te; */
} hio_datainfo_t;

/* status information */
typedef struct{
  int32_t rwmode;
  int32_t opened;
  /* FILE *fp; */
  phid_t hfid; /**< file_id for poh5 file */
  fpos_t eoh;
} hio_statusinfo_t;

/* file structure */
typedef struct{
  hio_headerinfo_t header;
  hio_datainfo_t *dinfo;
  hio_statusinfo_t status;
} hio_fileinfo_t;


/** endian change *****************************************************/
extern void hio_ednchg( void* avp_pointer,
                        const int32_t ai_size,
                        const int32_t ai_num   );

/** filename generator ************************************************/
extern void hio_mk_fname( char *fname,
                          const char *base,
                          const char *ext,
                          int32_t i,
                          int32_t y         );

/** string preprocess *************************************************/
extern void hio_set_str( char *_str,
                         const char *str,
                         const int str_len );

/** string postprocess *************************************************/
extern void hio_trim_str( char *_str,
                          const char *str,
                          const int str_len );

/** check system & initialze ******************************************/
extern int32_t hio_syscheck( void );

/** store common informtation *****************************************/
extern int32_t hio_put_commoninfo( int32_t fmode,
                                   int32_t endiantype,
                                   int32_t grid_topology,
                                   int32_t glevel,
                                   int32_t rlevel,
                                   int32_t num_of_rgn,
                                   int32_t rgnid[]        );

/** put common informtation from file *********************************/
extern int32_t hio_put_commoninfo_fromfile( int32_t fid,
                                            int32_t endiantype );

/** add new file structure ********************************************/
static int32_t hio_new_finfo( void );                     /*<internal>*/

/** add new file structure ********************************************/
static int32_t hio_new_datainfo( int32_t fid );           /*<internal>*/

/** put package information (full) ************************************/
extern int32_t hio_put_pkginfo( int32_t fid,
                                hio_headerinfo_t hinfo );

/** get package information (full) ************************************/
extern hio_headerinfo_t hio_get_pkginfo( int32_t fid );

/** put data information (full) ***************************************/
extern int32_t hio_put_datainfo( int32_t fid,
                                 int32_t did,
                                 hio_datainfo_t ditem );

/** get data information (full) ***************************************/
extern hio_datainfo_t hio_get_datainfo( int32_t fid,
                                    int32_t did  );

/** get time information(ts,te) ***************************************/
void hio_get_timeinfo( int32_t fid,
                       int32_t did,
                       int64_t ts[],
                       int64_t te[]);

/** seek data id by varname and step **********************************/
extern int32_t hio_seek_datainfo( int32_t fid,
                                  char *varname,
                                  int32_t step     );

/** return num_of_var in file *****************************************/
int32_t hio_get_num_of_var( int32_t fid );

/** open file IO stream ***********************************************/
extern int32_t hio_fopen( int32_t fid, int32_t mode );

/** close file IO stream **********************************************/
extern int32_t hio_fclose( int32_t fid );

/** write package information *****************************************/
extern int32_t hio_write_pkginfo( int32_t fid );

/** read package information ******************************************/
extern int32_t hio_read_pkginfo( int32_t fid );

/** write data information ********************************************/
extern int32_t hio_write_datainfo( int32_t fid,
                                   int32_t did  );

/** read data information *********************************************/
extern int32_t hio_read_datainfo( int32_t fid );

/** write data array **************************************************/
extern int32_t hio_write_data( int32_t fid,
                               int32_t did,
                        int32_t step,
                        int64_t ts,
                        int64_t te,
                               void *data   );

/** write data array (1 region) ***************************************/
extern int32_t hio_write_data_1rgn( int32_t fid,
                                    int32_t did,
                        int32_t step,
                                    int32_t ll,
                        int64_t ts,
                        int64_t te,
                                    void *data   );

/** read data array (full size) ***************************************/
extern int32_t hio_read_data( int32_t fid,
                              int32_t did,
                       const int32_t step,
                       int64_t *ts,
                       int64_t *te,  
                              void *data   );

/** register new file *************************************************/
extern int32_t hio_register_file( char *fname );

/** put & write package information (quick put) ***********************/
extern int32_t hio_put_write_pkginfo( int32_t fid,
                                      char *description,
                                      char *note         );

/** validate package information with common **************************/
extern int32_t hio_valid_pkginfo( int32_t fid );

/** validate package information with common (except rgnid) ***********/
extern int32_t hio_valid_pkginfo_validrgn( int32_t fid,
                                           int32_t rgnid[] );

/** put & write data information and write data ***********************/
extern int32_t hio_put_write_datainfo_data( int32_t fid,
                                            hio_datainfo_t ditem,
                                     int32_t step,  /* \todo tentative. How to set these ?? */
                                     int64_t ts,
                                     int64_t te,
                                            void *data        );

/** put & write data information **************************************/
extern int32_t hio_put_write_datainfo( int32_t fid,
                                       hio_datainfo_t ditem );

/** read pkginfo and datainfo *****************************************/
extern int32_t hio_read_allinfo( int32_t fid );

/** read pkginfo and datainfo, with validating rgnid ******************/
extern int32_t hio_read_allinfo_validrgn( int32_t fid,
                                          int32_t rgnid[] );

/** allocate and copy datainfo ****************************************/
extern int32_t hio_copy_datainfo( int32_t fid,
                                  int32_t fid_org );

/** dump package summary of all finfo *********************************/
extern int32_t hio_dump_finfolist( void );

/** dump package summary of all finfo *********************************/
extern int32_t hio_dump_finfo( int32_t fid,
                               int32_t endiantype,
                               int32_t dumptype    );


#endif
