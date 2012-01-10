/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (CORE)                                   *
 *                                                                    *
 *  HISTORY                                                           *
 *    0.80      11-07-27  H.Tomita  : [NEW]                           *
 *    0.90      11-08-19  H.Yashiro : Incorporate into NICAM          *
 *    1.00      11-08-25  H.Yashiro : Complete format specification   *
 *                                                                    *
 **********************************************************************/
#include "fio.h"

/** endian change *****************************************************/
extern void fio_ednchg_( void* avp_pointer,
                         const int32_t *ai_size,
                         const int32_t *ai_num   );

/** filename generator ************************************************/
extern void fio_mk_fname_( char *fname,
                           const char *base,
                           const char *ext,
                           int32_t *i,
                           int32_t *y,
                           int32_t fname_len,
                           int32_t base_len,
                           int32_t ext_len    );

/** check system & initialze ******************************************/
extern void fio_syscheck_( void );

/** put common informtation *******************************************/
extern void fio_put_commoninfo_( int32_t *use_mpiio,
                                 int32_t *fmode,
                                 int32_t *endiantype,
                                 int32_t *grid_topology,
                                 int32_t *glevel,
                                 int32_t *rlevel,
                                 int32_t *num_of_rgn,
                                 int32_t *rgnid          );

/** put package information (full) ************************************/
extern void fio_put_pkginfo_( int32_t *fid,
                              headerinfo_t *hinfo );

/** get package information (full) ************************************/
extern void fio_get_pkginfo_( int32_t *fid, 
                              headerinfo_t *hinfo );

/** put data information (full) ***************************************/
extern void fio_put_datainfo_( int32_t *fid,
                               int32_t *did,
                               datainfo_t *ditem );

/** get data information (full) ***************************************/
extern void fio_get_datainfo_( int32_t *fid, 
                               int32_t *did, 
                               datainfo_t *ditem );

/** seek data id by varname and step **********************************/
extern void fio_seek_datainfo_( int32_t *did, 
                                int32_t *fid, 
                                char *varname,
                                int32_t *step,
                                int32_t varname_len );


/** open file IO stream ***********************************************/
extern void fio_fopen_( int32_t *fid, int32_t *mode );

/** close file IO stream **********************************************/
extern void  fio_fclose_( int32_t *fid );

/** write package information *****************************************/
extern void fio_write_pkginfo_( int32_t *fid );

/** read package information ******************************************/
extern void fio_read_pkginfo_( int32_t *fid );

/** write data information *****************************************/
extern void fio_write_datainfo_( int32_t *fid,
                                 int32_t *did  );

/** read data information ******************************************/
extern void fio_read_datainfo_( int32_t *fid );


/** write data array **************************************************/
extern void fio_write_data_( int32_t *fid,
                             int32_t *did,
                             void *data   );

/** read data array (full size) ***************************************/
extern void fio_read_data_( int32_t *fid,
                            int32_t *did,
                            void *data   );


/** register new file *************************************************/
extern void fio_register_file_( int32_t *fid,
                                char *fname,
                                int32_t fname_len );

/** put & write package information (quick put) ***********************/
extern void fio_put_write_pkginfo_( int32_t *fid,
                                    char *description,
                                    char *note,
                                    int32_t description_len,
                                    int32_t note_len         );


/** validate package information with common **************************/
extern void fio_valid_pkginfo_( int32_t *fid );


/** put & write data information and write data ***********************/
extern void fio_put_write_datainfo_data_( int32_t *did,
                                          int32_t *fid,
                                          datainfo_t *ditem,
                                          void *data         );

/** put & write data information **************************************/
extern void fio_put_write_datainfo_( int32_t *did,
                                     int32_t *fid,
                                     datainfo_t *ditem  );

/** read pkginfo and datainfo *****************************************/
extern void fio_read_allinfo_( int32_t *fid );

/** read pkginfo and datainfo *****************************************/
extern void fio_read_allinfo_novalid_( int32_t *fid );

/** read pkginfo and datainfo and get pkginfo *************************/
extern void fio_read_allinfo_get_pkginfo_( int32_t *fid,
                                           headerinfo_t *hinfo );

/** dump package summary of all finfo *********************************/
extern void fio_dump_finfolist_( void );

/** dump package summary of all finfo *********************************/
extern void fio_dump_finfo_( int32_t *fid,
                             int32_t *endiantype,
                             int32_t *dumptype    );
