/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (CORE)                                   *
 *                                                                    *
 *  HISTORY                                                           *
 *    0.80      11-07-27  H.Tomita  : [NEW]                           *
 *    0.90      11-08-19  H.Yashiro : Incorporate into NICAM          *
 *    1.00      11-08-25  H.Yashiro : Complete format specification   *
 *    1.23      13-04-18  C.Kodama  : [add] fio_read_datainfo_tmpdata,*
 *                                      fio_register_vname_tmpdata,   *
 *                                      fio_read_data_tmpdata,        *
 *                                      fio_read_allinfo_tmpdata,     *
 *                                      fio_copy_datainfo             *
 *                                                                    *
 **********************************************************************/
#include "hio.h"

/** endian change *****************************************************/
extern void hio_ednchg_( void* avp_pointer,
                         const int32_t *ai_size,
                         const int32_t *ai_num   );

/** filename generator ************************************************/
extern void hio_mk_fname_( char *fname,
                           const char *base,
                           const char *ext,
                           int32_t *i,
                           int32_t *y,
                           int32_t fname_len,
                           int32_t base_len,
                           int32_t ext_len    );

/** check system & initialze ******************************************/
extern void hio_syscheck_( void );

/** put common informtation *******************************************/
extern void hio_put_commoninfo_( int32_t *fmode,
                                 int32_t *endiantype,
                                 int32_t *grid_topology,
                                 int32_t *glevel,
                                 int32_t *rlevel,
                                 int32_t *num_of_rgn,
                                 int32_t *rgnid          );

/** put common informtation from file *********************************/
extern void hio_put_commoninfo_fromfile_( int32_t *fid,
                                          int32_t *endiantype );

/** put package information (full) ************************************/
extern void hio_put_pkginfo_( int32_t *fid,
                              hio_headerinfo_t *hinfo );

/** get package information (full) ************************************/
extern void hio_get_pkginfo_( int32_t *fid,
                              hio_headerinfo_t *hinfo );

/** put data information (full) ***************************************/
extern void hio_put_datainfo_( int32_t *fid,
                               int32_t *did,
                               hio_datainfo_t *ditem );

/** get data information (full) ***************************************/
extern void hio_get_datainfo_( int32_t *fid,
                               int32_t *did,
                               hio_datainfo_t *ditem );

/** get time information(ts,te) ***************************************/
void hio_get_timeinfo_( int32_t *fid,
                        int32_t *did,
                        int64_t *ts,
                        int64_t *te);

/** seek data id by varname and step **********************************/
extern void hio_seek_datainfo_( int32_t *did,
                                int32_t *fid,
                                char *varname,
                                int32_t *step,
                                int32_t varname_len );


/** open file IO stream ***********************************************/
extern void hio_fopen_( int32_t *fid, int32_t *mode );

/** close file IO stream **********************************************/
extern void  hio_fclose_( int32_t *fid );

/** write package information *****************************************/
extern void hio_write_pkginfo_( int32_t *fid );

/** read package information ******************************************/
extern void hio_read_pkginfo_( int32_t *fid );

/** write data information *****************************************/
extern void hio_write_datainfo_( int32_t *fid,
                                 int32_t *did  );

/** read data information ******************************************/
extern void hio_read_datainfo_( int32_t *fid );

/** write data array **************************************************/
extern void hio_write_data_( int32_t *fid,
                             int32_t *did,
                        int32_t *step,
                        int64_t *ts,
                        int64_t *te,
                             void *data   );

/** read data array (full size) ***************************************/
extern void hio_read_data_( int32_t *fid,
                            int32_t *did,
                        int32_t *step,
                        int64_t *ts,
                        int64_t *te,
                            void *data   );

/** register new file *************************************************/
extern void hio_register_file_( int32_t *fid,
                                char *fname,
                                int32_t fname_len );

/** put & write package information (quick put) ***********************/
extern void hio_put_write_pkginfo_( int32_t *fid,
                                    char *description,
                                    char *note,
                                    int32_t description_len,
                                    int32_t note_len         );


/** validate package information with common **************************/
extern void hio_valid_pkginfo_( int32_t *fid );

/** validate package information with common (except rgnid) ***********/
extern void hio_valid_pkginfo_validrgn_( int32_t *fid,
                                         int32_t *rgnid );

/** put & write data information and write data ***********************/
extern void hio_put_write_datainfo_data_( int32_t *did,
                                          int32_t *fid,
                        int32_t *step,
                        int64_t *ts,
                        int64_t *te,
                                          hio_datainfo_t *ditem,
                                          void *data         );

/** put & write data information **************************************/
extern void hio_put_write_datainfo_( int32_t *did,
                                     int32_t *fid,
                                     hio_datainfo_t *ditem  );

/** read pkginfo and datainfo *****************************************/
extern void hio_read_allinfo_( int32_t *fid );

/** read pkginfo and datainfo, with validating rgnid ******************/
extern void hio_read_allinfo_validrgn_( int32_t *fid,
                                        int32_t *rgnid );

/** allocate and copy datainfo ****************************************/
/* [add] C.Kodama 13-04-18 */
extern void hio_copy_datainfo_( int32_t *fid,
                                int32_t *fid_org );

/** dump package summary of all finfo *********************************/
extern void hio_dump_finfolist_( void );

/** dump package summary of all finfo *********************************/
extern void hio_dump_finfo_( int32_t *fid,
                             int32_t *endiantype,
                             int32_t *dumptype    );

/** return num_of_var in file *****************************************/
void hio_get_num_of_var_( int32_t *fid,
                          int32_t *num_of_var );
