/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (CORE)                                   *
 *                                                                    *
 *  HISTORY                                                           *
 *    0.80      11-07-27  H.Tomita  : [NEW]                           *
 *    0.90      11-08-19  H.Yashiro : Incorporate into NICAM          *
 *    1.00      11-08-25  H.Yashiro : Complete format specification   *
 *    1.22      12-06-21  H.Yashiro : [fix] bad char tail treatment   *
 *                                    thanks to ohno-san@riken        *
 *    1.23      13-04-18  C.Kodama  : [add] fio_read_datainfo_tmpdata,*
 *                                      fio_register_vname_tmpdata,   *
 *                                      fio_read_data_tmpdata,        *
 *                                      fio_read_allinfo_tmpdata,     *
 *                                      fio_copy_datainfo             *
 *                                                                    *
 **********************************************************************/
#include "fiof.h"


/** endian change *****************************************************/
void fio_ednchg_( void* avp_pointer,
                  const int32_t *ai_size,
                  const int32_t *ai_num   )
{
  fio_ednchg( avp_pointer,
              *ai_size,
              *ai_num      );
}

/** filename generator ************************************************/
void fio_mk_fname_( char *fname,
                    const char *base,
                    const char *ext,
                    int32_t *i,
                    int32_t *y,
                    int32_t fname_len,
                    int32_t base_len,
                    int32_t ext_len    )
{
  char _base[FIO_HLONG];
  char _ext[FIO_HSHORT];
  int32_t _base_len = base_len;
  int32_t _ext_len  = ext_len;
  if( _base_len > FIO_HLONG-1  ){ _base_len = FIO_HLONG-1;  } /* [fix] H.Yashiro 20120621 */
  if( _ext_len  > FIO_HSHORT-1 ){ _ext_len  = FIO_HSHORT-1; } /* [fix] H.Yashiro 20120621 */

  fio_set_str( _base, base, _base_len );
  fio_set_str( _ext,  ext,  _ext_len  );

  fio_mk_fname( fname,
                _base,
                _ext,
                *i,
                *y      );
}

/** check system & initialze ******************************************/
void fio_syscheck_( void )
{
  int32_t ierr;

  ierr = fio_syscheck();
}

/** put common informtation *******************************************/
void fio_put_commoninfo_( int32_t *fmode,
                          int32_t *endiantype,
                          int32_t *grid_topology,
                          int32_t *glevel,
                          int32_t *rlevel,
                          int32_t *num_of_rgn,
                          int32_t *rgnid          )
{
  int32_t ierr;

  ierr = fio_put_commoninfo( *fmode,
                             *endiantype,
                             *grid_topology,
                             *glevel,
                             *rlevel,
                             *num_of_rgn,
                             rgnid           );
}

/** put common informtation from file *********************************/
void fio_put_commoninfo_fromfile_( int32_t *fid,
                                   int32_t *endiantype )
{
  int32_t ierr;

  ierr = fio_put_commoninfo_fromfile( *fid,
                                      *endiantype );
}

/** put package information (full) ************************************/
void fio_put_pkginfo_( int32_t *fid,
                       headerinfo_t *hinfo )
{
  int32_t ierr;

  ierr = fio_put_pkginfo( *fid,
                          *hinfo );
}

/** get package information (full) ************************************/
void fio_get_pkginfo_( int32_t *fid, 
                       headerinfo_t *hinfo )
{
  *hinfo = fio_get_pkginfo( *fid );
}

/** put data information (full) ***************************************/
void fio_put_datainfo_( int32_t *fid,
                        int32_t *did,
                        datainfo_t *ditem )
{
  int32_t ierr;

  ierr = fio_put_datainfo( *fid,
                           *did,
                           *ditem );
}

/** get data information (full) ***************************************/
void fio_get_datainfo_( int32_t *fid, 
                        int32_t *did,
                        datainfo_t *ditem )
{
  *ditem = fio_get_datainfo( *fid,
                             *did );
}

/** seek data id by varname and step **********************************/
void fio_seek_datainfo_( int32_t *did, 
                         int32_t *fid, 
                         char *varname,
                         int32_t *step,
                         int32_t varname_len )
{
  char _varname[FIO_HSHORT];
  int32_t _varname_len = varname_len;
  if( _varname_len > FIO_HSHORT-1 ){ _varname_len = FIO_HSHORT-1; } /* [fix] H.Yashiro 20120621 */

  fio_set_str( _varname, varname, _varname_len );

  *did = fio_seek_datainfo( *fid, 
                            _varname,
                            *step     );
}


/** open file IO stream ***********************************************/
void fio_fopen_( int32_t *fid, int32_t *mode )
{
  int32_t ierr;

  ierr = fio_fopen( *fid, *mode );
}

/** close file IO stream **********************************************/
void  fio_fclose_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_fclose( *fid );
}

/** write package information *****************************************/
void fio_write_pkginfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_write_pkginfo( *fid );
}

/** read package information ******************************************/
void fio_read_pkginfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_read_pkginfo( *fid );
}

/** write data information *****************************************/
void fio_write_datainfo_( int32_t *fid,
                          int32_t *did  )
{
  int32_t ierr;

  ierr = fio_write_datainfo( *fid,
                             *did  );
}

/** read data information *********************************************/
void fio_read_datainfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_read_datainfo( *fid );
}

/** read data information and store data as tmpdata *******************/
/* [add] C.Kodama 13-04-18 */
void fio_read_datainfo_tmpdata_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_read_datainfo_tmpdata( *fid );
}


/** write data array **************************************************/
void fio_write_data_( int32_t *fid,
                      int32_t *did,
                      void *data   )
{
  int32_t ierr;

  ierr = fio_write_data( *fid,
                         *did,
                         data  );
}

/** write data array (1 region) ***************************************/
void fio_write_data_1rgn_( int32_t *fid,
                           int32_t *did,
                           void *data   )
{
  int32_t ierr;

  ierr = fio_write_data_1rgn( *fid,
                              *did,
                              data  );
}

/** read data array (full size) ***************************************/
void fio_read_data_( int32_t *fid,
                     int32_t *did,
                     void *data   )
{
  int32_t ierr;

  ierr = fio_read_data( *fid,
                        *did,
                        data  );
}

/* [add] C.Kodama 13-04-18 */
void fio_read_data_tmpdata_( int32_t *fid,
                             int32_t *did,
                             void *data   )
{
  int32_t ierr;

  ierr = fio_read_data_tmpdata( *fid,
                                *did,
                                data  );
}


/** register new file *************************************************/
void fio_register_file_( int32_t *fid,
                         char *fname,
                         int32_t fname_len )
{
  char _fname[FIO_HLONG];
  int32_t _fname_len = fname_len;
  if( _fname_len > FIO_HLONG-1 ){ _fname_len = FIO_HLONG-1; } /* [fix] H.Yashiro 20120621 */

  fio_set_str( _fname, fname, _fname_len );

  *fid = fio_register_file( _fname );
}

/** put & write package information (quick put) ***********************/
void fio_put_write_pkginfo_( int32_t *fid,
                             char *description,
                             char *note,
                             int32_t description_len,
                             int32_t note_len         )
{
  int32_t ierr;

  char _description[FIO_HMID];
  char _note[FIO_HLONG];
  int32_t _description_len = description_len;
  int32_t _note_len        = note_len;
  if( _description_len > FIO_HMID-1  ){ _description_len = FIO_HMID-1;  } /* [fix] H.Yashiro 20120621 */
  if( _note_len        > FIO_HLONG-1 ){ _note_len        = FIO_HLONG-1; } /* [fix] H.Yashiro 20120621 */

  fio_set_str( _description, description, _description_len );
  fio_set_str( _note,        note,        _note_len        );

  ierr = fio_put_write_pkginfo( *fid,
                                _description,
                                _note         );
}

/** validate package information with common **************************/
void fio_valid_pkginfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_valid_pkginfo( *fid );
}

/** validate package information with common (except rgnid) ***********/
void fio_valid_pkginfo_validrgn_( int32_t *fid,
                                  int32_t *rgnid )
{
  int32_t ierr;

  ierr = fio_valid_pkginfo_validrgn( *fid,
                                     rgnid );
}

/** put & write data information and write data ***********************/
void fio_put_write_datainfo_data_( int32_t *did,
                                   int32_t *fid,
                                   datainfo_t *ditem,
                                   void *data         )
{
  *did = fio_put_write_datainfo_data( *fid,
                                      *ditem,
                                      data    );
}

/** put & write data information **************************************/
void fio_put_write_datainfo_( int32_t *did,
                              int32_t *fid,
                              datainfo_t *ditem  )
{
  *did = fio_put_write_datainfo( *fid,
                                 *ditem );
}

/** read pkginfo and datainfo *****************************************/
void fio_read_allinfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_read_allinfo( *fid );
}

/** read pkginfo and datainfo, with validating rgnid ******************/
void fio_read_allinfo_validrgn_( int32_t *fid,
                                 int32_t *rgnid )
{
  int32_t ierr;

  ierr = fio_read_allinfo_validrgn( *fid,
                                    rgnid );
}

/** read pkginfo and datainfo and store data as tmpdata ***************/
/* [add] C.Kodama 13-04-18 */
void fio_read_allinfo_tmpdata_( int32_t *fid )
{
  int32_t ierr;

  ierr = fio_read_allinfo_tmpdata( *fid );
}

/* [add] C.Kodama 13-04-18 */
void fio_register_vname_tmpdata_( const char *vname_in, 
                                  int32_t    *vname_len )
{
  char _vname_in[FIO_HSHORT];
  int32_t _vname_len  = *vname_len;
  if( _vname_len  > FIO_HSHORT-1 ){ _vname_len  = FIO_HSHORT-1; }
  fio_set_str( _vname_in,  vname_in,  _vname_len  );
  fio_register_vname_tmpdata( _vname_in );
}

/** allocate and copy datainfo ****************************************/
/* [add] C.Kodama 13-04-18 */
void fio_copy_datainfo_( int32_t *fid,
                         int32_t *fid_org )
{
  int32_t ierr;

  ierr = fio_copy_datainfo( *fid, *fid_org );
}

/** dump package summary of all finfo *********************************/
void fio_dump_finfolist_( void )
{
  int32_t ierr;

  ierr = fio_dump_finfolist();
}

/** dump package summary of all finfo *********************************/
void fio_dump_finfo_( int32_t *fid,
                      int32_t *endiantype,
                      int32_t *dumptype    )
{
  int32_t ierr;

  ierr = fio_dump_finfo( *fid,
                         *endiantype,
                         *dumptype    );
}


