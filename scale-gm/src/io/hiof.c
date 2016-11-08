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
#include "hiof.h"


/** endian change *****************************************************/
void hio_ednchg_( void* avp_pointer,
                  const int32_t *ai_size,
                  const int32_t *ai_num   )
{
  hio_ednchg( avp_pointer,
              *ai_size,
              *ai_num      );
}

/** filename generator ************************************************/
void hio_mk_fname_( char *fname,
                    const char *base,
                    const char *ext,
                    int32_t *i,
                    int32_t *y,
                    int32_t fname_len,
                    int32_t base_len,
                    int32_t ext_len    )
{
  char _base[HIO_HLONG];
  char _ext[HIO_HSHORT];
  int32_t _base_len = base_len;
  int32_t _ext_len  = ext_len;
  if( _base_len > HIO_HLONG-1  ){ _base_len = HIO_HLONG-1;  } /* [fix] H.Yashiro 20120621 */
  if( _ext_len  > HIO_HSHORT-1 ){ _ext_len  = HIO_HSHORT-1; } /* [fix] H.Yashiro 20120621 */

  hio_set_str( _base, base, _base_len );
  hio_set_str( _ext,  ext,  _ext_len  );

  hio_mk_fname( fname,
                _base,
                _ext,
                *i,
                *y      );
}

/** check system & initialze ******************************************/
void hio_syscheck_( void )
{
  int32_t ierr;

  ierr = hio_syscheck();
}

/** put common informtation *******************************************/
void hio_put_commoninfo_( int32_t *fmode,
                          int32_t *endiantype,
                          int32_t *grid_topology,
                          int32_t *glevel,
                          int32_t *rlevel,
                          int32_t *num_of_rgn,
                          int32_t *rgnid          )
{
  int32_t ierr;

  ierr = hio_put_commoninfo( *fmode,
                             *endiantype,
                             *grid_topology,
                             *glevel,
                             *rlevel,
                             *num_of_rgn,
                             rgnid           );
}

/** put common informtation from file *********************************/
void hio_put_commoninfo_fromfile_( int32_t *fid,
                                   int32_t *endiantype )
{
  int32_t ierr;

  ierr = hio_put_commoninfo_fromfile( *fid,
                                      *endiantype );
}

/** put package information (full) ************************************/
void hio_put_pkginfo_( int32_t *fid,
                       hio_headerinfo_t *hinfo )
{
  int32_t ierr;

  ierr = hio_put_pkginfo( *fid,
                          *hinfo );
}

/** get package information (full) ************************************/
void hio_get_pkginfo_( int32_t *fid, 
                       hio_headerinfo_t *hinfo )
{
  *hinfo = hio_get_pkginfo( *fid );
}

/** put data information (full) ***************************************/
void hio_put_datainfo_( int32_t *fid,
                        int32_t *did,
                        hio_datainfo_t *ditem )
{
  int32_t ierr;

  ierr = hio_put_datainfo( *fid,
                           *did,
                           *ditem );
}

/** get data information (full) ***************************************/
void hio_get_datainfo_( int32_t *fid, 
                        int32_t *did,
                        hio_datainfo_t *ditem )
{
  *ditem = hio_get_datainfo( *fid, *did );

}

/** get time information(ts,te) ***************************************/
void hio_get_timeinfo_( int32_t *fid,
                        int32_t *did,
                        int64_t *ts,
                        int64_t *te)
{
  hio_get_timeinfo( *fid, *did, ts, te );

#ifdef DEBUG
  fprintf(stdout,"dbg:hio_get_timeinfo_:fid=%d,did=%d\n",*fid,*did);
  for ( int n=0; n<2; n++ ){
    fprintf(stdout,"n=%d,ts[n]=%ld,te[n]=%ld\n",
            n, ts[n], te[n]);
  }
#endif
}



/** seek data id by varname and step **********************************/
void hio_seek_datainfo_( int32_t *did, 
                         int32_t *fid, 
                         char *varname,
                         int32_t *step,
                         int32_t varname_len )
{
  char _varname[HIO_HSHORT];
  int32_t _varname_len = varname_len;
  if( _varname_len > HIO_HSHORT-1 ){ _varname_len = HIO_HSHORT-1; } /* [fix] H.Yashiro 20120621 */

  hio_set_str( _varname, varname, _varname_len );

  *did = hio_seek_datainfo( *fid, 
                            _varname,
                            *step     );
}


/** open file IO stream ***********************************************/
void hio_fopen_( int32_t *fid, int32_t *mode )
{
  int32_t ierr;

  ierr = hio_fopen( *fid, *mode );
}

/** close file IO stream **********************************************/
void  hio_fclose_( int32_t *fid )
{
  int32_t ierr;

  ierr = hio_fclose( *fid );
}

/** write package information *****************************************/
void hio_write_pkginfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = hio_write_pkginfo( *fid );
}

/** read package information ******************************************/
void hio_read_pkginfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = hio_read_pkginfo( *fid );
}

/** write data information *****************************************/
void hio_write_datainfo_( int32_t *fid,
                          int32_t *did  )
{
  int32_t ierr;

  ierr = hio_write_datainfo( *fid,
                             *did  );
}

/** read data information *********************************************/
void hio_read_datainfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = hio_read_datainfo( *fid );
}

/** write data array **************************************************/
void hio_write_data_( int32_t *fid,
                      int32_t *did,
                        int32_t *step,
                        int64_t *ts,
                        int64_t *te,
                      void *data   )
{
  int32_t ierr;

  ierr = hio_write_data( *fid,
                         *did,
                        *step,
                        *ts,
                        *te,
                         data  );
}

/** write data array (1 region) ***************************************/
void hio_write_data_1rgn_( int32_t *fid,
                           int32_t *did,
                           int32_t *step,
                           int32_t *ll,
                           int64_t *ts,
                           int64_t *te,
                           void *data   )
{
  int32_t ierr;

  ierr = hio_write_data_1rgn( *fid,
                              *did,
                        *step,
                              *ll,
                        *ts,
                        *te,
                              data  );
}

/** read data array (full size) ***************************************/
void hio_read_data_( int32_t *fid,
                     int32_t *did,
                        int32_t *step,
                        int64_t *ts,
                        int64_t *te,
                     void *data   )
{
  int32_t ierr;

  ierr = hio_read_data( *fid,
                        *did,
                        *step,
                        ts,
                        te,
                        data  );
}

/** register new file *************************************************/
void hio_register_file_( int32_t *fid,
                         char *fname,
                         int32_t fname_len )
{
  char _fname[HIO_HLONG];
  int32_t _fname_len = fname_len;
  if( _fname_len > HIO_HLONG-1 ){ _fname_len = HIO_HLONG-1; } /* [fix] H.Yashiro 20120621 */

  hio_set_str( _fname, fname, _fname_len );

  *fid = hio_register_file( _fname );
}

/** put & write package information (quick put) ***********************/
void hio_put_write_pkginfo_( int32_t *fid,
                             char *description,
                             char *note,
                             int32_t description_len,
                             int32_t note_len         )
{
  int32_t ierr;

  char _description[HIO_HMID];
  char _note[HIO_HLONG];
  int32_t _description_len = description_len;
  int32_t _note_len        = note_len;
  if( _description_len > HIO_HMID-1  ){ _description_len = HIO_HMID-1;  } /* [fix] H.Yashiro 20120621 */
  if( _note_len        > HIO_HLONG-1 ){ _note_len        = HIO_HLONG-1; } /* [fix] H.Yashiro 20120621 */

  hio_set_str( _description, description, _description_len );
  hio_set_str( _note,        note,        _note_len        );

  ierr = hio_put_write_pkginfo( *fid,
                                _description,
                                _note         );
}

/** validate package information with common **************************/
void hio_valid_pkginfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = hio_valid_pkginfo( *fid );
}

/** validate package information with common (except rgnid) ***********/
void hio_valid_pkginfo_validrgn_( int32_t *fid,
                                  int32_t *rgnid )
{
  int32_t ierr;

  ierr = hio_valid_pkginfo_validrgn( *fid,
                                     rgnid );
}

/** put & write data information and write data ***********************/
void hio_put_write_datainfo_data_( int32_t *did,
                                   int32_t *fid,
                        int32_t *step,
                        int64_t *ts,
                        int64_t *te,
                                   hio_datainfo_t *ditem,
                                   void *data         )
{
  *did = hio_put_write_datainfo_data( *fid,
                                      *ditem,
                                      *step,
                                      *ts,
                                      *te,
                                      data    );
}

/** put & write data information **************************************/
void hio_put_write_datainfo_( int32_t *did,
                              int32_t *fid,
                              hio_datainfo_t *ditem  )
{
  *did = hio_put_write_datainfo( *fid,
                                 *ditem );
}

/** read pkginfo and datainfo *****************************************/
void hio_read_allinfo_( int32_t *fid )
{
  int32_t ierr;

  ierr = hio_read_allinfo( *fid );
}

/** read pkginfo and datainfo, with validating rgnid ******************/
void hio_read_allinfo_validrgn_( int32_t *fid,
                                 int32_t *rgnid )
{
  int32_t ierr;

  ierr = hio_read_allinfo_validrgn( *fid,
                                    rgnid );
}

/** allocate and copy datainfo ****************************************/
/* [add] C.Kodama 13-04-18 */
void hio_copy_datainfo_( int32_t *fid,
                         int32_t *fid_org )
{
  int32_t ierr;

  ierr = hio_copy_datainfo( *fid, *fid_org );
}

/** dump package summary of all finfo *********************************/
void hio_dump_finfolist_( void )
{
  int32_t ierr;

  ierr = hio_dump_finfolist();
}

/** dump package summary of all finfo *********************************/
void hio_dump_finfo_( int32_t *fid,
                      int32_t *endiantype,
                      int32_t *dumptype    )
{
  int32_t ierr;

  ierr = hio_dump_finfo( *fid,
                         *endiantype,
                         *dumptype    );
}


/** return num_of_var in file *****************************************/
void hio_get_num_of_var_( int32_t *fid,
                          int32_t *num_of_var )
{
  *num_of_var = hio_get_num_of_var( *fid );

}
