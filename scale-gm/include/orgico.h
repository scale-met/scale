#ifndef __ORGICO_H__
#define __ORGICO_H__

#include <stdint.h>

extern int32_t orgico_glevel;
extern int32_t orgico_rlevel;

extern char **orgico_fname_;
extern int32_t orgico_num_of_rgn;
extern int64_t orgico_block1R1L;

extern void orgico_setup( int32_t gl, int32_t rl );
extern void orgico_mk_fname( char *fname, char *base, int32_t i, int32_t y );
extern void orgico_readgriddata( char *fname, double *data, int i );
extern void orgico_setbasename( char *basename );
extern void orgico_readdata_seq( char *fname, int did, int nol, void *data, int32_t esize );
extern void orgico_readdata_dir( char *fname, int did, int nol, void *data, int32_t esize );


extern void orgico_setup_( int32_t *gl, int32_t *rl );
extern void orgico_mk_fname_( char *fname, char *base, int32_t *i, int32_t *y, int32_t fname_len, int32_t base_len );
extern void orgico_readgriddata_( char *fname, double *data, int32_t *i, int32_t fname_len );
extern void orgico_setbasename_( char *basename, int32_t basename_len );
extern void orgico_readdata_seq_( char *fname, int *did, int *nol, void *data, int32_t *esize, int32_t fname_len );
extern void orgico_readdata_dir_( char *fname, int *did, int *nol, void *data, int32_t *esize, int32_t fname_len );



#endif
