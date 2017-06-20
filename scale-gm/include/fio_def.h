/******************************************************************//**
 *                                                                    *
 *  Advanced file I/O module (defined vars)                           *
 *                                                                    *
 *  HISTORY                                                           *
 *    0.80      11-07-27  H.Tomita  : [NEW]                           *
 *    0.90      11-08-19  H.Yashiro : Incorporate into NICAM          *
 *    1.00      11-08-25  H.Yashiro : Complete format specification   *
 *    1.10      11-09-07  H.Yashiro : sepalate READ/APPEND mode       *
 *                                                                    *
 **********************************************************************/
#ifndef __FIO_DEF_H__
#define __FIO_DEF_H__

/* These parameters are also defined in mod_fio */

/* character length */
#define FIO_HSHORT  16
#define FIO_HMID    64
#define FIO_HLONG  256

/* data type */
#define FIO_REAL4     0
#define FIO_REAL8     1
#define FIO_INTEGER4  2
#define FIO_INTEGER8  3

/* data endian */
#define FIO_UNKNOWN_ENDIAN  0
#define FIO_LITTLE_ENDIAN   1
#define FIO_BIG_ENDIAN      2

/* topology */
#define FIO_ICOSAHEDRON  0
#define FIO_IGA_LCP      1
#define FIO_IGA_MLCP     2

/* file mode (partial or complete) */
#define FIO_SPLIT_FILE 0
#define FIO_INTEG_FILE 1


/* proccessor type */
#define FIO_SINGLE_PROC 0
#define FIO_MULTI_PROC  1


/* action type */
#define FIO_FREAD   0
#define FIO_FWRITE  1
#define FIO_FAPPEND 2 /* [add] H.Yashiro 20110907 overwrite mode */

/* data dump type */
#define FIO_DUMP_OFF      0
#define FIO_DUMP_HEADER   1
#define FIO_DUMP_ALL      2
#define FIO_DUMP_ALL_MORE 3

/* return type */
#define ERROR_CODE   -1
#define SUCCESS_CODE  1

#endif
