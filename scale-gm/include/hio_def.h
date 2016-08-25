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
#ifndef __HIO_DEF_H__
#define __HIO_DEF_H__

/* These parameters are also defined in mod_fio */

/* character length */
#define HIO_HSHORT  16
#define HIO_HMID    64
#define HIO_HLONG  256

/* data type */
#define HIO_REAL4     0
#define HIO_REAL8     1
#define HIO_INTEGER4  2
#define HIO_INTEGER8  3

/* data endian */
#define HIO_UNKNOWN_ENDIAN  0
#define HIO_LITTLE_ENDIAN   1
#define HIO_BIG_ENDIAN      2

/* topology */
#define HIO_ICOSAHEDRON  0
#define HIO_IGA_LCP      1
#define HIO_IGA_MLCP     2

/* file mode (partial or complete) */
#define HIO_SPLIT_FILE 0
#define HIO_INTEG_FILE 1


/* proccessor type */
#define HIO_SINGLE_PROC 0
#define HIO_MULTI_PROC  1


/* action type */
#define HIO_FREAD   0
#define HIO_FWRITE  1
#define HIO_FAPPEND 2 /* [add] H.Yashiro 20110907 overwrite mode */

/* data dump type */
#define HIO_DUMP_OFF      0
#define HIO_DUMP_HEADER   1
#define HIO_DUMP_ALL      2
#define HIO_DUMP_ALL_MORE 3

/* return type */
#define ERROR_CODE   -1
#define SUCCESS_CODE  1

#endif
