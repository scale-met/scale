#ifndef __POH5_H__
#define __POH5_H__


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hdf5.h"

/* #include "hdf5.h" */

#define POH5_VERSION 93

#ifndef POH5_DO_SHUF
#define POH5_DO_SHUF
#endif

#ifndef POH5_DO_GZIP
#define POH5_DO_GZIP 6
#endif

/*
 * \note below must be consistent with the one onf hio.h
 *
 */

/* character length */
#define POH5_HLNG 256
#define POH5_HMID  64
#define POH5_HSHT  16

/* data type */
#define POH5_REAL4     0
#define POH5_REAL8     1
#define POH5_INTEGER4  2
#define POH5_INTEGER8  3

/* flie open mode */
#define POH5_FREAD   0
#define POH5_FWRITE  1
#define POH5_FAPPEND 2

typedef float (real32_t);
typedef double(real64_t);

typedef int phid_t; /* must be consistent with hid_t defined in H5Ipoublic.h */

/*========================================================================*/
/** \brief Setup poh5 library.
 */
extern int poh5_setup();




/*========================================================================*/
/** \brief Open poh5 file.
 *
 * Open poh5 file with given `fname`.
 *
 * `mode` must be one of below:
 * - POH5_FREAD: read only (no modification allowed)
 * - POH5_FWRITE: newly created (existing file will be discarded)
 * - POH5_FAPPEND: open to append/modify
 * .
 * \return hdf5 file id.
 *
 */
extern hid_t poh5_open_file(
                     const char *fname, /**< [in] file name */
                     const int mode     /**< [in] file open mode */
                      );

/*========================================================================*/
/** \brief Open poh5 file.
 *
 * Close poh5 file opened bu poh5_close_file
 *
 */
extern int poh5_close_file(
                    const hid_t fid /**< [in] file id */
                    ) ;


/*========================================================================*/
/** \brief Write global attributes.
 *
 * Global attributes written to given `file_id`.
 *
 * \return 0 on success
 *
 */
extern int poh5_write_global_attr(
                      const hid_t    file_id,       /**< [in] poh5 file_id */
                      const int32_t  glevel,        /**< [in] glevel */
                      const int32_t  rlevel,        /**< [in] rlevel */
                      const int32_t  grid_topology, /**< [in] grid_topology */
                      const int      is_complete,   /**< [in] 1 if complete file */
                      const int32_t  num_of_rgn,    /**< [in] num of region this file contains. */
                      const int32_t *rgnid,         /**< [in] array of region ids */
                      const char    *description,   /**< [in] description of this file */
                      const char    *note,          /**< [in] longer note of this file */
                      const int32_t  num_of_var   /**< [in] number of data in this file. */
                                  );


/*========================================================================*/
/**  \brief Read global attributes.
 *
 * Global attributes read from given `file_id`.
 *
 * \note rgnid[] is allocated in this routine, previous memories will be free.
 *
 * \return 0 on success.
 *
 */
extern int poh5_read_global_attr(
                      const hid_t file_id,    /**< [in] poh5 file_id */
                      int32_t *glevel,        /**< [out] glevel */
                      int32_t *rlevel,        /**< [out] rlevel */
                      int32_t *grid_topology, /**< [out] grid_topology */
                      int     *is_complete,   /**< [out] 1 if complete file */
                      int32_t *num_of_rgn,    /**< [out] num of region this file contains. */
                      int32_t *rgnid[],       /**< [out] array of region id's */
                      char    description[64],   /**< [out] description of this file */
                      char    note[256] ,         /**< [out] longer note of this file */
                      int32_t *num_of_var  /**< [out] number of data in this file. */
                                 );




/*========================================================================*/
/** \brief Create new variable for write.
 *
 * Create new dataset for one variable.
 *
 * Open entry for one variable, writeout attributes,
 *
 * - length of name < POH5_HSHT
 * - length of desc < POH5_HMID
 * - length of note < POH5_HLNG
 * - length of unit < POH5_HSHT
 * - length of lname < POH5_HSHT
 *
 * example of layername:  "ZSSFC1","ZSDEF40","ZSALL42",etc.
 *
 * \return group_id for this variable.
 *
 */

extern hid_t poh5_create_variable(
                    const hid_t file_id,  /**< [in] poh5 file id*/
                    const char *name,     /**< [in] variable name */
                    const char *dscr,     /**< [in] variable description */
                    const char *note,     /**< [in] variable long note */
                    const char *unit,     /**< [in] variable unit */
                    const char *lname,    /**< [in] name of vertical layer */
                    const int32_t nlayer, /**< [in] number of layer */
                    const int gall1d,     /**< [in] gall1d */
                    const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                    const int num_of_rgn  /**< [in] number of regions */
                                  );



/*========================================================================*/
/** \brief Write one variable data.
 *
 * Write one variable data to given `var_gid`.
 *
 * Origina rank of `var_data` is 4 and denoted as var_data(i,j,k,l) in Fortran.
 * Dimension of var_data is [ gall1d x gall1d x nlayer x num_of_rgn ].
 *
 *
 * It is necessary to specify type of data as `dtype`, which is one of
 * - POH5_REAL4
 * - POH5_REAL8
 * - POH5_INTEGER4
 * - POH5_INTEGER8
 * .
 *
 */
extern void poh5_write_variable_data(
                           const hid_t var_gid,     /**< gid for variable opened by open_variable().*/
                           const int32_t step,      /**< current step???? */
                           const int64_t time_start,/**< start time this data represents */
                           const int64_t time_end,  /**< end time this data represents */
                           const int dtype,         /**< [in] POH5_{REAL|INTEGER}{4|8} */
                           const void *var_data     /**< variable data at current step */
                                     );


/*========================================================================*/
/** \brief Write one variable data of 1rgn
 *
 * Write one variable data to given `var_gid`.
 *
 * Original rank of `var_data` is 3 and denoted as var_data(i,j,k) in Fortran.
 * Dimension of var_data is [ gall1d x gall1d x nlayer ].
 *
 *
 * It is necessary to specify type of data as dtype, which is one of
 * - POH5_REAL4
 * - POH5_REAL8
 * - POH5_INTEGER4
 * - POH5_INTEGER8
 * .
 *
 */
extern void poh5_write_variable_data_1rgn(
                              const hid_t v_gid, /**<[in] poh5 file id .*/
                              const int32_t step,      /**< step count */
                              const int32_t ll, /**< rgn number, 0-based */
                              const int64_t time_start, /**< start time this data represents */
                              const int64_t time_end,   /**< end time this data represents */
                              const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                              const void *var_data   /**< variable data at current step */
                                          );

/*========================================================================*/
/** \brief Open variable for read.
 *
 *  Open variable with spacified `vname`.
 *
 * \return group id of its variable.
 *
 */
extern hid_t poh5_open_variable(
                    const hid_t file_id, /**<[in] poh5 file id .*/
                    const char *vname   /**<[in] variable name */
                    );


/*========================================================================*/
/** \brief Open variable for read by index.
 *
 *  Open 'idx'th variable.
 *
 * \note `idx` is 0-based.
 *
 * \return group id of the variable.
 *
 */
extern hid_t poh5_open_variable_by_idx(
                    const hid_t file_id, /**<[in] poh5 file id .*/
                    const int idx      /**<[in] variable index */
                    );


/*========================================================================*/
/** \brief Close variable.
 *
 * Close variable opened by poh5_open_variable() or poh5_open_variable_by_idx().
 *
 * \return returned value of H5Gclose();
 */
extern int poh5_close_variable(
                               const hid_t v_gid); /**<[in] group id .*/


/*========================================================================*/
/** \brief Read attributes of variable.
 *
 * Read attributes of given variable specified by `var_gid` and `step`.
 *
 * \return 0 on success.
 *
 */
extern int poh5_read_variable_attr(
                   const hid_t var_gid, /**< [in] group id of variable */
                   char *name,       /**< [out] variable name  */
                   char *dscr,       /**< [out] variable description */
                   char *note,          /**< [out] variable long note */
                   char *unit,          /**< [out] variable unit */
                   char *lname,         /**< [out] name of vertical layer */
                   int32_t *nlayer,     /**< [out] number of layers */
                   int32_t *nsteps,     /**< [out] number of steps */
                   int *dtype)          /**< [out] POH5_{REAL4,REAL8,INTEGER4,INTEGER8} */
  ;


/*========================================================================*/
/** \brief Read out variable time.
 *
  * Read out time_start/end of variable data specfied var_gid and step.
  *
  * \return 0 on success.
  *
  */
int poh5_read_variable_time(
                   const hid_t v_gid, /**< [in] group id of variable */
                   int64_t ts[],         /**< [out] start time of this step */
                   int64_t te[]);         /**< [out] end time of this step */



/*========================================================================*/
/** \brief Read out variable data.
 *
  * Read out variable data specfied var_gid and step.
  *
  * Also returns time_start and time_end of given step.
  *
  * \return 0 on success.
  *
  */
extern int poh5_read_variable_data(
                   const hid_t var_gid, /**< [in] group id of variable */
                   const int32_t step,        /**<  [in] step counter */
                   int64_t *time_start, /**< [out] start time of this step */
                   int64_t *time_end,   /**< [out] end time of this step */
                   const int dtype,     /**< [in] POH5_{REAL4,REAL8,INTEGER4,INTEGER8} */
                   void *vdata)          /**< [out] data */
  ;


#endif /* __POH5_H__ */
