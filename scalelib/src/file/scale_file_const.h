/* character length */
#define File_HSHORT   32
#define File_HMID    128
#define File_HLONG  1024

/* data type */
#define File_REAL4     0
#define File_REAL8     1
#define File_INTEGER2  2
#define File_INTEGER4  3
#define File_INTEGER8  4

/* action type */
#define File_FREAD   0
#define File_FWRITE  1
#define File_FAPPEND 2

/* return type */
#define ERROR_CODE          -1
#define SUCCESS_CODE         0
#define ALREADY_CLOSED_CODE  1
#define ALREADY_EXISTED_CODE 2

/* limit of files */
#define FILE_MAX 512
/* limit of variables */
#define VAR_MAX 40960
/* limit of dimensions */
#define RANK_MAX 10

/* missing value */
#define RMISS -9.9999e+30
