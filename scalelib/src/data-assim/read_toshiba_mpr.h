/*----------------------------------------------------------------------
     FILE NAME: grads_MPR.h
----------------------------------------------------------------------*/
/* constant fixed parameter */
#define  RDIM     800      /* max number of Range bins */
#define  AZDIM    320      /* max number of AZ angles */
#define  ELDIM    120      /* max number of EL angles */
#define  NFNAME   128      /* number of characters in file names */
#define  DMISS   -327.68   /* missing data in output (& input data offset) */
#define  DNOISE  -327.00   /* noise level data in output */

#define  RAD   (3.141592653589/180.0)

// added by Otsuka
#define LEN_BUFRDHD 512
#define LEN_BUFPBHD 10368
#define LEN_BUFELHD 432

/* Dimension of output XYZ(CAPPI), XY(PPI) and XZ(RHI) */
#define  XYRESO     0.25   /* horizontal resolution in km */
#define  ZRESO      0.25   /* vertical resolution in km */

#define  CXDIM    481      /* X-dim of CAPPI output (include center) */
#define  CYDIM    481      /* Y-dim of CAPPI output (include center) */
#define  CZDIM     65      /* Z-dim of CAPPI output (include surface) */
#define  CXKM0    -60.0    /* start X coords from the radar (in km) */
#define  CYKM0    -60.0    /* start Y coords from the radar (in km) */
#define  CZKM0      0.0    /* start Z coords from the radar (in km) */

#define  PXDIM    481      /* X-dim of PPI output (include center) */
#define  PYDIM    481      /* Y-dim of PPI output (include center) */
#define  PELDIM    20      /* EL-dim of output EL number (< ELDIM) */
#define  PXKM0    -60.0    /* start X coords from the radar (in km) */
#define  PYKM0    -60.0    /* start Y coords from the radar (in km) */

#define  RXDIM    241      /* X-dim of RHI output (include center) */
#define  RZDIM     65      /* Z-dim of RHI output (include surface) */
#define  RAZDIM   300      /* AZ-dim of RHI output (< AZDIM) */
#define  RZKM0      0.0    /* start Z coords from the radar (in km) */

typedef struct {
  char   data_name[32], site_name[32], sq_name[16];
  int    s_yr, s_mn, s_dy, s_hr, s_mi, s_sc;
  int    e_yr, e_mn, e_dy, e_hr, e_mi, e_sc;
  int    el_num, ray_num, range_num;
  int    range_res;
  double latitude, longitude, altitude;
  float  start_az, start_el, end_az, end_el;
  float  mesh_offset;
} mppawr_header;

int16_t char2int16(void *input);
uint16_t char2uint16(void *input);
int32_t char2int32(void *input);

int read_toshiba_mpr(char *in_file, 
                     int opt_verbose,
                     mppawr_header *hd,
                     float az[ELDIM][AZDIM], 
                     float el[ELDIM][AZDIM],
                     float rtdat[ELDIM][AZDIM][RDIM]);

int decode_toshiba_mpr(size_t bufsize, unsigned char *buf,
                       int opt_verbose,
                       mppawr_header *hd,
                       float az[ELDIM][AZDIM],
                       float el[ELDIM][AZDIM],
                       float rtdat[ELDIM][AZDIM][RDIM]);

int ggrads_ctl(char *c_name, char *g_name, int nx, int ny, int nz, int nt,
               double xmin, double ymin, double zmin, double tmin,
               double dx, double dy, double dz, double dt, double no_data,
               int ihour, int imin, int iday, int imonth, int iyear,
               char *c_var, char *c_memo);
