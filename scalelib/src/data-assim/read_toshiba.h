/*----------------------------------------------------------------------
     FILE NAME: read_toshiba.h
     created by Shinsuke Satoh and modified by Shigenori Otsuka
----------------------------------------------------------------------*/

#define  RDIM     600      /* max number of Range bins */
#define  AZDIM    320      /* max number of AZ angles */
#define  ELDIM    121      /* max number of AZ angles */
#define  NFNAME   128      /* number of characters for file names */
#define  DMISS   -327.68   /* missing data in output (& input data offset) */
#define  DNOISE  -327.00   /* noise level data in output */

#define  RAD   (3.141592653589/180.0)

typedef struct {
  int        s_yr, s_mn, s_dy, s_hr, s_mi, s_sc;
  int        e_yr, e_mn, e_dy, e_hr, e_mi, e_sc;
  int        data_size;
  int        total_step_num, el_num, total_el_num;
  int        hit_num, sector_num, range_num, range_res, mesh_size;
  double     latitude, longitude, altitude;
  float      start_angle, end_angle, mesh_lsb, mesh_offset;
  float      tx_freq, tx_power, pulse_len_l, pulse_len_s;
  float      ant_gain, beam_wid_h, beam_wid_v;
  float      tx_loss, rx_loss, smin_h, smin_l;
  float      prf_l, prf_h, zr_b, zr_beta;
} pawr_header;

/*
int jitdt_read_toshiba(int n_type, char *jitdt_place, pawr_header hd[n_type],
                       float az[n_type][ELDIM][AZDIM], float el[n_type][ELDIM][AZDIM],
                       float rtdat[n_type][ELDIM][AZDIM][RDIM]);
*/

int decode_toshiba(size_t bufsize, unsigned char *buf, pawr_header *hd,
                   float az[ELDIM][AZDIM], float el[ELDIM][AZDIM],
                   float rtdat[ELDIM][AZDIM][RDIM]);
