/*---------------------------------------------------------------
                                        06 Mar. 2013  Shin Satoh
FUNCTION: int read_toshiba
  to read PAWR data in Toshiba original format
----------------------------------------------------------------*/

/*
  modified by Shigenori Otsuka
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#ifdef MacOSX
#include <machine/endian.h>
#else
#include <endian.h>
#endif
#include "read_toshiba.h"

int16_t char2int16(void *input)
{
  int16_t output;

#if __BYTE_ORDER == __LITTLE_ENDIAN
  output = *((int16_t*)input);
#elif __BYTE_ORDER == __BIG_ENDIAN
  output = ((int16_t)(((unsigned char*)input)[1]) << 8) | (int16_t)(((unsigned char*)input)[0]);
#else
  output = ((char*)input)[1] * 0x100 + ((unsigned char*)input)[0];
#endif
  return output;
}

uint16_t char2uint16(void *input)
{
  uint16_t output;

#if __BYTE_ORDER == __LITTLE_ENDIAN
  output = *((uint16_t*)input);
#elif __BYTE_ORDER == __BIG_ENDIAN
  output = ((uint16_t)(((unsigned char*)input)[1]) << 8) | (uint16_t)(((unsigned char*)input)[0]);
#else
  output = ((unsigned char*)input)[1] * 0x100 + ((unsigned char*)input)[0];
#endif
  return output;
}

int32_t char2int32(void *input)
{
  int32_t output;
#if __BYTE_ORDER == __LITTLE_ENDIAN
  output = *((int32_t*)input);
#elif __BYTE_ORDER == __BIG_ENDIAN
  output = ((int32_t)(((unsigned char*)input)[3]) << 24) | ((int32_t)(((unsigned char*)input)[2]) << 16) | ((int32_t)(((unsigned char*)input)[1]) << 8) | (int32_t)(((unsigned char*)input)[0]);
#else
  output = ((char*)input)[3] * 0x1000000 + ((unsigned char*)input)[2] * 0x10000 + ((unsigned char*)input)[1] * 0x100 + ((unsigned char*)input)[0];
#endif
  return output;
}

int read_toshiba(char *fname, pawr_header *hd,
                 float az[ELDIM][AZDIM], float el[ELDIM][AZDIM],
                 float rtdat[ELDIM][AZDIM][RDIM])
{
  const size_t bufsize = 40 * 1024 * 1024; // fixed size
  size_t bsize; // actual data size
  int ierr;
  unsigned char *buf;
  FILE *fp;

  buf = malloc(bufsize);
  if(buf == NULL){
    printf("failed to allocate memory in read_toshiba");
    return -99;
  }
  if((fp = fopen(fname, "r")) == NULL){
    //printf("file not found: %s\n", fname);
    return -9;
  }
  bsize = fread(buf, 1, bufsize, fp);
  if(bsize == 0){
    printf("file size is 0: %s\n", fname);
    return -9;
  }

  ierr = decode_toshiba(bsize, buf, hd, az, el, rtdat);
  if(ierr != 0) return ierr;

  free(buf);

  return 0;
}


int decode_toshiba(size_t bufsize, unsigned char *buf, pawr_header *hd,
                   float az[ELDIM][AZDIM], float el[ELDIM][AZDIM],
                   float rtdat[ELDIM][AZDIM][RDIM])
{
  const size_t len_bufchd = 96;
  const size_t len_bufpcd = 320;
  const size_t len_bufblk = 64;
  unsigned char *bufchd, *bufpcd, *bufblk, *bufdat;
  float *flp;
  float s_az, e_az, s_el, e_el;
  float tx_pilot,return_pilot;
  int   i, j, k, iss, i0, i1, i2, i3, rnum, aznum, elnum;
  int   r_byte, irb, irb1, irb2;
  int   ios_obs, ios_spec, ios_loc, ios_rt, ios_mesh;
  int   ideg, ibunshi, ibunbo;
  int   prf, hit_num, long_pulse;
  const float ideg2deg = (180.0 / 8192.0);
  size_t buf_offset;

/*----- open RT file -----*/
  buf_offset = 0;

/*----- SCAN data (PPI) loop -----*/
  elnum = ELDIM;
  /* elnum=86; */
  for(k = 0; k < elnum; k++) {

/*----- read common header (96 byte)-----*/
    if(buf_offset == bufsize) {
      //printf("# find EOF at elnum=%d\n", k);
      goto READEND;
    } else if(buf_offset + len_bufchd > bufsize) {
      printf("# abnormal record length in common header\n");
      goto READERR;
    }
    bufchd = buf + buf_offset;
    buf_offset += len_bufchd;

    hd->s_yr = 1000 * (bufchd[20] >> 4) + 100 * (bufchd[20] & 15) + 10 * (bufchd[21] >> 4) + (bufchd[21] & 15);
    hd->s_mn = 10 * (bufchd[22] >> 4) + (bufchd[22] & 15);
    hd->s_dy = 10 * (bufchd[23] >> 4) + (bufchd[23] & 15);
    hd->s_hr = 10 * (bufchd[24] >> 4) + (bufchd[24] & 15);
    hd->s_mi = 10 * (bufchd[25] >> 4) + (bufchd[25] & 15);
    hd->s_sc = 10 * (bufchd[26] >> 4) + (bufchd[26] & 15);
    hd->e_yr = 1000 * (bufchd[28] >> 4) + 100 * (bufchd[28] & 15) + 10 * (bufchd[29] >> 4) + (bufchd[29] & 15);
    hd->e_mn = 10 * (bufchd[30] >> 4) + (bufchd[30] & 15);
    hd->e_dy = 10 * (bufchd[31] >> 4) + (bufchd[31] & 15);
    hd->e_hr = 10 * (bufchd[32] >> 4) + (bufchd[32] & 15);
    hd->e_mi = 10 * (bufchd[33] >> 4) + (bufchd[33] & 15);
    hd->e_sc = 10 * (bufchd[34] >> 4) + (bufchd[34] & 15);
    hd->data_size = 0x10000 * bufchd[10] + 0x100 * bufchd[9] + bufchd[8]; // ensure little endian

/*----- read polar-coords-data header (320 byte)-----*/
    if(buf_offset + len_bufpcd > bufsize) {
      printf("# abnormal record length in polar-coords-data header\n");
      goto READERR;
    }
    bufpcd = buf + buf_offset;
    buf_offset += len_bufpcd;

    ios_obs  = 56;   /* offset of obs info */
    ios_spec = 88;   /* offset of spec info */
    ios_loc  = 216;  /* offset of location info */
    ios_rt   = 232;  /* offset of RT coords info */
    ios_mesh = 264;  /* offset of mesh data info */

    hd->total_step_num = bufpcd[ios_obs +  4];
    hd->el_num         = bufpcd[ios_obs + 10];
    hd->total_el_num   = bufpcd[ios_obs + 11];

    // SPEC part
    hd->tx_freq     = 0.0001 * char2int32(bufpcd + ios_spec +   0);
    hd->tx_power    = 0.0001 * char2int32(bufpcd + ios_spec +   8);
    hd->pulse_len_l = 0.0001 * char2int32(bufpcd + ios_spec +  16);
    hd->pulse_len_s = 0.0001 * char2int32(bufpcd + ios_spec +  20);
    hd->ant_gain    = 0.0001 * char2int32(bufpcd + ios_spec +  32);
    hd->beam_wid_h  = 0.0001 * char2int32(bufpcd + ios_spec +  40);
    hd->beam_wid_v  = 0.0001 * char2int32(bufpcd + ios_spec +  48);
    hd->tx_loss     = 0.0001 * char2int32(bufpcd + ios_spec +  56);
    hd->rx_loss     = 0.0001 * char2int32(bufpcd + ios_spec +  64);
    hd->smin_h      = 0.0001 * char2int32(bufpcd + ios_spec +  72);
    hd->smin_l      = 0.0001 * char2int32(bufpcd + ios_spec +  76);
    hd->prf_l       =          char2int32(bufpcd + ios_spec +  88);
    hd->prf_h       =          char2int32(bufpcd + ios_spec +  92);
    hd->zr_b        = 0.0001 * char2int32(bufpcd + ios_spec +  96);
    hd->zr_beta     = 0.0001 * char2int32(bufpcd + ios_spec + 100);

    // LOC part
    hd->latitude  = 0.001 * char2int32(bufpcd + ios_loc + 0);
    hd->longitude = 0.001 * char2int32(bufpcd + ios_loc + 4);
    hd->altitude  = 0.1   * char2int32(bufpcd + ios_loc + 8);

    // RT part
    hd->start_angle = char2int16(bufpcd + ios_rt +  0) * ideg2deg;
    hd->end_angle   = char2int16(bufpcd + ios_rt +  2) * ideg2deg;
    hd->hit_num     = char2uint16(bufpcd + ios_rt +  4);
    hd->sector_num  = char2uint16(bufpcd + ios_rt +  8);
    hd->range_num   = char2uint16(bufpcd + ios_rt + 16);
    hd->range_res   = char2uint16(bufpcd + ios_rt + 20);

    // MESH part
    hd->mesh_size = bufpcd[ios_mesh + 0];
    ibunshi = char2int16(bufpcd + ios_mesh + 4);
    ibunbo  = char2int16(bufpcd + ios_mesh + 8);
    hd->mesh_lsb = (float)ibunshi / (float)ibunbo;
    ibunshi = char2int16(bufpcd + ios_mesh + 12);
    ibunbo  = char2int16(bufpcd + ios_mesh + 16);
    hd->mesh_offset = (float)ibunshi / (float)ibunbo;

    /*--- print header info ---*/
    /*
    if(k == 0){
      printf("## date & time : %d/%d/%d, %d:%d:%d  -",
             hd->s_yr, hd->s_mn, hd->s_dy, hd->s_hr, hd->s_mi, hd->s_sc); 
      printf("  %d/%d/%d, %d:%d:%d\n",
             hd->e_yr, hd->e_mn, hd->e_dy, hd->e_hr, hd->e_mi, hd->e_sc); 
      printf("## data_size=%d, total step num=%d, el num=%d/%d,",
             hd->data_size, hd->total_step_num, hd->el_num, hd->total_el_num); 
      printf(" lat=%8.3f, lon=%8.3f, alt=%8.1f\n",
             hd->latitude, hd->longitude, hd->altitude); 
      printf("## freq=%6.1f, power=%5.3f, pulse_len(L/S)=%4.1f/%3.1f,",
             hd->tx_freq, hd->tx_power, hd->pulse_len_l, hd->pulse_len_s);
      printf(" ant_gain=%5.2f, beam_wid(H/V)=%4.2f/%4.2f, loss(Tx/Rx)=%4.2f/%4.2f\n",
             hd->ant_gain, hd->beam_wid_h, hd->beam_wid_v, hd->tx_loss, hd->rx_loss);
      printf("## smin(H/L)=%6.1f/%6.1f, PRI(L/H)=%6.4f/%6.4f, ZR(B/beta)=%6.2f/%6.4f,",
             hd->smin_h, hd->smin_l, hd->prf_l, hd->prf_h, hd->zr_b, hd->zr_beta);
      printf(" start angle=%5.1f, end angle=%5.1f\n",
             hd->start_angle, hd->end_angle);
      printf("## hit num=%d, sector num=%d, range num=%d,",
             hd->hit_num, hd->sector_num, hd->range_num); 
      printf(" mesh size=%d, LSB=%4.2f, offset=%4.2f\n",
             hd->mesh_size, hd->mesh_lsb, hd->mesh_offset);
    }
    */

    /*----- read data block -----*/
    aznum = hd->sector_num;
    rnum = hd->range_num;

    for(j = 0; j < aznum; j++) {

      /*--- read block header (64 byte) ---*/ 
      if(buf_offset + len_bufblk > bufsize) {
        printf("# abnormal record length in data block\n");
        goto READERR;
      }
      bufblk = buf + buf_offset;
      buf_offset += len_bufblk;

      s_el       = char2int16(bufblk +  4) * ideg2deg;
      s_az       = char2int16(bufblk +  6) * ideg2deg;
      e_el       = char2int16(bufblk +  8) * ideg2deg;
      e_az       = char2int16(bufblk + 10) * ideg2deg;
      prf        = char2int16(bufblk + 14);
      long_pulse = bufblk[34];
      hit_num    = char2int16(bufblk + 36);
      /*
        tx_pilot = -327.68+0.01*(65536*bufpcd[46]+256*bufpcd[45]+bufpcd[44]);
        return_pilot = -327.68+0.01*(65536*bufpcd[50]+256*bufpcd[49]+bufpcd[48]);
      */    

      /* ------------------ <check later: hit num etc.>
         if(j<2){
         printf(" K=%d J=%d: EL=%5.1f-%5.1f, AZ=%5.1f-%5.1f,",
         k,j,s_el,e_el,s_az,e_az); 
         printf(" PRF=%d, hit=%d, long_pulse=%d\n",
         prf,hit_num,long_pulse); 
         }
         -------------------- */

      /*--- read data block ---*/
      irb1 = (hd->data_size / hd->sector_num) - 64;
      irb2 = hd->range_num * (hd->mesh_size / 8);
      irb = irb2;
      bufdat = buf + buf_offset;
      buf_offset += irb;

      /*--- store data ---*/
      el[k][j] = s_el;
      az[k][j] = s_az;
      if(hd->mesh_size == 4){
      } else if(hd->mesh_size == 8){
        for(i = 0; i < rnum; i++) {
          rtdat[k][j][i] = hd->mesh_offset + hd->mesh_lsb * (*((unsigned char*)(bufdat + i)));
        }
      } else if(hd->mesh_size == 16){
        for(i = 0; i < rnum; i++) {
          rtdat[k][j][i] = hd->mesh_offset + hd->mesh_lsb * (char2uint16(bufdat + i * 2));
        }
      }

    } /* end of j-loop */
  } /* end of k-loop */

    /*--- print el & az info ---
      for(k=0; k<elnum; k++) {
      printf("EL No=%3d: EL=%5.2f %5.2f %5.2f | AZ=%5.1f %5.1f %5.1f %5.1f %5.1f --- ",
      k+1,el[k][0],el[k][10],el[k][aznum-1],az[k][0],az[k][1],az[k][2],az[k][3],az[k][4]);
      printf("%5.1f %5.1f %5.1f %5.1f %5.1f\n",
      az[k][aznum-5],az[k][aznum-4],az[k][aznum-3],az[k][aznum-2],az[k][aznum-1]);
      }
      --- */

  goto READEND;

 READERR:
  return -8;

 READEND:
  return 0;
}
