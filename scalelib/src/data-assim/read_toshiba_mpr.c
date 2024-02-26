/*---------------------------------------------------------------
                                        15 Oct 2019  Shin Satoh
FUNCTION: int read_toshiba_mpr
  to read MP-PAWR RAW data in Toshiba original format
  opt_verbose = 0: no printout (silent)
                1: print common (RAW DATA) header info
                2: print Polar Block (RAY) header info
                3: print EL header info
----------------------------------------------------------------*/
/*
  modified by Shigenori Otsuka
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#ifdef MacOSX
#include <machine/endian.h>
#else
#include <endian.h>
#endif
#include <zlib.h>
#include "read_toshiba_mpr.h"

/*
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
*/

size_t ungzip_toshiba_mpr(size_t outbufsize, size_t bufsize, unsigned char *buf){
  unsigned char *outbuf;
  size_t datsize;
  z_stream strm;
  int ret;

#ifdef DA
  outbuf = (unsigned char*)malloc(outbufsize);
  if(outbuf == NULL){
    printf("malloc failed in ungzip_toshiba_mpr\n");
    datsize = 0;
  } else {
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    ret = inflateInit2(&strm, 47);

    strm.next_in = buf;
    strm.avail_in = bufsize;
    strm.next_out = outbuf;
    strm.avail_out = outbufsize;
    ret = inflate(&strm, Z_NO_FLUSH);

    datsize = outbufsize - strm.avail_out;

    memcpy(buf, outbuf, datsize);
    printf("maxbuf: %d, inbuf: %d, outbuf: %d\n",
           outbufsize, bufsize, datsize);

    ret = inflateEnd(&strm);
    free(outbuf);
  }
#endif
  return datsize;
}


int read_toshiba_mpr(char *in_file,
                     int opt_verbose,
                     mppawr_header *hd,
                     float az[ELDIM][AZDIM],
                     float el[ELDIM][AZDIM],
                     float rtdat[ELDIM][AZDIM][RDIM])
{
  const size_t bufsize = AZDIM * (ELDIM * (RDIM * 2 + LEN_BUFELHD) + LEN_BUFPBHD) + LEN_BUFRDHD; // max size
  size_t bsize; // actual data size
  int ierr;
  unsigned char *buf;
  FILE *fp;
  char *is_gzip;

  buf = malloc(bufsize);
  if(buf == NULL){
    printf("failed to allocate memory in read_toshiba");
    return -99;
  }
  if((fp = fopen(in_file, "r")) == NULL){
    is_gzip = "true";
    if((fp = fopen(strcat(in_file,".gz"), "r")) == NULL){
      //printf("file not found: %s\n", in_file);
      return -9;
    }
  }
  bsize = fread(buf, 1, bufsize, fp);
  if(bsize == 0){
    printf("file size is 0: %s\n", in_file);
    return -9;
  }

  if(is_gzip == "true") bsize = ungzip_toshiba_mpr(bufsize, bsize, buf);

  ierr = decode_toshiba_mpr(bsize, buf, opt_verbose, hd, az, el, rtdat);
  if(ierr != 0) return ierr;

  free(buf);

  return 0;
}


int decode_toshiba_mpr(size_t bufsize, unsigned char *buf,
                       int opt_verbose,
                       mppawr_header *hd,
                       float az[ELDIM][AZDIM],
                       float el[ELDIM][AZDIM],
                       float rtdat[ELDIM][AZDIM][RDIM])
{
  unsigned char *bufrdhd; /* RAW DATA hdr */
  unsigned char *bufpbhd; /* POLAR BLOCK hdr */
  unsigned char *bufelhd; /* EL hdr */
  unsigned char *bufdat;  /* 2 byte x 800 range (max) */
  int   i, j, k, iss, i0, i1, i2, i3;
  int   aznum, elnum, rnum, sp_rnum, lp_rnum;
  int   r_byte, irb, mesh_size, data_size;
  int   ideg, ideg_tc, inum, ibunshi, ibunbo;
  int   s1_yr, s1_mn, s1_dy, s1_hr, s1_mi, s1_sc;
  int   s2_yr, s2_mn, s2_dy, s2_hr, s2_mi, s2_sc;
  int   s3_yr, s3_mn, s3_dy, s3_hr, s3_mi, s3_sc;
  int   el_num_sec, cpi_num, sec_num, rnum1;
  int   hit_num, spa_num, lpa_num, spv_num, lpv_num;
  int   pw_hsp, pw_hlp, pw_vsp, pw_vlp, prf;
  float mesh_lsb, mesh_offset;
  float start_az_sec, end_az_sec, start_el_sec, end_el_sec;
  float start_az_cpi, end_az_cpi, start_el_cpi, end_el_cpi;
  float lower_el_angle, upper_el_angle;
  float tx_freq_h, tx_freq_v, tx_power_h, tx_power_v;
  float tx_plen_hs, tx_plen_hl, tx_plen_vs, tx_plen_vl;
  float tx_antgain_h, tx_antgain_v, rx_antgain_h, rx_antgain_v;
  float hori_beamwidth_h, hori_beamwidth_v, vert_beamwidth_h, vert_beamwidth_v;
  float start_obs_az, start_obs_el, end_obs_az, end_obs_el;
  // added by Otsuka
  size_t buf_offset;
  const float ideg2deg = 360.0 / 65536.0;
  const float ideg_tc2deg = 180.0 / 32768.0;
  int tmp_range_res;

  // init param
  hd->range_res = -1;

  /*----- open RAW DATA file -----*/
  buf_offset = 0;

  /*----- read RAW DATA header (512 byte)-----*/
  if(buf_offset == bufsize) {
    goto READEND;
  } else if(buf_offset + LEN_BUFRDHD > bufsize) {
    printf("# abnormal record length (RDhd)\n");
    goto READERR;
  }
  bufrdhd = buf + buf_offset;
  buf_offset += LEN_BUFRDHD;

  strncpy(hd->data_name, &bufrdhd[16], 32);
  hd->s_yr = 1000 * (bufrdhd[48] >> 4) + 100 * (bufrdhd[48] & 15) + 10 * (bufrdhd[49] >> 4) + (bufrdhd[49] & 15);
  hd->s_mn = 10 * (bufrdhd[50] >> 4) + (bufrdhd[50] & 15);
  hd->s_dy = 10 * (bufrdhd[51] >> 4) + (bufrdhd[51] & 15);
  hd->s_hr = 10 * (bufrdhd[52] >> 4) + (bufrdhd[52] & 15);
  hd->s_mi = 10 * (bufrdhd[53] >> 4) + (bufrdhd[53] & 15);
  hd->s_sc = 10 * (bufrdhd[54] >> 4) + (bufrdhd[54] & 15);
  hd->e_yr = 1000 * (bufrdhd[56] >> 4) + 100 * (bufrdhd[56] & 15) + 10 * (bufrdhd[57] >> 4) + (bufrdhd[57] & 15);
  hd->e_mn = 10 * (bufrdhd[58] >> 4) + (bufrdhd[58] & 15);
  hd->e_dy = 10 * (bufrdhd[59] >> 4) + (bufrdhd[59] & 15);
  hd->e_hr = 10 * (bufrdhd[60] >> 4) + (bufrdhd[60] & 15);
  hd->e_mi = 10 * (bufrdhd[61] >> 4) + (bufrdhd[61] & 15);
  hd->e_sc = 10 * (bufrdhd[62] >> 4) + (bufrdhd[62] & 15);
  /*--- Scan-start(s1), Sequence-start(s2), Mode-start(s3) time-- */
  s1_yr = 1000 * (bufrdhd[64] >> 4) + 100 * (bufrdhd[64] & 15) + 10 * (bufrdhd[65] >> 4) + (bufrdhd[65] & 15);
  s1_mn = 10 * (bufrdhd[66] >> 4) + (bufrdhd[66] & 15);
  s1_dy = 10 * (bufrdhd[67] >> 4) + (bufrdhd[67] & 15);
  s1_hr = 10 * (bufrdhd[68] >> 4) + (bufrdhd[68] & 15);
  s1_mi = 10 * (bufrdhd[69] >> 4) + (bufrdhd[69] & 15);
  s1_sc = 10 * (bufrdhd[70] >> 4) + (bufrdhd[70] & 15);
  s2_yr = 1000 * (bufrdhd[72] >> 4) + 100 * (bufrdhd[72] & 15) + 10 * (bufrdhd[73] >> 4) + (bufrdhd[73] & 15);
  s2_mn = 10 * (bufrdhd[74] >> 4) + (bufrdhd[74] & 15);
  s2_dy = 10 * (bufrdhd[75] >> 4) + (bufrdhd[75] & 15);
  s2_hr = 10 * (bufrdhd[76] >> 4) + (bufrdhd[76] & 15);
  s2_mi = 10 * (bufrdhd[77] >> 4) + (bufrdhd[77] & 15);
  s2_sc = 10 * (bufrdhd[78] >> 4) + (bufrdhd[78] & 15);
  s3_yr = 1000 * (bufrdhd[80] >> 4) + 100 * (bufrdhd[80] & 15) + 10 * (bufrdhd[81] >> 4) + (bufrdhd[81] & 15);
  s3_mn = 10 * (bufrdhd[82] >> 4) + (bufrdhd[82] & 15);
  s3_dy = 10 * (bufrdhd[83] >> 4) + (bufrdhd[83] & 15);
  s3_hr = 10 * (bufrdhd[84] >> 4) + (bufrdhd[84] & 15);
  s3_mi = 10 * (bufrdhd[85] >> 4) + (bufrdhd[85] & 15);
  s3_sc = 10 * (bufrdhd[86] >> 4) + (bufrdhd[86] & 15);
  strncpy(hd->site_name, &bufrdhd[136], 32);
  hd->latitude  = 0.000001 * char2int32(bufrdhd + 168 + 0);
  hd->longitude = 0.000001 * char2int32(bufrdhd + 172 + 0);
  hd->altitude  = 0.01 *     char2int32(bufrdhd + 180 + 0);
  strncpy(hd->sq_name, &bufrdhd[216], 16);

  hd->el_num    = bufrdhd[236 + 1];
  hd->ray_num   = char2uint16(bufrdhd + 238); //unsigned
  hd->start_az  = char2uint16(bufrdhd + 240) * ideg2deg; //unsigned
  hd->start_el  = char2int16(bufrdhd + 242) * ideg_tc2deg;
  hd->end_az    = char2int16(bufrdhd + 244) * ideg2deg; //unsigned
  hd->end_el    = char2int16(bufrdhd + 246) * ideg_tc2deg;
  mesh_size     = char2int16(bufrdhd + 384) * 8; /* mesh_size stored in byte? */
  ibunshi       = char2int32(bufrdhd + 388);
  ibunbo        = char2int32(bufrdhd + 392);
  mesh_lsb      = (float)ibunshi / (float)ibunbo;
  ibunshi       = char2int32(bufrdhd + 396);
  ibunbo        = char2int32(bufrdhd + 400);
  mesh_offset   = (float)ibunshi / (float)ibunbo;
  hd->range_num = RDIM; /* <== shoud be fixed num, becasue of rtdat's base dim */

  hd->mesh_offset = mesh_offset; //save offset

  /*--- print RAW DATA header ---*/
  if (opt_verbose >= 1){
    //printf("## input file: %s\n", in_file); 
    printf("## data name: %s\n", hd->data_name); 
    printf("## scan collection start/end time: ");
    printf("%d/%d/%d, %d:%d:%d - %d/%d/%d, %d:%d:%d\n",
           hd->s_yr, hd->s_mn, hd->s_dy, hd->s_hr, hd->s_mi, hd->s_sc, 
           hd->e_yr, hd->e_mn, hd->e_dy, hd->e_hr, hd->e_mi, hd->e_sc);
    printf("## scan start time (s1): %d/%d/%d, %d:%d:%d\n",
           s1_yr, s1_mn, s1_dy, s1_hr, s1_mi, s1_sc);
    printf("## sequence start time (s2): %d/%d/%d, %d:%d:%d\n",
           s2_yr, s2_mn, s2_dy, s2_hr, s2_mi, s2_sc);
    printf("##     mode start time (s3): %d/%d/%d, %d:%d:%d\n",
           s3_yr, s3_mn, s3_dy, s3_hr, s3_mi, s3_sc);
    printf("## site name: %s", hd->site_name); 
    printf("   lat=%7.4f, lon=%8.4f, alt=%6.2f\n",
           hd->latitude, hd->longitude, hd->altitude); 
    printf("## sequence name: %s", hd->sq_name); 
    printf("   el_num=%d, ray_num=%d\n", hd->el_num, hd->ray_num);
    /*--- 
      printf("## start/end AZ angle=%5.1f-%5.1f ",hd->start_az, hd->end_az);
      printf("   start/end EL angle=%5.2f-%5.2f\n",hd->start_el, hd->end_el);
      ---*/
    printf("## mesh size=%d, LSB=%7.3f, offset=%7.3f\n",
           mesh_size, mesh_lsb, mesh_offset);
  }

  if(mesh_size != 16){  /* 2-byte data */
    printf("# abnormal mesh size: %d bit\n", mesh_size);
    printf("# this program (read_toshiba_mpr) is not supported !\n");
    goto READERR;
  }

  /*----- RAY data (PPI) loop -----*/
  aznum = hd->ray_num;
  if(aznum > AZDIM){
    printf("# abnormal AZ dimension: %d (> %d)\n", aznum, AZDIM);
    goto READERR;
  }
  for(j = 0; j < aznum; j++) {    /*= AZ loop uses j =*/
    /*----- read POLAR BLOCK (RAY) header (10368 byte)-----*/
    if(buf_offset == bufsize) {
      goto READEND;
    } else if(buf_offset + LEN_BUFPBHD > bufsize) {
      printf("# abnormal record length (PB header)\n");
      goto READERR;
    }
    bufpbhd = buf + buf_offset;
    buf_offset += LEN_BUFPBHD;

    /*--- get PB header info for print ---*/
    if (opt_verbose >= 2){
      el_num_sec     = bufpbhd[9]; //unsigned
      cpi_num        = char2uint16(bufpbhd + 12); //unsigned
      sec_num        = char2uint16(bufpbhd + 14); //unsigned
      start_az_sec   = char2uint16(bufpbhd + 16) * ideg2deg; //unsigned
      start_el_sec   = char2int16(bufpbhd + 18) * ideg_tc2deg;
      end_az_sec     = char2uint16(bufpbhd + 20) * ideg2deg; //unsigned
      end_el_sec     = char2int16(bufpbhd + 22) * ideg_tc2deg;
      lower_el_angle = char2int16(bufpbhd + 26) * ideg_tc2deg;
      upper_el_angle = char2int16(bufpbhd + 30) * ideg_tc2deg;
      if(lower_el_angle > 180.0) lower_el_angle -= 360.0;

      /*--- print PB header ---*/
      printf("## Polar Block (RAY) header: j=%d\n", j); 
      printf("    el_num_sec=%d, cpi_num=%d, sec_num=%d  | rnum=%d\n", 
             el_num_sec, cpi_num, sec_num, rnum); 
      printf("    start/end AZ in sector=%6.2f-%6.2f ",
             start_az_sec, end_az_sec);
      printf(" start/end EL of pedestal=%6.2f-%6.2f ",
             start_el_sec, end_el_sec);
      printf(" lower/upper EL angle=%6.2f-%6.2f\n",
             lower_el_angle, upper_el_angle);
    }

    /*----- EL data (RHI) loop -----*/
    elnum = hd->el_num;
    if(elnum > ELDIM){
      printf("# abnormal EL dimension: %d (> %d)\n", elnum, ELDIM);
      goto READERR;
    }
    for(k = 0; k < elnum; k++) {   /*= EL loop uses k =*/

      /*----- read EL header (432 byte)-----*/
      if(buf_offset == bufsize) {
        goto READEND;
      } else if(buf_offset + LEN_BUFELHD > bufsize) {
        printf("# abnormal record length (EL header)\n");
        goto READERR;
      }
      bufelhd = buf + buf_offset;
      buf_offset += LEN_BUFELHD;

      /***   if (opt_verbose >= 3 && (j==0 || j==1) && (k==3 || k==4)){  ***/
      if (opt_verbose >= 3){

        /*--- get EL header info for print ---*/
        data_size    = char2uint16(bufelhd + 0); //unsigned
        hit_num      = char2uint16(bufelhd + 8); //unsigned
        cpi_num      = char2uint16(bufelhd + 10); //unsigned
        spa_num      = char2uint16(bufelhd + 12); /*short pulse align range num*/ //unsigned
        lpa_num      = char2uint16(bufelhd + 14); /*long pulse align range num*/ //unsigned
        spv_num      = char2uint16(bufelhd + 16); /*short pulse valid range num*/ //unsigned
        lpv_num      = char2uint16(bufelhd + 18); /*long pulse valid range num*/ //unsigned
        start_az_cpi = char2uint16(bufelhd + 128) * ideg2deg; //unsigned
        start_el_cpi = char2int16(bufelhd + 130) * ideg_tc2deg;
        end_az_cpi   = char2uint16(bufelhd + 132) * ideg2deg; //unsigned
        end_el_cpi   = char2int16(bufelhd + 134) * ideg_tc2deg;
        /*short pulse: 0=0.5us, 1=1.0us, 2=2.0us, 3=4.0us*/
        /*long pulse: 0=24us, 1=32us, 2=48us, 3=72us*/
        pw_hsp = bufelhd[144]; /*pulse width of H-pol short pulse*/
        pw_hlp = bufelhd[145]; /*pulse width of H-pol long pulse*/
        pw_vsp = bufelhd[146]; /*pulse width of V-pol short pulse*/
        pw_vlp = bufelhd[147]; /*pulse width of V-pol long pulse*/
        prf = char2uint16(bufelhd + 148); //unsigned
        tx_freq_h        = 0.0001 * char2int32(bufelhd + 176);
        tx_freq_v        = 0.0001 * char2int32(bufelhd + 180);
        tx_power_h       = 0.0001 * char2int32(bufelhd + 184);
        tx_power_v       = 0.0001 * char2int32(bufelhd + 188);
        tx_plen_hs       = 0.0001 * char2int32(bufelhd + 192);
        tx_plen_hl       = 0.0001 * char2int32(bufelhd + 196);
        tx_plen_vs       = 0.0001 * char2int32(bufelhd + 200);
        tx_plen_vl       = 0.0001 * char2int32(bufelhd + 204);
        tx_antgain_h     = 0.0001 * char2int32(bufelhd + 208);
        tx_antgain_v     = 0.0001 * char2int32(bufelhd + 212);
        rx_antgain_h     = 0.0001 * char2int32(bufelhd + 216);
        rx_antgain_v     = 0.0001 * char2int32(bufelhd + 220);
        hori_beamwidth_h = 0.0001 * char2int32(bufelhd + 224);
        hori_beamwidth_v = 0.0001 * char2int32(bufelhd + 228);
        vert_beamwidth_h = 0.0001 * char2int32(bufelhd + 232);
        vert_beamwidth_v = 0.0001 * char2int32(bufelhd + 236);

         /*--- print EL header ---*/
        printf("## EL header (%d byte): j=%d  k=%d | ", data_size, j, k); 
        printf(" hit_num=%d, cpi_num=%d, prf=%d,", 
               hit_num, cpi_num, prf); 
        printf(" spv/a=%d/%d, lpv/a=%d/%d\n", 
               spv_num, spa_num, lpv_num, lpa_num); 

        if (opt_verbose >= 4){
          printf("    start/end AZ in cpi=%6.2f-%6.2f\n",
                 start_az_sec, end_az_sec);
          printf("    start/end EL in cpi=%6.2f-%6.2f\n",
                 start_el_sec, end_el_sec);
          printf("  < pw_short pulse: 0=0.5us, 1=1.0us, 2=2.0us, 3=4.0us >\n");
          printf("  < pw_long pulse: 0=24us, 1=32us, 2=48us, 3=72us >\n");
          printf("    pw_hsp=%d, pw_hlp=%d, pw_vsp=%d, pw_vlp=%d\n", 
                 pw_hsp, pw_hlp, pw_vsp, pw_vlp); 
          printf("    tx_freq_h=%.1f, tx_freq_v=%.1f, ", 
                 tx_freq_h, tx_freq_v);
          printf("tx_power_h=%.1f, tx_power_v=%.1f\n", 
                 tx_power_h, tx_power_v);
          printf("    tx_plen_hs=%.1f, tx_plen_hl=%.1f, ", 
                 tx_plen_hs, tx_plen_hl);
          printf("tx_plen_vs=%.1f, tx_plen_vl=%.1f\n", 
                 tx_plen_vs, tx_plen_vl);
          printf("    tx_antgain_h=%.2f, tx_antgain_v=%.2f, ", 
                 tx_antgain_h, tx_antgain_v);
          printf("rx_antgain_h=%.2f, rx_antgain_v=%.2f\n", 
                 rx_antgain_h, rx_antgain_v);
          printf("    hori_beamwidth_h=%.2f, hori_beamwidth_v=%.2f, ", 
                 hori_beamwidth_h, hori_beamwidth_v);
          printf("vert_beamwidth_h=%.2f, vert_beamwidth_v=%.2f\n", 
                 vert_beamwidth_h, vert_beamwidth_v);
        }
      } /* end of if (opt_verbose >= 3) */ 

      /* range resolution (RAY) */
      tmp_range_res = char2uint16(bufelhd + 372); //unsigned
      if((hd->range_res > 0) && (hd->range_res != tmp_range_res)){
        printf("!!! causion: non-uniform range resolution is not supported !!! %d, %d\n", tmp_range_res, hd->range_res);
        goto READERR;
      }
      hd->range_res = tmp_range_res; //overwrite

      /*--- get range num and angle data---*/
      sp_rnum = char2uint16(bufelhd + 376); //unsigned
      lp_rnum = char2uint16(bufelhd + 378); //unsigned
      rnum = sp_rnum + lp_rnum;
      /* hd->range_num = rnum; */
      start_obs_az = char2uint16(bufelhd + 380) * ideg2deg; //unsigned
      start_obs_el = char2int16(bufelhd + 382) * ideg_tc2deg;
      end_obs_az   = char2uint16(bufelhd + 384) * ideg2deg; //unsigned
      end_obs_el   = char2int16(bufelhd + 386) * ideg_tc2deg;

      if(end_obs_az < start_obs_az){
        end_obs_az += 360.0;
      } 
      if(start_obs_el > 180.0){
        start_obs_el -= 360.0;
      }

      /*--- store AZ and EL ---*/
      az[k][j] = 0.5 * (start_obs_az + end_obs_az);
      el[k][j] = 0.5 * (start_obs_el + end_obs_el);

      if (opt_verbose >= 3){
        printf("=== short/long rnum=%d %d (%d)",sp_rnum, lp_rnum, rnum);
        printf("  start/end AZ=%6.2f-%6.2f (%6.2f) ",start_obs_az, end_obs_az, az[k][j]);
        printf("  start/end EL=%5.2f-%5.2f (%5.2f)\n",start_obs_el, end_obs_el, el[k][j]);
      }

      /*----- read data block -----*/
      irb = rnum * (mesh_size / 8);
      bufdat = buf + buf_offset;
      buf_offset += irb;

      /*--- store raw data ---*/
      for(i = 0; i < rnum; i++)
        rtdat[k][j][i] = mesh_offset + mesh_lsb * (float)(char2uint16(bufdat + 2 * i));

    } /* end of k-loop */
  } /* end of j-loop */

  goto READEND;

 READERR:
  return -8;

 READEND:
  return 0;
}
