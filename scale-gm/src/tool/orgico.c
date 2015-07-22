#include "orgico.h"
#include "fio_def.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int32_t orgico_glevel;
int32_t orgico_rlevel;
int64_t orgico_block1R1L;

static int32_t halo=1;
int32_t orgico_num_of_rgn;
char **orgico_fname_;

static int32_t system_endian;

static void cm_ednchg( void* avp_pointer,
		       int ai_size,
		       int ai_num)
{
  int ai_csize, ai_cnum;   
  char ac_buf[16];         
  char *acp_tmp;
  char *acp_local;         
  memset(ac_buf,'\0',sizeof(ac_buf));
  acp_tmp = avp_pointer;
  acp_local = avp_pointer;
  for ( ai_cnum=0; ai_cnum< ai_num; ai_cnum++) {
    memcpy(ac_buf, acp_local, ai_size);	
    for ( ai_csize=0; ai_csize< ai_size; ai_csize++){
      *acp_local = ac_buf[ai_size -1 -ai_csize];	
      acp_local++;
    }
    acp_tmp+=ai_size;	
  }
}




void orgico_setup( int32_t gl, int32_t rl )
{
  int32_t i=1;

  orgico_glevel=gl;
  orgico_rlevel=rl;
  orgico_num_of_rgn=pow(4,rl)*10;
  orgico_block1R1L
    = (pow(2,orgico_glevel-orgico_rlevel)+halo+halo)
    * (pow(2,orgico_glevel-orgico_rlevel)+halo+halo);

  if(*(char*)&i ){
    system_endian = FIO_LITTLE_ENDIAN;
  } else {
    system_endian = FIO_BIG_ENDIAN;
  }


}



void orgico_mk_fname( char *fname, char *base, int32_t i, int32_t y )
{
  switch (y){
  case 4 :
    sprintf(fname,"%s.rgn%04d",base,i);
    break;
  case 5 :
    sprintf(fname,"%s.rgn%05d",base,i);
    break;
  case 6 :
    sprintf(fname,"%s.rgn%06d",base,i);
    break;
  default :
    break;
  }

}

void orgico_setbasename( char *basename )
{
  int i;
  orgico_fname_ = (char **)malloc(sizeof(char **)*orgico_num_of_rgn);
  for(i=0;i<orgico_num_of_rgn;i++){
    orgico_fname_[i]=malloc(sizeof(char)*128);
    orgico_mk_fname( orgico_fname_[i],basename,i,5);
  }
}


void orgico_readgriddata( char *fname, double *data, int i )
{
  FILE *fp;
  int64_t offset[9];
  if( (fp=fopen(fname,"rb"))==NULL ){
	fprintf(stderr,"Can not open file : %s!\n",fname);
	exit(1);
  }

  /* GRD_x(XDIR) */
  offset[0]=(4+4+4)+4; 
  /* GRD_x(YDIR) */
  offset[1]=offset[0]+(8*orgico_block1R1L+4)+4;
  /* GRD_x(ZDIR) */
  offset[2]=offset[1]+(8*orgico_block1R1L+4)+4;
  /* GRD_xt(TI,XDIR) */
  offset[3]=offset[2]+(8*orgico_block1R1L+4)+4;
  /* GRD_xt(TJ,XDIR) */
  offset[4]=offset[3]+(8*orgico_block1R1L);
  /* GRD_xt(TI,YDIR) */
  offset[5]=offset[4]+(8*orgico_block1R1L+4)+4;
  /* GRD_xt(TJ,YDIR) */
  offset[6]=offset[5]+(8*orgico_block1R1L);
  /* GRD_xt(TI,ZDIR) */
  offset[7]=offset[6]+(8*orgico_block1R1L+4)+4;
  /* GRD_xt(TJ,ZDIR) */
  offset[8]=offset[7]+(8*orgico_block1R1L);
 

  fseek(fp,(long)offset[i],SEEK_SET);
  fread(data,sizeof(double),orgico_block1R1L,fp);
  /* orginal file is by BIG ENDIAN */
  if(system_endian==FIO_LITTLE_ENDIAN){
    cm_ednchg(data,sizeof(double),orgico_block1R1L);
  }

  fclose(fp);

  return;
}


void orgico_readdata_seq( char *fname, int did, int nol, void *data, int32_t esize )
{
  FILE *fp;
  int64_t offset;

  int i;

  if( (fp=fopen(fname,"rb"))==NULL ){
	fprintf(stderr,"Can not open file : %s!\n",fname);    
	exit(1);
  }

  offset=4;
  for(i=0;i<did;i++){
    offset = offset+(esize*orgico_block1R1L*nol+4)+4;
  }

  fseek(fp,(long)offset,SEEK_SET);
  fread(data,esize,orgico_block1R1L*nol,fp);
  /* orginal file is by BIG ENDIAN */
  if(system_endian==FIO_LITTLE_ENDIAN){
    cm_ednchg(data,esize,orgico_block1R1L*nol);
  }

  fclose(fp);

  return;
}

void orgico_readdata_dir( char *fname, int did, int nol, void *data, int32_t esize )
{
  FILE *fp;
  int64_t offset;

  int i;

  if( (fp=fopen(fname,"rb"))==NULL ){
	fprintf(stderr,"Can not open file : %s!\n",fname);    
	exit(1);
  }

  offset=0;
  for(i=0;i<did;i++){
    offset = offset+(esize*orgico_block1R1L*nol);
  }

  fseek(fp,(long)offset,SEEK_SET);
  fread(data,esize,orgico_block1R1L*nol,fp);
  /* orginal file is by BIG ENDIAN */
  if(system_endian==FIO_LITTLE_ENDIAN){
    cm_ednchg(data,esize,orgico_block1R1L*nol);
  }

  fclose(fp);

  return;
}


static void set_str( char *_str, char *str, int str_len )
{
  int i;

  strncpy(_str,str,str_len);

  if(_str[str_len-1]!=' '){
    _str[str_len]='\0';
  } else {
    for(i=str_len-1;i>=0;i--){
      if(_str[i]==' '){
	_str[i]='\0';
      } else {
	break;
      }
    }
  }
}



void orgico_setup_( int32_t *gl, int32_t *rl )
{
  orgico_setup(*gl,*rl);
}

void orgico_mk_fname_( char *fname,  char *base, int32_t *i, int32_t *y, 
		       int32_t fname_len, int32_t base_len )
{
  int n;
  char _base[FIO_HLONG];
  set_str( _base, base, base_len );

  orgico_mk_fname(fname,_base,*i,*y);
  for(n=strlen(fname);n<fname_len;n++){
    fname[n]=' ';
  }
  /*  printf("%d %c\n",strlen(fname),fname[strlen(fname)-1]); */

}

void orgico_readgriddata_( char *fname, double *data, int32_t *i, 
			   int32_t fname_len )
{
  char _fname[FIO_HLONG];
  set_str( _fname, fname, fname_len );
  orgico_readgriddata( _fname, data, *i );
}

void orgico_setbasename_( char *basename, 
			  int32_t basename_len )
{
  char _basename[FIO_HLONG];
  set_str( _basename, basename, basename_len );
  orgico_setbasename( _basename );

}

void orgico_readdata_seq_( char *fname, int *did, int *nol, void *data, int32_t *esize, 
			   int32_t fname_len )
{
  char _fname[FIO_HLONG];
  set_str( _fname, fname, fname_len );
  orgico_readdata_seq( _fname, *did, *nol, data, *esize );
}

void orgico_readdata_dir_( char *fname, int *did, int *nol, void *data, int32_t *esize, 
			   int32_t fname_len )
{
  char _fname[FIO_HLONG];
  set_str( _fname, fname, fname_len );
  orgico_readdata_dir( _fname, *did, *nol, data, *esize );
}
