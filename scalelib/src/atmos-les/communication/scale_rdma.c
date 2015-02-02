/*
 * rdma.c
 *
 *  Created on: 2012/01/16
 *      Author: ohno
 */

#include <stdlib.h>
#include "rdma.h"
#include <mpi.h>
#include <mpi-ext.h>

/** Private parameters & variables **/
#define WEST			0
#define NORTH			1
#define EAST			2
#define SOUTH			3
#define BEARING_CNT	4

int32_t COMM_vsize_max ;

int32_t	IA, JA, KA ;
int32_t	IHALO, JHALO;
int32_t	IS, IE, JS, JE ;
#define	offset(Z,X,Y)	(sizeof(var_t)*((Z)+(X)*(KA)+(Y)*(KA)*(IA)))

size_t		datasize_NS, datasize_WE;

int32_t	RANK_W, RANK_N, RANK_E, RANK_S;

#define	RDMA_TAG_NUM_MAX	15
#define	RDMA_TAG		0
#define	RDMA_TAG_TAIL		1

int32_t	memid_cnt;

int32_t	*memid ;
uint64_t	*lvar ;
uint64_t	**rvar ;


typedef uint64_t	sf_t ;
volatile sf_t	*status_flag;
int32_t		memid_sf;
uint64_t	local_sf;
uint64_t	*remote_sf;
#define LOCAL_RECV_READY	0
#define REMOTE_RECV_READY	1
#define LOCAL_PUT_DONE		2
#define REMOTE_PUT_DONE		3
#define LOCAL_RECV_DONE		3
#define FLAG_CNT				4
#define sidx(DIR,FLAG)		((DIR)+(BEARING_CNT)*(FLAG))
#define soffset(DIR,FLAG)	(sizeof(sf_t)*sidx(DIR,FLAG))

int64_t		*sending; 
#define FALSE	0
#define TRUE	1


volatile sf_t	rdma_put_id ;
#define	RDMA_PUT_ID_MAX	0xFFFFFFFF


#define FJMPI_RDMA_PUT_FLAGS_TO_NORTH	(FJMPI_RDMA_LOCAL_NIC0 | FJMPI_RDMA_REMOTE_NIC2 | FJMPI_RDMA_PATH0 )
#define FJMPI_RDMA_PUT_FLAGS_TO_SOUTH	(FJMPI_RDMA_LOCAL_NIC2 | FJMPI_RDMA_REMOTE_NIC0 | FJMPI_RDMA_PATH0 )
#define FJMPI_RDMA_PUT_FLAGS_TO_EAST	(FJMPI_RDMA_LOCAL_NIC3 | FJMPI_RDMA_REMOTE_NIC1 | FJMPI_RDMA_PATH0 )
#define FJMPI_RDMA_PUT_FLAGS_TO_WEST	(FJMPI_RDMA_LOCAL_NIC1 | FJMPI_RDMA_REMOTE_NIC3 | FJMPI_RDMA_PATH0 )

#include <stdio.h>

void rdma_setup_(
		const int32_t *COMM_vsize_max_in,
		const int32_t *IA_in,
		const int32_t *JA_in,
		const int32_t *KA_in,
		const int32_t *IHALO_in,
		const int32_t *JHALO_in,
		const int32_t *IS_in,
		const int32_t *IE_in,
		const int32_t *JS_in,
		const int32_t *JE_in,
		const int32_t *RANK_W_in,
		const int32_t *RANK_N_in,
		const int32_t *RANK_E_in,
		const int32_t *RANK_S_in)
{
	int v;

	COMM_vsize_max = *COMM_vsize_max_in ;
	IA = *IA_in ;
	JA = *JA_in ;
	KA = *KA_in ;
	IHALO = *IHALO_in ;
	JHALO = *JHALO_in ;
	IS = *IS_in ;
	IE = *IE_in ;
	JS = *JS_in ;
	JE = *JE_in ;
	RANK_W = *RANK_W_in ;
	RANK_N = *RANK_N_in ;
	RANK_E = *RANK_E_in ;
	RANK_S = *RANK_S_in ;

	datasize_NS = sizeof(var_t) * (IE-IS+1) * KA ;
	datasize_WE = sizeof(var_t) * IHALO * KA ;

	memid_cnt = 0;

	memid = (int32_t *) malloc(sizeof(int32_t) * COMM_vsize_max);
	lvar = (uint64_t *) malloc(sizeof(uint64_t ) * COMM_vsize_max) ;
	rvar = (uint64_t **) malloc(sizeof(uint64_t *) * COMM_vsize_max) ;
	for(v=0; v<COMM_vsize_max; v++) rvar[v] = (uint64_t *) malloc(sizeof(uint64_t) * BEARING_CNT) ;

	status_flag = (sf_t *) calloc(BEARING_CNT*FLAG_CNT, sizeof(sf_t)) ;
	remote_sf = (uint64_t *) malloc(sizeof(uint64_t) * BEARING_CNT) ;

	sending  = (int64_t *) calloc(BEARING_CNT, sizeof(int64_t)) ;

	rdma_put_id = 0 ;

	FJMPI_Rdma_init() ;

	memid_sf = memid_cnt ;
	local_sf = FJMPI_Rdma_reg_mem(memid_cnt, status_flag, sizeof(sf_t)*BEARING_CNT*FLAG_CNT) ;
	memid_cnt++ ;

	MPI_Barrier(MPI_COMM_WORLD) ;

	if( RANK_W != MPI_PROC_NULL )
		remote_sf[WEST]  = FJMPI_Rdma_get_remote_addr(RANK_W, memid_sf) ;
	else	
		remote_sf[WEST]  = 0 ;
	
	if( RANK_N != MPI_PROC_NULL )
		remote_sf[NORTH] = FJMPI_Rdma_get_remote_addr(RANK_N, memid_sf) ;
	else
		remote_sf[NORTH] = 0 ;

	if( RANK_E != MPI_PROC_NULL )
		remote_sf[EAST]  = FJMPI_Rdma_get_remote_addr(RANK_E, memid_sf) ;
	else
		remote_sf[EAST]  = 0 ;

	if( RANK_S != MPI_PROC_NULL )
		remote_sf[SOUTH] = FJMPI_Rdma_get_remote_addr(RANK_S, memid_sf) ;
	else
		remote_sf[SOUTH] = 0 ;

}


void set_rdma_variable_(
		const var_t		*var,
		const int32_t		*vid )
{

	memid[*vid] = memid_cnt ;
	lvar[*vid] = FJMPI_Rdma_reg_mem(memid_cnt, var, sizeof(var_t)*IA*JA*KA) ;
	memid_cnt++ ;


	MPI_Barrier(MPI_COMM_WORLD) ;

	if( RANK_W != MPI_PROC_NULL ) rvar[*vid][WEST]  = FJMPI_Rdma_get_remote_addr(RANK_W, memid[*vid]) ;
	if( RANK_N != MPI_PROC_NULL ) rvar[*vid][NORTH] = FJMPI_Rdma_get_remote_addr(RANK_N, memid[*vid]) ;
	if( RANK_E != MPI_PROC_NULL ) rvar[*vid][EAST]  = FJMPI_Rdma_get_remote_addr(RANK_E, memid[*vid]) ;
	if( RANK_S != MPI_PROC_NULL ) rvar[*vid][SOUTH] = FJMPI_Rdma_get_remote_addr(RANK_S, memid[*vid]) ;

}


void rdma_put_(const int32_t *vid, const int32_t *num)
{
	struct FJMPI_Rdma_cq cq ;
	int	j , v;

	rdma_put_id = (rdma_put_id % RDMA_PUT_ID_MAX) + 1 ;

	/* set status_flag(recv) */
	status_flag[sidx(WEST ,LOCAL_RECV_READY)] = rdma_put_id ;
	status_flag[sidx(NORTH,LOCAL_RECV_READY)] = rdma_put_id ;
	status_flag[sidx(EAST ,LOCAL_RECV_READY)] = rdma_put_id ;
	status_flag[sidx(SOUTH,LOCAL_RECV_READY)] = rdma_put_id ;

	if( RANK_S != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_S, RDMA_TAG,
				remote_sf[SOUTH]+soffset(NORTH,REMOTE_RECV_READY),
				local_sf+soffset(SOUTH,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_SOUTH | FJMPI_RDMA_STRONG_ORDER ) ;

	if( RANK_N != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_N, RDMA_TAG,
				remote_sf[NORTH]+soffset(SOUTH,REMOTE_RECV_READY),
				local_sf+soffset(NORTH,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_NORTH | FJMPI_RDMA_STRONG_ORDER ) ;

	if( RANK_E != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_E, RDMA_TAG,
				remote_sf[EAST]+soffset(WEST,REMOTE_RECV_READY),
				local_sf+soffset(EAST,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_EAST | FJMPI_RDMA_STRONG_ORDER ) ;

	if( RANK_W != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_W, RDMA_TAG,
				remote_sf[WEST]+soffset(EAST,REMOTE_RECV_READY),
				local_sf+soffset(WEST,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_WEST | FJMPI_RDMA_STRONG_ORDER ) ;

	/* send data */
	do {
		// to north
		if( status_flag[sidx(NORTH,LOCAL_PUT_DONE)] != rdma_put_id && RANK_N == MPI_PROC_NULL )
		{
			status_flag[sidx(NORTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			status_flag[sidx(NORTH,LOCAL_RECV_DONE)] = rdma_put_id ;
			sending[NORTH] = FALSE ;
		}
		
		if( status_flag[sidx(NORTH,LOCAL_PUT_DONE)] != rdma_put_id &&
				status_flag[sidx(NORTH,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[NORTH] = TRUE ;
			
			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_N, RDMA_TAG,
							rvar[*vid+v][NORTH]+offset(0,IS-1,j+JE),
							lvar[*vid+v]+offset(0,IS-1,j+JS-1),
							datasize_NS,
							FJMPI_RDMA_PUT_FLAGS_TO_NORTH ) ;
				}
			}

			// set status_flag(send)
			status_flag[sidx(NORTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			FJMPI_Rdma_put(RANK_N, RDMA_TAG,
					remote_sf[NORTH]+soffset(SOUTH,REMOTE_PUT_DONE),
					local_sf+soffset(NORTH,LOCAL_PUT_DONE),
					sizeof(sf_t),
					FJMPI_RDMA_PUT_FLAGS_TO_NORTH | FJMPI_RDMA_STRONG_ORDER ) ;
		}

		// to south
		if( status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] != rdma_put_id && RANK_S == MPI_PROC_NULL )
		{
			status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			status_flag[sidx(SOUTH,LOCAL_RECV_DONE)] = rdma_put_id ;
			sending[SOUTH] = FALSE;
		}
	
		if( status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] != rdma_put_id  &&
				status_flag[sidx(SOUTH,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[SOUTH] = TRUE;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_S, RDMA_TAG,
							rvar[*vid+v][SOUTH]+offset(0,IS-1,j+JS-JHALO-1),
							lvar[*vid+v]+offset(0,IS-1,j+JE-JHALO),
							datasize_NS,
							FJMPI_RDMA_PUT_FLAGS_TO_SOUTH ) ;
				}
			}

			// set status_flag(send)
			status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			FJMPI_Rdma_put(RANK_S, RDMA_TAG,
					remote_sf[SOUTH]+soffset(NORTH,REMOTE_PUT_DONE),
					local_sf+soffset(SOUTH,LOCAL_PUT_DONE),
					sizeof(sf_t),
					FJMPI_RDMA_PUT_FLAGS_TO_SOUTH | FJMPI_RDMA_STRONG_ORDER ) ;
		}

		// to west
		if( status_flag[sidx(WEST,LOCAL_PUT_DONE)] != rdma_put_id && RANK_W == MPI_PROC_NULL )
		{
			status_flag[sidx(WEST,LOCAL_PUT_DONE)] = rdma_put_id ;
			status_flag[sidx(WEST,LOCAL_RECV_DONE)] = rdma_put_id ;
			sending[WEST] = FALSE;
		}

		if( status_flag[sidx(WEST,LOCAL_PUT_DONE)] != rdma_put_id  &&
				status_flag[sidx(WEST,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[WEST] = TRUE;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=JS-1; j<JE; j++)
				{
					FJMPI_Rdma_put(RANK_W, RDMA_TAG,
							rvar[*vid+v][WEST]+offset(0,IE,j),
							lvar[*vid+v]+offset(0,IS-1,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_WEST ) ;
				}
			}
	
			// set status_flag(send)
			status_flag[sidx(WEST,LOCAL_PUT_DONE)] = rdma_put_id ;
			FJMPI_Rdma_put(RANK_W, RDMA_TAG,
					remote_sf[WEST]+soffset(EAST,REMOTE_PUT_DONE),
					local_sf+soffset(WEST,LOCAL_PUT_DONE),
					sizeof(sf_t),
					FJMPI_RDMA_PUT_FLAGS_TO_WEST | FJMPI_RDMA_STRONG_ORDER ) ;
		}

		// to east
		if( status_flag[sidx(EAST,LOCAL_PUT_DONE)] != rdma_put_id && RANK_E == MPI_PROC_NULL )
		{
			status_flag[sidx(EAST,LOCAL_PUT_DONE)] = rdma_put_id ;
			status_flag[sidx(EAST,LOCAL_RECV_DONE)] = rdma_put_id ;
			sending[EAST] = FALSE;
		}

		if( status_flag[sidx(EAST,LOCAL_PUT_DONE)] != rdma_put_id &&
				status_flag[sidx(EAST,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[TRUE] = FALSE;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=JS-1; j<JE; j++)
				{
					FJMPI_Rdma_put(RANK_E, RDMA_TAG,
							rvar[*vid+v][EAST]+offset(0,IS-IHALO-1,j),
							lvar[*vid+v]+offset(0,IE-IHALO,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_EAST ) ;
				}
			}

			// set status_flag(send)
			status_flag[sidx(EAST,LOCAL_PUT_DONE)] = rdma_put_id ;
			FJMPI_Rdma_put(RANK_E, RDMA_TAG,
					remote_sf[EAST]+soffset(WEST,REMOTE_PUT_DONE),
					local_sf+soffset(EAST,LOCAL_PUT_DONE),
					sizeof(sf_t),
					FJMPI_RDMA_PUT_FLAGS_TO_EAST | FJMPI_RDMA_STRONG_ORDER ) ;
		}
	} while( status_flag[sidx(WEST, LOCAL_PUT_DONE)] != rdma_put_id ||
			status_flag[sidx(NORTH,LOCAL_PUT_DONE)] != rdma_put_id ||
			status_flag[sidx(EAST, LOCAL_PUT_DONE)] != rdma_put_id ||
			status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] != rdma_put_id ) ;

	/* completion check (put) */
	do {
		// to west
		if( RANK_W != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC1, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[WEST] = FALSE ;
			}
		}
		// to east
		if( RANK_E != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC3, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[EAST] = FALSE ;
			}
		}
		// to north
		if( RANK_N != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC0, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[NORTH] = FALSE ;
			}
		}
		// to south
		if( RANK_S != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC2, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[SOUTH] = FALSE ;
			}
		}
	 } while(sending[WEST]  || sending[NORTH] || sending[EAST]  || sending[SOUTH] ) ;

	/* completion check (recv) */
	while( status_flag[sidx(WEST,LOCAL_RECV_DONE)] != rdma_put_id ||
			status_flag[sidx(NORTH,LOCAL_RECV_DONE)] != rdma_put_id ||
			status_flag[sidx(EAST,LOCAL_RECV_DONE)]  != rdma_put_id ||
			status_flag[sidx(SOUTH,LOCAL_RECV_DONE)] != rdma_put_id ) ;

}

void rdma_put8_(const int32_t *vid, const int32_t *num)
{
	struct FJMPI_Rdma_cq cq ;
	int	j , v;

	rdma_put_id = (rdma_put_id % RDMA_PUT_ID_MAX) + 1 ;

	/* set status_flag(recv) */
	status_flag[sidx(WEST ,LOCAL_RECV_READY)] = rdma_put_id ;
	status_flag[sidx(NORTH,LOCAL_RECV_READY)] = rdma_put_id ;
	status_flag[sidx(EAST ,LOCAL_RECV_READY)] = rdma_put_id ;
	status_flag[sidx(SOUTH,LOCAL_RECV_READY)] = rdma_put_id ;

	if( RANK_S != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_S, RDMA_TAG,
				remote_sf[SOUTH]+soffset(NORTH,REMOTE_RECV_READY),
				local_sf+soffset(SOUTH,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_SOUTH | FJMPI_RDMA_STRONG_ORDER ) ;

	if( RANK_N != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_N, RDMA_TAG,
				remote_sf[NORTH]+soffset(SOUTH,REMOTE_RECV_READY),
				local_sf+soffset(NORTH,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_NORTH | FJMPI_RDMA_STRONG_ORDER ) ;

	if( RANK_E != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_E, RDMA_TAG,
				remote_sf[EAST]+soffset(WEST,REMOTE_RECV_READY),
				local_sf+soffset(EAST,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_EAST | FJMPI_RDMA_STRONG_ORDER ) ;

	if( RANK_W != MPI_PROC_NULL )
		FJMPI_Rdma_put(RANK_W, RDMA_TAG,
				remote_sf[WEST]+soffset(EAST,REMOTE_RECV_READY),
				local_sf+soffset(WEST,LOCAL_RECV_READY),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_WEST | FJMPI_RDMA_STRONG_ORDER ) ;

	/* send data */
	do {
		// to north
		if( status_flag[sidx(NORTH,LOCAL_PUT_DONE)] != rdma_put_id && RANK_N == MPI_PROC_NULL )
		{
			status_flag[sidx(NORTH,LOCAL_PUT_DONE)] = rdma_put_id ;	
			sending[NORTH] =  FALSE ;
		}

		if( status_flag[sidx(NORTH,LOCAL_PUT_DONE)] != rdma_put_id &&
				status_flag[sidx(NORTH,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[NORTH] =  TRUE ;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_N, RDMA_TAG,
							rvar[*vid+v][NORTH]+offset(0,IS-1,j+JE),
							lvar[*vid+v]+offset(0,IS-1,j+JS-1),
							datasize_NS,
							FJMPI_RDMA_PUT_FLAGS_TO_NORTH ) ;

				}
			}

			// set status_flag(send)
			status_flag[sidx(NORTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			FJMPI_Rdma_put(RANK_N, RDMA_TAG_TAIL,
					remote_sf[NORTH]+soffset(SOUTH,REMOTE_PUT_DONE),
					local_sf+soffset(NORTH,LOCAL_PUT_DONE),
					sizeof(sf_t),
					FJMPI_RDMA_PUT_FLAGS_TO_NORTH | FJMPI_RDMA_STRONG_ORDER ) ;	
		}

		// to south	
		if( status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] != rdma_put_id && RANK_S == MPI_PROC_NULL )
		{
			status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			sending[SOUTH] = FALSE;
		}

		if( status_flag[sidx(SOUTH,LOCAL_PUT_DONE)]  != rdma_put_id &&
				status_flag[sidx(SOUTH,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[SOUTH] = TRUE;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_S, RDMA_TAG,
							rvar[*vid+v][SOUTH]+offset(0,IS-1,j+JS-JHALO-1),
							lvar[*vid+v]+offset(0,IS-1,j+JE-JHALO),
							datasize_NS,
							FJMPI_RDMA_PUT_FLAGS_TO_SOUTH ) ;
				}
			}

			// set status_flag(send)
			status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] = rdma_put_id ;
			FJMPI_Rdma_put(RANK_S, RDMA_TAG_TAIL,
					remote_sf[SOUTH]+soffset(NORTH,REMOTE_PUT_DONE),
					local_sf+soffset(SOUTH,LOCAL_PUT_DONE),
					sizeof(sf_t),
					FJMPI_RDMA_PUT_FLAGS_TO_SOUTH | FJMPI_RDMA_STRONG_ORDER ) ;
		}

		// to west
		if( status_flag[sidx(WEST,LOCAL_PUT_DONE)] != rdma_put_id && RANK_W == MPI_PROC_NULL )
		{
			status_flag[sidx(WEST,LOCAL_PUT_DONE)] = rdma_put_id ;
			sending[WEST] = FALSE ;
		}
		
		if( status_flag[sidx(WEST,LOCAL_PUT_DONE)] != rdma_put_id &&
				status_flag[sidx(WEST,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[WEST] = TRUE ;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=JS-1; j<JE; j++)
				{
					FJMPI_Rdma_put(RANK_W, RDMA_TAG,
							rvar[*vid+v][WEST]+offset(0,IE,j),
							lvar[*vid+v]+offset(0,IS-1,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_WEST ) ;
				}
			}

			// set status_flag(send)
			status_flag[sidx(WEST,LOCAL_PUT_DONE)] = rdma_put_id ;
		}

		// to east
		if( status_flag[sidx(EAST,LOCAL_PUT_DONE)] != rdma_put_id && RANK_E == MPI_PROC_NULL )
		{
			status_flag[sidx(EAST,LOCAL_PUT_DONE)] = rdma_put_id ;
			sending[EAST] = FALSE ;
		}

		if( status_flag[sidx(EAST,LOCAL_PUT_DONE)] != rdma_put_id &&
				status_flag[sidx(EAST,REMOTE_RECV_READY)] == rdma_put_id )
		{
			sending[EAST] = TRUE ;

			// put_data
			for(v=0; v<*num; v++)
			{
				for(j=JS-1; j<JE; j++)
				{
					FJMPI_Rdma_put(RANK_E, RDMA_TAG,
							rvar[*vid+v][EAST]+offset(0,IS-IHALO-1,j),
							lvar[*vid+v]+offset(0,IE-IHALO,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_EAST ) ;
				}
			}

			// set status_flag(send)
			status_flag[sidx(EAST,LOCAL_PUT_DONE)] = rdma_put_id ;
		}

	} while( status_flag[sidx(WEST,LOCAL_PUT_DONE)]  != rdma_put_id ||
			status_flag[sidx(NORTH,LOCAL_PUT_DONE)] != rdma_put_id ||
			status_flag[sidx(EAST,LOCAL_PUT_DONE)]  != rdma_put_id ||
			status_flag[sidx(SOUTH,LOCAL_PUT_DONE)] != rdma_put_id ) ;	

	/* completion check (recv from NORTH and SOUTH) */
	if( RANK_N == MPI_PROC_NULL ) status_flag[sidx(NORTH, LOCAL_RECV_DONE)] = rdma_put_id ;	
	if( RANK_S == MPI_PROC_NULL ) status_flag[sidx(SOUTH, LOCAL_RECV_DONE)] = rdma_put_id ;	

	while( status_flag[sidx(NORTH,LOCAL_RECV_DONE)] != rdma_put_id ||
			status_flag[sidx(SOUTH,LOCAL_RECV_DONE)] != rdma_put_id ) ;

	/* send data from North and South to West and East */
	for(v=0; v<*num; v++)
	{
		// data form north	
		if( RANK_N != MPI_PROC_NULL ) {
			if( RANK_W != MPI_PROC_NULL && RANK_E != MPI_PROC_NULL ) {
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_W, RDMA_TAG,
							rvar[*vid+v][WEST]+offset(0,IE,j),
							lvar[*vid+v]+offset(0,IS-1,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_WEST ) ;
					FJMPI_Rdma_put(RANK_E, RDMA_TAG,
							rvar[*vid+v][EAST]+offset(0,IS-IHALO-1,j),
							lvar[*vid+v]+offset(0,IE-IHALO,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_EAST ) ;
				}
			}
			else if( RANK_W != MPI_PROC_NULL ) {
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_W, RDMA_TAG,
							rvar[*vid+v][WEST]+offset(0,IE,j),
							lvar[*vid+v]+offset(0,IS-1,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_WEST ) ;
				}
			}
			else if( RANK_E != MPI_PROC_NULL ) {
				for(j=0; j<JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_E, RDMA_TAG,
							rvar[*vid+v][EAST]+offset(0,IS-IHALO-1,j),
							lvar[*vid+v]+offset(0,IE-IHALO,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_EAST ) ;
				}
			}
		}

		// data from south
		if( RANK_S != MPI_PROC_NULL ) {
			if( RANK_W != MPI_PROC_NULL && RANK_E != MPI_PROC_NULL ) {
				for(j=JE; j<JE+JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_W, RDMA_TAG,
							rvar[*vid+v][WEST]+offset(0,IE,j),
							lvar[*vid+v]+offset(0,IS-1,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_WEST ) ;
					FJMPI_Rdma_put(RANK_E, RDMA_TAG,
							rvar[*vid+v][EAST]+offset(0,IS-IHALO-1,j),
							lvar[*vid+v]+offset(0,IE-IHALO,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_EAST ) ;
				}
			}
			else if( RANK_W != MPI_PROC_NULL ) {
				for(j=JE; j<JE+JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_W, RDMA_TAG,
							rvar[*vid+v][WEST]+offset(0,IE,j),
							lvar[*vid+v]+offset(0,IS-1,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_WEST ) ;
				}
			}
			else if( RANK_E != MPI_PROC_NULL ) {
				for(j=JE; j<JE+JHALO; j++)
				{
					FJMPI_Rdma_put(RANK_E, RDMA_TAG,
							rvar[*vid+v][EAST]+offset(0,IS-IHALO-1,j),
							lvar[*vid+v]+offset(0,IE-IHALO,j),
							datasize_WE,
							FJMPI_RDMA_PUT_FLAGS_TO_EAST ) ;
				}
			}
		}
	}
	
	if( RANK_W != MPI_PROC_NULL ) {
		FJMPI_Rdma_put(RANK_W, RDMA_TAG_TAIL,
				remote_sf[WEST]+soffset(EAST,REMOTE_PUT_DONE),
				local_sf+soffset(WEST,LOCAL_PUT_DONE),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_WEST | FJMPI_RDMA_STRONG_ORDER ) ;
	}

	if( RANK_E != MPI_PROC_NULL ) {
		FJMPI_Rdma_put(RANK_E, RDMA_TAG_TAIL,
				remote_sf[EAST]+soffset(WEST,REMOTE_PUT_DONE),
				local_sf+soffset(EAST,LOCAL_PUT_DONE),
				sizeof(sf_t),
				FJMPI_RDMA_PUT_FLAGS_TO_EAST | FJMPI_RDMA_STRONG_ORDER ) ;
	}
	
	/* completion check (put) */
	do {
		// to west
		if( RANK_W != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC1, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[WEST] = FALSE ;
			}
		}
		// to east
		if( RANK_E != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC3, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[EAST] = FALSE ;
			}
		}
		// to north
		if( RANK_N != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC0, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[NORTH] = FALSE ;
			}
		}
		// to south
		if( RANK_S != MPI_PROC_NULL ) {
			while(FJMPI_Rdma_poll_cq(FJMPI_RDMA_NIC2, &cq) == FJMPI_RDMA_NOTICE ) {
				if(cq.tag == RDMA_TAG_TAIL) sending[SOUTH] = FALSE ;
			}
		}
	} while(sending[WEST]  || sending[NORTH] || sending[EAST]  || sending[SOUTH] ) ;

	/* completion check (recv from WEST and EAST) */
	if( RANK_W == MPI_PROC_NULL ) status_flag[sidx(WEST, LOCAL_RECV_DONE)] = rdma_put_id  ;
	if( RANK_E == MPI_PROC_NULL ) status_flag[sidx(EAST, LOCAL_RECV_DONE)] = rdma_put_id  ;


	while( status_flag[sidx(WEST,LOCAL_RECV_DONE)] != rdma_put_id ||
			status_flag[sidx(EAST,LOCAL_RECV_DONE)] != rdma_put_id ) ;
}


