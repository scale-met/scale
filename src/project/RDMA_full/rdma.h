/*
 * rdma.h
 *
 *  Created on: 2012/01/16
 *      Author: ohno
 */

#ifndef RDMA_H_
#define RDMA_H_

#include <stdint.h>

typedef double real64_t;
typedef float  real32_t;

typedef real64_t	var_t ;

extern void rdma_setup_(
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
		const int32_t *RANK_S_in );

extern void set_rdma_variable_(
		const var_t		*var,
		const int32_t		*vid );

extern void rdma_put_(
		const int32_t	*vid,
		const int32_t *num );

extern void rdma_put8_(
		const int32_t	*vid,
		const int32_t *num );


#endif /* RDMA_H_ */
