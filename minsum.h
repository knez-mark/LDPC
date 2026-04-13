#ifndef MINSUM_H_
#define MINSUM_H_

#include <stdint.h>
#include <float.h>
#include <math.h>
#include "matrix.h"
#include "LDPC_impl.h"

void reset_decode_state (lifting_struct* R_mj, uint16_t size);
void layered_normalized_minsum (ldpc_decoder_t* ldpc, ldpc_quantized_t * Lq, uint16_t len, uint8_t lifting_size);
void layered_sum_product (ldpc_decoder_t* ldpc, float * Lq, uint16_t len, uint8_t lifting_size);

#endif /* MINSUM_H_ */