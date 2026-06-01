#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "LDPC_impl.h"
#include "base_graph.h"
#if USE_SUM_PRODUCT
#include <float.h>
#include <math.h>
#endif
#include "matrix.h"

#define TYPE_MAX_T(T) _Generic((T)0, \
    int8_t: INT8_MAX, \
    int16_t: INT16_MAX, \
    int32_t: INT32_MAX, \
    int64_t: INT64_MAX, \
    uint8_t: UINT8_MAX, \
    uint16_t: UINT16_MAX, \
    uint32_t: UINT32_MAX, \
    uint64_t: UINT64_MAX \
)

static void reset_decode_state (ldpc_decoder_t* ldpc, void *Lq_void) {
    memset (ldpc->R_mj, 0, sizeof(ldpc->R_mj));
    memset (ldpc->Lq, 0, sizeof (ldpc->Lq));

    uint16_t Z = ldpc->cfg.lifting_size;
    uint16_t msg_len = ldpc->cfg.msg_len;
    uint16_t code_len = ldpc->cfg.target_code_len;
    uint16_t parity_len = code_len - msg_len;

    //uint16_t actual_msg_len = ((msg_len + (Z-1))/Z)*Z;
    uint16_t actual_msg_len = ldpc->msg_size*Z;
    #if USE_SUM_PRODUCT
        memcpy (ldpc->Lq, Lq_void, msg_len*sizeof(float));
        for (int i = msg_len; i < actual_msg_len; i++) {
            ldpc->Lq [i] = 10.0f;
        }
        memcpy (ldpc->Lq + actual_msg_len, Lq_void + msg_len, parity_len*sizeof(float));
    #else
        memcpy (ldpc->Lq, Lq_void, msg_len*sizeof(ldpc_quantized_t));
        for (int i = msg_len; i < actual_msg_len; i++) {
            ldpc->Lq [i] = TYPE_MAX_T (ldpc_quantized_t);
        }
        memcpy (ldpc->Lq + actual_msg_len, Lq_void + msg_len, parity_len*sizeof(ldpc_quantized_t));
    #endif
}

#if USE_SUM_PRODUCT
static inline float boxplus(float a, float b)
{
    float sign = ((a >= 0) == (b >= 0)) ? 1.0f : -1.0f;

    float abs_a = fabsf(a);
    float abs_b = fabsf(b);

    float min_ab = (abs_a < abs_b) ? abs_a : abs_b;

    float term1 = log1pf(expf(-fabsf(a + b)));
    float term2 = log1pf(expf(-fabsf(a - b)));

    return sign * (min_ab + term1 - term2);
}
#endif

#define CLAMP(x) _Generic((ldpc_quantized_t)0, \
    int8_t:  clamp_int8, \
    int16_t: clamp_int16, \
    int32_t: clamp_int32 \
)(x)

static inline int8_t clamp_int8(int x)
{
    if (x > 127) return 127;
    if (x < -127) return -127;
    return (int8_t)x;
}

static inline int16_t clamp_int16(int x)
{
    if (x > 32767) return 32767;
    if (x < -32767) return -32767;
    return (int16_t)x;
}

static inline int32_t clamp_int32(int x)
{
    if (x > INT32_MAX) return INT32_MAX;
    if (x < -INT32_MAX) return -INT32_MAX;
    return (int32_t)x;
}

static void min2_sign (ldpc_quantized_t * arr, uint8_t size, ldpc_quantized_t * v1, ldpc_quantized_t * v2, uint8_t * i1, int8_t * sign) {
    
    ldpc_quantized_t min1 = TYPE_MAX_T(ldpc_quantized_t);
    ldpc_quantized_t min2 = TYPE_MAX_T(ldpc_quantized_t);
    uint8_t min_idx = 0;
    int8_t sign_prod = 1;

    for (int i = 0; i < size; i++) {
        ldpc_quantized_t x = arr[i];

        //Treat 0 as +1
        int8_t s = (x < 0) ? -1 : 1;
        sign_prod *= s;

        //Absolute value
        ldpc_quantized_t ax = (x < 0) ? -x : x;

        if (ax < min1) {
            min2 = min1;
            min1 = ax;
            min_idx = i;
        }
        else if (ax < min2) {
            min2 = ax;
        }
    }

    *v1 = min1;
    *v2 = min2;
    *i1 = min_idx;
    *sign = sign_prod;
}

static inline uint8_t is_less_than_zero (void * val, int idx) {
    #if USE_SUM_PRODUCT
        return ((float *)val)[idx] < 0;
    #else
        return ((ldpc_quantized_t *)val)[idx] < 0;
    #endif
}

static void get_hard_decision_codeword (void * Lq, uint16_t len, uint8_t * codeword) {
    uint8_t byte = 0;
    uint8_t bit = 0;
    uint16_t out_idx = 0;

    for (int i = 0; i < len; i++) {
        byte |= is_less_than_zero (Lq, i) << bit;

        bit++;
        if (bit == 8) {
            codeword[out_idx] = byte;
            out_idx++;
            
            byte = 0;
            bit = 0;
        }
    }

    if (bit != 0) {
        codeword[out_idx] = byte;
    }
}

static uint16_t check_syndrome (ldpc_decoder_t * ldpc, uint8_t * codeword) {

    uint16_t lifting_size = ldpc->cfg.lifting_size;
    uint16_t len = lifting_size*(ldpc->msg_size + ldpc->parity_size);

    uint8_t* syndrome = ldpc->syndrome;
    memset (syndrome, 0, sizeof (ldpc->syndrome));

    get_hard_decision_codeword (ldpc->Lq, len, codeword);

    uint8_t* codeword_aligned = ldpc->codeword;
    get_byte_aligned (codeword, codeword_aligned, len, lifting_size);
    uint16_t block_size = DIV_CEIL(lifting_size, NUM_BITS_PER_BYTE);
    uint8_t* temp = ldpc->temp;

    circular_matrix_multiply (&ldpc->Hm, codeword_aligned, syndrome, temp, lifting_size);
    circular_matrix_multiply (&ldpc->Hp, codeword_aligned + ldpc->Hm.cols*block_size, syndrome, temp, lifting_size);

    return find_vector_weight (syndrome, ldpc->Hm.rows*block_size*NUM_BITS_PER_BYTE); //Returns number of parity check equation failures
}

static void LDPC_decode_one_iter (ldpc_decoder_t* ldpc) {

    uint16_t lifting_size = ldpc->cfg.lifting_size;

    #if USE_SUM_PRODUCT
        float * Lq = ldpc->Lq;
        float * Lq_mj = ldpc->Lq_mj;

        float * prefix = ldpc->prefix; 
        float * suffix = ldpc->suffix;
    #else
        ldpc_quantized_t * Lq = ldpc->Lq;
        ldpc_quantized_t * Lq_mj = ldpc->Lq_mj;

        float alpha = ldpc->cfg.alpha;
    #endif

    quasi_cyclic_matrix_t* Hm = &ldpc->Hm;
    quasi_cyclic_matrix_t* Hp = &ldpc->Hp;

    uint16_t* columnOffset = ldpc->columnOffset;
	uint16_t* circular_shift = ldpc->circular_shift;

	for (int m = 0; m < Hm->rows; m++) {
        lifting_struct* R_mj = &ldpc->R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]];

        uint8_t rowWeight = 0;
        const uint8_t *col = &Hm->columnIndexMap[Hm->rowOffset[m]];
        const uint16_t *base_graph = &Hm->base_graph[Hm->rowOffset[m]];

		for (int k = 0; k < Hm->rowWeight[m]; k++) {
            if (col[k] >= Hm->cols) {
                break;
            }

		    columnOffset[rowWeight] = lifting_size*col[k];
		    circular_shift[rowWeight] = base_graph[k]%lifting_size;
        
            rowWeight++;
		}

        col = &Hp->columnIndexMap[Hp->rowOffset[m]];
        base_graph = &Hp->base_graph[Hp->rowOffset[m]];

		for (int k = 0; k < Hp->rowWeight[m]; k++) {
            if (col[k] >= Hp->cols) {
                break;
            }

		    columnOffset[rowWeight] = lifting_size*(Hm->cols + col[k]);
		    circular_shift[rowWeight] = base_graph[k]%lifting_size;
        
            rowWeight++;
		}

    	for (int j = 0; j < lifting_size; j++) {
            for (int i = 0; i < rowWeight; i++) {
                uint16_t shift = circular_shift[i] + j;
                if (shift >= lifting_size) shift -= lifting_size;

                #if USE_SUM_PRODUCT
                Lq_mj[i] = Lq[columnOffset[i]+shift] - R_mj[i].element[j];
                #else
                Lq_mj[i] = CLAMP (Lq[columnOffset[i]+shift] - R_mj[i].element[j]);
                #endif
            }

            #if USE_SUM_PRODUCT
            prefix[0] = Lq_mj[0];
            for (int i = 1; i < rowWeight; i++) {
                prefix[i] = boxplus(prefix[i-1], Lq_mj[i]);
            }
            suffix[rowWeight-1] = Lq_mj[rowWeight-1];
            for (int i = rowWeight-2; i >= 0; i--) {
                suffix[i] = boxplus (Lq_mj[i], suffix[i+1]);
            }
            #else
            ldpc_quantized_t v1, v2;
            uint8_t i1 = 0;
            int8_t sign_prod;
            min2_sign (Lq_mj, rowWeight, &v1, &v2, &i1, &sign_prod);
            #endif

            for (int i = 0; i < rowWeight; i++) {
                #if USE_SUM_PRODUCT
                if (i == 0) {
                    R_mj[i].element[j] = suffix[1];
                }
                else if (i == rowWeight-1) {
                    R_mj[i].element[j] = prefix[rowWeight-2];
                }
                else {
                    R_mj[i].element[j] = boxplus(prefix[i-1], suffix[i+1]);
                }
                #else
                int8_t s = (Lq_mj[i] < 0) ? -1 : 1;

                ldpc_quantized_t val = (i == i1) ? v2 : v1;

                R_mj[i].element[j] = CLAMP (alpha * sign_prod * s * val);
                #endif

                uint16_t shift = circular_shift[i] + j;
                if (shift >= lifting_size) shift -= lifting_size;

                #if USE_SUM_PRODUCT
            	Lq[columnOffset[i]+shift] = Lq_mj[i] + R_mj[i].element[j];
                #else
                Lq[columnOffset[i]+shift] = CLAMP (Lq_mj[i] + R_mj[i].element[j]);
                #endif
            }
        }
    }
}

uint16_t LDPC_decode (ldpc_decoder_t* ldpc, void * Lq, uint8_t * decoded, uint16_t * num_iters) {

    uint8_t lifting_size = ldpc->cfg.lifting_size;

    if (!is_valid_lifting_size(lifting_size)) {
        return 0;
    }

    uint16_t iters = 0;
    uint16_t parity_check_errors = 0;

    reset_decode_state (ldpc, Lq);

    for (iters = 0; iters < ldpc->cfg.max_iters; iters++) {

        parity_check_errors = check_syndrome (ldpc, decoded);
        if (parity_check_errors == 0) {
            break;
        }

        LDPC_decode_one_iter (ldpc);
    }
    
    *num_iters = iters;
    return parity_check_errors;
}