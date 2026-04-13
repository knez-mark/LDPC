#include "minsum.h"
#include <string.h>

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

void reset_decode_state (lifting_struct* R_mj, uint16_t size) {
    memset (R_mj, 0, size*sizeof(lifting_struct));
}

static float sign (float in) {
    if (in > 0) {
        return 1;
    }
    else if (in < 0) {
        return -1;
    }
    else {
        return 0;
    }
}

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

void layered_normalized_minsum (ldpc_decoder_t* ldpc, ldpc_quantized_t * Lq, uint16_t len, uint8_t lifting_size) {

    quasi_cyclic_matrix_t* Hm = ldpc->Hm;
    quasi_cyclic_matrix_t* Hp = ldpc->Hp;
    #if !USE_SUM_PRODUCT
    float alpha = ldpc->cfg.alpha;
    #else
    float alpha = 0.75f;
    #endif

    uint16_t columnOffset [MAX_ROW_WEIGHT]= {0};
	uint16_t circular_shift [MAX_ROW_WEIGHT]= {0};

	for (int m = 0; m < Hm->rows; m++) {
        lifting_struct* R_mj = &ldpc->R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]];

        uint8_t rowWeight = 0;
        const uint8_t *col = &Hm->columnIndexMap[Hm->rowOffset[m]];
        const uint8_t *base_graph = &Hm->base_graph[Hm->rowOffset[m]];

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
            ldpc_quantized_t Lq_mj [MAX_ROW_WEIGHT];
            for (int i = 0; i < rowWeight; i++) {
                uint16_t shift = circular_shift[i] + j;
                if (shift >= lifting_size) shift -= lifting_size;

                Lq_mj[i] = CLAMP (Lq[columnOffset[i]+shift] - R_mj[i].element[j]);
            }

            ldpc_quantized_t v1, v2;
            uint8_t i1 = 0;
            int8_t sign_prod;
            min2_sign (Lq_mj, rowWeight, &v1, &v2, &i1, &sign_prod);

            for (int i = 0; i < rowWeight; i++) {
                int8_t s = (Lq_mj[i] < 0) ? -1 : 1;

                ldpc_quantized_t val = (i == i1) ? v2 : v1;

                R_mj[i].element[j] = CLAMP (alpha * sign_prod * s * val);

                uint16_t shift = circular_shift[i] + j;
                if (shift >= lifting_size) shift -= lifting_size;

                Lq[columnOffset[i]+shift] = CLAMP (Lq_mj[i] + R_mj[i].element[j]);
            }
        }
    }
}

void layered_sum_product (ldpc_decoder_t* ldpc, float * Lq, uint16_t len, uint8_t lifting_size) {

    quasi_cyclic_matrix_t* Hm = ldpc->Hm;
    quasi_cyclic_matrix_t* Hp = ldpc->Hp;

    uint16_t columnOffset [MAX_ROW_WEIGHT]= {0};
	uint16_t circular_shift [MAX_ROW_WEIGHT]= {0};

	for (int m = 0; m < Hm->rows; m++) {
        lifting_struct* R_mj = &ldpc->R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]];

        uint8_t rowWeight = 0;
        const uint8_t *col = &Hm->columnIndexMap[Hm->rowOffset[m]];
        const uint8_t *base_graph = &Hm->base_graph[Hm->rowOffset[m]];

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
            float Lq_mj [MAX_ROW_WEIGHT];
            for (int i = 0; i < rowWeight; i++) {
                uint16_t shift = circular_shift[i] + j;
                if (shift >= lifting_size) shift -= lifting_size;

                Lq_mj[i] = Lq[columnOffset[i]+shift] - R_mj[i].element[j];
            }

            float prefix[MAX_ROW_WEIGHT];
            float suffix[MAX_ROW_WEIGHT];

            prefix[0] = Lq_mj[0];
            for (int i = 1; i < rowWeight; i++) {
                prefix[i] = boxplus(prefix[i-1], Lq_mj[i]);
            }
            suffix[rowWeight-1] = Lq_mj[rowWeight-1];
            for (int i = rowWeight-2; i >= 0; i--) {
                suffix[i] = boxplus (Lq_mj[i], suffix[i+1]);
            }

            for (int i = 0; i < rowWeight; i++) {
                if (i == 0) {
                    R_mj[i].element[j] = suffix[1];
                }
                else if (i == rowWeight-1) {
                    R_mj[i].element[j] = prefix[rowWeight-2];
                }
                else {
                    R_mj[i].element[j] = boxplus(prefix[i-1], suffix[i+1]);
                }
            }

            for (int i = 0; i < rowWeight; i++) {
                uint16_t shift = circular_shift[i] + j;
                if (shift >= lifting_size) shift -= lifting_size;

            	Lq[columnOffset[i]+shift] = Lq_mj[i] + R_mj[i].element[j];
            }
        }
    }
}