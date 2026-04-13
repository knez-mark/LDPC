#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "LDPC_impl.h"
#include "matrix.h"
#include "minsum.h"

static inline uint8_t is_less_than_zero (void * val, int idx) {
    #if USE_SUM_PRODUCT
        return ((float *)val)[idx] < 0;
    #else
        return ((ldpc_quantized_t *)val)[idx] < 0;
    #endif
}

static void get_hard_decision_codeword (void * Lq, uint16_t len, uint8_t lifting_size, vector_t codeword) {
    int8_t bit_ptr = 0;
    uint8_t byte_ptr = 0;

    for (int i = 0; i < len; i++) {
        if (lifting_size == 32) {
            if (is_less_than_zero (Lq, i)) { 
                codeword.data32 [byte_ptr] |= 1 << bit_ptr;
            }
        }
        else if (lifting_size == 16) {
            if (is_less_than_zero (Lq, i)) { 
                codeword.data16 [byte_ptr] |= 1 << bit_ptr;
            }
        }
        else { //lifting size of 8
            if (is_less_than_zero (Lq, i)) { 
                codeword.data8 [byte_ptr] |= 1 << bit_ptr;
            }
        }

        bit_ptr++;
        if (bit_ptr == lifting_size) {
            bit_ptr = 0;
            byte_ptr++;
        }
    }
}

static uint16_t check_syndrome (void * Lq, uint16_t len, quasi_cyclic_matrix_t * Hm, quasi_cyclic_matrix_t * Hp, uint8_t * codeword) {

    uint8_t lifting_size = len/(Hm->cols + Hp->cols);

    uint8_t syndrome [MAX_CODE_LEN] = {0};
    memset (syndrome, 0, MAX_CODE_LEN);
    memset (codeword, 0, MAX_CODE_LEN);

    get_hard_decision_codeword (Lq, len, lifting_size, (vector_t) codeword);

    circular_matrix_multiply (Hm, (vector_t) codeword, (vector_t) syndrome, lifting_size);
    circular_matrix_multiply (Hp, (vector_t) (codeword + Hm->cols*lifting_size/EIGHT_BITS_PER_BYTE), (vector_t) syndrome, lifting_size);

    return find_vector_weight ((vector_t) syndrome, len); //Returns number of parity check equation failures
}

uint16_t LDPC_decode (ldpc_decoder_t* ldpc, void * Lq, uint16_t len, uint8_t * decoded, uint16_t * num_iters) {

    uint8_t lifting_size = len/(ldpc->cfg.msg_size + ldpc->cfg.parity_size);

    if (!(lifting_size == 32 || lifting_size == 16 || lifting_size == 8)) {
        return 0;
    }

    uint16_t iters = 0;
    uint16_t parity_check_errors = 0;

    reset_decode_state (ldpc->R_mj, NUM_EDGES);

    for (iters = 0; iters < ldpc->cfg.max_iters; iters++) {

        parity_check_errors = check_syndrome (Lq, len, ldpc->Hm, ldpc->Hp, decoded);
        if (parity_check_errors == 0) {
            break;
        }

        #if USE_SUM_PRODUCT
            layered_sum_product (ldpc, Lq, len, lifting_size);
        #else
            layered_normalized_minsum (ldpc, Lq, len, lifting_size);
        #endif
    }
    
    *num_iters = iters;
    return parity_check_errors;
}