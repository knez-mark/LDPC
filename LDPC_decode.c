#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "LDPC.h"
#include "matrix.h"
#include "minsum.h"

static void get_hard_decision_codeword (float * Lq, uint16_t len, uint8_t lifting_size, vector_t codeword) {
    int8_t bit_ptr = 0;
    uint8_t byte_ptr = 0;

    for (int i = 0; i < len; i++) {
        if (lifting_size == 32) {
            if (Lq [i] < 0) { 
                codeword.data32 [byte_ptr] |= 1 << bit_ptr;
            }
        }
        else if (lifting_size == 16) {
            if (Lq [i] < 0) { 
                codeword.data16 [byte_ptr] |= 1 << bit_ptr;
            }
        }
        else { //lifting size of 8
            if (Lq [i] < 0) { 
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

static uint16_t check_syndrome (float * Lq, uint16_t len, quasi_cyclic_matrix_t * Hm, quasi_cyclic_matrix_t * Hp, uint8_t * codeword) {

    uint8_t lifting_size = len/(Hm->cols + Hp->cols);

    uint8_t syndrome [MAX_CODE_LEN] = {0};
    memset (syndrome, 0, MAX_CODE_LEN);
    memset (codeword, 0, MAX_CODE_LEN);

    get_hard_decision_codeword (Lq, len, lifting_size, (vector_t) codeword);

    circular_matrix_multiply (Hm, (vector_t) codeword, (vector_t) syndrome, lifting_size);
    circular_matrix_multiply (Hp, (vector_t) (codeword + Hm->cols*lifting_size/MIN_LIFTING_SIZE), (vector_t) syndrome, lifting_size);

    return find_vector_weight ((vector_t) syndrome, len); //Returns number of parity check equation failures
}

uint8_t LDPC_decode (float * Lq, uint16_t len, uint8_t * decoded, uint16_t max_iters, uint16_t * num_iters) {

    quasi_cyclic_matrix_t * Hm = get_Hm ();
    quasi_cyclic_matrix_t * Hp = get_Hp ();

    uint8_t lifting_size = len/(Hm->cols + Hp->cols);

    if (!(lifting_size == 32 || lifting_size == 16 || lifting_size == 8)) {
        return 0;
    }

    uint16_t iters = 0;
    uint16_t parity_check_errors = 0;

    for (iters = 0; iters < max_iters; iters++) {

        parity_check_errors = check_syndrome (Lq, len, Hm, Hp, decoded);
        if (parity_check_errors == 0) {
            break;
        }

        #if USE_SUM_PRODUCT
            layered_sum_product (Lq, len, Hm, Hp);
        #else
            layered_normalized_minsum (Lq, len, Hm, Hp);
        #endif
    }

    reset_minsum ();
    
    *num_iters = iters;
    return parity_check_errors;
}