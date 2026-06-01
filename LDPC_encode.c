#include "LDPC_impl.h"
#include "base_graph.h"
#include <string.h>

static void mult_by_B_inv (B_matrix_t B, uint8_t * input, uint8_t * result, uint16_t lifting_size) {
    uint16_t block_size = (lifting_size + (NUM_BITS_PER_BYTE - 1)) / NUM_BITS_PER_BYTE;

    if (B == BG1_B1) {
        vector_add4 (input, 
                    input + 1*block_size, 
                    input + 2*block_size, 
                    input + 3*block_size, 
                    result, 
                    lifting_size);
        
        circular_shift (result, 
                        result + 3*block_size,
                        lifting_size,
                        1);

        vector_add2 (input,
                    result + 3*block_size,
                    result + 1*block_size,
                    lifting_size);

        vector_add2 (input + 3*block_size,
                    result + 3*block_size,
                    result + 3*block_size,
                    lifting_size);

        vector_add2 (input + 2*block_size,
                    result + 3*block_size,
                    result + 2*block_size,
                    lifting_size);
    }
    else if (B == BG1_B2) {
        vector_add4 (input, 
                    input + 1*block_size, 
                    input + 2*block_size, 
                    input + 3*block_size, 
                    result + 1*block_size, 
                    lifting_size);
        
        circular_shift (result + 1*block_size,  
                        result,
                        lifting_size,
                        lifting_size - 105 % lifting_size);

        vector_add2 (input,
                    result,
                    result + 1*block_size,
                    lifting_size);

        vector_add2 (input + 3*block_size,
                    result,
                    result + 3*block_size,
                    lifting_size);   

        vector_add2 (input + 2*block_size,
                    result + 3*block_size,
                    result + 2*block_size,
                    lifting_size);   
    }
    else if (B == BG2_B1) {
        vector_add4 (input, 
                    input + 1*block_size, 
                    input + 2*block_size, 
                    input + 3*block_size, 
                    result + 1*block_size, 
                    lifting_size);
        
        circular_shift (result + 1*block_size,  
                        result,
                        lifting_size,
                        lifting_size - 1);
        
        vector_add2 (input,
                    result,
                    result + 1*block_size,
                    lifting_size);

        vector_add2 (input + 1*block_size,
                    result + 1*block_size,
                    result + 2*block_size,
                    lifting_size);

        vector_add2 (input + 3*block_size,
                    result,
                    result + 3*block_size,
                    lifting_size);
    }
    else {
        vector_add4 (input, 
                    input + 1*block_size, 
                    input + 2*block_size, 
                    input + 3*block_size, 
                    result, 
                    lifting_size);

        circular_shift (result,  
                        result + 1*block_size,
                        lifting_size,
                        1);

        vector_add2 (input + 3*block_size,
                    result + 1*block_size,
                    result + 3*block_size,
                    lifting_size);

        vector_add2 (input,
                    result + 1*block_size,
                    result + 1*block_size,
                    lifting_size);

        vector_add2 (input + 1*block_size,
                    result + 1*block_size,
                    result + 2*block_size,
                    lifting_size);
    }

}

uint8_t LDPC_encode (ldpc_encoder_t* ldpc, uint8_t * data, uint8_t * parity) {
    //Encoding algorithm taken from "Low-Latency QC-LDPC Encoder Design for 5G NR"
    //by Tian et al.

    uint16_t lifting_size = ldpc->cfg.lifting_size;
    uint16_t len = ldpc->cfg.msg_len;

    if (!is_valid_lifting_size(lifting_size)) {
        return 0;
    }

    uint8_t * A_mult_S = ldpc->A_mult_S;
    memset (A_mult_S, 0, sizeof (ldpc->A_mult_S));
    uint8_t * P1 = ldpc->parity;

    uint8_t * D_mult_P1 = ldpc->D_mult_P1;
    memset (D_mult_P1, 0, sizeof (ldpc->D_mult_P1));
    uint8_t * C_mult_S = ldpc->C_mult_S;
    memset (C_mult_S, 0, sizeof (ldpc->C_mult_S));
    
    uint16_t block_size = DIV_CEIL(lifting_size, NUM_BITS_PER_BYTE);
    uint8_t * P2 = P1 + 4*block_size;

    uint8_t *temp = ldpc->temp;

    // 1) Multiply A with S
    uint8_t * data_aligned = ldpc->data;
    get_byte_aligned(data, data_aligned, len, lifting_size);

    circular_matrix_multiply (&ldpc->A, data_aligned, A_mult_S, temp, lifting_size);
    // 2) Find B^-1 * (A * S) to obtain P1
    mult_by_B_inv (ldpc->B, A_mult_S, P1, lifting_size);
    // 3) Multiply D with P1
    circular_matrix_multiply (&ldpc->D, P1, D_mult_P1, temp, lifting_size);
    // 4) Multiply C with S
    circular_matrix_multiply (&ldpc->C, data_aligned, C_mult_S, temp, lifting_size);
    // 5) Add the results from 3) and 4)
    vector_add2 (D_mult_P1, C_mult_S, P2, NUM_BITS_PER_BYTE*ldpc->C.rows*block_size);

    get_packed (P1, parity, (ldpc->A.rows + ldpc->C.rows)*lifting_size, lifting_size);

    return 1;
}
