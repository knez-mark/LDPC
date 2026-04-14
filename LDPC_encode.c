#include "LDPC_impl.h"
#include <string.h>

void mult_by_B_inv (vector_t input, vector_t result, uint8_t lifting_size) {
    if (lifting_size == 32) {
        result.data32 [0] = input.data32 [0] ^ input.data32 [1] ^ input.data32 [2] ^ input.data32 [3];
        result.data32 [1] = input.data32 [0] ^ ((result.data32 [0] >> 1) | (result.data32 [0] << 31));
        result.data32 [3] = input.data32 [3] ^ ((result.data32 [0] >> 1) | (result.data32 [0] << 31));
        result.data32 [2] = input.data32 [2] ^ result.data32 [3];
    }
    else if (lifting_size == 16) {
        result.data16 [0] = input.data16 [0] ^ input.data16 [1] ^ input.data16 [2] ^ input.data16 [3];
        result.data16 [1] = input.data16 [0] ^ ((result.data16 [0] >> 1) | (result.data16 [0] << 15));
        result.data16 [3] = input.data16 [3] ^ ((result.data16 [0] >> 1) | (result.data16 [0] << 15));
        result.data16 [2] = input.data16 [2] ^ result.data16 [3];
    }
    else { //lifting size of 8
        result.data8 [0] = input.data8 [0] ^ input.data8 [1] ^ input.data8 [2] ^ input.data8 [3];
        result.data8 [1] = input.data8 [0] ^ ((result.data8 [0] >> 1) | (result.data8 [0] << 7));
        result.data8 [3] = input.data8 [3] ^ ((result.data8 [0] >> 1) | (result.data8 [0] << 7));
        result.data8 [2] = input.data8 [2] ^ result.data8 [3];
    }
}

uint8_t LDPC_encode (ldpc_encoder_t* ldpc, uint8_t * data, uint16_t len, uint8_t * parity) {
    //Encoding algorithm taken from "Low-Latency QC-LDPC Encoder Design for 5G NR"
    //by Tian et al.

    uint8_t lifting_size = len/(ldpc->cfg.msg_size + ldpc->cfg.parity_size);

    if (!is_valid_lifting_size(lifting_size)) {
        return 0;
    }

    //Maximum size is 4 * 32 bits
    uint8_t * A_mult_S = ldpc->A_mult_S;
    memset (A_mult_S, 0, sizeof (ldpc->A_mult_S));
    uint8_t * P1 = parity;

    uint8_t * D_mult_P1 = ldpc->D_mult_P1;
    memset (D_mult_P1, 0, sizeof (ldpc->D_mult_P1));
    uint8_t * C_mult_S = ldpc->C_mult_S;
    memset (C_mult_S, 0, sizeof (ldpc->C_mult_S));
    uint8_t * P2 = parity + 4*lifting_size/EIGHT_BITS_PER_BYTE;

    // 1) Multiply A with S
    circular_matrix_multiply (&ldpc->A, (vector_t) data, (vector_t) A_mult_S, lifting_size);
    // 2) Find B^-1 * (A * S) to obtain P1
    mult_by_B_inv ((vector_t) A_mult_S, (vector_t) P1, lifting_size);
    // 3) Multiply D with P1
    circular_matrix_multiply (&ldpc->D, (vector_t) P1, (vector_t) D_mult_P1, lifting_size);
    // 4) Multiply C with S
    circular_matrix_multiply (&ldpc->C, (vector_t) data, (vector_t) C_mult_S, lifting_size);
    // 5) Add the results from 3) and 4)
    vector_add ((vector_t) D_mult_P1, (vector_t) C_mult_S, (vector_t) P2, ldpc->C.rows*lifting_size/EIGHT_BITS_PER_BYTE);

    return 1;
}
