#include "LDPC.h"
#include <stdio.h>

void mult_by_B_inv (vector_t input, vector_t result) {
    if (get_A()->lifting_size == 32) {
        result.data32 [0] = input.data32 [0] ^ input.data32 [1] ^ input.data32 [2] ^ input.data32 [3];
        result.data32 [1] = input.data32 [0] ^ ((result.data32 [0] >> 1) | (result.data32 [0] << 31));
        result.data32 [3] = input.data32 [3] ^ ((result.data32 [0] >> 1) | (result.data32 [0] << 31));
        result.data32 [2] = input.data32 [2] ^ result.data32 [3];
    }
    else if (get_A()->lifting_size == 16) {
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

void LDPC_encode (uint8_t * data, uint16_t len, uint8_t * parity) {
    //Encoding algorithm taken from "Low-Latency QC-LDPC Encoder Design for 5G NR"
    //by Tian et al.

    uint8_t lifting_size = len/get_H()->cols;

    get_H()->lifting_size = lifting_size;
    get_A()->lifting_size = lifting_size;
    get_C()->lifting_size = lifting_size;
    get_D()->lifting_size = lifting_size;

    //Maximum size is 4 * 32 bits
    uint8_t A_mult_S [4*4] = {0};
    uint8_t * P1 = parity;

    uint8_t D_mult_P1 [8*4] = {0};
    uint8_t C_mult_S [8*4] = {0};
    uint8_t * P2 = parity + 4*4*lifting_size/MAX_LIFTING_SIZE;

    // 1) Multiply A with S
    circular_matrix_multiply (get_A(), (vector_t) data, (vector_t) A_mult_S);
    // 2) Find B^-1 * (A * S) to obtain P1
    mult_by_B_inv ((vector_t) A_mult_S, (vector_t) P1);
    // 3) Multiply D with P1
    circular_matrix_multiply (get_D(), (vector_t) P1, (vector_t) D_mult_P1);
    // 4) Multiply C with S
    circular_matrix_multiply (get_C(), (vector_t) data, (vector_t) C_mult_S);
    // 5) Add the results from 3) and 4)
    vector_add ((vector_t) D_mult_P1, (vector_t) C_mult_S, (vector_t) P2, lifting_size, 8);
}
