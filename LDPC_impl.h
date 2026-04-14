#ifndef LDPC_IMPL_H_
#define LDPC_IMPL_H_

#include "LDPC.h"
#include "matrix.h"

#define EIGHT_BITS_PER_BYTE 8

#define NUM_EDGES_MSG 245
#define NUM_EDGES_PARITY 71
#define NUM_EDGES (NUM_EDGES_MSG + NUM_EDGES_PARITY)
#define BG1_ROWS 46
#define BG1_COLS 68
#define MAX_CODE_LEN BG1_COLS*MAX_LIFTING_SIZE/EIGHT_BITS_PER_BYTE //In bytes

#define MAX_ROW_WEIGHT 19

typedef struct{
    #if USE_SUM_PRODUCT
    float element [MAX_LIFTING_SIZE];
    #else
    ldpc_quantized_t element [MAX_LIFTING_SIZE];
    #endif
}lifting_struct;

struct ldpc_encoder {
    quasi_cyclic_matrix_t A;
    quasi_cyclic_matrix_t C;
    quasi_cyclic_matrix_t D;

    uint8_t A_mult_S [(4*MAX_LIFTING_SIZE+(EIGHT_BITS_PER_BYTE-1))/EIGHT_BITS_PER_BYTE];
    uint8_t D_mult_P1 [((BG1_COLS-4)*MAX_LIFTING_SIZE+(EIGHT_BITS_PER_BYTE-1))/EIGHT_BITS_PER_BYTE];
    uint8_t C_mult_S [((BG1_COLS-4)*MAX_LIFTING_SIZE+(EIGHT_BITS_PER_BYTE-1))/EIGHT_BITS_PER_BYTE];

    ldpc_encoder_cfg_t cfg;
};

struct ldpc_decoder {
    quasi_cyclic_matrix_t Hm;
    quasi_cyclic_matrix_t Hp;

    uint8_t syndrome [(BG1_ROWS*MAX_LIFTING_SIZE+(EIGHT_BITS_PER_BYTE-1))/EIGHT_BITS_PER_BYTE];

    lifting_struct R_mj [NUM_EDGES];

    uint16_t columnOffset [MAX_ROW_WEIGHT];
	uint16_t circular_shift [MAX_ROW_WEIGHT];

    #if USE_SUM_PRODUCT
    float Lq_mj [MAX_ROW_WEIGHT];
    float prefix[MAX_ROW_WEIGHT];
    float suffix[MAX_ROW_WEIGHT];
    #else
    ldpc_quantized_t Lq_mj [MAX_ROW_WEIGHT];
    #endif

    ldpc_decoder_cfg_t cfg;
};

uint8_t is_valid_lifting_size (uint16_t lifting_size);

#endif /* LDPC_IMPL_H_ */