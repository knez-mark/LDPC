#ifndef LDPC_IMPL_H_
#define LDPC_IMPL_H_

#include "LDPC.h"
#include "matrix.h"

#define DIV_CEIL(a, b) (((a) + (b) - 1) / (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

#define NUM_BITS_PER_BYTE 8

#define BG1_NUM_EDGES_MSG 245
#define BG1_NUM_EDGES_PARITY 71
#define BG1_NUM_EDGES (BG1_NUM_EDGES_MSG + BG1_NUM_EDGES_PARITY)
#define BG1_ROWS 46
#define BG1_COLS 68

#define BG2_NUM_EDGES_MSG 113
#define BG2_NUM_EDGES_PARITY 84
#define BG2_NUM_EDGES (BG2_NUM_EDGES_MSG + BG2_NUM_EDGES_PARITY)
#define BG2_ROWS 42
#define BG2_COLS 52

#define BG1_MAX_ROW_WEIGHT 19
#define BG2_MAX_ROW_WEIGHT 10
#define MAX_ROW_WEIGHT MAX(BG1_MAX_ROW_WEIGHT, BG2_MAX_ROW_WEIGHT)

typedef struct{
    #if USE_SUM_PRODUCT
    float element [MAX_LIFTING_SIZE_DECODE];
    #else
    ldpc_quantized_t element [MAX_LIFTING_SIZE_DECODE];
    #endif
}lifting_struct;

typedef enum {
    BG1_B1,
    BG1_B2,
    BG2_B1,
    BG2_B2
} B_matrix_t;

struct ldpc_encoder {
    quasi_cyclic_matrix_t A;
    B_matrix_t B;
    quasi_cyclic_matrix_t C;
    quasi_cyclic_matrix_t D;

    uint8_t A_mult_S [4*DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];
    uint8_t D_mult_P1 [MAX(BG1_COLS-4, BG2_COLS-4)*DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];
    uint8_t C_mult_S [MAX(BG1_COLS-4, BG2_COLS-4)*DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];

    uint8_t data [MAX(BG1_COLS-BG1_ROWS, BG2_COLS-BG2_ROWS)*DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];
    uint8_t parity [MAX(BG1_ROWS, BG2_ROWS)*DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];

    uint8_t temp [DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];

    ldpc_encoder_cfg_t cfg;

    uint8_t msg_size;
    uint8_t parity_size;
};

struct ldpc_decoder {
    quasi_cyclic_matrix_t Hm;
    quasi_cyclic_matrix_t Hp;

    uint8_t syndrome [MAX(BG1_ROWS, BG2_ROWS)*DIV_CEIL(MAX_LIFTING_SIZE_DECODE, NUM_BITS_PER_BYTE)];
    uint8_t codeword [MAX(BG1_COLS, BG2_COLS)*DIV_CEIL(MAX_LIFTING_SIZE_DECODE, NUM_BITS_PER_BYTE)];

    uint8_t temp [DIV_CEIL(MAX_LIFTING_SIZE_DECODE, NUM_BITS_PER_BYTE)];

    lifting_struct R_mj [MAX(BG1_NUM_EDGES, BG2_NUM_EDGES)];

    uint16_t columnOffset [MAX_ROW_WEIGHT];
	uint16_t circular_shift [MAX_ROW_WEIGHT];

    #if USE_SUM_PRODUCT
    float Lq_mj [MAX_ROW_WEIGHT];
    float prefix[MAX_ROW_WEIGHT];
    float suffix[MAX_ROW_WEIGHT];

    float Lq [BG1_COLS*MAX_LIFTING_SIZE];
    #else
    ldpc_quantized_t Lq_mj [MAX_ROW_WEIGHT];

    ldpc_quantized_t Lq [MAX(BG1_COLS, BG2_COLS)*MAX_LIFTING_SIZE_DECODE];
    #endif

    ldpc_decoder_cfg_t cfg;

    uint8_t msg_size;
    uint8_t parity_size;
};

#endif /* LDPC_IMPL_H_ */