#ifndef LDPC_IMPL_H_
#define LDPC_IMPL_H_

#include "LDPC.h"
#include "matrix.h"

#define EIGHT_BITS_PER_BYTE 8

typedef struct{
    #if USE_SUM_PRODUCT
    float element [MAX_LIFTING_SIZE];
    #else
    int8_t element [MAX_LIFTING_SIZE];
    #endif
}lifting_struct;

struct ldpc_encoder {
    quasi_cyclic_matrix_t* A;
    quasi_cyclic_matrix_t* C;
    quasi_cyclic_matrix_t* D;

    ldpc_encoder_cfg_t cfg;
};

struct ldpc_decoder {
    quasi_cyclic_matrix_t* Hm;
    quasi_cyclic_matrix_t* Hp;

    lifting_struct* R_mj;

    ldpc_decoder_cfg_t cfg;
};

#define NUM_EDGES_MSG 245
#define NUM_EDGES_PARITY 71
#define NUM_EDGES (NUM_EDGES_MSG + NUM_EDGES_PARITY)
#define BG1_ROWS 46
#define BG1_COLS 68
#define MAX_CODE_LEN BG1_COLS*MAX_LIFTING_SIZE/EIGHT_BITS_PER_BYTE //In bytes

#define MAX_ROW_WEIGHT 19

#endif /* LDPC_IMPL_H_ */