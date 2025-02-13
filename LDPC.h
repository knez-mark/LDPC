#ifndef LDPC_H_
#define LDPC_H_

#include <stdint.h>
#include "matrix.h"

#define ALPHA 0.75

#define MAX_LIFTING_SIZE 32

#define NUM_EDGES 316
#define BG1_ROWS 46
#define BG1_COLS 68
#define MAX_CODE_LEN BG1_COLS*MAX_LIFTING_SIZE/8 //In bytes

#define MAX_ROW_WEIGHT 19

typedef struct{
    float element [MAX_LIFTING_SIZE];
}float_struct;

quasi_cyclic_matrix_t* get_H ();

void LDPC_init ();

uint8_t LDPC_decode (float * Lq, uint16_t len, uint8_t * decoded, uint16_t max_iters, uint16_t * num_iters);

#endif /* LDPC_H_ */