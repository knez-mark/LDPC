#ifndef LDPC_H_
#define LDPC_H_

#include <stdint.h>
#include "matrix.h"

#define ALPHA 0.75

#define MAX_LIFTING_SIZE 32
#define CODEWORD_SIZE 32

//"Condensed" base graph
#define NUM_EDGES 137
#define BG1_ROWS 12
#define BG1_COLS 34
#define MAX_CODE_LEN BG1_COLS*MAX_LIFTING_SIZE/8 //In bytes

//Original base graph
//#define NUM_EDGES 316
//#define BG1_ROWS 46
//#define BG1_COLS 68
//#define MAX_CODE_LEN BG1_COLS*MAX_LIFTING_SIZE/8 //In bytes

#define MAX_ROW_WEIGHT 19

typedef struct{
    float element [MAX_LIFTING_SIZE];
}float_struct;

quasi_cyclic_matrix_t* get_H ();
quasi_cyclic_matrix_t* get_A ();
quasi_cyclic_matrix_t* get_C ();
quasi_cyclic_matrix_t* get_D ();

void LDPC_init ();

void LDPC_encode (uint8_t * data, uint16_t len, uint8_t * parity);
uint8_t LDPC_decode (float * Lq, uint16_t len, uint8_t * decoded, uint16_t max_iters, uint16_t * num_iters);

#endif /* LDPC_H_ */