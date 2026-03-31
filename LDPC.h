#ifndef LDPC_H_
#define LDPC_H_

#include <stdint.h>
#include "matrix.h"

#define ALPHA 0.75

#define MIN_LIFTING_SIZE 8
#define MAX_LIFTING_SIZE 32
#define CODEWORD_SIZE 32

#define NUM_EDGES_MSG 245
#define NUM_EDGES_PARITY 71
#define NUM_EDGES (NUM_EDGES_MSG + NUM_EDGES_PARITY)
#define BG1_ROWS 46
#define BG1_COLS 68
#define MAX_CODE_LEN BG1_COLS*MAX_LIFTING_SIZE/8 //In bytes

#define MAX_ROW_WEIGHT 19

typedef struct{
    float element [MAX_LIFTING_SIZE];
}float_struct;

quasi_cyclic_matrix_t* get_Hm ();
quasi_cyclic_matrix_t* get_Hp ();
quasi_cyclic_matrix_t* get_A ();
quasi_cyclic_matrix_t* get_C ();
quasi_cyclic_matrix_t* get_D ();

void LDPC_init (uint8_t msg_size, uint8_t parity_size);
void LDPC_encode_init (uint8_t msg_size, uint8_t parity_size);
void LDPC_decode_init (uint8_t msg_size, uint8_t parity_size);

uint8_t LDPC_encode (uint8_t * data, uint16_t len, uint8_t * parity);
uint8_t LDPC_decode (float * Lq, uint16_t len, uint8_t * decoded, uint16_t max_iters, uint16_t * num_iters);

#endif /* LDPC_H_ */