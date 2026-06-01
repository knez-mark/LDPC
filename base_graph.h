#ifndef BASE_GRAPH_H_
#define BASE_GRAPH_H_

#include "matrix.h"

const quasi_cyclic_matrix_t * get_BG1_msg (uint8_t idx);
const quasi_cyclic_matrix_t * get_BG1_parity (uint8_t idx);
const quasi_cyclic_matrix_t * get_BG2_msg (uint8_t idx);
const quasi_cyclic_matrix_t * get_BG2_parity (uint8_t idx);

uint8_t is_valid_lifting_size (uint16_t lifting_size);
uint16_t get_nearest_lifting_size (uint16_t msg_len, uint16_t code_len, uint8_t bgn);
int8_t get_set_index (uint16_t lifting_size);

#endif /* BASE_GRAPH_H_ */