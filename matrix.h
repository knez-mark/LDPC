#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdint.h>

typedef struct {
    uint8_t rows;
    uint8_t cols;

    const uint8_t *base_graph;
    struct {
        const uint16_t *rowOffset;
        const uint8_t *rowWeight;
        const uint8_t *columnIndexMap;
    };
} quasi_cyclic_matrix_t;

typedef union{
    uint8_t *data8;
    uint16_t *data16;
    uint32_t *data32;
} vector_t;

uint16_t find_vector_weight (vector_t x, uint16_t len);
void vector_add2 (uint8_t *x, uint8_t *y, uint8_t *result, uint16_t length);
void vector_add4 (uint8_t *w, uint8_t *x, uint8_t *y, uint8_t *z, uint8_t *result, uint16_t length);
void get_byte_aligned (uint8_t *input, uint8_t *output, uint16_t len_bits, uint16_t lifting_size);
void get_packed (uint8_t *input, uint8_t *output, uint16_t len_bits, uint16_t lifting_size);
void circular_shift (uint8_t *input, uint8_t *output, uint16_t lifting_size, uint16_t shift);
void circular_matrix_multiply (quasi_cyclic_matrix_t * H, vector_t x, vector_t b, uint8_t lifting_size);
void circular_matrix_multiply_general (quasi_cyclic_matrix_t * H, uint8_t *x, uint8_t *b, uint8_t *temp, uint8_t lifting_size);

#endif /* MATRIX_H_ */