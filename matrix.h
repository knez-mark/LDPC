#ifndef MATRIX_H_
#define MATRIX_H_

#include <stdint.h>

typedef struct {
    uint8_t rows;
    uint8_t cols;
    uint8_t lifting_size;

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

uint8_t find_vector_weight (vector_t x, uint8_t len);
void circular_matrix_multiply (quasi_cyclic_matrix_t * H, vector_t x, vector_t b);

#endif /* MATRIX_H_ */