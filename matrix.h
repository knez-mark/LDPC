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
void vector_add (vector_t x, vector_t y, vector_t result, uint16_t length);
void circular_matrix_multiply (quasi_cyclic_matrix_t * H, vector_t x, vector_t b, uint8_t lifting_size);

#endif /* MATRIX_H_ */