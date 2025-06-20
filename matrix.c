#include "matrix.h"

uint8_t find_vector_weight (vector_t x, uint8_t len) {
    uint8_t weight = 0;
    for (int j = 0; j < len/8; j++) {
        for (int i = 0; i < 8; i++) {
            weight += 0x01 & (x.data8[j] >> i);
        }
    }
    return weight;
}

void vector_add (vector_t x, vector_t y, vector_t result, uint8_t lifting_size, uint16_t length) {
    for (int i = 0; i < length; i++) {
        if (lifting_size == 32) {
            result.data32 [i] = x.data32 [i] ^ y.data32 [i];
        }
        else if (lifting_size == 16) {
            result.data16 [i] = x.data16 [i] ^ y.data16 [i];
        }
        else {
            result.data8 [i] = x.data8 [i] ^ y.data8 [i];
        }
    }
}

void circular_matrix_multiply (quasi_cyclic_matrix_t * H, vector_t x, vector_t b) {

    for (int j = 0; j < H->rows; j++) {

        for (int i = 0; i < H->rowWeight[j]; i++) {
            
            if (H->columnIndexMap[H->rowOffset[j]+i] >= H->cols) {
                break;
            }

            if (H->lifting_size == 32) {
                b.data32 [j] ^= (x.data32 [H->columnIndexMap[H->rowOffset[j]+i]] >> H->base_graph[H->rowOffset[j]+i]%H->lifting_size) | 
                (x.data32 [H->columnIndexMap[H->rowOffset[j]+i]] << (H->lifting_size - H->base_graph[H->rowOffset[j]+i]%H->lifting_size));
            }
            else if (H->lifting_size == 16) {
                b.data16 [j] ^= (x.data16 [H->columnIndexMap[H->rowOffset[j]+i]] >> H->base_graph[H->rowOffset[j]+i]%H->lifting_size) | 
                (x.data16 [H->columnIndexMap[H->rowOffset[j]+i]] << (H->lifting_size - H->base_graph[H->rowOffset[j]+i]%H->lifting_size));
            }
            else {
                b.data8 [j] ^= (x.data8 [H->columnIndexMap[H->rowOffset[j]+i]] >> H->base_graph[H->rowOffset[j]+i]%H->lifting_size) | 
                (x.data8 [H->columnIndexMap[H->rowOffset[j]+i]] << (H->lifting_size - H->base_graph[H->rowOffset[j]+i]%H->lifting_size));
            }
        }
    }
}