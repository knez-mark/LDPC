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

void circular_matrix_multiply (quasi_cyclic_matrix_t * H, vector_t x, vector_t b) {

    for (int j = 0; j < H->rows; j++) {

        for (int i = 0; i < H->rowWeight[j]; i++) {
            
            if (H->lifting_size == 32) {
                b.data32 [j] ^= (x.data32 [H->columnIndexMap[H->rowOffset[j]+i]] << H->base_graph[H->rowOffset[j]+i]%H->lifting_size) | 
                (x.data32 [H->columnIndexMap[H->rowOffset[j]+i]] >> (H->lifting_size - H->base_graph[H->rowOffset[j]+i]%H->lifting_size));
            }
            else if (H->lifting_size == 16) {
                b.data16 [j] ^= (x.data16 [H->columnIndexMap[H->rowOffset[j]+i]] << H->base_graph[H->rowOffset[j]+i]%H->lifting_size) | 
                (x.data16 [H->columnIndexMap[H->rowOffset[j]+i]] >> (H->lifting_size - H->base_graph[H->rowOffset[j]+i]%H->lifting_size));
            }
            else {
                b.data8 [j] ^= (x.data8 [H->columnIndexMap[H->rowOffset[j]+i]] >> H->base_graph[H->rowOffset[j]+i]%H->lifting_size) | 
                (x.data8 [H->columnIndexMap[H->rowOffset[j]+i]] << (H->lifting_size - H->base_graph[H->rowOffset[j]+i]%H->lifting_size));
            }
        }
    }
}