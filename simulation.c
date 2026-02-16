#include "simulation.h"
#include <stdlib.h>

uint8_t rand_u8_range(uint8_t min, uint8_t max)
{
    if (min > max) {
        return min;
    }

    return (uint8_t)(min + (rand() % (max - min + 1)));
}

uint8_t generate_random_bytes(uint8_t *buffer, uint16_t size) {
    if (!buffer || size == 0) {
        return 0;
    }

    for (size_t i = 0; i < size; i++) {
        buffer[i] = rand_u8_range(0, 255);
    }

    return size;
}