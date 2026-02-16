#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <stdint.h>

uint8_t rand_u8_range(uint8_t min, uint8_t max);
uint8_t generate_random_bytes(uint8_t *buffer, uint16_t size);

uint16_t bytes_to_bpsk(const uint8_t *bytes, uint16_t num_bytes, float* symbols);
void add_awgn(float *symbols, uint16_t len, float noise_stddev);

#endif /* SIMULATION_H_ */