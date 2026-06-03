#ifndef SIMULATION_H_
#define SIMULATION_H_

#include <stdint.h>
#include "LDPC_impl.h"

typedef struct {
    uint8_t  msg[MAX(BG1_COLS-BG1_ROWS, BG2_COLS-BG2_ROWS) *
                 DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)];
} message_t;

typedef struct {
    uint32_t length_field;
    uint8_t  ldpc[MAX(BG1_COLS, BG2_COLS) *
                  DIV_CEIL(MAX_LIFTING_SIZE_DECODE, NUM_BITS_PER_BYTE)];
} encoded_t;

#define BGN 1
#define USE_CONST_RATE false
#define USE_EXPLICIT_LIFTING_SIZE true
#define USE_FULL_PCM false

#define MAX_TRANSMISSIONS 100000000
#define TARGET_ERRORS 1000
#define LDPC_NUM_ITERS 5

uint8_t rand_u8_range(uint8_t min, uint8_t max);
uint8_t generate_random_bytes(uint8_t *buffer, uint16_t size);

uint16_t bytes_to_bpsk(const uint8_t *bytes, uint16_t num_bytes, float* symbols);
void add_awgn(float *symbols, uint16_t len, float noise_stddev);

float find_bit_error_rate (uint8_t msg_len, float noise_stddev);

#endif /* SIMULATION_H_ */