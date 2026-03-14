#include "simulation.h"
#include "frame.h"
#include "LDPC.h"
#include "length_field.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

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

uint16_t bytes_to_bpsk(const uint8_t *bytes, uint16_t num_bytes, float* symbols) {
    uint16_t total_bits = num_bytes * 8;

    uint16_t idx = 0;
    for (uint16_t i = 0; i < num_bytes; i++) {
        for (uint8_t bit = 0; bit < 8; bit++) {
            uint8_t b = (bytes[i] >> bit) & 1;
            symbols[idx++] = b ? -1.0f : 1.0f;
        }
    }

    return total_bits;
}

static float randn(void) {
    float u1 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);
    float u2 = ((float)rand() + 1.0f) / ((float)RAND_MAX + 2.0f);

    return sqrtf(-2.0f * logf(u1)) * cosf(2.0f * (float)M_PI * u2);
}

void add_awgn(float *symbols, uint16_t len, float noise_stddev) {
    for (uint16_t i = 0; i < len; i++) {
        symbols[i] += noise_stddev*randn();
        symbols[i] *= 2/(noise_stddev*noise_stddev)/500;
    }
}

static uint16_t calculate_bit_errors (uint8_t *msg, uint8_t *decoded, uint8_t msg_len) {
    uint8_t diff_vec;
    uint16_t bit_errors = 0;

    for (uint8_t i = 0; i < msg_len; i++) {
        diff_vec = msg[i] ^ decoded[i];
        for (uint8_t j = 0; j < 8; j++) {
            bit_errors += 0x01 & (diff_vec >> j);
        }
    }   
    return bit_errors;
}

float find_bit_error_rate (uint8_t msg_len, float noise_stddev) {

    uint32_t num_transmissions = 0;
    uint32_t bit_errors = 0;
    uint16_t num_iters;

    while (num_transmissions < MAX_TRANSMISSIONS && bit_errors < TARGET_ERRORS) {
        static message_t message = {0};

        static encoded_t encoding = {0};
        static encoded_t decoding = {0};

        static uint8_t parity [272] = {0};

        static float bpsk_symbols [272*8] = {0};
        
        memset (&message, 0, sizeof (message_t));
        memset (&encoding, 0, sizeof (encoded_t));
        memset (&decoding, 0, sizeof (encoded_t));
        memset (parity, 0, sizeof (parity));
        memset (bpsk_symbols, 0, sizeof (bpsk_symbols));

        //Generate message
        message.msg_size = msg_len;
        generate_random_bytes(message.msg, message.msg_size);

        uint8_t LDPC_msg_size;
        uint8_t LDPC_parity_size;
        uint16_t num_bytes;
        uint8_t lifting_size;

        if (message.msg_size == 0) {
            LDPC_msg_size = 2;
            num_bytes = 68;
            lifting_size = 8;
        }
        else if (message.msg_size < 22) {
            LDPC_msg_size = message.msg_size + 1;
            num_bytes = 68;
            lifting_size = 8;
        }
        else if (message.msg_size < 44) {
            LDPC_msg_size = message.msg_size + 1;
            num_bytes = 68*2;
            lifting_size = 16;
        }
        else {
            LDPC_msg_size = message.msg_size + 1;
            num_bytes = 68*4;
            lifting_size = 32;
        }
        LDPC_parity_size = 34*lifting_size/MIN_LIFTING_SIZE - LDPC_msg_size;

        //Encode data
        LDPC_encode ((uint8_t *) &message, num_bytes*8, parity);

        encoding.ldpc[0] = message.msg_size;
        for (int i = 0; i < 22*lifting_size/MIN_LIFTING_SIZE - 1; i++) {
            encoding.ldpc[1+i] = message.msg[i];
        }
        for (int i = 0; i < LDPC_parity_size; i++) {
            encoding.ldpc[22*lifting_size/MIN_LIFTING_SIZE+i] = parity [i];
        }

        //Transmit through channel
        uint16_t num_bits = bytes_to_bpsk(encoding.ldpc, num_bytes, bpsk_symbols);
        add_awgn(bpsk_symbols, num_bits, noise_stddev);

        for (int i = 0; i < 2*lifting_size; i++) {
            bpsk_symbols [i] = 0;
        }

        for (int i = 8*(message.msg_size + 1); i < lifting_size*(22); i++) {
            bpsk_symbols [i] = 10.0f;
        }

        for (int i = 8*num_bytes - lifting_size*(BG1_COLS - 22) + 8*LDPC_parity_size; i < 8*num_bytes; i++) {
            bpsk_symbols [i] = 0;
        }

        //Decode message
        LDPC_decode (bpsk_symbols, num_bytes*8, decoding.ldpc, LDPC_NUM_ITERS, &num_iters);

        //Calculate BER
        bit_errors += calculate_bit_errors ((uint8_t*) &message, decoding.ldpc, LDPC_msg_size);
        num_transmissions += 8*(LDPC_msg_size);
    }

    return (float)bit_errors/num_transmissions;
}