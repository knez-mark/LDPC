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
        symbols[i] += noise_stddev * randn();
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

        static uint8_t parity [48] = {0};

        static float bpsk_symbols [32+128*8] = {0};
        
        static float ldpc_llr [136*8] = {0};

        memset (&message, 0, sizeof (message_t));
        memset (&encoding, 0, sizeof (encoded_t));
        memset (&decoding, 0, sizeof (encoded_t));
        memset (parity, 0, sizeof (parity));
        memset (bpsk_symbols, 0, sizeof (bpsk_symbols));
        memset (ldpc_llr, 0, sizeof (ldpc_llr));

        //Generate message
        message.msg_size = msg_len;
        generate_random_bytes(message.msg, message.msg_size);
        encoding.length_field = generate_length_field(message.msg_size + 1);

        //Encode data
        LDPC_encode ((uint8_t *) &message, get_LDPC_len(encoding.length_field), parity);
        format_encoded_data((uint8_t *) &message, parity, encoding.ldpc, encoding.length_field);

        //Transmit through channel
        uint16_t num_bytes = get_frame_len (encoding.length_field);
        uint16_t num_bits = bytes_to_bpsk((uint8_t*) &encoding, num_bytes, bpsk_symbols);
        add_awgn(bpsk_symbols, num_bits, noise_stddev);

        //Decode message
        decoding.length_field = len_field_decode_soft (bpsk_symbols);
        format_decoded_data(bpsk_symbols + 32, ldpc_llr, decoding.length_field);
        LDPC_decode (ldpc_llr, get_LDPC_len(decoding.length_field), decoding.ldpc, LDPC_NUM_ITERS, &num_iters);

        //Calculate BER
        bit_errors += calculate_bit_errors ((uint8_t*) &message, decoding.ldpc, message.msg_size + 1);
        num_transmissions += 8*(message.msg_size + 1);
    }

    return (float)bit_errors/num_transmissions;
}