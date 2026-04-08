#include "simulation.h"
#include "frame.h"
#include "LDPC.h"
#include "length_field.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <correct.h>

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

uint16_t bytes_to_bpsk_reverse_endianness(const uint8_t *bytes, uint16_t num_bytes, float* symbols) {
    uint16_t total_bits = num_bytes * 8;

    uint16_t idx = 0;
    for (uint16_t i = 0; i < num_bytes; i++) {
        for (int8_t bit = 7; bit >= 0; bit--) {
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
        symbols[i] *= 2/(noise_stddev*noise_stddev);
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

//Use reduced PCM
float find_bit_error_rate2 (uint8_t msg_len, float noise_stddev) {

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
            num_bytes = 34;
            lifting_size = 8;
        }
        else if (message.msg_size < 22) {
            LDPC_msg_size = message.msg_size + 1;
            num_bytes = 34;
            lifting_size = 8;
        }
        else if (message.msg_size < 44) {
            LDPC_msg_size = message.msg_size + 1;
            num_bytes = 34*2;
            lifting_size = 16;
        }
        else {
            LDPC_msg_size = message.msg_size + 1;
            num_bytes = 34*4;
            lifting_size = 32;
        }

        uint8_t LDPC_msg_size_for_encoder;
        LDPC_msg_size_for_encoder = (LDPC_msg_size*MIN_LIFTING_SIZE + lifting_size - 1)/lifting_size; //Implements ceiling operation
        LDPC_parity_size = 34 - LDPC_msg_size_for_encoder;

        LDPC_init(LDPC_msg_size_for_encoder, LDPC_parity_size);

        //Encode data
        LDPC_encode ((uint8_t *) &message, num_bytes*8, parity);

        encoding.ldpc[0] = message.msg_size;
        for (int i = 0; i < LDPC_msg_size_for_encoder*lifting_size/MIN_LIFTING_SIZE - 1; i++) {
            encoding.ldpc[1+i] = message.msg[i];
        }
        for (int i = 0; i < LDPC_parity_size*lifting_size/MIN_LIFTING_SIZE; i++) {
            encoding.ldpc[LDPC_msg_size_for_encoder*lifting_size/MIN_LIFTING_SIZE+i] = parity [i];
        }

        //Transmit through channel
        uint16_t num_bits = bytes_to_bpsk(encoding.ldpc, num_bytes, bpsk_symbols);
        add_awgn(bpsk_symbols, num_bits, noise_stddev);

        for (int i = 0; i < 2*lifting_size; i++) {
            bpsk_symbols [i] = 0;
        }

        //for (int i = 8*LDPC_msg_size; i < LDPC_msg_size_for_encoder*lifting_size; i++) {
        //    bpsk_symbols [i] = 10.0f;
        //}

        //Decode message
        LDPC_decode (bpsk_symbols, num_bytes*8, decoding.ldpc, LDPC_NUM_ITERS, &num_iters);

        //Calculate BER
        bit_errors += calculate_bit_errors ((uint8_t*) &message, decoding.ldpc, LDPC_msg_size);
        num_transmissions += 8*(LDPC_msg_size);
    }

    return (float)bit_errors/num_transmissions;
}

static void cc_puncture (float* symbols, uint16_t len) {
    for (int i = 0; i < len; i++) {
        if ((i - 3) % 6 == 0 || (i - 4) % 6 == 0) {
            symbols [i] = 0.0f;
        }
    }
}

#define CC_SOFT_ONE 255
#define CC_SOFT_ZERO 0
#define CC_SOFT_ERASURE 128

static void cc_quantize_hard (float* symbols, uint16_t len, uint8_t* quantized) {
    for (uint16_t i = 0; i < len; i++) {
        if (symbols[i] > 0.0f) {
            quantized [i] = CC_SOFT_ZERO;
        }
        else if (symbols[i] < 0.0f) {
            quantized [i] = CC_SOFT_ONE;
        }
        else {
            quantized [i] = CC_SOFT_ERASURE;
        }
    }
}

#define L_MAX 8.0f

static void cc_quantize_soft (float* symbols, uint16_t len, uint8_t* quantized) {
    for (uint16_t i = 0; i < len; i++) {
        float L = symbols [i];
        if (L > L_MAX) L = L_MAX;
        if (L < -L_MAX) L = -L_MAX;

        float scale = -127.0f / L_MAX;
        int q = (int)roundf(scale*L);

        if (q > 127) q = 127;
        if (q < -127) q = -127;
        q = q + 128;

        quantized [i] = (uint8_t)q;
    }
}

//Use reduced PCM
float find_bit_error_rate_rs_cc (uint8_t msg_len, float noise_stddev, bool is_soft) {

    uint32_t num_transmissions = 0;
    uint32_t bit_errors = 0;
    uint16_t num_iters;

    while (num_transmissions < MAX_TRANSMISSIONS && bit_errors < TARGET_ERRORS) {
        static message_t message = {0};

        static encoded_t encoding_rs = {0};
        static encoded_t encoding_rs_cc = {0};
        static encoded_t decoding_rs_cc = {0};
        static encoded_t decoding_rs = {0};

        static float bpsk_symbols [8*(18*4*2+3)] = {0};
        static uint8_t quantized_bpsk_symbols [8*(18*4*2+3)] = {0};
        
        memset (&message, 0, sizeof (message_t));
        memset (&encoding_rs, 0, sizeof (encoded_t));
        memset (&encoding_rs_cc, 0, sizeof (encoded_t));
        memset (&decoding_rs_cc, 0, sizeof (encoded_t));
        memset (&decoding_rs, 0, sizeof (encoded_t));
        memset (bpsk_symbols, 0, sizeof (bpsk_symbols));
        memset (quantized_bpsk_symbols, 0, sizeof (quantized_bpsk_symbols));

        //Generate message
        message.msg_size = msg_len;
        generate_random_bytes(message.msg, message.msg_size);

        uint8_t rs_msg_size, cc_msg_size, rs_num_roots;

        if (message.msg_size < 16) {
            rs_msg_size = 16;
            cc_msg_size = 18;
            rs_num_roots = 2;
        }
        else if (message.msg_size < 16*2) {
            rs_msg_size = 16*2;
            cc_msg_size = 18*2;
            rs_num_roots = 2*2;
        }
        else if (message.msg_size < 16*3) {
            rs_msg_size = 16*3;
            cc_msg_size = 18*3;
            rs_num_roots = 2*3;
        }
        else {
            rs_msg_size = 16*4;
            cc_msg_size = 18*4;
            rs_num_roots = 2*4;
        }

        static correct_convolutional *conv;
        //static const correct_convolutional_polynomial_t correct_conv_r12_7_polynomial [] = {0133, 0171};
        conv = correct_convolutional_create(2, 7, correct_conv_r12_7_polynomial);
        //static const correct_convolutional_polynomial_t correct_conv_r12_4_polynomial [] = {05, 07};
        //conv = correct_convolutional_create(2, 4, correct_conv_r12_4_polynomial);

        static correct_reed_solomon *rs;
        rs = correct_reed_solomon_create(correct_rs_primitive_polynomial_8_4_3_2_0, 0, 1, rs_num_roots);

        //Encode data
        correct_reed_solomon_encode (rs, (uint8_t *)&message, rs_msg_size, (uint8_t *)&encoding_rs);
        correct_convolutional_encode (conv, (uint8_t *)&encoding_rs, cc_msg_size, (uint8_t *)&encoding_rs_cc);

        //Transmit through channel
        uint16_t num_bits = bytes_to_bpsk_reverse_endianness((uint8_t *)&encoding_rs_cc, (correct_convolutional_encode_len(conv, cc_msg_size) + 8)/8, bpsk_symbols);
        add_awgn(bpsk_symbols, num_bits, noise_stddev);

        //Puncture and quantize
        cc_puncture (bpsk_symbols, num_bits);
        if (is_soft) {
            cc_quantize_soft (bpsk_symbols, num_bits, quantized_bpsk_symbols);
        } else {
            cc_quantize_hard (bpsk_symbols, num_bits, quantized_bpsk_symbols);
        }

        //Decode message
        correct_convolutional_decode_soft (conv, quantized_bpsk_symbols, correct_convolutional_encode_len(conv, cc_msg_size), (uint8_t *)&decoding_rs_cc);
        correct_reed_solomon_decode (rs, (uint8_t *)&decoding_rs_cc, cc_msg_size, (uint8_t *)&decoding_rs);

        correct_convolutional_destroy (conv);
        correct_reed_solomon_destroy (rs);

        //Calculate BER
        bit_errors += calculate_bit_errors ((uint8_t*) &message, (uint8_t *)&decoding_rs, rs_msg_size);
        num_transmissions += 8*(rs_msg_size);
    }

    return (float)bit_errors/num_transmissions;
}