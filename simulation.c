#include "simulation.h"
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

#define TYPE_MAX_T(T) _Generic((T)0, \
    int8_t: INT8_MAX, \
    int16_t: INT16_MAX, \
    int32_t: INT32_MAX, \
    int64_t: INT64_MAX, \
    uint8_t: UINT8_MAX, \
    uint16_t: UINT16_MAX, \
    uint32_t: UINT32_MAX, \
    uint64_t: UINT64_MAX \
)

#define L_MAX 8.0f

static void quantize (float* symbols, uint16_t len, ldpc_quantized_t* quantized) {
    for (uint16_t i = 0; i < len; i++) {
        float L = symbols [i];
        if (L > L_MAX) L = L_MAX;
        if (L < -L_MAX) L = -L_MAX;

        float scale = TYPE_MAX_T(ldpc_quantized_t) / L_MAX /4; //Using full range causes issues (divide by 4)
        int q = (int)roundf(scale*L);

        if (q > TYPE_MAX_T(ldpc_quantized_t)) q = TYPE_MAX_T(ldpc_quantized_t);
        if (q < -TYPE_MAX_T(ldpc_quantized_t)) q = -TYPE_MAX_T(ldpc_quantized_t);

        quantized [i] = (ldpc_quantized_t)q;
    }
}

typedef struct {
    uint16_t num_bytes;
    uint16_t lifting_size;
} ldpc_settings_t;

static ldpc_settings_t get_codeword_settings (uint16_t msg_len) {
    uint16_t num_bytes;
    uint16_t lifting_size;

    #if (BGN == 1)
        if (msg_len < 22) {
            num_bytes = 34;
            lifting_size = 8;
        }
        else if (msg_len < 44) {
            num_bytes = 34*2;
            lifting_size = 16;
        }
        else if (msg_len < 66) {
            num_bytes = 34*3;
            lifting_size = 24;
        }
        else {
            num_bytes = 34*4;
            lifting_size = 32;
        }
    #else 
        if (msg_len < 20) {
            num_bytes = 34;
            lifting_size = 16;
        }
        else if (msg_len < 40) {
            num_bytes = 34*2;
            lifting_size = 32;
        }
        else if (msg_len < 60) {
            num_bytes = 34*3;
            lifting_size = 48;
        }
        else if (msg_len < 80) {
            num_bytes = 34*4;
            lifting_size = 64;
        }
        else {
            num_bytes = 34*5;
            lifting_size = 80;
        }
    #endif

    #if USE_CONST_RATE
        num_bytes = (8*(msg_len + 1))/((float)22/34);
        num_bytes = (num_bytes + 7)/8;
    #endif

    ldpc_settings_t settings = {
        .num_bytes = num_bytes,
        .lifting_size = lifting_size
    };

    return settings;
}

float find_bit_error_rate (uint8_t msg_len, float noise_stddev) {

    uint32_t num_transmissions = 0;
    uint32_t bit_errors = 0;
    uint16_t num_iters;

    while (num_transmissions < MAX_TRANSMISSIONS && bit_errors < TARGET_ERRORS) {
        static message_t message = {0};

        static encoded_t encoding = {0};
        static encoded_t decoding = {0};

        static uint8_t parity [MAX(BG1_ROWS, BG2_ROWS) *
                                DIV_CEIL(MAX_LIFTING_SIZE_ENCODE, NUM_BITS_PER_BYTE)] = {0};

        static float bpsk_symbols [MAX(BG1_COLS, BG2_COLS)*MAX_LIFTING_SIZE_DECODE] = {0};
        static ldpc_quantized_t quantized_symbols 
                                [MAX(BG1_COLS, BG2_COLS)*MAX_LIFTING_SIZE_DECODE] = {0};
        
        memset (&message, 0, sizeof (message_t));
        memset (&encoding, 0, sizeof (encoded_t));
        memset (&decoding, 0, sizeof (encoded_t));
        memset (parity, 0, sizeof (parity));
        memset (bpsk_symbols, 0, sizeof (bpsk_symbols));
        memset (quantized_symbols, 0, sizeof (quantized_symbols));

        //Generate message
        generate_random_bytes(message.msg, msg_len + 1);

        ldpc_settings_t settings = get_codeword_settings (msg_len);

        uint16_t num_bytes = settings.num_bytes;
        uint16_t lifting_size = settings.lifting_size;

        ldpc_encoder_cfg_t ldpc_encoder_cfg = {0};
        ldpc_decoder_cfg_t ldpc_decoder_cfg = {0};

        ldpc_encoder_cfg.msg_len = 8*(msg_len + 1);
        ldpc_encoder_cfg.target_code_len = 8*num_bytes;
        ldpc_encoder_cfg.bgn = BGN;
        #if USE_EXPLICIT_LIFTING_SIZE
            ldpc_encoder_cfg.lifting_mode = LDPC_LIFTING_EXPLICIT;
            ldpc_encoder_cfg.lifting_size = lifting_size; 
        #else
            ldpc_encoder_cfg.lifting_mode = LDPC_LIFTING_AUTO;
        #endif

        ldpc_decoder_cfg.msg_len = 8*(msg_len + 1);
        ldpc_decoder_cfg.target_code_len = 8*num_bytes;
        ldpc_decoder_cfg.bgn = BGN;
        ldpc_decoder_cfg.use_full_PCM = USE_FULL_PCM;
        #if USE_EXPLICIT_LIFTING_SIZE
            ldpc_decoder_cfg.lifting_mode = LDPC_LIFTING_EXPLICIT;
            ldpc_decoder_cfg.lifting_size = lifting_size; 
        #else
            ldpc_decoder_cfg.lifting_mode = LDPC_LIFTING_AUTO;
        #endif

        #if !USE_SUM_PRODUCT
        ldpc_decoder_cfg.alpha = 0.75f;
        #endif
        ldpc_decoder_cfg.max_iters = LDPC_NUM_ITERS;

        ldpc_encoder_t* ldpc_e = LDPC_encoder_create (ldpc_encoder_cfg);
        ldpc_decoder_t* ldpc_d = LDPC_decoder_create (ldpc_decoder_cfg);

        //Encode data
        LDPC_encode (ldpc_e, (uint8_t *) &message, parity);

        for (int i = 0; i < msg_len + 1; i++) {
            encoding.ldpc[i] = message.msg[i];
        }
        for (int i = msg_len + 1; i < (ldpc_decoder_cfg.target_code_len + 7)/8; i++) {
            encoding.ldpc[i] = parity [i - (msg_len + 1)];
        }

        //Transmit through channel
        uint16_t num_bits = bytes_to_bpsk(encoding.ldpc, (ldpc_decoder_cfg.target_code_len + 7)/8, bpsk_symbols);
        add_awgn(bpsk_symbols, ldpc_decoder_cfg.target_code_len, noise_stddev);

        ldpc_encoder_cfg = LDPC_get_encoder_config (ldpc_e);
        ldpc_decoder_cfg = LDPC_get_decoder_config (ldpc_d);

        if (num_transmissions == 0) {
            printf ("lifting size: %d\n", ldpc_encoder_cfg.lifting_size);
            printf ("codeword len: %d\n", ldpc_encoder_cfg.target_code_len);
        }

        for (int i = 0; i < 2*ldpc_decoder_cfg.lifting_size; i++) {
            bpsk_symbols [i] = 0;
        }

        quantize (bpsk_symbols, ldpc_decoder_cfg.target_code_len, quantized_symbols);

        //Decode message
        #if USE_SUM_PRODUCT
        LDPC_decode (ldpc_d, bpsk_symbols, decoding.ldpc, &num_iters);
        #else
        LDPC_decode (ldpc_d, quantized_symbols, decoding.ldpc, &num_iters);
        #endif
        //Calculate BER
        bit_errors += calculate_bit_errors ((uint8_t*) &message, decoding.ldpc, msg_len + 1);
        num_transmissions += 8*(msg_len + 1);

        LDPC_encoder_destroy (ldpc_e);
        LDPC_decoder_destroy (ldpc_d);
    }

    return (float)bit_errors/num_transmissions;
}

float find_packet_error_rate_lf (float noise_stddev) {
    uint32_t num_transmissions = 0;
    uint32_t packet_errors = 0;
    
    while (num_transmissions < MAX_TRANSMISSIONS_LF && packet_errors < TARGET_ERRORS) {
        static float bpsk_symbols[32];
        static ldpc_quantized_t quantized_symbols[32];

        uint8_t msg = rand_u8_range(1, 0x3F);
        uint32_t lf_e = len_field_encode(msg);

        bytes_to_bpsk((uint8_t *)&lf_e, 4, bpsk_symbols);
        add_awgn(bpsk_symbols, 32, noise_stddev);
        
        uint32_t lf_d = 0;
        for (int i = 0; i < 32; i++) {
            lf_d |= (bpsk_symbols[i] < 0) << i;
        }

        lf_d = len_field_decode(lf_d);

        if ((lf_e ^ lf_d) != 0) {
            packet_errors++;
        }

        num_transmissions++;
    }

    return (float)packet_errors/num_transmissions;
}

float find_packet_error_rate_lf_soft (float noise_stddev) {
    uint32_t num_transmissions = 0;
    uint32_t packet_errors = 0;
    
    while (num_transmissions < MAX_TRANSMISSIONS_LF && packet_errors < TARGET_ERRORS) {
        static float bpsk_symbols[32];
        static ldpc_quantized_t quantized_symbols[32];

        uint8_t msg = rand_u8_range(1, 0x3F);
        uint32_t lf_e = len_field_encode(msg);

        bytes_to_bpsk((uint8_t *)&lf_e, 4, bpsk_symbols);
        add_awgn(bpsk_symbols, 32, noise_stddev);
        quantize(bpsk_symbols, 32, quantized_symbols);

        uint32_t lf_d = len_field_decode_soft(quantized_symbols);

        if ((lf_e ^ lf_d) != 0) {
            packet_errors++;
        }

        num_transmissions++;
    }

    return (float)packet_errors/num_transmissions;
}