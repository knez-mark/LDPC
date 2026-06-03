#ifndef LDPC_H_
#define LDPC_H_

#include <stdint.h>
#include <stdbool.h>

#define USE_SUM_PRODUCT 0

#define MAX_LIFTING_SIZE_ENCODE 384
#define MAX_LIFTING_SIZE_DECODE 384

#define NUM_ENCODERS 1
#define NUM_DECODERS 1

typedef enum {
    LDPC_LIFTING_AUTO,
    LDPC_LIFTING_EXPLICIT
} lifting_mode_t;

typedef enum {
    LDPC_CONFIG_FAILURE,
    LDPC_CONFIG_SUCCESS
} ldpc_config_status_t;

typedef struct {
    uint16_t msg_len;
    uint16_t target_code_len;

    uint8_t bgn;

    lifting_mode_t lifting_mode;
    uint16_t lifting_size;
} ldpc_encoder_cfg_t;

typedef struct {
    uint16_t msg_len;
    uint16_t target_code_len;

    uint8_t bgn;

    #if (!USE_SUM_PRODUCT)
        float alpha;
    #endif
    uint16_t max_iters;
    bool use_full_PCM;

    lifting_mode_t lifting_mode;
    uint16_t lifting_size;
} ldpc_decoder_cfg_t;

typedef int8_t ldpc_quantized_t;

typedef struct ldpc_encoder ldpc_encoder_t;
typedef struct ldpc_decoder ldpc_decoder_t;

ldpc_encoder_t* LDPC_encoder_create (ldpc_encoder_cfg_t cfg);
ldpc_decoder_t* LDPC_decoder_create (ldpc_decoder_cfg_t cfg);

void LDPC_encoder_destroy (ldpc_encoder_t* ldpc);
void LDPC_decoder_destroy (ldpc_decoder_t* ldpc);

ldpc_encoder_cfg_t LDPC_get_encoder_config (ldpc_encoder_t* ldpc);
ldpc_decoder_cfg_t LDPC_get_decoder_config (ldpc_decoder_t* ldpc);

ldpc_config_status_t LDPC_set_encoder_config (ldpc_encoder_t* ldpc, ldpc_encoder_cfg_t cfg);
ldpc_config_status_t LDPC_set_decoder_config (ldpc_decoder_t* ldpc, ldpc_decoder_cfg_t cfg);

uint8_t LDPC_encode (ldpc_encoder_t* ldpc, uint8_t * data, uint8_t * parity);
uint16_t LDPC_decode (ldpc_decoder_t* ldpc, void * Lq, uint8_t * decoded, uint16_t* num_iters);

uint32_t len_field_encode (uint8_t msg);
uint32_t len_field_decode (uint32_t codeword);
uint32_t len_field_decode_soft (void* llr);

#endif /* LDPC_H_ */