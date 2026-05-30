#ifndef LDPC_H_
#define LDPC_H_

#include <stdint.h>

#define LDPC_ENCODE_ENABLE 1
#define LDPC_DECODE_ENABLE 1
#define LDPC_BG1_ENABLE 1
#define LDPC_BG2_ENABLE 0

#define USE_SUM_PRODUCT 0

#define MAX_LIFTING_SIZE 384

//#define NUM_ENCODERS 0
//#define NUM_DECODERS 1

typedef enum {
    LDPC_LIFTING_AUTO,
    LDPC_LIFTING_EXPLICIT
} lifting_mode_t;

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

    lifting_mode_t lifting_mode;
    uint16_t lifting_size;
} ldpc_decoder_cfg_t;

typedef int8_t ldpc_quantized_t;

typedef struct ldpc_encoder ldpc_encoder_t;
typedef struct ldpc_decoder ldpc_decoder_t;

void LDPC_memory_init (void);

ldpc_encoder_t* LDPC_encoder_create (ldpc_encoder_cfg_t cfg);
ldpc_decoder_t* LDPC_decoder_create (ldpc_decoder_cfg_t cfg);

void LDPC_encoder_destroy (ldpc_encoder_t* ldpc);
void LDPC_decoder_destroy (ldpc_decoder_t* ldpc);

ldpc_encoder_cfg_t LDPC_get_encoder_config (ldpc_encoder_t* ldpc);
ldpc_decoder_cfg_t LDPC_get_decoder_config (ldpc_decoder_t* ldpc);

uint8_t LDPC_set_encoder_config (ldpc_encoder_t* ldpc, ldpc_encoder_cfg_t cfg);
uint8_t LDPC_set_decoder_config (ldpc_encoder_t* ldpc, ldpc_decoder_cfg_t cfg);

uint8_t LDPC_encode (ldpc_encoder_t* ldpc, uint8_t * data, uint8_t * parity);
uint16_t LDPC_decode (ldpc_decoder_t* ldpc, void * Lq, uint8_t * decoded, uint16_t* num_iters);

#endif /* LDPC_H_ */