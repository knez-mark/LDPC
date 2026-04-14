#include "LDPC_memory.h"
#include "matrix.h"
#include <stdbool.h>
#include <stddef.h>

static ldpc_encoder_t encoders;
static ldpc_decoder_t decoders;

static bool enc_in_use;
static bool dec_in_use;

void LDPC_memory_init (void) {
    enc_in_use = false;

    dec_in_use = false;
}

ldpc_encoder_t* LDPC_encoder_alloc () {
    if (enc_in_use) return NULL;
    enc_in_use = 1;
    return &encoders;
}

ldpc_decoder_t* LDPC_decoder_alloc () {
    if (dec_in_use) return NULL;
    dec_in_use = 1;
    return &decoders;
}

void LDPC_encoder_free (ldpc_encoder_t* encoder) {
    if (encoder == &encoders) {
        enc_in_use = 0;
    }
}

void LDPC_decoder_free (ldpc_decoder_t* decoder) {
    if (decoder == &decoders) {
        dec_in_use = 0;
    }
}


