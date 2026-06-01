#include "LDPC_memory.h"
#include "matrix.h"
#include <stdbool.h>
#include <stddef.h>

static ldpc_encoder_t encoders [NUM_ENCODERS];
static ldpc_decoder_t decoders [NUM_DECODERS];

static bool enc_in_use [NUM_ENCODERS] = {false};
static bool dec_in_use [NUM_DECODERS] = {false};

ldpc_encoder_t* LDPC_encoder_alloc () {
    for (int i = 0; i < NUM_ENCODERS; i++) {
        if (!enc_in_use[i]) {
            enc_in_use[i] = true;
            return &encoders[i];
        }
    }
    return NULL;
}

ldpc_decoder_t* LDPC_decoder_alloc () {
    for (int i = 0; i < NUM_DECODERS; i++) {
        if (!dec_in_use[i]) {
            dec_in_use [i] = true;
            return &decoders[i];
        }
    }
    return NULL;
}

void LDPC_encoder_free (ldpc_encoder_t* encoder) {
    for (int i = 0; i < NUM_ENCODERS; i++) {
        if (encoder == &encoders[i]) {
            enc_in_use[i] = false;
            return;
        }
    }
}

void LDPC_decoder_free (ldpc_decoder_t* decoder) {
    for (int i = 0; i < NUM_DECODERS; i++) {
        if (decoder == &decoders[i]) {
            dec_in_use[i] = false;
            return;
        }
    }
}
