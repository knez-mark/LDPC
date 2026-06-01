#ifndef LDPC_MEMORY_H_
#define LDPC_MEMORY_H_

#include "LDPC_impl.h"

ldpc_encoder_t* LDPC_encoder_alloc ();
ldpc_decoder_t* LDPC_decoder_alloc ();

void LDPC_encoder_free (ldpc_encoder_t* encoder);
void LDPC_decoder_free (ldpc_decoder_t* decoder);

#endif /* LDPC_MEMORY_H_ */