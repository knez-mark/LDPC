#ifndef LENGTH_FIELD_H_
#define LENGTH_FIELD_H_

#include <stdint.h>

uint32_t len_field_encode(uint8_t msg);
uint32_t len_field_decode(uint32_t codeword);
uint32_t len_field_decode_soft(float* llr);

#endif /* LENGTH_FIELD_H_ */