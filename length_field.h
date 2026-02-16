#ifndef LENGTH_FIELD_H_
#define LENGTH_FIELD_H_

#include <stdint.h>

uint32_t len_field_encode(uint8_t msg);
uint8_t len_field_decode(uint32_t codeword);

#endif /* LENGTH_FIELD_H_ */