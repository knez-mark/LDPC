#ifndef FRAME_H_
#define FRAME_H_

#include <stdint.h>

#define MSG_SIZE_MAX 88

uint16_t get_LDPC_len (uint32_t len_field);
uint32_t generate_length_field(uint8_t msg_size);
void format_encoded_data(
    uint8_t *msg,
    uint8_t *parity,
    uint8_t *encoded_msg,
    uint32_t length_field
);

#endif /* FRAME_H_ */