#include "LDPC.h"
#include "frame.h"
#include <float.h>

uint16_t get_LDPC_len (uint32_t len_field) {
    uint8_t bits = (len_field >> 3) & 0x03;
    switch (bits) {
        case 0: return 544;
        case 1: return 544*2;
        case 2:
        case 3:
            return 544*4;
        default: return 544*4;
    }
}

static uint16_t get_frame_LDPC_len (uint32_t len_field) {
    uint8_t bits = (len_field >> 3) & 0x03;
    switch (bits) {
        case 0: return 272;
        case 1: return 272*2;
        case 2:
        case 3:
            return 272*4;
        default: return 272*4;
    }
}

uint16_t get_frame_len (uint32_t len_field) {
    return 4 + 32*(1 + ((len_field >> 3) & 0x03));
}

static uint8_t find_last_three_bits(uint8_t remainder)
{
    if (remainder <= 1)  return 0;
    if (remainder <= 4)  return 1;
    if (remainder <= 7)  return 2;
    if (remainder <= 10) return 3;
    if (remainder <= 13) return 4;
    if (remainder <= 16) return 5;
    if (remainder <= 19) return 6;
    return 7;
}

static int map_offset(uint8_t bits) {
    switch (bits) {
        case 22: return 2;
        case 44: return 4;
        case 66:
        case 88:
            return 8;
        default: return 2;
    }
}

static int map_first_two_bits(uint8_t bits) {
    switch (bits) {
        case 0: return 22;
        case 1: return 44;
        case 2: return 66;
        case 3: return 88;
        default: return 0;
    }
}

static int map_last_three_bits(uint8_t bits) {
    switch (bits) {
        case 0: return 1;
        case 1: return 4;
        case 2: return 7;
        case 3: return 10;
        case 4: return 13;
        case 5: return 16;
        case 6: return 19;
        case 7: return 22;
        default: return 0;
    }
}

uint32_t generate_length_field(uint8_t msg_size)
{
    if (msg_size > MSG_SIZE_MAX) {
        return 0;
    }

    uint8_t first_two_bits = (msg_size - 1) / 22;
    if (first_two_bits > 3) {
        first_two_bits = 3;
    }

    uint8_t remainder = msg_size % 22;
    if (remainder == 0) {
        remainder = 22;
    }

    uint8_t last_three_bits = find_last_three_bits(remainder);

    uint32_t length_field_data =
        ((uint32_t)first_two_bits << 3) |
        ((uint32_t)last_three_bits);

    return len_field_encode(length_field_data);
}

void format_encoded_data(
    uint8_t *msg,
    uint8_t *parity,
    uint8_t *encoded_msg,
    uint32_t length_field
) {

    uint8_t base_len = map_first_two_bits((length_field >> 3) & 0x03);
    uint8_t offset   = map_offset(base_len);
    uint8_t end_range = base_len - 22 + map_last_three_bits(length_field & 0x07);

    if ((map_first_two_bits((length_field >> 3) & 0x03) == 22) && (map_last_three_bits(length_field & 0x07) <= 13)) {
        offset = 0;
    }

    /* Determine encoded length (bytes) */
    uint16_t encoded_bytes;
    uint8_t parity_bytes;

    switch (base_len) {
        case 22:
            encoded_bytes = 32;
            parity_bytes = 12;
            break;
        case 44:
            encoded_bytes = 64;
            parity_bytes = 24;
            break;
        case 66:
            encoded_bytes = 96;
            parity_bytes = 36;
            break;
        case 88:
            encoded_bytes = 128;
            parity_bytes = 48;
            break;
        default:
            return;
    }

    /* Clear output buffer */
    for (uint16_t i = 0; i < encoded_bytes; i++) {
        encoded_msg[i] = 0;
    }

    uint16_t msg_idx = offset;
    uint16_t parity_idx = 0;
    uint8_t use_msg = 1;

    for (uint16_t i = 0; i < encoded_bytes; i++) {
        if (use_msg) {
            encoded_msg[i] = msg[msg_idx];

            if (msg_idx == end_range - 1) {
                msg_idx = offset;
                use_msg = 0;
            }
            else {
                msg_idx++;
            }
        } else {
            encoded_msg[i] = parity[parity_idx];

            if (parity_idx == parity_bytes - 1) {
                parity_idx = 0;
                use_msg = 1;
            }
            else {
                parity_idx++;
            }
        }
    }
}

void format_decoded_data(
    float *symbols,
    float *llr,
    uint32_t length_field
) {

    uint8_t base_len = map_first_two_bits((length_field >> 3) & 0x03);
    uint8_t offset   = map_offset(base_len);
    uint8_t end_range = base_len - 22 + map_last_three_bits(length_field & 0x07);

    if ((map_first_two_bits((length_field >> 3) & 0x03) == 22) && (map_last_three_bits(length_field & 0x07) <= 13)) {
        offset = 0;
    }

    /* Determine encoded length (bytes) */
    uint16_t encoded_bytes;
    uint8_t parity_bytes_f; //parity bytes in frame
    uint8_t parity_bytes_l; //parity bytes in LDPC

    switch (base_len) {
        case 22:
            encoded_bytes = 32;
            parity_bytes_f = 12;
            parity_bytes_l = 12;
            break;
        case 44:
            encoded_bytes = 64;
            parity_bytes_f = 24;
            parity_bytes_l = 24;
            break;
        case 66:
            encoded_bytes = 96;
            parity_bytes_f = 36;
            parity_bytes_l = 48;
            break;
        case 88:
            encoded_bytes = 128;
            parity_bytes_f = 48;
            parity_bytes_l = 48;
            break;
        default:
            return;
    }

    uint16_t LPDC_len = get_frame_LDPC_len (length_field);

    /* Clear output buffer */
    for (uint16_t i = 0; i < LPDC_len; i++) {
        llr[i] = 0;
    }

    for (uint16_t i = end_range * 8; i < LPDC_len - parity_bytes_l*8; i++) {
        llr[i] = 10.0f; //Set to large positive value
    }

    uint16_t msg_idx = offset*8;
    uint16_t parity_idx = LPDC_len - parity_bytes_l*8;
    uint8_t use_msg = 1;

    for (uint16_t i = 0; i < encoded_bytes*8; i++) {
        if (use_msg) {
            llr[msg_idx] += symbols[i];

            if (msg_idx == end_range*8 - 1) {
                msg_idx = offset*8;
                use_msg = 0;
            }
            else {
                msg_idx++;
            }
        } else {
            llr[parity_idx] += symbols[i];

            if (parity_idx == LPDC_len - (parity_bytes_l - parity_bytes_f)*8 - 1) {
                parity_idx = LPDC_len - parity_bytes_l*8;
                use_msg = 1;
            }
            else {
                parity_idx++;
            }
        }
    }
}