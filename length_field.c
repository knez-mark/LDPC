#include <stdint.h>
#include "matrix.h"

//All 64 possible codewords for a [32, 6] linear block code
uint32_t codewords[64] = {
    0x00000000, 0x07925AB1, 0x095566A9, 0x0EC73C18, 0x12A8B162, 0x153AEBD3, 0x1BFDD7CB, 0x1C6F8D7A,
    0x21CC4D92, 0x265E1723, 0x28992B3B, 0x2F0B718A, 0x3364FCF0, 0x34F6A641, 0x3A319A59, 0x3DA3C0E8,
    0x432A9364, 0x44B8C9D5, 0x4A7FF5CD, 0x4DEDAF7C, 0x51822206, 0x561078B7, 0x58D744AF, 0x5F451E1E,
    0x62E6DEF6, 0x65748447, 0x6BB3B85F, 0x6C21E2EE, 0x704E6F94, 0x77DC3525, 0x791B093D, 0x7E89538C,
    0x82954695, 0x85071C24, 0x8BC0203C, 0x8C527A8D, 0x903DF7F7, 0x97AFAD46, 0x9968915E, 0x9EFACBEF,
    0xA3590B07, 0xA4CB51B6, 0xAA0C6DAE, 0xAD9E371F, 0xB1F1BA65, 0xB663E0D4, 0xB8A4DCCC, 0xBF36867D,
    0xC1BFD5F1, 0xC62D8F40, 0xC8EAB358, 0xCF78E9E9, 0xD3176493, 0xD4853E22, 0xDA42023A, 0xDDD0588B,
    0xE0739863, 0xE7E1C2D2, 0xE926FECA, 0xEEB4A47B, 0xF2DB2901, 0xF54973B0, 0xFB8E4FA8, 0xFC1C1519,
};

uint32_t len_field_encode(uint8_t msg)
{
    if (msg == 0 || msg >= 64) return 0;
    return codewords[msg];
}

uint8_t len_field_decode(uint32_t codeword)
{
    uint8_t best = 1;
    uint8_t best_dist = 32;

    for (uint8_t msg = 1; msg < 64; msg++) {
        uint32_t diff = codeword ^ codewords[msg];
        uint8_t d = find_vector_weight((vector_t)&diff, 32);

        if (d < best_dist) {
            best_dist = d;
            best = msg;
            if (d == 0) break;
        }
    }
    return best;
}