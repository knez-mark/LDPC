#include "matrix.h"
#include <float.h>

//All 64 possible codewords for a [32, 6] linear block code
static uint32_t codewords[64] = {
    0x00000000, 0xA962A941, 0x26C954C2, 0x8FABFD83, 0x49B23384, 0xE0D09AC5, 0x6F7B6746, 0xC619CE07,
    0x468D1548, 0xEFEFBC09, 0x6044418A, 0xC926E8CB, 0x0F3F26CC, 0xA65D8F8D, 0x29F6720E, 0x8094DB4F,
    0x9566AA90, 0x3C0403D1, 0xB3AFFE52, 0x1ACD5713, 0xDCD49914, 0x75B63055, 0xFA1DCDD6, 0x537F6497,
    0xD3EBBFD8, 0x7A891699, 0xF522EB1A, 0x5C40425B, 0x9A598C5C, 0x333B251D, 0xBC90D89E, 0x15F271DF,
    0x8D5A49E0, 0x2438E0A1, 0xAB931D22, 0x02F1B463, 0xC4E87A64, 0x6D8AD325, 0xE2212EA6, 0x4B4387E7,
    0xCBD75CA8, 0x62B5F5E9, 0xED1E086A, 0x447CA12B, 0x82656F2C, 0x2B07C66D, 0xA4AC3BEE, 0x0DCE92AF,
    0x183CE370, 0xB15E4A31, 0x3EF5B7B2, 0x97971EF3, 0x518ED0F4, 0xF8EC79B5, 0x77478436, 0xDE252D77,
    0x5EB1F638, 0xF7D35F79, 0x7878A2FA, 0xD11A0BBB, 0x1703C5BC, 0xBE616CFD, 0x31CA917E, 0x98A8383F,
};

uint32_t len_field_encode(uint8_t msg)
{
    if (msg == 0 || msg >= 64) return 0;
    return codewords[msg];
}

uint32_t len_field_decode(uint32_t codeword)
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
    return codewords[best];
}

uint32_t len_field_decode_soft(float* llr)
{
    uint8_t best = 1;
    float best_dist = FLT_MAX;

    for (uint8_t msg = 1; msg < 64; msg++) {
        float d = 0;
        for (uint8_t bit = 0; bit < 32; bit++) {
            uint8_t b = (codewords[msg] >> bit) & 1;
            d += b ? llr[bit]: -llr[bit];
        }

        if (d < best_dist) {
            best_dist = d;
            best = msg;
        }
    }
    return codewords[best];
}