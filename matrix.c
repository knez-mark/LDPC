#include "matrix.h"
#include <string.h>

uint16_t find_vector_weight (vector_t x, uint16_t len) {
    uint16_t weight = 0;

    uint16_t bytes = len / 8;
    uint16_t rem = len % 8;

    for (uint16_t i = 0; i < bytes; i++) {
        weight += __builtin_popcount(x.data8[i]);
    }

    if (rem) {
        uint8_t last = x.data8[bytes] & ((1 << rem) - 1);
        weight += __builtin_popcount(last);
    }

    return weight;
}

void vector_add2 (uint8_t *x, uint8_t *y, uint8_t *result, uint16_t length) {
    uint16_t bytes = length / 8;
    uint16_t rem = length % 8;
    
    for (int i = 0; i < bytes; i++) {
        result[i] = x[i] ^ y[i];
    }

    if (rem) {
        result[bytes] = (x[bytes] ^ y[bytes]) & ((1 << rem) - 1);
    }
}

void vector_add4 (uint8_t *w, uint8_t *x, uint8_t *y, uint8_t *z, uint8_t *result, uint16_t length) {
    uint16_t bytes = length / 8;
    uint16_t rem = length % 8;
    
    for (int i = 0; i < bytes; i++) {
        result[i] = w[i] ^ x[i] ^ y[i] ^ z [i];
    }

    if (rem) {
        result[bytes] = (w[bytes] ^ x[bytes] ^ y[bytes] ^ z[bytes]) & ((1 << rem) - 1);
    }
}

void get_byte_aligned(
    uint8_t *input,
    uint8_t *output,
    uint16_t len_bits,
    uint16_t lifting_size
)
{
    uint16_t bytes_per_block =
        (lifting_size + 7) >> 3;

    uint16_t n_blocks =
        (len_bits + lifting_size - 1)
        / lifting_size;

    uint32_t output_size =
        n_blocks * bytes_per_block;

    memset(output, 0, output_size);

    for (uint32_t bit_idx = 0;
         bit_idx < len_bits;
         bit_idx++)
    {
        //Read input bit

        uint32_t in_byte = bit_idx >> 3;
        uint32_t in_bit  = (bit_idx & 7);

        uint8_t bit =
            (input[in_byte] >> in_bit) & 0x01;

        //Compute output location

        uint32_t block =
            bit_idx / lifting_size;

        uint32_t block_bit =
            bit_idx % lifting_size;

        uint32_t out_byte =
            block * bytes_per_block
            + (block_bit >> 3);

        uint32_t out_bit = (block_bit & 7);

        //Write bit

        output[out_byte] |=
            bit << out_bit;
    }
}

void get_packed(
    uint8_t *input,
    uint8_t *output,
    uint16_t len_bits,
    uint16_t lifting_size
)
{
    uint16_t bytes_per_block =
        (lifting_size + 7) >> 3;

    uint16_t n_blocks =
        (len_bits + lifting_size - 1)
        / lifting_size;

    uint32_t output_size =
        (len_bits + 7) >> 3;

    memset(output, 0, output_size);

    for (uint32_t bit_idx = 0;
         bit_idx < len_bits;
         bit_idx++)
    {
        //Find source location

        uint32_t block =
            bit_idx / lifting_size;

        uint32_t block_bit =
            bit_idx % lifting_size;

        uint32_t in_byte =
            block * bytes_per_block
            + (block_bit >> 3);

        uint32_t in_bit = (block_bit & 7);

        uint8_t bit =
            (input[in_byte] >> in_bit) & 0x01;

        //Write packed output

        uint32_t out_byte =
            bit_idx >> 3;

        uint32_t out_bit = (bit_idx & 7);

        output[out_byte] |=
            bit << out_bit;
    }
}

static void circular_shift_32 (
    uint32_t *input,
    uint32_t *output,
    uint16_t shift
)
{
    shift &= 31;
    *output = (*input >> shift) | 
                (*input << (32 - shift));
}

static void circular_shift_16 (
    uint16_t *input,
    uint16_t *output,
    uint16_t shift
)
{
    shift &= 15;
    *output = (*input >> shift) | 
                (*input << (16 - shift));
}

static void circular_shift_8 (
    uint8_t *input,
    uint8_t *output,
    uint16_t shift
)
{
    shift &= 7;
    *output = (*input >> shift) | 
                (*input << (8 - shift));
}

void circular_shift_mult_of_byte(
    uint8_t *input,
    uint8_t *output,
    uint16_t lifting_size,
    uint16_t shift
)
{
    shift %= lifting_size;

    const uint16_t num_bytes =
        (lifting_size + 7) >> 3;

    const uint16_t byte_shift =
        shift >> 3;

    const uint8_t bit_shift =
        shift & 7;


    for (uint16_t i = 0;
         i < num_bytes;
         i++)
    {
        uint16_t src0 =
            (i + byte_shift)
            % num_bytes;

        uint16_t src1 =
            (src0 + 1)
            % num_bytes;

        output[i] =
            (input[src0]
                >> bit_shift)
            |
            (input[src1]
                << (8 - bit_shift));
    }
}

static void circular_shift_generic(
    uint8_t *input,
    uint8_t *output,
    uint16_t lifting_size,
    uint16_t shift
)
{
    shift %= lifting_size;

    const uint16_t num_bytes =
        (lifting_size + 7) >> 3;

    const uint16_t byte_shift =
        shift >> 3;

    const uint8_t bit_shift =
        shift & 0x07;

    memset(output, 0, num_bytes);

    for (uint16_t dst_bit = 0;
         dst_bit < lifting_size;
         dst_bit++)
    {
        //Left circular shift
        uint16_t src_bit =
            (dst_bit
             + lifting_size
             + shift)
            % lifting_size;

        //Read source bit
        uint8_t bit =
            (input[src_bit >> 3]
             >> (src_bit & 7))
            & 1;

        //Write destination bit
        output[dst_bit >> 3]
            |= bit
               << (dst_bit & 7);
    }
}

void circular_shift(
    uint8_t *input,
    uint8_t *output,
    uint16_t lifting_size,
    uint16_t shift
)
{
    if (lifting_size == 32) {
        circular_shift_32 ((uint32_t*)input, (uint32_t*)output, shift);
    }
    else if (lifting_size == 16) {
        circular_shift_16 ((uint16_t*)input, (uint16_t*)output, shift);
    }
    else if (lifting_size == 8) {
        circular_shift_8 ((uint8_t*)input, (uint8_t*)output, shift);
    }
    else if ((lifting_size % 8) == 0) {
        circular_shift_mult_of_byte (input, output, lifting_size, shift);
    }
    else {
        circular_shift_generic (input, output, lifting_size, shift);
    }
}

void circular_matrix_multiply (quasi_cyclic_matrix_t * H, uint8_t *x, uint8_t *b, uint8_t *temp, uint8_t lifting_size) {

    uint16_t block_size = (lifting_size + 7) / 8;

    for (int j = 0; j < H->rows; j++) {

        const uint8_t *col = &H->columnIndexMap[H->rowOffset[j]];
        const uint16_t *base_graph = &H->base_graph[H->rowOffset[j]];

        for (int i = 0; i < H->rowWeight[j]; i++) {
            
            if (col[i] >= H->cols) {
                break;
            }

            circular_shift (x + col[i]*block_size, temp, lifting_size, base_graph[i]);
            vector_add2 (temp, b + j*block_size, b + j*block_size, lifting_size);
        }
    }
}