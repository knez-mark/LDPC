#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "matrix.h"
#include "simulation.h"
#include "length_field.h"
#include "frame.h"
#include "LDPC.h"

typedef struct {
    uint8_t  msg_size;
    uint8_t  msg[255];
} message_t;

typedef struct {
    uint32_t length_field;
    uint8_t  ldpc[128];
} encoded_t;

static void print_5bits(uint32_t value)
{
    for (int i = 4; i >= 0; i--) {
        printf("%u", (value >> i) & 1);
    }
}

int main() {
    LDPC_init ();
    
    srand((unsigned)time(NULL));

    /* Static, zero-initialized struct */
    static message_t message = {0};

    static encoded_t encoding = {0};

    static uint8_t parity [48] = {0};

    /* Generate random message size: 0–87 */
    message.msg_size = rand_u8_range(0, 87);

    /* Fill message payload */
    generate_random_bytes(message.msg, message.msg_size);

    /* Generate length field */
    encoding.length_field = generate_length_field(message.msg_size + 1);

    LDPC_encode ((uint8_t *) &message, get_LDPC_len(encoding.length_field), parity);

    format_encoded_data((uint8_t *) &message, parity, encoding.ldpc, encoding.length_field);

    printf("msg_size      = %d\n", message.msg_size);
    printf("length_field  = %x\n", encoding.length_field);
    print_5bits(encoding.length_field);
    printf("\n");

    for (int i = 0; i < message.msg_size + 1; i++) {
        if (i == 0) {
            printf ("%x, ", message.msg_size);
        }
        else {
            printf ("%x, ", message.msg [i-1]);
        }
        if ((i+1) % 8 == 0) {
            printf ("\n");
        }
    }

    printf ("\n\n");

    for (int i = 0; i < 32*(1 + ((encoding.length_field >> 3) & 0x03)); i++) {
        printf ("%x, ", encoding.ldpc [i]);
        if ((i+1) % 8 == 0) {
            printf ("\n");
        }
    }

    return 0;
}