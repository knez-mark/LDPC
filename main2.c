#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include "matrix.h"
#include "simulation.h"
#include "length_field.h"
#include "frame.h"
#include "LDPC.h"

static void print_5bits(uint32_t value)
{
    for (int i = 4; i >= 0; i--) {
        printf("%u", (value >> i) & 1);
    }
}

int main() {
    LDPC_memory_init ();
    
    float noise_stddev = 0.6f;

    srand((unsigned)time(NULL));

    /*
    //Static, zero-initialized struct
    static message_t message = {0};

    static encoded_t encoding = {0};
    static encoded_t decoding = {0};

    static uint8_t parity [272] = {0};

    static float bpsk_symbols [32+128*8] = {0};
    
    static float ldpc_llr [68*32] = {0};

    //Generate random message size: 0–87
    message.msg_size = 21;//rand_u8_range(0, 87);

    //Fill message payload
    generate_random_bytes(message.msg, message.msg_size);

    //Generate length field
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

    for (int i = 0; i < get_frame_len (encoding.length_field) - 4; i++) {
        printf ("%x, ", encoding.ldpc [i]);
        if ((i+1) % 8 == 0) {
            printf ("\n");
        }
    }

    uint16_t num_bytes = get_frame_len (encoding.length_field);
    uint16_t num_bits = bytes_to_bpsk((uint8_t*) &encoding, num_bytes, bpsk_symbols);
    add_awgn(bpsk_symbols, num_bits, noise_stddev);

    for (int i = 32; i < num_bits; i++) {
        printf ("%5.2f, ", bpsk_symbols [i]);
        if ((i+1) % 8 == 0) {
            printf ("\n");
        }
    }

    decoding.length_field = len_field_decode_soft (bpsk_symbols);

    printf("length_field  = %x\n", decoding.length_field);
    printf ("\n");

    format_decoded_data(bpsk_symbols + 32, ldpc_llr, decoding.length_field);

    
    for (int i = 0; i < get_LDPC_len(encoding.length_field); i++) {
        printf ("%5.2f, ", ldpc_llr [i]);
        if ((i+1) % 8 == 0) {
            printf ("\n");
        }
    }

    uint16_t num_iters;
    printf ("parity_check_errors: %d\n", LDPC_decode (ldpc_llr, get_LDPC_len(decoding.length_field), decoding.ldpc, 50, &num_iters));
    printf ("num_iters: %d\n", num_iters);

    printf ("LDPC decoded: \n");

    for (int i = 0; i < message.msg_size + 1; i++) {
        printf ("%x, ", decoding.ldpc[i]);
        if ((i+1) % 8 == 0) {
            printf ("\n");
        }
    }
    */

    struct timespec start, end;
    FILE *fp = fopen("sim_data.csv", "w"); 

    for (uint8_t msg_len = 0; msg_len < 88; msg_len ++) {

        float BER_vals [65] = {0};
        int i = 0;
        printf ("MESSAGE LENGTH = %d\n\n", msg_len);

        for (float EbN0_dB = 0; EbN0_dB <= 16; EbN0_dB += 0.25, i++) {
            clock_gettime(CLOCK_MONOTONIC, &start);

            float sigma = sqrt(1.0 / (2.0 * pow(10.0, EbN0_dB / 10.0))); 
            float BER = find_bit_error_rate3 (msg_len, sigma);

            clock_gettime(CLOCK_MONOTONIC, &end);

            double elapsed =
                (end.tv_sec - start.tv_sec) +
                (end.tv_nsec - start.tv_nsec) / 1e9;

            printf ("Elapsed time: %.2f seconds\n", elapsed);
            printf ("BER = %.2e at Eb/N0 = %.2f dB (sigma = %.5f)\n", BER, EbN0_dB, sigma);
            printf ("\n");

            BER_vals [i] = BER;

            if (BER == (float)0.0){
                break;
            }
        }
        fprintf(fp, "%d,", msg_len);
        for (i = 0; i < 65; i++) {
            fprintf(fp, "%.6e", BER_vals [i]);
            if (i < 64) {
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}