#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"
#include "simulation.h"
#include "LDPC.h"

int main(int argc, char **argv) {
    
    float noise_stddev = 0.6f;

    srand((unsigned)time(NULL));

    struct timespec start, end;
    FILE *fp = fopen("sim_data.csv", "w"); 

    if (argc >= 2 && (strcmp(argv[1], "lf") == 0)) {
        float PER_vals [65] = {0};
        int i = 0;

        for (float EbN0_dB = -6; EbN0_dB <= 10; EbN0_dB += 0.25, i++) {
            clock_gettime(CLOCK_MONOTONIC, &start);

            float sigma = sqrt(1.0 / (2.0 * pow(10.0, EbN0_dB / 10.0))); 
            float PER = find_packet_error_rate_lf_soft (sigma);

            clock_gettime(CLOCK_MONOTONIC, &end);

            double elapsed =
                (end.tv_sec - start.tv_sec) +
                (end.tv_nsec - start.tv_nsec) / 1e9;

            printf ("Elapsed time: %.2f seconds\n", elapsed);
            printf ("PER = %.2e at Eb/N0 = %.2f dB (sigma = %.5f)\n", PER, EbN0_dB, sigma);
            printf ("\n");

            PER_vals [i] = PER;

            if (PER == (float)0.0){
                break;
            }
        
        }

        for (i = 0; i < 65; i++) {
            fprintf(fp, "%.6e", PER_vals [i]);
            if (i < 64) {
                fprintf(fp, ",");
            }
        }
        fprintf(fp, "\n");
        fclose(fp);

        return 0;
    }

    for (uint8_t msg_len = 0; msg_len < 88; msg_len ++) {

        float BER_vals [65] = {0};
        int i = 0;
        printf ("MESSAGE LENGTH = %d\n\n", msg_len);

        for (float EbN0_dB = 0; EbN0_dB <= 16; EbN0_dB += 0.25, i++) {
            clock_gettime(CLOCK_MONOTONIC, &start);

            float sigma = sqrt(1.0 / (2.0 * pow(10.0, EbN0_dB / 10.0))); 
            float BER = find_bit_error_rate (msg_len, sigma);

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