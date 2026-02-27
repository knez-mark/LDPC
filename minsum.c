#include "minsum.h"
#include <string.h>

static float_struct R_mj [NUM_EDGES] = {0};

void reset_minsum () {
    memset (R_mj, 0, sizeof(R_mj));
}

static float sign (float in) {
    if (in > 0) {
        return 1;
    }
    else if (in < 0) {
        return -1;
    }
    else {
        return 0;
    }
}

static void min2 (float * arr, uint8_t size, float * v1, float * v2, uint8_t * i1, uint8_t * i2) {
    
    *v1 = FLT_MAX;
    *v2 = FLT_MAX;

    for (int i = 0; i < size; i++) {
        if (fabs(arr[i])< *v1) {
            *v2 = *v1;
            *i2 = *i1;
            *v1 = fabs(arr[i]);
            *i1 = i;
        }
        else if (fabs(arr[i])< *v2) {
            *v2 = fabs(arr[i]);
            *i2 = i;
        }
    }
}

void layered_normalized_minsum (float * Lq, uint16_t len, quasi_cyclic_matrix_t * H) {

	uint8_t lifting_size = len/H->cols;

	for (int m = 0; m < BG1_ROWS; m++) {
		static uint16_t temp_arr [MAX_ROW_WEIGHT]= {0};
		static uint16_t temp_arr_2 [MAX_ROW_WEIGHT]= {0};
		for (int k = 0; k < H->rowWeight[m]; k++) {
			temp_arr[k] = lifting_size*H->columnIndexMap[H->rowOffset[m]+k];
		  temp_arr_2[k] = H->base_graph[H->rowOffset[m]+k];
		}
    	for (int j = 0; j < lifting_size; j++) {
            float Lq_mj [MAX_ROW_WEIGHT];
            for (int i = 0; i < H->rowWeight[m]; i++) {
                //Lq_mj[i] = Lq[lifting_size*H->columnIndexMap[H->rowOffset[m]+i]+(H->base_graph[H->rowOffset[m]+i]+j)%lifting_size] - R_mj[H->rowOffset[m]+i].element[j];
                Lq_mj[i] = Lq[temp_arr[i]+((temp_arr_2[i]+j)%lifting_size)] - R_mj[H->rowOffset[m]+i].element[j];
            }

            float f1, f2;
            uint8_t i1 = 0, i2 = 0;
            min2 (Lq_mj, H->rowWeight[m], &f1, &f2, &i1, &i2);

            if (f1 == 0 && f2 == 0) {
                for (int i = 0; i < H->rowWeight[m]; i++) {
                    R_mj[H->rowOffset[m]+i].element[j] = 0;
                }
            }
            else {
                float temp [MAX_ROW_WEIGHT];
                float prod = 1;
                if (f1 == 0) {
                    for (int i = 0; i < H->rowWeight[m]; i++) {
                    	temp [i] = 0.0f;
                    	if (i == i1) {
                            continue;
                        }
                        prod *= sign(Lq_mj [i]);
                    }
                    temp [i1] = prod * f2;
                }
                else {
                    for (int i = 0; i < H->rowWeight[m]; i++) {
                        prod *= sign(Lq_mj [i]);
                    }
                    for (int i = 0; i < H->rowWeight[m]; i++) {
                        if (i == i1) {
                            temp [i] = prod * sign (Lq_mj [i]) * f2;
                        }
                        else {
                            temp [i] *= prod * sign (Lq_mj [i]) * f1;
                        }
                    }
                }
                for (int i = 0; i < H->rowWeight[m]; i++) {
                    R_mj[H->rowOffset[m]+i].element[j] = ALPHA * temp [i];
                }
            }

            for (int i = 0; i < H->rowWeight[m]; i++) {
            	//Lq[lifting_size*H->columnIndexMap[H->rowOffset[m]+i]+(H->base_graph[H->rowOffset[m]+i]+j)%lifting_size] = Lq_mj[i] + R_mj[H->rowOffset[m]+i].element[j];
            	Lq[temp_arr[i]+((temp_arr_2[i]+j)%lifting_size)] = Lq_mj[i] + R_mj[H->rowOffset[m]+i].element[j];
            }
        }
    }

    //for (int i = 0; i < NUM_EDGES; i++) {
    //    for (int j = 0; j < lifting_size; j++) {
    //        if (R_mj[i].element[j] != 0) {
                //printf ("i=%d, j=%d, %f\n", i, j, R_mj[i].element[j]);
    //        }
    //    }
    //}
    //printf ("\n");

}