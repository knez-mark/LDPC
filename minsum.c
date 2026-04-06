#include "minsum.h"
#include <string.h>
#include <stdio.h>

static float_struct R_mj [NUM_EDGES] = {0};

#include <math.h>

static inline float boxplus(float a, float b)
{
    float sign = ((a >= 0) == (b >= 0)) ? 1.0f : -1.0f;

    float abs_a = fabsf(a);
    float abs_b = fabsf(b);

    float min_ab = (abs_a < abs_b) ? abs_a : abs_b;

    float term1 = log1pf(expf(-fabsf(a + b)));
    float term2 = log1pf(expf(-fabsf(a - b)));

    return sign * (min_ab + term1 - term2);
}

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

void layered_normalized_minsum (float * Lq, uint16_t len, quasi_cyclic_matrix_t * Hm, quasi_cyclic_matrix_t * Hp) {

	uint8_t lifting_size = len/(Hm->cols + Hp->cols);

    //for (int i = NUM_EDGES - 3; i < NUM_EDGES; i++) {
        //for (int j = 0; j < lifting_size; j++) {
            //if (R_mj[i].element[j] != 0) {
                //printf ("i=%d, j=%d, %f\n", i, j, R_mj[i].element[j]);
            //}
        //}
    //}
    //printf ("\n");

	for (int m = 0; m < Hm->rows; m++) {
		static uint16_t columnOffset [MAX_ROW_WEIGHT]= {0};
		static uint16_t circular_shift [MAX_ROW_WEIGHT]= {0};

        uint8_t rowWeight_msg = 0, rowWeight_parity = 0, total_rowWeight = 0;

		for (int k = 0; k < Hm->rowWeight[m]; k++) {
            if (Hm->columnIndexMap[Hm->rowOffset[m]+k] >= Hm->cols) {
                break;
            }

		    columnOffset[k] = lifting_size*Hm->columnIndexMap[Hm->rowOffset[m]+k];
		    circular_shift[k] = Hm->base_graph[Hm->rowOffset[m]+k];
        
            rowWeight_msg++;
		}
		for (int k = 0; k < Hp->rowWeight[m]; k++) {
            if (Hp->columnIndexMap[Hp->rowOffset[m]+k] >= Hp->cols) {
                break;
            }

		    columnOffset[rowWeight_msg+k] = lifting_size*(Hm->cols + Hp->columnIndexMap[Hp->rowOffset[m]+k]);
		    circular_shift[rowWeight_msg+k] = Hp->base_graph[Hp->rowOffset[m]+k];
        
            rowWeight_parity++;
		}        
        total_rowWeight = rowWeight_msg + rowWeight_parity;

    	for (int j = 0; j < lifting_size; j++) {
            float Lq_mj [MAX_ROW_WEIGHT];
            for (int i = 0; i < total_rowWeight; i++) {
                //Lq_mj[i] = Lq[lifting_size*H->columnIndexMap[H->rowOffset[m]+i]+(H->base_graph[H->rowOffset[m]+i]+j)%lifting_size] - R_mj[H->rowOffset[m]+i].element[j];
                Lq_mj[i] = Lq[columnOffset[i]+((circular_shift[i]+j)%lifting_size)] - R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j];
            }

            float f1, f2;
            uint8_t i1 = 0, i2 = 0;
            min2 (Lq_mj, total_rowWeight, &f1, &f2, &i1, &i2);

            if (f1 == 0 && f2 == 0) {
                for (int i = 0; i < total_rowWeight; i++) {
                    R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j] = 0;
                }
            }
            else {
                float temp [MAX_ROW_WEIGHT];
                float prod = 1;
                if (f1 == 0) {
                    for (int i = 0; i < total_rowWeight; i++) {
                    	temp [i] = 0.0f;
                    	if (i == i1) {
                            continue;
                        }
                        prod *= sign(Lq_mj [i]);
                    }
                    temp [i1] = prod * f2;
                }
                else {
                    for (int i = 0; i < total_rowWeight; i++) {
                        prod *= sign(Lq_mj [i]);
                    }
                    for (int i = 0; i < total_rowWeight; i++) {
                        if (i == i1) {
                            temp [i] = prod * sign (Lq_mj [i]) * f2;
                        }
                        else {
                            temp [i] = prod * sign (Lq_mj [i]) * f1;
                        }
                    }
                }
                for (int i = 0; i < total_rowWeight; i++) {
                    R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j] = ALPHA * temp [i];
                }
            }

            for (int i = 0; i < total_rowWeight; i++) {
            	//Lq[lifting_size*H->columnIndexMap[H->rowOffset[m]+i]+(H->base_graph[H->rowOffset[m]+i]+j)%lifting_size] = Lq_mj[i] + R_mj[H->rowOffset[m]+i].element[j];
            	Lq[columnOffset[i]+((circular_shift[i]+j)%lifting_size)] = Lq_mj[i] + R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j];
            }
        }
    }
}

void layered_sum_product (float * Lq, uint16_t len, quasi_cyclic_matrix_t * Hm, quasi_cyclic_matrix_t * Hp) {

	uint8_t lifting_size = len/(Hm->cols + Hp->cols);

	for (int m = 0; m < Hm->rows; m++) {
		static uint16_t columnOffset [MAX_ROW_WEIGHT]= {0};
		static uint16_t circular_shift [MAX_ROW_WEIGHT]= {0};

        uint8_t rowWeight_msg = 0, rowWeight_parity = 0, total_rowWeight = 0;

		for (int k = 0; k < Hm->rowWeight[m]; k++) {
            if (Hm->columnIndexMap[Hm->rowOffset[m]+k] >= Hm->cols) {
                break;
            }

		    columnOffset[k] = lifting_size*Hm->columnIndexMap[Hm->rowOffset[m]+k];
		    circular_shift[k] = Hm->base_graph[Hm->rowOffset[m]+k];
        
            rowWeight_msg++;
		}
		for (int k = 0; k < Hp->rowWeight[m]; k++) {
            if (Hp->columnIndexMap[Hp->rowOffset[m]+k] >= Hp->cols) {
                break;
            }

		    columnOffset[rowWeight_msg+k] = lifting_size*(Hm->cols + Hp->columnIndexMap[Hp->rowOffset[m]+k]);
		    circular_shift[rowWeight_msg+k] = Hp->base_graph[Hp->rowOffset[m]+k];
        
            rowWeight_parity++;
		}        
        total_rowWeight = rowWeight_msg + rowWeight_parity;

    	for (int j = 0; j < lifting_size; j++) {
            float Lq_mj [MAX_ROW_WEIGHT];
            for (int i = 0; i < total_rowWeight; i++) {
                //Lq_mj[i] = Lq[lifting_size*H->columnIndexMap[H->rowOffset[m]+i]+(H->base_graph[H->rowOffset[m]+i]+j)%lifting_size] - R_mj[H->rowOffset[m]+i].element[j];
                Lq_mj[i] = Lq[columnOffset[i]+((circular_shift[i]+j)%lifting_size)] - R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j];
            }

            float prefix[MAX_ROW_WEIGHT];
            float suffix[MAX_ROW_WEIGHT];

            prefix[0] = Lq_mj[0];
            for (int i = 1; i < total_rowWeight; i++) {
                prefix[i] = boxplus(prefix[i-1], Lq_mj[i]);
            }
            suffix[total_rowWeight-1] = Lq_mj[total_rowWeight-1];
            for (int i = total_rowWeight-2; i >= 0; i--) {
                suffix[i] = boxplus (Lq_mj[i], suffix[i+1]);
            }

            for (int i = 0; i < total_rowWeight; i++) {
                if (i == 0) {
                    R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j] = suffix[1];
                }
                else if (i == total_rowWeight-1) {
                    R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j] = prefix[total_rowWeight-2];
                }
                else {
                    R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j] = boxplus(prefix[i-1], suffix[i+1]);
                }
            }

            for (int i = 0; i < total_rowWeight; i++) {
            	//Lq[lifting_size*H->columnIndexMap[H->rowOffset[m]+i]+(H->base_graph[H->rowOffset[m]+i]+j)%lifting_size] = Lq_mj[i] + R_mj[H->rowOffset[m]+i].element[j];
            	Lq[columnOffset[i]+((circular_shift[i]+j)%lifting_size)] = Lq_mj[i] + R_mj[Hm->rowOffset[m]+Hp->rowOffset[m]+i].element[j];
            }
        }
    }
}