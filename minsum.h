#ifndef MINSUM_H_
#define MINSUM_H_

#include <stdint.h>
#include <float.h>
#include <math.h>
#include "matrix.h"
#include "LDPC.h"

void layered_normalized_minsum (float * Lq, uint16_t len, quasi_cyclic_matrix_t * H);

#endif /* MINSUM_H_ */