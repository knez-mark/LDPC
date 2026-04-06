#include "LDPC.h"
#include "matrix.h"

//Original base graph

static const uint8_t base_graph_msg [NUM_EDGES_MSG] = {250, 69, 226, 159, 100, 10, 59, 229, 110, 191, 9, 195, 23, 190, 35, 239, 31, 2, 239, 117, 124, 71, 222, 
                                                       104, 173, 220, 102, 109, 132, 142, 155, 255, 28, 106, 111, 185, 63, 117, 93, 229, 177, 95, 39, 142, 225, 
                                                       225, 245, 205, 251, 117, 121, 89, 84, 20, 150, 131, 243, 136, 86, 246, 219, 211, 240, 76, 244, 144, 12, 
                                                       157, 102, 205, 236, 194, 231, 28, 123, 183, 22, 28, 67, 244, 11, 157, 211, 220, 44, 159, 31, 167, 104, 
                                                       112, 4, 7, 211, 102, 164, 109, 103, 182, 109, 21, 142, 14, 61, 216, 98, 149, 167, 160, 49, 58, 77, 41, 
                                                       83, 182, 78, 160, 42, 21, 32, 234, 7, 177, 248, 151, 185, 206, 55, 206, 127, 16, 229, 40, 96, 65, 63, 75, 
                                                       64, 49, 49, 51, 7, 164, 59, 1, 144, 42, 233, 8, 155, 147, 60, 73, 72, 127, 224, 151, 186, 217, 47, 249, 
                                                       121, 109, 131, 171, 64, 142, 188, 158, 156, 147, 170, 152, 112, 86, 236, 116, 23, 136, 116, 182, 195, 
                                                       243, 215, 61, 25, 104, 194, 128, 165, 181, 63, 86, 236, 84, 216, 73, 120, 95, 177, 221, 112, 199, 2, 187, 
                                                       41, 211, 127, 167, 164, 159, 161, 197, 207, 37, 105, 51, 120, 198, 220, 167, 151, 157, 163, 173, 139, 149, 
                                                       0, 157, 137, 149, 167, 173, 139, 151, 149, 157, 151, 163, 173, 139, 157, 163, 149, 151, 167};

static const uint8_t base_graph_parity [NUM_EDGES_PARITY] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 115, 0, 0, 0, 241, 90, 0, 0, 0, 252, 22, 0, 0, 62, 0, 0, 179, 0, 154, 0, 0, 
                                                             0, 0, 160, 0, 0, 0, 0, 222, 0, 0, 0, 0, 0, 6, 0, 9, 0, 172, 61, 0, 121, 0, 0, 0, 103, 0, 0, 122, 0, 0, 0, 
                                                             0, 0, 137, 0, 139, 0, 173, 0, 0};

static const uint16_t rowOffset_msg [BG1_ROWS] = {0, 17, 33, 50, 67, 69, 75, 83, 89, 96, 104, 110, 115, 121, 125, 131, 136, 140, 145, 150, 155, 159, 164, 168, 172, 176, 
                                                  180, 184, 187, 191, 194, 197, 199, 202, 206, 210, 213, 217, 219, 223, 227, 230, 234, 236, 239, 242};

static const uint16_t rowOffset_parity [BG1_ROWS] = {0, 2, 5, 7, 9, 10, 12, 13, 14, 17, 18, 19, 22, 23, 25, 26, 28, 30, 31, 32, 33, 35, 36, 37, 38, 40, 41, 42, 43, 44, 
                                                     46, 48, 51, 53, 54, 55, 57, 58, 60, 61, 62, 63, 64, 66, 68, 70};

static const uint8_t rowWeight_msg [BG1_ROWS] = {17, 16, 17, 17, 2, 6, 8, 6, 7, 8, 6, 5, 6, 4, 6, 5, 4, 5, 5, 5, 4, 5, 4, 4, 4, 4, 4, 3, 4, 3, 3, 2, 3, 4, 4, 3, 4, 2, 
                                                 4, 4, 3, 4, 2, 3, 3, 3};

static const uint8_t rowWeight_parity [BG1_ROWS] = {2, 3, 2, 2, 1, 2, 1, 1, 3, 1, 1, 3, 1, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 3, 2, 1, 1, 2, 1, 2, 1, 
                                                    1, 1, 1, 2, 2, 2, 1};

static const uint8_t columnIndexMap_msg [NUM_EDGES_MSG] = {0, 1, 2, 3, 5, 6, 9, 10, 11, 12, 13, 15, 16, 18, 19, 20, 21, 0, 2, 3, 4, 5, 7, 8, 9, 11, 12, 14, 15, 16, 
                                                           17, 19, 21, 0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 17, 18, 19, 20, 0, 1, 3, 4, 6, 7, 8, 10, 11, 12, 
                                                           13, 14, 16, 17, 18, 20, 21, 0, 1, 0, 1, 3, 12, 16, 21, 0, 6, 10, 11, 13, 17, 18, 20, 0, 1, 4, 7, 8, 14, 
                                                           0, 1, 3, 12, 16, 19, 21, 0, 1, 10, 11, 13, 17, 18, 20, 1, 2, 4, 7, 8, 14, 0, 1, 12, 16, 21, 0, 1, 10, 11, 
                                                           13, 18, 0, 3, 7, 20, 0, 12, 15, 16, 17, 21, 0, 1, 10, 13, 18, 1, 3, 11, 20, 0, 14, 16, 17, 21, 1, 12, 13, 
                                                           18, 19, 0, 1, 7, 8, 10, 0, 3, 9, 11, 1, 5, 16, 20, 21, 0, 12, 13, 17, 1, 2, 10, 18, 0, 3, 4, 11, 1, 6, 7, 
                                                           14, 0, 2, 4, 15, 1, 6, 8, 0, 4, 19, 21, 1, 14, 18, 0, 10, 13, 1, 7, 0, 12, 14, 1, 2, 11, 21, 0, 7, 15, 17, 
                                                           1, 6, 12, 0, 14, 15, 18, 1, 13, 0, 9, 10, 12, 1, 3, 7, 19, 0, 8, 17, 1, 3, 9, 18, 0, 4, 1, 16, 18, 0, 7, 
                                                           9, 1, 6, 10};

static const uint8_t columnIndexMap_parity [NUM_EDGES_PARITY] = {0, 1, 0, 1, 2, 2, 3, 0, 3, 4, 0, 5, 6, 7, 0, 2, 8, 9, 10, 0, 1, 11, 12, 1, 13, 14, 3, 15, 0, 16, 17, 
                                                                 18, 19, 0, 20, 21, 22, 23, 0, 24, 25, 26, 27, 28, 3, 29, 2, 30, 0, 3, 31, 2, 32, 33, 34, 0, 35, 36, 1, 
                                                                 37, 38, 39, 40, 41, 2, 42, 3, 43, 0, 44, 45};

static quasi_cyclic_matrix_t Hm;
static quasi_cyclic_matrix_t Hp;
static quasi_cyclic_matrix_t A;
static quasi_cyclic_matrix_t C;
static quasi_cyclic_matrix_t D;

quasi_cyclic_matrix_t* get_Hm () {
    return &Hm;
}

quasi_cyclic_matrix_t* get_Hp () {
    return &Hp;
}

quasi_cyclic_matrix_t* get_A () {
    return &A;
}

quasi_cyclic_matrix_t* get_C () {
    return &C;
}

quasi_cyclic_matrix_t* get_D () {
    return &D;
}

void LDPC_init (uint8_t msg_size, uint8_t parity_size) {
    LDPC_encode_init (msg_size, parity_size);
    LDPC_decode_init (msg_size, parity_size);
}

void LDPC_encode_init (uint8_t msg_size, uint8_t parity_size) {
    if (msg_size > 22 || parity_size > 46) {
        return;
    }
    
    A.rows = 4;
    A.cols = msg_size;
    
    A.base_graph = base_graph_msg;
    A.rowOffset = rowOffset_msg;
    A.rowWeight = rowWeight_msg;
    A.columnIndexMap = columnIndexMap_msg;

    C.rows = parity_size - 4;
    C.cols = msg_size;
    
    C.base_graph = base_graph_msg;
    C.rowOffset = rowOffset_msg + 4;
    C.rowWeight = rowWeight_msg + 4;
    C.columnIndexMap = columnIndexMap_msg;

    D.rows = parity_size - 4;
    D.cols = 4;
    
    D.base_graph = base_graph_parity;
    D.rowOffset = rowOffset_parity + 4;
    D.rowWeight = rowWeight_parity + 4;
    D.columnIndexMap = columnIndexMap_parity;
}

void LDPC_decode_init (uint8_t msg_size, uint8_t parity_size) {
    if (msg_size > 22 || parity_size > 46) {
        return;
    }

    Hm.rows = parity_size;
    Hm.cols = msg_size;
    
    Hm.base_graph = base_graph_msg;
    Hm.rowOffset = rowOffset_msg;
    Hm.rowWeight = rowWeight_msg;
    Hm.columnIndexMap = columnIndexMap_msg;

    Hp.rows = parity_size;
    Hp.cols = parity_size;
    
    Hp.base_graph = base_graph_parity;
    Hp.rowOffset = rowOffset_parity;
    Hp.rowWeight = rowWeight_parity;
    Hp.columnIndexMap = columnIndexMap_parity;    
}