#include "LDPC_impl.h"
#include "matrix.h"
#include "base_graph.h"
#include "LDPC_memory.h"
#include <stddef.h>

ldpc_encoder_cfg_t LDPC_get_encoder_config (ldpc_encoder_t* ldpc) {
    if (ldpc != NULL) {
        return ldpc->cfg;
    }
    ldpc_encoder_cfg_t cfg = {0};
    return cfg;
}
ldpc_decoder_cfg_t LDPC_get_decoder_config (ldpc_decoder_t* ldpc) {
    if (ldpc != NULL) {
        return ldpc->cfg;
    }
    ldpc_decoder_cfg_t cfg = {0};
    return cfg;
}

ldpc_encoder_t* LDPC_encoder_create (ldpc_encoder_cfg_t cfg) {
    if (cfg.msg_len > (BG1_COLS - BG1_ROWS)*MAX_LIFTING_SIZE_ENCODE || cfg.target_code_len > BG1_COLS*MAX_LIFTING_SIZE_ENCODE || cfg.msg_len > cfg.target_code_len || cfg.bgn != 1) {
        return NULL;
    }
    
    ldpc_encoder_t* ldpc = LDPC_encoder_alloc ();
    if (ldpc == NULL) return NULL;

    if (cfg.lifting_mode == LDPC_LIFTING_AUTO) {
        cfg.lifting_size = get_nearest_lifting_size (cfg.msg_len, cfg.target_code_len);
        if (cfg.lifting_size == 0 || cfg.lifting_size > MAX_LIFTING_SIZE_ENCODE) {
            LDPC_encoder_free (ldpc);
            return NULL;
        }
    }
    else if (cfg.lifting_mode == LDPC_LIFTING_EXPLICIT) {
        if (!is_valid_lifting_size (cfg.lifting_size) || cfg.lifting_size > MAX_LIFTING_SIZE_ENCODE) {
            LDPC_encoder_free (ldpc);
            return NULL;
        }
    }
    else {
        LDPC_encoder_free (ldpc);
        return NULL;
    }

    int8_t set_index = get_set_index (cfg.lifting_size);
    if (set_index == -1) {
        LDPC_encoder_free (ldpc);
        return NULL;
    }

    const quasi_cyclic_matrix_t * BG1_msg = get_BG1_msg (set_index);
    const quasi_cyclic_matrix_t * BG1_parity = get_BG1_parity (set_index);

    uint16_t msg_size = DIV_CEIL(cfg.msg_len, cfg.lifting_size);
    uint16_t parity_size = DIV_CEIL(cfg.target_code_len - cfg.msg_len, cfg.lifting_size);

    if (msg_size > 22 || parity_size > 46) {
        LDPC_encoder_free (ldpc);
        return NULL;
    }

    ldpc->msg_size = msg_size;
    ldpc->parity_size = parity_size;


    ldpc->A.rows = 4;
    ldpc->A.cols = msg_size;
    
    ldpc->A.base_graph = BG1_msg->base_graph;
    ldpc->A.rowOffset = BG1_msg->rowOffset;
    ldpc->A.rowWeight = BG1_msg->rowWeight;
    ldpc->A.columnIndexMap = BG1_msg->columnIndexMap;

    ldpc->C.rows = parity_size - 4;
    ldpc->C.cols = msg_size;
    
    ldpc->C.base_graph = BG1_msg->base_graph;
    ldpc->C.rowOffset = BG1_msg->rowOffset + 4;
    ldpc->C.rowWeight = BG1_msg->rowWeight + 4;
    ldpc->C.columnIndexMap = BG1_msg->columnIndexMap;

    ldpc->D.rows = parity_size - 4;
    ldpc->D.cols = 4;
    
    ldpc->D.base_graph = BG1_parity->base_graph;
    ldpc->D.rowOffset = BG1_parity->rowOffset + 4;
    ldpc->D.rowWeight = BG1_parity->rowWeight + 4;
    ldpc->D.columnIndexMap = BG1_parity->columnIndexMap;

    ldpc->cfg = cfg;
    
    return ldpc;
}

ldpc_decoder_t* LDPC_decoder_create (ldpc_decoder_cfg_t cfg) {
    if (cfg.msg_len > (BG1_COLS - BG1_ROWS)*MAX_LIFTING_SIZE_DECODE || cfg.target_code_len > BG1_COLS*MAX_LIFTING_SIZE_DECODE || cfg.msg_len > cfg.target_code_len ||cfg.bgn != 1) {
        return NULL;
    }

    ldpc_decoder_t* ldpc = LDPC_decoder_alloc ();
    if (ldpc == NULL) return NULL;

    if (cfg.lifting_mode == LDPC_LIFTING_AUTO) {
        cfg.lifting_size = get_nearest_lifting_size (cfg.msg_len, cfg.target_code_len);
        if (cfg.lifting_size == 0 || cfg.lifting_size > MAX_LIFTING_SIZE_DECODE) {
            LDPC_decoder_free (ldpc);
            return NULL;
        }
    }
    else if (cfg.lifting_mode == LDPC_LIFTING_EXPLICIT){
        if (!is_valid_lifting_size (cfg.lifting_size) || cfg.lifting_size > MAX_LIFTING_SIZE_DECODE) {
            LDPC_decoder_free (ldpc);
            return NULL;
        }
    }
    else {
        LDPC_decoder_free (ldpc);
        return NULL;
    }
    
    int8_t set_index = get_set_index (cfg.lifting_size);
    if (set_index == -1) {
        LDPC_decoder_free (ldpc);
        return NULL;
    }

    const quasi_cyclic_matrix_t * BG1_msg = get_BG1_msg (set_index);
    const quasi_cyclic_matrix_t * BG1_parity = get_BG1_parity (set_index);

    uint16_t msg_size = DIV_CEIL(cfg.msg_len, cfg.lifting_size);
    uint16_t parity_size = DIV_CEIL(cfg.target_code_len - cfg.msg_len, cfg.lifting_size);

    if (msg_size > 22 || parity_size > 46) {
        LDPC_decoder_free (ldpc);
        return NULL;
    }

    ldpc->msg_size = msg_size;
    ldpc->parity_size = parity_size;


    ldpc->Hm.rows = parity_size;
    ldpc->Hm.cols = msg_size;
    
    ldpc->Hm.base_graph = BG1_msg->base_graph;
    ldpc->Hm.rowOffset = BG1_msg->rowOffset;
    ldpc->Hm.rowWeight = BG1_msg->rowWeight;
    ldpc->Hm.columnIndexMap = BG1_msg->columnIndexMap;

    ldpc->Hp.rows = parity_size;
    ldpc->Hp.cols = parity_size;
    
    ldpc->Hp.base_graph = BG1_parity->base_graph;
    ldpc->Hp.rowOffset = BG1_parity->rowOffset;
    ldpc->Hp.rowWeight = BG1_parity->rowWeight;
    ldpc->Hp.columnIndexMap = BG1_parity->columnIndexMap;

    ldpc->cfg = cfg;

    return ldpc;
}

void LDPC_encoder_destroy (ldpc_encoder_t* ldpc) {
    if (ldpc == NULL) return;
    LDPC_encoder_free (ldpc);
}

void LDPC_decoder_destroy (ldpc_decoder_t* ldpc) {
    if (ldpc == NULL) return;
    LDPC_decoder_free (ldpc);
}