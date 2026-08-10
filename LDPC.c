#include "LDPC_impl.h"
#include "matrix.h"
#include "base_graph.h"
#include "LDPC_memory.h"
#include <stddef.h>

typedef struct {
    uint16_t msg_len;
    uint16_t target_code_len;

    uint8_t bgn;

    lifting_mode_t lifting_mode;
    uint16_t lifting_size;
    uint16_t max_lifting_size;

} ldpc_common_cfg_t;

static ldpc_config_status_t LDPC_validate_common_config (ldpc_common_cfg_t common) {

    if (common.msg_len > common.target_code_len) {
        return LDPC_CONFIG_FAILURE;
    }
    
    if (common.bgn == 1) {
        if (common.msg_len > (BG1_COLS - BG1_ROWS)*common.max_lifting_size 
            || common.target_code_len > BG1_COLS*common.max_lifting_size) {
            return LDPC_CONFIG_FAILURE;
        }
    }
    else if (common.bgn == 2) {
        if (common.msg_len > (BG2_COLS - BG2_ROWS)*common.max_lifting_size 
            || common.target_code_len > BG2_COLS*common.max_lifting_size) {
            return LDPC_CONFIG_FAILURE;
        }
    }
    else {
        return LDPC_CONFIG_FAILURE;
    }

    if (common.lifting_mode == LDPC_LIFTING_AUTO) {
        common.lifting_size = get_nearest_lifting_size (common.msg_len, common.target_code_len, common.bgn);
        if (common.lifting_size == 0 
            || common.lifting_size > common.max_lifting_size) {
            return LDPC_CONFIG_FAILURE;
        }
    }
    else if (common.lifting_mode == LDPC_LIFTING_EXPLICIT) {
        if (!is_valid_lifting_size (common.lifting_size) 
            || common.lifting_size > common.max_lifting_size) {
            return LDPC_CONFIG_FAILURE;
        }
    }
    else {
        return LDPC_CONFIG_FAILURE;
    }

    uint16_t msg_size = DIV_CEIL(common.msg_len, common.lifting_size);
    uint16_t parity_size = DIV_CEIL(common.target_code_len - common.msg_len, common.lifting_size);

    if (parity_size <= 4) {
        return LDPC_CONFIG_FAILURE;
    }
    if (common.bgn == 1 && (msg_size > (BG1_COLS - BG1_ROWS) 
        || parity_size > BG1_ROWS)) {
        return LDPC_CONFIG_FAILURE;
    }
    if (common.bgn == 2 && (msg_size > (BG2_COLS - BG2_ROWS) 
        || parity_size > BG2_ROWS)) {
        return LDPC_CONFIG_FAILURE;
    }

    return LDPC_CONFIG_SUCCESS;
}

ldpc_config_status_t LDPC_set_encoder_config (ldpc_encoder_t* ldpc, ldpc_encoder_cfg_t cfg) {

    if (ldpc == NULL) return LDPC_CONFIG_FAILURE;

    ldpc_common_cfg_t common = {
        .msg_len = cfg.msg_len,
        .target_code_len = cfg.target_code_len,

        .bgn = cfg.bgn,

        .lifting_mode = cfg.lifting_mode,
        .lifting_size = cfg.lifting_size,
        .max_lifting_size = MAX_LIFTING_SIZE_ENCODE
    };
    
    ldpc_config_status_t status = LDPC_validate_common_config (common);
    if (status == LDPC_CONFIG_FAILURE) return LDPC_CONFIG_FAILURE;
    
    if (common.lifting_mode == LDPC_LIFTING_AUTO) {
        cfg.lifting_size = get_nearest_lifting_size (cfg.msg_len, cfg.target_code_len, cfg.bgn);
    }
    ldpc->cfg = cfg;

    int8_t set_index = get_set_index (cfg.lifting_size);

    if (cfg.bgn == 1) {
        if (set_index == 6) {
            ldpc->B = BG1_B2;
        }
        else {
            ldpc->B = BG1_B1;
        }
    }
    else {
        if (set_index == 3 || set_index == 7) {
            ldpc->B = BG2_B2;
        }
        else {
            ldpc->B = BG2_B1;
        }
    }

    const quasi_cyclic_matrix_t * BG_msg = (cfg.bgn == 1) ? get_BG1_msg (set_index) : get_BG2_msg (set_index);
    const quasi_cyclic_matrix_t * BG_parity = (cfg.bgn == 1) ? get_BG1_parity (set_index) : get_BG2_parity (set_index);

    uint16_t msg_size = DIV_CEIL(cfg.msg_len, cfg.lifting_size);
    uint16_t parity_size = DIV_CEIL(cfg.target_code_len - cfg.msg_len, cfg.lifting_size);

    ldpc->msg_size = msg_size;
    ldpc->parity_size = parity_size;

    ldpc->A.rows = 4;
    ldpc->A.cols = msg_size;
    
    ldpc->A.base_graph = BG_msg->base_graph;
    ldpc->A.rowOffset = BG_msg->rowOffset;
    ldpc->A.rowWeight = BG_msg->rowWeight;
    ldpc->A.columnIndexMap = BG_msg->columnIndexMap;

    ldpc->C.rows = parity_size - 4;
    ldpc->C.cols = msg_size;
    
    ldpc->C.base_graph = BG_msg->base_graph;
    ldpc->C.rowOffset = BG_msg->rowOffset + 4;
    ldpc->C.rowWeight = BG_msg->rowWeight + 4;
    ldpc->C.columnIndexMap = BG_msg->columnIndexMap;

    ldpc->D.rows = parity_size - 4;
    ldpc->D.cols = 4;
    
    ldpc->D.base_graph = BG_parity->base_graph;
    ldpc->D.rowOffset = BG_parity->rowOffset + 4;
    ldpc->D.rowWeight = BG_parity->rowWeight + 4;
    ldpc->D.columnIndexMap = BG_parity->columnIndexMap;

    return LDPC_CONFIG_SUCCESS;

}

ldpc_config_status_t LDPC_set_decoder_config (ldpc_decoder_t* ldpc, ldpc_decoder_cfg_t cfg) {
    
    if (ldpc == NULL) return LDPC_CONFIG_FAILURE;

    ldpc_common_cfg_t common = {
        .msg_len = cfg.msg_len,
        .target_code_len = cfg.target_code_len,

        .bgn = cfg.bgn,

        .lifting_mode = cfg.lifting_mode,
        .lifting_size = cfg.lifting_size,
        .max_lifting_size = MAX_LIFTING_SIZE_DECODE
    };
    
    ldpc_config_status_t status = LDPC_validate_common_config (common);
    if (status == LDPC_CONFIG_FAILURE) return LDPC_CONFIG_FAILURE;

    #if !USE_SUM_PRODUCT
        if (cfg.alpha > 1.0f || cfg.alpha <= 0.0f) {
            return LDPC_CONFIG_FAILURE;
        }
    #endif

    if (common.lifting_mode == LDPC_LIFTING_AUTO) {
        cfg.lifting_size = get_nearest_lifting_size (cfg.msg_len, cfg.target_code_len, cfg.bgn);
    }
    ldpc->cfg = cfg;

    int8_t set_index = get_set_index (cfg.lifting_size);

    const quasi_cyclic_matrix_t * BG_msg = (cfg.bgn == 1) ? get_BG1_msg (set_index) : get_BG2_msg (set_index);
    const quasi_cyclic_matrix_t * BG_parity = (cfg.bgn == 1) ? get_BG1_parity (set_index) : get_BG2_parity (set_index);

    uint16_t msg_size = DIV_CEIL(cfg.msg_len, cfg.lifting_size);
    uint16_t parity_size;
    if (cfg.use_full_PCM) {
        if (cfg.bgn == 1) {
            parity_size = BG1_ROWS;
        }
        else {
            parity_size = BG2_ROWS;
        }
    }
    else {
        parity_size = DIV_CEIL(cfg.target_code_len - cfg.msg_len, cfg.lifting_size);
    }

    ldpc->msg_size = msg_size;
    ldpc->parity_size = parity_size;

    ldpc->Hm.rows = parity_size;
    ldpc->Hm.cols = msg_size;
    
    ldpc->Hm.base_graph = BG_msg->base_graph;
    ldpc->Hm.rowOffset = BG_msg->rowOffset;
    ldpc->Hm.rowWeight = BG_msg->rowWeight;
    ldpc->Hm.columnIndexMap = BG_msg->columnIndexMap;

    ldpc->Hp.rows = parity_size;
    ldpc->Hp.cols = parity_size;
    
    ldpc->Hp.base_graph = BG_parity->base_graph;
    ldpc->Hp.rowOffset = BG_parity->rowOffset;
    ldpc->Hp.rowWeight = BG_parity->rowWeight;
    ldpc->Hp.columnIndexMap = BG_parity->columnIndexMap;

    return LDPC_CONFIG_SUCCESS;

}


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
    
    ldpc_encoder_t* ldpc = LDPC_encoder_alloc ();
    if (ldpc == NULL) return NULL;

    ldpc_config_status_t status = LDPC_set_encoder_config (ldpc, cfg);

    if (status == LDPC_CONFIG_FAILURE) {
        LDPC_encoder_free (ldpc);
        return NULL;
    }
    
    return ldpc;
}

ldpc_decoder_t* LDPC_decoder_create (ldpc_decoder_cfg_t cfg) {
    
    ldpc_decoder_t* ldpc = LDPC_decoder_alloc ();
    if (ldpc == NULL) return NULL;  

    ldpc_config_status_t status = LDPC_set_decoder_config (ldpc, cfg);

    if (status == LDPC_CONFIG_FAILURE) {
        LDPC_decoder_free (ldpc);
        return NULL;
    }

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