#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

void transpose_columns(code_t *H);
void transpose_rows(code_t *H);

void compute_codeword(cw_t codeword, code_t *H, msg_t message);
void compute_syndrome(syndrome_t *syndrome, code_t *H, e_t *e_dense);
void compute_counters(counters_t counters, bit_t *syndrome, code_t *H);

void error_sparse_to_dense(e_t *e_dense, const sparse_t e_sparse,
                           index_t weight, index_t blocklength);
void syndrome_add_sparse_error(syndrome_t *syndrome, const sparse_t e_sparse,
                               index_t weight);
