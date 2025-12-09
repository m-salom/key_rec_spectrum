#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

index_t insert_sorted(sparse_t array, index_t value, index_t max_i);
void insert_sorted_noinc(sparse_t array, index_t value, index_t max_i);

sparse_t sparse_new(index_t weight);
void sparse_free(sparse_t array);
void sparse_rand(sparse_t array, index_t weight, index_t length, prng_t prng);

void transpose(sparse_t dst, const sparse_t src, index_t weight,
               index_t length);
void transpose_columns(code_t *H);
void transpose_rows(code_t *H);

void multiply_xor_mod2(dense_t z, const sparse_t x, const dense_t y,
                       index_t block_weight, index_t block_length);
void multiply_add(dense_t z, const sparse_t x, const dense_t y,
                  index_t block_weight, index_t block_length);
#ifdef AVX
void multiply_xor_mod2_avx2(dense_t z, const sparse_t x, const dense_t y,
                            index_t block_weight, index_t block_length);
void multiply_avx2(dense_t z, const sparse_t x, const dense_t y,
                   index_t block_weight, index_t block_length);
#endif
