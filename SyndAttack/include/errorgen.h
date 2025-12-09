#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

void generate_random_error(sparse_t e_block, index_t weight,
                           index_t blocklength, prng_t prng);
void generate_around_word(sparse_t e_dst, index_t weight_dst, sparse_t e_src,
                          index_t weight_src, index_t intersections,
                          index_t blocklength, prng_t prng);
void generate_near_codeword(sparse_t e_block, code_t *H, int error_weight, int error_floor_p, prng_t prng);
void generate_near_codeword2(sparse_t e_block, code_t *H, int error_weight, int error_floor_p, prng_t prng);
void generate_codeword(sparse_t e_block, code_t *H, int error_weight, int error_floor_p, prng_t prng);

void errgen_init(int weight, int param);
int errgen_weight();
void errgen_sparse(index_t *error_sparse, code_t *H, long error_seed, long key_seed);
void errgen(e_t *e, code_t *H, long error_seed, long key_seed);
