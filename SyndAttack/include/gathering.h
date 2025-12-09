#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

void gathering_error(sparse_t e_block, index_t blocklength, index_t error_weight, index_t m, index_t epsilon, prng_t prng);
void gathering_error_morphism(sparse_t e_block, index_t blocklength, index_t error_weight, index_t m, index_t epsilon, prng_t prng, index_t l);
void gathering_key(code_t *H, index_t blocklength, index_t m, index_t epsilon, prng_t prng);
