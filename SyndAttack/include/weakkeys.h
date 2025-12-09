#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

void generate_weak_type1(code_t *H, index_t blocklength, int weak_p, prng_t prng);
void generate_weak_type2(code_t *H, index_t blocklength, int weak_p, prng_t prng);
void generate_weak_type3(code_t *H, index_t blocklength, int weak_p, prng_t prng);
