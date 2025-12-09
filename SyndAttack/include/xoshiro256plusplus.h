#pragma once
#include <stdint.h>

uint64_t random_uint64_t(uint64_t *s);
uint64_t random_lim(uint64_t limit, uint64_t *s);

struct PRNG {
    uint64_t s[4];
    uint64_t (*random_lim)(uint64_t limit, uint64_t *s);
    uint64_t (*random_uint64_t)(uint64_t *s);
};

typedef struct PRNG *prng_t;

void prng_from_seed(prng_t prng, uint64_t seed);
void prng_from_seed2(prng_t prng, uint64_t seed1, uint64_t seed2);
void prng_init(prng_t prng);
