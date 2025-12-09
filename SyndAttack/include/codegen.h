#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

#define NONE 0
#define TYPE_I 1
#define TYPE_II 2
#define TYPE_III 3
#define M_GATHERING 4

void generate_random_code(code_t *H, prng_t prng);
void keygen_init(int length, int type, int p1, int p2);
int keygen_blocklength();
void keygen(code_t *H, long seed);
