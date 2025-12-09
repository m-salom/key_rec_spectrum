#include "codegen.h"
#include "gathering.h"
#include "param.h"
#include "sparse_cyclic.h"
#include "weakkeys.h"

static int blocklength = BLOCK_LENGTH;
static int weak_key_type = NONE;
static int weak_param_1 = 0, weak_param_2 = 0;

void generate_random_code(code_t *H, prng_t prng) {
    for (index_t i = 0; i < INDEX; ++i) {
        sparse_rand(H->columns[i], BLOCK_WEIGHT, blocklength, prng);
    }
    H->blocklength = blocklength;
    transpose_columns(H);
}

void keygen_init(int length, int type, int p1, int p2) {
    blocklength = length;
    weak_key_type = type;
    weak_param_1 = p1;
    weak_param_2 = p2;
}

int keygen_blocklength() {
    return blocklength;
}

void keygen(code_t *H, long seed) {
    struct PRNG prng;

    prng_from_seed(&prng, seed);

    if (weak_key_type == M_GATHERING) {
        gathering_key(H, blocklength, weak_param_1, weak_param_2, &prng);
    } else if (weak_key_type == TYPE_I) {
        generate_weak_type1(H, blocklength, weak_param_1, &prng);
    } else if (weak_key_type == TYPE_II) {
        generate_weak_type2(H, blocklength, weak_param_1, &prng);
    } else if (weak_key_type == TYPE_III) {
        generate_weak_type3(H, blocklength, weak_param_1, &prng);
    } else {
        generate_random_code(H, &prng);
    }
}


