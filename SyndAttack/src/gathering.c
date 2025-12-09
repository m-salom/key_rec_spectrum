#include <stdlib.h>
#include <string.h>

#include "param.h"
#include "sparse_cyclic.h"

void gathering_error(sparse_t e_block, index_t blocklength, index_t error_weight, index_t m, index_t epsilon, prng_t prng) {
    index_t i, t0 = error_weight / 2;
    if (t0 < epsilon) // just a safeguard, should not happen
        epsilon = t0;
    if (m > blocklength)
        m = blocklength;
    sparse_rand(e_block, t0 - epsilon, m, prng);
    sparse_rand(e_block + t0 - epsilon, epsilon, blocklength - m, prng);
    for (i = t0 - epsilon; i < t0; ++i)
        e_block[i] += m;
    sparse_rand(e_block + t0, error_weight - t0, blocklength, prng);
    for (i = t0; i < error_weight; ++i)
        e_block[i] += blocklength;
}

int cmp_int(const void *x, const void *y) {
    return *((index_t *)x) - *((index_t *)y);
}

void gathering_error_morphism(sparse_t e_block, index_t blocklength, index_t error_weight, index_t m, index_t epsilon, prng_t prng, index_t l) {
    index_t i, t0 = error_weight / 2;
    if (t0 < epsilon) // just a safeguard, should not happen
        epsilon = t0;
    if (m > blocklength)
        m = blocklength;
    sparse_rand(e_block, t0 - epsilon, m, prng);
    sparse_rand(e_block + t0 - epsilon, epsilon, blocklength - m, prng);
    for (i = t0 - epsilon; i < t0; ++i)
        e_block[i] += m;
    /* apply morphism */
    for (i = 0; i < t0; ++i)
        e_block[i] = (e_block[i] * l) % blocklength;
    qsort(e_block, t0, sizeof(index_t), cmp_int);
    sparse_rand(e_block + t0, error_weight - t0, blocklength, prng);
    for (i = t0; i < error_weight; ++i)
        e_block[i] += blocklength;
}

void gathering_key(code_t *H, index_t blocklength, index_t m, index_t epsilon, prng_t prng) {
    index_t i;
    sparse_rand(H->columns[0], BLOCK_WEIGHT - epsilon, m, prng);
    sparse_rand(H->columns[0] + BLOCK_WEIGHT - epsilon, epsilon, blocklength - m, prng);
    for (i = BLOCK_WEIGHT - epsilon; i < BLOCK_WEIGHT; ++i)
        H->columns[0][i] += m;
    for (index_t i = 1; i < INDEX; ++i) {
        sparse_rand(H->columns[i], BLOCK_WEIGHT, blocklength, prng);
    }
    H->blocklength = blocklength;
    transpose_columns(H);
}
