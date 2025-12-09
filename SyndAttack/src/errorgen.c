#include "errorgen.h"
#include "code.h"
#include "param.h"
#include "sparse_cyclic.h"

void generate_random_error(sparse_t e_block, index_t weight, index_t blocklength, prng_t prng) {
    sparse_rand(e_block, weight, INDEX * blocklength, prng);
}

/* Generate an error pattern with a fixed number intersections and specific
 * total weight. */
void generate_around_word(sparse_t e_dst, index_t weight_dst, sparse_t e_src,
                          index_t weight_src, index_t intersections,
                          index_t blocklength, prng_t prng) {
    bit_t e[INDEX * BLOCK_LENGTH] = {0};
    bit_t h[INDEX * BLOCK_LENGTH] = {0};

    /* Pick a near-codeword. */
    index_t shift = prng->random_lim(blocklength, prng->s);
    for (index_t l = 0; l < weight_src; ++l) {
        index_t k = e_src[l] >= blocklength;
        index_t i =
            k ? (e_src[l] - blocklength + shift) % blocklength + blocklength
              : (e_src[l] + shift) % blocklength;
        h[i] = 1;
    }

    /* Pick `intersections` common positions with the near-codeword. */
    index_t error_weight = 0;
    while (error_weight < intersections) {
        index_t lrand = prng->random_lim(weight_src, prng->s);
        index_t k = e_src[lrand] >= blocklength;
        index_t i = k ? (e_src[lrand] - blocklength + shift) % blocklength +
                            blocklength
                      : (e_src[lrand] + shift) % blocklength;
        if (!e[i]) {
            e[i] = 1;
            insert_sorted_noinc(e_dst, i, error_weight++);
        }
    }

    /* Complete the error pattern. */
    while (error_weight < weight_dst) {
        index_t jrand = prng->random_lim(INDEX * blocklength, prng->s);
        if (!e[jrand] && !h[jrand]) {
            e[jrand] = 1;
            insert_sorted_noinc(e_dst, jrand, error_weight++);
        }
    }
}

/* Generate an error pattern with error_floor_p intersection with a
 * (BLOCK_WEIGHT, BLOCK_WEIGHT) near-codeword. */
void generate_near_codeword(sparse_t e_block, code_t *H, int error_weight, int error_floor_p, prng_t prng) {
    index_t e_src[BLOCK_WEIGHT];

    /* Pick a near-codeword. */
    index_t k = prng->random_lim(INDEX, prng->s);

    for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
        e_src[l] = k * H->blocklength + H->columns[k][l];
    }

    generate_around_word(e_block, error_weight, e_src, BLOCK_WEIGHT,
                         error_floor_p, H->blocklength, prng);
}

/* Generate an error pattern with error_floor_p intersection with a
 * (2 * BLOCK_WEIGHT, ~2 * BLOCK_WEIGHT) near-codeword. */
void generate_near_codeword2(sparse_t e_block, code_t *H, int error_weight, int error_floor_p, prng_t prng) {
    index_t e_src[2 * BLOCK_WEIGHT];

    for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
        e_src[l] = H->columns[0][l];
    }

    /* Pick a near-codeword. */
    index_t shift = prng->random_lim(H->blocklength, prng->s);
    for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
        index_t i = H->columns[1][l] + shift;
        i = (i <= H->blocklength) ? i : (i - H->blocklength);
        e_src[BLOCK_WEIGHT + l] = H->blocklength + i;
    }

    generate_around_word(e_block, error_weight, e_src, 2 * BLOCK_WEIGHT,
                         error_floor_p, H->blocklength, prng);
}

/* Generate an error pattern with error_floor_p intersection with a
 * codeword of weight 2 * BLOCK_WEIGHT. */
void generate_codeword(sparse_t e_block, code_t *H, int error_weight, int error_floor_p, prng_t prng) {
    index_t e_src[INDEX * BLOCK_WEIGHT];

    for (index_t k = 0; k < INDEX; ++k) {
        index_t l;
        for (l = 0; l < BLOCK_WEIGHT; ++l) {
            index_t i = H->columns[INDEX - 1 - k][l];
            if (i >= H->blocklength)
                break;
            e_src[k * BLOCK_WEIGHT + l] = i + k * H->blocklength;
        }
        for (; l < BLOCK_WEIGHT; ++l) {
            index_t i = H->columns[INDEX - 1 - k][l];
            e_src[k * BLOCK_WEIGHT + l] = i + (k - 1) * H->blocklength;
        }
    }

    generate_around_word(e_block, error_weight, e_src, 2 * BLOCK_WEIGHT,
                         error_floor_p, H->blocklength, prng);
}

static int error_weight = ERROR_WEIGHT;
static int error_param = -1;

void errgen_init(int weight, int param) {
    error_weight = weight;
    error_param = param;
}

int errgen_weight() {
    return error_weight;
}

void errgen_sparse(index_t *error_sparse, code_t *H, long error_seed, long key_seed) {
    struct PRNG prng;

    prng_from_seed2(&prng, error_seed, key_seed);
    if (error_param >= 0)
        generate_near_codeword(error_sparse, H, error_weight, error_param, &prng);
    else
        generate_random_error(error_sparse, error_weight, H->blocklength, &prng);
}

void errgen(e_t *e, code_t *H, long error_seed, long key_seed) {
    index_t error_sparse[error_weight];
    errgen_sparse(error_sparse, H, error_seed, key_seed);
    error_sparse_to_dense(e, error_sparse, error_weight, H->blocklength);
}
