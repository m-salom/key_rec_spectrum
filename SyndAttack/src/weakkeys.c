#include "weakkeys.h"
#include "param.h"
#include "sparse_cyclic.h"

void generate_weak_type1(code_t *H, index_t blocklength, int weak_p, prng_t prng) {
    index_t k = prng->random_lim(INDEX, prng->s);
    index_t delta = 1 + prng->random_lim(blocklength / 2, prng->s);
    index_t shift = prng->random_lim(blocklength, prng->s);

    index_t length_left = blocklength;
    for (index_t i = 0; i < weak_p; i++) {
        uint32_t a = (delta * (i + shift)) % blocklength;
        insert_sorted_noinc(H->columns[k], a, i);
    }
    length_left -= weak_p;
    for (index_t i = weak_p; i < BLOCK_WEIGHT; i++) {
        uint32_t rand = prng->random_lim(length_left--, prng->s);
        insert_sorted(H->columns[k], rand, i);
    }

    sparse_rand(H->columns[INDEX - 1 - k], BLOCK_WEIGHT, blocklength, prng);
    H->blocklength = blocklength;
    transpose_columns(H);
}

/* Generate a polynomial with a multiplicity of weak_p using the stars and bars
 * principle. */
void generate_weak_type2(code_t *H, index_t blocklength, int weak_p, prng_t prng) {
    index_t k = prng->random_lim(INDEX, prng->s);
    index_t delta = 1 + prng->random_lim(blocklength / 2, prng->s);

    const index_t s = BLOCK_WEIGHT - weak_p;
    /* First the ois and zis represent the "bars". */
    index_t ois[s + 1];
    index_t zis[s + 1];
    ois[0] = 0;
    ois[s] = BLOCK_WEIGHT;
    index_t left = BLOCK_WEIGHT - 1;
    for (index_t i = 1; i < s; i++) {
        uint32_t rand = prng->random_lim(left--, prng->s);
        insert_sorted(ois, rand, i);
    }
    zis[0] = 0;
    zis[s] = blocklength - BLOCK_WEIGHT;
    left = blocklength - BLOCK_WEIGHT - 1;
    for (index_t i = 1; i < s; i++) {
        uint32_t rand = prng->random_lim(left--, prng->s);
        insert_sorted(zis, rand, i);
    }
    /* Convert ois and zis to run-length. */
    for (index_t i = 0; i < s; ++i) {
        ois[i] = ois[i + 1] - ois[i];
        zis[i] = zis[i + 1] - zis[i];
    }

    index_t shift = prng->random_lim(ois[0] + zis[0], prng->s);
    index_t current_pos = (blocklength - shift) % blocklength;
    index_t i = 0;
    for (index_t l1 = 0; l1 < s; ++l1) {
        current_pos = (current_pos + zis[l1]) % blocklength;
        for (index_t l2 = 0; l2 < ois[l1]; ++l2) {
            uint32_t a = (delta * (current_pos + l2)) % blocklength;
            insert_sorted_noinc(H->columns[k], a, i++);
        }
        current_pos = (current_pos + ois[l1]) % blocklength;
    }

    sparse_rand(H->columns[INDEX - 1 - k], BLOCK_WEIGHT, blocklength, prng);
    H->blocklength = blocklength;
    transpose_columns(H);
}

void generate_weak_type3(code_t *H, index_t blocklength, int weak_p, prng_t prng) {
    index_t length = blocklength;
    index_t shift = prng->random_lim(blocklength, prng->s);

    /* Choose weak_p common values. */
    for (index_t i = 0; i < weak_p; i++) {
        index_t rand = prng->random_lim(length--, prng->s);
        rand = insert_sorted(H->columns[0], rand, i);
        uint32_t a = (rand + shift) % blocklength;
        insert_sorted_noinc(H->columns[1], a, i);
    }

    /* Complete H->columns[0]. */
    for (index_t i = weak_p; i < BLOCK_WEIGHT; i++) {
        index_t rand = prng->random_lim(length--, prng->s);
        insert_sorted(H->columns[0], rand, i);
    }
    length += BLOCK_WEIGHT - weak_p;

    /* Complete H->columns[1] without increasing the number of intersections
     * with H->columns[0].
     */
    for (index_t i = weak_p; i < BLOCK_WEIGHT; i++) {
    gen : {
        index_t rand = prng->random_lim(length, prng->s);

        for (index_t j = 0; j < i && H->columns[1][j] <= rand; j++, rand++)
            ;
        for (index_t k = 0; k < BLOCK_WEIGHT; k++) {
            uint32_t a = (H->columns[0][k] + shift) % blocklength;
            if (a == rand)
                goto gen;
        }
        insert_sorted_noinc(H->columns[1], rand, i);
        --length;
    }
    }
    H->blocklength = blocklength;
    transpose_columns(H);
}
