#include <stdlib.h>
#include <string.h>

#include "poly.h"

poly_t poly_alloc(int length, int weight) {
    poly_t h = malloc(sizeof (struct poly));
    h->length = length;
    h->weight = weight;
    h->index = (uint16_t *) calloc(weight, sizeof (uint16_t));
    h->alloc = 0;
    h->extended = NULL;
    return h;
}

int rand_int(int n) {
    return random() % n;
}

void rand_sparse(uint16_t * index, int length, int weight) {
    for (int i = weight - 1; i >= 0; --i) {
        int j = i + rand_int(length - i);
        for (int l = weight - 1; l > i; --l) {
            if (index[l] == j) {
                j = i;
                break;
            }
        }
        index[i] = j;
    }
}

void poly_rand(poly_t h) {
    rand_sparse(h->index, h->length, h->weight);
}

void poly_free(poly_t h) {
    free(h->index);
    if (h->extended != NULL)
        free(h->extended);
    free(h);
}

void extended_alloc(poly_t f) {
    f->alloc = 1 + (f->length - 1) / (8 * sizeof (uint64_t));
    f->extended = (uint64_t *) calloc(f->alloc, sizeof (uint64_t));
}

void extended_set_all_to_one(poly_t f) { // all coeffs in extension get 1
    memset(f->extended, -1, f->alloc * sizeof (uint64_t));
    int l = f->length % (8 * sizeof (uint64_t));
    if (l > 0) { // pad the end of the final block with zeroes
        f->extended[f->alloc - 1] &= (1UL << l) - 1;
    }
}

void extended_set_to_zero(poly_t f, int i) {
    int wordsize = 8 * sizeof (uint64_t);
    int pos = i / wordsize;
    uint64_t bit = (1UL << (i % wordsize));
    f->extended[pos] &= ~bit;
}

void extended_set_to_one(poly_t f, int i) {
    int wordsize = 8 * sizeof (uint64_t);
    int pos = i / wordsize;
    uint64_t bit = (1UL << (i % wordsize));
    f->extended[pos] |= bit;
}

int extended_weight(poly_t h) {
    int i, w;
    if (h->extended == NULL)
        return -1;
    for (i = 0, w = 0; i < h->alloc; ++i)
        w += __builtin_popcountl(h->extended[i]);
    return w;
}

bool extended_coeff_is_zero(poly_t h, int i) {
    int wordsize = 8 * sizeof (uint64_t);
    int pos = i / wordsize;
    int j = i % wordsize;
    return (h->extended[pos] & (1UL << j)) == 0;
}

/*
   ZER0 in positions that cannot be ONE with the given spectrum multiplicity S
   assumes length is odd
   assumes 0 < 2 * j < length
*/
poly_t extend_binomial_init(int j, int8_t * S, int length) {
    poly_t f = poly_alloc(length, 2);
    f->index[0] = 0;
    f->index[1] = j;
    extended_alloc(f);
    extended_set_all_to_one(f);

#ifdef FULL_MULTIPLICITY
    if (S[j] == 1) {
        extended_set_to_zero(f, length - j);
        extended_set_to_zero(f, 2 * j);
    }
#endif

    // general case
    for (int i = 1; i <= f->length / 2; ++i) {
        if (S[i] == 0) {
            extended_set_to_zero(f, i);
            extended_set_to_zero(f, f->length - i);
            extended_set_to_zero(f, j + i);
            if (j - i >= 0) {
                extended_set_to_zero(f, j - i);
            }
            else {
                extended_set_to_zero(f, f->length + j - i);
            }
        }
    }

    extended_set_to_one(f, 0);
    extended_set_to_one(f, j);

    return f;
}

poly_t poly_shift_union(poly_t g, poly_t f, int l) {
    int i, j, k, wordsize, r = g->length;
    poly_t h = poly_alloc(r, g->weight + f->weight);

    /* sparse polynomial */
    memcpy(h->index, g->index, g->weight * sizeof (g->index[0]));
    /* append coeffs of f, avoiding duplicates */
    for (i = 0, k = g->weight; i < f->weight; ++i) {
        h->index[k] = (f->index[i] + l) % r;
        for (j = 0; j < g->weight; ++j)
            if (h->index[j] == h->index[k])
                break;
        if (j == g->weight)
            ++k;
    }
    h->weight = k;
    

    /* extended dense polynomial */
    if ((f->extended == 0) || (g->extended == 0))
        return h;
    extended_alloc(h);
    /* first set h_ext <- x^l * f_ext */
    /* first l coordinates (copy end of f at the beginning of h) */
    /* assumes that the end of f is correctly padded with zeroes */
    wordsize = 8 * sizeof (uint64_t);
    j = (r - l) / wordsize;
    k = (r - l) % wordsize;
    if (k == 0) {
        memcpy(h->extended, f->extended + j, (f->alloc - j) * sizeof (uint64_t));
    } else {
        for (i = 0; i + j + 1 < f->alloc; ++i) {
            h->extended[i] = (f->extended[i + j] >> k) | (f->extended[i + j + 1] << (wordsize - k));
        }
        h->extended[i] = f->extended[i + j] >> k;
    }
    /* last r-l coordinates (copy beginning of f at the end of h) */
    j = l / wordsize;
    k = l % wordsize;
    if (k == 0) {
        memcpy(h->extended + j, f->extended, (f->alloc - j) * sizeof (uint64_t));
    } else {
        h->extended[j] |= f->extended[0] << k;
        for (i = 0; i + j + 1 < h->alloc; ++i) {
            h->extended[i + j + 1] = (f->extended[i] >> (wordsize - k)) | (f->extended[i + 1] << k);
        }
        /* padding the end of h last block with zeroes is done next
           (assuming g is correctly padded with 0) */
    }

    /* second h_ext <- h_ext & g_ext */
    for (i = 0; i < g->alloc; ++i) {
        h->extended[i] &= g->extended[i];
    }
    return h;
}
