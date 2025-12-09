#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "poly.h"
#include "spectrum.h"

/* put in the array res[] all the offsets s such that
   Spectrum(g union x^s+x^(s+delta)) is included in S
   returns the number of such s */
int valid_offset(uint16_t * res, poly_t g, uint16_t delta, int8_t * S) {
    int i, j, l, k = 0, r = g->length;
    for (l = 0; l < r; ++l) {
        j = (l < r - l) ? l : (r - l);
        if ((l == 0) || (S[j] > 0)) {
            j = (l + delta) % r;
            for (i = 0; i < g->weight; ++i) {
                if (check_distance(l, g->index[i], S, r) &&
                    check_distance(j, g->index[i], S, r))
                    continue;
                break;
            }
            if (i == g->weight) {
                res[k++] = l;
            }
        }
    }
    return k;
}

/* extract the extension of h into g
   returns g if its spectrum is exactly S
   else returns NULL */
poly_t check_key(poly_t h, int8_t * S) {
    int i, j, w, r = h->length;
    uint16_t index[r]; // oversized
    int8_t mult[(r + 1) / 2];
    memset(mult, 0, ((r + 1) / 2) * sizeof (int8_t));

    for (i = 0, w = 0; i < r; ++i) 
        if (!extended_coeff_is_zero(h, i))
            index[w++] = i;

    for (i = 0; i < w; ++i)
        for (j = 0; j < i; ++j)
#ifndef FULL_MULTIPLICITY
            if (mult[distance(index[i], index[j], r)] == 0)
#endif
                mult[distance(index[i], index[j], r)] += 1;

    if (memcmp(S, mult, ((r + 1) / 2) * sizeof (int8_t)) != 0) {
        return NULL;
    }

    poly_t g = poly_alloc(r, w);
    memcpy(g->index, index, w * sizeof (uint16_t));

    return g;
}

// extract the extension of h 
poly_t extract_ext(poly_t h) {
    int i, w, r = h->length;
    uint16_t index[r];
  
    for (i = 0, w = 0; i < r; ++i) 
        if (!extended_coeff_is_zero(h, i))
            index[w++] = i;
  
    poly_t g = poly_alloc(r, w);
    memcpy(g->index, index, w * sizeof (uint16_t));

    return g;
}


int count_check_mult = 0, fail_check_mult = 0, count_guess_key = 0, count_nb_ext = 0;

/* return true if the spectrum of the extension of h contains the
   spectrum S */
bool check_mult(poly_t h, int8_t * S) {
    int i, j, w, r = h->length;
    uint16_t index[r];
    int16_t mult[(r + 1) / 2];

    count_check_mult++;

    for (i = 0, w = 0; i < r; ++i)
        if (!extended_coeff_is_zero(h, i))
            index[w++] = i;
    
    memset(mult, 0, ((r + 1) / 2) * sizeof (mult[0]));
    
    for (i = 0; i < w; ++i)
        for (j = 0; j < i; ++j) 
            mult[distance(index[i], index[j], r)] += 1;

    for (i = 0; i < (r + 1) / 2; ++i) {
#ifdef FULL_MULTIPLICITY
        if (mult[i] < S[i])
#else
        if ((S[i] > 0) && (mult[i] == 0))
#endif
            {
                fail_check_mult++;
                return false;
            }
    }
    return true;
}

/* return true if the spectrum of the extension of h is contained
   in the spectrum S 
*/
bool check_spectrum(poly_t h, int8_t * S) {
    int i, j, w, r = h->length;
    uint16_t index[r];
    int8_t mult[(r + 1) / 2];

    count_check_mult++;

    for (i = 0, w = 0; i < r; ++i)
        if (!extended_coeff_is_zero(h, i))
            index[w++] = i;

    memset(mult, 0, ((r + 1) / 2) * sizeof (mult[0]));
    for (i = 0; i < w; ++i)
        for (j = 0; j < i; ++j)
            mult[distance(index[i], index[j], r)] += 1;

    for (i = 0; i < (r + 1) / 2; ++i) {
        if ((S[i] == 0) && (mult[i] > 0)) {
                fail_check_mult++;
                return false;
            }
    }

    return true;
}

poly_t guess_key_rec(uint16_t * delta, poly_t * f, int depth, int max_depth, poly_t g, int8_t * S, int weight) {
    int i, nb_offset, r = g->length;
    poly_t h, key = NULL;
    uint16_t V[r];
    count_guess_key++;
    
    /*
      the dense vector g depends on (delta[i], offset[i]), 0 <= i < depth
    */

    if (depth >= max_depth)
        return NULL;

    if (f[depth] == NULL) { // initialize, only once per depth
        f[depth] = extend_binomial_init(delta[depth], S, r);
    }
    
    nb_offset = valid_offset(V, g, delta[depth], S);

    /*
      the first nb_offset entries of V are valid offsets for the distance delta[depth]
      compatible simultaneously with all delta[i] shifted by offset[i], 0 <= i < depth
    */
    
    for (i = 0; i < nb_offset; ++i) {
        count_nb_ext++;
        h = poly_shift_union(g, f[depth], V[i]);

        // weight of the extension of h
        int w = extended_weight(h);

        // cannot lead to a word of target weight -> quit and try next entry of V
        if (w < weight) {
            poly_free(h);
            continue;
        }

        // no word deriving from h can reach the target spectrum multiplicity
        if (check_mult(h, S) == false) {
            poly_free(h);
            continue;
        }

        // if w equals the target weight we check the only possible solution
        if (w == weight) {
            key = check_key(h, S);
           
            poly_free(h);
            if (key != NULL) { // success
                key->alloc = depth; // hack to store success depth
                break;
            }
            else {
                continue;
            }
        }
        
        key = guess_key_rec(delta, f, depth + 1, max_depth, h, S, weight);
        poly_free(h);
        if (key != NULL) {
            break;
        }
    }
    return key;
}

// reconstruct key from partial spectrum : high bound and low bound
poly_t guess_key_rec_bounds(uint16_t * delta, poly_t * f, int depth, int max_depth, poly_t g, 
                            int8_t * S, int8_t * S_low, int weight) {
    int i, nb_offset, r = g->length;
    poly_t h, key = NULL;
    uint16_t V[r];

    count_guess_key++;
    
    /*
      the dense vector g depends on (delta[i], offset[i]), 0 <= i < depth
    */

    if (depth >= max_depth)
        return NULL;

    if (f[depth] == NULL) { // initialize, only once per depth
        f[depth] = extend_binomial_init(delta[depth], S, r);
    }
    
    nb_offset = valid_offset(V, g, delta[depth], S);

    /*
      the first nb_offset entries of V are valid offsets for the distance delta[depth]
      compatible simultaneously with all delta[i] shifted by offset[i], 0 <= i < depth
    */
    for (i = 0; i < nb_offset; ++i) {
        count_nb_ext++;
        // current polynomial h is built from g and f[depth] shifted by V[i]
        h = poly_shift_union(g, f[depth], V[i]);

        // weight of the extension of h
        int w = extended_weight(h);
    
        // cannot lead to a word of target weight -> quit and try next entry of V
        if (w < weight) {
            poly_free(h);
            continue;
        }

        // no word deriving from h can reach the target spectrum multiplicity
        if (check_mult(h, S_low) == false) {
        
            poly_free(h);
            continue;
        }

        // if w equals the target weight we check the only possible solution
        if (w == weight) {  
            key = extract_ext(h);
            poly_free(h);
            if (key != NULL) { // success
                key->alloc = depth; // hack to store success depth
                break;
            }
            else {
                continue;
            }
        }
        
        // no solution, no contradiction, we continue from h and the next delta
        key = guess_key_rec_bounds(delta, f, depth + 1, max_depth, h, S, S_low, weight);
        poly_free(h);
        if (key != NULL) {
            break;
        }
    }
    return key;
}

poly_t guess_key(uint16_t * delta, int8_t * S, int r, int max_depth, int weight) {
    int i;
    poly_t f[max_depth];
    poly_t h;

    for (i = 0; i < max_depth; ++i) {
        f[i] = NULL;
    }
    f[0] = extend_binomial_init(delta[0], S, r);
    
    h = guess_key_rec(delta, f, 1, max_depth, f[0], S, weight);
    for (i = 0; i < max_depth; ++i) {
        if (f[i] != NULL)
            poly_free(f[i]);
    }
    return h;
}

poly_t guess_key_bounds(uint16_t * delta, int8_t * S, int8_t * S_low, int r, int max_depth, int weight) {
    int i;
    poly_t f[max_depth];
    poly_t h;

    for (i = 0; i < max_depth; ++i) {
        f[i] = NULL;
    }
    f[0] = extend_binomial_init(delta[0], S, r);

    h = guess_key_rec_bounds(delta, f, 1, max_depth, f[0], S, S_low, weight);
    for (i = 0; i < max_depth; ++i) {
        if (f[i] != NULL)
            poly_free(f[i]);
    }
    return h;
}

void sample_delta(uint16_t * delta, int8_t * S, int size, int nb_pos) {
    int i, j, l;
    for (i = 0; i < nb_pos; ++i) {
        do {
            do {
                j = 1 + random() % (size - 1);
            } while (S[j] == 0);
            for (l = 0; l < i; ++l)
                if (j == delta[l])
                    break;
        } while (l < i);
        delta[i] = j;
    }
}

static int int16comp (const void * a, const void * b) {
    return *((uint16_t *) a) - *((uint16_t *) b);
}

poly_t normalize(poly_t h) {
    int i, j, l, offset = 0, r = h->length, d = h->weight;
    poly_t h_normalized = poly_alloc(r, d);
    int8_t S[(r + 1) / 2];

    spectrum(S, h->index, r, d);
    /* l = smallest distance with multiplicity 1 */
    /* conjecture: at least one distance has multiplicity 1 */
    for (l = 0; S[l] != 1; ++l);

    for (i = 0; i < d; ++i) {
        for (j = 0; j < i; ++j) {
            if (distance(h->index[i], h->index[j], r) == l)
                break;
        }
        if (i != j) {
            if (h->index[i] == (h->index[j] + l) % r)
                offset = h->index[j];
            else
                offset = h->index[i];
            break;
        }
    }

    for (i = 0; i < d; ++i)
        h_normalized->index[i] = (r - offset + h->index[i]) % r;

    qsort(h_normalized->index, d, sizeof (uint16_t), int16comp);

    return h_normalized;
}

// Return the mirror of a polynomial
poly_t mirror_poly(poly_t p) {
    int d = p->weight;
    int r = p->length;

    poly_t m = poly_alloc(r,d);
    if (!m) return NULL;

    for (int i = 0; i < d; i++) {
        m->index[i] = (r - p->index[i]) % r;
    }

    return m;
}

bool is_key(poly_t guess_key, poly_t key_normalized) {
    if (guess_key == NULL) return false;
    poly_t guess_key_normalized = normalize(guess_key);
    bool success = true;
           
    for (int i = 0; i < key_normalized->weight; i++) {
        if (guess_key_normalized->index[i] != key_normalized->index[i]) {
            success = false;
            break;
        }
    }
                
    if (!success) {
        success = true;
        poly_t key_mirror = mirror_poly(guess_key);
        poly_t key_mirror_normalized = normalize(key_mirror);
        for (int i = 0; i < key_normalized->weight; i++) {
            if (key_mirror_normalized->index[i] != key_normalized->index[i]) {
                success = false;
                break;
            }
        }
        poly_free(key_mirror);
        poly_free(key_mirror_normalized);
    }
    poly_free(guess_key_normalized);

    return success;
}

enum Mode {FULL, PARTIAL, BOTH};


void reconstruct_full(uint16_t *delta, int8_t *S, poly_t h_normalized, int LENGTH, int WEIGHT) {
    count_check_mult = fail_check_mult = count_guess_key = count_nb_ext = 0;

    clock_t start = clock();
    poly_t key = guess_key(delta, S, LENGTH, WEIGHT, WEIGHT);
    clock_t end = clock();

    if (is_key(key, h_normalized)) {
        printf("success full\t%d\t%d\t%d\t%d\t%d\t%g ms\n",
            key->alloc, // max depth
            count_guess_key, // number of recursive calls
            count_nb_ext, // number of extensions
            count_check_mult, // number of checks
            fail_check_mult, // number of failed checks
            ((double) end - start) / CLOCKS_PER_SEC * 1000);
    } 
    
    else {
        printf("fail full\t%d\t%d\t%d\t%d\t%g ms\n",
            count_guess_key,
            count_nb_ext,
            count_check_mult,
            fail_check_mult,
            ((double) end - start) / CLOCKS_PER_SEC * 1000);
    }

    poly_free(key);
}

void reconstruct_partial(uint16_t *delta, int8_t *S, int8_t *S_high, int8_t *S_low, 
                         int cpt_zero, int cpt_ones, poly_t h_normalized, float percent_out, 
                         float percent_in) {
    int r = h_normalized->length;
    int d = h_normalized->weight;
    count_check_mult = fail_check_mult = count_guess_key = count_nb_ext = 0;
    int nb_out = (float)cpt_zero * (percent_out / 100.0);
    int nb_in  = (float)cpt_ones  * (percent_in / 100.0);

    uint16_t S_out[cpt_zero];
    spectrum_high(S_high, S_out, S, cpt_zero, r, nb_out);

    uint16_t S_in[cpt_ones];
    spectrum_low(S_low, S_in, S, cpt_ones, r, nb_in);

    clock_t start = clock();
    poly_t key = guess_key_bounds(delta, S_high, S_low, r, d, d);
    clock_t end = clock();

    if (is_key(key, h_normalized)) {
        printf("\n");
        printf("success partial\t%d\t%d\t%d\t%d\t%d\t%g ms\n",
               key->alloc, // max depth
               count_guess_key, // number of recursive calls
               count_nb_ext, // number of recursive calls
               count_check_mult, // number of checks
               fail_check_mult, // number of failed checks
               ((double) end - start) / CLOCKS_PER_SEC * 1000);
    } 
    
    else {
        printf("fail partial\t%d\t%d\t%d\t%d\t%g ms\n",
               count_guess_key,
               count_nb_ext,
               count_check_mult,
               fail_check_mult,
               ((double) end - start) / CLOCKS_PER_SEC * 1000);
    }

    poly_free(key);
}

int main(int argc, char ** argv) {
    poly_t h;
    int nb_samples = 1;
    unsigned seed = 0;
    int LENGTH = 12323;
    int WEIGHT = 71;
    enum Mode mode = BOTH;
    float percent_out = 70.0;
    float percent_in  = 1.0;

    // Arguments : LENGTH WEIGHT nb_samples seed mode nb_dist_out_percent nb_dist_in_percent
    if (argc > 1) LENGTH = atoi(argv[1]);
    if (argc > 2) WEIGHT = atoi(argv[2]);
    if (argc > 3) nb_samples = atoi(argv[3]);
    if (argc > 4) seed = atoi(argv[4]);
    if (argc > 5) {
        if (strcmp(argv[5], "full") == 0) mode = FULL;
        else if (strcmp(argv[5], "partial") == 0) mode = PARTIAL;
        else if (strcmp(argv[5], "both") == 0) mode = BOTH;
        else {
            printf("Mode must be 'full', 'partial', or 'both'\n");
            exit(1);
        }
    }
    if (argc > 6) percent_out = atof(argv[6]);
    if (argc > 7) percent_in  = atof(argv[7]);

    if ((LENGTH & 1) == 0) {
        printf("Length must be odd\n");
        exit(1);
    }

    h = poly_alloc(LENGTH, WEIGHT);
    uint16_t delta[WEIGHT];
    int8_t S[(LENGTH + 1)/2];
    int8_t S_high[(LENGTH + 1)/2];
    int8_t S_low[(LENGTH + 1)/2];

    for (int count = 0; count < nb_samples; ++count) {
        srandom(seed + count);
        poly_rand(h);
      
        poly_t h_normalized = normalize(h);

        int cpt_zero = 0;
        int cpt_ones = 0;
        spectrum_wo_mult(S, &cpt_zero, &cpt_ones, h->index, LENGTH, WEIGHT);
        sample_delta(delta, S, (LENGTH + 1)/2, WEIGHT);
     
        printf("\n ********* Reconstruction %d *********\n\n", count);

        if (mode == FULL || mode == BOTH) {
            reconstruct_full(delta, S, h_normalized, LENGTH, WEIGHT);
        }

        if (mode == PARTIAL || mode == BOTH) {
            reconstruct_partial(delta, S, S_high, S_low, cpt_zero, cpt_ones,
                                h_normalized, percent_out, percent_in);
        }

        poly_free(h_normalized);
    }
    poly_free(h);

    return 0;
}