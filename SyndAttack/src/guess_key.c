#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <omp.h>
#include <sys/time.h>

#include "model.h"
#include "poly.h"
#include "spectrum.h"
#include "synd_io.h"
#include "guess_key.h"

double get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return t.tv_sec + t.tv_usec * 1e-6;
}

//Comparison functions
extern const double *global_arrays;

int compare(const void *a, const void *b) {
    int idx1 = *(const int *)a;
    int idx2 = *(const int *)b;
    return (global_arrays[idx1] > global_arrays[idx2]) - (global_arrays[idx1] < global_arrays[idx2]);
}

static int int16comp (const void * a, const void * b) {
    return *((uint16_t *) a) - *((uint16_t *) b);
}

typedef struct {
    poly_t poly;
    double ext_w;
    double ext_wm1;
    uint16_t ext_monomials[BLOCK_WEIGHT];
} offset_score_t;

static int compare_offset_ext_desc(const void *a, const void *b) {
    double sa = ((offset_score_t *)a)->ext_wm1;
    double sb = ((offset_score_t *)b)->ext_wm1;
    return (sb > sa) - (sb < sa);
}

extern volatile sig_atomic_t interrupted;
int count_check_mult = 0, fail_check_mult = 0, count_guess_key = 0, count_nb_ext = 0;

/* put in the array res[] all the offsets s such that
   Spectrum(g union x^s+x^(s+delta)) is included in S
   returns the number of such s */
int valid_offset(uint16_t * res, poly_t g, uint16_t delta, int8_t * S) {
    int i, j, l, k = 0, r = g->length;
    for (l = 0; l < r; ++l) {
        j = (l < r - l) ? l : (r - l);
        if ((l == 0) || (S[j-1] > 0)) {
            j = (l + delta) % r;
            for (i = 0; i < g->weight; ++i) {
                if (check_distance(l, g->index[i], S, r) &&
                    check_distance(j, g->index[i], S, r))
                    continue;
                break;
            }
            if (i == g->weight)
                res[k++] = l;
        }
    }
    return k;
}

int valid_offset_soft(uint16_t * V, poly_t g, uint16_t delta, double S[BLOCK_LENGTH/2][max_dist_h0], int depth, double *T) {
    int i, j, l, k = 0, r = g->length;
    int d = g->weight;
    double res[r];
    poly_t h = poly_alloc(r, d+2);
    memcpy(h->index, g->index, d * sizeof(uint16_t));
    
    for (l = 0; l < r; ++l) {
        j = (l < r - l) ? l : (r - l);

        if ((l == 0) || (S[j-1][0] != 1.)) {
            j = (l + delta) % r;
            h->index[d] = l;
            h->index[d+1] = j;
            res[l] = log_score_poly(h, S);
        }

        else {
            res[l] = -10e10;
        }
    }

    int index[BLOCK_LENGTH];
    global_arrays = res; 
    for (i = 0; i < BLOCK_LENGTH; i++) index[i] = i; 
    qsort(index, BLOCK_LENGTH, sizeof(int), compare_decrease); 

    for (i = 0; i < r; i++) V[i] = index[i];
    qsort(res, r, sizeof(double), compare_floats);

    /*
    Global threshold for the score of polynomials (below this score, it is not likely we are exploring a promising branch.)
    */
    if (res[0] < -85.) return 0;
    if (depth == 1) return 30;
    if (depth == 2) return 10;
    
    if (T != NULL && depth >= 3 && depth < 12) {
        if (res[0] < T[depth]) // thresholds for different sets of parameters in guess_key.h
            return 0;
    }

    for (i = 0; i < 10; i++) {
        if (res[i] - res[i+1] > 1.5) {
            k = i+1;
            break;
        }
    }

    if (k==0) {
        for (i = 0; i < r-1; i++) {
            if (res[i] < -2.) {
                k = i;
                break;
            }
        }
    }

    return k;
}

/* extract the extension of h into g
   returns g if its spectrum is exactly S (without multiplicity)
   else returns NULL */
poly_t check_key(poly_t h, int8_t * S) {
    int i, j, w, r = h->length;
    uint16_t index[r]; // oversized
    int8_t mult[r / 2];
    memset(mult, 0, (r / 2) * sizeof (int8_t));
    for (i = 0, w = 0; i < r; ++i) 
        if (!extended_coeff_is_zero(h, i))
            index[w++] = i;
    for (i = 0; i < w; ++i) {
        for (j = 0; j < i; ++j) {
            if (mult[distance(index[i], index[j], r)] == 0) {
                mult[distance(index[i], index[j], r)] += 1;
            }
        }
    }

    if (memcmp(S, mult, (r / 2) * sizeof (int8_t)) != 0) {
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

/* returns true if the spectrum of the extension of h contains the
   spectrum S */
bool check_mult(poly_t h, int8_t * S) {
    int i, j, w, r = h->length;
    uint16_t index[r];
    int8_t mult[r / 2];

    count_check_mult++;

    for (i = 0, w = 0; i < r; ++i)
        if (!extended_coeff_is_zero(h, i))
            index[w++] = i;
    memset(mult, 0, (r / 2) * sizeof (mult[0]));

    for (i = 0; i < w; ++i) {
        for (j = 0; j < i; ++j) {
            mult[distance(index[i], index[j], r)] += 1;
        }
    }

    for (i = 0; i < r / 2; ++i) {
        if ((S[i] > 0) && (mult[i] == 0)) {
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
        //current polynomial h is built from g and f[depth] shifted by V[i]
        h = poly_shift_union(g, f[depth], V[i]);
    
        // weight of the extension of h
        int w = extended_weight(h);
    
        // cannot lead to a word of target weight -> quit and try next entry of V
        if (w < weight) {
            poly_free(h);
            continue;
        }

        //no word deriving from h can reach the target spectrum multiplicity
        if (check_mult(h, S) == false) {
            poly_free(h);
            continue;
        }

        //if w equals the target weight we check the only possible solution
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
        
        // no solution, no contradiction, we continue from h and the next delta
        key = guess_key_rec(delta, f, depth+1, max_depth, h, S, weight);
    
        poly_free(h);
        if (key != NULL)
            break;
    }

    return key;
}

poly_t guess_key_rec_bounds(uint16_t * delta, poly_t * f, int depth, int max_depth, poly_t g, int8_t * S, int8_t * S_low, int weight) {
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

        //no word deriving from h can reach the target spectrum multiplicity
        if (check_mult(h, S_low) == false) {
            poly_free(h);
            continue;
        }

        //if w equals the target weight we check the only possible solution
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

poly_t guess_key_rec_soft(uint16_t * delta, poly_t * f, int depth, int max_depth, poly_t g, 
                          double S[BLOCK_LENGTH/2][max_dist_h0], int weight, double gap_threshold,
                          double *T) {

    int nb_offset, r = g->length;
    poly_t key = NULL;
    uint16_t V[r];
  
    if (depth >= max_depth)
        return NULL;

    if (f[depth] == NULL) { // initialize, only once per depth
        f[depth] = poly_alloc(r, 2);
        f[depth]->index[0] = 0;
        f[depth]->index[1] = delta[depth];
    }

    nb_offset = valid_offset_soft(V, g, delta[depth], S, depth, T);
 
    offset_score_t offset_scores[nb_offset];

    #pragma omp parallel for
    for (int i = 0; i < nb_offset; ++i) {
        double Ext[r];
        uint16_t Ext_monomials[r];

        poly_t p = poly_shift_union_soft(Ext, Ext_monomials, g, f[depth], V[i], S);

        offset_scores[i].poly = p;
        offset_scores[i].ext_w = (p != NULL) ? Ext[weight] : -1e11;
        offset_scores[i].ext_wm1 = (p != NULL) ? Ext[weight - 1] : -1e11;

        if (p != NULL) {
            memcpy(offset_scores[i].ext_monomials, Ext_monomials, sizeof(uint16_t) * weight);
        } 
        
        else {
            memset(offset_scores[i].ext_monomials, 0, sizeof(uint16_t) * weight);
        }
    }

    qsort(offset_scores, nb_offset, sizeof(offset_score_t), compare_offset_ext_desc);

    for (int i = 0; i < nb_offset; ++i) {
        offset_score_t *o = &offset_scores[i];
        poly_t p = o->poly;

        if (p == NULL)
            continue;

        if (o->ext_wm1 - o->ext_w > gap_threshold) {
            key = poly_alloc(r, weight);
            for (int j = 0; j < weight; j++) {
                key->index[j] = o->ext_monomials[j];
            }
            break;
        }

        if (p->weight == weight) {
            key = poly_alloc(r, weight);
            for (int j = 0; j < weight; j++) {
                key->index[j] = p->index[j];
            }
            break;
        }

        if (interrupted) return NULL;
        key = guess_key_rec_soft(delta, f, depth + 1, max_depth, p, S, weight, gap_threshold, T);
        poly_free(p);
        if (key != NULL) break;
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

poly_t guess_key_soft(uint16_t * delta, double S[BLOCK_LENGTH/2][max_dist_h0], int r, int max_depth, int weight, double gap_threshold, double *T) {
    int i;
    poly_t f[max_depth];
    poly_t h;

    for (i = 0; i < max_depth; ++i) {
        f[i] = NULL;
    }

    f[0] = poly_alloc(r, 2);
    f[0]->index[0] = 0;
    f[0]->index[1] = delta[0];
    
    h = guess_key_rec_soft(delta, f, 1, max_depth, f[0], S, weight, gap_threshold, T);
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
        delta[i] = j + 1;
    }
}

void sample_delta_soft(uint16_t *delta, double *S, int nb_pos) {
    
    int index[(BLOCK_LENGTH)/2];

    global_arrays = S; 
    for (int i = 0; i < (BLOCK_LENGTH) / 2; i++) index[i] = i; 
    qsort(index, (BLOCK_LENGTH) / 2, sizeof(int), compare); 

    for (int i = 0; i < nb_pos; i++) {
        delta[i] = index[i]+1;
    }
}

poly_t normalize(poly_t h) {
    int i, j, l, offset = 0, r = h->length, d = h->weight;
   
    poly_t h_normalized = poly_alloc(r, d);
    int8_t S[r / 2];

    spectrum_r(S, h->index, r, d);
    /* l = smallest distance with multiplicity 1 */
    /* conjecture: at least one distance has multiplicity 1 */
    for (l = 0; S[l] != 1; ++l);

    for (i = 0; i < d; ++i) {
        for (j = 0; j < i; ++j) {
            if (distance(h->index[i], h->index[j], r) == l) {
                break;
            }   
        }
      
        if (i != j) {
            if (h->index[i] == (h->index[j] + l + 1) % r)
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

poly_t mirror_poly(poly_t p) {
    int w = p->weight;
    int r = p->length;

    poly_t m = poly_alloc(r,w);
    if (!m) return NULL;

    for (int i = 0; i < w; i++) {
        m->index[i] = (r - p->index[i]) % r;
    }

    return m;
}

bool is_key(poly_t guess_key, poly_t key_normalized) {
    if (guess_key == NULL) return false;
    poly_t guess_key_normalized = normalize(guess_key);
    bool success = true;

    for (int i = 0; i < BLOCK_WEIGHT; i++) {
        if (guess_key_normalized->index[i] != key_normalized->index[i]){
            success = false;
            break;
        }
    }

    if (!success) {
        success = true;
        poly_t key_mirror = mirror_poly(guess_key);
        poly_t key_mirror_normalized = normalize(key_mirror);
        for (int i = 0; i < BLOCK_WEIGHT; i++) {
            if (key_mirror_normalized->index[i] != key_normalized->index[i]) {
                success = false;
                break;
            }
        }
        poly_free(key_mirror);
        poly_free(key_mirror_normalized);
    }
    
    poly_free(guess_key_normalized);
    poly_free(guess_key);

    return success;
}

void reconstruct_partial(poly_t h, poly_t h_normalized) {
    count_check_mult = fail_check_mult = count_guess_key = count_nb_ext = 0;
    int8_t S_high[BLOCK_LENGTH/2];
    int8_t S_low[BLOCK_LENGTH/2];
    uint16_t S_out[BLOCK_LENGTH/2];
    uint16_t S_in[BLOCK_LENGTH/2];
    uint16_t delta[BLOCK_WEIGHT];
    int cpt_avg = 0;
    int cpt_avg_l = 0;
    // These limits are used to recover the key with a partial spectrum with probability 0.95 using 750000 queries. 
    // They can be adapted depending on the probability and the amount of queries.
    double low_limit = 4867.9; 
    double high_limit = 4868.933; 
    int8_t S[BLOCK_LENGTH/2];
    int cpt_zero = 0;
    int cpt_ones = 0;
    
    spectrum_wo_mult(S, &cpt_zero, &cpt_ones, h->index, BLOCK_LENGTH, BLOCK_WEIGHT);
    compute_S_out_S_in(S_out, S_in, &cpt_avg, &cpt_avg_l, low_limit, high_limit);
    sample_delta(delta, S, BLOCK_LENGTH / 2, BLOCK_WEIGHT);
    
    printf("Number of 0 in S_high : %d = %f %% \n", cpt_avg, ((float)cpt_avg/cpt_zero)*100.);
    printf("Number of 1 in S_low : %d = %f %% \n", cpt_avg_l, ((float)cpt_avg_l/cpt_ones)*100.);
    spectrum_high(S_high, S_out, BLOCK_LENGTH, cpt_avg);
    spectrum_low(S_low, S_in, BLOCK_LENGTH, cpt_avg_l);
    printf("\n ********* Reconstruction from partial spectrum *********\n\n");
    double t0 = get_time();
    poly_t key = guess_key_bounds(delta, S_high, S_low, BLOCK_LENGTH, BLOCK_WEIGHT, BLOCK_WEIGHT);
    double t1 = get_time();

    if (is_key(key, h_normalized)) {
        printf("\n");
        printf("success partial\t%d\t%d\t%d\t%d\t%d\t%g ms\n",
               key->alloc, // max depth
               count_guess_key, // number of recursive calls
               count_nb_ext, // number of recursive calls
               count_check_mult, // number of checks
               fail_check_mult, // number of failed checks
               (t1-t0)*1000);
    } 
    
    else {
        printf("fail partial\t%d\t%d\t%d\t%d\t%g ms\n",
               count_guess_key,
               count_nb_ext,
               count_check_mult,
               fail_check_mult,
               (t1-t0)*1000);
    }
}

void reconstruct_soft(long int nb_dec, double sigma, poly_t h_normalized) {
    double gap_threshold = 6.0;
    double avg_dist_model[max_dist_h0];
    double std_dist_model[max_dist_h0];
    long double count_dist_model[max_dist_h0];
    double S_soft[BLOCK_LENGTH/2][max_dist_h0];
    uint16_t delta[BLOCK_WEIGHT];

    compute_model(avg_dist_model, std_dist_model, count_dist_model, sigma, nb_dec);
    spectrum_soft_info(S_soft, means, avg_dist_model, std_dist_model, count_dist_model, BLOCK_LENGTH);
    sample_delta_soft(delta, means, BLOCK_WEIGHT);
    printf("\n ********* Reconstruction from soft spectrum *********\n\n");

    double T_l[12];
    double *T = NULL;
 
    if (sigma == 0 && nb_dec == 200000) {
        memcpy(T_l, thresholds_200k, sizeof(T_l));
        T = T_l;
    } 
    else if (sigma == 0 && nb_dec == 300000){ 
        memcpy(T_l, thresholds_300k, sizeof(T_l));
        T = T_l;
    }
    else if (sigma == 10 && nb_dec == 200000) {
        memcpy(T_l, thresholds_200k_sigma10, sizeof(T_l));
        T = T_l;
    }
    else if (sigma == 30 && nb_dec == 280000) {
        memcpy(T_l, thresholds_280k_sigma30, sizeof(T_l));
        T = T_l;
    }
    else if (sigma == 50 && nb_dec == 400000) {
        memcpy(T_l, thresholds_400k_sigma50, sizeof(T_l));
        T = T_l;
    }
   
    double t0 = get_time();
    poly_t key = guess_key_soft(delta, S_soft, BLOCK_LENGTH, BLOCK_WEIGHT, BLOCK_WEIGHT, gap_threshold, T);
    double t1 = get_time();

    if (is_key(key, h_normalized)) {
        printf("\n");
        printf("success soft\t%g s\n",
               t1-t0);
    } 
    
    else {
        printf("fail soft\t%g s\n",
               t1-t0);
    }
}