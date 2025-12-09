#include <math.h>
#include <stdio.h>

#include "model.h"

/* log(n choose k) */
long double log_nCk_ld(int n, int k) {
    if (k < 0 || k > n) return -INFINITY;
    if (k > n - k) k = n - k;
    return lgammal((long double)n + 1.0L)
         - lgammal((long double)k + 1.0L)
         - lgammal((long double)(n - k) + 1.0L);
}

long double S_n(int n, int w, int t) {
    int limit = (w < t) ? w : t;
    long double max_log = -INFINITY;

    for (int i = 1; i <= limit; i += 2) {
        int ti = t - i;
        if (ti < 0 || ti > n - w) continue;
        long double lnum = log_nCk_ld(w, i) + log_nCk_ld(n - w, ti);
        if (lnum > max_log) max_log = lnum;
    }

    if (max_log == -INFINITY) return 0.0L;

    long double sum = 0.0L;
    for (int i = 1; i <= limit; i += 2) {
        int ti = t - i;
        if (ti < 0 || ti > n - w) continue;
        long double lnum = log_nCk_ld(w, i) + log_nCk_ld(n - w, ti);
        sum += expl(lnum - max_log);
    }

    long double logsumexp = max_log + logl(sum);
    long double logden = log_nCk_ld(n, t);
    long double result = expl(logsumexp - logden);

    return result;
}

double W(long int n, int w, int t, int nu) {
    int r = n / 2;
    int d =  w / 2;

    return ((r - 2 * (d - nu) - nu ) * S_n(n - 2, w, t - 2) + 
            2 * (d - nu) * (1 - S_n(n - 2, w - 1, t -2)) + 
            nu * S_n(n - 2, w - 2, t - 2));
}

/*
  Note: in the paper, the model is adjusted by centering the observations subtracting the observed 
  average syndrome weight and the model subtracting the average syndrome weight. Here, we only 
  center the observations subtracting both the average syndrome weight and the observed average 
  syndrome weight but it amounts to the same thing.
*/
void shifting_values(double *syndromes, long int nb_dec) {
    long double som = 0.0;
    double coeff = 1.0 / nb_dec;
    double rS_n = BLOCK_LENGTH * S_n(2 * BLOCK_LENGTH, 2 * BLOCK_WEIGHT, ERROR_WEIGHT);

    for (long int i = 0; i < nb_dec; i++)
        som += coeff * (rS_n - syndromes[i]);


    for (int i = 0; i < BLOCK_LENGTH / 2; i++)
        means[i] += som;
}

/*
  Compute parameters of the model: mean, standard deviation and multiplicity distribution.
*/
void compute_model(double avg_dist_model[max_dist_h0], double std_dist_model[max_dist_h0], 
                   long double count_dist_model[max_dist_h0], double noise, long int nb_dec) {
    double sigma0 = 83.3;
	double std_val = (sigma0*sigma0 + noise*noise) / nb_dec;
    long double logC_N1 = log_nCk_ld(BLOCK_LENGTH - 1, BLOCK_WEIGHT - 1);

    for (int m = 0; m < max_dist_h0; m++) {
        std_dist_model[m] = std_val;
        avg_dist_model[m] = W(2*BLOCK_LENGTH, 2*BLOCK_WEIGHT, ERROR_WEIGHT, m);
        count_dist_model[m] = expl(log_nCk_ld(BLOCK_WEIGHT, m) 
                            + log_nCk_ld(BLOCK_LENGTH - BLOCK_WEIGHT - 1, BLOCK_WEIGHT - m - 1) - logC_N1);
    }
}

/*
  Compute the set of distances believed to be outside the spectrum: distances whose syndrome weight is above the high limit and
  the set of distances believed to be inside the spectrum: distances whose syndrome weight is below the low limit. 
*/
void compute_S_out_S_in(uint16_t S_out[BLOCK_LENGTH/2], uint16_t S_in[BLOCK_LENGTH/2], int *cpt_avg, 
                        int *cpt_avg_l, double low_limit, double high_limit) {
    int cpt_wrong = 0;
    int cpt_wrong_l = 0;

    for (int j = 0; j < BLOCK_LENGTH / 2; ++j) {
        double x = means[j];

        if (x < low_limit) {
            S_in[(*cpt_avg_l)++] = j;
            if (S_h[0][j] == 0)
                cpt_wrong_l++;
        }

        if (x > high_limit) {
            S_out[(*cpt_avg)++] = j;
            if (S_h[0][j] != 0)
                cpt_wrong++;
        }    
    }

    if (cpt_wrong != 0 || cpt_wrong_l != 0) {
        printf("Warning ! \n");
        printf("%d distance(s) with mult > 0 above limit\n", cpt_wrong);
        printf("%d distance(s) with mult = 0 below limit\n", cpt_wrong_l);
    }
}