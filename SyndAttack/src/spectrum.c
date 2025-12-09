#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "spectrum.h"
#include "sparse_cyclic.h"

int split_array(sparse_t array0, sparse_t array1, sparse_t array, index_t weight, index_t length0) {
    index_t i, weight0;
    for (i = 0; (i < weight) && (array[i] < length0); ++i) {
        array0[i] = array[i];
    }
    weight0 = i;
    for (; i < weight; ++i) {
        array1[i - weight0] = array[i] - length0;
    }
    return weight0;
}

void spectrum(int *S, sparse_t array, index_t weight, index_t length) {
    memset(S, 0, sizeof(int) * (length / 2));
    for (index_t i = 0; i < weight; ++i) {
        for (index_t j = i + 1; j < weight; ++j) {
            int l = array[j] - array[i];
            if (l > length / 2)
                l = length - l;
            S[l - 1]++;
        }
    }
}

int spectrum_multiplicity_profile(int *M, sparse_t array, index_t weight, index_t length) {
    memset(M, 0, sizeof(int) * (weight + 1));
    int i, max = 0, S[length / 2];
    spectrum(S, array, weight, length);
    for (i = 0; i < length / 2; ++i) {
        M[S[i]]++;
        if (S[i] > max)
            max = S[i];
    }
    return max;
}

int spectrum_size(sparse_t array, index_t weight, index_t length) {
    int i, w, S[length / 2];
    spectrum(S, array, weight, length);
    for (i = 0, w = 0; i < length / 2; ++i) {
        if (S[i])
            w++;
    }
    return w;
}

int spectrum_L2(sparse_t array, index_t weight, index_t length) {
    int i, w, S[length / 2];
    spectrum(S, array, weight, length);
    for (i = 0, w = 0; i < length / 2; ++i) {
        w += S[i] * S[i];
    }
    return w;
}

void cross_spectrum(int *S, sparse_t array0, index_t weight0, sparse_t array1, index_t weight1, index_t length) {
    memset(S, 0, sizeof(int) * (length * 2));
    for (index_t i = 0; i < weight0; ++i) {
        for (index_t j = 0; j < weight1; ++j) {
            int l = length + array1[j] - array0[i];
            S[l - 1]++;
        }
    }
}

int spectrum_distance(int *D, int *S0, int *S1, index_t size) {
    int w = 0;
    for (index_t i = 0; i < size; ++i) {
        int l = S0[i] * S1[i];
        if (D != NULL)
            D[i] = l;
        w += l;
    }
    return w;
}

int error_key_distance(sparse_t e, index_t t, sparse_t h0, sparse_t h1, index_t d, index_t r) {
    sparse_t e0 = sparse_new(t);
    sparse_t e1 = sparse_new(t);
    int t0 = split_array(e0, e1, e, t, r);
    int t1 = t - t0;
    int w = 0;
    int Sp_e[2 * r];
    int Sp_h[2 * r];

    if (h0) {
        spectrum(Sp_e, e0, t0, r);
        spectrum(Sp_h, h0, d, r);
        w += spectrum_distance(NULL, Sp_e, Sp_h, r / 2);
    }
    if (h1) {
        spectrum(Sp_e, e1, t1, r);
        spectrum(Sp_h, h1, d, r);
        w += spectrum_distance(NULL, Sp_e, Sp_h, r / 2);
    }
    /*
      cross_spectrum(Sp_h, h0, d, h1, d, r);
      cross_spectrum(Sp_e, e0, t0, e1, t1, r);
      w += spectrum_distance(NULL, Sp_e, Sp_h, r * 2);
    */
    sparse_free(e0);
    sparse_free(e1);
    return w;
}

void spectral_product(int *S, sparse_t e, int t, sparse_t h, int w, int r) {
    int i, j, l;

    memset(S, 0, sizeof(int) * (1 + (r - 1) / 2));
    for (i = 0; i < t; ++i) {
        for (j = 0; j < w; ++j) {
            l = (r + e[i] - h[j]) % r;
            if (2 * l < r) {
                S[l]++;
            }
        }
    }
}

int error_pscw_distance(sparse_t e, index_t t, sparse_t h0, sparse_t h1, index_t d, index_t r) {
    sparse_t e0 = sparse_new(t);
    sparse_t e1 = sparse_new(t);
    int t0 = split_array(e0, e1, e, t, r);
    int t1 = t - t0;
    int w = 0;
    int S0[1 + (r - 1) / 2], S1[1 + (r - 1) / 2];

    spectral_product(S0, e0, t0, h0, d, r);
    spectral_product(S1, e1, t1, h1, d, r);

    for (int l = 0; 2 * l < r; ++l) {
        if (S0[l] > w)
            w = S0[l];
        if (S1[l] > w)
            w = S1[l];
    }

    return w;
}

/* min(|i - j|, r - |i - j|)
   assumes 0 <= i,j < r (no need to reduce abs(i - j) modulo r) */
int distance(int i, int j, int r) {
    return ((r - abs(r - 2 * abs(i - j))) / 2) - 1;
}

/* checks whether distance(i,j,r) appears in the spectrum S
   assumes 0 <= i,j < r */
bool check_distance(int i, int j, int8_t * S, int r) {
    if (i == j) {
        return true;
    }

    else { 
        int d = distance(i, j, r);
        return (S[d] != 0);
    }
}

void spectrum_r(int8_t * S, uint16_t * index, int r, int d) {
    memset(S, 0, (r / 2) * sizeof (int8_t));
    
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < i; ++j) {
            S[distance(index[i], index[j], r)]++;
        }
    }  
}

void spectrum_wo_mult(int8_t * S, int * cpt_zero, int * cpt_ones, uint16_t * index, int r, int d) {
    memset(S, 0, (r / 2) * sizeof (int8_t));

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < i; ++j)
            if (S[distance(index[i], index[j], r)] == 0) {
                S[distance(index[i], index[j], r)]++;
            }
    
    for (int i = 0; i < r / 2; i++) {
        if (S[i] == 0) {
            (*cpt_zero)++;
        }
        else {
            (*cpt_ones)++;
        }
    }
}   

void spectrum_high(int8_t *S_high, uint16_t * S_out, int r, int nb_out) {
    memset(S_high, 1, (r / 2) * sizeof (int8_t));

    for (int i = 0; i < nb_out; i++)
        S_high[S_out[i]] = 0;

}

void spectrum_low(int8_t *S_low, uint16_t * S_in, int r, int nb_in) {
    memset(S_low, 0, (r / 2) * sizeof (int8_t));

    for (int i = 0; i < nb_in; i++)
        S_low[S_in[i]] = 1;

}

double normal_pdf(double x, double mu, double sigma) {
    return (1.0 / (sigma * sqrt(2 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2));
}

//Compute the probability of a distance to belong to each class (multiplicity)
void bayesian_inference(double x, double *means, double *std_dev, long double *count_dist, 
                        int num_classes, double *posteriors) {
    double denominator = 0.0;
    double f_k[num_classes];

    for (int i = 0; i < num_classes; i++) { 
        f_k[i] = normal_pdf(x, means[i], sqrt(std_dev[i]));
        denominator += f_k[i] * count_dist[i];
    }

    for (int i = 0; i < num_classes; i++)
        posteriors[i] = (f_k[i] * count_dist[i]) / denominator;  
    
}

/*
    Compute spectrum using soft information: 
    Each line i corresponds to a distance 
    Each column j corresponds to the probability of the multiplicity j
*/
void spectrum_soft_info(double S[BLOCK_LENGTH/2][max_dist_h0], double * means, double * avg_dist, 
                        double * std_devs, long double *count_dist, int r) {
    double posteriors[max_dist_h0];

    for (int i = 0; i < r/2; i++) {
        bayesian_inference(means[i], avg_dist, std_devs, count_dist, max_dist_h0, posteriors);
        for (int j = 0; j < max_dist_h0; j++) {
            S[i][j] = posteriors[j];
        }
    }
}







