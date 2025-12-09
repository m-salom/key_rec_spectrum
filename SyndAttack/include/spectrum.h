#pragma once
#include <stdbool.h>

#include "types.h"

int split_array(sparse_t array0, sparse_t array1, sparse_t array, index_t weight, index_t length0);
void spectrum(int *S, sparse_t array, index_t weight, index_t length);
int spectrum_multiplicity_profile(int *M, sparse_t array, index_t weight, index_t length);
int spectrum_size(sparse_t array, index_t weight, index_t length);
int spectrum_L2(sparse_t array, index_t weight, index_t length);
int error_key_distance(sparse_t e, index_t t, sparse_t h0, sparse_t h1, index_t d, index_t r);
int distance(int i, int j, int r);
bool check_distance(int i, int j, int8_t * S, int r);
void spectrum_r(int8_t * S, uint16_t * index, int r, int d);
void spectrum_wo_mult(int8_t * S, int * cpt_zero, int * cpt_ones, uint16_t * index, int r, int d);
void spectrum_high(int8_t *S_high, uint16_t * S_out, int r, int nb_out); 
void spectrum_low(int8_t *S_low, uint16_t * S_in, int r, int nb_in);
void spectrum_soft_info(double S[BLOCK_LENGTH/2][max_dist_h0], double * means, double * avg_dist, 
                        double * std_devs, long double *count_dist, int r);