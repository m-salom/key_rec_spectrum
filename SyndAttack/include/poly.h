#pragma once
#include <stdbool.h>

#include "param.h"
#include "types.h"

int compare_floats(const void *a, const void *b);
int compare_decrease(const void *a, const void *b);
poly_t poly_alloc(int length, int weight);
void poly_rand(poly_t h);
void poly_free(poly_t h);
int extended_weight(poly_t h);
bool extended_coeff_is_zero(poly_t h, int i);
poly_t extend_binomial_init(int j, int8_t * S, int length);
poly_t poly_shift_union(poly_t g, poly_t f, int l);
poly_t poly_shift_union_soft(double * res, uint16_t* V, poly_t g, poly_t f, 
                             int l, double S[BLOCK_LENGTH/2][max_dist_h0]);
double log_score_poly(poly_t h, double S[BLOCK_LENGTH/2][max_dist_h0]);
