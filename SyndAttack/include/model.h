#pragma once
#include <stdint.h>

#include "param.h"

extern double means[BLOCK_LENGTH/2];
extern int S_h[INDEX][BLOCK_LENGTH / 2];

void shifting_values(double *syndromes, long int nb_dec);
void compute_model(double avg_dist_model[max_dist_h0], double std_dist_model[max_dist_h0], 
                   long double count_dist_model[max_dist_h0], double noise, long int nb_dec);
void compute_S_out_S_in(uint16_t S_out[BLOCK_LENGTH/2], uint16_t S_in[BLOCK_LENGTH/2], 
                        int *cpt_avg, int *cpt_avg_l, double low_limit, double high_limit);