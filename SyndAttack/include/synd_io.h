#pragma once
#include <stdbool.h>

#include "types.h"

enum Mode {
    SOFT,
    PARTIAL,
    BOTH
};

bool catch_interrupt();

void init_synd_attack(code_t * H, int length, int e_weight, long key_seed, long error_seed);
void close_synd_attack();
void log_synd_attack(double * syndromes, sparse_t error_sparse, int syndrome_weight, long int nb_dec, double sigma);
void display_synd_attack(bool final, double sigma);
void parse_arguments_and_init_synd(int argc, char *argv[], long *seed, long *error_seed, long *nb_dec, bool *once, double *sigma, enum Mode *mode);
