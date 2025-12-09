#include <stdint.h>
#include <stdbool.h>

int distance(int i, int j, int r);
bool check_distance(int i, int j, int8_t *S, int r);
void spectrum(int8_t *S, uint16_t *index, int r, int d);
void spectrum_wo_mult(int8_t *S, int *cpt_zero, int *cpt_ones, uint16_t *index, int r, int d);
void spectrum_high(int8_t *S_high, uint16_t *S_out, int8_t *S, int cpt_zero, int r, int nb_out); 
void spectrum_low(int8_t *S_low, uint16_t *S_in, int8_t *S, int cpt_ones, int r, int nb_in);