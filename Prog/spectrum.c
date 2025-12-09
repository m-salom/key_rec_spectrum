#include <stdlib.h>
#include <string.h>

#include "spectrum.h"

/* min(|i - j|, r - |i - j|)
   assumes 0 <= i,j < r (no need to reduce abs(i - j) modulo r) */
int distance(int i, int j, int r) {
    return (r - abs(r - 2 * abs(i - j))) / 2;
}

/* checks whether distance(i,j,r) appears in the spectrum S
   assumes 0 <= i,j < r */
bool check_distance(int i, int j, int8_t * S, int r) {
    int d = distance(i, j, r);
    return (d == 0) || (S[d] != 0);
}

void spectrum(int8_t * S, uint16_t * index, int r, int d) {
    memset(S, 0, ((r+1) / 2) * sizeof (int8_t));
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < i; ++j)
            S[distance(index[i], index[j], r)] += 1;
}

//Computes spectrum without multiplicity and counts number of distances with multiplicity 0 and > 0
void spectrum_wo_mult(int8_t * S, int * cpt_zero, int * cpt_ones, uint16_t * index, int r, int d) {
    memset(S, 0, ((r+1) / 2) * sizeof (int8_t));
    
    for (int i = 0; i < d; ++i) {
        for (int j = 0; j < i; ++j) {
#ifndef FULL_MULTIPLICITY
            if (S[distance(index[i], index[j], r)] == 0)
#endif
                S[distance(index[i], index[j], r)] += 1;
        }
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

void shuffle(uint16_t *array, size_t n) { 
    if (n > 1) { 
        for (size_t i = 0; i < n - 1; i++) { 
            size_t j = i + rand() / (RAND_MAX / (n - i) + 1); 
            uint16_t temp = array[j]; 
            array[j] = array[i]; 
            array[i] = temp; 
        } 
    } 
} 

void spectrum_high(int8_t *S_high, uint16_t * S_out, int8_t * S, int cpt_zero, int r, int nb_out) {
    memset(S_high, 1, ((r+1) / 2) * sizeof (int8_t));
    int i,j;
    uint16_t null_positions[cpt_zero];
    uint16_t numbers[cpt_zero]; 

    for (i = 0, j = 0; i < (r+1) / 2; i++) {
        if (S[i] == 0)
            null_positions[j++] = i;
    }

    for (i = 0; i < cpt_zero; i++)
        numbers[i] = null_positions[i]; 
    
    shuffle(numbers, cpt_zero); 
    
    for (i = 0; i < nb_out; i++)
        S_out[i] =  numbers[i];

    for (i = 0; i < nb_out; i++)
        S_high[S_out[i]] = 0;
}

void spectrum_low(int8_t *S_low, uint16_t * S_in, int8_t * S, int cpt_ones, int r, int nb_in) {
    memset(S_low, 0, ((r+1) / 2) * sizeof (int8_t));
    int i,j;
    uint16_t non_null_positions[cpt_ones];
    uint16_t numbers[cpt_ones]; 

    for (i = 0, j = 0; i < (r + 1) / 2; i++) {
        if (S[i] != 0)
            non_null_positions[j++] = i;
    }

    for (i = 0; i < cpt_ones; i++)
        numbers[i] = non_null_positions[i]; 

    shuffle(numbers, cpt_ones); 
    
    for (i = 0; i < nb_in; i++)
        S_in[i] =  numbers[i];

    for (i = 0; i < nb_in; i++)
        S_low[S_in[i]] = 1;
}



