#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>

#include "model.h"
#include "code.h"
#include "codegen.h"
#include "errorgen.h"
#include "param.h"
#include "types.h"
#include "sparse_cyclic.h"
#include "poly.h"
#include "guess_key.h"
#include "synd_io.h"
#include "spectrum.h"


volatile sig_atomic_t interrupted = 0;

void handle_sigint(int sig) {
    (void)sig;
    interrupted = 1;
}

int synd_attack(double * syndromes, uint16_t * index, long key_seed, long error_seed, long nb_dec, double sigma) {
    long int i;
    bool abort = false;
    int error_weight = errgen_weight();
    int blocklength = keygen_blocklength();
    code_t H;
    e_t e;
    syndrome_t syndrome;
    index_t error_sparse[error_weight];
        
    keygen(&H, key_seed);
    
    for (int i = 0; i < BLOCK_WEIGHT; i++)
        index[i] =  H.rows[0][i];

    init_synd_attack(&H, blocklength, error_weight, key_seed, error_seed);

    for (i = 0; ((nb_dec < 0) ||  (i < nb_dec)) && (!abort); ++i) {
        errgen_sparse(error_sparse, &H, error_seed + i, key_seed);
        error_sparse_to_dense(&e, error_sparse, error_weight, BLOCK_LENGTH);
        compute_syndrome(&syndrome, &H, &e);
        
        log_synd_attack(syndromes, error_sparse, syndrome.weight, i, sigma);

        abort = catch_interrupt();
    }

    display_synd_attack(true, sigma);

    close_synd_attack();
    
    return abort;
}


int main(int argc, char *argv[]) {
        
    bool abort;
    long seed = -1;
    long error_seed = -1;
    double sigma = 0;
    long int nb_dec = 1000000; // default number of decoding samples before log or exit
    bool once = true;         // if false continue until interrupted
    enum Mode mode = BOTH;
    uint16_t index[BLOCK_WEIGHT];
    poly_t h = poly_alloc(BLOCK_LENGTH, BLOCK_WEIGHT);
    poly_t h_normalized = poly_alloc(BLOCK_LENGTH, BLOCK_WEIGHT);

    parse_arguments_and_init_synd(argc, argv, &seed, &error_seed, &nb_dec, &once, &sigma, &mode);

    double *syndromes = malloc(nb_dec * sizeof(double));
    
    do {
        abort = synd_attack(syndromes, index, seed, error_seed, nb_dec, sigma);
        seed += 1;
    } while (!(once || abort));

    for (int i = 0; i < BLOCK_WEIGHT; i++) {
        h->index[i] = index[i];
    }

    h_normalized = normalize(h);
    shifting_values(syndromes, nb_dec);
 
    signal(SIGINT, handle_sigint);
    

    if (mode == PARTIAL || mode == BOTH) {
        reconstruct_partial(h, h_normalized); 
    }

    if (mode == SOFT || mode == BOTH) {
        reconstruct_soft(nb_dec, sigma, h_normalized); 
    }

    poly_free(h);
    poly_free(h_normalized);

    exit(EXIT_SUCCESS);
}
