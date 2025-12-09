#include <stdlib.h>
#include <signal.h>
#include <stdint.h>
#include <getopt.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "synd_io.h"
#include "types.h"
#include "spectrum.h"
#include "sparse_cyclic.h"
#include "codegen.h"
#include "errorgen.h"

unsigned verbose_synd = 1;

static bool log_file = false; // if log_file=true output into a file (else stdout)

long random_seed_synd() {
    FILE *urandom_fp;
    unsigned long s;

    urandom_fp = fopen("/dev/urandom", "r");
    if (urandom_fp == NULL)
        return 0;
    if (fread(&s, 8, 1, urandom_fp) != 1)
        return 0;
    fclose(urandom_fp);

    return s & ((1UL << 48) - 1); // 48 bits
}

#define CONT 0
#define STOP 1
#define SHOW 2

static int sig = CONT;

static void inthandler(int signo) {
    (void)signo;
    sig = STOP;
}

static void huphandler(int signo) {
    (void)signo;
    sig = SHOW;
}

static int blocklength, error_weight;
static int max_dist;
static long nb_samples;
static double total_weight;
static long key_seed, error_seed;
static long count_synd[INDEX][BLOCK_LENGTH / 2];
static double total_synd[INDEX][BLOCK_LENGTH / 2];
double total_synd2[INDEX][BLOCK_LENGTH / 2];
static int S_e[INDEX][BLOCK_LENGTH / 2];
int S_h[INDEX][BLOCK_LENGTH / 2];
static sparse_t e[INDEX];
static sparse_t h[INDEX];
double means[BLOCK_LENGTH/2];

void init_synd_attack(code_t * H, int length, int e_weight, long k_seed, long e_seed) {
    int i, j;
    key_seed = k_seed;
    error_seed = e_seed;
    blocklength = length;
    error_weight = e_weight;
    nb_samples = 0;
    total_weight = 0;

    for (i = 0; i < INDEX; ++i) {
        h[i] = H->rows[i];
     
        spectrum(S_h[i], h[i], BLOCK_WEIGHT, blocklength);
        memset(count_synd[i], 0, (BLOCK_LENGTH / 2) * sizeof (long));
        memset(total_synd[i], 0, (BLOCK_LENGTH / 2) * sizeof (double));
        memset(total_synd2[i], 0, (BLOCK_LENGTH / 2) * sizeof (double));

        e[i] = sparse_new(errgen_weight());
    }

    for (i = 0; i < INDEX; ++i)
        for (j = 0; j < blocklength / 2; ++j)
            if (S_h[i][j] > max_dist)
                max_dist = S_h[i][j];

    if (max_dist >= max_dist_h0)
        exit(0);

}

void close_synd_attack() {
    int i;
    for (i = 0; i < INDEX; ++i) {
        sparse_free(e[i]);
    }
}

bool catch_interrupt() {
    if (sig == STOP)
        return true;
    sig = CONT;
    return false;
}

double normal_random() {
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
}

// N(0, sigma^2)
double normal_random_sigma(double sigma) {
    return sigma * normal_random();
}

#define WEIGHTED_SYND_W

void log_synd_attack(double * syndromes, sparse_t error_sparse, int syndrome_weight, long int nb_dec, double sigma) {
    int i, j, t[INDEX];
    nb_samples++;

    if (sigma > 0) {
        double noise = normal_random_sigma(sigma);
        syndrome_weight += noise;
    }

    total_weight += syndrome_weight;
    // INDEX = 2
    t[0] = split_array(e[0], e[1], error_sparse, error_weight, blocklength);
    t[1] = error_weight - t[0];
    for (i = 0; i < INDEX; ++i) {
        spectrum(S_e[i], e[i], t[i], blocklength);
    }

    syndromes[nb_dec] = syndrome_weight;

    for (i = 0; i < INDEX; ++i) {
        for (j = 0; j < blocklength / 2; ++j) {
            if (S_e[i][j] > 0) {
#ifdef WEIGHTED_SYND_W // weighted average syndrome weight
                count_synd[i][j] += S_e[i][j];
                total_synd[i][j] += S_e[i][j] * syndrome_weight;
                total_synd2[i][j] += S_e[i][j] * syndrome_weight * syndrome_weight; 
#else // average syndrome weight
                count_synd[i][j] ++;
                total_synd[i][j] += syndrome_weight;
                total_synd2[i][j] += syndrome_weight * syndrome_weight;
#endif
            }
        }
    }
}

/* given (h_0,h_1) of weight w = 2*d, |h_0| = |h_1| = d
   average syndrome weight for a given distance (b,delta) b \in {0,1}, 0 < delta < (r + 1) / 2:
   sum(|e_0*h_0 + e_1*h_1|, |(e_0,e_1)| = t, delta \in Sp(e_b)) / {(e_0,e_1), |(e_0,e_1)| = t, delta \in Sp(e_b))
   weighted average syndrome weight:
   sum(mu(delta,e_b) * |e_0*h_0 + e_1*h_1|, |(e_0,e_1)| = t) / sum(mu(delta,e_b), |(e_0,e_1)| = t) 
*/
void display_synd_attack(bool final, double sigma) {
    int i, j;
    char logfilename[100];
    long count_dist[INDEX][max_dist_h0];
    double std_dist[INDEX][max_dist_h0];
    double avg_dist[INDEX][max_dist_h0];
    unsigned long N[2] = {0};
    double T[2] = {0};
    FILE *file = stdout;

    sprintf(logfilename, "synd_attack%ld-%ld", nb_samples, key_seed);
    
    if (log_file) {
        file = fopen(logfilename, "a");
        while (lockf(fileno(file), F_LOCK, 0) < 0)
            ;
    } else if (!final) {
        file = stderr;
    }

    for (i = 0; i < INDEX; ++i) {
        memset(count_dist[i], 0, max_dist_h0 * sizeof (long));
        memset(avg_dist[i], 0, max_dist_h0 * sizeof (double));
        memset(std_dist[i], 0, max_dist_h0 * sizeof (double));
    }

    fprintf(file, "# %d %d %d %ld\n", blocklength, BLOCK_WEIGHT, error_weight, nb_samples);
    fprintf(file, "# SIGMA : %f\n", sigma);
    fprintf(file, "# key seed: %lu\t error seed: %lu\n", key_seed, error_seed);

    if (verbose_synd & 1) {
        fprintf(file, "# block 0 spectrum size/L2: %d / %d\tblock 1 spectrum size/L2: %d / %d\n",
                spectrum_size(h[0], BLOCK_WEIGHT, blocklength),
                spectrum_L2(h[0], BLOCK_WEIGHT, blocklength),
                spectrum_size(h[1], BLOCK_WEIGHT, blocklength),
                spectrum_L2(h[1], BLOCK_WEIGHT, blocklength)
                );
#ifdef WEIGHTED_SYND_W // weighted average syndrome weight
    fprintf(file, "# weighted average syndrome weight (weighted with multiplicity in error)\n");
#else
    fprintf(file, "# average syndrome weight (when multiplicity in error > 0)\n");
#endif
    }

    for (i = 0; i < INDEX; ++i) {
        for (j = 0; j < blocklength / 2; ++j) {
            if (S_h[i][j] < max_dist_h0) {
                double x = total_synd[i][j] / count_synd[i][j];
                avg_dist[i][S_h[i][j]] += x;
                std_dist[i][S_h[i][j]] += x * x;
                count_dist[i][S_h[i][j]]++;
                N[i] += 1;
                T[i] += x;
            }
        }
        
        T[i] /= N[i];
    }

    if (verbose_synd & 2)
        for (i = 0; i < INDEX; ++i) {
            for (j = 0; j <= max_dist; ++j) {
                if (count_dist[i][j] > 0) {
                    avg_dist[i][j] /= count_dist[i][j];
                    std_dist[i][j] /= count_dist[i][j];
                    std_dist[i][j] -= avg_dist[i][j] * avg_dist[i][j];
                    if (count_dist[i][j] == 1)
                        std_dist[i][j] = 0;
                    fprintf(file, "# %d\t%d\t%-16.16g\t%5ld\t%-16.16g\n", i, j, avg_dist[i][j], count_dist[i][j], sqrt(std_dist[i][j]));
                }
            }
            fprintf(file, "# %d\t\t%-16.16g\t%5ld\n", i, T[i], N[i]);
        }
    fprintf(file, "# \t\t%-16.16g\n", total_weight / nb_samples);

    for (i = 0; i < INDEX; ++i) {
        for (j = 0; j < blocklength / 2; ++j) {
            double x = total_synd[i][j] / count_synd[i][j];
            double y = total_synd2[i][j] / count_synd[i][j];
            y -= x * x;

            if (i == 0)
                means[j] = x;

            if ((verbose_synd & 4) && (file != stderr))
                fprintf(file, "%-16.16g\t%-16.16g\t%d\t%d\t%d\n", x, sqrt(y), i, j, S_h[i][j]);
        }
    }

    if (log_file) {
        while (lockf(fileno(file), F_ULOCK, 0) < 0)
            ;
        fclose(file);
    }
}

static void print_usage(FILE *f, char *arg0) {
    fprintf(f,
            "usage: %s [OPTIONS]\n"
            "\n"
            "-D, --default          display simulation features and exit\n"
            "-r, --block_length     block length (smaller or equal to default)\n"
            "-t, --error_weight     error weight (else default value)\n"
            "-s, --seed             key seed\n"
            "-S, --error_seed       error seed (single decoding and exit)\n"
            "-N, --rounds           number of decoding to perform (default 10^6, infinite if < 0)\n"
            "-v, --verbose_synd     print level\n"
            "-o, --once             toggle run once/run forever (default is once)\n"
            "-l, --log              toggle output in stdout/file (default on stdout)\n"
            "-n, --sigma            sigma, noise (default 0)\n"
            "-m, --mode             reconstruction mode: partial, soft or both (default is both)\n",
            arg0);
    exit(2);
}

/*** parse command line arguments ***/
void parse_arguments_and_init_synd(int argc, char *argv[], long *seed, long *error_seed, long *nb_dec, bool *once, double *sigma, enum Mode *mode) {
    const char *options = "Dr:t:i:N:s:S:lov:n:m:";
    static struct option longopts[] = {{"default", no_argument, 0, 'D'},
                                       {"block_length", required_argument, 0, 'r'},
                                       {"error_weight", required_argument, 0, 't'},
                                       {"rounds", required_argument, 0, 'N'},
                                       {"seed", required_argument, 0, 's'},
                                       {"error_seed", required_argument, 0, 'S'},
                                       {"log", no_argument, 0, 'l'},
                                       {"once", no_argument, 0, 'o'},
                                       {"verbose_synd", required_argument, 0, 'v'},
                                       {"noise", required_argument, 0, 'n'},
                                       {"mode", required_argument, 0, 'm'},
                                       {NULL, 0, 0, 0}};
    /* default values */
    int blocklength = BLOCK_LENGTH;
    int error_weight = ERROR_WEIGHT;
    int weak_key_type = 0;
    int weak_param_1 = 0, weak_param_2 = 0;
    int error_param = -1;
    bool show_default = false;


    int ch;
    while ((ch = getopt_long(argc, argv, options, longopts, NULL)) != -1) {
        switch (ch) {
        case 'D':
            show_default = true;
            break;
            /* blocklength = BLOCK_LENGTH by default
               user may select blocklength < BLOCK_LENGTH
            */
        case 'r':
            blocklength = atoi(optarg);
            if (blocklength > BLOCK_LENGTH)
                print_usage(stderr, argv[0]);
            break;
            /* use an error weight different from ERROR_WEIGHT

               unlike -r, with the -t option the threshold schedule is
               modified to match with the actual error weight
            */
        case 't':
            error_weight = atoi(optarg);
            break;
        case 's':
            *seed = atol(optarg);
            break;
        case 'S':
            *error_seed = atol(optarg);
            break;

            /************************ begin simulation parameters */
            /* nb_dec decoding before exit or log */
        case 'N':
            *nb_dec = atol(optarg);
            break;
            /* never stop simulation */
        case 'o':
            *once = false;
            break;
            /************************ end simulation parameters */

            /* verbose_synd level, meaningful up to 7
               print extra info on stdout */
        case 'v':
            verbose_synd = atoi(optarg); /* global variable */
            break;
            /* log simulation data in a file */
        case 'l':
            log_file = !log_file; /* global variable, -l toggles default value */
            break;
            /* iter > allowed_iter is logged as a fail */
        case 'n':
            *sigma = atol(optarg);
            break;
        case 'm':
            if (strcmp(optarg, "partial") == 0) {*mode = PARTIAL;}
            else if (strcmp(optarg, "soft") == 0) {*mode = SOFT;}
            else if (strcmp(optarg, "both") == 0) {*mode = BOTH;}
            break;
           
        default:
            print_usage(stderr, argv[0]);
            break;
        }
    }

    if ((weak_key_type > 0) && (weak_param_1 <= 0))
        print_usage(stderr, argv[0]); // exit

    /* At the moment, special (i.e. non-uniform) keys and special
     * errors are generated up to isomorphism (cycle and power) so
     * using both may give incorrect outcome -> warning */
    if ((weak_key_type != NONE) && (error_param >= 0))
        fprintf(stderr, "warning: simultaneous use of special keys and errors\n");

    keygen_init(blocklength, weak_key_type, weak_param_1, weak_param_2);
    errgen_init(error_weight, error_param);

    if (show_default) {
    }

    if (*seed < 0)
        *seed = random_seed_synd();

    if (*error_seed < 0)
        *error_seed = random_seed_synd();

    /*** interuption handling ***/
    struct sigaction action_hup;
    action_hup.sa_handler = huphandler;
    sigemptyset(&action_hup.sa_mask);
    action_hup.sa_flags = 0;

    struct sigaction action_int;
    action_int.sa_handler = inthandler;
    sigemptyset(&action_int.sa_mask);
    action_int.sa_flags = 0;

    sigaction(SIGINT, &action_int, NULL);
    sigaction(SIGTERM, &action_int, NULL);
    sigaction(SIGHUP, &action_hup, NULL);
}
