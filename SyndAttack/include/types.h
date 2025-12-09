#pragma once
#include <stdbool.h>
#include <stdint.h>

#include "param.h"

/* Round relevant arrays size to the next multiple of 16 * 256 bits (to use the
 * 16 ymm AVX registers). */
#define ROUND_UP(N, S) ((((N) + (S)-1) / (S)) * (S))
#define SIZE_AVX (ROUND_UP(BLOCK_LENGTH * 8 * sizeof(bit_t), 256 * 16) / 8)

typedef int64_t index_t;
typedef index_t *sparse_t;

typedef uint8_t bit_t;
typedef bit_t *dense_t;

typedef struct decoder *decoder_t;

typedef struct {
    index_t blocklength; // effective blocklength <= BLOCK_LENGTH
    index_t columns[INDEX][BLOCK_WEIGHT];
    index_t rows[INDEX][BLOCK_WEIGHT];
} code_t;

typedef struct {
    bit_t vec[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
    index_t weight;
    index_t initial_weight; // for monitoring
} e_t;

typedef struct {
    bit_t vec[2 * SIZE_AVX] __attribute__((aligned(32)));
    index_t weight;
} syndrome_t;

typedef struct poly {
    int length, weight, alloc;
    uint16_t * index;
    uint64_t * extended;
} * poly_t;

typedef bit_t msg_t[2 * SIZE_AVX] __attribute__((aligned(32)));
typedef bit_t cw_t[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
typedef bit_t bits_t[INDEX][BLOCK_LENGTH];
typedef bit_t counters_t[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));

/* State of the decoder */
struct decoder {
    code_t *H;
    syndrome_t *syndrome;
    e_t *e;
    bits_t bits;
    counters_t counters;
    index_t iter, max_iter;
    bool blocked;
    double gap;
    double threshold;
    long key_seed, error_seed;
};
