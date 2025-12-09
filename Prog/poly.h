#include <stdint.h>
#include <stdbool.h>

typedef struct poly {
    int length, weight, alloc;
    uint16_t * index;
    uint64_t * extended;
} * poly_t;

poly_t poly_alloc(int length, int weight);
void poly_rand(poly_t h);
void poly_free(poly_t h);
int extended_weight(poly_t h);
bool extended_coeff_is_zero(poly_t h, int i);
poly_t extend_binomial_init(int j, int8_t *S, int length);
poly_t poly_shift_union(poly_t g, poly_t f, int l);
