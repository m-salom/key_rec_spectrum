/*
  Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

  To the extent possible under law, the author has dedicated all copyright and
  related and neighboring rights to this software to the public domain
  worldwide. This software is distributed without any warranty.

  See <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#include <stdint.h>
#include <stdio.h>

#include "xoshiro256plusplus.h"

/*
  This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators. It
  has excellent (sub-ns) speed, a state (256 bits) that is large enough for any
  parallel application, and it passes all tests we are aware of.

  For generating just floating-point numbers, xoshiro256+ is even faster.

  The state must be seeded so that it is not everywhere zero. If you have a
  64-bit seed, we suggest to seed a splitmix64 generator and use its output to
  fill s.
*/

/* splitmix64 */
uint64_t splitmix64(uint64_t x) {
    uint64_t z = x + 0x9e3779b97f4a7c15;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

void seed_expand2(uint64_t *s, uint64_t x, uint64_t y) {
    s[0] = splitmix64(x);
    s[1] = splitmix64(s[0] ^ y);
    s[2] = splitmix64(s[1]);
    s[3] = splitmix64(s[2]);
}

void seed_expand(uint64_t *s, uint64_t x) {
    s[0] = splitmix64(x);
    s[1] = splitmix64(s[0]);
    s[2] = splitmix64(s[1]);
    s[3] = splitmix64(s[2]);
}

int seed_random(uint64_t *s) {
    FILE *urandom_fp;

    urandom_fp = fopen("/dev/urandom", "r");
    if (urandom_fp == NULL)
        return 0;
    if (fread(s, 8, 4, urandom_fp) != 4)
        return 0;
    fclose(urandom_fp);

    return 1;
}

void prng_from_seed(prng_t prng, uint64_t seed) {
    seed_expand(prng->s, seed); // initialize the 256 bits prng state from a 64 bits seed
    prng->random_lim = random_lim;
    prng->random_uint64_t = random_uint64_t;
}

void prng_from_seed2(prng_t prng, uint64_t seed1, uint64_t seed2) {
    seed_expand2(prng->s, seed1, seed2); // initialize the 256 bits prng state from two 64 bits seeds
    prng->random_lim = random_lim;
    prng->random_uint64_t = random_uint64_t;
}

void prng_init(prng_t prng) {
    seed_random(prng->s); // initialize the 256 bits prng state from /dev/urandom
    prng->random_lim = random_lim;
    prng->random_uint64_t = random_uint64_t;
}

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t random_uint64_t(uint64_t *s) {
    const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = rotl(s[3], 45);

    return result;
}

/* See
 * <https://lemire.me/blog/2019/06/06/nearly-divisionless-random-integer-generation-on-various-systems/>.
 */
uint64_t random_lim(uint64_t limit, uint64_t *s) {
    uint64_t x = random_uint64_t(s);
    __uint128_t m = (__uint128_t)x * (__uint128_t)limit;
    uint64_t l = (uint64_t)m;
    if (l < limit) {
        uint64_t t = -limit % limit;
        while (l < t) {
            x = random_uint64_t(s);
            m = (__uint128_t)x * (__uint128_t)limit;
            l = (uint64_t)m;
        }
    }
    return m >> 64;
}
