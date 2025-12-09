#pragma once

#define INDEX 2
#define max_dist_h0 20

#ifdef NIST

#if NIST == 1
#define BLOCK_LENGTH 12323
#define BLOCK_WEIGHT 71
#define ERROR_WEIGHT 134
#define THRESHOLD_C1 0.006254868353074983
#define THRESHOLD_C0 11.101432337243956
#define BIKE_MAX_ITER 7
#define DELTA 3
#elif NIST == 3
#define BLOCK_LENGTH 24659
#define BLOCK_WEIGHT 103
#define ERROR_WEIGHT 199
#define THRESHOLD_C1 0.004533882596007288
#define THRESHOLD_C0 13.282669604666431
#define BIKE_MAX_ITER 7
#define DELTA 5
#elif NIST == 5
#define BLOCK_LENGTH 40973
#define BLOCK_WEIGHT 137
#define ERROR_WEIGHT 264
#define THRESHOLD_C1 0.0036083738659016262
#define THRESHOLD_C0 15.430866686308178
#define BIKE_MAX_ITER 7
#define DELTA 6
#endif

#endif // NIST

#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 12323
#define BLOCK_WEIGHT 71
#define ERROR_WEIGHT 134
#endif

#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0 // will force the threshold computation
#define THRESHOLD_C0 0
#endif

#ifndef BIKE_MAX_ITER
#define BIKE_MAX_ITER 7
#endif

#ifndef DELTA
#define DELTA 3
#endif
