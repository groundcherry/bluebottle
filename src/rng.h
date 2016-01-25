/* Written by Adam Sierakowski 2015
 * for JHU HEART: The Computer Laboratory
 *
 * Compile with: 'make'
 */

#ifndef _RNG_H
#define _RNG_H

#include <stdlib.h>

extern unsigned long long int rng_j; // seed value
extern unsigned long long int rng_v; // random number sequence

// set seed value
void rng_init(unsigned long long int seed);

// random unsigned 64-bit integer
unsigned long long int rng_uint64(void);

// random unsigned 32-bit integer
unsigned int rng_uint32(void);

// random 64-bit integer
long long int rng_int64(void);

// random 32-bit integer
int rng_int32(void);

// random float between 0 and 1
float rng_flt(void);

// random double between 0 and 1
double rng_dbl(void);

#endif
