/* Written by Adam Sierakowski 2015
 * for JHU HEART: The Computer Laboratory
 *
 * Compile with: 'make'
 */

#include "rng.h"
/* Implementation of Numerical Recipes Third Edition Ranq1 (p. 351) */

unsigned long long int rng_j = 1;
unsigned long long int rng_v = 4101842887655102017ULL;

void rng_init(unsigned long long seed)
{
  // set seed
  rng_j = seed;
  // XOR assignment
  rng_v ^= rng_j;
  // run random number generator
  rng_v = rng_uint64();
}

// random unsigned 64-bit integer
unsigned long long int rng_uint64(void)
{
  // XOR assignments of bitshifted values
  rng_v ^= rng_v >> 21;
  rng_v ^= rng_v << 35;
  rng_v ^= rng_v >> 4;
  return rng_v * 2685821657736338717ULL;
}

// random unsigned 32-bit integer
unsigned int rng_uint32(void)
{
  return (unsigned int) rng_uint64();
}

// random 64-bit integer
long long int rng_int64(void)
{
  return (long long int) rng_uint64();
}

// random 32-bit integer
int rng_int32(void)
{
  return (int) rng_uint64();
}

// random float between 0 and 1
float rng_flt(void)
{
  return (float) rng_dbl();
}

// random double between 0 and 1
double rng_dbl(void)
{
  return 5.42101086242752217e-20 * rng_uint64();
}
