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

// reference:http://c-faq.com/lib/gaussian.html
// generate double random gaussian distribution
double gaussrand(void)
{
        static double V1, V2, S;
        static int phase = 0;
        double X;

        if(phase == 0){
                do{
                        double U1 = rng_dbl();
                        double U2 = rng_dbl();

                        V1 = 2*U1 - 1;
                        V2 = 2*U2 - 1;
                        S = V1 * V1 + V2 * V2;
                        } while(S >=1 || S == 0);
                X = V1 * sqrt(-2 * log(S)/S);
        }else
                X = V2 * sqrt(-2 * log(S)/S);

        phase = 1 - phase;

        return X;
}
