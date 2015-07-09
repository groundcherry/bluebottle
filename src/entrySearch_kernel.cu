/*******************************************************************************
 ********************************* BLUEBOTTLE **********************************
 *******************************************************************************
 *
 *  Copyright 2012 - 2015 Adam Sierakowski, The Johns Hopkins University
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *  Please contact the Johns Hopkins University to use Bluebottle for
 *  commercial and/or for-profit applications.
 ******************************************************************************/

/* entrySearch
 *
 * An algorithm to find the minimum or maximum value of an array using a GPU.
 * This program searches an array of reals of arbitrary length and returns the
 * minimum or maximum value of the array.
 *
 * Adam J. Sierakowski, JHU/APL 2011
 */

#include "entrySearch.h"
#include <stdio.h>
#include <float.h>

/* The GPU kernel that performs the power-of-two minimum value search
 * algorithm.
 */
__global__ void entrySearch_min_kernel(real *g_iarr, real *g_minarr,
  int size)
{
    // create shared memory
    extern __shared__ real sarr[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
 
    if(i + blockDim.x < size) {
      if(g_iarr[i] < g_iarr[i + blockDim.x]) {
        sarr[tid] = g_iarr[i];
      } else {
        sarr[tid] = g_iarr[i + blockDim.x];
      }
    } else if (i < size) {
      sarr[tid] = g_iarr[i];
    } else {
#ifdef DOUBLE
        sarr[tid] = DBL_MAX;
#else
        sarr[tid] = FLT_MAX;
#endif
    }

    __syncthreads();

    // do comparison in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
      if(tid < s) {
        if(sarr[tid] > sarr[tid + s]) {
          sarr[tid] = sarr[tid + s];
        }
      }
      __syncthreads();
    }
  
    // write result for this block to global mem
    if(tid == 0) {
      g_minarr[blockIdx.x] = sarr[0];
    }
}

/* The GPU kernel that performs the power-of-two maximum value search
 * algorithm.
 */
__global__ void entrySearch_max_kernel(real *g_iarr, real *g_maxarr,
  int size)
{
    // create shared memory
    extern __shared__ real sarr[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
 
    if(i + blockDim.x < size) {
      if(g_iarr[i] > g_iarr[i + blockDim.x]) {
        sarr[tid] = g_iarr[i];
      } else {
        sarr[tid] = g_iarr[i + blockDim.x];
      }
    } else if (i < size) {
      sarr[tid] = g_iarr[i];
    } else {
#ifdef DOUBLE
        sarr[tid] = DBL_MIN;
#else
        sarr[tid] = FLT_MIN;
#endif
    }

    __syncthreads();

    // do comparison in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
      if(tid < s) {
        if(sarr[tid] < sarr[tid + s]) {
          sarr[tid] = sarr[tid + s];
        }
      }
      __syncthreads();
    }
  
    // write result for this block to global mem
    if(tid == 0) {
      g_maxarr[blockIdx.x] = sarr[0];
    }
}

/* The GPU kernel that performs the power-of-two maximum value search
 * algorithm.
 */
__global__ void entrySearch_max_int_kernel(int *g_iarr, int *g_maxarr,
  int size)
{
    // create shared memory
    extern __shared__ int sarr_int[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
 
    if(i + blockDim.x < size) {
      if(g_iarr[i] > g_iarr[i + blockDim.x]) {
        sarr_int[tid] = g_iarr[i];
      } else {
        sarr_int[tid] = g_iarr[i + blockDim.x];
      }
    } else if (i < size) {
      sarr_int[tid] = g_iarr[i];
    } else {
        sarr_int[tid] = INT_MIN;
    }

    __syncthreads();

    // do comparison in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
      if(tid < s) {
        if(sarr_int[tid] < sarr_int[tid + s]) {
          sarr_int[tid] = sarr_int[tid + s];
        }
      }
      __syncthreads();
    }
  
    // write result for this block to global mem
    if(tid == 0) {
      g_maxarr[blockIdx.x] = sarr_int[0];
    }
}

/* The GPU kernel that performs the power-of-two maximum magnitude value search
 * algorithm.
 */
__global__ void entrySearch_max_mag_kernel(real *g_iarr, real *g_maxarr,
  int size)
{
    // create shared memory
    extern __shared__ real sarr[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
 
    if(i + blockDim.x < size) {
      if(g_iarr[i] > fabs(g_iarr[i + blockDim.x])) {
        sarr[tid] = fabs(g_iarr[i]);
      } else {
        sarr[tid] = fabs(g_iarr[i + blockDim.x]);
      }
    } else if (i < size) {
      sarr[tid] = fabs(g_iarr[i]);
    } else {
#ifdef DOUBLE
        sarr[tid] = DBL_MIN;
#else
        sarr[tid] = FLT_MIN;
#endif
    }

    __syncthreads();

    // do comparison in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
      if(tid < s) {
        if(sarr[tid] < sarr[tid + s]) {
          sarr[tid] = sarr[tid + s];
        }
      }
      __syncthreads();
    }
  
    // write result for this block to global mem
    if(tid == 0) {
      g_maxarr[blockIdx.x] = sarr[0];
    }
}

/* The GPU kernel that performs the power-of-two average entries
 * algorithm.
 */
__global__ void entrySearch_avg_entries_kernel(real *g_iarr, real *g_maxarr,
  int size)
{
    // create shared memory
    extern __shared__ real sarr[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
 
    if(i + blockDim.x < size) {
      sarr[tid] = (g_iarr[i] + g_iarr[i + blockDim.x]);
    } else if (i < size) {
      sarr[tid] = g_iarr[i];
    } else {
      sarr[tid] = 0.;
    }

    __syncthreads();

    // do sum in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
      if(tid < s) {
        sarr[tid] = (sarr[tid] + sarr[tid + s]);
      }
      __syncthreads();
    }
  
    // write result for this block to global mem
    if(tid == 0) {
      g_maxarr[blockIdx.x] = sarr[0];
    }
}
