/*******************************************************************************
 ********************************* BLUEBOTTLE **********************************
 *******************************************************************************
 *
 *  Copyright 2012 - 2016 Adam Sierakowski, The Johns Hopkins University
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

#include <cuda.h>

#define MAXTHREADS 128
#define MAXBLOCKS 64

/* A bitwise function to determine the maximum exponent x that satisfies the
 * inequality 2^x < n.
 */
int floorLog2(unsigned int n) {
  int pos = 0;
  if (n >= 1<<16) { n >>= 16; pos += 16; }
  if (n >= 1<< 8) { n >>=  8; pos +=  8; }
  if (n >= 1<< 4) { n >>=  4; pos +=  4; }
  if (n >= 1<< 2) { n >>=  2; pos +=  2; }
  if (n >= 1<< 1) {           pos +=  1; }
  return ((n == 0) ? (-1) : pos);
}

/* A bitwise function to determine the minimum number n that satisfies the
 * inequality n > x, where n = 2^a for arbitrary a.
 */
unsigned int nextPow2(unsigned int x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

/* A function to determine the proper number of blocks and threads into which
 * the array should be split when parallelized on the GPU.
 */
void getNumBlocksAndThreads(int n, int &blocks, int &threads)
{
  threads = (n < MAXTHREADS * 2) ? nextPow2((n + 1) / 2): MAXTHREADS;
  blocks = (n + threads * 2 - 1) / (threads * 2);
}

/* A function to create random input data on CPU for program testing.  During
 * generation, the minimum value is recorded and returned for use in verifying
 * the GPU test result.
 */
void randArrGen(int size, real *arr, real* minmax)
{
  srand(time(NULL));
  for(int i=0; i<size; i++) {
    arr[i] = (rand() % size) - size / 2;
    if (arr[i] < minmax[0]) {
      minmax[0] = arr[i];
    }
    if (arr[i] > minmax[1]) {
      minmax[1] = arr[i];
    }
  }
}

/* The base function of the minimum search algorithm. */
real find_min(int size, real *d_iarr)
{
  int blocks = 0;
  int threads = 0;

  getNumBlocksAndThreads(size, blocks, threads);

  // create minarr on device
  int h_bytes = blocks * sizeof(real);
  real *d_minarr = NULL;
  (cudaMalloc((void**)&d_minarr, h_bytes));
  gpumem += h_bytes;

  cudaThreadSynchronize();

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(real);

  // run kernel
  entrySearch_min_kernel<<<dimGrid, dimBlock, smemSize>>>(d_iarr,
    d_minarr, size);

  cudaThreadSynchronize();

  // if there was more than one block, re-run the kernel on the minimum values 
  // from each of the blocks, which now reside in the first block_number indices
  // in d_minarr
  while(blocks > 1) {
    // use only the first block_number indices in min_arr
    size = blocks;
    getNumBlocksAndThreads(size, blocks, threads);

    entrySearch_min_kernel<<<dimGrid, dimBlock, smemSize>>>(d_minarr,
      d_minarr, size);

    cudaThreadSynchronize();
  }

  // grab final answer
  real min;
  (cudaMemcpy(&min, d_minarr, sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaFree(d_minarr));
  return min;
}

/* The base function of the maximum search algorithm. */
real find_max(int size, real *d_iarr)
{
  int blocks = 0;
  int threads = 0;

  getNumBlocksAndThreads(size, blocks, threads);

  // create minarr on device
  int h_bytes = blocks * sizeof(real);
  real *d_maxarr = NULL;
  (cudaMalloc((void**)&d_maxarr, h_bytes));
  gpumem += h_bytes;

  cudaThreadSynchronize();

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(real);

  // run kernel
  entrySearch_max_kernel<<<dimGrid, dimBlock, smemSize>>>(d_iarr,
    d_maxarr, size);

  cudaThreadSynchronize();

  // if there was more than one block, re-run the kernel on the maximum values 
  // from each of the blocks, which now reside in the first block_number indices
  // in d_minarr
  while(blocks > 1) {
    // use only the first block_number indices in min_arr
    size = blocks;
    getNumBlocksAndThreads(size, blocks, threads);

    entrySearch_max_kernel<<<dimGrid, dimBlock, smemSize>>>(d_maxarr,
      d_maxarr, size);

    cudaThreadSynchronize();
  }

  // grab final answer
  real max;
  (cudaMemcpy(&max, d_maxarr, sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaFree(d_maxarr));
  return max;
}

/* The base function of the maximum search algorithm. */
int find_max_int(int size, int *d_iarr)
{
  int blocks = 0;
  int threads = 0;

  getNumBlocksAndThreads(size, blocks, threads);

  // create minarr on device
  int h_bytes = blocks * sizeof(int);
  int *d_maxarr = NULL;
  (cudaMalloc((void**)&d_maxarr, h_bytes));
  gpumem += h_bytes;

  cudaThreadSynchronize();

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(int);

  // run kernel
  entrySearch_max_int_kernel<<<dimGrid, dimBlock, smemSize>>>(d_iarr,
    d_maxarr, size);

  cudaThreadSynchronize();

  // if there was more than one block, re-run the kernel on the maximum values 
  // from each of the blocks, which now reside in the first block_number indices
  // in d_minarr
  while(blocks > 1) {
    // use only the first block_number indices in min_arr
    size = blocks;
    getNumBlocksAndThreads(size, blocks, threads);

    entrySearch_max_int_kernel<<<dimGrid, dimBlock, smemSize>>>(d_maxarr,
      d_maxarr, size);

    cudaThreadSynchronize();
  }

  // grab final answer
  int max;
  (cudaMemcpy(&max, d_maxarr, sizeof(int),
    cudaMemcpyDeviceToHost));
  (cudaFree(d_maxarr));
  return max;
}

/* The base function of the maximum search algorithm. */
real find_max_mag(int size, real *d_iarr)
{
  int blocks = 0;
  int threads = 0;

  getNumBlocksAndThreads(size, blocks, threads);

  // create minarr on device
  int h_bytes = blocks * sizeof(real);
  real *d_maxarr = NULL;
  (cudaMalloc((void**)&d_maxarr, h_bytes));
  gpumem += h_bytes;

  cudaThreadSynchronize();

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(real);

  // run kernel
  entrySearch_max_mag_kernel<<<dimGrid, dimBlock, smemSize>>>(d_iarr,
    d_maxarr, size);

  cudaThreadSynchronize();

  // if there was more than one block, re-run the kernel on the maximum values 
  // from each of the blocks, which now reside in the first block_number indices
  // in d_minarr
  while(blocks > 1) {
    // use only the first block_number indices in min_arr
    size = blocks;
    getNumBlocksAndThreads(size, blocks, threads);

    entrySearch_max_mag_kernel<<<dimGrid, dimBlock, smemSize>>>(d_maxarr,
      d_maxarr, size);

    cudaThreadSynchronize();
  }

  // grab final answer
  real max;
  (cudaMemcpy(&max, d_maxarr, sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaFree(d_maxarr));
  return max;
}

/* The base function of the average algorithm. */
real avg_entries(int size, real *d_iarr)
{
  int blocks = 0;
  int threads = 0;

  int size_in = size;

  getNumBlocksAndThreads(size, blocks, threads);

  // create minarr on device
  int h_bytes = blocks * sizeof(real);
  real *d_maxarr = NULL;
  (cudaMalloc((void**)&d_maxarr, h_bytes));
  gpumem += h_bytes;

  cudaThreadSynchronize();

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(real);

  // run kernel
  entrySearch_avg_entries_kernel<<<dimGrid, dimBlock, smemSize>>>(d_iarr,
    d_maxarr, size);

  cudaThreadSynchronize();

  // if there was more than one block, re-run the kernel on the maximum values 
  // from each of the blocks, which now reside in the first block_number indices
  // in d_minarr
  while(blocks > 1) {
    // use only the first block_number indices in min_arr
    size = blocks;
    getNumBlocksAndThreads(size, blocks, threads);

    entrySearch_avg_entries_kernel<<<dimGrid, dimBlock, smemSize>>>(d_maxarr,
      d_maxarr, size);

    cudaThreadSynchronize();
  }

  // grab final answer
  real max;
  (cudaMemcpy(&max, d_maxarr, sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaFree(d_maxarr));
  return max / size_in;
}

/* The base function of the sum algorithm. */
real sum_entries(int size, real *d_iarr)
{
  int blocks = 0;
  int threads = 0;

  getNumBlocksAndThreads(size, blocks, threads);

  // create minarr on device
  int h_bytes = blocks * sizeof(real);
  real *d_maxarr = NULL;
  (cudaMalloc((void**)&d_maxarr, h_bytes));
  gpumem += h_bytes;

  cudaThreadSynchronize();

  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  int smemSize = threads * sizeof(real);

  // run kernel
  entrySearch_avg_entries_kernel<<<dimGrid, dimBlock, smemSize>>>(d_iarr,
    d_maxarr, size);

  cudaThreadSynchronize();

  // if there was more than one block, re-run the kernel on the maximum values 
  // from each of the blocks, which now reside in the first block_number indices
  // in d_minarr
  while(blocks > 1) {
    // use only the first block_number indices in min_arr
    size = blocks;
    getNumBlocksAndThreads(size, blocks, threads);

    entrySearch_avg_entries_kernel<<<dimGrid, dimBlock, smemSize>>>(d_maxarr,
      d_maxarr, size);

    cudaThreadSynchronize();
  }

  // grab final answer
  real max;
  (cudaMemcpy(&max, d_maxarr, sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaFree(d_maxarr));
  return max;
}

/* The main test function that creates a test array of random values and calls
 * find_min(...).  It displays both the known result as maintained through the
 * CPU-generated array and the GPU test result.
 */
/*int main(int argc, char** argv) 
{
  cudaDeviceProp deviceProp;
  deviceProp.major = 1;
  deviceProp.minor = 0;

  // force use of device number zero
  int dev = 0;

  (cudaSetDevice(dev));
  (cudaGetDeviceProperties(&deviceProp, dev));

  printf("\nUsing device %d: \"%s\"\n", dev, deviceProp.name);
  (cudaSetDevice(dev));

  // number of elements to reduce
  int size = pow(2, 23);

  printf("\nSearching %d randomly-generated elements", size);
  printf(" for the minimum value...\n");

  // create random input data on CPU
  real* h_arr = (real*) malloc(size * sizeof(real));
  cpumem += size * sizeof(real);
  real* minmax = (real*) malloc(2 * sizeof(real));
  cpumem += 2 * sizeof(real);
  randArrGen(size, h_arr, minmax);

  // load host data to device
  int numBlocks = 0;
  int numThreads = 0;
  getNumBlocksAndThreads(size, numBlocks, numThreads);
  unsigned int inbytes = size * sizeof(real);
  real* d_iarr = NULL;
  (cudaMalloc((void**)&d_iarr, inbytes));
  gpumem += in_bytes;
  (cudaMemcpy(d_iarr, h_arr, inbytes,
    cudaMemcpyHostToDevice));

  // run test
  real gpu_result = 0;
  int numcount = 100; // the number of iterations to test

  // timing stuff
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  real elapsedTime;

  // run GPU test
  printf("\nComputing GPU result %d times...\n", numcount);
  cudaEventRecord(start, 0);
  for(int count = 0; count < numcount; count++) {
    gpu_result = find_min(size, d_iarr);
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("...completed in %0.0f ms.\n", elapsedTime);

  // run CPU test
  printf("\nComputing CPU result %d times...\n", numcount);
  cudaEventRecord(start, 0);
  real cpu_result = size * 2;
  for(int count = 0; count < numcount; count++) {
    for(int z = 0; z < size; z++) {
      if(h_arr[z] < cpu_result) {
        cpu_result = h_arr[z];
      }
    }
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("...completed in %0.0f ms.\n", elapsedTime);

  // final minimum values
  printf("\nKnown result = %0.0f\n", minmax[0]);
  printf("CPU result = %0.0f\n", cpu_result);
  printf("GPU result = %0.0f\n", gpu_result);

  printf("\nSearching %d randomly-generated elements", size);
  printf(" for the maximum value...\n");

  // run GPU test
  printf("\nComputing GPU result %d times...\n", numcount);
  cudaEventRecord(start, 0);
  for(int count = 0; count < numcount; count++) {
    gpu_result = find_max(size, d_iarr);
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("...completed in %0.0f ms.\n", elapsedTime);

  // run CPU test
  printf("\nComputing CPU result %d times...\n", numcount);
  cudaEventRecord(start, 0);
  cpu_result = -size * 2;
  for(int count = 0; count < numcount; count++) {
    for(int z = 0; z < size; z++) {
      if(h_arr[z] > cpu_result) {
        cpu_result = h_arr[z];
      }
    }
  }
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("...completed in %0.0f ms.\n", elapsedTime);

  // final maximum values
  printf("\nKnown result = %0.0f\n", minmax[1]);
  printf("CPU result = %0.0f\n", cpu_result);
  printf("GPU result = %0.0f\n", gpu_result);

  // clean up
  (cudaFree(d_iarr));
  free(h_arr);
  free(minmax);

  cudaThreadExit();
  cutilExit(argc, argv);
}
*/
