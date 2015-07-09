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

/****h* Bluebottle/entrySearch
 * NAME
 *  entrySearch
 * FUNCTION
 *  An algorithm to find the minimum or maximum value of an array using a GPU.
 *  This program searches an array of reals of arbitrary length and returns the
 *  minimum or maximum value of the array.
 ******
 */

/****h* Bluebottle/entrySearch_kernel
 * NAME
 *  entrySearch_kernel
 * FUNCTION
 *  EntrySearch CUDA kernel functions.
 ******
 */

#ifndef __ENTRYSEARCH_H__
#define __ENTRYSEARCH_H__

extern "C"
{
#include "bluebottle.h"
}


/****f* entrySearch/find_min()
 * NAME
 *  find_min()
 * USAGE
 */
real find_min(int size, real *d_iarr);
/*
 * FUNCTION
 *  Return the minimum value of array d_iarr.  Note that d_iarr will be
 *  destroyed by the reduction algorithm.
 * ARGUMENTS
 *  * size -- the size of the array
 *  * d_iarr -- the array to search (must be allocated on the device); it will
 *    be destroyed by the reduction algorithm
 ******
 */

/****f* entrySearch/find_max()
 * NAME 
 *  find_max()
 * USAGE
 */
real find_max(int size, real *d_iarr);
/*
 * FUNCTION
 *  Return the maximum value of array d_iarr.  Note that d_iarr will be
 *  destroyed by the reduction algorithm.
 * ARGUMENTS
 *  * size -- the size of the array
 *  * d_iarr -- the array to search (must be allocated on the device); it will
 *    be destroyed by the reduction algorithm
 ******
 */

/****f* entrySearch/find_max_int()
 * NAME 
 *  find_max_int()
 * USAGE
 */
int find_max_int(int size, int *d_iarr);
/*
 * FUNCTION
 *  Return the maximum value of array d_iarr.  Note that d_iarr will be
 *  destroyed by the reduction algorithm.
 * ARGUMENTS
 *  * size -- the size of the array
 *  * d_iarr -- the array to search (must be allocated on the device); it will
 *    be destroyed by the reduction algorithm
 ******
 */

/****f* entrySearch/find_max_mag()
 * NAME 
 *  find_max_mag()
 * USAGE
 */
real find_max_mag(int size, real *d_iarr);
/*
 * FUNCTION
 *  Return the maximum magnitude value of array d_iarr.  Note that d_iarr will
 *  be destroyed by the reduction algorithm.
 * ARGUMENTS
 *  * size -- the size of the array
 *  * d_iarr -- the array to search (must be allocated on the device); it will
 *    be destroyed by the reduction algorithm
 ******
 */

/****f* entrySearch/avg_entries()
 * NAME 
 *  avg_entries()
 * USAGE
 */
real avg_entries(int size, real *d_iarr);
/*
 * FUNCTION
 *  Return the average of the entries of array d_iarr.  Note that d_iarr will
 *  be destroyed by the reduction algorithm.
 * ARGUMENTS
 *  * size -- the size of the array
 *  * d_iarr -- the array to average (must be allocated on the device); it will
 *    be destroyed by the reduction algorithm
 ******
 */

/****f* entrySearch/sum_entries()
 * NAME 
 *  sum_entries()
 * USAGE
 */
real sum_entries(int size, real *d_iarr);
/*
 * FUNCTION
 *  Return the sum of the entries of array d_iarr.  Note that d_iarr will
 *  be destroyed by the reduction algorithm.
 * ARGUMENTS
 *  * size -- the size of the array
 *  * d_iarr -- the array to average (must be allocated on the device); it will
 *    be destroyed by the reduction algorithm
 ******
 */

/****f* entrySearch_kernel/entrySearch_min_kernel<<<>>>()
 * NAME
 *  entrySearch_min_kernel<<<>>>()
 * USAGE
 */
__global__ void entrySearch_min_kernel(real *iarr, real *minarr,
  int size);
/*
 * FUNCTION
 *  Kernel function called by find_min().
 * ARGUMENTS
 *  * iarr -- the input array
 *  * minarr -- the ouput array
 *  * size -- the size of the input array
 ******
 */

/****f* entrySearch_kernel/entrySearch_max_kernel<<<>>>()
 * NAME
 *  entrySearch_max_kernel<<<>>>()
 * USAGE
 */
__global__ void entrySearch_max_kernel(real *iarr, real *maxarr,
  int size);
/*
 * FUNCTION
 *  Kernel function called by find_max().
 * ARGUMENTS
 *  * iarr -- the input array
 *  * maxarr -- the ouput array
 *  * size -- the size of the input array
 ******
 */

/****f* entrySearch_kernel/entrySearch_max_int_kernel<<<>>>()
 * NAME
 *  entrySearch_max_int_kernel<<<>>>()
 * USAGE
 */
__global__ void entrySearch_max_int_kernel(int *iarr, int *maxarr,
  int size);
/*
 * FUNCTION
 *  Kernel function called by find_max().
 * ARGUMENTS
 *  * iarr -- the input array
 *  * maxarr -- the ouput array
 *  * size -- the size of the input array
 ******
 */

/****f* entrySearch_kernel/entrySearch_max_mag_kernel<<<>>>()
 * NAME
 *  entrySearch_max_mag_kernel<<<>>>()
 * USAGE
 */
__global__ void entrySearch_max_mag_kernel(real *iarr, real *maxarr,
  int size);
/*
 * FUNCTION
 *  Kernel function called by find_max_mag().
 * ARGUMENTS
 *  * iarr -- the input array
 *  * maxarr -- the ouput array
 *  * size -- the size of the input array
 ******
 */

/****f* entrySearch_kernel/entrySearch_avg_entries_kernel<<<>>>()
 * NAME
 *  entrySearch_avg_entries_kernel<<<>>>()
 * USAGE
 */
__global__ void entrySearch_avg_entries_kernel(real *iarr, real *maxarr,
  int size);
/*
 * FUNCTION
 *  Kernel function called by avg_entries().
 * ARGUMENTS
 *  * iarr -- the input array
 *  * maxarr -- the ouput array
 *  * size -- the size of the input array
 ******
 */

#endif
