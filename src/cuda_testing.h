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

/****h* Bluebottle/cuda_testing_kernel
 * NAME
 *  cuda_testing_kernel
 * FUNCTION
 *  Bluebottle testing CUDA kernel functions.
 ******
 */

#ifndef _CUDA_TESTING_H
#define _CUDA_TESTING_H

extern "C"
{
#include "bluebottle.h"
}

/****f* cuda_testing_kernel/PP_memcpy_p_test<<<>>>()
 * NAME
 *  PP_memcpy_p_test<<<>>>()
 * USAGE
 */
__global__ void PP_memcpy_p_test(real *dst, real *src, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the computed _rhs_p to _p for testing output.
 * ARGUMENTS
 *  dst -- the destination _p
 *  src -- the source _rhs_p
 *  dom -- the dom_struct associated with this device
 ******
 */

#endif
