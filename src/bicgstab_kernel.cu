/*******************************************************************************
 ******************************* BLUEBOTTLE-1.0 ********************************
 *******************************************************************************
 *
 *  Copyright 2012 - 2014 Adam Sierakowski, The Johns Hopkins University
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

#include "cuda_bicgstab.h"

#include <cusp/dia_matrix.h>

__global__ void PP_rhs(real rho_f, real *u_star, real *v_star, real *w_star,
  real *rhs, dom_struct *dom, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
  __shared__ real s_rhs[MAX_THREADS_DIM * MAX_THREADS_DIM];     // solution

  // working constants
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  int C;

  // loop over z-planes
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_w0[ti + tj*blockDim.x] = w_star[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w_star[i + j*dom->Gfz._s1b
        + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u_star[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v_star[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

    s_rhs[ti + tj*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
      s_rhs[ti + tj*blockDim.x] = (s_u[(ti+1) + tj*blockDim.x]
        - s_u[ti + tj*blockDim.x]) * ddx;
      s_rhs[ti + tj*blockDim.x] += (s_v[ti + (tj+1)*blockDim.x]
        - s_v[ti + tj*blockDim.x]) * ddy;
      s_rhs[ti + tj*blockDim.x] += (s_w1[ti + tj*blockDim.x]
        - s_w0[ti + tj*blockDim.x]) * ddz;

      s_rhs[ti + tj*blockDim.x] *= rho_f / dt;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
      && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = (i-DOM_BUF) + (j-DOM_BUF)*dom->Gcc._s1 + (k-DOM_BUF)*dom->Gcc._s2;
      rhs[C] = s_rhs[ti + tj*blockDim.x];
    }
  }
}

__global__ void div_U(real *u, real *v, real *w,
  real *out, dom_struct *dom)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
  __shared__ real s_rhs[MAX_THREADS_DIM * MAX_THREADS_DIM];     // solution

  // working constants
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  int C;

  // loop over z-planes
  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b
        + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

    s_rhs[ti + tj*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
      s_rhs[ti + tj*blockDim.x] = (s_u[(ti+1) + tj*blockDim.x]
        - s_u[ti + tj*blockDim.x]) * ddx;
      s_rhs[ti + tj*blockDim.x] += (s_v[ti + (tj+1)*blockDim.x]
        - s_v[ti + tj*blockDim.x]) * ddy;
      s_rhs[ti + tj*blockDim.x] += (s_w1[ti + tj*blockDim.x]
        - s_w0[ti + tj*blockDim.x]) * ddz;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
      && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      out[C] = s_rhs[ti + tj*blockDim.x];
    }
  }
}

__global__ void coeffs_init(dom_struct *dom, int pitch, real *values)
{
  int i;  // iterator
  int C;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(tj < dom->Gcc.jn && tk < dom->Gcc.kn) {
    for(i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      values[C + 0*pitch]  = 0.;
      values[C + 1*pitch]  = 0.;
      values[C + 2*pitch]  = 0.;
      values[C + 3*pitch]  = 0.;
      values[C + 4*pitch]  = 0.;
      values[C + 5*pitch]  = 0.;
      values[C + 6*pitch]  = 0.;
      values[C + 7*pitch]  = 0.;
      values[C + 8*pitch]  = 0.;
      values[C + 9*pitch]  = 0.;
      values[C + 10*pitch] = 0.;
      values[C + 11*pitch] = 0.;
      values[C + 12*pitch] = 0.;
    }
  }
}

__global__ void coeffs(dom_struct *dom, int *flag_u, int *flag_v, int *flag_w,
  int pitch, real *values)
{
  int i;  // iterator
  int C, W, E, S, N, B, T;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  // loop over slices to set values
  if(tj < dom->Gcc.jn + DOM_BUF && tk < dom->Gcc.kn + DOM_BUF) {
    for(i = dom->Gcc.is; i < dom->Gcc.ie; i++) {
      W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      S = i + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      N = i + (tj+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      B = i + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
      T = i + tj*dom->Gfz.s1b + (tk+1)*dom->Gfz.s2b;
      C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
      values[C + pitch * 1]  += (real)abs(flag_w[B]) * ddz;
      values[C + pitch * 3]  += (real)abs(flag_v[S]) * ddy;
      values[C + pitch * 5]  += (real)abs(flag_u[W]) * ddx;
      values[C + pitch * 6]  -= (real)(abs(flag_u[W]) + abs(flag_u[E])) * ddx;
      values[C + pitch * 6]  -= (real)(abs(flag_v[S]) + abs(flag_v[N])) * ddy;
      values[C + pitch * 6]  -= (real)(abs(flag_w[B]) + abs(flag_w[T])) * ddz;
      values[C + pitch * 7]  += (real)abs(flag_u[E]) * ddx;
      values[C + pitch * 9]  += (real)abs(flag_v[N]) * ddy;
      values[C + pitch * 11] += (real)abs(flag_w[T]) * ddz;
    }
  }
}

__global__ void coeffs_periodic_W(dom_struct *dom, int pitch, real *values)
{
  int i = dom->Gcc.is;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 5] -= ddx;
    values[C + pitch * 8] += ddx;
  }
}

__global__ void coeffs_periodic_E(dom_struct *dom, int pitch, real *values)
{
  int i = dom->Gcc.ie-1;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 4] += ddx;
    values[C + pitch * 7] -= ddx;
  }
}

__global__ void coeffs_periodic_S(dom_struct *dom, int pitch, real *values)
{
  int j = dom->Gcc.js;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 3]  -= ddy;
    values[C + pitch * 10] += ddy;
  }
}

__global__ void coeffs_periodic_N(dom_struct *dom, int pitch, real *values)
{
  int j = dom->Gcc.je-1;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 2] += ddy;
    values[C + pitch * 9] -= ddy;
  }
}

__global__ void coeffs_periodic_B(dom_struct *dom, int pitch, real *values)
{
  int k = dom->Gcc.ks;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 1]  -= ddz;
    values[C + pitch * 12] += ddz;
  }
}

__global__ void coeffs_periodic_T(dom_struct *dom, int pitch, real *values)
{
  int k = dom->Gcc.ke-1;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 0]  += ddz;
    values[C + pitch * 11] -= ddz;
  }
}

__global__ void coeffs_particle(dom_struct *dom, int pitch, real *values,
  int *phase_shell)
{
  int i;  // iterator
  int C, CC;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  real denom = -2. / (dom->dx * dom->dx);
  denom -= 2./ (dom->dy * dom->dy);
  denom -= 2./ (dom->dz * dom->dz);

  if(tj < dom->Gcc.jn && tk < dom->Gcc.kn) {
    for(i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc.s1b + (tk+DOM_BUF)*dom->Gcc.s2b;
      if(phase_shell[CC] < 1) {
        values[C + 0*pitch] =  0.;
        values[C + 1*pitch] =  0.;
        values[C + 2*pitch] =  0.;
        values[C + 3*pitch] =  0.;
        values[C + 4*pitch] =  0.;
        values[C + 5*pitch] =  0.;
        values[C + 6*pitch] =  1.;
        values[C + 7*pitch] =  0.;
        values[C + 8*pitch] =  0.;
        values[C + 9*pitch] =  0.;
        values[C + 10*pitch] = 0.;
        values[C + 11*pitch] = 0.;
        values[C + 12*pitch] = 0.;
      }
    }
  }
}

__global__ void coeffs_zeros(dom_struct *dom, int pitch, real *values,
  real *rhs)
{
  int i;  // iterator
  int C;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  real denom = -2. / (dom->dx * dom->dx);
  denom -= 2./ (dom->dy * dom->dy);
  denom -= 2./ (dom->dz * dom->dz);

  // loop over all slices to find rows with all zeros
  if(tj < dom->Gcc.jn && tk < dom->Gcc.kn) {
    for(i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      real val = 0.;
      val += values[C + 0*pitch]*values[C + 0*pitch];
      val += values[C + 1*pitch]*values[C + 1*pitch];
      val += values[C + 2*pitch]*values[C + 2*pitch];
      val += values[C + 3*pitch]*values[C + 3*pitch];
      val += values[C + 4*pitch]*values[C + 4*pitch];
      val += values[C + 5*pitch]*values[C + 5*pitch];
      val += values[C + 6*pitch]*values[C + 6*pitch];
      val += values[C + 7*pitch]*values[C + 7*pitch];
      val += values[C + 8*pitch]*values[C + 8*pitch];
      val += values[C + 9*pitch]*values[C + 9*pitch];
      val += values[C + 10*pitch]*values[C + 10*pitch];
      val += values[C + 11*pitch]*values[C + 11*pitch];
      val += values[C + 12*pitch]*values[C + 12*pitch];
      // set rows with all zeros equal to one on the main diagonal
      if(val == 0.) values[C + 6*pitch] = 1.;//denom;
    }
  }
}

/*
// p; west
__global__ void PP_BC_p_W(real *A, int pitch, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if((tj < dom->Gcc._jn) && (tk < dom->Gcc._kn)) {
    int row = (dom->Gcc._is-DOM_BUF) + tj*dom->Gcc._s1 + tk*dom->Gcc._s2;
    A[row + 3 * pitch] += A[row + 2 * pitch];
    A[row + 2 * pitch] = 0.;
  }
}

// p; east
__global__ void PP_BC_p_E(real *A, int pitch, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if((tj < dom->Gcc._jn) && (tk < dom->Gcc._kn)) {
    int row = (dom->Gcc._ie-1-DOM_BUF) + tj*dom->Gcc._s1 + tk*dom->Gcc._s2;
    A[row + 3 * pitch] += A[row + 4 * pitch];
    A[row + 4 * pitch] = 0.;
  }
}

// p; south
__global__ void PP_BC_p_S(real *A, int pitch, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if((tk < dom->Gcc._kn) && (ti < dom->Gcc._in)) {
    int row = ti + (dom->Gcc._js-DOM_BUF)*dom->Gcc._s1 + tk*dom->Gcc._s2;
    A[row + 3 * pitch] += A[row + 1 * pitch];
    A[row + 1 * pitch] = 0.;
  }
}

// p; north
__global__ void PP_BC_p_N(real *A, int pitch, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if((tk < dom->Gcc._kn) && (ti < dom->Gcc._in)) {
    int row = ti + (dom->Gcc._je-1-DOM_BUF)*dom->Gcc._s1 + tk*dom->Gcc._s2;
    A[row + 3 * pitch] += A[row + 5 * pitch];
    A[row + 5 * pitch] = 0.;
  }
}

// p; bottom
__global__ void PP_BC_p_B(real *A, int pitch, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if((ti < dom->Gcc._in) && (tj < dom->Gcc._jn)) {
    int row = ti + tj*dom->Gcc._s1 + (dom->Gcc._ks-DOM_BUF)*dom->Gcc._s2;
    A[row + 3 * pitch] += A[row + 0 * pitch];
    A[row + 0 * pitch] = 0.;
  }
}

// p; top
__global__ void PP_BC_p_T(real *A, int pitch, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if((ti < dom->Gcc._in) && (tj < dom->Gcc._jn)) {
    int row = ti + tj*dom->Gcc._s1 + (dom->Gcc._ke-1-DOM_BUF)*dom->Gcc._s2;
    A[row + 3 * pitch] += A[row + 6 * pitch];
    A[row + 6 * pitch] = 0.;
  }
}
*/
