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

#include "cuda_bicgstab.h"

//#include <cusp/dia_matrix.h>

__global__ void ustar_rhs(real rho_f, real nu, real *u, real *v, real *w,
  real *p, real *f, real *conv0, real *conv, real *u_star, dom_struct *dom,
  real dt, real dt0)
{
  // create shared memory
  __shared__ real s_u0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u back
  __shared__ real s_u1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u center
  __shared__ real s_u2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
  __shared__ real s_w01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w back
  __shared__ real s_w12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w forward
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_u_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution

  // working constants
  real ab0 = 0.5 * dt / dt0;   // for Adams-Bashforth stepping
  real ab = 1. + ab0;          // for Adams-Bashforth stepping
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over u-planes
  for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int j = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int k = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tj = threadIdx.x;
    int tk = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    // TODO: THIS CAN BE FIXED BY PADDING ALL OF THESE ARRAYS WHEN COPYING FROM
    // HOST TO DEVICE
    if((k >= dom->Gfx._ksb && k < dom->Gfx._keb)
      && (j >= dom->Gfx._jsb && j < dom->Gfx._jeb)) {
      s_u0[tj + tk*blockDim.x] = u[(i-1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u1[tj + tk*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u2[tj + tk*blockDim.x] = u[(i+1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((k >= dom->Gfy._ksb && k < dom->Gfy._keb)
      && (j >= dom->Gfy._jsb && j < dom->Gfy._jeb)) {
      s_v01[tj + tk*blockDim.x] = v[(i-1) + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v12[tj + tk*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((k >= dom->Gfz._ksb && k < dom->Gfz._keb)
      && (j >= dom->Gfz._jsb && j < dom->Gfz._jeb)) {
      s_w01[tj + tk*blockDim.x] = w[(i-1) + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tj + tk*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }

    s_u_star[tj + tk*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((tj > 0 && tj < blockDim.x-1) && (tk > 0 && tk < blockDim.y-1)
      && j < dom->Gfx.jeb && k < dom->Gfx.keb) {
      // pressure gradient
      s_u_star[tj + tk*blockDim.x] =
        (p[(i-1) + j*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[i + j*dom->Gcc._s1b + k*dom->Gcc._s2b]) * ddx / rho_f;

      // grab the required data points for calculations
      real u011 = s_u0[tj + tk*blockDim.x];
      real u111 = s_u1[tj + tk*blockDim.x];
      real u211 = s_u2[tj + tk*blockDim.x];

      real u101 = s_u1[(tj-1) + tk*blockDim.x];
      real u121 = s_u1[(tj+1) + tk*blockDim.x];
      real v011 = s_v01[tj + tk*blockDim.x];
      real v111 = s_v12[tj + tk*blockDim.x];
      real v021 = s_v01[(tj+1) + tk*blockDim.x];
      real v121 = s_v12[(tj+1) + tk*blockDim.x];

      real u110 = s_u1[tj + (tk-1)*blockDim.x];
      real u112 = s_u1[tj + (tk+1)*blockDim.x];
      real w011 = s_w01[tj + tk*blockDim.x];
      real w111 = s_w12[tj + tk*blockDim.x];
      real w012 = s_w01[tj + (tk+1)*blockDim.x];
      real w112 = s_w12[tj + (tk+1)*blockDim.x];

      // compute convection term (Adams-Bashforth stepping)
      real duudx = (u211 + u111)*(u211 + u111) - (u111 + u011)*(u111 + u011);
      duudx *= 0.25 * ddx;

      real duvdy = (u121 + u111)*(v121 + v021) - (u111 + u101)*(v111 + v011);
      duvdy *= 0.25 * ddy;

      real duwdz = (u112 + u111)*(w112 + w012) - (u111 + u110)*(w111 + w011);
      duwdz *= 0.25 * ddz;

      s_c[tj + tk*blockDim.x] = duudx + duvdy + duwdz;

      // convection term sums into right-hand side
      if(dt0 > 0) // Adams-Bashforth
        s_u_star[tj + tk*blockDim.x] += (-ab * s_c[tj + tk*blockDim.x]
          + ab0 * conv0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]);
      else        // forward Euler
        s_u_star[tj + tk*blockDim.x] += -s_c[tj + tk*blockDim.x];

      // compute diffusion term (Adams-Bashforth stepping)
      real dud1 = (u211 - u111) * ddx;
      real dud0 = (u111 - u011) * ddx;
      real ddudxx = (dud1 - dud0) * ddx;

      dud1 = (u121 - u111) * ddy;
      dud0 = (u111 - u101) * ddy;
      real ddudyy = (dud1 - dud0) * ddy;

      dud1 = (u112 - u111) * ddz;
      dud0 = (u111 - u110) * ddz;
      real ddudzz = (dud1 - dud0) * ddz;

      s_d[tj + tk*blockDim.x] = nu * (ddudxx + ddudyy + ddudzz);

      // diffusive term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += 0.5*s_d[tj + tk*blockDim.x];

      // add on imposed pressure gradient
      s_u_star[tj + tk*blockDim.x] += f[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];

      // multiply by dt
      s_u_star[tj + tk*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += u111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((k >= dom->Gfx._ks && k < dom->Gfx._ke)
      && (j >= dom->Gfx._js && j < dom->Gfx._je)
      && (tj > 0 && tj < (blockDim.x-1))
      && (tk > 0 && tk < (blockDim.y-1))) {
      u_star[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]
        = s_u_star[tj + tk*blockDim.x];
      conv[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = s_c[tj + tk*blockDim.x];
    }
  }
}

__global__ void vstar_rhs(real rho_f, real nu, real *u, real *v, real *w,
  real *p, real *f, real *conv0, real *conv, real *v_star, dom_struct *dom,
  real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_v0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v back
  __shared__ real s_v1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v center
  __shared__ real s_v2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v forward
  __shared__ real s_w01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w back
  __shared__ real s_w12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_v_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution

  // working constants
  real ab0 = 0.5 * dt / dt0;   // for Adams-Bashforth stepping
  real ab = 1. + ab0;          // for Adams-Bashforth stepping
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over v-planes
  for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int k = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int i = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tk = threadIdx.x;
    int ti = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    if((i >= dom->Gfy._isb && i < dom->Gfy._ieb)
      && (k >= dom->Gfy._ksb && k < dom->Gfy._keb)) {
      s_v0[tk + ti*blockDim.x] = v[i + (j-1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v1[tk + ti*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v2[tk + ti*blockDim.x] = v[i + (j+1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((i >= dom->Gfz._isb && i < dom->Gfz._ieb)
      && (k >= dom->Gfz._ksb && k < dom->Gfz._keb)) {
      s_w01[tk + ti*blockDim.x] = w[i + (j-1)*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tk + ti*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
    if((i >= dom->Gfx._isb && i < dom->Gfx._ieb)
      && (k >= dom->Gfx._ksb && k < dom->Gfx._keb)) {
      s_u01[tk + ti*blockDim.x] = u[i + (j-1)*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u12[tk + ti*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }

    s_v_star[tk + ti*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((tk > 0 && tk < blockDim.x-1) && (ti > 0 && ti < blockDim.y-1)
      && k < dom->Gfy.keb && i < dom->Gfy.ieb) {
      // pressure gradient
      s_v_star[tk + ti*blockDim.x] =
        (p[i + (j-1)*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[i + j*dom->Gcc._s1b + k*dom->Gcc._s2b]) * ddy / rho_f;

      // grab the required data points for calculations
      real v101 = s_v0[tk + ti*blockDim.x];
      real v111 = s_v1[tk + ti*blockDim.x];
      real v121 = s_v2[tk + ti*blockDim.x];

      real v110 = s_v1[(tk-1) + ti*blockDim.x];
      real v112 = s_v1[(tk+1) + ti*blockDim.x];
      real w101 = s_w01[tk + ti*blockDim.x];
      real w111 = s_w12[tk + ti*blockDim.x];
      real w102 = s_w01[(tk+1) + ti*blockDim.x];
      real w112 = s_w12[(tk+1) + ti*blockDim.x];

      real v011 = s_v1[tk + (ti-1)*blockDim.x];
      real v211 = s_v1[tk + (ti+1)*blockDim.x];
      real u101 = s_u01[tk + ti*blockDim.x];
      real u111 = s_u12[tk + ti*blockDim.x];
      real u201 = s_u01[tk + (ti+1)*blockDim.x];
      real u211 = s_u12[tk + (ti+1)*blockDim.x];

      // compute convection term (Adams-Bashforth stepping)
      real dvudx = (v211 + v111)*(u211 + u201) - (v111 + v011)*(u111 + u101);
      dvudx *= 0.25 * ddx;

      real dvvdy = (v121 + v111)*(v121 + v111) - (v111 + v101)*(v111 + v101);
      dvvdy *= 0.25 * ddy;

      real dvwdz = (v112 + v111)*(w112 + w102) - (v111 + v110)*(w111 + w101);
      dvwdz *= 0.25 * ddz;

      s_c[tk + ti*blockDim.x] = dvudx + dvvdy + dvwdz;

      // convection term sums into right-hand side
      if(dt0 > 0) // Adams-Bashforth
        s_v_star[tk + ti*blockDim.x] += (-ab * s_c[tk + ti*blockDim.x]
          + ab0 * conv0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b]);
      else        // forward Euler
        s_v_star[tk + ti*blockDim.x] += -s_c[tk + ti*blockDim.x];

      // compute diffusive term
      real dvd1 = (v211 - v111) * ddx;
      real dvd0 = (v111 - v011) * ddx;
      real ddvdxx = (dvd1 - dvd0) * ddx;

      dvd1 = (v121 - v111) * ddy;
      dvd0 = (v111 - v101) * ddy;
      real ddvdyy = (dvd1 - dvd0) * ddy;

      dvd1 = (v112 - v111) * ddz;
      dvd0 = (v111 - v110) * ddz;
      real ddvdzz = (dvd1 - dvd0) * ddz;

      s_d[tk + ti*blockDim.x] = nu * (ddvdxx + ddvdyy + ddvdzz);

      // diffusive term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += 0.5*s_d[tk + ti*blockDim.x];

      // add on imposed pressure gradient
      s_v_star[tk + ti*blockDim.x] += f[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];

      // multiply by dt
      s_v_star[tk + ti*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += v111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((i >= dom->Gfy._is && i < dom->Gfy._ie)
      && (k >= dom->Gfy._ks && k < dom->Gfy._ke)
      && (tk > 0 && tk < (blockDim.x-1))
      && (ti > 0 && ti < (blockDim.y-1))) {
      v_star[i+ j*dom->Gfy._s1b + k*dom->Gfy._s2b]
        = s_v_star[tk + ti*blockDim.x];
      conv[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = s_c[tk + ti*blockDim.x];
    }
  }
}

__global__ void wstar_rhs(real rho_f, real nu, real *u, real *v, real *w,
  real *p, real *f, real *conv0, real *conv, real *w_star, dom_struct *dom,
  real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_w_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution

  // working constants
  real ab0 = 0.5 * dt / dt0;   // for Adams-Bashforth stepping
  real ab = 1. + ab0;          // for Adams-Bashforth stepping
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over w-planes
  for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
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
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u01[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + (k-1)*dom->Gfx._s2b];
      s_u12[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v01[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + (k-1)*dom->Gfy._s2b];
      s_v12[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

    s_w_star[ti + tj*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)
      && i < dom->Gfz.ieb && j < dom->Gfz.jeb) {
      // pressure gradient
      s_w_star[ti + tj*blockDim.x] =
        (p[i + j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]
        - p[i + j*dom->Gcc._s1b + k*dom->Gcc._s2b]) * ddz / rho_f;

      // grab the required data points for calculations
      real w110 = s_w0[ti + tj*blockDim.x];
      real w111 = s_w1[ti + tj*blockDim.x];
      real w112 = s_w2[ti + tj*blockDim.x];

      real w011 = s_w1[(ti-1) + tj*blockDim.x];
      real w211 = s_w1[(ti+1) + tj*blockDim.x];
      real u110 = s_u01[ti + tj*blockDim.x];
      real u111 = s_u12[ti + tj*blockDim.x];
      real u210 = s_u01[(ti+1) + tj*blockDim.x];
      real u211 = s_u12[(ti+1) + tj*blockDim.x];

      real w101 = s_w1[ti + (tj-1)*blockDim.x];
      real w121 = s_w1[ti + (tj+1)*blockDim.x];
      real v110 = s_v01[ti + tj*blockDim.x];
      real v111 = s_v12[ti + tj*blockDim.x];
      real v120 = s_v01[ti + (tj+1)*blockDim.x];
      real v121 = s_v12[ti + (tj+1)*blockDim.x];

      // compute convection term (Adams-Bashforth stepping)
      real dwudx = (w211 + w111)*(u211 + u210) - (w111 + w011)*(u111 + u110);
      dwudx *= 0.25 * ddx;

      real dwvdy = (w121 + w111)*(v121 + v120) - (w111 + w101)*(v111 + v110);
      dwvdy *= 0.25 * ddy;

      real dwwdz = (w112 + w111)*(w112 + w111) - (w111 + w110)*(w111 + w110);
      dwwdz *= 0.25 * ddz;

      s_c[ti + tj*blockDim.x] = dwudx + dwvdy + dwwdz;

      // convection term sums into right-hand side
      if(dt0 > 0) // Adams-Bashforth
        s_w_star[ti + tj*blockDim.x] += (-ab * s_c[ti + tj*blockDim.x]
          + ab0 * conv0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]);
      else        // forward Euler
        s_w_star[ti + tj*blockDim.x] += -s_c[ti + tj*blockDim.x];

      // compute diffusive term
      real dwd1 = (w211 - w111) * ddx;
      real dwd0 = (w111 - w011) * ddx;
      real ddwdxx = (dwd1 - dwd0) * ddx;

      dwd1 = (w121 - w111) * ddy;
      dwd0 = (w111 - w101) * ddy;
      real ddwdyy = (dwd1 - dwd0) * ddy;

      dwd1 = (w112 - w111) * ddz;
      dwd0 = (w111 - w110) * ddz;
      real ddwdzz = (dwd1 - dwd0) * ddz;

      s_d[ti + tj*blockDim.x] = nu * (ddwdxx + ddwdyy + ddwdzz);

      // diffusive term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += 0.5*s_d[ti + tj*blockDim.x];

      // add on imposed pressure gradient
      s_w_star[ti + tj*blockDim.x] += f[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];

      // multiply by dt
      s_w_star[ti + tj*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += w111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      w_star[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]
        = s_w_star[ti + tj*blockDim.x];
      conv[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_c[ti + tj*blockDim.x];
    }
  }
}

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

__global__ void ustar_coeffs_init(dom_struct *dom, int pitch, real *values)
{
  int i;  // iterator
  int C;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    for(i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gfx.s1 + tk*dom->Gfx.s2;
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

__global__ void vstar_coeffs_init(dom_struct *dom, int pitch, real *values)
{
  int j;  // iterator
  int C;  // cell

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    for(j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      C = ti + j*dom->Gfy.s1 + tk*dom->Gfy.s2;
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

__global__ void wstar_coeffs_init(dom_struct *dom, int pitch, real *values)
{
  int k;  // iterator
  int C;  // cell

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    for(k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      C = ti + tj*dom->Gfz.s1 + k*dom->Gfz.s2;
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

__global__ void ustar_coeffs_particles(dom_struct *dom, int pitch, real *values,
  int *flag_u)
{
  int i;  // iterator
  int C, CC;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    for(i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gfx.s1 + tk*dom->Gfx.s2;
      CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
      if(flag_u[CC] < 0) {
        values[C + 0*pitch]  = 0.;
        values[C + 1*pitch]  = 0.;
        values[C + 2*pitch]  = 0.;
        values[C + 3*pitch]  = 0.;
        values[C + 4*pitch]  = 0.;
        values[C + 5*pitch]  = 0.;
        values[C + 6*pitch]  = 1.;
        values[C + 7*pitch]  = 0.;
        values[C + 8*pitch]  = 0.;
        values[C + 9*pitch]  = 0.;
        values[C + 10*pitch] = 0.;
        values[C + 11*pitch] = 0.;
        values[C + 12*pitch] = 0.;
      }
    }
  }
}

__global__ void vstar_coeffs_particles(dom_struct *dom, int pitch, real *values,
  int *flag_v)
{
  int j;  // iterator
  int C, CC;  // cell

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    for(j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      C = ti + j*dom->Gfy.s1 + tk*dom->Gfy.s2;
      CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
      if(flag_v[CC] < 0) {
        values[C + 0*pitch]  = 0.;
        values[C + 1*pitch]  = 0.;
        values[C + 2*pitch]  = 0.;
        values[C + 3*pitch]  = 0.;
        values[C + 4*pitch]  = 0.;
        values[C + 5*pitch]  = 0.;
        values[C + 6*pitch]  = 1.;
        values[C + 7*pitch]  = 0.;
        values[C + 8*pitch]  = 0.;
        values[C + 9*pitch]  = 0.;
        values[C + 10*pitch] = 0.;
        values[C + 11*pitch] = 0.;
        values[C + 12*pitch] = 0.;
      }
    }
  }
}

__global__ void wstar_coeffs_particles(dom_struct *dom, int pitch, real *values,
  int *flag_w)
{
  int k;  // iterator
  int C, CC;  // cell

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    for(k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      C = ti + tj*dom->Gfz.s1 + k*dom->Gfz.s2;
      CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (k+DOM_BUF)*dom->Gfz.s2b;
      if(flag_w[CC] < 0) {
        values[C + 0*pitch]  = 0.;
        values[C + 1*pitch]  = 0.;
        values[C + 2*pitch]  = 0.;
        values[C + 3*pitch]  = 0.;
        values[C + 4*pitch]  = 0.;
        values[C + 5*pitch]  = 0.;
        values[C + 6*pitch]  = 1.;
        values[C + 7*pitch]  = 0.;
        values[C + 8*pitch]  = 0.;
        values[C + 9*pitch]  = 0.;
        values[C + 10*pitch] = 0.;
        values[C + 11*pitch] = 0.;
        values[C + 12*pitch] = 0.;
      }
    }
  }
}

__global__ void ustar_coeffs(real nu, real dt, dom_struct *dom, int pitch,
  real *values, int *flag_u, int *flag_v, int *flag_w)
{
  int i;  // iterator
  int C, W, E, S, N, B, T;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  // loop over slices to set values
  if(tj < dom->Gfx.jn + DOM_BUF && tk < dom->Gfx.kn + DOM_BUF) {
    for(i = dom->Gfx.is; i < dom->Gfx.ie; i++) {
      W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      S = i + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      N = i + (tj+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      B = i + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
      T = i + tj*dom->Gfz.s1b + (tk+1)*dom->Gfz.s2b;
      C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gfx.s1 + (tk-DOM_BUF)*dom->Gfx.s2;
      values[C + pitch * 1]  -= (real)abs(flag_w[B])*0.5*nu*dt*ddz;
      values[C + pitch * 3]  -= (real)abs(flag_v[S])*0.5*nu*dt*ddy;
      values[C + pitch * 5]  -= (real)abs(flag_u[W])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += 1.;
      values[C + pitch * 6]  += (real)abs(flag_u[W])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += (real)abs(flag_u[E])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += (real)abs(flag_v[S])*0.5*nu*dt*ddy;
      values[C + pitch * 6]  += (real)abs(flag_v[N])*0.5*nu*dt*ddy;
      values[C + pitch * 6]  += (real)abs(flag_w[B])*0.5*nu*dt*ddz;
      values[C + pitch * 6]  += (real)abs(flag_w[T])*0.5*nu*dt*ddz;
      values[C + pitch * 7]  -= (real)abs(flag_u[E])*0.5*nu*dt*ddx;
      values[C + pitch * 9]  -= (real)abs(flag_v[N])*0.5*nu*dt*ddy;
      values[C + pitch * 11] -= (real)abs(flag_w[T])*0.5*nu*dt*ddz;
    }
  }
}

__global__ void vstar_coeffs(real nu, real dt, dom_struct *dom, int pitch,
  real *values, int *flag_u, int *flag_v, int *flag_w)
{
  int j;  // iterator
  int C, W, E, S, N, B, T;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  // loop over slices to set values
  if(tk < dom->Gfy.kn + DOM_BUF && ti < dom->Gfy.in + DOM_BUF) {
    for(j = dom->Gfy.js; j < dom->Gfy.je; j++) {
      W = ti + j*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      E = (ti+1) + j*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      S = ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      N = ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      B = ti + j*dom->Gfz.s1b + tk*dom->Gfz.s2b;
      T = ti + j*dom->Gfz.s1b + (tk+1)*dom->Gfz.s2b;
      C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gfy.s1 + (tk-DOM_BUF)*dom->Gfy.s2;
      values[C + pitch * 1]  -= (real)abs(flag_w[B])*0.5*nu*dt*ddz;
      values[C + pitch * 3]  -= (real)abs(flag_v[S])*0.5*nu*dt*ddy;
      values[C + pitch * 5]  -= (real)abs(flag_u[W])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += 1.;
      values[C + pitch * 6]  += (real)abs(flag_u[W])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += (real)abs(flag_u[E])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += (real)abs(flag_v[S])*0.5*nu*dt*ddy;
      values[C + pitch * 6]  += (real)abs(flag_v[N])*0.5*nu*dt*ddy;
      values[C + pitch * 6]  += (real)abs(flag_w[B])*0.5*nu*dt*ddz;
      values[C + pitch * 6]  += (real)abs(flag_w[T])*0.5*nu*dt*ddz;
      values[C + pitch * 7]  -= (real)abs(flag_u[E])*0.5*nu*dt*ddx;
      values[C + pitch * 9]  -= (real)abs(flag_v[N])*0.5*nu*dt*ddy;
      values[C + pitch * 11] -= (real)abs(flag_w[T])*0.5*nu*dt*ddz;
    }
  }
}

__global__ void wstar_coeffs(real nu, real dt, dom_struct *dom, int pitch,
  real *values, int *flag_u, int *flag_v, int *flag_w)
{
  int k;  // iterator
  int C, W, E, S, N, B, T;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  // loop over slices to set values
  if(ti < dom->Gfz.in + DOM_BUF && tj < dom->Gfz.jn + DOM_BUF) {
    for(k = dom->Gfz.ks; k < dom->Gfz.ke; k++) {
      W = ti + tj*dom->Gfx.s1b + k*dom->Gfx.s2b;
      E = (ti+1) + tj*dom->Gfx.s1b + k*dom->Gfx.s2b;
      S = ti + tj*dom->Gfy.s1b + k*dom->Gfy.s2b;
      N = ti + (tj+1)*dom->Gfy.s1b + k*dom->Gfy.s2b;
      B = ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b;
      T = ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b;
      C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gfz.s1 + (k-DOM_BUF)*dom->Gfz.s2;
      values[C + pitch * 1]  -= (real)abs(flag_w[B])*0.5*nu*dt*ddz;
      values[C + pitch * 3]  -= (real)abs(flag_v[S])*0.5*nu*dt*ddy;
      values[C + pitch * 5]  -= (real)abs(flag_u[W])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += 1.;
      values[C + pitch * 6]  += (real)abs(flag_u[W])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += (real)abs(flag_u[E])*0.5*nu*dt*ddx;
      values[C + pitch * 6]  += (real)abs(flag_v[S])*0.5*nu*dt*ddy;
      values[C + pitch * 6]  += (real)abs(flag_v[N])*0.5*nu*dt*ddy;
      values[C + pitch * 6]  += (real)abs(flag_w[B])*0.5*nu*dt*ddz;
      values[C + pitch * 6]  += (real)abs(flag_w[T])*0.5*nu*dt*ddz;
      values[C + pitch * 7]  -= (real)abs(flag_u[E])*0.5*nu*dt*ddx;
      values[C + pitch * 9]  -= (real)abs(flag_v[N])*0.5*nu*dt*ddy;
      values[C + pitch * 11] -= (real)abs(flag_w[T])*0.5*nu*dt*ddz;
    }
  }
}

__global__ void ustar_coeffs_dirichlet_W(real bc, real *u_star, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfx.is-DOM_BUF;  // iterator
  int C, CC;  // cell location

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;
  
  if((tj < dom->Gfx.jn) && (tk < dom->Gfx.kn)) {
    C = i + tj*dom->Gfx.s1 + tk*dom->Gfx.s2;
    CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
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
    u_star[CC] = bc;
  }
}

__global__ void ustar_coeffs_dirichlet_E(real bc, real *u_star, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfx.ie-1-DOM_BUF;  // iterator
  int C, CC;  // cell location

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;
  
  if((tj < dom->Gfx.jn) && (tk < dom->Gfx.kn)) {
    C = i + tj*dom->Gfx.s1 + tk*dom->Gfx.s2;
    CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
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
    u_star[CC] = bc;
  }
}

__global__ void ustar_coeffs_dirichlet_S(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfx.js-DOM_BUF;  // iterator
  int C, CC, S;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfx.kn) && (ti < dom->Gfx.in)) {
    C = ti + j*dom->Gfx.s1 + tk*dom->Gfx.s2;
    CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
    S = (ti+DOM_BUF) + (j-1+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
    values[C + 3*pitch] =  0.;
    u_star[CC] += u_star[S] * 0.5 * nu * ddy;
  }
}

__global__ void ustar_coeffs_dirichlet_N(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfx.je-1-DOM_BUF;  // iterator
  int C, CC, N;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfx.kn) && (ti < dom->Gfx.in)) {
    C = ti + j*dom->Gfx.s1 + tk*dom->Gfx.s2;
    CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
    N = (ti+DOM_BUF) + (j+1+DOM_BUF)*dom->Gfx.s1b + (tk+DOM_BUF)*dom->Gfx.s2b;
    values[C + 9*pitch] =  0.;
    u_star[CC] += u_star[N] * 0.5 * nu * ddy;
  }
}

__global__ void ustar_coeffs_dirichlet_B(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfx.ks-DOM_BUF;  // iterator
  int C, CC, B;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfx.in) && (tj < dom->Gfx.jn)) {
    C = ti + tj*dom->Gfx.s1 + k*dom->Gfx.s2;
    CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (k+DOM_BUF)*dom->Gfx.s2b;
    B = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (k-1+DOM_BUF)*dom->Gfx.s2b;
    values[C + 1*pitch] =  0.;
    u_star[CC] += u_star[B] * 0.5 * nu * ddz;
  }
}

__global__ void ustar_coeffs_dirichlet_T(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfx.ke-1-DOM_BUF;  // iterator
  int C, CC, T;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfx.in) && (tj < dom->Gfx.jn)) {
    C = ti + tj*dom->Gfx.s1 + k*dom->Gfx.s2;
    CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (k+DOM_BUF)*dom->Gfx.s2b;
    T = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx.s1b + (k+1+DOM_BUF)*dom->Gfx.s2b;
    values[C + 11*pitch] =  0.;
    u_star[CC] += u_star[T] * 0.5 * nu * ddz;
  }
}

__global__ void vstar_coeffs_dirichlet_W(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfy.is-DOM_BUF;  // iterator
  int C, CC, W;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfy.jn) && (tk < dom->Gfy.kn)) {
    C = i + tj*dom->Gfy.s1 + tk*dom->Gfy.s2;
    CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
    W = (i-1+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
    values[C + 5*pitch] =  0.;
    v_star[CC] += v_star[W] * 0.5 * nu * ddx;
  }
}

__global__ void vstar_coeffs_dirichlet_E(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfy.ie-1-DOM_BUF;  // iterator
  int C, CC, E;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfy.jn) && (tk < dom->Gfy.kn)) {
    C = i + tj*dom->Gfy.s1 + tk*dom->Gfy.s2;
    CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
    E = (i+1+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
    values[C + 7*pitch] =  0.;
    v_star[CC] += v_star[E] * 0.5 * nu * ddx;
  }
}

__global__ void vstar_coeffs_dirichlet_S(real bc, real *v_star, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfy.js-DOM_BUF;  // iterator
  int C, CC;  // cell locations

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfy.kn) && (ti < dom->Gfy.in)) {
    C = ti + j*dom->Gfy.s1 + tk*dom->Gfy.s2;
    CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
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
    v_star[CC] = bc;
  }
}

__global__ void vstar_coeffs_dirichlet_N(real bc, real *v_star, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfy.je-1-DOM_BUF;  // iterator
  int C, CC;  // cell locations

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfy.kn) && (ti < dom->Gfy.in)) {
    C = ti + j*dom->Gfy.s1 + tk*dom->Gfy.s2;
    CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy.s1b + (tk+DOM_BUF)*dom->Gfy.s2b;
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
    v_star[CC] = bc;
  }
}

__global__ void vstar_coeffs_dirichlet_B(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfy.ks-DOM_BUF;  // iterator
  int C, CC, B;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfy.in) && (tj < dom->Gfy.jn)) {
    C = ti + tj*dom->Gfy.s1 + k*dom->Gfy.s2;
    CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (k+DOM_BUF)*dom->Gfy.s2b;
    B = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (k-1+DOM_BUF)*dom->Gfy.s2b;
    values[C + 1*pitch] =  0.;
    v_star[CC] += v_star[B] * 0.5 * nu * ddz;
  }
}

__global__ void vstar_coeffs_dirichlet_T(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfy.ke-1-DOM_BUF;  // iterator
  int C, CC, T;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfy.in) && (tj < dom->Gfy.jn)) {
    C = ti + tj*dom->Gfy.s1 + k*dom->Gfy.s2;
    CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (k+DOM_BUF)*dom->Gfy.s2b;
    T = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfy.s1b + (k+1+DOM_BUF)*dom->Gfy.s2b;
    values[C + 11*pitch] = 0.;
    v_star[CC] += v_star[T] * 0.5 * nu * ddz;
  }
}

__global__ void wstar_coeffs_dirichlet_W(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfz.is-DOM_BUF;  // iterator
  int C, CC, W;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfz.jn) && (tk < dom->Gfz.kn)) {
    C = i + tj*dom->Gfz.s1 + tk*dom->Gfz.s2;
    CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    W = (i-1+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    values[C + 5*pitch] =  0.;
    w_star[CC] += w_star[W] * 0.5 * nu * ddx;
  }
}

__global__ void wstar_coeffs_dirichlet_E(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfz.ie-1-DOM_BUF;  // iterator
  int C, CC, E;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfz.jn) && (tk < dom->Gfz.kn)) {
    C = i + tj*dom->Gfz.s1 + tk*dom->Gfz.s2;
    CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    E = (i+1+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    values[C + 7*pitch] =  0.;
    w_star[CC] += w_star[E] * 0.5 * nu * ddx;
  }
}

__global__ void wstar_coeffs_dirichlet_S(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfz.js-DOM_BUF;  // iterator
  int C, CC, S;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfz.kn) && (ti < dom->Gfz.in)) {
    C = ti + j*dom->Gfz.s1 + tk*dom->Gfz.s2;
    CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    S = (ti+DOM_BUF) + (j-1+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    values[C + 3*pitch] =  0.;
    w_star[CC] += w_star[S] * 0.5 * nu * ddy;
  }
}

__global__ void wstar_coeffs_dirichlet_N(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfz.je-1-DOM_BUF;  // iterator
  int C, CC, N;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfz.kn) && (ti < dom->Gfz.in)) {
    C = ti + j*dom->Gfz.s1 + tk*dom->Gfz.s2;
    CC = (ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    N = (ti+DOM_BUF) + (j+1+DOM_BUF)*dom->Gfz.s1b + (tk+DOM_BUF)*dom->Gfz.s2b;
    values[C + 9*pitch] =  0.;
    w_star[CC] += w_star[N] * 0.5 * nu * ddy;
  }
}

__global__ void wstar_coeffs_dirichlet_B(real bc, real *w_star, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfz.ks-DOM_BUF;  // iterator
  int C, CC;  // cell locations

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfz.in) && (tj < dom->Gfz.jn)) {
    C = ti + tj*dom->Gfz.s1 + k*dom->Gfz.s2;
    CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (k+DOM_BUF)*dom->Gfz.s2b;
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
    w_star[CC] = bc;
  }
}

__global__ void wstar_coeffs_dirichlet_T(real bc, real *w_star, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfz.ke-1-DOM_BUF;  // iterator
  int C, CC;  // cell locations

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfz.in) && (tj < dom->Gfz.jn)) {
    C = ti + tj*dom->Gfz.s1 + k*dom->Gfz.s2;
    CC = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz.s1b + (k+DOM_BUF)*dom->Gfz.s2b;
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
    w_star[CC] = bc;
  }
}

__global__ void ustar_coeffs_periodic_W(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfx.is-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfx.jn) && (tk < dom->Gfx.kn)) {
    C = i + tj*dom->Gfx.s1 + tk*dom->Gfx.s2;
    values[C + pitch * 5] += 0.5*nu*dt*ddx;
    values[C + pitch * 8] -= 0.5*nu*dt*ddx;
  }
}

__global__ void ustar_coeffs_periodic_E(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfx.ie-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfx.jn) && (tk < dom->Gfx.kn)) {
    C = i + tj*dom->Gfx.s1 + tk*dom->Gfx.s2;
    values[C + pitch * 4] -= 0.5*nu*dt*ddx;
    values[C + pitch * 7] += 0.5*nu*dt*ddx;
  }
}

__global__ void ustar_coeffs_periodic_S(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfx.js-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfx.kn) && (ti < dom->Gfx.in)) {
    C = ti + j*dom->Gfx.s1 + tk*dom->Gfx.s2;
    values[C + pitch * 3]  += 0.5*nu*dt*ddy;
    values[C + pitch * 10] -= 0.5*nu*dt*ddy;
  }
}

__global__ void ustar_coeffs_periodic_N(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfx.je-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfx.kn) && (ti < dom->Gfx.in)) {
    C = ti + j*dom->Gfx.s1 + tk*dom->Gfx.s2;
    values[C + pitch * 2] -= 0.5*nu*dt*ddy;
    values[C + pitch * 9] += 0.5*nu*dt*ddy;
  }
}

__global__ void ustar_coeffs_periodic_B(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfx.ks-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfx.in) && (tj < dom->Gfx.jn)) {
    C = ti + tj*dom->Gfx.s1 + k*dom->Gfx.s2;
    values[C + pitch * 1]  += 0.5*nu*dt*ddz;
    values[C + pitch * 12] -= 0.5*nu*dt*ddz;
  }
}

__global__ void ustar_coeffs_periodic_T(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfx.ke-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfx.in) && tj < (dom->Gfx.jn)) {
    C = ti + tj*dom->Gfx.s1 + k*dom->Gfx.s2;
    values[C + pitch * 0]  -= 0.5*nu*dt*ddz;
    values[C + pitch * 11] += 0.5*nu*dt*ddz;
  }
}

__global__ void vstar_coeffs_periodic_W(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfy.is-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfy.jn) && (tk < dom->Gfy.kn)) {
    C = i + tj*dom->Gfy.s1 + tk*dom->Gfy.s2;
    values[C + pitch * 5] += 0.5*nu*dt*ddx;
    values[C + pitch * 8] -= 0.5*nu*dt*ddx;
  }
}

__global__ void vstar_coeffs_periodic_E(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfy.ie-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfy.jn) && (tk < dom->Gfy.kn)) {
    C = i + tj*dom->Gfy.s1 + tk*dom->Gfy.s2;
    values[C + pitch * 4] -= 0.5*nu*dt*ddx;
    values[C + pitch * 7] += 0.5*nu*dt*ddx;
  }
}

__global__ void vstar_coeffs_periodic_S(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfy.js-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfy.kn) && (ti < dom->Gfy.in)) {
    C = ti + j*dom->Gfy.s1 + tk*dom->Gfy.s2;
    values[C + pitch * 3]  += 0.5*nu*dt*ddy;
    values[C + pitch * 10] -= 0.5*nu*dt*ddy;
  }
}

__global__ void vstar_coeffs_periodic_N(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfy.je-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfy.kn) && (ti < dom->Gfy.in)) {
    C = ti + j*dom->Gfy.s1 + tk*dom->Gfy.s2;
    values[C + pitch * 2] -= 0.5*nu*dt*ddy;
    values[C + pitch * 9] += 0.5*nu*dt*ddy;
  }
}

__global__ void vstar_coeffs_periodic_B(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfy.ks-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfy.in) && (tj < dom->Gfy.jn)) {
    C = ti + tj*dom->Gfy.s1 + k*dom->Gfy.s2;
    values[C + pitch * 1]  += 0.5*nu*dt*ddz;
    values[C + pitch * 12] -= 0.5*nu*dt*ddz;
  }
}

__global__ void vstar_coeffs_periodic_T(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfy.ke-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfy.in) && (tj < dom->Gfy.jn)) {
    C = ti + tj*dom->Gfy.s1 + k*dom->Gfy.s2;
    values[C + pitch * 0]  -= 0.5*nu*dt*ddz;
    values[C + pitch * 11] += 0.5*nu*dt*ddz;
  }
}

__global__ void wstar_coeffs_periodic_W(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfz.is-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfz.jn) && (tk < dom->Gfz.kn)) {
    C = i + tj*dom->Gfz.s1 + tk*dom->Gfz.s2;
    values[C + pitch * 5] += 0.5*nu*dt*ddx;
    values[C + pitch * 8] -= 0.5*nu*dt*ddx;
  }
}

__global__ void wstar_coeffs_periodic_E(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gfz.ie-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if((tj < dom->Gfz.jn) && (tk < dom->Gfz.kn)) {
    C = i + tj*dom->Gfz.s1 + tk*dom->Gfz.s2;
    values[C + pitch * 4] -= 0.5*nu*dt*ddx;
    values[C + pitch * 7] += 0.5*nu*dt*ddx;
  }
}

__global__ void wstar_coeffs_periodic_S(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfz.js-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfz.kn) && (ti < dom->Gfz.in)) {
    C = ti + j*dom->Gfz.s1 + tk*dom->Gfz.s2;
    values[C + pitch * 3]  += 0.5*nu*dt*ddy;
    values[C + pitch * 10] -= 0.5*nu*dt*ddy;
  }
}

__global__ void wstar_coeffs_periodic_N(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gfz.je-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if((tk < dom->Gfz.kn) && (ti < dom->Gfz.in)) {
    C = ti + j*dom->Gfz.s1 + tk*dom->Gfz.s2;
    values[C + pitch * 2] -= 0.5*nu*dt*ddy;
    values[C + pitch * 9] += 0.5*nu*dt*ddy;
  }
}

__global__ void wstar_coeffs_periodic_B(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfz.ks-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfz.in) && (tj < dom->Gfz.jn)) {
    C = ti + tj*dom->Gfz.s1 + k*dom->Gfz.s2;
    values[C + pitch * 1]  += 0.5*nu*dt*ddz;
    values[C + pitch * 12] -= 0.5*nu*dt*ddz;
  }
}

__global__ void wstar_coeffs_periodic_T(real nu, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gfz.ke-1-DOM_BUF;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if((ti < dom->Gfz.in) && (tj < dom->Gfz.jn)) {
    C = ti + tj*dom->Gfz.s1 + k*dom->Gfz.s2;
    values[C + pitch * 0]  -= 0.5*nu*dt*ddz;
    values[C + pitch * 11] += 0.5*nu*dt*ddz;
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

__global__ void coeffs_periodic_W(dom_struct *dom, int pitch, real *values,
  int *flag_u)
{
  int i = dom->Gcc.is;  // iterator
  int C, W;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    values[C + pitch * 5] -= (real)abs(flag_u[W])*ddx;
    values[C + pitch * 8] += (real)abs(flag_u[W])*ddx;
  }
}

__global__ void coeffs_periodic_E(dom_struct *dom, int pitch, real *values,
  int *flag_u)
{
  int i = dom->Gcc.ie-1;  // iterator
  int C, E;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    values[C + pitch * 4] += (real)abs(flag_u[E])*ddx;
    values[C + pitch * 7] -= (real)abs(flag_u[E])*ddx;
  }
}

__global__ void coeffs_periodic_S(dom_struct *dom, int pitch, real *values,
  int *flag_v)
{
  int j = dom->Gcc.js;  // iterator
  int C, S;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    S = ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    values[C + pitch * 3]  -= (real)abs(flag_v[S])*ddy;
    values[C + pitch * 10] += (real)abs(flag_v[S])*ddy;
  }
}

__global__ void coeffs_periodic_N(dom_struct *dom, int pitch, real *values,
  int *flag_v)
{
  int j = dom->Gcc.je-1;  // iterator
  int C, N;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    N = ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    values[C + pitch * 2] += (real)abs(flag_v[N])*ddy;
    values[C + pitch * 9] -= (real)abs(flag_v[N])*ddy;
  }
}

__global__ void coeffs_periodic_B(dom_struct *dom, int pitch, real *values,
  int *flag_w)
{
  int k = dom->Gcc.ks;  // iterator
  int C, B;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    B = ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b;
    values[C + pitch * 1]  -= (real)abs(flag_w[B])*ddz;
    values[C + pitch * 12] += (real)abs(flag_w[B])*ddz;
  }
}

__global__ void coeffs_periodic_T(dom_struct *dom, int pitch, real *values,
  int *flag_w)
{
  int k = dom->Gcc.ke-1;  // iterator
  int C, T;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    T = ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b;
    values[C + pitch * 0]  += (real)abs(flag_w[T])*ddz;
    values[C + pitch * 11] -= (real)abs(flag_w[T])*ddz;
  }
}

__global__ void coeffs_particle(dom_struct *dom, int pitch, real *values,
  int *phase)
{
  int i;  // iterator
  int C, CC;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.jn && tk < dom->Gcc.kn) {
    for(i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      CC = (i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc.s1b + (tk+DOM_BUF)*dom->Gcc.s2b;
      if(phase[CC] > -1) { // decouple all nodes inside particles
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
