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

#include "cuda_bluebottle.h"

// pressure; west; periodic
__global__ void BC_p_W_P(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[dom->Gcc._isb + tj*s1b + tk*s2b] = p[(dom->Gcc._ie-1) + tj*s1b + tk*s2b];
}

// pressure; west; Neumann
__global__ void BC_p_W_N(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[dom->Gcc._isb + tj*s1b + tk*s2b] = p[dom->Gcc._is + tj*s1b + tk*s2b];
}

// pressure; east; periodic
__global__ void BC_p_E_P(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = p[dom->Gcc._is + tj*s1b + tk*s2b];
}

// pressure; east; Neumann
__global__ void BC_p_E_N(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = p[(dom->Gcc._ie-1)
      + tj*s1b + tk*s2b];
}

// pressure; south; periodic
__global__ void BC_p_S_P(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + dom->Gcc._jsb*s1b + tk*s2b] = p[ti + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// pressure; south; Neumann
__global__ void BC_p_S_N(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + dom->Gcc._jsb*s1b + tk*s2b] = p[ti + dom->Gcc._js*s1b + tk*s2b];
}

// pressure; north; periodic
__global__ void BC_p_N_P(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = p[ti + dom->Gcc._js*s1b + tk*s2b];
}

// pressure; north; Neumann
__global__ void BC_p_N_N(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = p[ti
      + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// pressure; bottom; periodic
__global__ void BC_p_B_P(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + dom->Gcc._ksb*s2b] = p[ti + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

// pressure; bottom; Neumann
__global__ void BC_p_B_N(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b; 
  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + dom->Gcc._ksb*s2b] = p[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// pressure; top; periodic
__global__ void BC_p_T_P(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = p[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// pressure; top; Neumann
__global__ void BC_p_T_N(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = p[ti
      + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

// u-velocity; west; periodic
__global__ void BC_u_W_P(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = u[(dom->Gfx._ie-2) + tj*s1b + tk*s2b];
    u[dom->Gfx._is + tj*s1b + tk*s2b] = u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b];
  }
}

// u-velocity; west; Dirichlet
__global__ void BC_u_W_D(real *u, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = 2. * bc
      - u[(dom->Gfx._is+1) + tj*s1b + tk*s2b];
    u[dom->Gfx._is + tj*s1b + tk*s2b] = bc;
  }
}

// u-velocity; west; Neumann
__global__ void BC_u_W_N(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb))
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = u[dom->Gfx._is + tj*s1b + tk*s2b];
}

// u-velocity; west; Turbulent precursor
__global__ void BC_u_W_T(real *u, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = 2. * bc[tj + tk*dom->Gfx.jnb]
      - u[(dom->Gfx._is+1) + tj*s1b + tk*s2b];
    u[dom->Gfx._is + tj*s1b + tk*s2b] = bc[tj + tk*dom->Gfx.jnb];
  }
}

// u-velocity; east; periodic
__global__ void BC_u_E_P(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = u[(dom->Gfx._is+1) + tj*s1b + tk*s2b];
    u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b] = u[dom->Gfx._is + tj*s1b + tk*s2b];
  }
}

// u-velocity; east; Dirichlet
__global__ void BC_u_E_D(real *u, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = 2. * bc - u[(dom->Gfx._ie-2)
      + tj*s1b + tk*s2b];
    u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b] = bc;
  }
}

// u-velocity; east; Neumann
__global__ void BC_u_E_N(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb))
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = u[(dom->Gfx._ie-1)
      + tj*s1b + tk*s2b];
}

// u-velocity; east; Turbulent precursor
__global__ void BC_u_E_T(real *u, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = 2. * bc[tj + tk*dom->Gfx.jnb]
      - u[(dom->Gfx._ie-2) + tj*s1b + tk*s2b];
    u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b] = bc[tj + tk*dom->Gfx.jnb];
  }
}

// u-velocity; south; periodic
__global__ void BC_u_S_P(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = u[ti + (dom->Gfx._je-1)*s1b + tk*s2b];
  }
}

// u-velocity; south; Dirichlet
__global__ void BC_u_S_D(real *u, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = 8./3. * bc
      - 2. * u[ti + dom->Gfx._js*s1b + tk*s2b]
      + 1./3. * u[ti + (dom->Gfx._js+1)*s1b + tk*s2b];
  }
}

// u-velocity; south; Neumann
__global__ void BC_u_S_N(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = u[ti + dom->Gfx._js*s1b + tk*s2b];
}

// u-velocity; south; Turbulent precursor
__global__ void BC_u_S_T(real *u, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = 2. * bc[tk + ti*dom->Gfx.knb]
      - u[ti + dom->Gfx._js*s1b + tk*s2b];
  }
}

// u-velocity; north; periodic
__global__ void BC_u_N_P(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = u[ti
      + dom->Gfx._js*s1b + tk*s2b];
}

// u-velocity; north; Dirichlet
__global__ void BC_u_N_D(real *u, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = 8./3. * bc
      - 2. * u[ti + (dom->Gfx._je-1)*s1b + tk*s2b]
      + 1./3. * u[ti + (dom->Gfx._je-2)*s1b + tk*s2b];
  }
}

// u-velocity; north; Neumann
__global__ void BC_u_N_N(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = u[ti
      + (dom->Gfx._je-1)*s1b + tk*s2b];
}

// u-velocity; north; Turbulent precursor
__global__ void BC_u_N_T(real *u, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = 2. * bc[tk + ti*dom->Gfx.knb]
      - u[ti + (dom->Gfx._je-1)*s1b + tk*s2b];
  }
}

// u-velocity; bottom; periodic
__global__ void BC_u_B_P(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = u[ti + tj*s1b + (dom->Gfx._ke-1)*s2b];
}

// u-velocity; bottom; Dirichlet
__global__ void BC_u_B_D(real *u, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = 8./3. * bc
      - 2. * u[ti + tj*s1b + dom->Gfx._ks*s2b]
      + 1./3. * u[ti + tj*s1b + (dom->Gfx._ks+1)*s2b];
}

// u-velocity; bottom; Neumann
__global__ void BC_u_B_N(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = u[ti + tj*s1b + dom->Gfx._ks*s2b];
}

// u-velocity; bottom; Turbulent precursor
__global__ void BC_u_B_T(real *u, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb)) {
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = 2. * bc[ti + tj*dom->Gfx.inb]
      - u[ti + tj*s1b + dom->Gfx._ks*s2b];
  }
}

// u-velocity; top; periodic
__global__ void BC_u_T_P(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = u[ti
      + tj*s1b + dom->Gfx._ks*s2b];
}

// u-velocity; top; Dirichlet
__global__ void BC_u_T_D(real *u, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = 8./3. * bc
      - 2. * u[ti + tj*s1b + (dom->Gfx._ke-1)*s2b]
      + 1./3. * u[ti + tj*s1b + (dom->Gfx._ke-2)*s2b];
}

// u-velocity; top; Neumann
__global__ void BC_u_T_N(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = u[ti
      + tj*s1b + (dom->Gfx._ke-1)*s2b];
}

// u-velocity; top; Turbulent precursor
__global__ void BC_u_T_T(real *u, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb)) {
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = 2. * bc[ti + tj*dom->Gfx.inb]
      - u[ti + tj*s1b + (dom->Gfx._ke-1)*s2b];
  }
}

// v-velocity; west; periodic
__global__ void BC_v_W_P(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = v[(dom->Gfy._ie-1) + tj*s1b + tk*s2b];
}

// v-velocity; west; Dirichlet
__global__ void BC_v_W_D(real *v, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = 8./3. * bc
      - 2. * v[dom->Gfy._is + tj*s1b + tk*s2b]
      + 1./3. * v[(dom->Gfy._is+1) + tj*s1b + tk*s2b];
  }
}

// v-velocity; west; Neumann
__global__ void BC_v_W_N(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = v[dom->Gfy._is + tj*s1b + tk*s2b];
}

// v-velocity; west; Turbulent precursor
__global__ void BC_v_W_T(real *v, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = 2. * bc[tj + tk*dom->Gfy.jnb]
      - v[(dom->Gfy._is) + tj*s1b + tk*s2b];
  }
}

// v-velocity; east; periodic
__global__ void BC_v_E_P(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = v[dom->Gfy._is
      + tj*s1b + tk*s2b];
}

// v-velocity; east; Dirichlet
__global__ void BC_v_E_D(real *v, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = 8./3. * bc
      - 2. * v[(dom->Gfy._ie-1) + tj*s1b + tk*s2b]
      + 1./3. * v[(dom->Gfy._ie-2) + tj*s1b + tk*s2b];
  }
}

// v-velocity; east; Neumann
__global__ void BC_v_E_N(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = v[(dom->Gfy._ie-1)
      + tj*s1b + tk*s2b];
}

// v-velocity; east; Turbulent precursor
__global__ void BC_v_E_T(real *v, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = 2. * bc[tj + tk*dom->Gfy.jnb]
      - v[(dom->Gfy._ie-1) + tj*s1b + tk*s2b];
  }
}

// v-velocity; south; periodic
__global__ void BC_v_S_P(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = v[ti + (dom->Gfy._je-2)*s1b + tk*s2b];
    v[ti + dom->Gfy._js*s1b + tk*s2b] = v[ti + (dom->Gfy._je-1)*s1b + tk*s2b];
  }
}

// v-velocity; south; Dirichlet
__global__ void BC_v_S_D(real *v, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = 2. * bc - v[ti
      + (dom->Gfy._js+1)*s1b + tk*s2b];
    v[ti + dom->Gfy._js*s1b + tk*s2b] = bc;
  }
}

// v-velocity; south; Neumann
__global__ void BC_v_S_N(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb))
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = v[ti + dom->Gfy._js*s1b + tk*s2b];
}

// v-velocity; south; Turbulent precursor
__global__ void BC_v_S_T(real *v, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = 2. * bc[tk + ti*dom->Gfy.knb]
      - v[ti + (dom->Gfy._js+1)*s1b + tk*s2b];
    v[ti + dom->Gfy._js*s1b + tk*s2b] = bc[tk + ti*dom->Gfy.knb];
  }
}

// v-velocity; north; periodic
__global__ void BC_v_N_P(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = v[ti
      + (dom->Gfy._js+1)*s1b + tk*s2b];
    v[ti + (dom->Gfy._je-1)*s1b + tk*s2b] = v[ti
      + dom->Gfy._js*s1b + tk*s2b];
  }
}

// v-velocity; north; Dirichlet
__global__ void BC_v_N_D(real *v, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = 2. * bc - v[ti +
      (dom->Gfy._je-2)*s1b + tk*s2b];
    v[ti + (dom->Gfy._je-1)*s1b + tk*s2b] = bc;
  }
}

// v-velocity; north; Neumann
__global__ void BC_v_N_N(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb))
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = v[ti
      + (dom->Gfy._je-1)*s1b + tk*s2b];
}

// v-velocity; north; Turbulent precursor
__global__ void BC_v_N_T(real *v, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = 2. * bc[tk + ti*dom->Gfy.knb]
      - v[ti + (dom->Gfy._je-2)*s1b + tk*s2b];
    v[ti + (dom->Gfy._je-1)*s1b + tk*s2b] = bc[tk + ti*dom->Gfy.knb];
  }
}

// v-velocity; bottom; periodic
__global__ void BC_v_B_P(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = v[ti
      + tj*s1b + (dom->Gfy._ke-1)*s2b];
}

// v-velocity; bottom; Dirichlet
__global__ void BC_v_B_D(real *v, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = 8./3. * bc
      - 2. * v[ti + tj*s1b + dom->Gfy._ks*s2b]
      + 1./3. * v[ti + tj*s1b + (dom->Gfy._ks+1)*s2b];
}

// v-velocity; bottom; Neumann
__global__ void BC_v_B_N(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = v[ti + tj*s1b + dom->Gfy._ks*s2b];
}

// v-velocity; bottom; Turbulent precursor
__global__ void BC_v_B_T(real *v, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb)) {
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = 2. * bc[ti + tj*dom->Gfy.inb]
      - v[ti + tj*s1b + dom->Gfy._ks*s2b];
  }
}

// v-velocity; top; periodic
__global__ void BC_v_T_P(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = v[ti
      + tj*s1b + dom->Gfy._ks*s2b];
}

// v-velocity; top; Dirichlet
__global__ void BC_v_T_D(real *v, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = 8./3. * bc
      - 2. * v[ti + tj*s1b + (dom->Gfy._ke-1)*s2b]
      + 1./3. * v[ti + tj*s1b + (dom->Gfy._ke-2)*s2b];
}

// v-velocity; top; Neumann
__global__ void BC_v_T_N(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = v[ti
      + tj*s1b + (dom->Gfy._ke-1)*s2b];
}

// v-velocity; top; Turbulent precursor
__global__ void BC_v_T_T(real *v, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb)) {
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = 2. * bc[ti + tj*dom->Gfy.inb]
      - v[ti + tj*s1b + (dom->Gfy._ke-1)*s2b];
  }
}

// w-velocity; west; periodic
__global__ void BC_w_W_P(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = w[(dom->Gfz._ie-1) + tj*s1b + tk*s2b];
}

// w-velocity; west; Dirichlet
__global__ void BC_w_W_D(real *w, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = 8./3. * bc
      - 2. * w[dom->Gfz._is + tj*s1b + tk*s2b]
      + 1./3. * w[(dom->Gfz._is+1) + tj*s1b + tk*s2b];
}

// w-velocity; west; Neumann
__global__ void BC_w_W_N(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = w[dom->Gfz._is + tj*s1b + tk*s2b];
}

// w-velocity; west; Turbulent precursor
__global__ void BC_w_W_T(real *w, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb)) {
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = 2. * bc[tj + tk*dom->Gfz.jnb]
      - w[(dom->Gfz._is) + tj*s1b + tk*s2b];
  }
}

// w-velocity; east; periodic
__global__ void BC_w_E_P(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = w[dom->Gfz._is
      + tj*s1b + tk*s2b];
}

// w-velocity; east; Dirichlet
__global__ void BC_w_E_D(real *w, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = 8./3. * bc
      - 2. * w[(dom->Gfz._ie-1) + tj*s1b + tk*s2b]
      + 1./3. * w[(dom->Gfz._ie-2) + tj*s1b + tk*s2b];
}

// w-velocity; east; Neumann
__global__ void BC_w_E_N(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = w[(dom->Gfz._ie-1)
      + tj*s1b + tk*s2b];
}

// w-velocity; east; Turbulent precursor
__global__ void BC_w_E_T(real *w, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb)) {
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = 2. * bc[tj + tk*dom->Gfz.jnb]
      - w[(dom->Gfz._ie-1) + tj*s1b + tk*s2b];
  }
}

// w-velocity; south; periodic
__global__ void BC_w_S_P(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = w[ti + (dom->Gfz._je-1)*s1b + tk*s2b];
  }
}

// w-velocity; south; Dirichlet
__global__ void BC_w_S_D(real *w, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = 8./3. * bc
      - 2. * w[ti + dom->Gfz._js*s1b + tk*s2b]
      + 1./3. * w[ti + (dom->Gfz._js+1)*s1b + tk*s2b];
}

// w-velocity; south; Neumann
__global__ void BC_w_S_N(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = w[ti + dom->Gfz._js*s1b + tk*s2b];
}

// w-velocity; south; Turbulent precursor
__global__ void BC_w_S_T(real *w, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = 2. * bc[tk + ti*dom->Gfz.knb]
      - w[ti + dom->Gfz._js*s1b + tk*s2b];
  }
}

// w-velocity; north; periodic
__global__ void BC_w_N_P(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = w[ti
      + dom->Gfz._js*s1b + tk*s2b];
}

// w-velocity; north; Dirichlet
__global__ void BC_w_N_D(real *w, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = 8./3. * bc
      - 2. * w[ti + (dom->Gfz._je-1)*s1b + tk*s2b]
      + 1./3. * w[ti + (dom->Gfz._je-2)*s1b + tk*s2b];
}

// w-velocity; north; Neumann
__global__ void BC_w_N_N(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = w[ti
      + (dom->Gfz._je-1)*s1b + tk*s2b];
}

// w-velocity; north; Turbulent precursor
__global__ void BC_w_N_T(real *w, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = 2. * bc[tk + ti*dom->Gfz.knb]
      - w[ti + (dom->Gfz._je-1)*s1b + tk*s2b];
  }
}

// w-velocity; bottom; periodic
__global__ void BC_w_B_P(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = w[ti + tj*s1b + (dom->Gfz._ke-2)*s2b];
    w[ti + tj*s1b + dom->Gfz._ks*s2b] = w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b];
  }
}

// w-velocity; bottom; Dirichlet
__global__ void BC_w_B_D(real *w, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = 2. * bc - w[ti
      + tj*s1b + (dom->Gfz._ks+1)*s2b];
    w[ti + tj*s1b + dom->Gfz._ks*s2b] = bc;
  }
}

// w-velocity; bottom; Neumann
__global__ void BC_w_B_N(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb))
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = w[ti + tj*s1b + dom->Gfz._ks*s2b];
}

// w-velocity; bottom; Turbulent precursor
__global__ void BC_w_B_T(real *w, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = 2. * bc[ti + tj*dom->Gfz.inb]
      - w[ti + tj*s1b + (dom->Gfz._ks+1)*s2b];
    w[ti + tj*s1b + dom->Gfz._ks*s2b] = bc[ti + tj*dom->Gfz.inb];
  }
}

// w-velocity; top; periodic
__global__ void BC_w_T_P(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = w[ti
      + tj*s1b + (dom->Gfz._ks+1)*s2b];
    w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b] = w[ti
      + tj*s1b + dom->Gfz._ks*s2b];
  }
}

// w-velocity; top; Dirichlet
__global__ void BC_w_T_D(real *w, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = 2. * bc - w[ti + tj*s1b +
      (dom->Gfz._ke-2)*s2b];
    w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b] = bc;
  }
}

// w-velocity; top; Neumann
__global__ void BC_w_T_N(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb))
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = w[ti
      + tj*s1b + (dom->Gfz._ke-1)*s2b];
}

// w-velocity; top; Turbulent precursor
__global__ void BC_w_T_T(real *w, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = 2. * bc[ti + tj*dom->Gfz.inb]
      - w[ti + tj*s1b + (dom->Gfz._ke-2)*s2b];
    w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b] = bc[ti + tj*dom->Gfz.inb];
  }
}

__global__ void project_u(real *u_star, real *p, real rho_f, real dt,
  real *u, dom_struct *dom, real ddx, int *flag_u, int *phase)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tj < dom->Gfx._je && tk < dom->Gfx._ke) {
    for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
      real gradPhi = abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
      u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = (u_star[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b] - dt / rho_f * gradPhi);
    }
  }
}

__global__ void project_v(real *v_star, real *p, real rho_f, real dt,
  real *v, dom_struct *dom, real ddy, int *flag_v, int *phase)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tk < dom->Gfy._ke && ti < dom->Gfy._ie) {
    for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
      real gradPhi = abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
      v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = (v_star[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b] - dt / rho_f * gradPhi);
    }
  }
}

__global__ void project_w(real *w_star, real *p, real rho_f, real dt,
  real *w, dom_struct *dom, real ddz, int *flag_w, int *phase)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(ti < dom->Gfz._ie && tj < dom->Gfz._je) {
    for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
      real gradPhi = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);
      w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = (w_star[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b] - dt / rho_f * gradPhi);
    }
  }
}

__global__ void update_p_laplacian(real *Lp, real *p, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
    for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
      int C = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int W = (i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int E = (i+1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int S = i + (tj-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int N = i + (tj+1)*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int B = i + tj*dom->Gcc._s1b + (tk-1)*dom->Gcc._s2b;
      int T = i + tj*dom->Gcc._s1b + (tk+1)*dom->Gcc._s2b;
      real ddpdxx = (p[E]-2.*p[C]+p[W])/dom->dx/dom->dx;
      real ddpdyy = (p[N]-2.*p[C]+p[S])/dom->dy/dom->dy;
      real ddpdzz = (p[T]-2.*p[C]+p[B])/dom->dz/dom->dz;
      Lp[C] = ddpdxx+ddpdyy+ddpdzz;
    }
  }
}

__global__ void update_p(real *Lp, real *p0, real *p, real *phi,
  dom_struct *dom, real nu, real dt, int *phase)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
    for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
      int C = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      p[C] = (phase[C] < 0) * (p0[C] + phi[C] - 0.5*nu*dt*Lp[C]);
    }
  }
}

__global__ void copy_p_ghost(real *p, real *p_tmp, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.je-DOM_BUF && tk < dom->Gcc.ke-DOM_BUF) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      p[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc.s1b
        + (tk+DOM_BUF)*dom->Gcc.s2b] = p_tmp[i + tj*dom->Gcc.s1
        + tk*dom->Gcc.s2];
    }
  }
}

__global__ void copy_p_noghost(real *p_noghost, real *p_ghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.je-DOM_BUF && tk < dom->Gcc.ke-DOM_BUF) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      p_noghost[i + tj*dom->Gcc._s1 + tk*dom->Gcc._s2]
        = p_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b
        + (tk+DOM_BUF)*dom->Gcc._s2b];
    }
  }
}

__global__ void copy_u_ghost(real *u_ghost, real *u_noghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.je-DOM_BUF && tk < dom->Gfx.ke-DOM_BUF) {
    for(int i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      u_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx._s1b
        + (tk+DOM_BUF)*dom->Gfx._s2b] = u_noghost[i + tj*dom->Gfx._s1
        + tk*dom->Gfx._s2];
    }
  }
}

__global__ void copy_u_noghost(real *u_noghost, real *u_ghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.je-DOM_BUF && tk < dom->Gfx.ke-DOM_BUF) {
    for(int i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      u_noghost[i + tj*dom->Gfx._s1 + tk*dom->Gfx._s2]
        = u_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx._s1b
        + (tk+DOM_BUF)*dom->Gfx._s2b];
    }
  }
}

__global__ void copy_v_ghost(real *v_ghost, real *v_noghost, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.ke-DOM_BUF && ti < dom->Gfy.ie-DOM_BUF) {
    for(int j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      v_ghost[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy._s1b
        + (tk+DOM_BUF)*dom->Gfy._s2b] = v_noghost[ti + j*dom->Gfy._s1
        + tk*dom->Gfy._s2];
    }
  }
}

__global__ void copy_v_noghost(real *v_noghost, real *v_ghost, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.ke-DOM_BUF && ti < dom->Gfy.ie-DOM_BUF) {
    for(int j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      v_noghost[ti + j*dom->Gfy._s1 + tk*dom->Gfy._s2]
        = v_ghost[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy._s1b
        + (tk+DOM_BUF)*dom->Gfy._s2b];
    }
  }
}

__global__ void copy_w_ghost(real *w_ghost, real *w_noghost, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.ie-DOM_BUF && tj < dom->Gfz.je-DOM_BUF) {
    for(int k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      w_ghost[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b
        + (k+DOM_BUF)*dom->Gfz._s2b] = w_noghost[ti + tj*dom->Gfz._s1
        + k*dom->Gfz._s2];
    }
  }
}

__global__ void copy_w_noghost(real *w_noghost, real *w_ghost, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.ie-DOM_BUF && tj < dom->Gfz.je-DOM_BUF) {
    for(int k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      w_noghost[ti + tj*dom->Gfz._s1 + k*dom->Gfz._s2]
        = w_ghost[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b
        + (k+DOM_BUF)*dom->Gfz._s2b];
    }
  }
}

__global__ void copy_u_fluid(real *u_noghost, real *u_ghost, int *phase, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.je-DOM_BUF && tk < dom->Gfx.ke-DOM_BUF) {
    for(int i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      int boo = 1;
      if(phase[(i+DOM_BUF-1) + (tj+DOM_BUF)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      else if(phase[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      u_noghost[i + tj*dom->Gfx._s1 + tk*dom->Gfx._s2]
        = boo * u_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx._s1b
        + (tk+DOM_BUF)*dom->Gfx._s2b];
    }
  }
}

__global__ void copy_v_fluid(real *v_noghost, real *v_ghost, int *phase, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.ke-DOM_BUF && ti < dom->Gfy.ie-DOM_BUF) {
    for(int j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      int boo = 1;
      if(phase[(ti+DOM_BUF) + (j+DOM_BUF-1)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      else if(phase[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      v_noghost[ti + j*dom->Gfy._s1 + tk*dom->Gfy._s2]
        = boo * v_ghost[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy._s1b
        + (tk+DOM_BUF)*dom->Gfy._s2b];
    }
  }
}

__global__ void copy_w_fluid(real *w_noghost, real *w_ghost, int *phase, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.ie-DOM_BUF && tj < dom->Gfz.je-DOM_BUF) {
    for(int k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      int boo = 1;
      if(phase[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b + (k+DOM_BUF-1)*dom->Gcc._s2b] > -1) boo = 0;
      else if(phase[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b + (k+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      w_noghost[ti + tj*dom->Gfz._s1 + k*dom->Gfz._s2]
        = boo * w_ghost[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b
        + (k+DOM_BUF)*dom->Gfz._s2b];
    }
  }
}

#ifndef IMPLICIT
__global__ void u_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *p, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *u_star,
  dom_struct *dom, real dt0, real dt, int *phase)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
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
      s_u0[tj + tk*blockDim.x] = u0[(i-1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u1[tj + tk*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u2[tj + tk*blockDim.x] = u0[(i+1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((k >= dom->Gfy._ksb && k < dom->Gfy._keb)
      && (j >= dom->Gfy._jsb && j < dom->Gfy._jeb)) {
      s_v01[tj + tk*blockDim.x] = v0[(i-1) + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v12[tj + tk*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((k >= dom->Gfz._ksb && k < dom->Gfz._keb)
      && (j >= dom->Gfz._jsb && j < dom->Gfz._jeb)) {
      s_w01[tj + tk*blockDim.x] = w0[(i-1) + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tj + tk*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
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
#ifndef STOKESFLOW
      if(dt0 > 0) // Adams-Bashforth
        s_u_star[tj + tk*blockDim.x] += (-ab * s_c[tj + tk*blockDim.x]
          + ab0 * conv0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]);
      else        // forward Euler
        s_u_star[tj + tk*blockDim.x] += -s_c[tj + tk*blockDim.x];
#endif

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
      if(dt0 > 0) // Adams-Bashforth
        s_u_star[tj + tk*blockDim.x] += (ab * s_d[tj + tk*blockDim.x]
          - ab0 * diff0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]);
      else
        s_u_star[tj + tk*blockDim.x] += s_d[tj + tk*blockDim.x];

      // add on imposed pressure gradient
      s_u_star[tj + tk*blockDim.x] += f[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];

      // multiply by dt
      s_u_star[tj + tk*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += u111;

      // zero contribution inside particles
      s_u_star[tj + tk*blockDim.x] *=
        (phase[(i-1) + j*dom->Gcc._s1b + k*dom->Gcc._s2b] < 0
        && phase[i + j*dom->Gcc._s1b + k*dom->Gcc._s2b] < 0);
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
      diff[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = s_d[tj + tk*blockDim.x];
    }
  }
}
#endif

#ifndef IMPLICIT
__global__ void v_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *p, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *v_star,
  dom_struct *dom, real dt0, real dt, int *phase)
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
      s_v0[tk + ti*blockDim.x] = v0[i + (j-1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v1[tk + ti*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v2[tk + ti*blockDim.x] = v0[i + (j+1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((i >= dom->Gfz._isb && i < dom->Gfz._ieb)
      && (k >= dom->Gfz._ksb && k < dom->Gfz._keb)) {
      s_w01[tk + ti*blockDim.x] = w0[i + (j-1)*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tk + ti*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
    if((i >= dom->Gfx._isb && i < dom->Gfx._ieb)
      && (k >= dom->Gfx._ksb && k < dom->Gfx._keb)) {
      s_u01[tk + ti*blockDim.x] = u0[i + (j-1)*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u12[tk + ti*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
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
#ifndef STOKESFLOW
      if(dt0 > 0) // Adams-Bashforth
        s_v_star[tk + ti*blockDim.x] += (-ab * s_c[tk + ti*blockDim.x]
          + ab0 * conv0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b]);
      else
        s_v_star[tk + ti*blockDim.x] += -s_c[tk + ti*blockDim.x];
#endif

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
      if(dt0 > 0) // Adams-Bashforth
        s_v_star[tk + ti*blockDim.x] += (ab * s_d[tk + ti*blockDim.x]
          - ab0 * diff0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b]);
      else
        s_v_star[tk + ti*blockDim.x] += s_d[tk + ti*blockDim.x];

      // add on imposed pressure gradient
      s_v_star[tk + ti*blockDim.x] += f[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];

      // multiply by dt
      s_v_star[tk + ti*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += v111;

      // zero contribution inside particles
      s_v_star[tk + ti*blockDim.x] *=
        (phase[i + (j-1)*dom->Gcc._s1b + k*dom->Gcc._s2b] < 0
        && phase[i + j*dom->Gcc._s1b + k*dom->Gcc._s2b] < 0);
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
      diff[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = s_d[tk + ti*blockDim.x];
    }
  }
}
#endif

#ifndef IMPLICIT
__global__ void w_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *p, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *w_star,
  dom_struct *dom, real dt0, real dt, int *phase)
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
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv0
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
      s_w0[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u01[ti + tj*blockDim.x] = u0[i + j*dom->Gfx._s1b + (k-1)*dom->Gfx._s2b];
      s_u12[ti + tj*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v01[ti + tj*blockDim.x] = v0[i + j*dom->Gfy._s1b + (k-1)*dom->Gfy._s2b];
      s_v12[ti + tj*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
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
#ifndef STOKESFLOW
      if(dt0 > 0) // Adams-Bashforth
        s_w_star[ti + tj*blockDim.x] += (-ab * s_c[ti + tj*blockDim.x]
          + ab0 * conv0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]);
      else        // forward Euler
        s_w_star[ti + tj*blockDim.x] += -s_c[ti + tj*blockDim.x];
#endif

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
      if(dt0 > 0) // Adams-Bashforth
        s_w_star[ti + tj*blockDim.x] += (ab * s_d[ti + tj*blockDim.x]
          - ab0 * diff0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]);
      else        // forward Euler
        s_w_star[ti + tj*blockDim.x] += s_d[ti + tj*blockDim.x];

      // add on imposed pressure gradient
      s_w_star[ti + tj*blockDim.x] += f[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];

      // multiply by dt
      s_w_star[ti + tj*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += w111;

      // zero contribution inside particles
      s_w_star[ti + tj*blockDim.x] *=
        (phase[i + j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b] < 0
        && phase[i + j*dom->Gcc._s1b + k*dom->Gcc._s2b] < 0);
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      w_star[i+ j*dom->Gfz._s1b + k*dom->Gfz._s2b]
        = s_w_star[ti + tj*blockDim.x];
      conv[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_c[ti + tj*blockDim.x];
      diff[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_d[ti + tj*blockDim.x];
    }
  }
}
#endif

__global__ void forcing_reset_x(real *fx, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0.;
    }
  }
}

__global__ void forcing_reset_y(real *fy, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0.;
    }
  }
}

__global__ void forcing_reset_z(real *fz, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
    if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
      fz[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 0.;
    }
  }
}

__global__ void forcing_add_x_const(real val, real *fx, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] += val;
    }
  }
}

__global__ void forcing_add_y_const(real val, real *fy, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] += val;
    }
  }
}

__global__ void forcing_add_z_const(real val, real *fz, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
    if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
      fz[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] += val;
    }
  }
}


__global__ void forcing_add_x_field(real scale, real *val, real *fx,
  dom_struct *dom, int *phase)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b]
        += scale * val[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b];
    }
  }
}

__global__ void forcing_add_y_field(real scale, real *val, real *fy,
  dom_struct *dom, int *phase)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b]
        += scale * val[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b];
    }
  }
}

__global__ void forcing_add_z_field(real scale, real *val, real *fz,
  dom_struct *dom, int *phase)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
    if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
      fz[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b]
        += scale * val[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
  }
}

__global__ void surf_int_x_copy(real *u_star, real *u_star_tmp,
  dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    int C = dom->Gfx.is + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    int CC = 0 + tj + tk*dom->Gfx.jn;
    u_star_tmp[CC] = -u_star[C];
    C = dom->Gfx.ie-1 + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    CC = dom->Gfx.jn*dom->Gfx.kn + tj + tk*dom->Gfx.jn;
    u_star_tmp[CC] = u_star[C];
  }
}

__global__ void surf_int_y_copy(real *v_star, real *v_star_tmp,
  dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    int C = (ti+DOM_BUF) + dom->Gfy.js*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    int CC = ti + 0 + tk*dom->Gfy.in;
    v_star_tmp[CC] = -v_star[C];
    C = (ti+DOM_BUF) + (dom->Gfy.je-1)*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    CC = ti + dom->Gfy.in*dom->Gfy.kn + tk*dom->Gfy.in;
    v_star_tmp[CC] = v_star[C];
  }
}

__global__ void surf_int_z_copy(real *w_star, real *w_star_tmp,
  dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    int C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + dom->Gfz.ks*dom->Gfz._s2b;
    int CC = ti + tj*dom->Gfz.in + 0;
    w_star_tmp[CC] = -w_star[C];
    C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + (dom->Gfz.ke-1)*dom->Gfz._s2b;
    CC = ti + tj*dom->Gfz.in + dom->Gfz.in*dom->Gfz.jn;
    w_star_tmp[CC] = w_star[C];
  }
}

__global__ void plane_eps_x_W(real eps, real *u_star, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y; 

  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    int C = dom->Gfx.is + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    u_star[C] = u_star[C] + eps;
  }
}

__global__ void plane_eps_x_E(real eps, real *u_star, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y; 

  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    int C = dom->Gfx.ie-1 + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    u_star[C] = u_star[C] - eps;
  }
}

__global__ void plane_eps_y_S(real eps, real *v_star, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    int C = (ti+DOM_BUF) + (dom->Gfy.js)*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    v_star[C] = v_star[C] + eps;
  }
}

__global__ void plane_eps_y_N(real eps, real *v_star, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    int C = (ti+DOM_BUF) + (dom->Gfy.je-1)*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    v_star[C] = v_star[C] - eps;
  }
}

__global__ void plane_eps_z_B(real eps, real *w_star, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    int C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + (dom->Gfz.ks)*dom->Gfz._s2b;
    w_star[C] = w_star[C] + eps;
  }
}

__global__ void plane_eps_z_T(real eps, real *w_star, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    int C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + (dom->Gfz.ke-1)*dom->Gfz._s2b;
    w_star[C] = w_star[C] - eps;
  }
}

__global__ void move_parts_a(dom_struct *dom, part_struct *parts, int nparts,
  real dt, real dt0, g_struct g, gradP_struct gradP, real rho_f, real ttime)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x; // particle number
  real vol = 4./3. * PI * parts[pp].r*parts[pp].r*parts[pp].r;
  real m = vol * parts[pp].rho;

  if(pp < nparts) {
    if(parts[pp].translating) {
      // update linear accelerations
      parts[pp].udot = (parts[pp].Fx + parts[pp].kFx + parts[pp].iFx
        + parts[pp].aFx - vol*gradP.x) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.x;
      parts[pp].vdot = (parts[pp].Fy + parts[pp].kFy + parts[pp].iFy
        + parts[pp].aFy - vol*gradP.y) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.y;
      parts[pp].wdot = (parts[pp].Fz + parts[pp].kFz + parts[pp].iFz
        + parts[pp].aFz - vol*gradP.z) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.z;

      // update linear velocities
      parts[pp].u = parts[pp].u0 + 0.5*dt*(parts[pp].udot + parts[pp].udot0);
      parts[pp].v = parts[pp].v0 + 0.5*dt*(parts[pp].vdot + parts[pp].vdot0);
      parts[pp].w = parts[pp].w0 + 0.5*dt*(parts[pp].wdot + parts[pp].wdot0);

      // do not update position
    }
    if(parts[pp].rotating) {
      // update angular accelerations
      real I = 0.4 * m * parts[pp].r*parts[pp].r;
      parts[pp].oxdot = (parts[pp].Lx + parts[pp].iLx + parts[pp].aLx) / I;
      parts[pp].oydot = (parts[pp].Ly + parts[pp].iLy + parts[pp].aLy) / I;
      parts[pp].ozdot = (parts[pp].Lz + parts[pp].iLz + parts[pp].aLz) / I;

      // update angular velocities
      parts[pp].ox = parts[pp].ox0 + 0.5*dt*(parts[pp].oxdot + parts[pp].oxdot0);
      parts[pp].oy = parts[pp].oy0 + 0.5*dt*(parts[pp].oydot + parts[pp].oydot0);
      parts[pp].oz = parts[pp].oz0 + 0.5*dt*(parts[pp].ozdot + parts[pp].ozdot0);
    }
  }
}

__global__ void move_parts_b(dom_struct *dom, part_struct *parts, int nparts,
  real dt, real dt0, g_struct g, gradP_struct gradP, real rho_f, real ttime)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x; // particle number
  real vol = 4./3. * PI * parts[pp].r*parts[pp].r*parts[pp].r;
  real m = vol * parts[pp].rho;

  if(pp < nparts) {
    if(parts[pp].translating) {
      // update linear accelerations
      parts[pp].udot = (parts[pp].Fx + parts[pp].kFx + parts[pp].iFx
        + parts[pp].aFx - vol*gradP.x) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.x;
      parts[pp].vdot = (parts[pp].Fy + parts[pp].kFy + parts[pp].iFy
        + parts[pp].aFy - vol*gradP.y) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.y;
      parts[pp].wdot = (parts[pp].Fz + parts[pp].kFz + parts[pp].iFz
          + parts[pp].aFz - vol*gradP.z) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.z;

      // update linear velocities
      parts[pp].u = parts[pp].u0 + 0.5*dt*(parts[pp].udot + parts[pp].udot0);
      parts[pp].v = parts[pp].v0 + 0.5*dt*(parts[pp].vdot + parts[pp].vdot0);
      parts[pp].w = parts[pp].w0 + 0.5*dt*(parts[pp].wdot + parts[pp].wdot0);

      // update position (trapezoidal rule)
      parts[pp].x = parts[pp].x0 + 0.5*dt*(parts[pp].u + parts[pp].u0);
      if(parts[pp].x < dom->xs) parts[pp].x = parts[pp].x + dom->xl;
      else if(parts[pp].x > dom->xe) parts[pp].x = parts[pp].x - dom->xl;

      parts[pp].y = parts[pp].y0 + 0.5*dt*(parts[pp].v + parts[pp].v0);
      if(parts[pp].y < dom->ys) parts[pp].y = parts[pp].y + dom->yl;
      else if(parts[pp].y > dom->ye) parts[pp].y = parts[pp].y - dom->yl;

      parts[pp].z = parts[pp].z0 + 0.5*dt*(parts[pp].w + parts[pp].w0);
      if(parts[pp].z < dom->zs) parts[pp].z = parts[pp].z + dom->zl;
      else if(parts[pp].z > dom->ze) parts[pp].z = parts[pp].z - dom->zl;

      // store for next time step
      parts[pp].x0 = parts[pp].x;
      parts[pp].y0 = parts[pp].y;
      parts[pp].z0 = parts[pp].z;
      parts[pp].u0 = parts[pp].u;
      parts[pp].v0 = parts[pp].v;
      parts[pp].w0 = parts[pp].w;
      parts[pp].udot0 = parts[pp].udot;
      parts[pp].vdot0 = parts[pp].vdot;
      parts[pp].wdot0 = parts[pp].wdot;
    }
    if(parts[pp].rotating) {
      // update angular accelerations
      real I = 0.4 * m * parts[pp].r*parts[pp].r;
      parts[pp].oxdot = (parts[pp].Lx + parts[pp].iLx + parts[pp].aLx) / I;
      parts[pp].oydot = (parts[pp].Ly + parts[pp].iLy + parts[pp].aLy) / I;
      parts[pp].ozdot = (parts[pp].Lz + parts[pp].iLz + parts[pp].aLz) / I;

      // update angular velocities
      parts[pp].ox = parts[pp].ox0 + 0.5*dt*(parts[pp].oxdot + parts[pp].oxdot0);
      parts[pp].oy = parts[pp].oy0 + 0.5*dt*(parts[pp].oydot + parts[pp].oydot0);
      parts[pp].oz = parts[pp].oz0 + 0.5*dt*(parts[pp].ozdot + parts[pp].ozdot0);

      /* update basis vectors */
      // calculate rotation magnitude (trapezoidal rule)
      real mag = 0.5*sqrt(parts[pp].ox*parts[pp].ox + parts[pp].oy*parts[pp].oy
        + parts[pp].oz*parts[pp].oz);
      mag += 0.5*sqrt(parts[pp].ox0*parts[pp].ox0 + parts[pp].oy0*parts[pp].oy0
        + parts[pp].oz0*parts[pp].oz0);
      // calculate normalized rotation axis
      real X = 0;
      real Y = 0;
      real Z = 0;
      if(mag > 0) {
        X = 0.5 * (parts[pp].ox + parts[pp].ox0) / mag;
        Y = 0.5 * (parts[pp].oy + parts[pp].oy0) / mag;
        Z = 0.5 * (parts[pp].oz + parts[pp].oz0) / mag;
      }
      // calculate rotation quaternion
      real theta = mag * dt;
      real qr = cos(0.5*theta);
      real qi = X * sin(0.5*theta);
      real qj = Y * sin(0.5*theta);
      real qk = Z * sin(0.5*theta);
      // compute quaternion conjugation to apply rotation to basis vectors
      rotate(qr, qi, qj, qk, &parts[pp].axx, &parts[pp].axy, &parts[pp].axz);
      rotate(qr, qi, qj, qk, &parts[pp].ayx, &parts[pp].ayy, &parts[pp].ayz);
      rotate(qr, qi, qj, qk, &parts[pp].azx, &parts[pp].azy, &parts[pp].azz);

      // store for next time step
      parts[pp].ox0 = parts[pp].ox;
      parts[pp].oy0 = parts[pp].oy;
      parts[pp].oz0 = parts[pp].oz;
      parts[pp].oxdot0 = parts[pp].oxdot;
      parts[pp].oydot0 = parts[pp].oydot;
      parts[pp].ozdot0 = parts[pp].ozdot;
    }
  }
}

__device__ void rotate(real qr, real qi, real qj, real qk,
  real *pi, real *pj, real *pk)
{
  real Pr = *pi*qi + *pj*qj + *pk*qk;
  real Pi = *pi*qr - *pj*qk + *pk*qj;
  real Pj = *pi*qk + *pj*qr - *pk*qi;
  real Pk = -*pi*qj + *pj*qi + *pk*qr;

  *pi = qr*Pi + qi*Pr + qj*Pk - qk*Pj;
  *pj = qr*Pj - qi*Pk + qj*Pr + qk*Pi;
  *pk = qr*Pk + qi*Pj - qj*Pi + qk*Pr;
}

__global__ void collision_init(part_struct *parts, int nparts)
{
  int j = threadIdx.x + blockIdx.x*blockDim.x;
  if(j < nparts) {
    parts[j].iFx = 0.;
    parts[j].iFy = 0.;
    parts[j].iFz = 0.;
    parts[j].iLx = 0.;
    parts[j].iLy = 0.;
    parts[j].iLz = 0.;
  }
}

__global__ void init(int *vector, int N, int val)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i < N) {
    vector[i] = val;
  }
}

__global__ void bin_fill(int *partInd, int *partBin, int nparts,
                  part_struct *parts, dom_struct *binDom, BC bc)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  int c;
  int ibin, jbin, kbin;

  // find the correct bin index for each part and store it
  if (pp < nparts) {
    ibin = floor((parts[pp].x - binDom->xs)/binDom->dx);
    jbin = floor((parts[pp].y - binDom->ys)/binDom->dy);
    kbin = floor((parts[pp].z - binDom->zs)/binDom->dz);
    c = ibin + jbin*binDom->Gcc.s1 + kbin*binDom->Gcc.s2;

    partInd[pp] = pp;         // index of particle
    partBin[pp] = c;          // bin index
    parts[pp].bin = c;        // bin index (stored in particle)
  }
}

__global__ void bin_partCount(int *binCount, int *binStart, int *binEnd,
                              dom_struct *binDom, BC bc, int nBins)
{
  int bin = threadIdx.x + blockIdx.x*blockDim.x;

  // fill binCount
  if (bin < nBins) {
    binCount[bin] = binEnd[bin] - binStart[bin];
  }
}

__global__ void bin_start(int *binStart, int *binEnd, int *partBin, int nparts)
{
  // This kernel function was adapted from NVIDIA CUDA 5.5 Examples
  // This software contains source code provided by NVIDIA Corporation
  extern __shared__ int sharedBin[];    //blockSize + 1
  int index = threadIdx.x + blockIdx.x*blockDim.x;
  int bin;

  // for a given bin index, the previous bins's index is stored in sharedBin
  if (index < nparts) {
    bin = partBin[index]; 

    // Load bin data into shared memory so that we can look
    // at neighboring particle's hash value without loading
    // two bin values per thread
    sharedBin[threadIdx.x + 1] = bin;

    if (index > 0 && threadIdx.x == 0) {
      // first thread in block must load neighbor particle bin
      sharedBin[0] = partBin[index - 1];
    }
  }
  __syncthreads();

  if (index < nparts) {
    // If this particle has a different cell index to the previous
    // particle then it must be the first particle in the cell,
    // so store the index of this particle in the cell.
    // As it isn't the first particle, it must also be the cell end of
    // the previous particle's cell
    bin = partBin[index]; 

    if (index == 0 || bin != sharedBin[threadIdx.x]) {
    binStart[bin] = index;

        if (index > 0)
            binEnd[sharedBin[threadIdx.x]] = index;
    }

    if (index == nparts - 1)
    {
        binEnd[bin] = index + 1;
    }
  }
}

__global__ void collision_parts(part_struct *parts, int nparts,
  dom_struct *dom, real eps, real mu, BC bc, int *binStart, int *binEnd,
  int *partBin, int *partInd, dom_struct *binDom, int interactionLength)
{
  int index = threadIdx.x + blockIdx.x*blockDim.x;

  if (index < nparts) {

    int i = partInd[index];
    int bin = partBin[index];

    int kbin = floorf(bin/binDom->Gcc.s2);
    int jbin = floorf((bin - kbin*binDom->Gcc.s2)/binDom->Gcc.s1);
    int ibin = bin - kbin*binDom->Gcc.s2 - jbin*binDom->Gcc.s1;

    int l, m, n;                          // adjacent bin iterators
    int target, j;                        // target indices
    int adjBin, adjStart, adjEnd;         // adjacent bin stuff
    int iStride, kStride, jStride;        // how to get to adjacent bin

    // predefine face locations 
    // -1, -2 due to local vs global indexing and defiinition of dom_struct
    int fW = binDom->Gcc.is - 1;
    int fE = binDom->Gcc.ie - 2;
    int fS = binDom->Gcc.js - 1;
    int fN = binDom->Gcc.je - 2;
    int fB = binDom->Gcc.ks - 1;
    int fT = binDom->Gcc.ke - 2;

    // size checks
    int xnBin = (binDom->xn > 2);
    int ynBin = (binDom->yn > 2);
    int znBin = (binDom->zn > 2);

    // loop over adjacent bins and take care of periodic conditions 
    for (n = -1; n <= 1; n++) {
      // if on a face and not periodic, continue
      // if on a face and periodic but only 2 bins, continue
      if ((n == -1 && kbin == fB && bc.uB != PERIODIC) || 
          (n == 1 && kbin == fT && bc.uT != PERIODIC) ||
          (n == -1 && kbin == fB && bc.uB == PERIODIC && znBin == 0) ||
          (n == 1 && kbin == fT && bc.uT == PERIODIC && znBin == 0)) {
        continue;
      // if on a face and periodic, flip to other side
      } else if (n == -1 && kbin == fB && bc.uB == PERIODIC) {
        kStride = fT*binDom->Gcc.s2;
      } else if (n == 1 && kbin == fT && bc.uT == PERIODIC) {
        kStride = fB*binDom->Gcc.s2;
      // else, we are in the middle, do nothing special
      } else {
        kStride = (kbin + n)*binDom->Gcc.s2;
      }

      for (m = -1; m <= 1; m++) {
        if ((m == -1 && jbin == fS && bc.uS != PERIODIC) ||
            (m == 1 && jbin == fN && bc.uN != PERIODIC) ||
            (m == -1 && jbin == fS && bc.uS == PERIODIC && ynBin == 0) ||
            (m == 1 && jbin == fN && bc.uN == PERIODIC && ynBin == 0)) {
          continue;
        } else if (m == -1 && jbin == fS && bc.uS == PERIODIC) {
          jStride = fN*binDom->Gcc.s1;  
        } else if (m == 1 && jbin == fN && bc.uN == PERIODIC) {
          jStride = fS*binDom->Gcc.s1;
        } else {
          jStride = (jbin + m)*binDom->Gcc.s1;
        }

        for (l = -1; l <= 1; l++) {
          if ((l == -1 && ibin == fW && bc.uW != PERIODIC) ||
              (l == 1 && ibin == fE && bc.uE != PERIODIC) ||
              (l == -1 && ibin == fW && bc.uW == PERIODIC && xnBin == 0) ||
              (l == 1 && ibin == fE && bc.uE == PERIODIC && xnBin == 0)) {
            continue;
          } else if (l == -1 && ibin == fW && bc.uW == PERIODIC) {
            iStride = fE;
          } else if (l == 1 && ibin == fE && bc.uE == PERIODIC) {
            iStride = fW;
          } else {
            iStride = ibin + l;
          }

          adjBin = iStride + jStride + kStride; 
          adjStart = binStart[adjBin];        // find start and end of bins
          adjEnd = binEnd[adjBin];
          if (adjStart != -1) {               // if bin is not empty
            for (target = adjStart; target < adjEnd; target++) {
              j = partInd[target];
              if (j != i) {                   // if its not original part

                // calculate forces

                real ai = parts[i].r;
                real aj = parts[j].r;
                real B = aj / ai;
                real hN = interactionLength;

                real ux, uy, uz;
                real rx, rx1, rx2, ry, ry1, ry2, rz, rz1, rz2, r;
                real h, ah, lnah;
                real nx, ny, nz, udotn;
                real unx, uny, unz, utx, uty, utz, ut;
                real tx, ty, tz, t, bx, by, bz, b;
                real omegax, omegay, omegaz, omega;
                real ocrossnx, ocrossny, ocrossnz;
                real utcrossnx, utcrossny, utcrossnz;
                real opB;
                real Fnx, Fny, Fnz, Ftx, Fty, Ftz, Lox, Loy, Loz;

                real xi = parts[i].x;
                real xj = parts[j].x;
                // check for neighbors across the domain when using periodic
                // boundaries
                rx = xi - xj;
                rx1 = xi - (xj + dom->xl);
                rx2 = xi - (xj - dom->xl);
                if(rx1*rx1 < rx*rx) rx = rx1;
                if(rx2*rx2 < rx*rx) rx = rx2;
                rx = (bc.uW == PERIODIC) * rx + (bc.uW != PERIODIC) * (xi - xj);

                real yi = parts[i].y;
                real yj = parts[j].y;
                // check for neighbors across the domain when using periodic
                // boundaries
                ry = yi - yj;
                ry1 = yi - (yj + dom->yl);
                ry2 = yi - (yj - dom->yl);
                if(ry1*ry1 < ry*ry) ry = ry1;
                if(ry2*ry2 < ry*ry) ry = ry2;
                ry = (bc.vS == PERIODIC) * ry + (bc.vS != PERIODIC) * (yi - yj);

                real zi = parts[i].z;
                real zj = parts[j].z;
                // check for neighbors across the domain when using periodic
                // boundaries
                rz = zi - zj;
                rz1 = zi - (zj + dom->zl);
                rz2 = zi - (zj - dom->zl);
                if(rz1*rz1 < rz*rz) rz = rz1;
                if(rz2*rz2 < rz*rz) rz = rz2;
                rz = (bc.wB == PERIODIC) * rz + (bc.wB != PERIODIC) * (zi - zj);

                ux = parts[i].u - parts[j].u;
                uy = parts[i].v - parts[j].v;
                uz = parts[i].w - parts[j].w;

                r = sqrt(rx*rx + ry*ry + rz*rz);

                omegax = parts[i].ox - parts[j].ox;
                omegay = parts[i].oy - parts[j].oy;
                omegaz = parts[i].oz - parts[j].oz;

                omega = sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz);

                h = r - ai - aj;

                nx = rx / r;
                ny = ry / r;
                nz = rz / r;

                udotn = ux * nx + uy * ny + uz * nz;

                unx = udotn * nx;
                uny = udotn * ny;
                unz = udotn * nz;

                utx = ux - unx;
                uty = uy - uny;
                utz = uz - unz;

                ut = sqrt(utx*utx + uty*uty + utz*utz);

                if(ut > 0) {
                  tx = utx / ut;
                  ty = uty / ut;
                  tz = utz / ut;

                  bx = ny*tz - nz*ty;
                  by = -nx*tz + nz*tx;
                  bz = nx*ty - ny*tx;

                  b = sqrt(bx*bx + by*by + bz*bz);

                  bx = bx / b;
                  by = by / b;
                  bz = bz / b;

                } else if(omega > 0) {
                  bx = omegax / omega;
                  by = omegay / omega;
                  bz = omegaz / omega;

                  tx = by*nz - bz*ny;
                  ty = -bx*nz + bz*nx;
                  tz = bx*ny - by*nx;

                  t = sqrt(tx*tx + ty*ty + tz*tz);

                  tx = tx / t;
                  ty = ty / t;
                  tz = tz / t;
                } else {
                  tx = 1.;
                  ty = 0.;
                  tz = 0.;

                  bx = ny*tz - nz*ty;
                  by = -nx*tz + nz*tx;
                  bz = nx*ty - ny*tx;

                  b = sqrt(bx*bx + by*by + bz*bz);

                  bx = bx / b;
                  by = by / b;
                  bz = bz / b;
                }

                opB = 1 + B;

                ocrossnx = omegay*nz - omegaz*ny;
                ocrossny = -omegax*nz + omegaz*nx;
                ocrossnz = omegax*ny - omegay*nx;

                utcrossnx = uty*nz - utz*ny;
                utcrossny = -utx*nz + utz*nx;
                utcrossnz = utx*ny - uty*nx;

                if(h < hN && h > 0) {
                  if(h < eps*parts[i].r) h = eps*parts[i].r;
                  ah = ai/h - ai/hN;
                  lnah = log(hN/h);
                  Fnx = -1. * B*B / (opB*opB) * ah
                    - B*(1.+7.*B+B*B)/(5.*opB*opB*opB)*lnah;
                  Fny = Fnx;
                  Fnz = Fnx;
                  Fnx *= 6.*PI*mu*ai*unx;
                  Fny *= 6.*PI*mu*ai*uny;
                  Fnz *= 6.*PI*mu*ai*unz;

                  Ftx = -6.*PI*mu*ai*utx*4.*B*(2.+B+2.*B*B)
                    /(15.*opB*opB*opB)*lnah;
                  Fty = -6.*PI*mu*ai*uty*4.*B*(2.+B+2.*B*B)
                    /(15.*opB*opB*opB)*lnah;
                  Ftz = -6.*PI*mu*ai*utz*4.*B*(2.+B+2.*B*B)
                    /(15.*opB*opB*opB)*lnah;
                  Ftx += 8.*PI*mu*ai*ai*ocrossnx*B*(4.+B)/(10.*opB*opB)*lnah;
                  Fty += 8.*PI*mu*ai*ai*ocrossny*B*(4.+B)/(10.*opB*opB)*lnah;
                  Ftz += 8.*PI*mu*ai*ai*ocrossnz*B*(4.+B)/(10.*opB*opB)*lnah;

                  Lox = -8.*PI*mu*ai*ai*utcrossnx*B*(4.+B)/(10.*opB*opB)*lnah;
                  Loy = -8.*PI*mu*ai*ai*utcrossny*B*(4.+B)/(10.*opB*opB)*lnah;
                  Loz = -8.*PI*mu*ai*ai*utcrossnz*B*(4.+B)/(10.*opB*opB)*lnah;
                  Lox += -8.*PI*mu*ai*ai*ai*omegax*2.*B/(5.*opB)*lnah;
                  Loy += -8.*PI*mu*ai*ai*ai*omegay*2.*B/(5.*opB)*lnah;
                  Loz += -8.*PI*mu*ai*ai*ai*omegaz*2.*B/(5.*opB)*lnah;
                } else {
                  ah = 0;
                  lnah = 0;
                  Fnx = 0;
                  Fny = 0;
                  Fnz = 0;
                  Ftx = 0;
                  Fty = 0;
                  Ftz = 0;
                  Lox = 0;
                  Loy = 0;
                  Loz = 0;
                }

          /** TODO implement variable alpha damping as is done for the walls **/
                if(h < 0) {
                  ah = 0;
                  lnah = 0;
                  real denom = 0.75*((1-parts[i].sigma*parts[i].sigma)
                    /parts[i].E
                    + (1-parts[j].sigma*parts[j].sigma)/parts[j].E)
                    *sqrt(1./ai + 1./aj);
                    
                  Fnx = (sqrt(-h*h*h)/denom
                    - sqrt(4./3.*PI*ai*ai*ai*parts[i].rho/denom*sqrt(-h))*udotn)
                    *nx;
                  Fny = (sqrt(-h*h*h)/denom
                    - sqrt(4./3.*PI*ai*ai*ai*parts[i].rho/denom*sqrt(-h))*udotn)
                    *ny;
                  Fnz = (sqrt(-h*h*h)/denom
                    - sqrt(4./3.*PI*ai*ai*ai*parts[i].rho/denom*sqrt(-h))*udotn)
                    *nz;

                }

                // assign forces
                parts[i].iFx += Fnx + Ftx;
                parts[i].iFy += Fny + Fty;
                parts[i].iFz += Fnz + Ftz;
                parts[i].iLx += Lox;
                parts[i].iLy += Loy;
                parts[i].iLz += Loz;
              }
            }
          }
        }
      }
    }
  }
}

__global__ void collision_walls(dom_struct *dom, part_struct *parts,
  int nparts, BC bc, real eps, real mu, real rhof, real nu,
  int interactionLength, real dt)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  /**** parallelize this further by using a CUDA block for each wall ****/

  if(i < nparts) {
    real dx = 0;
    real dy = 0;
    real dz = 0;
    real Un, Utx, Uty, Utz;
    real omx, omy, omz;

    real ai = parts[i].r;
    real h = 0;
    real hN = interactionLength;
    real ah, lnah;

    real Fnx, Fny, Fnz, Ftx, Fty, Ftz;
    real Lox, Loy, Loz;

    int isTrue = 0;

    // west wall
    dx = fabs(parts[i].x - (dom->xs + bc.dsW));
    h = dx - ai;
    isTrue = (bc.pW == NEUMANN); // collision force applied ifTrue
    if(h < hN && h > 0) {
      Un = parts[i].u - bc.uWD;
      if(h < eps*parts[i].r) h = eps*parts[i].r;
      ah = ai/h - ai/hN;
      lnah = log(hN/h);

      Utx = 0.;
      Uty = parts[i].v - bc.vWD;
      Utz = parts[i].w - bc.wWD;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = -6.*PI*mu*ai*Un*ah;
      Fny = 0.;
      Fnz = 0.;

      Ftx = 0.;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += 0.;
      Fty += 8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Ftz += -8.*PI*mu*ai*ai*omy*1./10.*lnah;

      Lox = 0.;
      Loy = -8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loz = 8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Lox += 0.;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      parts[i].iFx += isTrue * (Fnx + Ftx);
      parts[i].iFy += isTrue * (Fny + Fty);
      parts[i].iFz += isTrue * (Fnz + Ftz);
      parts[i].iLx += isTrue * Lox;
      parts[i].iLy += isTrue * Loy;
      parts[i].iLz += isTrue * Loz;

      // if no contact, mark with a -1
      parts[i].St = -1.;
    }
    if(h < 0) {
      Un = parts[i].u - bc.uWD;

      // if this is the first time step with contact, set St
      // otherwise, use St that was set at the beginning of the contact
      if(parts[i].St == -1.) {
        parts[i].St = 1./9.*parts[i].rho/rhof*2.*parts[i].r*abs(Un)/nu;
      }
      lnah = 0;
      /** for now, use particle material as wall materal in second V term **/
      real k = 4./3./((1.-parts[i].sigma*parts[i].sigma)/parts[i].E
        + (1.-parts[i].sigma*parts[i].sigma)/parts[i].E)/sqrt(1./ai);

      // estimate alpha according to Tsuji (1991) p. 32 and Joseph (2000) p. 344
      real xcx0 = parts[i].l_rough/(2.*parts[i].r)*100.;
      real e = parts[i].e_dry + (1.+parts[i].e_dry)/parts[i].St*log(xcx0);
      if(e < 0) e = 0;
      real alpha = -2.263*pow(e,0.3948)+2.22;

      real eta = alpha*sqrt(4./3.*PI*ai*ai*ai*parts[i].rho*k*sqrt(-h));

      parts[i].iFx += isTrue * (sqrt(-h*h*h)*k - eta*Un);
    }

    // east wall
    dx = fabs(parts[i].x - (dom->xe - bc.dsE));
    h = dx - ai;
    isTrue = (bc.pE == NEUMANN);
    if(h < hN && h > 0) {
      if(h < eps*parts[i].r) h = eps*parts[i].r;
      ah = ai/h - ai/hN;
      lnah = log(hN/h);

      Un = parts[i].u - bc.uED;
      Utx = 0.;
      Uty = parts[i].v - bc.vED;
      Utz = parts[i].w - bc.wED;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = -6.*PI*mu*ai*Un*ah;
      Fny = 0.;
      Fnz = 0.;

      Ftx = 0.;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += 0.;
      Fty += -8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Ftz += 8.*PI*mu*ai*ai*omy*1./10.*lnah;

      Lox = 0.;
      Loy = 8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loz = -8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Lox += 0.;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      parts[i].iFx += isTrue * (Fnx + Ftx);
      parts[i].iFy += isTrue * (Fny + Fty);
      parts[i].iFz += isTrue * (Fnz + Ftz);
      parts[i].iLx += isTrue * Lox;
      parts[i].iLy += isTrue * Loy;
      parts[i].iLz += isTrue * Loz;

      // if no contact, mark with a -1
      parts[i].St = -1.;
    } 
    if(h < 0) {
      Un = -(parts[i].u - bc.uED);

      // if this is the first time step with contact, set St
      // otherwise, use St that was set at the beginning of the contact
      if(parts[i].St == -1.) {
        parts[i].St = 1./9.*parts[i].rho/rhof*2.*parts[i].r*abs(Un)/nu;
      }
      lnah = 0;
      /** for now, use particle material as wall materal in second V term **/
      real k = 4./3./((1.-parts[i].sigma*parts[i].sigma)/parts[i].E
        + (1.-parts[i].sigma*parts[i].sigma)/parts[i].E)/sqrt(1./ai);

      // estimate alpha according to Tsuji (1991) p. 32 and Joseph (2000) p. 344
      real xcx0 = parts[i].l_rough/(2.*parts[i].r)*100.;
      real e = parts[i].e_dry + (1.+parts[i].e_dry)/parts[i].St*log(xcx0);
      if(e < 0) e = 0;
      real alpha = -2.263*pow(e,0.3948)+2.22;

      real eta = alpha*sqrt(4./3.*PI*ai*ai*ai*parts[i].rho*k*sqrt(-h));

      parts[i].iFx -= isTrue * (sqrt(-h*h*h)*k - eta*Un);
    }

    // south wall
    dy = fabs(parts[i].y - (dom->ys + bc.dsS));
    h = dy - ai;
    isTrue = (bc.pS == NEUMANN);
    if(h < hN && h > 0) {
      if(h < eps*parts[i].r) h = eps*parts[i].r;
      ah = ai/h - ai/hN;
      lnah = log(hN/h);

      Un = parts[i].v - bc.vSD;
      Utx = parts[i].u - bc.uSD;
      Uty = 0.;
      Utz = parts[i].w - bc.wSD;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = -6.*PI*mu*ai*Un*ah;
      Fnz = 0.;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = 0.;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += -8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Fty += 0.;
      Ftz += 8.*PI*mu*ai*ai*omx*1./10.*lnah;

      Lox = 8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loy = 0.;
      Loz = -8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += 0.;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      parts[i].iFx += isTrue * (Fnx + Ftx);
      parts[i].iFy += isTrue * (Fny + Fty);
      parts[i].iFz += isTrue * (Fnz + Ftz);
      parts[i].iLx += isTrue * Lox;
      parts[i].iLy += isTrue * Loy;
      parts[i].iLz += isTrue * Loz;

      // if no contact, mark with a -1
      parts[i].St = -1.;
    }
    if(h < 0) {
      Un = parts[i].v - bc.vSD;

      // if this is the first time step with contact, set St
      // otherwise, use St that was set at the beginning of the contact
      if(parts[i].St == -1.) {
        parts[i].St = 1./9.*parts[i].rho/rhof*2.*parts[i].r*abs(Un)/nu;
      }
      lnah = 0;
      /** for now, use particle material as wall materal in second V term **/
      real k = 4./3./((1.-parts[i].sigma*parts[i].sigma)/parts[i].E
        + (1.-parts[i].sigma*parts[i].sigma)/parts[i].E)/sqrt(1./ai);

      // estimate alpha according to Tsuji (1991) p. 32 and Joseph (2000) p. 344
      real xcx0 = parts[i].l_rough/(2.*parts[i].r)*100.;
      real e = parts[i].e_dry + (1.+parts[i].e_dry)/parts[i].St*log(xcx0);
      if(e < 0) e = 0;
      real alpha = -2.263*pow(e,0.3948)+2.22;

      real eta = alpha*sqrt(4./3.*PI*ai*ai*ai*parts[i].rho*k*sqrt(-h));

      parts[i].iFy += isTrue * (sqrt(-h*h*h)*k - eta*Un);
    }

    // north wall
    dy = fabs(parts[i].y - (dom->ye - bc.dsN));
    h = dy - ai;
    isTrue = (bc.pN == NEUMANN);
    if(h < hN && h > 0) {
      if(h < eps*parts[i].r) h = eps*parts[i].r;
      ah = ai/h - ai/hN;
      lnah = log(hN/h);

      Un = parts[i].v - bc.vND;
      Utx = parts[i].u - bc.uND;
      Uty = 0.;
      Utz = parts[i].w - bc.wND;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = -6.*PI*mu*ai*Un*ah;
      Fnz = 0.;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = 0.;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += 8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Fty += 0.;
      Ftz += -8.*PI*mu*ai*ai*omx*1./10.*lnah;

      Lox = -8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loy = 0.;
      Loz = 8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += 0.;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      parts[i].iFx += isTrue * (Fnx + Ftx);
      parts[i].iFy += isTrue * (Fny + Fty);
      parts[i].iFz += isTrue * (Fnz + Ftz);
      parts[i].iLx += isTrue * Lox;
      parts[i].iLy += isTrue * Loy;
      parts[i].iLz += isTrue * Loz;

      // if no contact, mark with a -1
      parts[i].St = -1.;
    }
    if(h < 0) {
      Un = -(parts[i].v - bc.vND);

      // if this is the first time step with contact, set St
      // otherwise, use St that was set at the beginning of the contact
      if(parts[i].St == -1.) {
        parts[i].St = 1./9.*parts[i].rho/rhof*2.*parts[i].r*abs(Un)/nu;
      }
      lnah = 0;
      /** for now, use particle material as wall materal in second V term **/
      real k = 4./3./((1.-parts[i].sigma*parts[i].sigma)/parts[i].E
        + (1.-parts[i].sigma*parts[i].sigma)/parts[i].E)/sqrt(1./ai);

      // estimate alpha according to Tsuji (1991) p. 32 and Joseph (2000) p. 344
      real xcx0 = parts[i].l_rough/(2.*parts[i].r)*100.;
      real e = parts[i].e_dry + (1.+parts[i].e_dry)/parts[i].St*log(xcx0);
      if(e < 0) e = 0;
      real alpha = -2.263*pow(e,0.3948)+2.22;

      real eta = alpha*sqrt(4./3.*PI*ai*ai*ai*parts[i].rho*k*sqrt(-h));

      parts[i].iFy -= isTrue * (sqrt(-h*h*h)*k - eta*Un);
    }

    // bottom wall
    dz = fabs(parts[i].z - (dom->zs + bc.dsB));
    h = dz - ai;
    isTrue = (bc.pB == NEUMANN);
    if(h < hN && h > 0) {
      if(h < eps*parts[i].r) h = eps*parts[i].r;
      ah = ai/h - ai/hN;
      lnah = log(hN/h);

      Un = parts[i].w - bc.wBD;
      Utx = parts[i].u - bc.uBD;
      Uty = parts[i].v - bc.vBD;
      Utz = 0.;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = 0.;
      Fnz = -6.*PI*mu*ai*Un*ah;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = 0.;
      Ftx += 8.*PI*mu*ai*ai*omy*1./10.*lnah;
      Fty += -8.*PI*mu*ai*ai*omx*1./10.*lnah;
      Ftz += 0.;

      Lox = -8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Loy = 8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Loz = 0.;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += 0.;

      parts[i].iFx += isTrue * (Fnx + Ftx);
      parts[i].iFy += isTrue * (Fny + Fty);
      parts[i].iFz += isTrue * (Fnz + Ftz);
      parts[i].iLx += isTrue * Lox;
      parts[i].iLy += isTrue * Loy;
      parts[i].iLz += isTrue * Loz;

      // if no contact, mark with a -1
      parts[i].St = -1.;
    }
    if(h < 0) {
      Un = parts[i].w - bc.wBD;

      // if this is the first time step with contact, set St
      // otherwise, use St that was set at the beginning of the contact
      if(parts[i].St == -1.) {
        parts[i].St = 1./9.*parts[i].rho/rhof*2.*parts[i].r*abs(Un)/nu;
      }
      lnah = 0;
      // for now, use particle material as wall materal in second V term //
      real k = 4./3./((1.-parts[i].sigma*parts[i].sigma)/parts[i].E
        + (1.-parts[i].sigma*parts[i].sigma)/parts[i].E)/sqrt(1./ai);

      // estimate alpha according to Tsuji (1991) p. 32 and Joseph (2000) p. 344
      real xcx0 = parts[i].l_rough/(2.*parts[i].r)*100.;
      real e = parts[i].e_dry + (1.+parts[i].e_dry)/parts[i].St*log(xcx0);
      if(e < 0) e = 0;
      real alpha = -2.263*pow(e,0.3948)+2.22;

      real eta = alpha*sqrt(4./3.*PI*ai*ai*ai*parts[i].rho*k*sqrt(-h));

      parts[i].iFz += isTrue * (sqrt(-h*h*h)*k - eta*Un);
    }

    // top wall
    dz = fabs(parts[i].z - (dom->ze - bc.dsT));
    h = dz - ai;
    isTrue = (bc.pT == NEUMANN);
    if(h < hN && h > 0) {
      if(h < eps*parts[i].r) h = eps*parts[i].r;
      ah = ai/h - ai/hN;
      lnah = log(hN/h);

      Un = parts[i].w - bc.wTD;
      Utx = parts[i].u - bc.uTD;
      Uty = parts[i].v - bc.vTD;
      Utz = 0.;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = 0.;
      Fnz = -6.*PI*mu*ai*Un*ah;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = 0.;
      Ftx += -8.*PI*mu*ai*ai*omy*1./10.*lnah;
      Fty += 8.*PI*mu*ai*ai*omx*1./10.*lnah;
      Ftz += 0.;

      Lox = 8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Loy = -8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Loz = 0.;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += 0.;

      parts[i].iFx += isTrue * (Fnx + Ftx);
      parts[i].iFy += isTrue * (Fny + Fty);
      parts[i].iFz += isTrue * (Fnz + Ftz);
      parts[i].iLx += isTrue * Lox;
      parts[i].iLy += isTrue * Loy;
      parts[i].iLz += isTrue * Loz;

      // if no contact, mark with a -1
      parts[i].St = -1.;
    }
    if(h < 0) {
      Un = -(parts[i].w - bc.wTD);

      // if this is the first time step with contact, set St
      // otherwise, use St that was set at the beginning of the contact
      if(parts[i].St == -1.) {
        parts[i].St = 1./9.*parts[i].rho/rhof*2.*parts[i].r*abs(Un)/nu;
      }
      lnah = 0;
      // for now, use particle material as wall materal in second V term //
      real k = 4./3./((1.-parts[i].sigma*parts[i].sigma)/parts[i].E
        + (1.-parts[i].sigma*parts[i].sigma)/parts[i].E)/sqrt(1./ai);

      // estimate alpha according to Tsuji (1991) p. 32 and Joseph (2000) p. 344
      real xcx0 = parts[i].l_rough/(2.*parts[i].r)*100.;
      real e = parts[i].e_dry + (1.+parts[i].e_dry)/parts[i].St*log(xcx0);
      if(e < 0) e = 0;
      real alpha = -2.263*pow(e,0.3948)+2.22;

      real eta = alpha*sqrt(4./3.*PI*ai*ai*ai*parts[i].rho*k*sqrt(-h));

      parts[i].iFz -= isTrue * (sqrt(-h*h*h)*k - eta*Un);
    }
  }
}

__global__ void spring_parts(part_struct *parts, int nparts)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if(i < nparts && parts[i].spring_k > 0.) {
    real nx = parts[i].x-parts[i].spring_x;
    real ny = parts[i].y-parts[i].spring_y;
    real nz = parts[i].z-parts[i].spring_z;
    real n = sqrt(nx*nx+ny*ny+nz*nz);
    real nhatx = nx / n;
    real nhaty = ny / n;
    real nhatz = nz / n;
    real lx = parts[i].spring_l*nhatx;
    real ly = parts[i].spring_l*nhaty;
    real lz = parts[i].spring_l*nhatz;
    real l = sqrt(lx*lx+ly*ly+lz*lz);

    real dx = parts[i].x-parts[i].spring_x-lx;
    real dy = parts[i].y-parts[i].spring_y-ly;
    real dz = parts[i].z-parts[i].spring_z-lz;

    parts[i].kFx = - parts[i].spring_k * dx;
    parts[i].kFy = - parts[i].spring_k * dy;
    parts[i].kFz = - parts[i].spring_k * dz;
  }
}

__global__ void yank_u_WE(real *u, dom_struct *dom, real *plane, real xpos,
  real vel)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  real ddx = 1. / dom->dx;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int i = floor((xpos - dom->xs) * ddx) + DOM_BUF;
    if(i < dom->Gfx.is) i += dom->Gfx.inb;
    if(i > dom->Gfx.ie-1) i -= dom->Gfx.inb;
    real xx = (i-DOM_BUF) * dom->dx + dom->xs;
    int W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    int E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    real dudx = (u[E] - u[W]) * ddx;

    plane[tj + tk*dom->Gfx.jnb] = u[W] + dudx * (xpos - xx) + vel;
  }
}

__global__ void yank_v_WE(real *v, dom_struct *dom, real *plane, real xpos,
  real vel)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  real ddx = 1. / dom->dx;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int i = floor((xpos - dom->xs) * ddx - 0.5) + DOM_BUF;
    if(i < dom->Gfy.is) i += dom->Gfy.inb;
    if(i > dom->Gfy.ie-1) i -= dom->Gfy.inb;
    real xx = (i-DOM_BUF+0.5) * dom->dx + dom->xs;
    int W = i + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    int E = (i+1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    real dvdx = (v[E] - v[W]) * ddx;

    plane[tj + tk*dom->Gfy.jnb] = v[W] + dvdx * (xpos - xx);
  }
}

__global__ void yank_w_WE(real *w, dom_struct *dom, real *plane, real xpos,
  real vel)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  real ddx = 1. / dom->dx;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int i = floor((xpos - dom->xs) * ddx - 0.5) + DOM_BUF;
    if(i < dom->Gfz.is) i += dom->Gfz.inb;
    if(i > dom->Gfz.ie-1) i -= dom->Gfz.inb;
    real xx = (i-DOM_BUF + 0.5) * dom->dx + dom->xs;
    int W = i + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    int E = (i+1) + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    real dwdx = (w[E] - w[W]) * ddx;

    plane[tj + tk*dom->Gfz.jnb] = w[W] + dwdx * (xpos - xx);
  }
}

__global__ void yank_u_SN(real *u, dom_struct *dom, real *plane, real ypos,
  real vel)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  real ddy = 1. / dom->dy;

  if((tk < dom->Gfx._inb) && (ti < dom->Gfx._inb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int j = floor((ypos - dom->ys) * ddy - 0.5) + DOM_BUF;
    if(j < dom->Gfx.js) j += dom->Gfx.jnb;
    if(j > dom->Gfx.je-1) j -= dom->Gfx.jnb;
    real yy = (j-DOM_BUF + 0.5) * dom->dy + dom->ys;
    int S = ti + j*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    int N = ti + (j+1)*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    real dudy = (u[N] - u[S]) * ddy;

    plane[tk + ti*dom->Gfx.knb] = u[S] + dudy * (ypos - yy);
  }
}

__global__ void yank_v_SN(real *v, dom_struct *dom, real *plane, real ypos,
  real vel)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  real ddy = 1. / dom->dy;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int j = floor((ypos - dom->ys) * ddy) + DOM_BUF;
    if(j < dom->Gfy.js) j += dom->Gfy.jnb;
    if(j > dom->Gfy.je-1) j -= dom->Gfy.jnb;
    real yy = (j-DOM_BUF) * dom->dy + dom->ys;
    int S = ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    int N = ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    real dvdy = (v[N] - v[S]) * ddy;

    plane[tk + ti*dom->Gfy.knb] = v[S] + dvdy * (ypos - yy) + vel;
  }
}

__global__ void yank_w_SN(real *w, dom_struct *dom, real *plane, real ypos,
  real vel)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  real ddy = 1. / dom->dy;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int j = floor((ypos - dom->ys) * ddy - 0.5) + DOM_BUF;
    if(j < dom->Gfz.js) j += dom->Gfz.jnb;
    if(j > dom->Gfz.je-1) j -= dom->Gfz.jnb;
    real yy = (j-DOM_BUF + 0.5) * dom->dy + dom->ys;
    int S = ti + j*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    int N = ti + (j+1)*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    real dwdy = (w[N] - w[S]) * ddy;

    plane[tk + ti*dom->Gfz.knb] = w[S] + dwdy * (ypos - yy);
  }
}

__global__ void yank_u_BT(real *u, dom_struct *dom, real *plane, real zpos,
  real vel)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  real ddz = 1. / dom->dz;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int k = floor((zpos - dom->zs) * ddz - 0.5) + DOM_BUF;
    if(k < dom->Gfx.ks) k += dom->Gfx.knb;
    if(k > dom->Gfx.ke-1) k -= dom->Gfx.knb;
    real zz = (k-DOM_BUF + 0.5) * dom->dz + dom->zs;
    int B = ti + tj*dom->Gfx.s1b + k*dom->Gfx.s2b;
    int T = ti + tj*dom->Gfx.s1b + (k+1)*dom->Gfx.s2b;
    real dudz = (u[T] - u[B]) * ddz;

    plane[ti + tj*dom->Gfx.inb] = u[B] + dudz * (zpos - zz);
  }
}

__global__ void yank_v_BT(real *v, dom_struct *dom, real *plane, real zpos,
  real vel)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  real ddz = 1. / dom->dz;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int k = floor((zpos - dom->zs) * ddz - 0.5) + DOM_BUF;
    if(k < dom->Gfy.ks) k += dom->Gfy.knb;
    if(k > dom->Gfy.ke-1) k -= dom->Gfy.knb;
    real zz = (k-DOM_BUF + 0.5) * dom->dz + dom->zs;
    int B = ti + tj*dom->Gfy.s1b + k*dom->Gfy.s2b;
    int T = ti + tj*dom->Gfy.s1b + (k+1)*dom->Gfy.s2b;
    real dvdz = (v[T] - v[B]) * ddz;

    plane[ti + tj*dom->Gfy.inb] = v[B] + dvdz * (zpos - zz);
  }
}

__global__ void yank_w_BT(real *w, dom_struct *dom, real *plane, real zpos,
  real vel)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  real ddz = 1. / dom->dx;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int k = floor((zpos - dom->zs) * ddz) + DOM_BUF;
    if(k < dom->Gfz.ks) k += dom->Gfz.knb;
    if(k > dom->Gfz.ke-1) k -= dom->Gfz.knb;
    real zz = (k-DOM_BUF) * dom->dz + dom->zs;
    int B = ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b;
    int T = ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b;
    real dwdz = (w[T] - w[B]) * ddz;

    plane[ti + tj*dom->Gfz.inb] = w[B] + dwdz * (zpos - zz) + vel;
  }
}

__global__ void colocate_Gfx(real *u, real *u_co, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tj < dom->Gfx.jnb-1) && (tk < dom->Gfx.knb-1)) {
    for(int i = dom->Gfx.is; i < dom->Gfx.ie-1; i++) {
      u_co[(i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2] =
        0.5 * (u[i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
        + u[(i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]);
    }
  }
}

__global__ void colocate_Gfy(real *v, real *v_co, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int ti = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tk < dom->Gfy.knb-1) && (ti < dom->Gfy.inb-1)) {
    for(int j = dom->Gfy.js; j < dom->Gfy.je-1; j++) {
      v_co[(ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2] =
        0.5 * (v[ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b]
        + v[ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b]);
    }
  }
}

__global__ void colocate_Gfz(real *w, real *w_co, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tj = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((ti < dom->Gfz.inb-1) && (tj < dom->Gfz.jnb-1)) {
    for(int k = dom->Gfz.ks; k < dom->Gfz.ke-1; k++) {
      w_co[(ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2] =
        0.5 * (w[ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b]
        + w[ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b]);
    }
  }
}

__global__ void energy_multiply(real *u_co, real *v_co, real *w_co, real *co,
  dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int C;  // memory location

  if((tj < dom->Gcc.jn) && (tk < dom->Gcc.kn)) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      co[C] = u_co[C]*u_co[C] + v_co[C]*v_co[C] + w_co[C]*w_co[C];
    }
  }
}

__device__ real ab_int(real dt0, real dt, real f0, real df0, real df)
{
  real DT = dt/dt0;
  if(dt0 < 0) 
    return f0 + df*dt;
  else
    return f0 + ((1+0.5*DT)*df - 0.5*DT*df0)*dt;
}

__global__ void internal_u(real *u, part_struct *parts, dom_struct *dom,
  int *flag_u, int *phase)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tj < dom->Gfx._je && tk < dom->Gfx._ke) {
    for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
      int C = i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b;
      int W = (i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int E = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;

      int pw = phase[W];
      int pe = phase[E];
      int f = flag_u[C];

      int p = (pw > -1 && pe > -1) * phase[E];

      real rx = (i - DOM_BUF) * dom->dx + dom->xs - parts[p].x;
      if(rx < dom->xs - 0.5*dom->dx) rx += dom->xl;
      if(rx > dom->xe + 0.5*dom->dx) rx -= dom->xl;
      real ry = (tj - 0.5) * dom->dy + dom->ys - parts[p].y;
      if(ry < dom->ys - 0.5*dom->dy) ry += dom->yl;
      if(ry > dom->ye + 0.5*dom->dy) ry -= dom->yl;
      real rz = (tk - 0.5) * dom->dz + dom->zs - parts[p].z;
      if(rz < dom->zs - 0.5*dom->dz) rz += dom->zl;
      if(rz > dom->ze + 0.5*dom->dz) rz -= dom->zl;

      real ocrossr_x = parts[p].oy*rz - parts[p].oz*ry;

      u[C] = (pw == -1 || pe == -1 || f == -1) * u[C]
        + (pw > -1 && pe > -1 && f != -1) * (ocrossr_x + parts[p].u);
    }
  }
}

__global__ void internal_v(real *v, part_struct *parts, dom_struct *dom,
  int *flag_v, int *phase)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tk < dom->Gfy._ke && ti < dom->Gfy._ie) {
    for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
      int C = ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b;
      int S = ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      int N = ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b;

      int ps = phase[S];
      int pn = phase[N];
      int f = flag_v[C];

      int p = (ps > -1 && pn > -1) * phase[N];

      real rx = (ti - 0.5) * dom->dx + dom->xs - parts[p].x;
      if(rx < dom->xs - 0.5*dom->dx) rx += dom->xl;
      if(rx > dom->xe + 0.5*dom->dx) rx -= dom->xl;
      real ry = (j - DOM_BUF) * dom->dy + dom->ys - parts[p].y;
      if(ry < dom->ys - 0.5*dom->dy) ry += dom->yl;
      if(ry > dom->ye + 0.5*dom->dy) ry -= dom->yl;
      real rz = (tk - 0.5) * dom->dz + dom->zs - parts[p].z;
      if(rz < dom->zs - 0.5*dom->dz) rz += dom->zl;
      if(rz > dom->ze + 0.5*dom->dz) rz -= dom->zl;

      real ocrossr_y = parts[p].oz*rx - parts[p].ox*rz;

      v[C] = (ps == -1 || pn == -1 || f == -1) * v[C]
        + (ps > -1 && pn > -1 && f != -1) * (ocrossr_y + parts[p].v);
    }
  }
}

__global__ void internal_w(real *w, part_struct *parts, dom_struct *dom,
  int *flag_w, int *phase)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(ti < dom->Gfz._ie && tj < dom->Gfz._je) {
    for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
      int C = ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b;
      int B = ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b;
      int T = ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b;

      int pb = phase[B];
      int pt = phase[T];
      int f = flag_w[C];

      int p = (pb > -1 && pt > -1) * phase[T];

      real rx = (ti - 0.5) * dom->dx + dom->xs - parts[p].x;
      if(rx < dom->xs - 0.5*dom->dx) rx += dom->xl;
      if(rx > dom->xe + 0.5*dom->dx) rx -= dom->xl;
      real ry = (tj - 0.5) * dom->dy + dom->ys - parts[p].y;
      if(ry < dom->ys - 0.5*dom->dy) ry += dom->yl;
      if(ry > dom->ye + 0.5*dom->dy) ry -= dom->yl;
      real rz = (k - DOM_BUF) * dom->dz + dom->zs - parts[p].z;
      if(rz < dom->zs - 0.5*dom->dz) rz += dom->zl;
      if(rz > dom->ze + 0.5*dom->dz) rz -= dom->zl;

      real ocrossr_z = parts[p].ox*ry - parts[p].oy*rx;

      w[C] = (pb == -1 || pt == -1 || f == -1) * w[C]
        + (pb > -1 && pt > -1 && f != -1) * (ocrossr_z + parts[p].w);
    }
  }
}
