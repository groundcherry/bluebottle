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

#include "cuda_particle.h"

__global__ void reset_phase(int *phase, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb)) {
    for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) {
      phase[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b] = -1;
    }
  }
}

__global__ void reset_phase_shell(int *phase_shell, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb)) {
    for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) {
      phase_shell[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b] = 1;
    }
  }
}

__global__ void reset_flag_u(int *flag_u, dom_struct *dom)
{
  int i;    // iterator
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    for(i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
      flag_u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 1;
    }
  }
}

__global__ void reset_flag_v(int *flag_v, dom_struct *dom)
{
  int j;    // iterator
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tk < dom->Gfy._knb) && (ti < dom->Gfy._inb)) {
    for(j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
      flag_v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 1;
    }
  }
}

__global__ void reset_flag_w(int *flag_w, dom_struct *dom)
{
  int k;    // iterator
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    for(k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
      flag_w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 1;
    }
  }
}

__global__ void build_phase(int p, part_struct *parts, int *phase,
  dom_struct *dom, real X, real Y, real Z,
  int is, int ie, int js, int je, int ks, int ke)
{
  real xx, yy, zz;  // distance from cell center to particle center along
                    // Cartesian basis
  real d;           // distance form cell center to particle center
  int C;            // current cell
  int cutoff;       // cage cutoff constant

  // update phase (use center of cell containing particle center)
  int ti = blockDim.x*blockIdx.x + threadIdx.x + is;
  int tj = blockDim.y*blockIdx.y + threadIdx.y + js;
  int tk = blockDim.z*blockIdx.z + threadIdx.z + ks;
  if((ti < ie) && (tj < je) && (tk < ke)) {
    xx = (ti-0.5)*dom->dx - (X-dom->xs);
    yy = (tj-0.5)*dom->dy - (Y-dom->ys);
    zz = (tk-0.5)*dom->dz - (Z-dom->zs);
    d = sqrt(xx * xx + yy * yy + zz * zz);

    C = ti + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;

    cutoff = (floor(d / (parts[p].r))) < 1;

    phase[C] = cutoff*p + (1 - cutoff)*phase[C];
  }
}

__global__ void build_cage(int p, part_struct *parts, int *phase,
  int *phase_shell, dom_struct *dom, real Y, real Z,
  int js, int je, int ks, int ke)
{
  real xx, yy, zz;  // distance from cell center to particle center along
                    // Cartesian basis
  real d;           // distance form cell center to particle center
  int cutoff;       // cage cutoff constant
  real X;           // particle center location

  // update phase (use center of cell containing particle center)
  int tj = blockDim.x*blockIdx.x + threadIdx.x + js;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + ks;

  if((tj < je) && (tk < ke)) {
    X = parts[p].x;
    if(parts[p].x < dom->xs) X = parts[p].x + dom->xl;
    for(int i = parts[p].cage.is; i < parts[p].cage.ibs; i++) {
      xx = (i-0.5)*dom->dx - (X-dom->xs);
      yy = (tj-0.5)*dom->dy - (Y-dom->ys);
      zz = (tk-0.5)*dom->dz - (Z-dom->zs);
      d = sqrt(xx * xx + yy * yy + zz * zz);

      cutoff =
        (1 - floor(d / (1.0*parts[p].r
        - 0.50*(dom->dx + dom->dy + dom->dz)/3.)));

      if((cutoff * (p+1) - 1) > -1)
        phase[i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b] = cutoff * (p + 1) - 1;

      /*cutoff = (cutoff>0) &&
        (ceil(d / (1.0*parts[p].r - 2.*(dom->dx + dom->dy + dom->dz)/3.))-1);

      if((cutoff*(p+1)-1) > -1)
        phase_shell[i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b] = 1-cutoff;
*/
    }
    X = parts[p].x;
    if(parts[p].x > dom->xe) X = parts[p].x - dom->xl;
    for(int i = parts[p].cage.ibe; i < parts[p].cage.ie; i++) {
      xx = (i-0.5)*dom->dx - (X-dom->xs);
      yy = (tj-0.5)*dom->dy - (Y-dom->ys);
      zz = (tk-0.5)*dom->dz - (Z-dom->zs);
      d = sqrt(xx * xx + yy * yy + zz * zz);

      cutoff =
        (1 - floor(d / (1.0*parts[p].r
        - 0.50*(dom->dx + dom->dy + dom->dz)/3.)));

      if((cutoff * (p+1) - 1) > -1)
        phase[i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b] = cutoff * (p + 1) - 1;

      /*cutoff = (cutoff>0) &&
        (ceil(d / (1.0*parts[p].r - 2.*(dom->dx + dom->dy + dom->dz)/3.))-1);

      if((cutoff*(p+1)-1) > -1)
        phase_shell[i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b] = 1-cutoff;
*/
    }

  }
}

__global__ void cage_phases_periodic_x(int *phase_type, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gcc.jnb && tk < dom->Gcc.knb) {
    phase_type[dom->Gcc.isb + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b]
      = phase_type[(dom->Gcc.ie-1) + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b];
    phase_type[(dom->Gcc.ieb-1) + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b]
      = phase_type[dom->Gcc.is + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b];
  }
}

__global__ void cage_phases_periodic_y(int *phase_type, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gcc.knb && ti < dom->Gcc.inb) {
    phase_type[ti + dom->Gcc.jsb*dom->Gcc.s1b + tk*dom->Gcc.s2b]
      = phase_type[ti + (dom->Gcc.je-1)*dom->Gcc.s1b + tk*dom->Gcc.s2b];
    phase_type[ti + (dom->Gcc.jeb-1)*dom->Gcc.s1b + tk*dom->Gcc.s2b]
      = phase_type[ti + dom->Gcc.js*dom->Gcc.s1b + tk*dom->Gcc.s2b];
  }
}

__global__ void cage_phases_periodic_z(int *phase_type, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gcc.inb && tj < dom->Gcc.jnb) {
    phase_type[ti + tj*dom->Gcc.s1b + dom->Gcc.ksb*dom->Gcc.s2b]
      = phase_type[ti + tj*dom->Gcc.s1b + (dom->Gcc.ke-1)*dom->Gcc.s2b];
    phase_type[ti + tj*dom->Gcc.s1b + (dom->Gcc.keb-1)*dom->Gcc.s2b]
      = phase_type[ti + tj*dom->Gcc.s1b + dom->Gcc.ks*dom->Gcc.s2b];
  }
}

__global__ void phase_shell_x(part_struct *parts,
  dom_struct *dom, int *phase, int *phase_shell)
{
  int i;    // iterator

  int W, E; // flag locations

  int tj = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tj < dom->Gcc.je) && (tk < dom->Gcc.ke)) {
    for(i = dom->Gcc.is; i <= dom->Gcc.ie; i++) {
      W = (i-1) + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;
      E = i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;
      
      // if phase changes from fluid to solid
      phase_shell[E] = phase_shell[E]*(1 - (phase[W] < 0 && phase[E] > -1));
      // if phase changes from solid to fluid
      phase_shell[W] = phase_shell[W]*(1 - (phase[W] > -1 && phase[E] < 0));

    }
  }
}

__global__ void phase_shell_y(part_struct *parts,
  dom_struct *dom, int *phase, int *phase_shell)
{
  int j;    // iterator

  int N, S; // flag locations

  int tk = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int ti = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tk < dom->Gcc.ke) && (ti < dom->Gcc.ie)) {
    for(j = dom->Gcc.js; j <= dom->Gcc.je; j++) {
      S = ti + (j-1)*dom->Gcc.s1b + tk*dom->Gcc.s2b;
      N = ti + j*dom->Gcc.s1b + tk*dom->Gcc.s2b;

      // if phase changes from fluid to solid
      phase_shell[N] = phase_shell[N]*(1 - (phase[S] < 0 && phase[N] > -1));
      // if phase changes from solid to fluid
      phase_shell[S] = phase_shell[S]*(1 - (phase[S] > -1 && phase[N] < 0));

    }
  }
}

__global__ void phase_shell_z(part_struct *parts,
  dom_struct *dom, int *phase, int *phase_shell)
{
  int k;    // iterator

  int B, T; // flag locations

  int ti = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tj = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((ti < dom->Gcc.ie) && (tj < dom->Gcc.je)) {
    for(k = dom->Gcc.ks; k <= dom->Gcc.ke; k++) {
      B = ti + tj*dom->Gcc.s1b + (k-1)*dom->Gcc.s2b;
      T = ti + tj*dom->Gcc.s1b + k*dom->Gcc.s2b;
      // if phase changes from fluid to solid
      phase_shell[T] = phase_shell[T]*(1 - (phase[B] < 0 && phase[T] > -1));
      // if phase changes from solid to fluid
      phase_shell[B] = phase_shell[B]*(1 - (phase[B] > -1 && phase[T] < 0));
    }
  }
}

__global__ void cage_flag_u_periodic_x(int *flag_u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfx.jnb && tk < dom->Gfx.knb) {
    flag_u[dom->Gfx.isb + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[(dom->Gfx.ie-2) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b];
    flag_u[(dom->Gfx.ieb-1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[(dom->Gfx.is+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_y(int *flag_u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfx.knb && ti < dom->Gfx.inb) {
    flag_u[ti + dom->Gfx.jsb*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[ti + (dom->Gfx.je-1)*dom->Gfx.s1b + tk*dom->Gfx.s2b];
    flag_u[ti + (dom->Gfx.jeb-1)*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[ti + dom->Gfx.js*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_z(int *flag_u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfx.inb && tj < dom->Gfx.jnb) {
    flag_u[ti + tj*dom->Gfx.s1b + dom->Gfx.ksb*dom->Gfx.s2b]
      = flag_u[ti + tj*dom->Gfx.s1b + (dom->Gfx.ke-1)*dom->Gfx.s2b];
    flag_u[ti + tj*dom->Gfx.s1b + (dom->Gfx.keb-1)*dom->Gfx.s2b]
      = flag_u[ti + tj*dom->Gfx.s1b + dom->Gfx.ks*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_v_periodic_x(int *flag_v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfy.jnb && tk < dom->Gfy.knb) {
    flag_v[dom->Gfy.isb + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[(dom->Gfy.ie-1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b];
    flag_v[(dom->Gfy.ieb-1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[dom->Gfy.is + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_y(int *flag_v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfy.knb && ti < dom->Gfy.inb) {
    flag_v[ti + dom->Gfy.jsb*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[ti + (dom->Gfy.je-2)*dom->Gfy.s1b + tk*dom->Gfy.s2b];
    flag_v[ti + (dom->Gfy.jeb-1)*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[ti + (dom->Gfy.js+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_z(int *flag_v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfy.inb && tj < dom->Gfy.jnb) {
    flag_v[ti + tj*dom->Gfy.s1b + dom->Gfy.ksb*dom->Gfy.s2b]
      = flag_v[ti + tj*dom->Gfy.s1b + (dom->Gfy.ke-1)*dom->Gfy.s2b];
    flag_v[ti + tj*dom->Gfy.s1b + (dom->Gfy.keb-1)*dom->Gfy.s2b]
      = flag_v[ti + tj*dom->Gfy.s1b + dom->Gfy.ks*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_w_periodic_x(int *flag_w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfz.jnb && tk < dom->Gfz.knb) {
    flag_w[dom->Gfz.isb + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[(dom->Gfz.ie-1)+ tj*dom->Gfz.s1b + tk*dom->Gfz.s2b];
    flag_w[(dom->Gfz.ieb-1) + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[dom->Gfz.is + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_y(int *flag_w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfz.knb && ti < dom->Gfz.inb) {
    flag_w[ti + dom->Gfz.jsb*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[ti + (dom->Gfz.je-1)*dom->Gfz.s1b + tk*dom->Gfz.s2b];
    flag_w[ti + (dom->Gfz.jeb-1)*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[ti + dom->Gfz.js*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_z(int *flag_w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfz.inb && tj < dom->Gfz.jnb) {
    flag_w[ti + tj*dom->Gfz.s1b + dom->Gfz.ksb*dom->Gfz.s2b]
      = flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.ke-2)*dom->Gfz.s2b];
    flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.keb-1)*dom->Gfz.s2b]
      = flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.ks+1)*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_u(int *flag_u, part_struct *parts, dom_struct *dom,
  int *phase, int *phase_shell)
{
  int i;    // iterator

  int W, E; // flag locations

  int tj = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tj < dom->Gcc.je) && (tk < dom->Gcc.ke)) {
    for(i = dom->Gcc.is; i <= dom->Gcc.ie; i++) {
      W = (i-1) + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;
      E = i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;

#ifdef STEPS
      flag_u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 
        1 - ((phase[W] < 0 && phase[E] > -1) 
            || (phase[W] > -1 && phase[E] < 0)
            || ((phase_shell[W] < 1 && phase_shell[E] < 1)));
#else
      flag_u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 
        1 - 2*((phase[W] < 0 && phase[E] > -1) 
            || (phase[W] > -1 && phase[E] < 0)
            || ((phase_shell[W] < 1 && phase_shell[E] < 1)));
#endif
    }
  }
}

__global__ void cage_flag_v(int *flag_v, part_struct *parts, dom_struct *dom,
  int *phase, int *phase_shell)
{
  int j;    // iterator
  int S, N; // flag locations

  int tk = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int ti = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tk < dom->Gcc.ke) && (ti < dom->Gcc.ie)) {
    for(j = dom->Gcc.js; j <= dom->Gcc.je; j++) {
      S = ti + (j-1)*dom->Gcc.s1b + tk*dom->Gcc.s2b;
      N = ti + j*dom->Gcc.s1b + tk*dom->Gcc.s2b;

#ifdef STEPS
      flag_v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 
        1 - ((phase[S] < 0 && phase[N] > -1)
            || (phase[S] > -1 && phase[N] < 0)
            || ((phase_shell[S] < 1 && phase_shell[N] < 1))); 
#else
      flag_v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 
        1 - 2*((phase[S] < 0 && phase[N] > -1)
            || (phase[S] > -1 && phase[N] < 0)
            || ((phase_shell[S] < 1 && phase_shell[N] < 1))); 
#endif
    }
  }
}

__global__ void cage_flag_w(int *flag_w, part_struct *parts, dom_struct *dom, 
  int *phase, int *phase_shell)
{
  int k;    // iterator
  int B, T; // flag locations

  int ti = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tj = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((ti < dom->Gcc.ie) && (tj < dom->Gcc.je)) {
    for(k = dom->Gcc.ks; k <= dom->Gcc.ke; k++) {
      B = ti + tj*dom->Gcc.s1b + (k-1)*dom->Gcc.s2b;
      T = ti + tj*dom->Gcc.s1b + k*dom->Gcc.s2b;

#ifdef STEPS
      flag_w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 
        1 - ((phase[B] < 0 && phase[T] > -1)
            || (phase[B] > -1 && phase[T] < 0)
            || ((phase_shell[B] < 1 && phase_shell[T] < 1)));
#else
      flag_w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 
        1 - 2*((phase[B] < 0 && phase[T] > -1)
            || (phase[B] > -1 && phase[T] < 0)
            || ((phase_shell[B] < 1 && phase_shell[T] < 1)));
#endif

    }
  }
}

__global__ void flag_external_u(int *flag_u, dom_struct *dom)
{
  
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;
  
  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    flag_u[dom->Gfx._is + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0;
    flag_u[dom->Gfx._ie-1 + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0;
  }
}
      
__global__ void flag_external_v(int *flag_v, dom_struct *dom)
{
  
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;
  
  if((tk < dom->Gfy._knb) && (ti < dom->Gfy._inb)) {
    flag_v[ti + dom->Gfy._js*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0;
    flag_v[ti + (dom->Gfy._je-1)*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0;
  }
}

__global__ void flag_external_w(int *flag_w, dom_struct *dom)
{
  
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;
  
  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    flag_w[ti + tj*dom->Gfz._s1b + dom->Gfz._ks*dom->Gfz._s2b] = 0;
    flag_w[ti + tj*dom->Gfz._s1b + (dom->Gfz._ke-1)*dom->Gfz._s2b] = 0;
  }
}

__global__ void part_BC_u(real *u, int *phase, int *flag_u,
  part_struct *parts, dom_struct *dom,
  real nu, int stride,
  real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x + dom->Gfx._jsb;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + dom->Gfx._ksb;
  int C, CW, CE;
  real x, y, z;         // velocity node location Cartesian
  real X, Y, Z;         // particle position
  real r, theta, phi;   // velocity node location spherical
  real Ux, Uy, Uz;      // Cartesian-directed pressure gradients
  int P, PP, PW, PE;    // particle number
  real a;               // particle radius
  int order;            // particle order
  real oy, oz;          // particle angular velocity
  real oydot, ozdot;    // particle angular acceleration
  real uu;              // particle velocity

  if(tj < dom->Gfx._jeb && tk < dom->Gfx._keb) {
    for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
      C = i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b;
      CW = (i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      CE = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      // get particle number
      PW = phase[CW];
      PE = phase[CE];
      if(PW > -1) {
        P = PW;
        PP = PW;
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        oy = parts[P].oy;
        oz = parts[P].oz;
        oydot = parts[P].oydot;
        ozdot = parts[P].ozdot;
        uu = parts[P].u;
      } else if(PE > -1) {
        P = PE;
        PP = PE;
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        oy = parts[P].oy;
        oz = parts[P].oz;
        oydot = parts[P].oydot;
        ozdot = parts[P].ozdot;
        uu = parts[P].u;
      } else {
        P = 0;
        PP = -1;
        a = (dom->dx + dom->dy + dom->dz) / 3.;
        X = (i-DOM_BUF) * dom->dx + dom->xs + a;
        Y = (tj-0.5) * dom->dy + dom->ys + a;
        Z = (tk-0.5) * dom->dz + dom->zs + a;
        order = 0;
        oy = 0;
        oz = 0;
        oydot = 0;
        ozdot = 0;
        uu = 0;
      }
      x = (i-DOM_BUF) * dom->dx + dom->xs - X;
      if(x < dom->xs - 0.5*dom->dx) x += dom->xl;
      if(x > dom->xe + 0.5*dom->dx) x -= dom->xl;
      y = (tj-0.5) * dom->dy + dom->ys - Y;
      if(y < dom->ys - 0.5*dom->dy) y += dom->yl;
      if(y > dom->ye + 0.5*dom->dy) y -= dom->yl;
      z = (tk-0.5) * dom->dz + dom->zs - Z;
      if(z < dom->zs - 0.5*dom->dz) z += dom->zl;
      if(z > dom->ze + 0.5*dom->dz) z -= dom->zl;
      xyz2rtp(x, y, z, &r, &theta, &phi);

      // calculate analytic solution
#ifdef STEPS
      int check = (flag_u[C] < 1) && (PP > -1);
      u[C] = - (check - 1) * u[C];
#else
      lamb_vel(order, a, r, theta, phi,
        nu, pnm_re, pnm_im, phinm_re, phinm_im,
        chinm_re, chinm_im,
        P, stride, &Ux, &Uy, &Uz);


      real ocrossr_x = oy*z - oz*y;
      real odotcrossr_x = oydot*z - ozdot*y;
      Ux += uu + ocrossr_x;
      Ux += 0.1/nu *(r*r*r*r*r-a*a*a*a*a)/(r*r*r) * odotcrossr_x;
      // boolean check if this is an analytically-posed node
      int check = (flag_u[C] < 1) && (PP > -1);
      u[C] = check * Ux + (1 - check) * u[C];
#endif
    }
  }
}

__global__ void part_BC_v(real *v, int *phase, int *flag_v,
  part_struct *parts, dom_struct *dom,
  real nu, int stride,
  real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x + dom->Gfy._ksb;
  int ti = blockDim.y*blockIdx.y + threadIdx.y + dom->Gfy._isb;
  int C, CS, CN;
  real x, y, z;         // velocity node location Cartesian
  real X, Y, Z;         // particle position
  real r, theta, phi;   // velocity node location spherical
  real Ux, Uy, Uz;      // Cartesian-directed pressure gradients
  int P, PP, PS, PN;    // particle number
  real a;               // particle radius
  int order;            // particle order
  real oz, ox;          // particle angular velocity
  real ozdot, oxdot;    // particle angular acceleration
  real vv;              // particle velocity

  if(tk < dom->Gfy._keb && ti < dom->Gfy._ieb) {
    for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
      C = ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b;
      CS = ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      CN = ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      // get particle number
      PS = phase[CS];
      PN = phase[CN];
      if(PS > -1) {
        P = PS;
        PP = PS;
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        oz = parts[P].oz;
        ox = parts[P].ox;
        ozdot = parts[P].ozdot;
        oxdot = parts[P].oxdot;
        vv = parts[P].v;
      } else if(PN > -1) {
        P = PN;
        PP = PN;
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        oz = parts[P].oz;
        ox = parts[P].ox;
        ozdot = parts[P].ozdot;
        oxdot = parts[P].oxdot;
        vv = parts[P].v;
      } else {
        P = 0;
        PP = -1;
        a = (dom->dx + dom->dy + dom->dz) / 3.;
        X = (ti-0.5) * dom->dx + dom->xs + a;
        Y = (j-DOM_BUF) * dom->dy + dom->ys + a;
        Z = (tk-0.5) * dom->dz + dom->zs + a;
        order = 0;
        oz = 0;
        ox = 0;
        ozdot = 0;
        oxdot = 0;
        vv = 0;
      }
      x = (ti-0.5) * dom->dx + dom->xs - X;
      if(x < dom->xs - 0.5*dom->dx) x += dom->xl;
      if(x > dom->xe + 0.5*dom->dx) x -= dom->xl;
      y = (j-DOM_BUF) * dom->dy + dom->ys - Y;
      if(y < dom->ys - 0.5*dom->dy) y += dom->yl;
      if(y > dom->ye + 0.5*dom->dy) y -= dom->yl;
      z = (tk-0.5) * dom->dz + dom->zs - Z;
      if(z < dom->zs - 0.5*dom->dx) z += dom->zl;
      if(z > dom->ze + 0.5*dom->dz) z -= dom->zl;
      xyz2rtp(x, y, z, &r, &theta, &phi);

      // calculate analytic solution
#ifdef STEPS
      int check = (flag_v[C] < 1) && (PP > -1);
      v[C] = - (check - 1) * v[C];
#else
      lamb_vel(order, a, r, theta, phi,
        nu, pnm_re, pnm_im, phinm_re, phinm_im,
        chinm_re, chinm_im,
        P, stride, &Ux, &Uy, &Uz);

      // switch reference frame and set boundary condition
      real ocrossr_y = -(ox*z - oz*x);
      real odotcrossr_y = -(oxdot*z - ozdot*x);
      Uy += vv + ocrossr_y;
      Uy += 0.1/nu *(r*r*r*r*r-a*a*a*a*a)/(r*r*r) * odotcrossr_y;
      // boolean check if this is an analytically-posed node
      int check = (flag_v[C] < 1) && (PP > -1);
      v[C] = check * Uy + (1 - check) * v[C];
#endif
    }
  }
}

__global__ void part_BC_w(real *w, int *phase, int *flag_w,
  part_struct *parts, dom_struct *dom,
  real nu, int stride,
  real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x + dom->Gfz._isb;
  int tj = blockDim.y*blockIdx.y + threadIdx.y + dom->Gfz._jsb;
  int C, CB, CT;
  real x, y, z;         // velocity node location Cartesian
  real X, Y, Z;         // particle position
  real r, theta, phi;   // velocity node location spherical
  real Ux, Uy, Uz;      // Cartesian-directed pressure gradients
  int P, PP, PB, PT;    // particle number
  real a;               // particle radius
  int order;            // particle order
  real ox, oy;          // particle angular velocity
  real oxdot, oydot;    // particle angular acceleration
  real ww;              // particle velocity

  if(ti < dom->Gfz._ieb && tj < dom->Gfz._jeb) {
    for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
      C = ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b;
      CB = ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b;
      CT = ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b;
      // get particle number
      PB = phase[CB];
      PT = phase[CT];
      if(PB > -1) {
        P = PB;
        PP = PB;
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        ox = parts[P].ox;
        oy = parts[P].oy;
        oxdot = parts[P].oxdot;
        oydot = parts[P].oydot;
        ww = parts[P].w;
      } else if(PT > -1) {
        P = PT;
        PP = PT;
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        ox = parts[P].ox;
        oy = parts[P].oy;
        oxdot = parts[P].oxdot;
        oydot = parts[P].oydot;
        ww = parts[P].w;
      } else {
        P = 0;
        PP = -1;
        a = (dom->dx + dom->dy + dom->dz) / 3.;
        X = (ti-0.5) * dom->dx + dom->xs + a;
        Y = (tj-0.5) * dom->dy + dom->ys + a;
        Z = (k-DOM_BUF) * dom->dz + dom->zs + a;
        order = 0;
        ox = 0;
        oy = 0;
        oxdot = 0;
        oydot = 0;
        ww = 0;
      }
      x = (ti-0.5) * dom->dx + dom->xs - X;
      if(x < dom->xs - 0.5*dom->dx) x += dom->xl;
      if(x > dom->xe + 0.5*dom->dx) x -= dom->xl;
      y = (tj-0.5) * dom->dy + dom->ys - Y;
      if(y < dom->ys - 0.5*dom->dy) y += dom->yl;
      if(y > dom->ye + 0.5*dom->dy) y -= dom->yl;
      z = (k-DOM_BUF) * dom->dz + dom->zs - Z;
      if(z < dom->zs - 0.5*dom->dz) z += dom->zl;
      if(z > dom->ze + 0.5*dom->dz) z -= dom->zl;
      xyz2rtp(x, y, z, &r, &theta, &phi);

      // calculate analytic solution
#ifdef STEPS
      int check = (flag_w[C] < 1) && (PP > -1);
      w[C] = - (check - 1) * w[C];
#else
      lamb_vel(order, a, r, theta, phi,
        nu, pnm_re, pnm_im, phinm_re, phinm_im,
        chinm_re, chinm_im,
        P, stride, &Ux, &Uy, &Uz);

      // switch reference frame and set boundary condition
      real ocrossr_z = ox*y - oy*x;
      real odotcrossr_z = oxdot*y - oydot*x;
      Uz += ww + ocrossr_z;
      Uz += 0.1/nu *(r*r*r*r*r-a*a*a*a*a)/(r*r*r) * odotcrossr_z;
      // boolean check if this is an analytically-posed node
      int check = (flag_w[C] < 1) && (PP > -1);
      w[C] = check * Uz + (1 - check) * w[C];
#endif
    }
  }
}

__global__ void part_BC_p(real *p, real *p_rhs, int *phase, int *phase_shell,
  part_struct *parts, dom_struct *dom,
  real mu, real nu, real dt, real dt0, gradP_struct gradP, real rho_f, int stride,
  real *pnm_re00, real *pnm_im00, real *phinm_re00, real *phinm_im00,
  real *chinm_re00, real *chinm_im00,
  real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x + dom->Gcc._js;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + dom->Gcc._ks;
  int C, CC;
  real x, y, z;         // pressure node location Cartesian
  real X, Y, Z;         // particle position
  real r, theta, phi;   // velocity node location spherical
  real pp_tmp;//, pp_tmp00;// temporary pressure
  int P;            // particle number
  real a;               // particle radius
  int order;            // particle order
  real ox, oy, oz;      // particle angular velocity
  real udot, vdot, wdot;// particle acceleration

  if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
    for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
      CC = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc._s1 + (tk-DOM_BUF)*dom->Gcc._s2;
      // get particle number
      P = phase[CC];
      if(P > -1) {
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        order = parts[P].order;
        ox = parts[P].ox;
        oy = parts[P].oy;
        oz = parts[P].oz;
        udot = parts[P].udot;
        vdot = parts[P].vdot;
        wdot = parts[P].wdot;
      } else {
        P = 0;
        a = (dom->dx + dom->dy + dom->dz) / 3.;
        X = (i-0.5) * dom->dx + dom->xs + a;
        Y = (tj-0.5) * dom->dy + dom->ys + a;
        Z = (tk-0.5) * dom->dz + dom->zs + a;
        order = 0;
        ox = 0;
        oy = 0;
        oz = 0;
        udot = 0;
        vdot = 0;
        wdot = 0;
      }
      x = (i-0.5) * dom->dx + dom->xs - X;
      if(x < dom->xs) x += dom->xl;
      if(x > dom->xe) x -= dom->xl;
      y = (tj-0.5) * dom->dy + dom->ys - Y;
      if(y < dom->ys) y += dom->yl;
      if(y > dom->ye) y -= dom->yl;
      z = (tk-0.5) * dom->dz + dom->zs - Z;
      if(z < dom->zs) z += dom->zl;
      if(z > dom->ze) z -= dom->zl;
      xyz2rtp(x, y, z, &r, &theta, &phi);

      // calculate analytic solution
#ifdef STEPS
      p_rhs[C] = phase_shell[CC] * p_rhs[C];
#else
      real ar = a / r;
      real ra = r / a;
      pp_tmp = X_pn(0, theta, phi, pnm_re, pnm_im, P, stride);
      //pp_tmp00 = X_pn(0, theta, phi, pnm_re00, pnm_im00, P, stride);
      for(int n = 1; n <= order; n++) {
        pp_tmp += (1.-0.5*n*(2.*n-1.)/(n+1.)*pow(ar,2.*n+1.))*pow(ra,n)
          * X_pn(n, theta, phi, pnm_re, pnm_im, P, stride);
        //pp_tmp00 += (1.-0.5*n*(2.*n-1.)/(n+1.)*pow(ar,2.*n+1.))*pow(ra,n)
          //* X_pn(n, theta, phi, pnm_re00, pnm_im00, P, stride);
        pp_tmp -= n*(2.*n-1.)*(2.*n+1.)/(n+1.)*pow(ar,n+1.)
          * X_phin(n, theta, phi, phinm_re, phinm_im, P, stride);
        //pp_tmp00 -= n*(2.*n-1.)*(2.*n+1.)/(n+1.)*pow(ar,n+1.)
          //* X_phin(n, theta, phi, phinm_re00, phinm_im00, P, stride);
      }
      pp_tmp *= mu*nu/(a*a);
      //pp_tmp00 *= mu*nu/(a*a);
      real ocrossr2 = (oy*z - oz*y) * (oy*z - oz*y);
      ocrossr2 += (ox*z - oz*x) * (ox*z - oz*x);
      ocrossr2 += (ox*y - oy*x) * (ox*y - oy*x);
      real rhoV = rho_f;
      real accdotr = (-gradP.x/rhoV - udot)*x + (-gradP.y/rhoV - vdot)*y
        + (-gradP.z/rhoV - wdot)*z;
      pp_tmp += 0.5 * rho_f * ocrossr2 + rho_f * accdotr;
      //pp_tmp00 += 0.5 * rho_f * ocrossr2 + rho_f * accdotr;
      // write BC if flagged, otherwise leave alone
      p_rhs[C] = (real) phase_shell[CC] * p_rhs[C]
        + (real) (1 - phase_shell[CC])
        * (pp_tmp - p[CC]);// + 0.5*mu*p_rhs[C]);
      p_rhs[C] = (real) (phase[CC] < 0 && phase_shell[CC]) * p_rhs[C];
#endif
    }
  }
}

__global__ void part_BC_p_fill(real *p, int *phase,
  part_struct *parts, dom_struct *dom,
  real mu, real nu, real rho_f, gradP_struct gradP, int stride,
  real *pnm_re, real *pnm_im)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x + dom->Gcc._js;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + dom->Gcc._ks;
  int CC;
  real x, y, z;         // pressure node location Cartesian
  real X, Y, Z;         // particle position
  real r, theta, phi;   // velocity node location spherical
  real pp_tmp;//, pp_tmp00;// temporary pressure
  int P;            // particle number
  real a;               // particle radius
  real ox, oy, oz;      // particle angular velocity
  real udot, vdot, wdot;// particle acceleration

  if(tj < dom->Gcc._je && tk < dom->Gcc._ke) {
    for(int i = dom->Gcc._is; i < dom->Gcc._ie; i++) {
      CC = i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b;
      // get particle number
      P = phase[CC];
      if(P > -1) {
        a = parts[P].r;
        X = parts[P].x;
        Y = parts[P].y;
        Z = parts[P].z;
        ox = parts[P].ox;
        oy = parts[P].oy;
        oz = parts[P].oz;
        udot = parts[P].udot;
        vdot = parts[P].vdot;
        wdot = parts[P].wdot;
      } else {
        P = 0;
        a = (dom->dx + dom->dy + dom->dz) / 3.;
        X = (i-0.5) * dom->dx + dom->xs + a;
        Y = (tj-0.5) * dom->dy + dom->ys + a;
        Z = (tk-0.5) * dom->dz + dom->zs + a;
        ox = 0;
        oy = 0;
        oz = 0;
        udot = 0;
        vdot = 0;
        wdot = 0;
      }
      x = (i-0.5) * dom->dx + dom->xs - X;
      if(x < dom->xs) x += dom->xl;
      if(x > dom->xe) x -= dom->xl;
      y = (tj-0.5) * dom->dy + dom->ys - Y;
      if(y < dom->ys) y += dom->yl;
      if(y > dom->ye) y -= dom->yl;
      z = (tk-0.5) * dom->dz + dom->zs - Z;
      if(z < dom->zs) z += dom->zl;
      if(z > dom->ze) z -= dom->zl;
      xyz2rtp(x, y, z, &r, &theta, &phi);

      // calculate analytic solution
#ifdef STEPS
      p[CC] = phase_shell[CC] * p[CC];
#else
      pp_tmp = X_pn(0, theta, phi, pnm_re, pnm_im, P, stride);
      pp_tmp *= mu*nu/(a*a);
      real ocrossr2 = (oy*z - oz*y) * (oy*z - oz*y);
      ocrossr2 += (ox*z - oz*x) * (ox*z - oz*x);
      ocrossr2 += (ox*y - oy*x) * (ox*y - oy*x);
      real rhoV = rho_f;
      real accdotr = (-gradP.x/rhoV - udot)*x + (-gradP.y/rhoV - vdot)*y
        + (-gradP.z/rhoV - wdot)*z;
      pp_tmp += 0.5 * rho_f * ocrossr2 + rho_f * accdotr;
      // write BC if inside particle, otherwise leave alone
      p[CC] = (real) (phase[CC] > -1) * pp_tmp
        + (1 - (phase[CC] > -1)) * p[CC];
#endif
    }
  }
}

__device__ real Nnm(int n, int m)
{
  real fact_top = 1;
  real fact_bot = 1;

  for(int i = 1; i <= (n-m); i++) fact_top *= (real)i;
  for(int i = 1; i <= (n+m); i++) fact_bot *= (real)i;

  return sqrt((2.*n+1.) / 4. / PI * fact_top / fact_bot);
}

__device__ real Pnm(int n, int m, real theta)
{
  real x = cos(theta);
  real y = sin(theta);

  switch(n) {
    case 0: return 1;
    case 1:
      switch(m) {
        //case -1: return 0.5*y;
        case 0: return x;
        case 1: return -y;
      }
    case 2:
      switch(m) {
        //case -2: return 0.125*y*y;
        //case -1: return 0.5*x*y;
        case 0: return 0.5*(3.*x*x - 1.);
        case 1: return -3.*x*y;
        case 2: return 3.*y*y;
      }
    case 3:
      switch(m) {
        //case -3: return 0.02083333333333*y*y*y;
        //case -2: return 0.125*x*y*y;
        //case -1: return -0.125*(1. - 5.*x*x)*y;
        case 0: return 0.5*x*(5.*x*x - 3.);
        case 1: return -1.5*(5.*x*x - 1.)*y;
        case 2: return 15.*x*y*y;
        case 3: return -15.*y*y*y;
      }
    case 4:
      switch(m) {
        //case -4: return .002604166666667*y*y*y*y;
        //case -3: return 0.02083333333333*x*y*y*y*y;
        //case -2: return 0.02083333333333*(7.*x*x - 1.)*y*y;
        //case -1: return -0.125*x*(3. - 7.*x*x)*y;
        case 0: return 0.125*(35.*x*x*x*x - 30.*x*x + 3.);
        case 1: return -2.5*(7.*x*x - 3.)*x*y;
        case 2: return 7.5*(7.*x*x - 1.)*y*y;
        case 3: return -105.*x*y*y*y;
        case 4: return 105.*y*y*y*y;
      }
    case 5:
      switch(m) {
        //case -5: return 0.000260416666667*y*y*y*y*y;
        //case -4: return 0.002604166666667*x*y*y*y*y;
        //case -3: return -0.002604166666667*y*y*y*(9.*x*x - 1.);
        //case -2: return 0.0625*x*y*y*(3.*x*x - 1.);
        //case -1: return -0.0625*(21.*x*x*x*x - 14.*x*x + 1.);
        case 0: return 0.125*x*(63.*x*x*x*x - 70.*x*x + 15.);
        case 1: return -1.875*y*(21.*x*x*x*x - 14.*x*x + 1.);
        case 2: return 52.5*x*y*y*(3.*x*x - 1.);
        case 3: return -52.5*y*y*y*(9.*x*x - 1.);
        case 4: return 945.*x*y*y*y*y;
        case 5: return -945.*y*y*y*y*y;
      }
  }
  return 0; // this should never be reached
}

__device__ void xyz2rtp(real x, real y, real z, real *r, real *theta, real *phi)
{
  real XY = x*x + y*y;
  real XYZ = XY + z*z;
  // We calculate the coefficients everywhere in space. If a particle is
  // centered at the center of a cell, XYZ will be zero.
  if(XYZ >= 0 && XYZ < DIV_ST) XYZ = DIV_ST;
  *r = sqrt(XYZ);
  *theta = acos(z / *r);
  // Note that XY cannot be set equal to one, because the values are used.
  if(XY >= 0 && XY < DIV_ST) XY = DIV_ST;
  *phi = acos(x / sqrt(XY));
  if(y < 0.) *phi = 2.*PI - *phi;
}

__device__ real X_pn(int n, real theta, real phi,
  real *pnm_re, real *pnm_im, int pp, int stride)
{
  int coeff = 0;
  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = Nnm(n,m)*Pnm(n,m,theta)*pnm_re[coeff];

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += 2.*Nnm(n,m)*Pnm(n,m,theta)
      *(pnm_re[coeff]*cos(m*phi)
      - pnm_im[coeff]*sin(m*phi));
  }

  return sum;
}

__device__ real X_phin(int n, real theta, real phi,
  real *phinm_re, real *phinm_im, int pp, int stride)
{
  int coeff = 0;
  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = Nnm(n,m)*Pnm(n,m,theta)*phinm_re[coeff];

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += 2.*Nnm(n,m)*Pnm(n,m,theta)
      *(phinm_re[coeff]*cos(m*phi)
      - phinm_im[coeff]*sin(m*phi));
  }

  return sum;
}

__device__ real Y_pn(int n, real theta, real phi,
  real *pnm_re, real *pnm_im, int pp, int stride)
{
  int coeff = 0;
  real ct = cos(theta);
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;

  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = Nnm(n,m)
    *(-(n+1)*ct/st*Pnm(n,m,theta)+(n-m+1)/st*Pnm(n+1,m,theta))
    *pnm_re[coeff];

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += 2.*Nnm(n,m)
      *(-(n+1)*ct/st*Pnm(n,m,theta)+(n-m+1)/st*Pnm(n+1,m,theta))
      *(pnm_re[coeff]*cos(m*phi)
      - pnm_im[coeff]*sin(m*phi));
  }

  return sum;
}

__device__ real Y_phin(int n, real theta, real phi,
  real *phinm_re, real *phinm_im, int pp, int stride)
{
  int coeff = 0;
  real ct = cos(theta);
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;

  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = Nnm(n,m)
    *(-(n+1)*ct/st*Pnm(n,m,theta)+(n-m+1)/st*Pnm(n+1,m,theta))
    *phinm_re[coeff];

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += 2.*Nnm(n,m)
      *(-(n+1)*ct/st*Pnm(n,m,theta)+(n-m+1)/st*Pnm(n+1,m,theta))
      *(phinm_re[coeff]*cos(m*phi)
      - phinm_im[coeff]*sin(m*phi));
  }

  return sum;
}

__device__ real Y_chin(int n, real theta, real phi,
  real *chinm_re, real *chinm_im, int pp, int stride)
{
  int coeff = 0;
  real ct = cos(theta);
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;

  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = Nnm(n,m)
    *(-(n+1)*ct/st*Pnm(n,m,theta)+(n-m+1)/st*Pnm(n+1,m,theta))
    *chinm_re[coeff];

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += 2.*Nnm(n,m)
      *(-(n+1)*ct/st*Pnm(n,m,theta)+(n-m+1)/st*Pnm(n+1,m,theta))
      *(chinm_re[coeff]*cos(m*phi)
      - chinm_im[coeff]*sin(m*phi));
  }

  return sum;
}

__device__ real Z_pn(int n, real theta, real phi,
  real *pnm_re, real *pnm_im, int pp, int stride)
{
  int coeff = 0;
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;

  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = 0.;

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += -2.*m/st*Nnm(n,m)*Pnm(n,m,theta)
      *(pnm_re[coeff]*sin(m*phi)
      + pnm_im[coeff]*cos(m*phi));
  }

  return sum;
}

__device__ real Z_phin(int n, real theta, real phi,
  real *phinm_re, real *phinm_im, int pp, int stride)
{
  int coeff = 0;
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;

  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = 0.;

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += -2.*m/st*Nnm(n,m)*Pnm(n,m,theta)
      *(phinm_re[coeff]*sin(m*phi)
      + phinm_im[coeff]*cos(m*phi));
  }

  return sum;
}

__device__ real Z_chin(int n, real theta, real phi,
  real *chinm_re, real *chinm_im, int pp, int stride)
{
  int coeff = 0;
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;

  for(int j = 0; j < n; j++) coeff += j+1;

  coeff = coeff + pp*stride;

  int m = 0;
  real sum = 0.;

  for(m = 1; m <= n; m++) {
    coeff++;
    sum += -2.*m/st*Nnm(n,m)*Pnm(n,m,theta)
      *(chinm_re[coeff]*sin(m*phi)
      + chinm_im[coeff]*cos(m*phi));
  }

  return sum;
}

__device__ void lamb_vel(int order, real a, real r, real theta, real phi,
  real nu, real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im,
  int p_ind, int stride, real *Ux, real *Uy, real *Uz)
{
  real ar = a / r;
  real ra = r / a;

  real ur = 0.;
  real ut = 0.5*ra*Y_pn(0, theta, phi, pnm_re, pnm_im, p_ind, stride);
  real up = 0.5*ra*Z_pn(0, theta, phi, pnm_re, pnm_im, p_ind, stride);


  for(int n = 1; n <= order; n++) {
    real powranp1 = pow(ra,n+1.);
    real powranm1 = pow(ra,n-1.);

    real powarnp2 = pow(ar,n+2.);
    real powarnp1 = pow(ar,n+1.);
    real powarn = pow(ar,n);

    real od2np3 = 1./(2.*n+3.);
    real odnp1 = 1./(n+1.);

    ur += (0.5*n*od2np3*powranp1
      + 0.25*n*((2.*n+1.)*od2np3*ar*ar-1.)*powarn)
      * X_pn(n, theta, phi, pnm_re, pnm_im, p_ind, stride);
    ur += (n*powranm1
      + 0.5*n*(2.*n-1.-(2.*n+1.)*ra*ra)*pow(ar,n+2.))
      * X_phin(n, theta, phi, phinm_re, phinm_im, p_ind, stride);

    ut += (0.5*(n+3.)*odnp1*od2np3*powranp1
      + 0.25*odnp1*(n-2.-n*(2.*n+1.)*od2np3*ar*ar)*powarn)
      * Y_pn(n, theta, phi, pnm_re, pnm_im, p_ind, stride);
    ut += (powranm1
      + 0.5*odnp1*((n-2.)*(2.*n+1.)*ra*ra-n*(2.*n-1.))*powarnp2)
      * Y_phin(n, theta, phi, phinm_re, phinm_im, p_ind, stride);
    ut += (powranm1
      - powarnp1)
      * Z_chin(n, theta, phi, chinm_re, chinm_im, p_ind, stride);
      
    up += (0.5*(n+3.)*odnp1*od2np3*powranp1
      + 0.25*odnp1*(n-2.-n*(2.*n+1.)*od2np3*ar*ar)*powarn)
      * Z_pn(n, theta, phi, pnm_re, pnm_im, p_ind, stride);
    up += (powranm1
      + 0.5*odnp1*((n-2.)*(2.*n+1.)*ra*ra-n*(2.*n-1.))*powarnp2)
      * Z_phin(n, theta, phi, phinm_re, phinm_im, p_ind, stride);
    up += (-powranm1
      + powarnp1)
      * Y_chin(n, theta, phi, chinm_re, chinm_im, p_ind, stride);
  }
  ur *= nu / a;
  ut *= nu / a;
  up *= nu / a;

  *Ux = ur*sin(theta)*cos(phi)+ut*cos(theta)*cos(phi)-up*sin(phi);
  *Uy = ur*sin(theta)*sin(phi)+ut*cos(theta)*sin(phi)+up*cos(phi);
  *Uz = ur*cos(theta)-ut*sin(theta);
}

__device__ void lamb_gradP(int order, real a, real r, real theta, real phi,
  real mu, real nu, real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  int p_ind, int stride, real *gPx, real *gPy, real *gPz)
{
  real ar = a / r;
  real ra = r / a;
  real st = sin(theta);
  if(st >= 0 && st < DIV_ST) st = DIV_ST;
  else if(st < 0 && st > -DIV_ST) st = -DIV_ST;
  
  real pr = 0;
  real pt = 1./r * Y_pn(0, theta, phi, pnm_re, pnm_im, p_ind, stride);
  real pp = 1./r * Z_pn(0, theta, phi, pnm_re, pnm_im, p_ind, stride);
  for(int n = 1; n <= order; n++) {
    pr += n * (ar + 0.5*(2.*n-1.)*pow(ar,2.*n+2.))*pow(ra,n)
      * X_pn(n, theta, phi, pnm_re, pnm_im, p_ind, stride);
    pr += n * (2.*n-1.)*(2.*n+1.)*pow(ar,n+2.)
      * X_phin(n, theta, phi, phinm_re, phinm_im, p_ind, stride);

    pt += (1. - 0.5*n*(2.*n-1.)/(n+1.)*pow(ar,2.*n+1.))*pow(ra,n)
      * Y_pn(n, theta, phi, pnm_re, pnm_im, p_ind, stride);
    pt += -n*(2.*n-1.)*(2.*n+1.)/(n+1.)*pow(ar,n+1.)
      * Y_phin(n, theta, phi, phinm_re, phinm_im, p_ind, stride);

    pp += (1. - 0.5*n*(2.*n-1.)/(n+1.)*pow(ar,2.*n+1.))*pow(ra,n)
      * Z_pn(n, theta, phi, pnm_re, pnm_im, p_ind, stride);
    pp += -n*(2.*n-1.)*(2.*n+1.)/(n+1.)*pow(ar,n+1.)
      * Z_phin(n, theta, phi, phinm_re, phinm_im, p_ind, stride);
  }
  pr *= mu * nu / (a*a*a);
  pt *= mu * nu / (a*a);
  pp *= mu * nu / (a*a);

  *gPx = pr*sin(theta)*cos(phi)+pt*cos(theta)*cos(phi)/r-pp*sin(phi)/r/st;
  *gPy = pr*sin(theta)*sin(phi)+pt*cos(theta)*sin(phi)/r+pp*cos(phi)/r/st;
  *gPz = pr*cos(theta)-pt*sin(theta)/r;
}

__global__ void predict_coeffs(real dt0, real dt,
  real *pnm_re00, real *pnm_im00, real *phinm_re00, real *phinm_im00,
  real *chinm_re00, real *chinm_im00,
  real *pnm_re0, real *pnm_im0, real *phinm_re0, real *phinm_im0,
  real *chinm_re0, real *chinm_im0,
  real *pnm_re, real *pnm_im, real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im, int stride)
{
  int coeff = threadIdx.x;
  int part = blockIdx.x;

  // extrapolate coefficients
  int ti = stride*part + coeff;
  real a = dt/dt0;
  // store last result
  /*pnm_re0[ti] = pnm_re[ti];
  pnm_im0[ti] = pnm_im[ti];
  phinm_re0[ti] = phinm_re[ti];
  phinm_im0[ti] = phinm_im[ti];
  chinm_re0[ti] = chinm_re[ti];
  chinm_im0[ti] = chinm_im[ti];
*/
  // predict starting point for iterations at the next timestep
  pnm_re[ti] = pnm_re[ti]*(1. + a) - pnm_re00[ti]*a;
  pnm_im[ti] = pnm_im[ti]*(1. + a) - pnm_im00[ti]*a;
  phinm_re[ti] = phinm_re[ti]*(1. + a) - phinm_re00[ti]*a;
  phinm_im[ti] = phinm_im[ti]*(1. + a) - phinm_im00[ti]*a;
  chinm_re[ti] = chinm_re[ti]*(1. + a) - chinm_re00[ti]*a;
  chinm_im[ti] = chinm_im[ti]*(1. + a) - chinm_im00[ti]*a;
}
