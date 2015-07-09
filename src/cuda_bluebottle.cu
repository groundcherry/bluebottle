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

#include <cuda.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/sort.h>

#include "cuda_bicgstab.h"
#include "cuda_bluebottle.h"
#include "cuda_particle.h"
#include "entrySearch.h"

extern "C"
void cuda_dom_malloc(void)
{
  // allocate device memory on host
  _dom = (dom_struct**) malloc(nsubdom * sizeof(dom_struct*));
  cpumem += nsubdom * sizeof(dom_struct*);
  _p0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _p = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _phi = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  //_divU = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _u = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _v = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _w = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _u0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _v0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _w0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _f_x = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _f_y = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _f_z = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
#ifndef IMPLICIT
  _diff0_u = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff0_v = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff0_w = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
#endif
  _conv0_u = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv0_v = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv0_w = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff_u = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff_v = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff_w = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv_u = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv_v = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv_w = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _u_star = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _v_star = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _w_star = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _u_WE = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _u_SN = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _u_BT = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _v_WE = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _v_SN = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _v_BT = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _w_WE = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _w_SN = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _w_BT = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _rhs_p = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);

  // allocate device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    (cudaMalloc((void**) &(_dom[dev]),
      sizeof(dom_struct)));
    gpumem += sizeof(dom_struct);
    // copy domain info to _dom
    (cudaMemcpy(_dom[dev], &dom[dev], sizeof(dom_struct),
      cudaMemcpyHostToDevice));

    (cudaMalloc((void**) &(_p0[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);

    (cudaMalloc((void**) &(_p[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);

    (cudaMalloc((void**) &(_phi[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);

    //(cudaMalloc((void**) &(_divU[dev]),
      //sizeof(real) * dom[dev].Gcc.s3b));
    //gpumem += dom[dev].Gcc.s3b * sizeof(real);

    (cudaMalloc((void**) &(_u[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_v[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_w[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
    (cudaMalloc((void**) &(_u0[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_v0[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_w0[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
    (cudaMalloc((void**) &(_f_x[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_f_y[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_f_z[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
    
#ifndef IMPLICIT
    (cudaMalloc((void**) &(_diff0_u[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_diff0_v[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_diff0_w[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
#endif
    (cudaMalloc((void**) &(_conv0_u[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_conv0_v[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_conv0_w[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
    (cudaMalloc((void**) &(_diff_u[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_diff_v[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_diff_w[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
    (cudaMalloc((void**) &(_conv_u[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_conv_v[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_conv_w[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);
    (cudaMalloc((void**) &(_u_star[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    gpumem += dom[dev].Gfx.s3b * sizeof(real);
    (cudaMalloc((void**) &(_v_star[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    gpumem += dom[dev].Gfy.s3b * sizeof(real);
    (cudaMalloc((void**) &(_w_star[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
    gpumem += dom[dev].Gfz.s3b * sizeof(real);

    (cudaMalloc((void**) &(_u_WE[dev]),
      sizeof(real) * dom[dev].Gfx.jnb*dom[dev].Gfx.knb));
    gpumem += dom[dev].Gfx.jnb*dom[dev].Gfx.knb * sizeof(real);
    (cudaMalloc((void**) &(_u_SN[dev]),
      sizeof(real) * dom[dev].Gfx.inb*dom[dev].Gfx.knb));
    gpumem += dom[dev].Gfx.inb*dom[dev].Gfx.knb * sizeof(real);
    (cudaMalloc((void**) &(_u_BT[dev]),
      sizeof(real) * dom[dev].Gfx.inb*dom[dev].Gfx.jnb));
    gpumem += dom[dev].Gfx.inb*dom[dev].Gfx.jnb * sizeof(real);
    (cudaMalloc((void**) &(_v_WE[dev]),
      sizeof(real) * dom[dev].Gfy.jnb*dom[dev].Gfy.knb));
    gpumem += dom[dev].Gfy.jnb*dom[dev].Gfy.knb * sizeof(real);
    (cudaMalloc((void**) &(_v_SN[dev]),
      sizeof(real) * dom[dev].Gfy.inb*dom[dev].Gfy.knb));
    gpumem += dom[dev].Gfy.inb*dom[dev].Gfy.knb * sizeof(real);
    (cudaMalloc((void**) &(_v_BT[dev]),
      sizeof(real) * dom[dev].Gfy.inb*dom[dev].Gfy.jnb));
    gpumem += dom[dev].Gfy.inb*dom[dev].Gfy.jnb * sizeof(real);
    (cudaMalloc((void**) &(_w_WE[dev]),
      sizeof(real) * dom[dev].Gfz.jnb*dom[dev].Gfz.knb));
    gpumem += dom[dev].Gfz.jnb*dom[dev].Gfz.knb * sizeof(real);
    (cudaMalloc((void**) &(_w_SN[dev]),
      sizeof(real) * dom[dev].Gfz.inb*dom[dev].Gfz.knb));
    gpumem += dom[dev].Gfz.inb*dom[dev].Gfz.knb * sizeof(real);
    (cudaMalloc((void**) &(_w_BT[dev]),
      sizeof(real) * dom[dev].Gfz.inb*dom[dev].Gfz.jnb));
    gpumem += dom[dev].Gfz.inb*dom[dev].Gfz.jnb * sizeof(real);

    (cudaMalloc((void**) &(_rhs_p[dev]),
      sizeof(real) * (dom[dev].Gcc.s3)));
    gpumem += dom[dev].Gcc.s3 * sizeof(real);

    // TODO add CUSP solver data structures to memory usage count

    //printf("Device %d of %d using %f Mb global memory.\n", dev, nsubdom, mb);
  }
}

extern "C"
void cuda_dom_push(void)
{
  // copy host data to device
  #pragma omp parallel num_threads(nsubdom)
  {
    int i, j, k;          // iterators
    int ii, jj, kk;       // helper iterators
    int C, CC;            // cell references

    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    // set up host working arrays for subdomain copy from host to device
    real *pp0 = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // cpumem += dom[dev].Gcc.s3b * sizeof(real);
    real *pp = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *pphi = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // cpumem += dom[dev].Gcc.s3b * sizeof(real);
    //real *pdivU = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // cpumem += dom[dev].Gcc.s3b * sizeof(real);
    real *uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *uu0 = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *vv0 = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *ww0 = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *uu_star = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *vv_star = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *ww_star = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
#ifndef IMPLICIT
    real *diff0_uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *diff0_vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *diff0_ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
#endif
    real *conv0_uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *conv0_vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *conv0_ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *diff_uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *diff_vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *diff_ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *conv_uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *conv_vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *conv_ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);

    // select appropriate subdomain
    // p
    for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
      for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
        for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
          ii = i - dom[dev].Gcc.isb;
          jj = j - dom[dev].Gcc.jsb;
          kk = k - dom[dev].Gcc.ksb;
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
          pp0[CC] = p0[C];
          pp[CC] = p[C];
          pphi[CC] = phi[C];
          //pdivU[CC] = divU[C];
        }
      }
    }
    // u
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          uu[CC] = u[C];
          uu0[CC] = u0[C];
          uu_star[CC] = u_star[C];
#ifndef IMPLICIT
          diff0_uu[CC] = diff0_u[C];
#endif
          conv0_uu[CC] = conv0_u[C];
          diff_uu[CC] = diff_u[C];
          conv_uu[CC] = conv_u[C];
        }
      }
    }
    // v
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          vv[CC] = v[C];
          vv0[CC] = v0[C];
          vv_star[CC] = v_star[C];
#ifndef IMPLICIT
          diff0_vv[CC] = diff0_v[C];
#endif
          conv0_vv[CC] = conv0_v[C];
          diff_vv[CC] = diff_v[C];
          conv_vv[CC] = conv_v[C];
        }
      }
    }
    // w
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          ww[CC] = w[C];
          ww0[CC] = w0[C];
          ww_star[CC] = w_star[C];
#ifndef IMPLICIT
          diff0_ww[CC] = diff0_w[C];
#endif
          conv0_ww[CC] = conv0_w[C];
          diff_ww[CC] = diff_w[C];
          conv_ww[CC] = conv_w[C];
        }
      }
    }

    // copy from host to device
    (cudaMemcpy(_p0[dev], pp0, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_p[dev], pp, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_phi[dev], pphi, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    //(cudaMemcpy(_divU[dev], pdivU, sizeof(real) * dom[dev].Gcc.s3b,
      //cudaMemcpyHostToDevice));
    (cudaMemcpy(_u[dev], uu, sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_v[dev], vv, sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_w[dev], ww, sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_u0[dev], uu0, sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_v0[dev], vv0, sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_w0[dev], ww0, sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_u_star[dev], uu_star,
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_v_star[dev], vv_star,
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_w_star[dev], ww_star,
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
#ifndef IMPLICIT
    (cudaMemcpy(_diff0_u[dev], diff0_uu,
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_diff0_v[dev], diff0_vv,
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_diff0_w[dev], diff0_ww,
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
#endif
    (cudaMemcpy(_conv0_u[dev], conv0_uu,
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_conv0_v[dev], conv0_vv,
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_conv0_w[dev], conv0_ww,
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_diff_u[dev], diff_uu,
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_diff_v[dev], diff_vv,
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_diff_w[dev], diff_ww,
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_conv_u[dev], conv_uu,
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_conv_v[dev], conv_vv,
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyHostToDevice));
    (cudaMemcpy(_conv_w[dev], conv_ww,
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyHostToDevice));

    // free host subdomain working arrays
    free(pp0);
    free(pp);
    free(pphi);
    //free(pdivU);
    free(uu);
    free(vv);
    free(ww);
    free(uu0);
    free(vv0);
    free(ww0);
    free(uu_star);
    free(vv_star);
    free(ww_star);
#ifndef IMPLICIT
    free(diff0_uu);
    free(diff0_vv);
    free(diff0_ww);
#endif
    free(conv0_uu);
    free(conv0_vv);
    free(conv0_ww);
    free(diff_uu);
    free(diff_vv);
    free(diff_ww);
    free(conv_uu);
    free(conv_vv);
    free(conv_ww);
  }
}

extern "C"
void cuda_dom_turb_planes_push(int *bc_configs)
{
  // copy host data to device
  #pragma omp parallel num_threads(nsubdom)
  {
    int i, j, k;    // iterators
    int ii, jj, kk; // helper iterators
    int C, CC;      // cell references

    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    if(bc_configs[ 0] == PRECURSOR || bc_configs[ 1] == PRECURSOR) {
      // working array
      real *uu_WE = (real*) malloc(dom[dev].Gfx.jnb*dom[dev].Gfx.knb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
        for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = j + k * Dom.Gfx.jnb;
          CC = jj + kk * dom[dev].Gfx.jnb;
          uu_WE[CC] = u_WE[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_u_WE[dev], uu_WE, sizeof(real)
        * dom[dev].Gfx.jnb*dom[dev].Gfx.knb, cudaMemcpyHostToDevice));

      // clean up
      free(uu_WE);
    }
    if(bc_configs[ 2] == PRECURSOR || bc_configs[ 3] == PRECURSOR) {
      // working array
      real *uu_SN = (real*) malloc(dom[dev].Gfx.inb*dom[dev].Gfx.knb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + k * Dom.Gfx.inb;
          CC = ii + kk * dom[dev].Gfx.inb;
          uu_SN[CC] = u_SN[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_u_SN[dev], uu_SN, sizeof(real)
        * dom[dev].Gfx.inb*dom[dev].Gfx.knb, cudaMemcpyHostToDevice));

      // clean up
      free(uu_SN);
    }
    if(bc_configs[ 4] == PRECURSOR || bc_configs[ 5] == PRECURSOR) {
      // working array
      real *uu_BT = (real*) malloc(dom[dev].Gfx.inb*dom[dev].Gfx.jnb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          C = i + j * Dom.Gfx.inb;
          CC = ii + jj * dom[dev].Gfx.inb;
          uu_BT[CC] = u_BT[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_u_BT[dev], uu_BT, sizeof(real)
        * dom[dev].Gfx.inb*dom[dev].Gfx.jnb, cudaMemcpyHostToDevice));

      // clean up
      free(uu_BT);
    }

    if(bc_configs[ 6] == PRECURSOR || bc_configs[ 7] == PRECURSOR) {
      // working array
      real *vv_WE = (real*) malloc(dom[dev].Gfy.jnb*dom[dev].Gfy.knb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
        for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = j + k * Dom.Gfy.jnb;
          CC = jj + kk * dom[dev].Gfy.jnb;
          vv_WE[CC] = v_WE[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_v_WE[dev], vv_WE, sizeof(real)
        * dom[dev].Gfy.jnb*dom[dev].Gfy.knb, cudaMemcpyHostToDevice));

      // clean up
      free(vv_WE);
    }
    if(bc_configs[ 8] == PRECURSOR || bc_configs[ 9] == PRECURSOR) {
      // working array
      real *vv_SN = (real*) malloc(dom[dev].Gfy.inb*dom[dev].Gfy.knb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + k * Dom.Gfy.inb;
          CC = ii + kk * dom[dev].Gfy.inb;
          vv_SN[CC] = v_SN[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_v_SN[dev], vv_SN, sizeof(real)
        * dom[dev].Gfy.inb*dom[dev].Gfy.knb, cudaMemcpyHostToDevice));

      // clean up
      free(vv_SN);
    }
    if(bc_configs[10] == PRECURSOR || bc_configs[11] == PRECURSOR) {
      // working array
      real *vv_BT = (real*) malloc(dom[dev].Gfy.inb*dom[dev].Gfy.jnb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          C = i + j * Dom.Gfy.inb;
          CC = ii + jj * dom[dev].Gfy.inb;
          vv_BT[CC] = v_BT[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_v_BT[dev], vv_BT, sizeof(real)
        * dom[dev].Gfy.inb*dom[dev].Gfy.jnb, cudaMemcpyHostToDevice));

      // clean up
      free(vv_BT);
    }

    if(bc_configs[12] == PRECURSOR || bc_configs[13] == PRECURSOR) {
      // working array
      real *ww_WE = (real*) malloc(dom[dev].Gfz.jnb*dom[dev].Gfz.knb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
        for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = j + k * Dom.Gfz.jnb;
          CC = jj + kk * dom[dev].Gfz.jnb;
          ww_WE[CC] = w_WE[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_w_WE[dev], ww_WE, sizeof(real)
        * dom[dev].Gfz.jnb*dom[dev].Gfz.knb, cudaMemcpyHostToDevice));

      // clean up
      free(ww_WE);
    }
    if(bc_configs[14] == PRECURSOR || bc_configs[15] == PRECURSOR) {
      // working array
      real *ww_SN = (real*) malloc(dom[dev].Gfz.inb*dom[dev].Gfz.knb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + k * Dom.Gfz.inb;
          CC = ii + kk * dom[dev].Gfz.inb;
          ww_SN[CC] = w_SN[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_w_SN[dev], ww_SN, sizeof(real)
        * dom[dev].Gfz.inb*dom[dev].Gfz.knb, cudaMemcpyHostToDevice));

      // clean up
      free(ww_SN);
    }
    if(bc_configs[16] == PRECURSOR || bc_configs[17] == PRECURSOR) {
      // working array
      real *ww_BT = (real*) malloc(dom[dev].Gfz.inb*dom[dev].Gfz.jnb
        * sizeof(real));

      // select appropriate subdomain
      // set position
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          C = i + j * Dom.Gfz.inb;
          CC = ii + jj * dom[dev].Gfz.inb;
          ww_BT[CC] = w_BT[C];
        }
      }

      // copy from host to device
      (cudaMemcpy(_w_BT[dev], ww_BT, sizeof(real)
        * dom[dev].Gfz.inb*dom[dev].Gfz.jnb, cudaMemcpyHostToDevice));

      // clean up
      free(ww_BT);
    }

  }
}

extern "C"
void cuda_dom_pull(void)
{
  // copy device data to host
  #pragma omp parallel num_threads(nsubdom)
  {
    int i, j, k;          // iterators
    int ii, jj, kk;       // helper iterators
    int C, CC;            // cell references

    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    // host working arrays for subdomain copy from device to host
    real *pp0 = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // cpumem += dom[dev].Gcc.s3b * sizeof(real);
    real *pp = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *pphi = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // cpumem += dom[dev].Gcc.s3b * sizeof(real);
    //real *pdivU = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // cpumem += dom[dev].Gcc.s3b * sizeof(real);
    real *uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *uu0 = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *vv0 = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *ww0 = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *diffuu0 = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *diffvv0 = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *diffww0 = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *convuu0 = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *convvv0 = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *convww0 = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *diffuu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *diffvv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *diffww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    real *convuu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *convvv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *convww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));

    // copy from device to host
    (cudaMemcpy(pp0, _p0[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(pp, _p[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(pphi, _phi[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost)); 
    //(cudaMemcpy(pdivU , _divU [dev],
      //sizeof(real) * dom[dev].Gcc.s3b, cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(uu, _u[dev], sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(vv, _v[dev], sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(ww, _w[dev], sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(uu0, _u0[dev], sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(vv0, _v0[dev], sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(ww0, _w0[dev], sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyDeviceToHost)); 
#ifndef IMPLICIT
    (cudaMemcpy(diffuu0, _diff0_u[dev],
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(diffvv0, _diff0_v[dev],
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(diffww0, _diff0_w[dev],
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost));
#endif
    (cudaMemcpy(convuu0, _conv0_u[dev],
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(convvv0, _conv0_v[dev],
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(convww0, _conv0_w[dev],
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(diffuu, _diff_u[dev],
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(diffvv, _diff_v[dev],
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(diffww, _diff_w[dev],
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(convuu, _conv_u[dev],
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(convvv, _conv_v[dev],
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost));
    (cudaMemcpy(convww, _conv_w[dev],
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost));

    real *uu_star = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    // cpumem += dom[dev].Gfx.s3b * sizeof(real);
    real *vv_star = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    // cpumem += dom[dev].Gfy.s3b * sizeof(real);
    real *ww_star = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
    // cpumem += dom[dev].Gfz.s3b * sizeof(real);
    (cudaMemcpy(uu_star, _u_star[dev],
      sizeof(real) * dom[dev].Gfx.s3b, cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(vv_star, _v_star[dev],
      sizeof(real) * dom[dev].Gfy.s3b, cudaMemcpyDeviceToHost)); 
    (cudaMemcpy(ww_star, _w_star[dev],
      sizeof(real) * dom[dev].Gfz.s3b, cudaMemcpyDeviceToHost)); 

#ifdef DEBUG // run test code
    // fill in apropriate subdomain (copy back ghost cells)
    // p
    for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
      for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
        for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
          ii = i - dom[dev].Gcc.isb;
          jj = j - dom[dev].Gcc.jsb;
          kk = k - dom[dev].Gcc.ksb;
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
          p0[C] = pp0[CC];
          phi[C] = pphi[CC];
          p[C] = pp[CC];
          //divU[C] = pdivU[CC];
        }
      }
    }
    // u
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          u[C] = uu[CC];
          u0[C] = uu0[CC];
        }
      }
    }
    // v
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          v[C] = vv[CC];
          v0[C] = vv0[CC];
        }
      }
    }
    // w
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          w[C] = ww[CC];
          w0[C] = ww0[CC];
        }
      }
    }
    // u
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          u_star[C] = uu_star[CC];
        }
      }
    }
    // v
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          v_star[C] = vv_star[CC];
        }
      }
    }
    // w
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          w_star[C] = ww_star[CC];
        }
      }
    }

#else // run simulation
    // fill in apropriate subdomain
    // p
    for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
      for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
        for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
          ii = i - dom[dev].Gcc.isb;
          jj = j - dom[dev].Gcc.jsb;
          kk = k - dom[dev].Gcc.ksb;
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
          p0[C] = pp0[CC];
          p[C] = pp[CC];
          phi[C] = pphi[CC];
          //divU[C] = pdivU[CC];
        }
      }
    }
    // u
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          u[C] = uu[CC];
          u0[C] = uu0[CC];
        }
      }
    }
    // v
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          v[C] = vv[CC];
          v0[C] = vv0[CC];
        }
      }
    }
    // w
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          w[C] = ww[CC];
          w0[C] = ww0[CC];
        }
      }
    }
    // conv_u
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
#ifndef IMPLICIT
          diff0_u[C] = diffuu0[CC];
#endif
          conv0_u[C] = convuu0[CC];
          diff_u[C] = diffuu[CC];
          conv_u[C] = convuu[CC];
        }
      }
    }
    // conv_v
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
#ifndef IMPLICIT
          diff0_v[C] = diffvv0[CC];
#endif
          conv0_v[C] = convvv0[CC];
          diff_v[C] = diffvv[CC];
          conv_v[C] = convvv[CC];
        }
      }
    }
    // conv_w
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
#ifndef IMPLICIT
          diff0_w[C] = diffww0[CC];
#endif
          conv0_w[C] = convww0[CC];
          diff_w[C] = diffww[CC];
          conv_w[C] = convww[CC];
        }
      }
    }

    // u
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          u_star[C] = uu_star[CC];
        }
      }
    }
    // v
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          v_star[C] = vv_star[CC];
        }
      }
    }
    // w
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          w_star[C] = ww_star[CC];
        }
      }
    }

#endif

    // free host subdomain working arrays
    free(pp0);
    free(pp);
    free(pphi);
    //free(pdivU);
    free(uu);
    free(vv);
    free(ww);
    free(uu0);
    free(vv0);
    free(ww0);
    free(diffuu0);
    free(diffvv0);
    free(diffww0);
    free(convuu0);
    free(convvv0);
    free(convww0);
    free(diffuu);
    free(diffvv);
    free(diffww);
    free(convuu);
    free(convvv);
    free(convww);
    free(uu_star);
    free(vv_star);
    free(ww_star);
  }
}

extern "C"
void cuda_dom_turb_planes_pull(int *bc_configs)
{
  // copy device data to host
  #pragma omp parallel num_threads(nsubdom)
  {
    int i, j, k;    // iterators
    int ii, jj, kk; // helper iterators
    int C, CC;      // cell references

    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    if(bc_configs[ 0] == PRECURSOR || bc_configs[ 1] == PRECURSOR) {
      // working array
      real *uu_WE = (real*) malloc(dom[dev].Gfx.jnb*dom[dev].Gfx.knb
        * sizeof(real));

      // copy from device to host
      (cudaMemcpy(uu_WE, _u_WE[dev], sizeof(real)
        * dom[dev].Gfx.jnb*dom[dev].Gfx.knb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
        for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = j + k * Dom.Gfx.jnb;
          CC = jj + kk * dom[dev].Gfx.jnb;
          u_WE[C] = uu_WE[CC];
        }
      }

      // clean up
      free(uu_WE);
    }
    if(bc_configs[ 2] == PRECURSOR || bc_configs[ 3] == PRECURSOR) {
      // working array
      real *uu_SN = (real*) malloc(dom[dev].Gfx.inb*dom[dev].Gfx.knb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(uu_SN, _u_SN[dev], sizeof(real)
        * dom[dev].Gfx.inb*dom[dev].Gfx.knb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + k * Dom.Gfx.inb;
          CC = ii + kk * dom[dev].Gfx.inb;
          u_SN[C] = uu_SN[CC];
        }
      }

      // clean up
      free(uu_SN);
    }
    if(bc_configs[ 4] == PRECURSOR || bc_configs[ 5] == PRECURSOR) {
      // working array
      real *uu_BT = (real*) malloc(dom[dev].Gfx.inb*dom[dev].Gfx.jnb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(uu_BT, _u_BT[dev], sizeof(real)
        * dom[dev].Gfx.inb*dom[dev].Gfx.jnb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          C = i + j * Dom.Gfx.inb;
          CC = ii + jj * dom[dev].Gfx.inb;
          u_BT[C] = uu_BT[CC];
        }
      }

      // clean up
      free(uu_BT);
    }

    if(bc_configs[ 6] == PRECURSOR || bc_configs[ 7] == PRECURSOR) {
      // working array
      real *vv_WE = (real*) malloc(dom[dev].Gfy.jnb*dom[dev].Gfy.knb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(vv_WE, _v_WE[dev], sizeof(real)
        * dom[dev].Gfy.jnb*dom[dev].Gfy.knb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
        for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = j + k * Dom.Gfy.jnb;
          CC = jj + kk * dom[dev].Gfy.jnb;
          v_WE[C] = vv_WE[CC];
        }
      }

      // clean up
      free(vv_WE);
    }
    if(bc_configs[ 8] == PRECURSOR || bc_configs[ 9] == PRECURSOR) {
      // working array
      real *vv_SN = (real*) malloc(dom[dev].Gfy.inb*dom[dev].Gfy.knb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(vv_SN, _v_SN[dev], sizeof(real)
        * dom[dev].Gfy.inb*dom[dev].Gfy.knb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + k * Dom.Gfy.inb;
          CC = ii + kk * dom[dev].Gfy.inb;
          v_SN[C] = vv_SN[CC];
        }
      }

      // clean up
      free(vv_SN);
    }
    if(bc_configs[10] == PRECURSOR || bc_configs[11] == PRECURSOR) {
      // working array
      real *vv_BT = (real*) malloc(dom[dev].Gfy.inb*dom[dev].Gfy.jnb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(vv_BT, _v_BT[dev], sizeof(real)
        * dom[dev].Gfy.inb*dom[dev].Gfy.jnb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          C = i + j * Dom.Gfy.inb;
          CC = ii + jj * dom[dev].Gfy.inb;
          v_BT[C] = vv_BT[CC];
        }
      }

      // clean up
      free(vv_BT);
    }

    if(bc_configs[12] == PRECURSOR || bc_configs[13] == PRECURSOR) {
      // working array
      real *ww_WE = (real*) malloc(dom[dev].Gfz.jnb*dom[dev].Gfz.knb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(ww_WE, _w_WE[dev], sizeof(real)
        * dom[dev].Gfz.jnb*dom[dev].Gfz.knb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
        for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = j + k * Dom.Gfz.jnb;
          CC = jj + kk * dom[dev].Gfz.jnb;
          w_WE[C] = ww_WE[CC];
        }
      }

      // clean up
      free(ww_WE);
    }
    if(bc_configs[14] == PRECURSOR || bc_configs[15] == PRECURSOR) {
      // working array
      real *ww_SN = (real*) malloc(dom[dev].Gfz.inb*dom[dev].Gfz.knb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(ww_SN, _w_SN[dev], sizeof(real)
        * dom[dev].Gfz.inb*dom[dev].Gfz.knb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + k * Dom.Gfz.inb;
          CC = ii + kk * dom[dev].Gfz.inb;
          w_SN[C] = ww_SN[CC];
        }
      }

      // clean up
      free(ww_SN);
    }
    if(bc_configs[16] == PRECURSOR || bc_configs[17] == PRECURSOR) {
      // working array
      real *ww_BT = (real*) malloc(dom[dev].Gfz.inb*dom[dev].Gfz.jnb
        * sizeof(real));

      // copy from host to device
      (cudaMemcpy(ww_BT, _w_BT[dev], sizeof(real)
        * dom[dev].Gfz.inb*dom[dev].Gfz.jnb, cudaMemcpyDeviceToHost));

      // select appropriate subdomain
      // set position
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          C = i + j * Dom.Gfz.inb;
          CC = ii + jj * dom[dev].Gfz.inb;
          w_BT[C] = ww_BT[CC];
        }
      }

      // clean up
      free(ww_BT);
    }
  }
}

extern "C"
void cuda_dom_free(void)
{
  // free device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    (cudaFree(_dom[dev]));
    (cudaFree(_p0[dev]));
    (cudaFree(_p[dev]));
    (cudaFree(_phi[dev]));
    //(cudaFree(_divU[dev]));
    (cudaFree(_u[dev]));
    (cudaFree(_v[dev]));
    (cudaFree(_w[dev]));
    (cudaFree(_u0[dev]));
    (cudaFree(_v0[dev]));
    (cudaFree(_w0[dev]));
    (cudaFree(_f_x[dev]));
    (cudaFree(_f_y[dev]));
    (cudaFree(_f_z[dev]));
#ifndef IMPLICIT
    (cudaFree(_diff0_u[dev]));
    (cudaFree(_diff0_v[dev]));
    (cudaFree(_diff0_w[dev]));
#endif
    (cudaFree(_conv0_u[dev]));
    (cudaFree(_conv0_v[dev]));
    (cudaFree(_conv0_w[dev]));
    (cudaFree(_diff_u[dev]));
    (cudaFree(_diff_v[dev]));
    (cudaFree(_diff_w[dev]));
    (cudaFree(_conv_u[dev]));
    (cudaFree(_conv_v[dev]));
    (cudaFree(_conv_w[dev]));
    (cudaFree(_u_star[dev]));
    (cudaFree(_v_star[dev]));
    (cudaFree(_w_star[dev]));
    (cudaFree(_u_WE[dev]));
    (cudaFree(_u_SN[dev]));
    (cudaFree(_u_BT[dev]));
    (cudaFree(_v_WE[dev]));
    (cudaFree(_v_SN[dev]));
    (cudaFree(_v_BT[dev]));
    (cudaFree(_w_WE[dev]));
    (cudaFree(_w_SN[dev]));
    (cudaFree(_w_BT[dev]));
    (cudaFree(_rhs_p[dev]));
  }

  // free device memory on host
  free(_dom);
  free(_p0);
  free(_p);
  free(_phi);
  //free(_divU);
  free(_u);
  free(_v);
  free(_w);
  free(_u0);
  free(_v0);
  free(_w0);
  free(_f_x);
  free(_f_y);
  free(_f_z);
#ifndef IMPLICIT
  free(_diff0_u);
  free(_diff0_v);
  free(_diff0_w);
#endif
  free(_conv0_u);
  free(_conv0_v);
  free(_conv0_w);
  free(_diff_u);
  free(_diff_v);
  free(_diff_w);
  free(_conv_u);
  free(_conv_v);
  free(_conv_w);
  free(_u_star);
  free(_v_star);
  free(_w_star);
  free(_u_WE);
  free(_u_SN);
  free(_u_BT);
  free(_v_WE);
  free(_v_SN);
  free(_v_BT);
  free(_w_WE);
  free(_w_SN);
  free(_w_BT);
  free(_rhs_p);
}

extern "C"
void cuda_dom_BC(void)
{
  // CPU threading for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // check whether each subdomain boundary (E, W, N, S, T, B) is
    // an external boundary
    if(dom[dev].W == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // u-velocity
      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);

      dim3 dimBlocks_u(threads_y, threads_z);
      dim3 numBlocks_u(blocks_y, blocks_z);

      // v-velocity
      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);

      dim3 dimBlocks_v(threads_y, threads_z);
      dim3 numBlocks_v(blocks_y, blocks_z);

      // w-velocity
      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);

      dim3 dimBlocks_w(threads_y, threads_z);
      dim3 numBlocks_w(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.pW) {
        case PERIODIC:
          BC_p_W_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_W_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
      switch(bc.uW) {
        case PERIODIC:
          BC_u_W_P<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_W_D<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], bc.uWD);
          break;
        case NEUMANN:
          BC_u_W_N<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_W_T<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_WE[dev]);
          break;
      }
      switch(bc.vW) {
        case PERIODIC:
          BC_v_W_P<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_W_D<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], bc.vWD);
          break;
        case NEUMANN:
          BC_v_W_N<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_W_T<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_WE[dev]);
          break;
      }
      switch(bc.wW) {
        case PERIODIC:
          BC_w_W_P<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_W_D<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], bc.wWD);
          break;
        case NEUMANN:
          BC_w_W_N<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_W_T<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_WE[dev]);
          break;
      }
    }
    if(dom[dev].E == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // u-velocity
      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);

      dim3 dimBlocks_u(threads_y, threads_z);
      dim3 numBlocks_u(blocks_y, blocks_z);

      // v-velocity
      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);

      dim3 dimBlocks_v(threads_y, threads_z);
      dim3 numBlocks_v(blocks_y, blocks_z);

      // w-velocity
      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);

      dim3 dimBlocks_w(threads_y, threads_z);
      dim3 numBlocks_w(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.pE) {
        case PERIODIC:
          BC_p_E_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_E_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
      switch(bc.uE) {
        case PERIODIC:
          BC_u_E_P<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_E_D<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], bc.uED);
          break;
        case NEUMANN:
          BC_u_E_N<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_E_T<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_WE[dev]);
          break;
      }
      switch(bc.vE) {
        case PERIODIC:
          BC_v_E_P<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_E_D<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], bc.vED);
          break;
        case NEUMANN:
          BC_v_E_N<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_E_T<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_WE[dev]);
          break;
      }
      switch(bc.wE) {
        case PERIODIC:
          BC_w_E_P<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_E_D<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], bc.wED);
          break;
        case NEUMANN:
          BC_w_E_N<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_E_T<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_WE[dev]);
          break;
      }
    }
    if(dom[dev].S == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

      // u-velocity
      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);

      dim3 dimBlocks_u(threads_z, threads_x);
      dim3 numBlocks_u(blocks_z, blocks_x);

      // v-velocity
      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);

      dim3 dimBlocks_v(threads_z, threads_x);
      dim3 numBlocks_v(blocks_z, blocks_x);

      // w-velocity
      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);

      dim3 dimBlocks_w(threads_z, threads_x);
      dim3 numBlocks_w(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.pS) {
        case PERIODIC:
          BC_p_S_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_S_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
      switch(bc.uS) {
        case PERIODIC:
          BC_u_S_P<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_S_D<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], bc.uSD);
          break;
        case NEUMANN:
          BC_u_S_N<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_S_T<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_SN[dev]);
          break;
      }
      switch(bc.vS) {
        case PERIODIC:
          BC_v_S_P<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_S_D<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], bc.vSD);
          break;
        case NEUMANN:
          BC_v_S_N<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_S_T<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_SN[dev]);
          break;
      }
      switch(bc.wS) {
        case PERIODIC:
          BC_w_S_P<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_S_D<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], bc.wSD);
          break;
        case NEUMANN:
          BC_w_S_N<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_S_T<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_SN[dev]);
          break;
      }
    }
    if(dom[dev].N == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

      // u-velocity
      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);

      dim3 dimBlocks_u(threads_z, threads_x);
      dim3 numBlocks_u(blocks_z, blocks_x);

      // v-velocity
      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);

      dim3 dimBlocks_v(threads_z, threads_x);
      dim3 numBlocks_v(blocks_z, blocks_x);

      // w-velocity
      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);

      dim3 dimBlocks_w(threads_z, threads_x);
      dim3 numBlocks_w(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.pN) {
        case PERIODIC:
          BC_p_N_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_N_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
      switch(bc.uN) {
        case PERIODIC:
          BC_u_N_P<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_N_D<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], bc.uND);
          break;
        case NEUMANN:
          BC_u_N_N<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_N_T<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_SN[dev]);
          break;
      }
      switch(bc.vN) {
        case PERIODIC:
          BC_v_N_P<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_N_D<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], bc.vND);
          break;
        case NEUMANN:
          BC_v_N_N<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_N_T<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_SN[dev]);
          break;
      }
      switch(bc.wN) {
        case PERIODIC:
          BC_w_N_P<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_N_D<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], bc.wND);
          break;
        case NEUMANN:
          BC_w_N_N<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_N_T<<<numBlocks_w, dimBlocks_w>>>(_u[dev], _dom[dev], _w_SN[dev]);
          break;
      }
    }
    if(dom[dev].B == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);

      // u-velocity
      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);

      dim3 dimBlocks_u(threads_x, threads_y);
      dim3 numBlocks_u(blocks_x, blocks_y);

      // v-velocity
      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);

      dim3 dimBlocks_v(threads_x, threads_y);
      dim3 numBlocks_v(blocks_x, blocks_y);

      // w-velocity
      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);

      dim3 dimBlocks_w(threads_x, threads_y);
      dim3 numBlocks_w(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.pB) {
        case PERIODIC:
          BC_p_B_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_B_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
      switch(bc.uB) {
        case PERIODIC:
          BC_u_B_P<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_B_D<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], bc.uBD);
          break;
        case NEUMANN:
          BC_u_B_N<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_B_T<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_BT[dev]);
          break;
      }
      switch(bc.vB) {
        case PERIODIC:
          BC_v_B_P<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_B_D<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], bc.vBD);
          break;
        case NEUMANN:
          BC_v_B_N<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_B_T<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_BT[dev]);
          break;
      }
      switch(bc.wB) {
        case PERIODIC:
          BC_w_B_P<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_B_D<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], bc.wBD);
          break;
        case NEUMANN:
          BC_w_B_N<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_B_T<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_BT[dev]);
          break;
      }
    }
    if(dom[dev].T == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);

      // u-velocity
      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);

      dim3 dimBlocks_u(threads_x, threads_y);
      dim3 numBlocks_u(blocks_x, blocks_y);

      // v-velocity
      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);

      dim3 dimBlocks_v(threads_x, threads_y);
      dim3 numBlocks_v(blocks_x, blocks_y);

      // w-velocity
      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);

      dim3 dimBlocks_w(threads_x, threads_y);
      dim3 numBlocks_w(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.pT) {
        case PERIODIC:
          BC_p_T_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_T_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
      switch(bc.uT) {
        case PERIODIC:
          BC_u_T_P<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_T_D<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], bc.uTD);
          break;
        case NEUMANN:
          BC_u_T_N<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_T_T<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_BT[dev]);
          break;
      }
      switch(bc.vT) {
        case PERIODIC:
          BC_v_T_P<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_T_D<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], bc.vTD);
          break;
        case NEUMANN:
          BC_v_T_N<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_T_T<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_BT[dev]);
          break;
      }
      switch(bc.wT) {
        case PERIODIC:
          BC_w_T_P<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_T_D<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], bc.wTD);
          break;
        case NEUMANN:
          BC_w_T_N<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_T_T<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_BT[dev]);
          break;
      }
    }
  }
}

extern "C"
void cuda_dom_BC_star(void)
{
  // CPU threading for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // check whether each subdomain boundary (E, W, N, S, T, B) is
    // an external boundary
    if(dom[dev].W == -1) {
      // set up kernel call

      // u-velocity
      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);

      dim3 dimBlocks_u(threads_y, threads_z);
      dim3 numBlocks_u(blocks_y, blocks_z);

      // v-velocity
      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);

      dim3 dimBlocks_v(threads_y, threads_z);
      dim3 numBlocks_v(blocks_y, blocks_z);

      // w-velocity
      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);

      dim3 dimBlocks_w(threads_y, threads_z);
      dim3 numBlocks_w(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.uW) {
        case PERIODIC:
          BC_u_W_P<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_W_D<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], bc.uWD);
          break;
        case NEUMANN:
          BC_u_W_N<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_W_T<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev],
            _u_WE[dev]);
          break;
      }
      switch(bc.vW) {
        case PERIODIC:
          BC_v_W_P<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_W_D<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev], bc.vWD);
          break;
        case NEUMANN:
          BC_v_W_N<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_W_T<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev],
            _v_WE[dev]);
          break;
      }
      switch(bc.wW) {
        case PERIODIC:
          BC_w_W_P<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_W_D<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], bc.wWD);
          break;
        case NEUMANN:
          BC_w_W_N<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_W_T<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev],
            _w_WE[dev]);
          break;
      }
    }
    if(dom[dev].E == -1) {
      // set up kernel call

      // u-velocity
      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);

      dim3 dimBlocks_u(threads_y, threads_z);
      dim3 numBlocks_u(blocks_y, blocks_z);

      // v-velocity
      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);

      dim3 dimBlocks_v(threads_y, threads_z);
      dim3 numBlocks_v(blocks_y, blocks_z);

      // w-velocity
      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);

      dim3 dimBlocks_w(threads_y, threads_z);
      dim3 numBlocks_w(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.uE) {
        case PERIODIC:
          BC_u_E_P<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_E_D<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], bc.uED);
          break;
        case NEUMANN:
          BC_u_E_N<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_E_T<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev],
            _u_WE[dev]);
          break;
      }
      switch(bc.vE) {
        case PERIODIC:
          BC_v_E_P<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_E_D<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev], bc.vED);
          break;
        case NEUMANN:
          BC_v_E_N<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_E_T<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev],
            _v_WE[dev]);
          break;
      }
      switch(bc.wE) {
        case PERIODIC:
          BC_w_E_P<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_E_D<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], bc.wED);
          break;
        case NEUMANN:
          BC_w_E_N<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_E_T<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], 
            _w_WE[dev]);
          break;
      }
    }
    if(dom[dev].S == -1) {
      // set up kernel call

      // u-velocity
      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);

      dim3 dimBlocks_u(threads_z, threads_x);
      dim3 numBlocks_u(blocks_z, blocks_x);

      // v-velocity
      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);

      dim3 dimBlocks_v(threads_z, threads_x);
      dim3 numBlocks_v(blocks_z, blocks_x);

      // w-velocity
      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);

      dim3 dimBlocks_w(threads_z, threads_x);
      dim3 numBlocks_w(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.uS) {
        case PERIODIC:
          BC_u_S_P<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_S_D<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], bc.uSD);
          break;
        case NEUMANN:
          BC_u_S_N<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_S_T<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev],
            _u_SN[dev]);
          break;
      }
      switch(bc.vS) {
        case PERIODIC:
          BC_v_S_P<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_S_D<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev], bc.vSD);
          break;
        case NEUMANN:
          BC_v_S_N<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_S_T<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev],
            _v_SN[dev]);
          break;
      }
      switch(bc.wS) {
        case PERIODIC:
          BC_w_S_P<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_S_D<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], bc.wSD);
          break;
        case NEUMANN:
          BC_w_S_N<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_S_T<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev],
            _w_SN[dev]);
          break;
      }
    }
    if(dom[dev].N == -1) {
      // set up kernel call

      // u-velocity
      if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);

      dim3 dimBlocks_u(threads_z, threads_x);
      dim3 numBlocks_u(blocks_z, blocks_x);

      // v-velocity
      if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);

      dim3 dimBlocks_v(threads_z, threads_x);
      dim3 numBlocks_v(blocks_z, blocks_x);

      // w-velocity
      if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfz.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);

      dim3 dimBlocks_w(threads_z, threads_x);
      dim3 numBlocks_w(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.uN) {
        case PERIODIC:
          BC_u_N_P<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_N_D<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], bc.uND);
          break;
        case NEUMANN:
          BC_u_N_N<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_N_T<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], 
            _u_SN[dev]);
          break;
      }
      switch(bc.vN) {
        case PERIODIC:
          BC_v_N_P<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_N_D<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev], bc.vND);
          break;
        case NEUMANN:
          BC_v_N_N<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_S_T<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev],
            _v_SN[dev]);
          break;
      }
      switch(bc.wN) {
        case PERIODIC:
          BC_w_N_P<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_N_D<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], bc.wND);
          break;
        case NEUMANN:
          BC_w_N_N<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_N_T<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev],
            _w_SN[dev]);
          break;
      }
    }
    if(dom[dev].B == -1) {
      // set up kernel call

      // u-velocity
      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);

      dim3 dimBlocks_u(threads_x, threads_y);
      dim3 numBlocks_u(blocks_x, blocks_y);

      // v-velocity
      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);

      dim3 dimBlocks_v(threads_x, threads_y);
      dim3 numBlocks_v(blocks_x, blocks_y);

      // w-velocity
      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);

      dim3 dimBlocks_w(threads_x, threads_y);
      dim3 numBlocks_w(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.uB) {
        case PERIODIC:
          BC_u_B_P<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_B_D<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], bc.uBD);
          break;
        case NEUMANN:
          BC_u_B_N<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_B_T<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev],
            _u_BT[dev]);
          break;
      }
      switch(bc.vB) {
        case PERIODIC:
          BC_v_B_P<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_B_D<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev], bc.vBD);
          break;
        case NEUMANN:
          BC_v_B_N<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_B_T<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev],
            _v_BT[dev]);
          break;
      }
      switch(bc.wB) {
        case PERIODIC:
          BC_w_B_P<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_B_D<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], bc.wBD);
          break;
        case NEUMANN:
          BC_w_B_N<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_B_T<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev],
            _w_BT[dev]);
          break;
      }
    }
    if(dom[dev].T == -1) {
      // set up kernel call

      // u-velocity
      if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfx.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);

      dim3 dimBlocks_u(threads_x, threads_y);
      dim3 numBlocks_u(blocks_x, blocks_y);

      // v-velocity
      if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfy.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);

      dim3 dimBlocks_v(threads_x, threads_y);
      dim3 numBlocks_v(blocks_x, blocks_y);

      // w-velocity
      if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);

      dim3 dimBlocks_w(threads_x, threads_y);
      dim3 numBlocks_w(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.uT) {
        case PERIODIC:
          BC_u_T_P<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_u_T_D<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev], bc.uTD);
          break;
        case NEUMANN:
          BC_u_T_N<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_u_T_T<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _dom[dev],
            _u_BT[dev]);
          break;
      }
      switch(bc.vT) {
        case PERIODIC:
          BC_v_T_P<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_v_T_D<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev], bc.vTD);
          break;
        case NEUMANN:
          BC_v_T_N<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_v_T_T<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _dom[dev],
            _v_BT[dev]);
          break;
      }
      switch(bc.wT) {
        case PERIODIC:
          BC_w_T_P<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_w_T_D<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev], bc.wTD);
          break;
        case NEUMANN:
          BC_w_T_N<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev]);
          break;
        case PRECURSOR:
          BC_w_T_T<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _dom[dev],
            _w_BT[dev]);
          break;
      }
    }
  }
}

#ifndef IMPLICIT
extern "C"
void cuda_U_star_2(void)
{
  // CPU thread
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // u-component
    if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx.jnb + 2;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx.knb + 2;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) (threads_y-2));
    blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) (threads_z-2));

    dim3 dimBlocks_u(threads_y, threads_z);
    dim3 numBlocks_u(blocks_y, blocks_z);

    u_star_2<<<numBlocks_u, dimBlocks_u>>>(rho_f, nu,
      _u0[dev], _v0[dev], _w0[dev], _p0[dev], _f_x[dev],
      _diff0_u[dev], _conv0_u[dev], _diff_u[dev], _conv_u[dev],
      _u_star[dev], _dom[dev], dt0, dt, _phase[dev]);

    // v-component
    if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy.knb + 2;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy.inb + 2;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) (threads_z-2));
    blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) (threads_x-2));

    dim3 dimBlocks_v(threads_z, threads_x);
    dim3 numBlocks_v(blocks_z, blocks_x);

    v_star_2<<<numBlocks_v, dimBlocks_v>>>(rho_f, nu,
      _u0[dev], _v0[dev], _w0[dev], _p0[dev], _f_y[dev],
      _diff0_v[dev], _conv0_v[dev], _diff_v[dev], _conv_v[dev],
      _v_star[dev], _dom[dev], dt0, dt, _phase[dev]);

    // w-component
    if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz.inb + 2;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz.jnb + 2;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) (threads_x-2));
    blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) (threads_y-2));

    dim3 dimBlocks_w(threads_x, threads_y);
    dim3 numBlocks_w(blocks_x, blocks_y);

    w_star_2<<<numBlocks_w, dimBlocks_w>>>(rho_f, nu,
      _u0[dev], _v0[dev], _w0[dev], _p0[dev], _f_z[dev],
      _diff0_w[dev], _conv0_w[dev], _diff_w[dev], _conv_w[dev],
      _w_star[dev], _dom[dev], dt0, dt, _phase[dev]);

  #ifdef TEST
    // copy _u_star back to _u
    cudaMemcpy(_u[dev], _u_star[dev], dom[dev].Gfx.s3b * sizeof(real),
      cudaMemcpyDeviceToDevice);

    // copy _v_star back to _v
    cudaMemcpy(_v[dev], _v_star[dev], dom[dev].Gfy.s3b * sizeof(real),
      cudaMemcpyDeviceToDevice);

    // copy _w_star back to _w
    cudaMemcpy(_w[dev], _w_star[dev], dom[dev].Gfz.s3b * sizeof(real),
      cudaMemcpyDeviceToDevice);

  #endif
  }
}
#endif

extern "C"
void cuda_dom_BC_phi(void)
{
  // CPU threading for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // check whether each subdomain boundary (E, W, N, S, T, B) is
    // an external boundary
    if(dom[dev].W == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.pW) {
        case PERIODIC:
          BC_p_W_P<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_W_N<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].E == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.pE) {
        case PERIODIC:
          BC_p_E_P<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_E_N<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].S == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.pS) {
        case PERIODIC:
          BC_p_S_P<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_S_N<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].N == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.pN) {
        case PERIODIC:
          BC_p_N_P<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_N_N<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].B == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.pB) {
        case PERIODIC:
          BC_p_B_P<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_B_N<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
      }

    }
    if(dom[dev].T == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.pT) {
        case PERIODIC:
          BC_p_T_P<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_T_N<<<numBlocks_p, dimBlocks_p>>>(_phi[dev], _dom[dev]);
          break;
      }
    }
  }
}

extern "C"
void cuda_dom_BC_p(void)
{
  // CPU threading for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // check whether each subdomain boundary (E, W, N, S, T, B) is
    // an external boundary
    if(dom[dev].W == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.pW) {
        case PERIODIC:
          BC_p_W_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_W_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].E == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(bc.pE) {
        case PERIODIC:
          BC_p_E_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_E_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].S == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.pS) {
        case PERIODIC:
          BC_p_S_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_S_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].N == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

      // apply BC to all fields for this face
      switch(bc.pN) {
        case PERIODIC:
          BC_p_N_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_N_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
    }
    if(dom[dev].B == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.pB) {
        case PERIODIC:
          BC_p_B_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_B_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }

    }
    if(dom[dev].T == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);

      // apply BC to all fields for this face
      switch(bc.pT) {
        case PERIODIC:
          BC_p_T_P<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_p_T_N<<<numBlocks_p, dimBlocks_p>>>(_p[dev], _dom[dev]);
          break;
      }
    }
  }
}

extern "C"
void cuda_project(void)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // solve for u
    if(dom[dev].Gfx._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx._kn / (real) threads_z);

    dim3 dimBlocks_u(threads_y, threads_z);
    dim3 numBlocks_u(blocks_y, blocks_z);

    project_u<<<numBlocks_u, dimBlocks_u>>>(_u_star[dev], _phi[dev],
      rho_f, dt, _u[dev], _dom[dev], 1. / dom[dev].dx, _flag_u[dev], _phase[dev]);

    // solve for v
    if(dom[dev].Gfy._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy._kn;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy._in;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy._kn / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy._in / (real) threads_x);

    dim3 dimBlocks_v(threads_z, threads_x);
    dim3 numBlocks_v(blocks_z, blocks_x);

    project_v<<<numBlocks_v, dimBlocks_v>>>(_v_star[dev], _phi[dev],
      rho_f, dt, _v[dev], _dom[dev], 1. / dom[dev].dy, _flag_v[dev], _phase[dev]);

    // solve for w
    if(dom[dev].Gfz._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz._jn;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz._jn / (real) threads_y);

    dim3 dimBlocks_w(threads_x, threads_y);
    dim3 numBlocks_w(blocks_x, blocks_y);

    project_w<<<numBlocks_w, dimBlocks_w>>>(_w_star[dev], _phi[dev],
      rho_f, dt, _w[dev], _dom[dev], 1. / dom[dev].dz, _flag_w[dev], _phase[dev]);
  }
}

extern "C"
real cuda_find_dt(void)
{
  // results from all devices
  real *dts = (real*) malloc(nsubdom * sizeof(real));
    // cpumem += nsubdom * sizeof(real);

  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    // search
    real u_max = find_max_mag(dom[dev].Gfx.s3b, _u[dev]);
    real v_max = find_max_mag(dom[dev].Gfy.s3b, _v[dev]);
    real w_max = find_max_mag(dom[dev].Gfz.s3b, _w[dev]);

#ifndef IMPLICIT
    real tmp = u_max / dom[dev].dx + 2.*nu/dom[dev].dx/dom[dev].dx;
    tmp += v_max / dom[dev].dy + 2.*nu/dom[dev].dy/dom[dev].dy;
    tmp += w_max / dom[dev].dz + 2.*nu/dom[dev].dz/dom[dev].dz;

    dts[dev] = tmp;
#else
    real tmp = u_max / dom[dev].dx;
    tmp += v_max / dom[dev].dy;
    tmp += w_max / dom[dev].dz;

    dts[dev] = tmp;
#endif

/*
#ifndef IMPLICIT
    real max = u_max / dom[dev].dx + 2.*nu/dom[dev].dx/dom[dev].dx;
    real tmp = v_max / dom[dev].dy + 2.*nu/dom[dev].dy/dom[dev].dy;
    if(tmp > max) max = tmp;
    tmp = w_max / dom[dev].dz + 2.*nu/dom[dev].dz/dom[dev].dz;
    if(tmp > max) max = tmp;

    dts[dev] = max;
#else
    real max = u_max / dom[dev].dx;
    real tmp = v_max / dom[dev].dy;
    if(tmp > max) max = tmp;
    tmp = w_max / dom[dev].dz;
    if(tmp > max) max = tmp;

    dts[dev] = max;
    //dts[dev] += v_max / dom[dev].dy;
    //dts[dev] += w_max / dom[dev].dz;
#endif
*/
    dts[dev] = CFL / dts[dev];
  }

  // find min of all devices
  real min = FLT_MAX;
  for(int i = 0; i < nsubdom; i++)
    if(dts[i] < min) min = dts[i];

  // clean up
  free(dts);

#ifdef IMPLICIT
  if(min > 1.5*dt) min = 1.5*dt;
#endif

  return min;
}

extern "C"
void cuda_store_u(void)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    (cudaMemcpy(_conv0_u[dev], _conv_u[dev],
      dom[dev].Gfx.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    (cudaMemcpy(_conv0_v[dev], _conv_v[dev],
      dom[dev].Gfy.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    (cudaMemcpy(_conv0_w[dev], _conv_w[dev],
      dom[dev].Gfz.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
#ifndef IMPLICIT
    (cudaMemcpy(_diff0_u[dev], _diff_u[dev],
      dom[dev].Gfx.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    (cudaMemcpy(_diff0_v[dev], _diff_v[dev],
      dom[dev].Gfy.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    (cudaMemcpy(_diff0_w[dev], _diff_w[dev],
      dom[dev].Gfz.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
#endif

    (cudaMemcpy(_p0[dev], _p[dev],
      dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));

    (cudaMemcpy(_u0[dev], _u[dev],
      dom[dev].Gfx.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    (cudaMemcpy(_v0[dev], _v[dev],
      dom[dev].Gfy.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    (cudaMemcpy(_w0[dev], _w[dev],
      dom[dev].Gfz.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
  }
}

extern "C"
void cuda_update_p()
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    int threads_y = 0;
    int threads_z = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // update pressure
    if(dom[dev].Gcc._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gcc._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gcc._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gcc._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gcc._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gcc._kn / (real) threads_z);

    dim3 dimBlocks_p(threads_y, threads_z);
    dim3 numBlocks_p(blocks_y, blocks_z);

    // create temporary working array
    real *_Lp;
    (cudaMalloc((void**) &_Lp,
      sizeof(real)*dom[dev].Gcc.s3b));

    update_p_laplacian<<<numBlocks_p, dimBlocks_p>>>(_Lp, _phi[dev], _dom[dev]);

    update_p<<<numBlocks_p, dimBlocks_p>>>(_Lp, _p0[dev], _p[dev], _phi[dev],
      _dom[dev], nu, dt, _phase[dev]);

    // clean up temporary array
    (cudaFree(_Lp));
  }
}

extern "C"
void cuda_compute_forcing(real *pid_int, real *pid_back, real Kp, real Ki,
  real Kd)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // x-component
    if(dom[dev].Gfx._jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx._jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx._knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx._knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx._jnb / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx._knb / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);

    // y-component
    if(dom[dev].Gfy._knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy._knb;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy._inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy._inb;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy._knb / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy._inb / (real) threads_x);

    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);

    // z-component
    if(dom[dev].Gfz._inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz._inb;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz._jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz._jnb;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz._inb / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz._jnb / (real) threads_y);

    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // reset forcing arrays
    forcing_reset_x<<<numBlocks_x, dimBlocks_x>>>(_f_x[dev], _dom[dev]);
    forcing_reset_y<<<numBlocks_y, dimBlocks_y>>>(_f_y[dev], _dom[dev]);
    forcing_reset_z<<<numBlocks_z, dimBlocks_z>>>(_f_z[dev], _dom[dev]);

    // linearly accelerate pressure gradient from zero
    real delta = ttime - p_tDelay;
    if (delta >= 0 ) {
      if(gradP.xa == 0) gradP.x = gradP.xm;
      else if(fabs(delta*gradP.xa) > fabs(gradP.xm)) gradP.x = gradP.xm;
      else gradP.x = delta*gradP.xa;

      if(gradP.ya == 0) gradP.y = gradP.ym;
      else if(fabs(delta*gradP.ya) > fabs(gradP.ym)) gradP.y = gradP.ym;
      else gradP.y = delta*gradP.ya;

      // turn off if PID controller is being used
      if(!(Kp > 0 || Ki > 0 || Kd > 0)) {
        if(gradP.za == 0) gradP.z = gradP.zm;
        else if(fabs(delta*gradP.za) > fabs(gradP.zm)) gradP.z = gradP.zm;
        else gradP.z = delta*gradP.za;
      }
    }

    // linearly accelerate gravitational acceleration from zero
    delta = ttime - g_tDelay;
    if (delta >= 0) {
      if(g.xa == 0) g.x = g.xm;
      else if(fabs(delta*g.xa) > fabs(g.xm)) g.x = g.xm;
      else g.x = delta*g.xa;

      if(g.ya == 0) g.y = g.ym;
      else if(fabs(delta*g.ya) > fabs(g.ym)) g.y = g.ym;
      else g.y = delta*g.ya;

      if(g.za == 0) g.z = g.zm;
      else if(fabs(delta*g.za) > fabs(g.zm)) g.z = g.zm;
      else g.z = delta*g.za;
    }

    delta = ttime - p_tDelay;
    // PID controller  TODO: make this into a kernel function
    if (delta >= 0) {
      if(Kp > 0 || Ki > 0 || Kd > 0) {
        cuda_part_pull();
        real acc_z = 0;
        real volp = 0;
        real massp = 0;
        for(int i = 0; i < nparts; i++) {
          volp += parts[i].r*parts[i].r*parts[i].r;
          massp += parts[i].rho*parts[i].r*parts[i].r*parts[i].r;
          acc_z += parts[i].wdot;
        }
        volp *= 4./3.*PI;
        massp *= 4./3.*PI;
        real volfrac = volp / (Dom.xl * Dom.yl * Dom.zl);
        real rho_avg = massp/volp*volfrac + rho_f*(1.-volfrac);
        acc_z /= nparts;

        *pid_int = *pid_int + acc_z*dt;
        gradP.z = gradP.z
          + (Kp*acc_z + Ki*(*pid_int)/ttime + (Kd)*(acc_z-*pid_back))*rho_avg;
        *pid_back = acc_z;
      }
    }
    forcing_add_x_const<<<numBlocks_x, dimBlocks_x>>>(-gradP.x / rho_f,
      _f_x[dev], _dom[dev]);
    forcing_add_y_const<<<numBlocks_y, dimBlocks_y>>>(-gradP.y / rho_f,
      _f_y[dev], _dom[dev]);
    forcing_add_z_const<<<numBlocks_z, dimBlocks_z>>>(-gradP.z / rho_f,
      _f_z[dev], _dom[dev]);

  }
}

extern "C"
void cuda_compute_turb_forcing(void)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    // x-component
    if(dom[dev].Gfx._jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx._jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx._knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx._knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx._jnb / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx._knb / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);

    // y-component
    if(dom[dev].Gfy._knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy._knb;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy._inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy._inb;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy._knb / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy._inb / (real) threads_x);

    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);

    // z-component
    if(dom[dev].Gfz._inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz._inb;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz._jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz._jnb;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz._inb / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz._jnb / (real) threads_y);

    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    real k = cuda_compute_energy();
    real k0 = 13.5 * turbA*turbA * 0.0361 * turbl*turbl;
    //printf("\nk = %e k0 = %e k0/k = %e\n", k, k0, k0/k);
    //k = 1;
    //k0 = 1;

    // do not reset forcing arrays
    // turbulence linear forcing A*(k0/k)*u'
    // create temporary workspace
    real *_u_co;  // colocated u-velocity workspace
    real *_v_co;  // colocated v-velocity workspace
    real *_w_co;  // colocated w-velocity workspace
    (cudaMalloc((void**) &_u_co, sizeof(real)*dom[dev].Gcc.s3));
    (cudaMalloc((void**) &_v_co, sizeof(real)*dom[dev].Gcc.s3));
    (cudaMalloc((void**) &_w_co, sizeof(real)*dom[dev].Gcc.s3));

    // colocate the velocity variables for computation
    cuda_colocate_Gfx(_u[dev], _u_co, _dom[dev]);
    cuda_colocate_Gfy(_v[dev], _v_co, _dom[dev]);
    cuda_colocate_Gfz(_w[dev], _w_co, _dom[dev]);

    real umean = avg_entries(dom[dev].Gcc.s3, _u_co);
    real vmean = avg_entries(dom[dev].Gcc.s3, _v_co);
    real wmean = avg_entries(dom[dev].Gcc.s3, _w_co);

    // clean up workspace
    (cudaFree(_u_co));
    (cudaFree(_v_co));
    (cudaFree(_w_co));

    // now add in the forcing
    // add in U then remove the mean to apply the perturbation forcing
    forcing_add_x_field<<<numBlocks_x, dimBlocks_x>>>(turbA*k0/k, _u[dev],
      _f_x[dev], _dom[dev], _phase[dev]);
    forcing_add_y_field<<<numBlocks_y, dimBlocks_y>>>(turbA*k0/k, _v[dev],
      _f_y[dev], _dom[dev], _phase[dev]);
    forcing_add_z_field<<<numBlocks_z, dimBlocks_z>>>(turbA*k0/k, _w[dev],
      _f_z[dev], _dom[dev], _phase[dev]);

    forcing_add_x_const<<<numBlocks_x, dimBlocks_x>>>(-turbA*k0/k*umean,
      _f_x[dev], _dom[dev]);
    forcing_add_y_const<<<numBlocks_y, dimBlocks_y>>>(-turbA*k0/k*vmean,
      _f_y[dev], _dom[dev]);
    forcing_add_z_const<<<numBlocks_z, dimBlocks_z>>>(-turbA*k0/k*wmean,
      _f_z[dev], _dom[dev]);
  }
}

extern "C"
real cuda_compute_energy(void)
{
  real *k_dom = (real*) malloc(sizeof(real)*nsubdom);
  real k = 0; // energy to return

  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    // create temporary workspace
    real *_u_co;  // colocated u-velocity workspace
    real *_v_co;  // colocated v-velocity workspace
    real *_w_co;  // colocated w-velocity workspace
    (cudaMalloc((void**) &_u_co, sizeof(real)*dom[dev].Gcc.s3));
    (cudaMalloc((void**) &_v_co, sizeof(real)*dom[dev].Gcc.s3));
    (cudaMalloc((void**) &_w_co, sizeof(real)*dom[dev].Gcc.s3));

    // colocate the velocity variables for computation
    cuda_colocate_Gfx(_u[dev], _u_co, _dom[dev]);
    cuda_colocate_Gfy(_v[dev], _v_co, _dom[dev]);
    cuda_colocate_Gfz(_w[dev], _w_co, _dom[dev]);

    /**** MAY NEED TO SUBTRACT MEAN FLOW ****/

    // do computation
    int threads_y = 0;
    int threads_z = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gcc._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gcc._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gcc._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gcc._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gcc._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gcc._kn / (real) threads_z);

    dim3 dimBlocks(threads_y, threads_z);
    dim3 numBlocks(blocks_y, blocks_z);
   
    // do the u*u multiplication and drop the result back into u
    energy_multiply<<<numBlocks, dimBlocks>>>(_u_co, _v_co, _w_co, _u_co,
      _dom[dev]);
    // do the averaging
    k_dom[dev] = 0.5 * avg_entries(dom[dev].Gcc.s3, _u_co);

    // clean up workspace
    (cudaFree(_u_co));
    (cudaFree(_v_co));
    (cudaFree(_w_co));
  }

  // TODO: rework for multi-gpu. For now, just use the only value we have
  k = k_dom[0];

  // clean up
  free(k_dom);

  return k;
}

extern "C"
void cuda_colocate_Gfx(real *_u, real *_u_co, dom_struct *_dom)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_y = 0;
    int threads_z = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gfx._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx._kn / (real) threads_z);

    dim3 dimBlocks(threads_y, threads_z);
    dim3 numBlocks(blocks_y, blocks_z);
    
    colocate_Gfx<<<numBlocks, dimBlocks>>>(_u, _u_co, _dom);
  }
}

extern "C"
void cuda_colocate_Gfy(real *_v, real *_v_co, dom_struct *_dom)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_z = 0;
    int threads_x = 0;
    int blocks_z = 0;
    int blocks_x = 0;

    if(dom[dev].Gfy._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy._kn;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy._in;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy._kn / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy._in / (real) threads_x);

    dim3 dimBlocks(threads_z, threads_x);
    dim3 numBlocks(blocks_z, blocks_x);
    
    colocate_Gfy<<<numBlocks, dimBlocks>>>(_v, _v_co, _dom);
  }
}

extern "C"
void cuda_colocate_Gfz(real *_w, real *_w_co, dom_struct *_dom)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int blocks_x = 0;
    int blocks_y = 0;

    if(dom[dev].Gfz._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz._jn;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz._jn / (real) threads_y);

    dim3 dimBlocks(threads_x, threads_y);
    dim3 numBlocks(blocks_x, blocks_y);
    
    colocate_Gfz<<<numBlocks, dimBlocks>>>(_w, _w_co, _dom);
  }
}

// NOTE: this needs to be reworked for multi-gpu usage
extern "C"
void cuda_solvability(void)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int Nx = dom[dev].Gcc.jn * dom[dev].Gcc.kn;
    int Ny = dom[dev].Gcc.in * dom[dev].Gcc.kn;
    int Nz = dom[dev].Gcc.in * dom[dev].Gcc.jn;
    real eps_x = 0;
    real eps_y = 0;
    real eps_z = 0;

    int threads_x;
    int threads_y;
    int threads_z;
    int blocks_x;
    int blocks_y;
    int blocks_z;

    // x-component
    if(dom[dev].Gfx._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx._kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);

    // y-component
    if(dom[dev].Gfy._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy._kn;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy._in;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy._kn / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy._in / (real) threads_x);

    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);

    // z-component
    if(dom[dev].Gfz._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz._jn;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz._jn / (real) threads_y);

    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // create temporary storage for reduction algorithms
    real *u_star_red;
    real *v_star_red;
    real *w_star_red;
    (cudaMalloc((void**) &u_star_red, sizeof(real) * 2 * Nx));
    gpumem += sizeof(real) * 2 * Nx;
    (cudaMalloc((void**) &v_star_red, sizeof(real) * 2 * Ny));
    gpumem += sizeof(real) * 2 * Ny;
    (cudaMalloc((void**) &w_star_red, sizeof(real) * 2 * Nz));
    gpumem += sizeof(real) * 2 * Nz;

    // calculate x-face surface integrals
    surf_int_x_copy<<<numBlocks_x, dimBlocks_x>>>(_u_star[dev],
      u_star_red, _dom[dev]);
    eps_x = sum_entries(2 * Nx, u_star_red);
    eps_x = eps_x * dom[dev].dy*dom[dev].dz;
    // calculate y-face surface integrals
    surf_int_y_copy<<<numBlocks_y, dimBlocks_y>>>(_v_star[dev],
      v_star_red, _dom[dev]);
    eps_y = sum_entries(2 * Ny, v_star_red);
    eps_y = eps_y * dom[dev].dx*dom[dev].dz;
    // calculate z-face surface integrals
    surf_int_z_copy<<<numBlocks_z, dimBlocks_z>>>(_w_star[dev],
      w_star_red, _dom[dev]);
    eps_z = sum_entries(2 * Nz, w_star_red);
    eps_z = eps_z * dom[dev].dx*dom[dev].dy;

    // subtract eps from outflow plane
    switch(out_plane) {
      case WEST:
        // normalize eps by face surface area
        eps_x = eps_x/dom[dev].yl/dom[dev].zl;
        //eps_y = eps_y/dom[dev].yl/dom[dev].zl;
        //eps_z = eps_z/dom[dev].yl/dom[dev].zl;
        plane_eps_x_W<<<numBlocks_x, dimBlocks_x>>>
          //(eps_x+eps_y+eps_z, _u_star[dev], _dom[dev]);
          (eps_x, _u_star[dev], _dom[dev]);
        break;
      case EAST:
        // normalize eps by face surface area
        eps_x = eps_x/dom[dev].yl/dom[dev].zl;
        //eps_y = eps_y/dom[dev].yl/dom[dev].zl;
        //eps_z = eps_z/dom[dev].yl/dom[dev].zl;
        plane_eps_x_E<<<numBlocks_x, dimBlocks_x>>>
          //(eps_x+eps_y+eps_z, _u_star[dev], _dom[dev]);
          (eps_x, _u_star[dev], _dom[dev]);
        break;
      case SOUTH:
        // normalize eps by face surface area
        //eps_x = eps_x/dom[dev].zl/dom[dev].xl;
        eps_y = eps_y/dom[dev].zl/dom[dev].xl;
        //eps_z = eps_z/dom[dev].zl/dom[dev].xl;
        plane_eps_y_S<<<numBlocks_y, dimBlocks_y>>>
          //(eps_x+eps_y+eps_z, _v_star[dev], _dom[dev]);
          (eps_y, _v_star[dev], _dom[dev]);
        break;
      case NORTH:
        // normalize eps by face surface area
        //eps_x = eps_x/dom[dev].zl/dom[dev].xl;
        eps_y = eps_y/dom[dev].zl/dom[dev].xl;
        //eps_z = eps_z/dom[dev].zl/dom[dev].xl;
        plane_eps_y_N<<<numBlocks_y, dimBlocks_y>>>
          //(eps_x+eps_y+eps_z, _v_star[dev], _dom[dev]);
          (eps_y, _v_star[dev], _dom[dev]);
        break;
      case BOTTOM:
        // normalize eps by face surface area
        //eps_x = eps_x/dom[dev].xl/dom[dev].yl;
        //eps_y = eps_y/dom[dev].xl/dom[dev].yl;
        eps_z = eps_z/dom[dev].xl/dom[dev].yl;
        plane_eps_z_B<<<numBlocks_z, dimBlocks_z>>>
          //(eps_x+eps_y+eps_z, _w_star[dev], _dom[dev]);
          (eps_z, _w_star[dev], _dom[dev]);
        break;
      case TOP:
        // normalize eps by face surface area
        //eps_x = eps_x/dom[dev].xl/dom[dev].yl;
        //eps_y = eps_y/dom[dev].xl/dom[dev].yl;
        eps_z = eps_z/dom[dev].xl/dom[dev].yl;
        plane_eps_z_T<<<numBlocks_z, dimBlocks_z>>>
          //(eps_x+eps_y+eps_z, _w_star[dev], _dom[dev]);
          (eps_z, _w_star[dev], _dom[dev]);
        break;
      case HOMOGENEOUS:
        // spread the errors over the entire domain
        eps_x = 0.5*eps_x/dom[dev].yl/dom[dev].zl;
        eps_y = 0.5*eps_y/dom[dev].zl/dom[dev].xl;
        eps_z = 0.5*eps_z/dom[dev].xl/dom[dev].yl;
        plane_eps_x_W<<<numBlocks_x, dimBlocks_x>>>
          (eps_x, _u_star[dev], _dom[dev]);
        plane_eps_x_E<<<numBlocks_x, dimBlocks_x>>>
          (eps_x, _u_star[dev], _dom[dev]);
        plane_eps_y_S<<<numBlocks_y, dimBlocks_y>>>
          (eps_y, _v_star[dev], _dom[dev]);
        plane_eps_y_N<<<numBlocks_y, dimBlocks_y>>>
          (eps_y, _v_star[dev], _dom[dev]);
        plane_eps_z_B<<<numBlocks_z, dimBlocks_z>>>
          (eps_z, _w_star[dev], _dom[dev]);
        plane_eps_z_T<<<numBlocks_z, dimBlocks_z>>>
          (eps_z, _w_star[dev], _dom[dev]);
        break;
      // no default because it should already be guaranteed that these are
      // the only options as checked by the input file reader
    }

    // clean up
    (cudaFree(u_star_red));
    (cudaFree(v_star_red));
    (cudaFree(w_star_red));
  }
}

// TODO: MULTI-GPU capabilities must be reworked
void cuda_move_parts_sub()
{
   // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) nparts / (real) threads);
    if(threads > nparts) {
      threads = nparts;
      blocks = 1;
    }

    dim3 dimBlocks(threads);
    dim3 numBlocks(blocks);

    if(nparts > 0) {
      real eps = 0.01;

      if(nparts == 1) {
        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);
        move_parts_a<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev], nparts,
          dt, dt0, g, gradP, rho_f, ttime);
        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);
      } else if(nparts > 1) {
        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);

        int nBins = binDom.Gcc.s3;

        // initialize threads for nBin size
        int threads_nb = MAX_THREADS_1D;
        int blocks_nb = (int)ceil((real) nBins / (real) threads_nb);
        if(threads_nb > nBins) {
          threads_nb = nBins;
          blocks_nb = 1;
        }
        dim3 dimBlocks_nb(threads_nb);
        dim3 numBlocks_nb(blocks_nb);

        /* go to each particle and find its bin */
        int *_partInd;
        int *_partBin;
        (cudaMalloc((void**) &_partInd, nparts*sizeof(int)));
        (cudaMalloc((void**) &_partBin, nparts*sizeof(int)));
        gpumem += nparts*sizeof(int);
        gpumem += nparts*sizeof(int);
        bin_fill<<<numBlocks, dimBlocks>>>(_partInd, _partBin, nparts,
          _parts[dev], _binDom, bc);

        /* sort by bin */
        thrust::device_ptr<int> ptr_partBin(_partBin);
        thrust::device_ptr<int> ptr_partInd(_partInd);
        thrust::sort_by_key(ptr_partBin, ptr_partBin + nparts, ptr_partInd);
        _partBin = thrust::raw_pointer_cast(ptr_partBin);
        _partInd = thrust::raw_pointer_cast(ptr_partInd);

        /* calculate start and end index of each bin */
        int *_binStart;
        int *_binEnd;
        (cudaMalloc((void**) &_binStart, nBins*sizeof(int)));
        (cudaMalloc((void**) &_binEnd, nBins*sizeof(int)));
        init<<<numBlocks_nb, dimBlocks_nb>>>(_binStart, nBins, -1);
        init<<<numBlocks_nb, dimBlocks_nb>>>(_binEnd, nBins, -1);
        gpumem += nBins*sizeof(int);
        gpumem += nBins*sizeof(int);

        int smemSize = sizeof(int)*(threads + 1);
        bin_start<<<blocks, threads, smemSize>>>(_binStart, _binEnd, _partBin,
          nparts);

        /* count the number of particles in each bin */
        //int *_binCount;
        //(cudaMalloc((void**) &_binCount, nBins*sizeof(int)));
        //init<<<numBlocks_nb, dimBlocks_nb>>>(_binCount, nBins, 0);
        //gpumem += nBins*sizeof(int);
        //bin_partCount<<<dimBlocks_nb, numBlocks_nb>>>(_binCount,_binStart,
        //  _binEnd, _binDom, bc, nBins);

        // launch a thread per particle to calc collision
        collision_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts,
         _dom[dev], eps, mu, bc, _binStart, _binEnd, _partBin, _partInd, 
         _binDom, interactionLength);

        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);
        move_parts_a<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev], nparts,
          dt, dt0, g, gradP, rho_f, ttime);

        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        
        collision_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts,
         _dom[dev], eps, mu, bc, _binStart, _binEnd, _partBin, _partInd, 
         _binDom, interactionLength);

        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);

        // free variables
        (cudaFree(_partInd));
        (cudaFree(_partBin));
        //(cudaFree(_binCount));
        (cudaFree(_binStart));
        (cudaFree(_binEnd));
      }
    }
  }
}

void cuda_move_parts()
{
   // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) nparts / (real) threads);
    if(threads > nparts) {
      threads = nparts;
      blocks = 1;
    }

    dim3 dimBlocks(threads);
    dim3 numBlocks(blocks);

    if(nparts > 0) {
      real eps = 0.01;

      if(nparts == 1) {
        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);
        move_parts_a<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev], nparts,
          dt, dt0, g, gradP, rho_f, ttime);
        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);
      } else if(nparts > 1) {
        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);

        int nBins = binDom.Gcc.s3;

        // initialize threads for nBin size
        int threads_nb = MAX_THREADS_1D;
        int blocks_nb = (int)ceil((real) nBins / (real) threads_nb);
        if(threads_nb > nBins) {
          threads_nb = nBins;
          blocks_nb = 1;
        }
        dim3 dimBlocks_nb(threads_nb);
        dim3 numBlocks_nb(blocks_nb);

        /* go to each particle and find its bin */
        int *_partInd;
        int *_partBin;
        (cudaMalloc((void**) &_partInd, nparts*sizeof(int)));
        (cudaMalloc((void**) &_partBin, nparts*sizeof(int)));
        gpumem += nparts*sizeof(int);
        gpumem += nparts*sizeof(int);
        bin_fill<<<numBlocks, dimBlocks>>>(_partInd, _partBin, nparts,
          _parts[dev], _binDom, bc);

        /* sort by bin */
        thrust::device_ptr<int> ptr_partBin(_partBin);
        thrust::device_ptr<int> ptr_partInd(_partInd);
        thrust::sort_by_key(ptr_partBin, ptr_partBin + nparts, ptr_partInd);
        _partBin = thrust::raw_pointer_cast(ptr_partBin);
        _partInd = thrust::raw_pointer_cast(ptr_partInd);

        /* calculate start and end index of each bin */
        int *_binStart;
        int *_binEnd;
        (cudaMalloc((void**) &_binStart, nBins*sizeof(int)));
        (cudaMalloc((void**) &_binEnd, nBins*sizeof(int)));
        init<<<numBlocks_nb, dimBlocks_nb>>>(_binStart, nBins, -1);
        init<<<numBlocks_nb, dimBlocks_nb>>>(_binEnd, nBins, -1);
        gpumem += nBins*sizeof(int);
        gpumem += nBins*sizeof(int);

        int smemSize = sizeof(int)*(threads + 1);
        bin_start<<<blocks, threads, smemSize>>>(_binStart, _binEnd,_partBin,
          nparts);

        /* count the number of particles in each bin */
        //int *_binCount;
        //(cudaMalloc((void**) &_binCount, nBins*sizeof(int)));
        //init<<<numBlocks_nb, dimBlocks_nb>>>(_binCount, nBins, 0);
        //gpumem += nBins*sizeof(int);
        //bin_partCount<<<dimBlocks_nb, numBlocks_nb>>>(_binCount,_binStart,
        //  _binEnd, _binDom, bc, nBins);

        // launch a thread per particle to calc collision
        collision_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts,
         _dom[dev], eps, mu, bc, _binStart, _binEnd, _partBin, _partInd, 
         _binDom, interactionLength);

        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);
        move_parts_a<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev], nparts,
          dt, dt0, g, gradP, rho_f, ttime);

        collision_init<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        
        collision_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts,
         _dom[dev], eps, mu, bc, _binStart, _binEnd, _partBin, _partInd, 
         _binDom, interactionLength);

        spring_parts<<<numBlocks, dimBlocks>>>(_parts[dev], nparts);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev],
          nparts, bc, eps, mu, rho_f, nu, interactionLength, dt);

        // free variables
        (cudaFree(_partInd));
        (cudaFree(_partBin));
        //(cudaFree(_binCount));
        (cudaFree(_binStart));
        (cudaFree(_binEnd));
      }

      move_parts_b<<<numBlocks, dimBlocks>>>(_dom[dev], _parts[dev], nparts,
        dt, dt0, g, gradP, rho_f, ttime);
    }
  }
}

void cuda_yank_turb_planes(int *bc_flow_configs, real *pos, real *vel)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    /** X-direction flow: WE-plane **/
    // set up kernel call
    // u-velocity
    if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx.knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);

    dim3 dimBlocks_u(threads_y, threads_z);
    dim3 numBlocks_u(blocks_y, blocks_z);

    // yank plane from precursor domain
    if(bc_flow_configs[ 0] == PRECURSOR)
      yank_u_WE<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_WE[dev],
        pos[ 0], vel[ 0]);
    else if(bc_flow_configs[ 1] == PRECURSOR)
      yank_u_WE<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_WE[dev],
        pos[ 0], vel[ 1]);

    // set up kernel call
    // v-velocity
    if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfy.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy.knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);

    dim3 dimBlocks_v(threads_y, threads_z);
    dim3 numBlocks_v(blocks_y, blocks_z);

    // yank plane from precursor domain
    if(bc_flow_configs[ 6] == PRECURSOR)
      yank_v_WE<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_WE[dev],
        pos[ 0], vel[ 0]);
    if(bc_flow_configs[ 7] == PRECURSOR)
      yank_v_WE<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_WE[dev],
        pos[ 0], vel[ 1]);

    // set up kernel call
    // w-velocity
    if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfz.knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);

    dim3 dimBlocks_w(threads_y, threads_z);
    dim3 numBlocks_w(blocks_y, blocks_z);

    // yank plane from precursor domain
    if(bc_flow_configs[12] == PRECURSOR)
      yank_w_WE<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_WE[dev],
        pos[ 0], vel[ 0]);
    if(bc_flow_configs[13] == PRECURSOR)
      yank_w_WE<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_WE[dev],
        pos[ 0], vel[ 1]);

    /** Y-direction flow: SN-plane **/
    // set up kernel call
    // u-velocity
    if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx.knb;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfx.inb;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);

    dimBlocks_u.x = threads_z;
    dimBlocks_u.y = threads_x;
    numBlocks_u.x = blocks_z;
    numBlocks_u.y = blocks_x;

    // yank plane from precursor domain
    if(bc_flow_configs[ 2] == PRECURSOR)
      yank_u_SN<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_SN[dev],
        pos[ 4], vel[ 8]);
    else if(bc_flow_configs[ 3] == PRECURSOR)
      yank_u_SN<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_SN[dev],
        pos[ 4], vel[ 9]);

    // set up kernel call
    // v-velocity
    if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy.knb;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy.inb;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);

    dimBlocks_v.x = threads_z;
    dimBlocks_v.y = threads_x;
    numBlocks_v.x = blocks_z;
    numBlocks_v.y = blocks_x;

    // yank plane from precursor domain
    if(bc_flow_configs[ 8] == PRECURSOR)
      yank_v_SN<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_SN[dev],
        pos[ 4], vel[ 8]);
    if(bc_flow_configs[ 9] == PRECURSOR)
      yank_v_SN<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_SN[dev],
        pos[ 4], vel[ 9]);

    // set up kernel call
    // w-velocity
    if(dom[dev].Gfz.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfz.knb;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz.inb;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);

    dimBlocks_w.x = threads_z;
    dimBlocks_w.y = threads_x;
    numBlocks_w.x = blocks_z;
    numBlocks_w.y = blocks_x;

    // yank plane from precursor domain
    if(bc_flow_configs[14] == PRECURSOR)
      yank_w_SN<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_SN[dev],
        pos[ 4], vel[ 8]);
    if(bc_flow_configs[15] == PRECURSOR)
      yank_w_SN<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_SN[dev],
        pos[ 4], vel[ 9]);

    /** Z-direction flow: BT-plane **/
    // set up kernel call
    // u-velocity
    if(dom[dev].Gfx.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfx.inb;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);

    dimBlocks_u.x = threads_x;
    dimBlocks_u.y = threads_y;
    numBlocks_u.x = blocks_x;
    numBlocks_u.y = blocks_y;

    // yank plane from precursor domain
    if(bc_flow_configs[ 4] == PRECURSOR)
      yank_u_BT<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_BT[dev],
        pos[ 8], vel[16]);
    else if(bc_flow_configs[ 5] == PRECURSOR)
      yank_u_BT<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _dom[dev], _u_BT[dev],
        pos[ 8], vel[17]);

    // set up kernel call
    // v-velocity
    if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy.inb;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfy.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfy.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);

    dimBlocks_v.x = threads_x;
    dimBlocks_v.y = threads_y;
    numBlocks_v.x = blocks_x;
    numBlocks_v.y = blocks_y;

    // yank plane from precursor domain
    if(bc_flow_configs[10] == PRECURSOR)
      yank_v_BT<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_BT[dev],
        pos[ 8], vel[16]);
    if(bc_flow_configs[11] == PRECURSOR)
      yank_v_BT<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _dom[dev], _v_BT[dev],
        pos[ 8], vel[17]);

    // set up kernel call
    // w-velocity
    if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz.inb;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);

    dimBlocks_w.x = threads_x;
    dimBlocks_w.y = threads_y;
    numBlocks_w.x = blocks_x;
    numBlocks_w.y = blocks_y;

    // yank plane from precursor domain
    if(bc_flow_configs[16] == PRECURSOR)
      yank_w_BT<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_BT[dev],
        pos[ 8], vel[16]);
    if(bc_flow_configs[17] == PRECURSOR)
      yank_w_BT<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _dom[dev], _w_BT[dev],
        pos[ 8], vel[17]);
  }
}

void cuda_parts_internal(void)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(nparts > 0) {

      // solve for u
      if(dom[dev].Gfx._jn < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfx._jn;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gfx._kn < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfx._kn;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gfx._jn / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gfx._kn / (real) threads_z);

      dim3 dimBlocks_u(threads_y, threads_z);
      dim3 numBlocks_u(blocks_y, blocks_z);

      internal_u<<<numBlocks_u, dimBlocks_u>>>(_u[dev], _parts[dev], _dom[dev],
        _flag_u[dev], _phase[dev]);

      // solve for v
      if(dom[dev].Gfy._kn < MAX_THREADS_DIM)
        threads_z = dom[dev].Gfy._kn;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gfy._in < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfy._in;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gfy._kn / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gfy._in / (real) threads_x);

      dim3 dimBlocks_v(threads_z, threads_x);
      dim3 numBlocks_v(blocks_z, blocks_x);

      internal_v<<<numBlocks_v, dimBlocks_v>>>(_v[dev], _parts[dev], _dom[dev],
        _flag_v[dev], _phase[dev]);

      // solve for w
      if(dom[dev].Gfz._in < MAX_THREADS_DIM)
        threads_x = dom[dev].Gfz._in;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gfz._jn < MAX_THREADS_DIM)
        threads_y = dom[dev].Gfz._jn;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gfz._in / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gfz._jn / (real) threads_y);

      dim3 dimBlocks_w(threads_x, threads_y);
      dim3 numBlocks_w(blocks_x, blocks_y);

      internal_w<<<numBlocks_w, dimBlocks_w>>>(_w[dev], _parts[dev], _dom[dev],
        _flag_w[dev], _phase[dev]);
    }
  }
}
