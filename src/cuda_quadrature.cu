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
#include "cuda_particle.h"
#include "cuda_quadrature.h"

#include <cuda.h>

extern "C"
void cuda_quad_check_nodes(int dev,
  real *node_t, real *node_p, int nnodes)
{
  int threads = nnodes;
  int blocks = nparts;

  dim3 dimBlocks(threads);
  dim3 numBlocks(blocks);

  if(nparts > 0)
    check_nodes<<<numBlocks, dimBlocks>>>(nparts, _parts[dev], _dom[dev],
      node_t, node_p, nnodes, bc);

}

extern "C"
void cuda_quad_interp(int dev,
  real *node_t, real *node_p, int nnodes,
  real *pp, real *ur, real *ut, real *up)
{
  int threads = nnodes;
  int blocks = nparts;

  dim3 dimBlocks(threads);
  dim3 numBlocks(blocks);

  // TODO: interpolate using shared memory
  if(nparts > 0)
    interpolate_nodes<<<numBlocks, dimBlocks>>>(_p0[dev], _p[dev],
      _u[dev], _v[dev], _w[dev], rho_f, nu, gradP,
      _parts[dev], _dom[dev], node_t, node_p, nnodes, pp, ur, ut, up, dt0, dt, bc);
}

extern "C"
void cuda_Lamb(void)
{
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    (cudaSetDevice(dev + dev_start));

    int i;  // iterator

    // set up coefficient table
    int nn[15] = {0,
                  1, 1,
                  2, 2, 2,
                  3, 3, 3, 3,
                  4, 4, 4, 4, 4};
    int mm[15] = {0,
                  0, 1,
                  0, 1, 2,
                  0, 1, 2, 3,
                  0, 1, 2, 3, 4};

    // set up quadrature nodes for 7th-order Lebedev quadrature
    real PI14 = 0.25 * PI;
    real PI12 = 0.5 * PI;
    real PI34 = 0.75 * PI;
    real PI54 = 1.25 * PI;
    real PI32 = 1.5 * PI;
    real PI74 = 1.75 * PI;
    real alph1 = 0.955316618124509; //54.736
    real alph2 = 2.186276035465284; //125.264
  /*
    real alph3 = 0.440510663004698; //25.239
    real alph4 = 2.701081990585095; //154.761
    real alph5 = 1.264518957625227; //72.452
    real alph6 = 1.877073695964566; //107.548
    real alph7 = 1.249045772398254; //71.565
    real alph8 = 1.892546881191539; //108.435
    real alph9 = 0.321750554396642; //18.435
    real alph10 = 2.819842099193151; //161.565
  */

    // weights
    real A1 = 0.598398600683775;//0.159572960182342;//0.119679720136760;//
    real A2 = 0.478718880547015;//0.283685262546377;//0.403919055461543;//
    real A3 = 0.403919055461543;//0.265071880146639;//0.403919055461543;
    real B = 0.;//0.253505610897313;//0.359039160410267;

    // nodes TODO: find a more elegant way of fixing the divide by sin(0)
    // TODO: put this in GPU constant memory
    real a1_t[6] = {PI12, PI12, PI12, PI12, 0.+DIV_ST, PI-DIV_ST};
    real a1_p[6] = {0., PI12, PI, PI32, 0., 0.};
    real a2_t[12] = {PI12, PI12, PI12, PI12,
                     PI14, PI14, PI14, PI14,
                     PI34, PI34, PI34, PI34};
    real a2_p[12] = {PI14, PI34, PI54, PI74,
                     0., PI12, PI, PI32,
                     0., PI12, PI, PI32};
    real a3_t[8] = {alph1, alph1, alph1, alph1,
                    alph2, alph2, alph2, alph2};
    real a3_p[8] = {PI14, PI34, PI54, PI74,
                    PI14, PI34, PI54, PI74};
    /*real b_t[24] = {alph3, alph4, alph3, alph4,
                    alph3, alph4, alph3, alph4,
                    alph5, alph5, alph6, alph6,
                    alph5, alph5, alph6, alph6,
                    alph5, alph5, alph6, alph6,
                    alph5, alph5, alph6, alph6};
    real b_p[24] = {PI14, PI14, PI74, PI74,
                    PI34, PI32, PI54, PI54,
                    alph7, -alph7, alph7, -alph7,
                    alph8, -alph8, alph8, -alph8,
                    alph9, alph10, alph9, alph10,
                    -alph9, -alph10, -alph9, -alph10};
  */

    int nnodes = 26;

    // put all quadrature nodes together for interpolation
    real node_t[nnodes];
    real node_p[nnodes];
    for(i = 0; i < 6; i++) {
      node_t[i] = a1_t[i];
      node_p[i] = a1_p[i];
    }
    for(i = 0; i < 12; i++) {
      node_t[6+i] = a2_t[i];
      node_p[6+i] = a2_p[i];
    }
    for(i = 0; i < 8; i++) {
      node_t[18+i] = a3_t[i];
      node_p[18+i] = a3_p[i];
    }
    //for(i = 0; i < 24; i++) {
      //node_t[26+i] = b_t[i];
      //node_p[26+i] = b_p[i];
    //}

    // create a place to temporarily store field variables at quadrature nodes
    int *_nn;
    int *_mm;
    (cudaMalloc((void**) &_nn, nnodes * sizeof(int)));
    gpumem += nnodes * sizeof(int);
    (cudaMalloc((void**) &_mm, nnodes * sizeof(int)));
    gpumem += nnodes * sizeof(int);
    (cudaMemcpy(_nn, nn, nnodes * sizeof(int),
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_mm, mm, nnodes * sizeof(int),
      cudaMemcpyHostToDevice));
    real *_node_t;
    real *_node_p;
    (cudaMalloc((void**) &_node_t, nnodes * sizeof(real)));
    gpumem += nnodes * sizeof(real);
    (cudaMalloc((void**) &_node_p, nnodes * sizeof(real)));
    gpumem += nnodes * sizeof(real);
    (cudaMemcpy(_node_t, node_t, nnodes * sizeof(real),
      cudaMemcpyHostToDevice));
    (cudaMemcpy(_node_p, node_p, nnodes * sizeof(real),
      cudaMemcpyHostToDevice));
    real *_pp;
    real *_ur;
    real *_ut;
    real *_up;
    (cudaMalloc((void**) &_pp, nnodes * nparts * sizeof(real)));
    gpumem += nnodes * nparts * sizeof(real);
    (cudaMalloc((void**) &_ur, nnodes * nparts * sizeof(real)));
    gpumem += nnodes * nparts * sizeof(real);
    (cudaMalloc((void**) &_ut, nnodes * nparts * sizeof(real)));
    gpumem += nnodes * nparts * sizeof(real);
    (cudaMalloc((void**) &_up, nnodes * nparts * sizeof(real)));
    gpumem += nnodes * nparts * sizeof(real);

    // interpolate field variables to quadrature nodes
    cuda_quad_check_nodes(dev, _node_t, _node_p, nnodes);
    cuda_quad_interp(dev, _node_t, _node_p, nnodes, _pp, _ur, _ut, _up);

/*
    real *pp = (real*) malloc(nnodes*nparts*sizeof(real));
    cpumem += nnodes*nparts*sizeof(real);
    real *ur = (real*) malloc(nnodes*nparts*sizeof(real));
    cpumem += nnodes*nparts*sizeof(real);
    real *ut = (real*) malloc(nnodes*nparts*sizeof(real));
    cpumem += nnodes*nparts*sizeof(real);
    real *up = (real*) malloc(nnodes*nparts*sizeof(real));
    cpumem += nnodes*nparts*sizeof(real);
    (cudaMemcpy(pp, _pp, nnodes*nparts*sizeof(real),
      cudaMemcpyDeviceToHost));
    (cudaMemcpy(ur, _ur, nnodes*nparts*sizeof(real),
      cudaMemcpyDeviceToHost));
    (cudaMemcpy(ut, _ut, nnodes*nparts*sizeof(real),
      cudaMemcpyDeviceToHost));
    (cudaMemcpy(up, _up, nnodes*nparts*sizeof(real),
      cudaMemcpyDeviceToHost));

    printf("\n");
    for(i = 0; i < nnodes*nparts; i++) {
      printf("%d: pp = %f    ur = %f    ut = %f    up = %f\n", i, pp[i], ur[i], ut[i], up[i]);
    }

    free(pp);
    free(ur);
    free(ut);
    free(up);
*/

    // create storage using maximum particle coefficient size
    // TODO: make this a constant value instead of calculating it every time
    int ncoeffs = 0;
    for(i = 0; i < nparts; i++) {
      if(parts[i].ncoeff > ncoeffs)
        ncoeffs = parts[i].ncoeff;
    }

    // create temporary storage for inner product integrands
    real *int_Yp_re;
    real *int_Yp_im;
    real *int_rDYu_re;
    real *int_rDYu_im;
    real *int_xXDYu_re;
    real *int_xXDYu_im;
    (cudaMalloc((void**) &int_Yp_re,
      nparts * nnodes * ncoeffs * sizeof(real)));
    gpumem += nparts * nnodes * ncoeffs * sizeof(real);
    (cudaMalloc((void**) &int_Yp_im,
      nparts * nnodes * ncoeffs * sizeof(real)));
    gpumem += nparts * nnodes * ncoeffs * sizeof(real);
    (cudaMalloc((void**) &int_rDYu_re,
      nparts * nnodes * ncoeffs * sizeof(real)));
    gpumem += nparts * nnodes * ncoeffs * sizeof(real);
    (cudaMalloc((void**) &int_rDYu_im,
      nparts * nnodes * ncoeffs * sizeof(real)));
    gpumem += nparts * nnodes * ncoeffs * sizeof(real);
    (cudaMalloc((void**) &int_xXDYu_re,
      nparts * nnodes * ncoeffs * sizeof(real)));
    gpumem += nparts * nnodes * ncoeffs * sizeof(real);
    (cudaMalloc((void**) &int_xXDYu_im,
      nparts * nnodes * ncoeffs * sizeof(real)));
    gpumem += nparts * nnodes * ncoeffs * sizeof(real);

    // apply quadrature rules to compute the three (six real) inner products
    // and set values of Lamb's coefficients
    dim3 dimBlocks(nnodes);
    dim3 numBlocks(nparts, ncoeffs);
    if(nparts > 0) {
      cuda_get_coeffs<<<numBlocks, dimBlocks>>>(_parts[dev],
        _nn, _mm, _node_t, _node_p,
        _pp, _ur, _ut, _up, mu, nu,
        coeff_stride, _pnm_re[dev], _pnm_im[dev],
        _phinm_re[dev], _phinm_im[dev],
        _chinm_re[dev], _chinm_im[dev],
        int_Yp_re, int_Yp_im,
        int_rDYu_re, int_rDYu_im,
        int_xXDYu_re, int_xXDYu_im,
        nnodes, ncoeffs, A1, A2, A3, B,
        _pnm_re0[dev], _pnm_im0[dev],
        _phinm_re0[dev], _phinm_im0[dev],
        _chinm_re0[dev], _chinm_im0[dev],
        lamb_relax);

      int threads = MAX_THREADS_1D;
      int blocks = (int)ceil((real) nparts / (real) threads);

      dim3 dimBlocks_f(threads);
      dim3 numBlocks_f(blocks);
      cuda_calc_forces<<<numBlocks_f, dimBlocks_f>>>(_dom[dev], _parts[dev],
        nparts, gradP,
        rho_f, mu, nu, coeff_stride,
        _pnm_re[dev], _pnm_im[dev],
        _phinm_re[dev], _phinm_im[dev],
        _chinm_re[dev], _chinm_im[dev]);
    }

    // clean up temporary variables
    (cudaFree(_nn));
    (cudaFree(_mm));
    (cudaFree(_pp));
    (cudaFree(_ur));
    (cudaFree(_ut));
    (cudaFree(_up));
    (cudaFree(_node_t));
    (cudaFree(_node_p));
    (cudaFree(int_Yp_re));
    (cudaFree(int_Yp_im));
    (cudaFree(int_rDYu_re));
    (cudaFree(int_rDYu_im));
    (cudaFree(int_xXDYu_re));
    (cudaFree(int_xXDYu_im));
  }
}

extern "C"
real cuda_lamb_err(void)
{
  real *errors = (real*) malloc(nsubdom * sizeof(real));
  // cpumem += nsubdom * sizeof(real);
  real max = FLT_MIN;

  if(nparts > 0) {
    #pragma omp parallel num_threads(nsubdom)
    {
      int dev = omp_get_thread_num();
      (cudaSetDevice(dev + dev_start));

      // create a place to store sorted coefficients and errors
      real *part_errors = (real*) malloc(nparts * sizeof(real));
      // cpumem += nparts * sizeof(real); 
      real *_sorted_coeffs;
      real *_sorted_errors;
      real *_part_errors;
      (cudaMalloc((void**) &_sorted_coeffs,
        nparts*6*coeff_stride*sizeof(real)));
      gpumem += 6 * nparts * coeff_stride * sizeof(real);
      (cudaMalloc((void**) &_sorted_errors,
        nparts*6*coeff_stride*sizeof(real)));
      gpumem += 6 * nparts * coeff_stride * sizeof(real);
      (cudaMalloc((void**) &_part_errors,
        nparts*sizeof(real)));
      gpumem += nparts * sizeof(real);
      
      // sort the coefficients and calculate errors along the way
      dim3 dimBlocks(1);
      dim3 numBlocks(nparts);

      compute_error<<<numBlocks, dimBlocks>>>(lamb_cut, coeff_stride, nparts,
        _pnm_re[dev], _pnm_re0[dev], _pnm_im[dev], _pnm_im0[dev],
        _phinm_re[dev], _phinm_re0[dev], _phinm_im[dev], _phinm_im0[dev],
        _chinm_re[dev], _chinm_re0[dev], _chinm_im[dev], _chinm_im0[dev],
        _sorted_coeffs, _sorted_errors, _part_errors, _dom[dev], nu);

      // copy the errors back to device
      (cudaMemcpy(part_errors, _part_errors,
        nparts*sizeof(real), cudaMemcpyDeviceToHost));

      // find maximum error of all particles
      real tmp = FLT_MIN;
      for(int i = 0; i < nparts; i++) {
        if(part_errors[i] > tmp) tmp = part_errors[i];
      }
      errors[dev] = tmp;

      // clean up
      (cudaFree(_sorted_coeffs));
      (cudaFree(_sorted_errors));
      (cudaFree(_part_errors));

      // store copy of coefficients for future calculation
      (cudaMemcpy(_pnm_re0[dev], _pnm_re[dev],
        coeff_stride*nparts*sizeof(real), cudaMemcpyDeviceToDevice));
      (cudaMemcpy(_pnm_im0[dev], _pnm_im[dev],
        coeff_stride*nparts*sizeof(real), cudaMemcpyDeviceToDevice));
      (cudaMemcpy(_phinm_re0[dev], _phinm_re[dev],
        coeff_stride*nparts*sizeof(real), cudaMemcpyDeviceToDevice));
      (cudaMemcpy(_phinm_im0[dev], _phinm_im[dev],
        coeff_stride*nparts*sizeof(real), cudaMemcpyDeviceToDevice));
      (cudaMemcpy(_chinm_re0[dev], _chinm_re[dev],
        coeff_stride*nparts*sizeof(real), cudaMemcpyDeviceToDevice));
      (cudaMemcpy(_chinm_im0[dev], _chinm_im[dev],
        coeff_stride*nparts*sizeof(real), cudaMemcpyDeviceToDevice));

      // clean up
      free(part_errors);
    }

    // find maximum error of all subdomains
    for(int i = 0; i < nsubdom; i++) {
      if(errors[i] > max) max = errors[i];
    }

    // clean up
    free(errors);

    return max;
  } else return 0.;
}
