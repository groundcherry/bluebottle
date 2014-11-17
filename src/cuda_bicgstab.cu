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
#include "cuda_bluebottle.h"
#include "entrySearch.h"
#include "cuda_particle.h"

#ifdef TEST
#include "cuda_testing.h"
#endif

#include <cuda.h>
#include <cusp/array1d.h>
#include <cusp/blas.h>
#include <cusp/dia_matrix.h>
#include <cusp/monitor.h>
#include <cusp/precond/diagonal.h>
#include <cusp/krylov/bicgstab.h>
//#include <cusp/krylov/cg.h>
#include <cusp/print.h>

extern "C"
void cuda_PP_bicgstab(int rank)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    // write right-hand side
    cuda_PP_rhs(dev);
    if(nparts > 0) {
      cuda_part_BC_p(dev);
    }

    // create temporary U* without ghost cells for bicgstab result
    cusp::array1d<real, cusp::device_memory> p_tmp(dom[dev].Gcc.s3, 0.);

    // create CUSP diagonal matrix for p
    cusp::dia_matrix<int, real, cusp::device_memory> *_A_p;
    _A_p = new cusp::dia_matrix<int, real, cusp::device_memory>
      (dom[dev].Gcc._s3, dom[dev].Gcc._s3, 0, 13);

    // set up the coefficient matrix
    _A_p->diagonal_offsets[0]  = -dom[dev].Gcc._s3 + dom[dev].Gcc._s2;
    _A_p->diagonal_offsets[1]  = -dom[dev].Gcc._s2;
    _A_p->diagonal_offsets[2]  = -dom[dev].Gcc._s2 + dom[dev].Gcc._s1;
    _A_p->diagonal_offsets[3]  = -dom[dev].Gcc._s1;
    _A_p->diagonal_offsets[4]  = -dom[dev].Gcc._s1 + 1;
    _A_p->diagonal_offsets[5]  = -1;
    _A_p->diagonal_offsets[6]  = 0;
    _A_p->diagonal_offsets[7]  = 1;
    _A_p->diagonal_offsets[8]  = dom[dev].Gcc._s1 - 1;
    _A_p->diagonal_offsets[9]  = dom[dev].Gcc._s1;
    _A_p->diagonal_offsets[10] = dom[dev].Gcc._s2 - dom[dev].Gcc._s1;
    _A_p->diagonal_offsets[11] = dom[dev].Gcc._s2;
    _A_p->diagonal_offsets[12] = dom[dev].Gcc._s3 - dom[dev].Gcc._s2;

    // write coefficients using kernel
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gcc._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gcc._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gcc._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gcc._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gcc._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gcc._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gcc._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gcc._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gcc._kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);
    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);
    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // build pressure-Poisson coefficient matrix
    coeffs_init<<<numBlocks_x, dimBlocks_x>>>(_dom[dev], _A_p->values.pitch,
      thrust::raw_pointer_cast(&_A_p->values.values[0]));
    coeffs<<<numBlocks_x, dimBlocks_x>>>(_dom[dev], _flag_u[dev], _flag_v[dev],
      _flag_w[dev], _A_p->values.pitch,
      thrust::raw_pointer_cast(&_A_p->values.values[0]));
    if(bc.pW == PERIODIC)
      coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]));
    if(bc.pE == PERIODIC)
      coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]));
    if(bc.pS == PERIODIC)
      coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]));
    if(bc.pN == PERIODIC)
      coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]));
    if(bc.pB == PERIODIC)
      coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]));
    if(bc.pT == PERIODIC)
      coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]));

    coeffs_particle<<<numBlocks_x, dimBlocks_x>>>(_dom[dev], _A_p->values.pitch,
      thrust::raw_pointer_cast(&_A_p->values.values[0]), _phase_shell[dev]);

    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_p(_rhs_p[dev]);
    cusp::array1d<real, cusp::device_memory> *_pp;
    _pp = new cusp::array1d<real, cusp::device_memory>(_ptr_p,
      _ptr_p + dom[dev].Gcc._s3);

    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_pp);
    if(norm == 0)
      norm = 1.;
    cusp::blas::scal(*_pp, 1. / norm);

/*
cusp::array1d<real, cusp::host_memory> PP = *_pp;
cusp::blas::scal(PP, norm);
//cusp::print(PP);
real ppsum = 0.;
for(int s = 0; s < dom[dev].Gcc.s3; s++) {
  ppsum += PP[s]*dom[dev].dx*dom[dev].dy*dom[dev].dz;
}

printf("PPSUM_1 = %e\n", ppsum*norm*dt/rho_f);
*/

    // call BiCGSTAB to solve for p_tmp
    cusp::convergence_monitor<real> monitor(*_pp, pp_max_iter, pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_p);
    cusp::krylov::bicgstab(*_A_p, p_tmp, *_pp, monitor, M);
    //cusp::krylov::cg(*_A_p, p_tmp, *_pp, monitor, M);
    // write convergence data to file
    if(rank == 0)
      recorder_bicgstab("solver_flow.rec", monitor.iteration_count(),
        monitor.residual_norm());
    else
      recorder_bicgstab("solver_turb.rec", monitor.iteration_count(),
        monitor.residual_norm());
    if(!monitor.converged()) {
      printf("The pressure-Poisson equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(p_tmp, norm);

    // calculate average pressure
    real p_avg = avg_entries(dom[dev].Gcc.s3,
      thrust::raw_pointer_cast(p_tmp.data()));

    // subtract average value from pressure
    cusp::array1d<real, cusp::device_memory> ones(dom[dev].Gcc.s3, 1.);
    cusp::blas::axpy(ones, p_tmp, -p_avg);

    // copy solution back to pressure field
    copy_p_ghost<<<numBlocks_x, dimBlocks_x>>>(_p[dev],
      thrust::raw_pointer_cast(p_tmp.data()), _dom[dev]);

    // clean up
    delete(_pp);
    delete(_A_p);
  }
}

/*extern "C"
void cuda_div_U(void)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    // write right-hand side
    cuda_div_U_launch(dev, _u[dev], _v[dev], _w[dev], _divU[dev]);
  }
}
*/

extern "C"
void cuda_PP_rhs(int dev)
{
  int threads_x = 0;
  int threads_y = 0;
  int blocks_x = 0;
  int blocks_y = 0;

  if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
    threads_x = dom[dev].Gcc.inb + 2;
  else
    threads_x = MAX_THREADS_DIM;

  if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
    threads_y = dom[dev].Gcc.jnb + 2;
  else
    threads_y = MAX_THREADS_DIM;

  blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) (threads_x-2));
  blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) (threads_y-2));

  dim3 dimBlocks(threads_x, threads_y);
  dim3 numBlocks(blocks_x, blocks_y);

  PP_rhs<<<numBlocks, dimBlocks>>>(rho_f, _u_star[dev], _v_star[dev],
    _w_star[dev], _rhs_p[dev], _dom[dev], dt);
}

extern "C"
void cuda_div_U_launch(int dev, real *u, real *v, real *w, real *divU)
{
  int threads_x = 0;
  int threads_y = 0;
  int blocks_x = 0;
  int blocks_y = 0;

  if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
    threads_x = dom[dev].Gcc.inb + 2;
  else
    threads_x = MAX_THREADS_DIM;

  if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
    threads_y = dom[dev].Gcc.jnb + 2;
  else
    threads_y = MAX_THREADS_DIM;

  blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) (threads_x-2));
  blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) (threads_y-2));

  dim3 dimBlocks(threads_x, threads_y);
  dim3 numBlocks(blocks_x, blocks_y);

  div_U<<<numBlocks, dimBlocks>>>(u, v, w, divU, _dom[dev]);
}
