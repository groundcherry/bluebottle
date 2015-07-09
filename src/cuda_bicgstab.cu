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
#include <cusp/array1d.h>
#include <cusp/blas/blas.h>
#include <cusp/dia_matrix.h>
#include <cusp/monitor.h>
#include <cusp/precond/diagonal.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/cg.h>
#include <cusp/print.h>
#include <thrust/device_ptr.h>

#include "cuda_bicgstab.h"
#include "cuda_bluebottle.h"
#include "entrySearch.h"
#include "cuda_particle.h"

#ifdef TEST
#include "cuda_testing.h"
#endif

extern "C"
void cuda_ustar_helmholtz(int rank)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    // write right-hand side
    cuda_ustar_rhs(dev);
    cuda_dom_BC_star();
    if(nparts > 0) {
      cuda_part_BC_star();
    }

    // create temporary U* without ghost cells for bicgstab result
    cusp::array1d<real, cusp::device_memory> ustar_tmp(dom[dev].Gfx.s3, 0.);

    // create CUSP diagonal matrix for p
    cusp::dia_matrix<int, real, cusp::device_memory> *_A_ustar;
    _A_ustar = new cusp::dia_matrix<int, real, cusp::device_memory>
      (dom[dev].Gfx._s3, dom[dev].Gfx._s3, 0, 13);

    // set up the coefficient matrix
    _A_ustar->diagonal_offsets[0]  = -dom[dev].Gfx._s3 + dom[dev].Gfx._s2;
    _A_ustar->diagonal_offsets[1]  = -dom[dev].Gfx._s2;
    _A_ustar->diagonal_offsets[2]  = -dom[dev].Gfx._s2 + dom[dev].Gfx._s1;
    _A_ustar->diagonal_offsets[3]  = -dom[dev].Gfx._s1;
    _A_ustar->diagonal_offsets[4]  = -dom[dev].Gfx._s1 + 2;
    _A_ustar->diagonal_offsets[5]  = -1;
    _A_ustar->diagonal_offsets[6]  = 0;
    _A_ustar->diagonal_offsets[7]  = 1;
    _A_ustar->diagonal_offsets[8]  = dom[dev].Gfx._s1 - 2;
    _A_ustar->diagonal_offsets[9]  = dom[dev].Gfx._s1;
    _A_ustar->diagonal_offsets[10] = dom[dev].Gfx._s2 - dom[dev].Gfx._s1;
    _A_ustar->diagonal_offsets[11] = dom[dev].Gfx._s2;
    _A_ustar->diagonal_offsets[12] = dom[dev].Gfx._s3 - dom[dev].Gfx._s2;

    // write coefficients using kernel
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gfx._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfx._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfx._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfx._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfx._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx._kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);
    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);
    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // create temporary ustar without ghost cells
    real *_ustar_noghost;
    (cudaMalloc((void**) &_ustar_noghost,
      sizeof(real) * dom[dev].Gfx.s3));
    // copy u_star into noghost structure for Helmholtz right-hand side
    copy_u_noghost<<<numBlocks_x, dimBlocks_x>>>(_ustar_noghost, _u_star[dev],
      _dom[dev]);

    // build pressure-Poisson coefficient matrix
    ustar_coeffs_init<<<numBlocks_x, dimBlocks_x>>>(_dom[dev],
      _A_ustar->values.pitch,
      thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    ustar_coeffs<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
      _A_ustar->values.pitch,
      thrust::raw_pointer_cast(&_A_ustar->values.values[0]),
      _flag_u[dev], _flag_v[dev], _flag_w[dev]);

/*
cusp::dia_matrix<int, real, cusp::host_memory> AA = *_A_ustar;
printf("\n");
for(int i = 0; i < dom[dev].Gfx.s3; i++) {
  for(int j = 0; j < dom[dev].Gfx.s3; j++) {
    if(j == AA.diagonal_offsets[0] + i)
      if(AA.values(i, 0) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 0));
    else if(j == AA.diagonal_offsets[1] + i)
      if(AA.values(i, 1) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 1));
    else if(j == AA.diagonal_offsets[2] + i)
      if(AA.values(i, 2) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 2));
    else if(j == AA.diagonal_offsets[3] + i)
      if(AA.values(i, 3) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 3));
    else if(j == AA.diagonal_offsets[4] + i)
      if(AA.values(i, 4) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 4));
    else if(j == AA.diagonal_offsets[5] + i)
      if(AA.values(i, 5) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 5));
    else if(j == AA.diagonal_offsets[6] + i)
      if(AA.values(i, 6) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 6));
    else if(j == AA.diagonal_offsets[7] + i)
      if(AA.values(i, 7) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 7));
    else if(j == AA.diagonal_offsets[8] + i)
      if(AA.values(i, 8) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 8));
    else if(j == AA.diagonal_offsets[9] + i)
      if(AA.values(i, 9) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 9));
    else if(j == AA.diagonal_offsets[10] + i)
      if(AA.values(i, 10) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 10));
    else if(j == AA.diagonal_offsets[11] + i)
      if(AA.values(i, 11) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 11));
    else if(j == AA.diagonal_offsets[12] + i)
      if(AA.values(i, 12) == 0)
        printf("    ");
      else
        printf("%3.0f ", AA.values(i, 12));
    else
      printf("    ");
  }
  printf("\n");
}
*/

    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_ustar(_ustar_noghost);
    cusp::array1d<real, cusp::device_memory> *_ustar_rhs;
    _ustar_rhs = new cusp::array1d<real, cusp::device_memory>(_ptr_ustar,
      _ptr_ustar + dom[dev].Gfx._s3);

    // account for boundary conditions
    if(bc.uW == PERIODIC)
      ustar_coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
        _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    else if(bc.uW == DIRICHLET)
      ustar_coeffs_dirichlet_W<<<numBlocks_x, dimBlocks_x>>>(bc.uWD,
        _u_star[dev], _dom[dev], _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    if(bc.uE == PERIODIC)
      ustar_coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
        _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    else if(bc.uE == DIRICHLET)
      ustar_coeffs_dirichlet_E<<<numBlocks_x, dimBlocks_x>>>(bc.uED,
        _u_star[dev], _dom[dev], _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    if(bc.uS == PERIODIC)
      ustar_coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
        _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    else if(bc.uS == DIRICHLET)
      ustar_coeffs_dirichlet_S<<<numBlocks_y, dimBlocks_y>>>(bc.uSD,
        _u_star[dev], _dom[dev], _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    if(bc.uN == PERIODIC)
      ustar_coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
        _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    else if(bc.uN == DIRICHLET)
      ustar_coeffs_dirichlet_N<<<numBlocks_y, dimBlocks_y>>>(bc.uND,
        _u_star[dev], _dom[dev], _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    if(bc.uB == PERIODIC)
      ustar_coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
        _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    else if(bc.uB == DIRICHLET)
      ustar_coeffs_dirichlet_B<<<numBlocks_z, dimBlocks_z>>>(bc.uBD,
        _u_star[dev], _dom[dev], _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    if(bc.uT == PERIODIC)
      ustar_coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
        _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));
    else if(bc.uT == DIRICHLET)
      ustar_coeffs_dirichlet_T<<<numBlocks_z, dimBlocks_z>>>(bc.uTD,
        _u_star[dev], _dom[dev], _A_ustar->values.pitch,
        thrust::raw_pointer_cast(&_A_ustar->values.values[0]));

    ustar_coeffs_particles<<<numBlocks_x, dimBlocks_x>>>(_dom[dev],
      _A_ustar->values.pitch,
      thrust::raw_pointer_cast(&_A_ustar->values.values[0]),
      _flag_u[dev]);

    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_ustar_rhs);
    if(norm == 0)
      norm = 1.;
    cusp::blas::scal(*_ustar_rhs, 1. / norm);

    // call BiCGSTAB to solve for ustar_tmp
    cusp::monitor<real> monitor(*_ustar_rhs, pp_max_iter,
      pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_ustar);
    //cusp::krylov::bicgstab(*_A_ustar, ustar_tmp, *_ustar_rhs, monitor, M);
    cusp::krylov::cg(*_A_ustar, ustar_tmp, *_ustar_rhs, monitor, M);
    // write convergence data to file
    if(rank == 0) {
      char nam[FILE_NAME_SIZE] = "solver_helmholtz_expd.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    } else { 
      char nam[FILE_NAME_SIZE] = "solver_helmholtz_prec.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    }
    if(!monitor.converged()) {
      printf("The u_star Helmholtz equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(ustar_tmp, norm);

    // copy solution back to _u_star
    copy_u_ghost<<<numBlocks_x, dimBlocks_x>>>(_u_star[dev],
      thrust::raw_pointer_cast(ustar_tmp.data()), _dom[dev]);

    // clean up
    delete(_ustar_rhs);
    delete(_A_ustar);
    (cudaFree(_ustar_noghost));

#ifdef TEST
  // copy _u_star to _u
  cudaMemcpy(_u[dev], _u_star[dev], dom[dev].Gfx.s3b * sizeof(real),
    cudaMemcpyDeviceToDevice);
#endif
  }
}

extern "C"
void cuda_vstar_helmholtz(int rank)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    // write right-hand side
    cuda_vstar_rhs(dev);
    cuda_dom_BC_star();
    if(nparts > 0) {
      cuda_part_BC_star();
    }

    // create temporary U* without ghost cells for bicgstab result
    cusp::array1d<real, cusp::device_memory> vstar_tmp(dom[dev].Gfy.s3, 0.);

    // create CUSP diagonal matrix for p
    cusp::dia_matrix<int, real, cusp::device_memory> *_A_vstar;
    _A_vstar = new cusp::dia_matrix<int, real, cusp::device_memory>
      (dom[dev].Gfy._s3, dom[dev].Gfy._s3, 0, 13);

    // set up the coefficient matrix
    _A_vstar->diagonal_offsets[0]  = -dom[dev].Gfy._s3 + dom[dev].Gfy._s2;
    _A_vstar->diagonal_offsets[1]  = -dom[dev].Gfy._s2;
    _A_vstar->diagonal_offsets[2]  = -dom[dev].Gfy._s2 + 2*dom[dev].Gfy._s1;
    _A_vstar->diagonal_offsets[3]  = -dom[dev].Gfy._s1;
    _A_vstar->diagonal_offsets[4]  = -dom[dev].Gfy._s1 + 1;
    _A_vstar->diagonal_offsets[5]  = -1;
    _A_vstar->diagonal_offsets[6]  = 0;
    _A_vstar->diagonal_offsets[7]  = 1;
    _A_vstar->diagonal_offsets[8]  = dom[dev].Gfy._s1 - 1;
    _A_vstar->diagonal_offsets[9]  = dom[dev].Gfy._s1;
    _A_vstar->diagonal_offsets[10] = dom[dev].Gfy._s2 - 2*dom[dev].Gfy._s1;
    _A_vstar->diagonal_offsets[11] = dom[dev].Gfy._s2;
    _A_vstar->diagonal_offsets[12] = dom[dev].Gfy._s3 - dom[dev].Gfy._s2;

    // write coefficients using kernel
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gfy._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfy._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfy._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfy._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfy._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfy._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfy._kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);
    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);
    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // create temporary ustar without ghost cells
    real *_vstar_noghost;
    (cudaMalloc((void**) &_vstar_noghost,
      sizeof(real) * dom[dev].Gfy.s3));
    // copy v_star into noghost structure for Helmholtz right-hand side
    copy_v_noghost<<<numBlocks_y, dimBlocks_y>>>(_vstar_noghost, _v_star[dev],
      _dom[dev]);

    // build pressure-Poisson coefficient matrix
    vstar_coeffs_init<<<numBlocks_y, dimBlocks_y>>>(_dom[dev],
      _A_vstar->values.pitch,
      thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    vstar_coeffs<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
      _A_vstar->values.pitch,
      thrust::raw_pointer_cast(&_A_vstar->values.values[0]),
      _flag_u[dev], _flag_v[dev], _flag_w[dev]);

    // account for boundary conditions
    if(bc.vW == PERIODIC)
      vstar_coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
        _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    else if(bc.vW == DIRICHLET)
      vstar_coeffs_dirichlet_W<<<numBlocks_x, dimBlocks_x>>>(bc.vWD,
        _v_star[dev], _dom[dev], _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    if(bc.vE == PERIODIC)
      vstar_coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
        _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    else if(bc.vE == DIRICHLET)
      vstar_coeffs_dirichlet_E<<<numBlocks_x, dimBlocks_x>>>(bc.vED,
        _v_star[dev], _dom[dev], _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    if(bc.vS == PERIODIC)
      vstar_coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
        _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    else if(bc.vS == DIRICHLET)
      vstar_coeffs_dirichlet_S<<<numBlocks_y, dimBlocks_y>>>(bc.vSD,
        _v_star[dev], _dom[dev], _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    if(bc.vN == PERIODIC)
      vstar_coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
        _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    else if(bc.vN == DIRICHLET)
      vstar_coeffs_dirichlet_N<<<numBlocks_y, dimBlocks_y>>>(bc.vND,
        _v_star[dev], _dom[dev], _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    if(bc.vB == PERIODIC)
      vstar_coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
        _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    else if(bc.vB == DIRICHLET)
      vstar_coeffs_dirichlet_B<<<numBlocks_z, dimBlocks_z>>>(bc.vBD,
        _v_star[dev], _dom[dev], _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    if(bc.vT == PERIODIC)
      vstar_coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
        _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));
    else if(bc.vT == DIRICHLET)
      vstar_coeffs_dirichlet_T<<<numBlocks_z, dimBlocks_z>>>(bc.vTD,
        _v_star[dev], _dom[dev], _A_vstar->values.pitch,
        thrust::raw_pointer_cast(&_A_vstar->values.values[0]));

    vstar_coeffs_particles<<<numBlocks_y, dimBlocks_y>>>(_dom[dev],
      _A_vstar->values.pitch,
      thrust::raw_pointer_cast(&_A_vstar->values.values[0]),
      _flag_v[dev]);

    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_vstar(_vstar_noghost);
    cusp::array1d<real, cusp::device_memory> *_vstar_rhs;
    _vstar_rhs = new cusp::array1d<real, cusp::device_memory>(_ptr_vstar,
      _ptr_vstar + dom[dev].Gfy._s3);

    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_vstar_rhs);
    if(norm == 0)
      norm = 1.;
    cusp::blas::scal(*_vstar_rhs, 1. / norm);

    // call BiCGSTAB to solve for ustar_tmp
    cusp::monitor<real> monitor(*_vstar_rhs, pp_max_iter,
      pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_vstar);
    //cusp::krylov::bicgstab(*_A_vstar, vstar_tmp, *_vstar_rhs, monitor, M);
    cusp::krylov::cg(*_A_vstar, vstar_tmp, *_vstar_rhs, monitor, M);
    // write convergence data to file
    if(rank == 0) {
      char nam[FILE_NAME_SIZE] = "solver_helmholtz_expd.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    } else {
      char nam[FILE_NAME_SIZE] = "solver_helmholtz_prec.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    }
    if(!monitor.converged()) {
      printf("The v_star Helmholtz equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(vstar_tmp, norm);

    // copy solution back to _v_star field
    copy_v_ghost<<<numBlocks_y, dimBlocks_y>>>(_v_star[dev],
      thrust::raw_pointer_cast(vstar_tmp.data()), _dom[dev]);

    // clean up
    delete(_vstar_rhs);
    delete(_A_vstar);
    (cudaFree(_vstar_noghost));

#ifdef TEST
  // copy _v_star to _v
  cudaMemcpy(_v[dev], _v_star[dev], dom[dev].Gfy.s3b * sizeof(real),
    cudaMemcpyDeviceToDevice);
#endif
  }
}

extern "C"
void cuda_wstar_helmholtz(int rank)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);

    // write right-hand side
    cuda_wstar_rhs(dev);
    cuda_dom_BC_star();
    if(nparts > 0) {
      cuda_part_BC_star();
    }

    // create temporary U* without ghost cells for bicgstab result
    cusp::array1d<real, cusp::device_memory> wstar_tmp(dom[dev].Gfz.s3, 0.);

    // create CUSP diagonal matrix for p
    cusp::dia_matrix<int, real, cusp::device_memory> *_A_wstar;
    _A_wstar = new cusp::dia_matrix<int, real, cusp::device_memory>
      (dom[dev].Gfz._s3, dom[dev].Gfz._s3, 0, 13);

    // set up the coefficient matrix
    _A_wstar->diagonal_offsets[0]  = -dom[dev].Gfz._s3 + 2*dom[dev].Gfz._s2;
    _A_wstar->diagonal_offsets[1]  = -dom[dev].Gfz._s2;
    _A_wstar->diagonal_offsets[2]  = -dom[dev].Gfz._s2 + dom[dev].Gfz._s1;
    _A_wstar->diagonal_offsets[3]  = -dom[dev].Gfz._s1;
    _A_wstar->diagonal_offsets[4]  = -dom[dev].Gfz._s1 + 1;
    _A_wstar->diagonal_offsets[5]  = -1;
    _A_wstar->diagonal_offsets[6]  = 0;
    _A_wstar->diagonal_offsets[7]  = 1;
    _A_wstar->diagonal_offsets[8]  = dom[dev].Gfz._s1 - 1;
    _A_wstar->diagonal_offsets[9]  = dom[dev].Gfz._s1;
    _A_wstar->diagonal_offsets[10] = dom[dev].Gfz._s2 - dom[dev].Gfz._s1;
    _A_wstar->diagonal_offsets[11] = dom[dev].Gfz._s2;
    _A_wstar->diagonal_offsets[12] = dom[dev].Gfz._s3 - 2*dom[dev].Gfz._s2;

    // write coefficients using kernel
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gfz._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfz._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfz._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfz._kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);
    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);
    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // create temporary ustar without ghost cells
    real *_wstar_noghost;
    (cudaMalloc((void**) &_wstar_noghost,
      sizeof(real) * dom[dev].Gfz.s3));
    // copy w_star into noghost structure for Helmholtz right-hand side
    copy_w_noghost<<<numBlocks_z, dimBlocks_z>>>(_wstar_noghost, _w_star[dev],
      _dom[dev]);

    // build pressure-Poisson coefficient matrix
    wstar_coeffs_init<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
      _A_wstar->values.pitch,
      thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    wstar_coeffs<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
      _A_wstar->values.pitch,
      thrust::raw_pointer_cast(&_A_wstar->values.values[0]),
      _flag_u[dev], _flag_v[dev], _flag_w[dev]);

    // account for boundary conditions
    if(bc.wW == PERIODIC)
      wstar_coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
        _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    else if(bc.wW == DIRICHLET)
      wstar_coeffs_dirichlet_W<<<numBlocks_x, dimBlocks_x>>>(bc.wWD,
        _w_star[dev], _dom[dev], _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    if(bc.wE == PERIODIC)
      wstar_coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(nu, dt, _dom[dev],
        _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    else if(bc.wE == DIRICHLET)
      wstar_coeffs_dirichlet_E<<<numBlocks_x, dimBlocks_x>>>(bc.wED,
        _w_star[dev], _dom[dev], _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    if(bc.wS == PERIODIC)
      wstar_coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
        _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    else if(bc.wS == DIRICHLET)
      wstar_coeffs_dirichlet_S<<<numBlocks_y, dimBlocks_y>>>(bc.wSD,
        _w_star[dev], _dom[dev], _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    if(bc.wN == PERIODIC)
      wstar_coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(nu, dt, _dom[dev],
        _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    else if(bc.wN == DIRICHLET)
      wstar_coeffs_dirichlet_N<<<numBlocks_y, dimBlocks_y>>>(bc.wND,
        _w_star[dev], _dom[dev], _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    if(bc.wB == PERIODIC)
      wstar_coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
        _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    else if(bc.wB == DIRICHLET)
      wstar_coeffs_dirichlet_B<<<numBlocks_z, dimBlocks_z>>>(bc.wBD,
        _w_star[dev], _dom[dev], _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    if(bc.wT == PERIODIC)
      wstar_coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(nu, dt, _dom[dev],
        _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));
    else if(bc.wT == DIRICHLET)
      wstar_coeffs_dirichlet_T<<<numBlocks_z, dimBlocks_z>>>(bc.wTD,
        _w_star[dev], _dom[dev], _A_wstar->values.pitch,
        thrust::raw_pointer_cast(&_A_wstar->values.values[0]));

    wstar_coeffs_particles<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
      _A_wstar->values.pitch,
      thrust::raw_pointer_cast(&_A_wstar->values.values[0]),
      _flag_w[dev]);

    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_wstar(_wstar_noghost);
    cusp::array1d<real, cusp::device_memory> *_wstar_rhs;
    _wstar_rhs = new cusp::array1d<real, cusp::device_memory>(_ptr_wstar,
      _ptr_wstar + dom[dev].Gfz.s3);

    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_wstar_rhs);
    if(norm == 0)
      norm = 1.;
    cusp::blas::scal(*_wstar_rhs, 1. / norm);

    // call BiCGSTAB to solve for ustar_tmp
    cusp::monitor<real> monitor(*_wstar_rhs, pp_max_iter,
      pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_wstar);
    //cusp::krylov::bicgstab(*_A_wstar, wstar_tmp, *_wstar_rhs, monitor, M);
    cusp::krylov::cg(*_A_wstar, wstar_tmp, *_wstar_rhs, monitor, M);
    // write convergence data to file
    if(rank == 0) {
      char nam[FILE_NAME_SIZE] = "solver_helmholtz_expd.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    } else {
      char nam[FILE_NAME_SIZE] = "solver_helmholtz_prec.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    }
    if(!monitor.converged()) {
      printf("The w_star Helmholtz equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(wstar_tmp, norm);

    // copy solution back to _v_star field
    copy_w_ghost<<<numBlocks_z, dimBlocks_z>>>(_w_star[dev],
      thrust::raw_pointer_cast(wstar_tmp.data()), _dom[dev]);

    // clean up
    delete(_wstar_rhs);
    delete(_A_wstar);
    (cudaFree(_wstar_noghost));

#ifdef TEST
  // copy _w_star to _w
  cudaMemcpy(_w[dev], _w_star[dev], dom[dev].Gfz.s3b * sizeof(real),
    cudaMemcpyDeviceToDevice);
#endif
  }
}

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
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]),
        _flag_u[dev]);
    if(bc.pE == PERIODIC)
      coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]),
        _flag_u[dev]);
    if(bc.pS == PERIODIC)
      coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]),
        _flag_v[dev]);
    if(bc.pN == PERIODIC)
      coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]),
        _flag_v[dev]);
    if(bc.pB == PERIODIC)
      coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]),
        _flag_w[dev]);
    if(bc.pT == PERIODIC)
      coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(_dom[dev],
        _A_p->values.pitch, thrust::raw_pointer_cast(&_A_p->values.values[0]),
        _flag_w[dev]);

    coeffs_particle<<<numBlocks_x, dimBlocks_x>>>(_dom[dev], _A_p->values.pitch,
      thrust::raw_pointer_cast(&_A_p->values.values[0]), _phase[dev]);

/*    cusp::dia_matrix<int, real, cusp::host_memory> AA = *_A_p;
    FILE *mat = fopen("mat.txt", "w");
    for(int i = 0; i < dom[dev].Gcc.s3; i++) {
      for(int j = 0; j < dom[dev].Gcc.s3; j++) {
        if(j == AA.diagonal_offsets[0] + i)
          if(AA.values(i, 0) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 0));
        else if(j == AA.diagonal_offsets[1] + i)
          if(AA.values(i, 1) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 1));
        else if(j == AA.diagonal_offsets[2] + i)
          if(AA.values(i, 2) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 2));
        else if(j == AA.diagonal_offsets[3] + i)
          if(AA.values(i, 3) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 3));
        else if(j == AA.diagonal_offsets[4] + i)
          if(AA.values(i, 4) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 4));
        else if(j == AA.diagonal_offsets[5] + i)
          if(AA.values(i, 5) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 5));
        else if(j == AA.diagonal_offsets[6] + i)
          if(AA.values(i, 6) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 6));
        else if(j == AA.diagonal_offsets[7] + i)
          if(AA.values(i, 7) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 7));
        else if(j == AA.diagonal_offsets[8] + i)
          if(AA.values(i, 8) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 8));
        else if(j == AA.diagonal_offsets[9] + i)
          if(AA.values(i, 9) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 9));
        else if(j == AA.diagonal_offsets[10] + i)
          if(AA.values(i, 10) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 10));
        else if(j == AA.diagonal_offsets[11] + i)
          if(AA.values(i, 11) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 11));
        else if(j == AA.diagonal_offsets[12] + i)
          if(AA.values(i, 12) == 0)
            fprintf(mat,"%f ", 0.);
          else
            fprintf(mat,"%f ", AA.values(i, 12));
        else
            fprintf(mat,"%f ", 0.);
      }
    }
    fclose(mat);
*/

    // copy p0 to array without ghost cells and use it as an initial guess and solution space
    real *_phinoghost;
    (cudaMalloc((void**) &_phinoghost,
      sizeof(real)*dom[dev].Gcc.s3b));
    copy_p_noghost<<<numBlocks_x, dimBlocks_x>>>(_phinoghost, _phi[dev],
      _dom[dev]);
    thrust::device_ptr<real> _ptr_p_sol(_phinoghost);
    cusp::array1d<real, cusp::device_memory> *_p_sol;
    _p_sol = new cusp::array1d<real, cusp::device_memory>(_ptr_p_sol,
      _ptr_p_sol + dom[dev].Gcc._s3);

    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_p(_rhs_p[dev]);
    cusp::array1d<real, cusp::device_memory> *_pp;
    _pp = new cusp::array1d<real, cusp::device_memory>(_ptr_p,
      _ptr_p + dom[dev].Gcc._s3);

/*
printf("%e\n", dt);
printf("_p_sol_in\n");
cusp::print(*_p_sol);
printf("_pp\n");
cusp::print(*_pp);
*/

    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_pp);
    if(norm == 0)
      norm = 1.;
    cusp::blas::scal(*_pp, 1. / norm);
    cusp::blas::scal(*_p_sol, 1. / norm);

    // call BiCGSTAB to solve for p_sol
    cusp::monitor<real> monitor(*_pp, pp_max_iter, pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_p);
    cusp::krylov::bicgstab(*_A_p, *_p_sol, *_pp, monitor, M);
    // write convergence data to file
    if(rank == 0) {
      char nam[FILE_NAME_SIZE] = "solver_expd.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    } else {
      char nam[FILE_NAME_SIZE] = "solver_prec.rec";
      recorder_bicgstab(nam, monitor.iteration_count(),
        monitor.residual_norm());
    }
    if(!monitor.converged()) {
      printf("The pressure-Poisson equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(*_p_sol, norm);

    // calculate average pressure
    real p_avg = avg_entries(dom[dev].Gcc.s3,
      thrust::raw_pointer_cast(_p_sol->data()));

    // subtract average value from pressure
    cusp::array1d<real, cusp::device_memory> ones(dom[dev].Gcc.s3, 1.);
    cusp::blas::axpy(ones, *_p_sol, -p_avg);

/*
printf("_p_sol_out\n");
cusp::print(*_p_sol);
*/

    // copy solution back to pressure field
    copy_p_ghost<<<numBlocks_x, dimBlocks_x>>>(_phi[dev],
      thrust::raw_pointer_cast(_p_sol->data()), _dom[dev]);

    // clean up
    (cudaFree(_phinoghost));
    delete(_p_sol);
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
void cuda_ustar_rhs(int dev)
{
  int threads_y = 0;
  int threads_z = 0;
  int blocks_y = 0;
  int blocks_z = 0;

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

  dim3 dimBlocks(threads_y, threads_z);
  dim3 numBlocks(blocks_y, blocks_z);

  ustar_rhs<<<numBlocks, dimBlocks>>>(rho_f, nu, _u0[dev], _v0[dev], _w0[dev],
    _p0[dev], _f_x[dev], _conv0_u[dev], _conv_u[dev], _u_star[dev], _dom[dev],
    dt, dt0);
}

extern "C"
void cuda_vstar_rhs(int dev)
{
  int threads_z = 0;
  int threads_x = 0;
  int blocks_z = 0;
  int blocks_x = 0;

  if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
    threads_z = dom[dev].Gfx.knb + 2;
  else
    threads_z = MAX_THREADS_DIM;

  if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
    threads_x = dom[dev].Gfx.inb + 2;
  else
    threads_x = MAX_THREADS_DIM;

  blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) (threads_z-2));
  blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) (threads_x-2));

  dim3 dimBlocks(threads_z, threads_x);
  dim3 numBlocks(blocks_z, blocks_x);

  vstar_rhs<<<numBlocks, dimBlocks>>>(rho_f, nu, _u0[dev], _v0[dev], _w0[dev],
    _p0[dev], _f_y[dev], _conv0_v[dev], _conv_v[dev], _v_star[dev], _dom[dev],
    dt, dt0);
}

extern "C"
void cuda_wstar_rhs(int dev)
{
  int threads_x = 0;
  int threads_y = 0;
  int blocks_x = 0;
  int blocks_y = 0;

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

  dim3 dimBlocks(threads_x, threads_y);
  dim3 numBlocks(blocks_x, blocks_y);

  wstar_rhs<<<numBlocks, dimBlocks>>>(rho_f, nu, _u0[dev], _v0[dev], _w0[dev],
    _p0[dev], _f_z[dev], _conv0_w[dev], _conv_w[dev], _w_star[dev], _dom[dev],
    dt, dt0);
}

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
