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

/****h* Bluebottle/cuda_quadrature_kernel
 * NAME
 *  cuda_quadrature_kernel
 * FUNCTION
 *  Lebedev quadrature CUDA kernel functions.
 ******
 */

#ifndef _CUDA_QUADRATURE_H
#define _CUDA_QUADRATURE_H

extern "C"
{
#include "bluebottle.h"
#include "particle.h"
}

/****f* cuda_quadrature_kernel/rtp2xyz<<<>>>()
 * NAME
 *  rtp2xyz<<<>>>()
 * USAGE
 */
__device__ void rtp2xyz(real r, real t, real p, real *x, real *y, real *z);
/*
 * PURPOSE
 *  Compute (x, y, z) from (r, theta, phi).
 * ARGUMENTS
 *  r -- r-component in spherical basis
 *  t -- theta-component in spherical basis
 *  p -- phi-component in spherical basis
 *  x -- x-component in Cartesian basis
 *  y -- y-component in Cartesian basis
 *  z -- z-component in Cartesian basis
 ******
 */

/****f* cuda_quadrature_kernel/cart2sphere<<<>>>()
 * NAME
 *  cart2sphere<<<>>>()
 * USAGE
 */
__device__ void cart2sphere(real u, real v, real w, real t, real p,
  real *ur, real *ut, real *up);
/*
 * PURPOSE
 *  Compute (u_r, u_theta, u_phi) from (u, v, w).
 * ARGUMENTS
 *  u -- the x-component of velocity in Cartesian basis
 *  v -- the y-component of velocity in Cartesian basis
 *  w -- the z-component of velocity in Cartesian basis
 *  t -- theta-component in spherical basis
 *  p -- phi-component in spherical basis
 *  ur -- r-component of velocity in spherical basis
 *  ut -- theta-component of velocity in spherical basis
 *  up -- phi-component of velocity in spherical basis
 ******
 */

/****f* cuda_quadrature_kernel/check_nodes<<<>>>()
 * NAME
 *  check_nodes<<<>>>()
 * USAGE
 */
__global__ void check_nodes(int nparts, part_struct *parts, dom_struct *dom,
  real *theta, real *phi, int nnodes, BC bc);
/*
 * PURPOSE
 *  CUDA kernel to interpolate field variables to Lebedev quadrature nodes.
 * ARGUMENTS
 *  * nparts -- the number of particles
 *  * parts -- list of particles on this device
 *  * dom -- domain specification
 *  * theta -- theta component of list of Lebedev nodes
 *  * phi -- phi component of list of Lebedev nodes
 *  * nnodes -- the number of Lebedev nodes in the list
 *  * bc -- domain boundary conditions
 ******
 */

/****f* cuda_quadrature_kernel/interpolate_nodes<<<>>>()
 * NAME
 *  interpolate_nodes<<<>>>()
 * USAGE
 */
__global__ void interpolate_nodes(real *p0, real *p, real *u, real *v, real *w,
  real rho_f, real nu, gradP_struct gradP,
  part_struct *parts, dom_struct *dom, real *theta, real *phi, int nnodes,
  real *pp, real *ur, real *ut, real *up, real dt0, real dt, BC bc);
/*
 * PURPOSE
 *  CUDA kernel to interpolate field variables to Lebedev quadrature nodes.
 * ARGUMENTS
 *  * p -- device pressure field
 *  * u -- device x-component velocity field
 *  * v -- device y-component velocity field
 *  * w -- device z-component velocity field
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * gradP -- body force
 *  * parts -- list of particles on this device
 *  * dom -- domain specification
 *  * theta -- theta component of list of Lebedev nodes
 *  * phi -- phi component of list of Lebedev nodes
 *  * nnodes -- the number of Lebedev nodes in the list
 *  * pp -- the interpolated pressure field
 *  * ur -- the interpolated r-component of velocity field
 *  * ut -- the interpolated theta-component of velocity field
 *  * up -- the interpolated phi-component of velocity field
 ******
 */

/****f* cuda_quadrature_kernel/nnm<<<>>>()
 * NAME
 *  nnm<<<>>>()
 * USAGE
 */
__device__ real nnm(int n, int m);
/*
 * PURPOSE
 *  Compute spherical harmonic normalization N_nm.
 * ARGUMENTS
 *  * n -- degree
 *  * m -- order
 ******
 */

/****f* cuda_quadrature_kernel/pnm<<<>>>()
 * NAME
 *  pnm<<<>>>()
 * USAGE
 */
__device__ real pnm(int n, int m, real t);
/*
 * PURPOSE
 *  Compute associated Legendre function P_nm(theta).
 * ARGUMENTS
 *  * n -- degree
 *  * m -- order
 *  * t -- theta
 ******
 */

/****f* cuda_quadrature_kernel/cuda_get_coeffs<<<>>>()
 * NAME
 *  cuda_get_coeffs<<<>>>()
 * USAGE
 */
__global__ void cuda_get_coeffs(part_struct *parts,
  int *nn, int *mm, real *node_t, real *node_p,
  real *pp, real *ur, real *ut, real *up, real mu, real nu,
  int stride, real *pnm_re, real *pnm_im,
  real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im,
  real *int_Yp_re, real *int_Yp_im,
  real *int_rDYu_re, real *int_rDYu_im,
  real *int_xXDYu_re, real *int_xXDYu_im,
  int nnodes, int ncoeffs, real A1, real A2, real A3, real B,
  real *pnm_re0, real *pnm_im0,
  real *phinm_re0, real *phinm_im0,
  real *chinm_re0, real *chinm_im0,
  real lambrelax);
/*
 * PURPOSE
 *  Compute the Lamb's coefficients.
 * ARGUMENTS
 *  * TODO!
 ******
 */

/****f* cuda_quadrature_kernel/cuda_calc_forces<<<>>>()
 * NAME
 *  cuda_calc_forces<<<>>>()
 * USAGE
 */
__global__ void cuda_calc_forces(dom_struct *dom, part_struct *parts,
  int nparts, gradP_struct gradP,
  real rho_f, real mu, real nu, int stride,
  real *pnm_re, real *pnm_im,
  real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im);
/* 
 * PURPOSE
 *  Calculate the hydrodynamic forces on each particle.
 * ARGUMENTS
 *  * dom -- domain information
 *  * parts -- the particle data
 *  * nparts -- the number of particles
 *  * gradP -- body force
 *  * rho_f -- fluid density
 *  * mu -- fluid dynamic viscosity
 *  * nu -- fluid kinematic viscosity
 *  * stride -- the Lamb's coefficient array access stride length
 *  * pnm_re -- Lamb's coefficients
 *  * pnm_im -- Lamb's coefficients
 *  * phinm_re -- Lamb's coefficients
 *  * phinm_im -- Lamb's coefficients
 *  * chinm_re -- Lamb's coefficients
 *  * chinm_im -- Lamb's coefficients
 ******
 */

/****f* cuda_quadrature_kernel/compute_error<<<>>>()
 * NAME
 *  compute_error<<<>>>()
 * USAGE
 */
__global__ void compute_error(real lamb_cut, int stride, int nparts,
  real *pnm_re, real *pnm_re0, real *pnm_im, real *pnm_im0,
  real *phinm_re, real *phinm_re0, real *phinm_im, real *phinm_im0,
  real *chinm_re, real *chinm_re0, real *chinm_im, real *chinm_im0,
  real *coeffs, real *errors, real *part_errors, dom_struct *dom, real nu);
/*
 * PURPOSE
 *  Compute the error between the current and previous iteration Lamb's
 *  coefficient.
 * ARGUMENTS
 *  * lamb_cut -- the magnitude of errors below which to ignore, referenced to
 *    the error with the greatest magnitude
 *  * stride -- the Lamb's coefficient array access stride length
 *  * nparts -- the total number of particles
 *  * pnm_re -- Lamb's coefficients
 *  * pnm_re0 -- previous Lamb's coefficients
 *  * pnm_im -- Lamb's coefficients
 *  * pnm_im0 -- previous Lamb's coefficients
 *  * phinm_re -- Lamb's coefficients
 *  * phinm_re0 -- previous Lamb's coefficients
 *  * phinm_im -- Lamb's coefficients
 *  * phinm_im0 -- previous Lamb's coefficients
 *  * chinm_re -- Lamb's coefficients
 *  * chinm_re0 -- previous Lamb's coefficients
 *  * chinm_im -- Lamb's coefficients
 *  * chinm_im0 -- previous Lamb's coefficients
 *  * coeffs -- the sorted coefficients
 *  * errors -- the errors corresponding to the sorted coefficients
 *  * part_errors -- the maximum error for each particle
 ******
*/

#endif
