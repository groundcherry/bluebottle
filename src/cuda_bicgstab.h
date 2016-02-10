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

/****h* Bluebottle/cuda_bicgstab
 * NAME
 *  cuda_bicgstab
 * FUNCTION
 *  Bluebottle BiCGSTAB CUDA host functions.
 ******
 */
/****h* Bluebottle/cuda_bicgstab_kernel
 * NAME
 *  cuda_bicgstab_kernel
 * FUNCTION
 *  Bluebottle BiCGSTAB CUDA kernel functions.
 ******
 */

#ifndef _CUDA_BICGSTAB_H
#define _CUDA_BICGSTAB_H

extern "C"
{
#include "bluebottle.h"
#include "particle.h"
#include "recorder.h"
}

//#include <cusp/dia_matrix.h>

/****f* cuda_bicgstab/cuda_PP_rhs()
 * NAME
 *  cuda_PP_rhs()
 * USAGE
 */
extern "C"
void cuda_PP_rhs(int dev);
/*
 * FUNCTION
 *  Build the right-hand side of the pressure-Poisson problem.  Note that
 *  this is called from cuda_PP_bicgstab, which is already OMP threaded.  Do not
 *  OMP thread this function.
 * ARGUMENTS
 *  dev -- the device on which to operate
 ******
 */

/****f* cuda_bicgstab/cuda_ustar_rhs()
 * NAME
 *  cuda_ustar_rhs()
 * USAGE
 */
extern "C"
void cuda_ustar_rhs(int dev);
/*
 * FUNCTION
 *  Build the right-hand side of the u_star Helmholtz problem.  Note that
 *  this is called from cuda_ustar_helmholtz, which is already OMP threaded.  Do
 *  not OMP thread this function.
 * ARGUMENTS
 *  dev -- the device on which to operate
 ******
 */

/****f* cuda_bicgstab/cuda_vstar_rhs()
 * NAME
 *  cuda_vstar_rhs()
 * USAGE
 */
extern "C"
void cuda_vstar_rhs(int dev);
/*
 * FUNCTION
 *  Build the right-hand side of the v_star Helmholtz problem.  Note that
 *  this is called from cuda_vstar_helmholtz, which is already OMP threaded.  Do
 *  not OMP thread this function.
 * ARGUMENTS
 *  dev -- the device on which to operate
 ******
 */

/****f* cuda_bicgstab/cuda_wstar_rhs()
 * NAME
 *  cuda_wstar_rhs()
 * USAGE
 */
extern "C"
void cuda_wstar_rhs(int dev);
/*
 * FUNCTION
 *  Build the right-hand side of the w_star Helmholtz problem.  Note that
 *  this is called from cuda_wstar_helmholtz, which is already OMP threaded.  Do
 *  not OMP thread this function.
 * ARGUMENTS
 *  dev -- the device on which to operate
 ******
 */

/****f* cuda_bicgstab/cuda_wstar_rhs()
 * NAME
 *  cuda_wstar_rhs()
 * USAGE
 */
extern "C"
void cuda_wstar_rhs(int dev);
/*
 * FUNCTION
 *  Build the right-hand side of the w_star Helmholtz problem.  Note that
 *  this is called from cuda_wstar_helmholtz, which is already OMP threaded.  Do
 *  not OMP thread this function.
 * ARGUMENTS
 *  dev -- the device on which to operate
 ******
 */

/****f* cuda_bicgstab/cuda_div_U_launch()
 * NAME
 *  cuda_div_U_launch()
 * USAGE
 */
extern "C"
void cuda_div_U_launch(int dev, real *u, real *v, real *w, real *divU);
/*
 * FUNCTION
 *  Launch the divergence calculation.
 * ARGUMENTS
 *  * dev -- the current subdomain
 *  * u -- the device x-velocity field for computation
 *  * v -- the device y-velocity field for computation
 *  * w -- the device z-velocity field for computation
 *  * divU -- the device result
 ******
 */

/****f* cuda_bicgstab_kernel/div_U<<<>>>()
 * NAME
 *  div_U<<<>>>()
 * USAGE
 */
__global__ void div_U(real *u, real *v, real *w,
  real *out, dom_struct *dom);
/*
 * FUNCTION
 *  CUDA kernel for divergence calculation.
 * ARGUMENTS
 *  * u -- subdomain u-velocity
 *  * v -- subdomain v-velocity
 *  * w -- subdomain w-velocity
 *  * out -- the output divergence field
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bicgstab_kernel/PP_rhs<<<>>>()
 * NAME
 *  PP_rhs<<<>>>()
 * USAGE
 */
__global__ void PP_rhs(real rho_f, real *u_star, real *v_star, real *w_star,
  real *rhs, dom_struct *dom, real dt);
/*
 * FUNCTION
 *  Compute the right-hand side of the pressure Poisson problem.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * u_star -- device subdomain u-component intermediate velocity field
 *  * v_star -- device subdomain v-component intermediate velocity field
 *  * w_star -- device subdomain w-component intermediate velocity field
 *  * rhs -- the right-hand side array to build
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_rhs<<<>>>()
 * NAME
 *  ustar_rhs<<<>>>()
 * USAGE
 */
__global__ void ustar_rhs(real rho_f, real nu, real *u, real *v, real *w,
  real *p, real *f, real *conv0, real *conv, real *u_star, dom_struct *dom,
  real dt, real dt0);
/*
 * FUNCTION
 *  Compute the right-hand side of the u_star Helmholtz problem.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u -- device subdomain u-component velocity field
 *  * v -- device subdomain v-component velocity field
 *  * w -- device subdomain w-component velocity field
 *  * p -- device subdomain pressure field
 *  * f -- device subdomain body forcing field
 *  * conv0 -- device subdomain previous time step stored convective term
 *  * conv -- device subdomain time step stored convective term
 *  * u_star -- device subdomain u-component intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 *  * dt0 -- the previous timestep
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_rhs<<<>>>()
 * NAME
 *  vstar_rhs<<<>>>()
 * USAGE
 */
__global__ void vstar_rhs(real rho_f, real nu, real *u, real *v, real *w,
  real *p, real *f, real *conv0, real *conv, real *v_star, dom_struct *dom,
  real dt, real dt0);
/*
 * FUNCTION
 *  Compute the right-hand side of the v_star Helmholtz problem.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u -- device subdomain u-component velocity field
 *  * v -- device subdomain v-component velocity field
 *  * w -- device subdomain w-component velocity field
 *  * p -- device subdomain pressure field
 *  * f -- device subdomain body forcing field
 *  * conv0 -- device subdomain previous time step stored convective term
 *  * conv -- device subdomain time step stored convective term
 *  * v_star -- device subdomain v-component intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 *  * dt0 -- the previous timestep
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_rhs<<<>>>()
 * NAME
 *  wstar_rhs<<<>>>()
 * USAGE
 */
__global__ void wstar_rhs(real rho_f, real nu, real *u, real *v, real *w,
  real *p, real *f, real *conv0, real *conv, real *w_star, dom_struct *dom,
  real dt, real dt0);
/*
 * FUNCTION
 *  Compute the right-hand side of the w_star Helmholtz problem.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u -- device subdomain u-component velocity field
 *  * v -- device subdomain v-component velocity field
 *  * w -- device subdomain w-component velocity field
 *  * p -- device subdomain pressure field
 *  * f -- device subdomain body forcing field
 *  * conv0 -- device subdomain previous time step stored convective term
 *  * conv -- device subdomain time step stored convective term
 *  * w_star -- device subdomain w-component intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 *  * dt0 -- the previous timestep
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_init<<<>>>()
 * NAME
 *  coeffs_init<<<>>>()
 * USAGE
 */
__global__ void coeffs_init(dom_struct *dom, int pitch, real *values);
/*
 * FUNCTION
 *  Initialize the coefficient matrix to zeros before calling coeffs<<<>>>().
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_init<<<>>>()
 * NAME
 *  ustar_coeffs_init<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_init(dom_struct *dom, int pitch, real *values);
/*
 * FUNCTION
 *  Initialize the coefficient matrix to zeros before calling coeffs<<<>>>().
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_init<<<>>>()
 * NAME
 *  vstar_coeffs_init<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_init(dom_struct *dom, int pitch, real *values);
/*
 * FUNCTION
 *  Initialize the coefficient matrix to zeros before calling coeffs<<<>>>().
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_init<<<>>>()
 * NAME
 *  wstar_coeffs_init<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_init(dom_struct *dom, int pitch, real *values);
/*
 * FUNCTION
 *  Initialize the coefficient matrix to zeros before calling coeffs<<<>>>().
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs<<<>>>()
 * NAME
 *  coeffs<<<>>>()
 * USAGE
 */
__global__ void coeffs(dom_struct *dom, int *flag_u, int *flag_v, int *flag_w,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set up the coefficient matrix for the pressure-Poisson problem.  This
 *  routine applies zero-normal gradient boundary conditions to all boundaries
 *  (external and internal) while constructing the coefficient matrix.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * flag_u -- the x-direction velocity boundary flag
 *  * flag_v -- the y-direction velocity boundary flag
 *  * flag_w -- the z-direction velocity boundary flag
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs<<<>>>()
 * NAME
 *  ustar_coeffs<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs(real nu, real dt, dom_struct *dom, int pitch,
  real *values, int *flag_u, int *flag_v, int *flag_w);
/*
 * FUNCTION
 *  Set up the coefficient matrix for the u_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid kinematic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_u -- the x-direction velocity boundary flag
 *  * flag_v -- the y-direction velocity boundary flag
 *  * flag_w -- the z-direction velocity boundary flag
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs<<<>>>()
 * NAME
 *  vstar_coeffs<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs(real nu, real dt, dom_struct *dom, int pitch,
  real *values, int *flag_u, int *flag_v, int *flag_w);
/*
 * FUNCTION
 *  Set up the coefficient matrix for the v_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid kinematic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_u -- the x-direction velocity boundary flag
 *  * flag_v -- the y-direction velocity boundary flag
 *  * flag_w -- the z-direction velocity boundary flag
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs<<<>>>()
 * NAME
 *  wstar_coeffs<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs(real nu, real dt, dom_struct *dom, int pitch,
  real *values, int *flag_u, int *flag_v, int *flag_w);
/*
 * FUNCTION
 *  Set up the coefficient matrix for the w_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid kinematic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_u -- the x-direction velocity boundary flag
 *  * flag_v -- the y-direction velocity boundary flag
 *  * flag_w -- the z-direction velocity boundary flag
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_particles<<<>>>()
 * NAME
 *  ustar_coeffs_particles<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_particles(dom_struct *dom, int pitch, real *values,
  int *flag_u);
/*
 * FUNCTION
 *  Edit the coefficient matrix to set particle boundary conditions.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_u -- the x-direction velocity boundary flag
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_particles<<<>>>()
 * NAME
 *  vstar_coeffs_particles<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_particles(dom_struct *dom, int pitch, real *values,
  int *flag_v);
/*
 * FUNCTION
 *  Edit the coefficient matrix to set particle boundary conditions.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_v -- the y-direction velocity boundary flag
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_particles<<<>>>()
 * NAME
 *  wstar_coeffs_particles<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_particles(dom_struct *dom, int pitch, real *values,
  int *flag_u);
/*
 * FUNCTION
 *  Edit the coefficient matrix to set particle boundary conditions.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_w -- the z-direction velocity boundary flag
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_dirichlet_W<<<>>>()
 * NAME
 *  ustar_coeffs_dirichlet_W<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_dirichlet_W(real bc, real *u_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for u_star Helmholtz problem.
 * ARGUMENTS
 *  * bc -- the value of the Dirichlet boundary condition
 *  * u_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_dirichlet_E<<<>>>()
 * NAME
 *  ustar_coeffs_dirichlet_E<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_dirichlet_E(real bc, real *u_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for u_star Helmholtz problem.
 * ARGUMENTS
 *  * bc -- the value of the Dirichlet boundary condition
 *  * u_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_dirichlet_S<<<>>>()
 * NAME
 *  ustar_coeffs_dirichlet_S<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_dirichlet_S(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for u_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * u_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_dirichlet_N<<<>>>()
 * NAME
 *  ustar_coeffs_dirichlet_N<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_dirichlet_N(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for u_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * u_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_dirichlet_B<<<>>>()
 * NAME
 *  ustar_coeffs_dirichlet_B<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_dirichlet_B(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for u_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * u_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_dirichlet_T<<<>>>()
 * NAME
 *  ustar_coeffs_dirichlet_T<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_dirichlet_T(real nu, real *u_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for u_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * u_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_dirichlet_W<<<>>>()
 * NAME
 *  vstar_coeffs_dirichlet_W<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_dirichlet_W(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for v_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * v_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_dirichlet_E<<<>>>()
 * NAME
 *  vstar_coeffs_dirichlet_E<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_dirichlet_E(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for v_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * v_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_dirichlet_S<<<>>>()
 * NAME
 *  vstar_coeffs_dirichlet_S<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_dirichlet_S(real bc, real *v_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for v_star Helmholtz problem.
 * ARGUMENTS
 *  * bc -- the value of the Dirichlet boundary condition
 *  * v_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_dirichlet_N<<<>>>()
 * NAME
 *  vstar_coeffs_dirichlet_N<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_dirichlet_N(real bc, real *v_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for v_star Helmholtz problem.
 * ARGUMENTS
 *  * bc -- the value of the Dirichlet boundary condition
 *  * v_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_dirichlet_B<<<>>>()
 * NAME
 *  vstar_coeffs_dirichlet_B<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_dirichlet_B(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for v_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * v_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_dirichlet_T<<<>>>()
 * NAME
 *  vstar_coeffs_dirichlet_T<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_dirichlet_T(real nu, real *v_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for v_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * v_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_dirichlet_W<<<>>>()
 * NAME
 *  wstar_coeffs_dirichlet_W<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_dirichlet_W(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for w_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * w_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_dirichlet_E<<<>>>()
 * NAME
 *  wstar_coeffs_dirichlet_E<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_dirichlet_E(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for w_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * w_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_dirichlet_S<<<>>>()
 * NAME
 *  wstar_coeffs_dirichlet_S<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_dirichlet_S(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for w_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * w_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_dirichlet_N<<<>>>()
 * NAME
 *  wstar_coeffs_dirichlet_N<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_dirichlet_N(real nu, real *w_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for w_star Helmholtz problem.
 * ARGUMENTS
 *  * nu -- fluid viscosity
 *  * w_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_dirichlet_B<<<>>>()
 * NAME
 *  wstar_coeffs_dirichlet_B<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_dirichlet_B(real bc, real *w_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for w_star Helmholtz problem.
 * ARGUMENTS
 *  * bc -- the value of the Dirichlet boundary condition
 *  * w_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_dirichlet_T<<<>>>()
 * NAME
 *  wstar_coeffs_dirichlet_T<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_dirichlet_T(real bc, real *w_star, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set Dirichlet boundary conditions for w_star Helmholtz problem.
 * ARGUMENTS
 *  * bc -- the value of the Dirichlet boundary condition
 *  * w_star -- the right-hand side vector
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_periodic_W<<<>>>()
 * NAME
 *  coeffs_periodic_W<<<>>>()
 * USAGE
 */
__global__ void coeffs_periodic_W(dom_struct *dom, int pitch, real *values,
  int *flag_u);
/*
 * FUNCTION
 *  Set periodic boundary conditions for pressure for W face.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_u -- subdomain device flag
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_periodic_E<<<>>>()
 * NAME
 *  coeffs_periodic_E<<<>>>()
 * USAGE
 */
__global__ void coeffs_periodic_E(dom_struct *dom, int pitch, real *values,
  int *flag_u);
/*
 * FUNCTION
 *  Set periodic boundary conditions for pressure for E face.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_u -- subdomain device flag
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_periodic_S<<<>>>()
 * NAME
 *  coeffs_periodic_S<<<>>>()
 * USAGE
 */
__global__ void coeffs_periodic_S(dom_struct *dom, int pitch, real *values,
  int *flag_v);
/*
 * FUNCTION
 *  Set periodic boundary conditions for pressure for S face.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_v -- subdomain device flag
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_periodic_N<<<>>>()
 * NAME
 *  coeffs_periodic_N<<<>>>()
 * USAGE
 */
__global__ void coeffs_periodic_N(dom_struct *dom, int pitch, real *values,
  int *flag_v);
/*
 * FUNCTION
 *  Set periodic boundary conditions for pressure for N face.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_v -- subdomain device flag
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_periodic_B<<<>>>()
 * NAME
 *  coeffs_periodic_B<<<>>>()
 * USAGE
 */
__global__ void coeffs_periodic_B(dom_struct *dom, int pitch, real *values,
  int *flag_w);
/*
 * FUNCTION
 *  Set periodic boundary conditions for pressure for B face.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_w -- subdomain device flag
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_periodic_T<<<>>>()
 * NAME
 *  coeffs_periodic_T<<<>>>()
 * USAGE
 */
__global__ void coeffs_periodic_T(dom_struct *dom, int pitch, real *values,
  int *flag_w);
/*
 * FUNCTION
 *  Set periodic boundary conditions for pressure for T face.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 *  * flag_w -- subdomain device flag
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_periodic_W<<<>>>()
 * NAME
 *  ustar_coeffs_periodic_W<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_periodic_W(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for u_star for W face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_periodic_E<<<>>>()
 * NAME
 *  ustar_coeffs_periodic_E<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_periodic_E(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for u_star for E face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_periodic_S<<<>>>()
 * NAME
 *  ustar_coeffs_periodic_S<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_periodic_S(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for u_star for S face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_periodic_N<<<>>>()
 * NAME
 *  ustar_coeffs_periodic_N<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_periodic_N(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for u_star for N face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_periodic_B<<<>>>()
 * NAME
 *  ustar_coeffs_periodic_B<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_periodic_B(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for u_star for B face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/ustar_coeffs_periodic_T<<<>>>()
 * NAME
 *  ustar_coeffs_periodic_T<<<>>>()
 * USAGE
 */
__global__ void ustar_coeffs_periodic_T(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for u_star for T face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_periodic_W<<<>>>()
 * NAME
 *  vstar_coeffs_periodic_W<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_periodic_W(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for v_star for W face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_periodic_E<<<>>>()
 * NAME
 *  vstar_coeffs_periodic_E<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_periodic_E(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for v_star for E face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_periodic_S<<<>>>()
 * NAME
 *  vstar_coeffs_periodic_S<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_periodic_S(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for v_star for S face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_periodic_N<<<>>>()
 * NAME
 *  vstar_coeffs_periodic_N<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_periodic_N(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for v_star for N face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_periodic_B<<<>>>()
 * NAME
 *  vstar_coeffs_periodic_B<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_periodic_B(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for v_star for B face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/vstar_coeffs_periodic_T<<<>>>()
 * NAME
 *  vstar_coeffs_periodic_T<<<>>>()
 * USAGE
 */
__global__ void vstar_coeffs_periodic_T(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for v_star for T face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_periodic_W<<<>>>()
 * NAME
 *  wstar_coeffs_periodic_W<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_periodic_W(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for w_star for W face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_periodic_E<<<>>>()
 * NAME
 *  wstar_coeffs_periodic_E<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_periodic_E(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for w_star for E face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_periodic_S<<<>>>()
 * NAME
 *  wstar_coeffs_periodic_S<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_periodic_S(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for w_star for S face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_periodic_N<<<>>>()
 * NAME
 *  wstar_coeffs_periodic_N<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_periodic_N(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for w_star for N face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_periodic_B<<<>>>()
 * NAME
 *  wstar_coeffs_periodic_B<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_periodic_B(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for w_star for B face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/wstar_coeffs_periodic_T<<<>>>()
 * NAME
 *  wstar_coeffs_periodic_T<<<>>>()
 * USAGE
 */
__global__ void wstar_coeffs_periodic_T(real nu, real dt, dom_struct *dom,
  int pitch, real *values);
/*
 * FUNCTION
 *  Set periodic boundary conditions for w_star for T face.
 * ARGUMENTS
 *  * nu -- kinemtic viscosity
 *  * dt -- time step
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_particle<<<>>>()
 * NAME
 *  coeffs_particle<<<>>>()
 * USAGE
 */
__global__ void coeffs_particle(dom_struct *dom, int pitch, real *values,
  int *phase);
/*
 * FUNCTION
 *  Parse the matrix generated by coeffs<<<>>>() and place a one on the
 *  main diagonal of any row that has all zeros.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */

/****f* cuda_bicgstab_kernel/coeffs_zeros<<<>>>()
 * NAME
 *  coeffs_zeros<<<>>>()
 * USAGE
 */
__global__ void coeffs_zeros(dom_struct *dom, int pitch, real *values,
  real *rhs);
/*
 * FUNCTION
 *  Parse the matrix generated by coeffs<<<>>>() and place a one on the
 *  main diagonal of any row that has all zeros.
 * ARGUMENTS
 *  * dom -- the subdomain associated with this device
 *  * pitch -- the pitch associated with the CUSP dia_matrix values data
 *    structure
 *  * values -- the device pointer to the device-allocated CUSP dia_matrix
 *    values data structure (access it via:
 *    thrust::raw_pointer_cast(&A->values.values[0]))
 ******
 */


/****f* cuda_bicgstab_kernel/PP_BC_p_W<<<>>>()
 * NAME
 *  PP_BC_p_W<<<>>>()
 * USAGE
 */
__global__ void PP_BC_p_W(real *A, int pitch, dom_struct *dom);
/*
 * FUNCTION
 *  Apply zero-normal Neumann boundary conditions to the west face p for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * A -- the subdomain coefficient matrix
 *  * pitch -- the pitch associated with A
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bicgstab_kernel/PP_BC_p_E<<<>>>()
 * NAME
 *  PP_BC_p_E<<<>>>()
 * USAGE
 */
__global__ void PP_BC_p_E(real *A, int pitch, dom_struct *dom);
/*
 * FUNCTION
 *  Apply zero-normal Neumann boundary conditions to the east face p for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * A -- the subdomain coefficient matrix
 *  * pitch -- the pitch associated with A
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bicgstab_kernel/PP_BC_p_S<<<>>>()
 * NAME
 *  PP_BC_p_S<<<>>>()
 * USAGE
 */
__global__ void PP_BC_p_S(real *A, int pitch, dom_struct *dom);
/*
 * FUNCTION
 *  Apply zero-normal Neumann boundary conditions to the south face p for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * A -- the subdomain coefficient matrix
 *  * pitch -- the pitch associated with A
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bicgstab_kernel/PP_BC_p_N<<<>>>()
 * NAME
 *  PP_BC_p_N<<<>>>()
 * USAGE
 */
__global__ void PP_BC_p_N(real *A, int pitch, dom_struct *dom);
/*
 * FUNCTION
 *  Apply zero-normal Neumann boundary conditions to the north face p for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * A -- the subdomain coefficient matrix
 *  * pitch -- the pitch associated with A
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bicgstab_kernel/PP_BC_p_B<<<>>>()
 * NAME
 *  PP_BC_p_B<<<>>>()
 * USAGE
 */
__global__ void PP_BC_p_B(real *A, int pitch, dom_struct *dom);
/*
 * FUNCTION
 *  Apply zero-normal Neumann boundary conditions to the bottom face p for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * A -- the subdomain coefficient matrix
 *  * pitch -- the pitch associated with A
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bicgstab_kernel/PP_BC_p_T<<<>>>()
 * NAME
 *  PP_BC_p_T<<<>>>()
 * USAGE
 */
__global__ void PP_BC_p_T(real *A, int pitch, dom_struct *dom);
/*
 * FUNCTION
 *  Apply zero-normal Neumann boundary conditions to the top face p for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * A -- the subdomain coefficient matrix
 *  * pitch -- the pitch associated with A
 *  * dom -- the device subdomain
 ******
 */

#endif
