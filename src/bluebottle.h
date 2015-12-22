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

/****h* Bluebottle/bluebottle
 * NAME
 *  bluebottle
 * FUNCTION
 *  Bluebottle main execution code and global variable declarations.
 ******
 */

#ifndef _BLUEBOTTLE_H
#define _BLUEBOTTLE_H

#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef DOUBLE
  typedef double real;
#else
  typedef float real;
#endif

#include "domain.h"
#include "recorder.h"
#include "vtk.h"

/****d* bluebottle/PI
 * NAME
 *  PI
 * TYPE
 */
#define PI 3.14159265358979323846
/*
 * PURPOSE
 *  Define the constant pi.
 ******
 */

/****d* blueblottle/DIV_ST
 * NAME
 *  DIV_ST
 * TYPE
 */
#define DIV_ST 1e-10
/*
 * PURPOSE
 *  Define a value to use for fudging the value of theta when a division by
 *  sin(theta) occurs.
 ******
 */

/****d* bluebottle/FILE_NAME_SIZE
 * NAME
 *  FILE_NAME_SIZE
 * TYPE
 */
#define FILE_NAME_SIZE 256
/*
 * PURPOSE
 *  Define the maximum length of a file name.
 ******
 */

/****d* bluebottle/CHAR_BUF_SIZE
 * NAME
 *  CHAR_BUF_SIZE
 * TYPE
 */
#define CHAR_BUF_SIZE 256
/*
 * PURPOSE
 *  Define the maximum length of a character buffer read.
 ******
 */

/****d* bluebottle/ROOT_DIR
 * NAME
 *  ROOT_DIR
 * TYPE
 */
#define ROOT_DIR "."
/*
 * PURPOSE
 *  Define the root directory for the project.
 ******
 */

/****d* bluebottle/OUTPUT_DIR
 * NAME
 *  OUTPUT_DIR
 * TYPE
 */
#define OUTPUT_DIR "./output/"
/*
 * PURPOSE
 *  Define the output directory for the project.
 ******
 */

/****d* bluebottle/INPUT_DIR
 * NAME
 *  INPUT_DIR 
 * TYPE
 */
#define INPUT_DIR "./input/"
/*
 * PURPOSE
 *  Define the input directory for the project.
 ******
 */

/****d* bluebottle/DOM_BUF
 * NAME
 *  DOM_BUF
 * TYPE
 */
#define DOM_BUF 1
/*
 * PURPOSE
 *  Define the size of the domain boundary condition ghost cell buffer (the
 *  number of ghost cells on one side of a give domain direction).
 ******
 */

/****d* bluebottle/PERIODIC
 * NAME
 *  PERIODIC
 * TYPE
 */
#define PERIODIC 0
/*
 * PURPOSE
 *  Define the periodic boundary condition type.
 ******
 */

/****d* bluebottle/DIRICHLET
 * NAME
 *  DIRICHLET
 * TYPE
 */
#define DIRICHLET 1
/*
 * PURPOSE
 *  Define the Dirichlet boundary condition type.
 ******
 */

/****d* bluebottle/NEUMANN
 * NAME
 *  NEUMANN
 * TYPE
 */
#define NEUMANN 2
/*
 * PURPOSE
 *  Define the Neumann boundary condition type.
 ******
 */

/****d* bluebottle/PRECURSOR
 * NAME
 *  PRECURSOR
 * TYPE
 */
#define PRECURSOR 3
/*
 * PURPOSE
 *  Define the turbulence precursor domain boundary condition type. This
 *  type of boundary is treated like a Dirichlet boundary, but takes its value
 *  from the precursor domain.
 ******
 */

/****d* bluebottle/HOMOGENEOUS
 * NAME
 *  HOMOGENEOUS
 * TYPE
 */
#define HOMOGENEOUS 10
/*
 * PURPOSE
 *  Define the homogeneous outflow plane condition.
 ******
 */

/****d* bluebottle/WEST
 * NAME
 *  WEST
 * TYPE
 */
#define WEST 0
/*
 * PURPOSE
 *  Define the West boundary.
 ******
 */

/****d* bluebottle/EAST
 * NAME
 *  EAST
 * TYPE
 */
#define EAST 1
/*
 * PURPOSE
 *  Define the East boundary.
 ******
 */

/****d* bluebottle/SOUTH
 * NAME
 *  SOUTH
 * TYPE
 */
#define SOUTH 2
/*
 * PURPOSE
 *  Define the South boundary.
 ******
 */

/****d* bluebottle/NORTH
 * NAME
 *  NORTH
 * TYPE
 */
#define NORTH 3
/*
 * PURPOSE
 *  Define the North boundary.
 ******
 */

/****d* bluebottle/BOTTOM
 * NAME
 *  BOTTOM
 * TYPE
 */
#define BOTTOM 4
/*
 * PURPOSE
 *  Define the Bottom boundary.
 ******
 */

/****d* bluebottle/TOP
 * NAME
 *  TOP
 * TYPE
 */
#define TOP 5
/*
 * PURPOSE
 *  Define the Top boundary.
 ******
 */

/****d* bluebottle/MAX_THREADS_1D
 * NAME
 *  MAX_THREADS_1D
 * TYPE
 */
#define MAX_THREADS_1D 128
/*
 * PURPOSE
 *  Define the maximum number of threads per block on a CUDA device.  Must be
 *  hardcoded, but does depend on the particular device being used.
 ******
 */

/****d* bluebottle/MAX_THREADS_DIM
 * NAME
 *  MAX_THREADS_DIM
 * TYPE
 */
#define MAX_THREADS_DIM 16
/*
 * PURPOSE
 *  Define the maximum number of threads per dimension per block on a CUDA
 *  device.  Must be hardcoded, but does depend on the particular device being
 *  used.
 ******
 */

/****d* bluebottle/MASTER
 * NAME
 *  MASTER
 * TYPE
 */
#define MASTER 0
/*
 * PURPOSE
 *  Define the master process in MPI communication.
 ******
 */

/****d* bluebottle/QUIESCENT
 * NAME
 *  QUIESCENT
 * TYPE
 */
#define QUIESCENT 0
/*
 * PURPOSE
 *  Define the initial condition quiescent flow case.
 ******
 */

/****d* bluebottle/SHEAR
 * NAME
 *  SHEAR
 * TYPE
 */
#define SHEAR 1
/*
 * PURPOSE
 *  Define the initial condition shear case.
 ******
 */

/****d* bluebottle/CHANNEL
 * NAME
 *  CHANNEL
 * TYPE
 */
#define CHANNEL 2
/*
 * PURPOSE
 *  Define the initial condition channel case.
 ******
 */

/****v* bluebottle/dev_start
 * NAME
 *  dev_start
 * TYPE
 */
extern int dev_start;
/*
 * PURPOSE
 *  The GPU device index list start device (inclusive).  To use devices three
 *  and four, e.g., use dev_start = 3, dev_end = 4.
 ******
 */

/****v* bluebottle/dev_end
 * NAME
 *  dev_end
 * TYPE
 */
extern int dev_end;
/*
 * PURPOSE
 *  The GPU device index list end device (inclusive).  To use devices three
 *  and four, e.g., use dev_start = 3, dev_end = 4.
 ******
 */

/****v* bluebottle/dom
 * NAME
 *  dom
 * TYPE
 */
extern dom_struct *dom;
/*
 * PURPOSE
 *  Carry GPU domain decomposition subdomain information.  Contains an array
 *  of dom_struct structs that each represent one subdomain.
 ******
 */

/****v* bluebottle/_dom
 * NAME
 *  _dom
 * TYPE
 */
extern dom_struct **_dom;
/*
 * PURPOSE
 *  CUDA device analog for dom.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* bluebottle/Dom
 * NAME
 *  Dom
 * TYPE
 */
extern dom_struct Dom;
/*
 * PURPOSE
 *  Carry total domain information (instance of dom_struct).
 ******
 */

extern real *p0;
extern real *phi;

/****v* bluebottle/p
 * NAME
 *  p
 * TYPE
 */
extern real *p;
/*
 * PURPOSE
 *  Pressure field vector (grid type Gcc; x-component varies first, then
 *  y-component, then z-component).
 ******
 */

extern real **_p0;
extern real **_phi;

/****v* bluebottle/_p
 * NAME
 *  _p
 * TYPE
 */
extern real **_p;
/*
 * PURPOSE
 *  CUDA device analog for p.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* bluebottle/divU
 * NAME
 *  divU
 * TYPE
 */
extern real *divU;
/*
 * PURPOSE
 *  Velocity divergence field vector.
 ******
 */

/****v* bluebottle/_divU
 * NAME
 *  _divU
 * TYPE
 */
extern real **_divU;
/*
 * PURPOSE
 *  CUDA device analog for divU.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* bluebottle/u
 * NAME
 *  u
 * TYPE
 */
extern real *u;
/*
 * PURPOSE
 *  Velocity field vector u-component (grid type Gfx; x-component varies
 *  first, then y-component, then z-component).
 ******
 */

/****v* bluebottle/_u
 * NAME
 *  _u
 * TYPE
 */
extern real **_u;
/*
 * PURPOSE
 *  CUDA device analog for u.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* bluebottle/u0
 * NAME
 *  u0
 * TYPE
 */
extern real *u0;
/*
 * PURPOSE
 *  Host u stored from the previous timestep.
 ******
 */

/****v* bluebottle/_u0
 * NAME
 *  _u0
 * TYPE
 */
extern real **_u0;
/*
 * PURPOSE
 *  CUDA device analog for u stored from the previous timestep.  It contains
 *  pointers to arrays containing the subdomain fields on which each device
 *  operates.
 ******
 */

/****v* bluebottle/v
 * NAME
 *  v
 * TYPE
 */
extern real *v;
/*
 * PURPOSE
 *  Velocity field vector v-component (grid type Gfy; x-component varies
 *  first, then y-component, then z-component).
 ******
 */

/****v* bluebottle/_v
 * NAME
 *  _v
 * TYPE
 */
extern real **_v;
/*
 * PURPOSE
 *  CUDA device analog for v.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* bluebottle/v0
 * NAME
 *  v0
 * TYPE
 */
extern real *v0;
/*
 * PURPOSE
 *  Host v stored from the previous timestep.
 ******
 */

/****v* bluebottle/_v0
 * NAME
 *  _v0
 * TYPE
 */
extern real **_v0;
/*
 * PURPOSE
 *  CUDA device analog for v stored from the previous timestep.  It contains
 *  pointers to arrays containing the subdomain fields on which each device
 *  operates.
 ******
 */

/****v* bluebottle/w
 * NAME
 *  w
 * TYPE
 */
extern real *w;
/*
 * PURPOSE
 *  Velocity field vector w-component (grid type Gfz; x-component varies
 *  first, then y-component, then z-component).
 ******
 */

/****v* bluebottle/_w
 * NAME
 *  _w
 * TYPE
 */
extern real **_w;
/*
 * PURPOSE
 *  CUDA device analog for w.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* bluebottle/w0
 * NAME
 *  w0
 * TYPE
 */
extern real *w0;
/*
 * PURPOSE
 *  Host w stored from the previous timestep.
 ******
 */

/****v* bluebottle/_w0
 * NAME
 *  _w0
 * TYPE
 */
extern real **_w0;
/*
 * PURPOSE
 *  CUDA device analog for w stored from the previous timestep.  It contains
 *  pointers to arrays containing the subdomain fields on which each device
 *  operates.
 ******
 */

/****v* bluebottle/f_x
 * NAME
 *  f_x
 * TYPE
 */
extern real *f_x;
/*
 * PURPOSE
 *  The body forcing in the x-direction.
 ******
 */

/****v* bluebottle/_f_x
 * NAME
 *  _f_x
 * TYPE
 */
extern real **_f_x;
/*
 * PURPOSE
 *  CUDA device analog for f_x.  It contains pointers to arrays containing the
 *  subdomain fiels on which each device operates.
 ******
 */

/****v* bluebottle/f_y
 * NAME
 *  f_y
 * TYPE
 */
extern real *f_y;
/*
 * PURPOSE
 *  The body forcing in the y-direction.
 ******
 */

/****v* bluebottle/_f_y
 * NAME
 *  _f_y
 * TYPE
 */
extern real **_f_y;
/*
 * PURPOSE
 *  CUDA device analog for f_y.  It contains pointers to arrays containing the
 *  subdomain fiels on which each device operates.
 ******
 */

/****v* bluebottle/f_z
 * NAME
 *  f_z
 * TYPE
 */
extern real *f_z;
/*
 * PURPOSE
 *  The body forcing in the z-direction.
 ******
 */

/****v* bluebottle/_f_z
 * NAME
 *  _f_z
 * TYPE
 */
extern real **_f_z;
/*
 * PURPOSE
 *  CUDA device analog for f_z.  It contains pointers to arrays containing the
 *  subdomain fiels on which each device operates.
 ******
 */

/****v* bluebottle/diff0_u
 * NAME
 *  diff0_u
 * TYPE
 */
extern real *diff0_u;
/*
 * PURPOSE
 *  Host array to store the x-component diffusion solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/diff0_v
 * NAME
 *  diff0_v
 * TYPE
 */
extern real *diff0_v;
/*
 * PURPOSE
 *  Host array to store the y-component diffusion solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/diff0_w
 * NAME
 *  diff0_w
 * TYPE
 */
extern real *diff0_w;
/*
 * PURPOSE
 *  Host array to store the z-component diffusion solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/conv0_u
 * NAME
 *  conv0_u
 * TYPE
 */
extern real *conv0_u;
/*
 * PURPOSE
 *  Host array to store the previous x-component convection solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/conv0_v
 * NAME
 *  conv0_v
 * TYPE
 */
extern real *conv0_v;
/*
 * PURPOSE
 *  Host array to store the previous y-component convection solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/conv0_w
 * NAME
 *  conv0_w
 * TYPE
 */
extern real *conv0_w;
/*
 * PURPOSE
 *  Host array to store the previous z-component convection solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/_diff0_u
 * NAME
 *  _diff0_u
 * TYPE
 */
extern real **_diff0_u;
/*
 * PURPOSE
 *  CUDA device array to store the x-component diffusion solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_diff0_v
 * NAME
 *  _diff0_v
 * TYPE
 */
extern real **_diff0_v;
/*
 * PURPOSE
 *  CUDA device array to store the y-component diffusion solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_diff0_w
 * NAME
 *  _diff0_w
 * TYPE
 */
extern real **_diff0_w;
/*
 * PURPOSE
 *  CUDA device array to store the z-component diffusion solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_conv0_u
 * NAME
 *  _conv0_u
 * TYPE
 */
extern real **_conv0_u;
/*
 * PURPOSE
 *  CUDA device array to store the previous x-component convection solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_conv0_v
 * NAME
 *  _conv0_v
 * TYPE
 */
extern real **_conv0_v;
/*
 * PURPOSE
 *  CUDA device array to store the previous y-component convection solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_conv0_w
 * NAME
 *  _conv0_w
 * TYPE
 */
extern real **_conv0_w;
/*
 * PURPOSE
 *  CUDA device array to store the previous z-component convection solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/diff_u
 * NAME
 *  diff_u
 * TYPE
 */
extern real *diff_u;
/*
 * PURPOSE
 *  Host array to store the x-component diffusion solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/diff_v
 * NAME
 *  diff_v
 * TYPE
 */
extern real *diff_v;
/*
 * PURPOSE
 *  Host array to store the y-component diffusion solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/diff_w
 * NAME
 *  diff_w
 * TYPE
 */
extern real *diff_w;
/*
 * PURPOSE
 *  Host array to store the z-component diffusion solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/conv_u
 * NAME
 *  conv_u
 * TYPE
 */
extern real *conv_u;
/*
 * PURPOSE
 *  Host array to store the previous x-component convection solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/conv_v
 * NAME
 *  conv_v
 * TYPE
 */
extern real *conv_v;
/*
 * PURPOSE
 *  Host array to store the previous y-component convection solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/conv_w
 * NAME
 *  conv_w
 * TYPE
 */
extern real *conv_w;
/*
 * PURPOSE
 *  Host array to store the previous z-component convection solution
 *  for use in the next Adams-Bashforth step.
 ******
 */

/****v* bluebottle/_diff_u
 * NAME
 *  _diff_u
 * TYPE
 */
extern real **_diff_u;
/*
 * PURPOSE
 *  CUDA device array to store the x-component diffusion solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_diff_v
 * NAME
 *  _diff_v
 * TYPE
 */
extern real **_diff_v;
/*
 * PURPOSE
 *  CUDA device array to store the y-component diffusion solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_diff_w
 * NAME
 *  _diff_w
 * TYPE
 */
extern real **_diff_w;
/*
 * PURPOSE
 *  CUDA device array to store the z-component diffusion solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_conv_u
 * NAME
 *  _conv_u
 * TYPE
 */
extern real **_conv_u;
/*
 * PURPOSE
 *  CUDA device array to store the previous x-component convection solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_conv_v
 * NAME
 *  _conv_v
 * TYPE
 */
extern real **_conv_v;
/*
 * PURPOSE
 *  CUDA device array to store the previous y-component convection solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_conv_w
 * NAME
 *  _conv_w
 * TYPE
 */
extern real **_conv_w;
/*
 * PURPOSE
 *  CUDA device array to store the previous z-component convection solution
 *  for use in the next Adams-Bashforth step.  It contains pointers to arrays
 *  containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_u_star
 * NAME
 *  _u_star
 * TYPE
 */
extern real **_u_star;
/*
 * PURPOSE
 *  CUDA device array for the x-direction intermediate velocity.  It contains
 *  pointers to arrays containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_v_star
 * NAME
 *  _v_star
 * TYPE
 */
extern real **_v_star;
/*
 * PURPOSE
 *  CUDA device array for the y-direction intermediate velocity.  It contains
 *  pointers to arrays containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/_w_star
 * NAME
 *  _w_star
 * TYPE
 */
extern real **_w_star;
/*
 * PURPOSE
 *  CUDA device array for z-direction intermediate velocity.  It contains
 *  pointers to arrays containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/u_star
 * NAME
 *  u_star
 * TYPE
 */
extern real *u_star;
/*
 * PURPOSE
 *  X-direction intermediate velocity.  It contains
 *  pointers to arrays containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/v_star
 * NAME
 *  v_star
 * TYPE
 */
extern real *v_star;
/*
 * PURPOSE
 *  Y-direction intermediate velocity.  It contains
 *  pointers to arrays containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/w_star
 * NAME
 *  w_star
 * TYPE
 */
extern real *w_star;
/*
 * PURPOSE
 *  Z-direction intermediate velocity.  It contains
 *  pointers to arrays containing the subdomain fields to be stored.
 ******
 */

/****v* bluebottle/u_WE
 * NAME
 *  u_WE
 * TYPE
 */
extern real *u_WE;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity W or E face.
 ******
 */

/****v* bluebottle/u_SN
 * NAME
 *  u_SN
 * TYPE
 */
extern real *u_SN_S;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity S or N face, southern plane
 ******
 */

/****v* bluebottle/u_SN
 * NAME
 *  u_SN
 * TYPE
 */
extern real *u_SN_N;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity S or N face, northern plane
 ******
 */

/****v* bluebottle/u_BT_B
 * NAME
 *  u_BT_B
 * TYPE
 */
extern real *u_BT_B;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity B or T face on ghost plane.
 ******
 */

/****v* bluebottle/u_BT_T
 * NAME
 *  u_BT_T
 * TYPE
 */
extern real *u_BT_T;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity B or T face on first layer plane.
 ******
 */


/****v* bluebottle/v_WE
 * NAME
 *  v_WE
 * TYPE
 */
extern real *v_WE_W;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity W or E face.
 ******
 */

/****v* bluebottle/v_WE
 * NAME
 *  v_WE
 * TYPE
 */
extern real *v_WE_E;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity W or E face.
 ******
 */

/****v* bluebottle/v_SN
 * NAME
 *  v_SN
 * TYPE
 */
extern real *v_SN;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity S or N face.
 ******
 */

/****v* bluebottle/v_BT_B
 * NAME
 *  v_BT_B
 * TYPE
 */
extern real *v_BT_B;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity B or T face on ghost cell.
 ******
 */

/****v* bluebottle/v_BT_T
 * NAME
 *  v_BT_T
 * TYPE
 */
extern real *v_BT_T;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity B or T face on the first layer.
 ******
 */


/****v* bluebottle/w_WE
 * NAME
 *  w_WE
 * TYPE
 */
extern real *w_WE_W;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity W or E face.
 ******
 */

/****v* bluebottle/w_WE
 * NAME
 *  w_WE
 * TYPE
 */
extern real *w_WE_E;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity W or E face.
 ******
 */

/****v* bluebottle/w_SN
 * NAME
 *  w_SN
 * TYPE
 */
extern real *w_SN_S;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity S or N face.
 ******
 */

/****v* bluebottle/w_SN
 * NAME
 *  w_SN
 * TYPE
 */
extern real *w_SN_N;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity S or N face.
 ******
 */

/****v* bluebottle/w_BT
 * NAME
 *  w_BT
 * TYPE
 */
extern real *w_BT;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity B or T face.
 ******
 */

/****v* bluebottle/_u_WE
 * NAME
 *  _u_WE
 * TYPE
 */
extern real **_u_WE;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity W or E face.
 *  Device version.
 ******
 */

/****v* bluebottle/_u_SN
 * NAME
 *  _u_SN
 * TYPE
 */
extern real **_u_SN_S;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity S or N face.
 *  Device version.
 ******
 */
/****v* bluebottle/_u_SN
 * NAME
 *  _u_SN
 * TYPE
 */
extern real **_u_SN_N;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity S or N face.
 *  Device version.
 ******
 */

/****v* bluebottle/_u_BT_B
 * NAME
 *  _u_BT_B
 * TYPE
 */
extern real **_u_BT_B;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity B or T face on ghost cell.
 *  Device version.
 ******
 */

/****v* bluebottle/_u_BT_T
 * NAME
 *  _u_BT_T
 * TYPE
 */
extern real **_u_BT_T;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the u-velocity B or T face on the first layer.
 *  Device version.
 ******
 */


/****v* bluebottle/_v_WE
 * NAME
 *  _v_WE
 * TYPE
 */
extern real **_v_WE_W;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity W or E face.
 *  Device version.
 ******
 */

/****v* bluebottle/_v_WE
 * NAME
 *  _v_WE
 * TYPE
 */
extern real **_v_WE_E;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity W or E face.
 *  Device version.
 ******
 */

/****v* bluebottle/_v_SN
 * NAME
 *  _v_SN
 * TYPE
 */
extern real **_v_SN;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity S or N face.
 *  Device version.
 ******
 */

/****v* bluebottle/_v_BT_B
 * NAME
 *  _v_BT_B
 * TYPE
 */
extern real **_v_BT_B;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity B or T face on ghost cell.
 *  Device version.
 ******
 */

/****v* bluebottle/_v_BT_T
 * NAME
 *  _v_BT_T
 * TYPE
 */
extern real **_v_BT_T;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the v-velocity B or T face on the first layer cell.
 *  Device version.
 ******
 */

/****v* bluebottle/_w_WE
 * NAME
 *  _w_WE
 * TYPE
 */
extern real **_w_WE_W;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity W or E face.
 *  Device version.
 ******
 */

/****v* bluebottle/_w_WE
 * NAME
 *  _w_WE
 * TYPE
 */
extern real **_w_WE_E;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity W or E face.
 *  Device version.
 ******
 */

/****v* bluebottle/_w_SN
 * NAME
 *  _w_SN
 * TYPE
 */
extern real **_w_SN_S;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity S or N face.
 *  Device version.
 ******
 */

/****v* bluebottle/_w_SN
 * NAME
 *  _w_SN
 * TYPE
 */
extern real **_w_SN_N;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity S or N face.
 *  Device version.
 ******
 */

/****v* bluebottle/_w_BT
 * NAME
 *  _w_BT
 * TYPE
 */
extern real **_w_BT;
/*
 * PURPOSE
 *  Hold the turbulence precursor plane for the w-velocity B or T face.
 *  Device version.
 ******
 */

/****v* bluebottle/_rhs_p
 * NAME
 *  _rhs_p
 * TYPE
 */
extern real **_rhs_p;
/*
 * PURPOSE
 *  CUDA device array for storing the right-hand side of the pressure-Poisson
 *  problem.  It contains pointers to arrays containing the subdomain fields
 *  on which each device operates.
 ******
 */

/****v* bluebottle/duration
 * NAME
 *  duration
 * TYPE
 */
extern real duration;
/*
 * PURPOSE
 *  The duration of the simulation (i.e., stop time).
 ******
 */

/****v* bluebottle/ttime
 * NAME
 *  ttime
 * TYPE
 */
extern real ttime;
/*
 * PURPOSE
 *  The accumulated time since the simulation began.
 ******
 */

/****v* bluebottle/vel_tDelay
  * NAME
  * vel_tDelay
  * TYPE
  */
extern real vel_tDelay;
/*
  * PURPOSE
  *  The time when the dirichlet velocity boundary conditions should be applied
  ******
  */

/****v* bluebottle/p_tDelay
  * NAME
  * p_tDelay
  * TYPE
  */
extern real p_tDelay;
/*
  * PURPOSE
  *  The time when the applied pressure should be applied
  ******
  */
 
/****v* bluebottle/g_tDelay
  * NAME
  * g_tDelay
  * TYPE
  */
extern real g_tDelay;
/*
  * PURPOSE
  *  The time when the applied gravity should be applied
  ******
  */

/****v* bluebottle/dt
 * NAME
 *  dt
 * TYPE
 */
extern real dt;
/*
 * PURPOSE
 *  The current timestep size.  Upon initialization via input file, dt is the
 *  timestep size used for the first timestep.  It is subsequently updated
 *  dynamically for the remainder of the simulation.
 ******
 */

/****v* bluebottle/dt0
 * NAME
 *  dt0
 * TYPE
 */
extern real dt0;
/*
 * PURPOSE
 *  The previous timestep size.
 ******
 */

/****v* bluebottle/CFL
 * NAME
 *  CFL
 * TYPE
 */
extern real CFL;
/*
 * PURPOSE
 *  Define the CFL condition number for adaptive timestep calculation.
 ******
 */

/****v* bluebottle/pp_max_iter
 * NAME
 *  pp_max_iter
 * TYPE
 */
extern int pp_max_iter;
/*
 * PURPOSE
 *  The maximum number of iterations for the pressure-Poisson problem solver.
 ******
 */

/****v* bluebottle/pp_residual
 * NAME
 *  pp_residual
 * TYPE
 */
extern real pp_residual;
/*
 * PURPOSE
 *  The maximum desired residual for the pressure-Poisson problem solver.
 ******
 */

/****v* bluebottle/lamb_max_iter
 * NAME
 *  lamb_max_iter
 * TYPE
 */
extern int lamb_max_iter;
/*
 * PURPOSE
 *  The maximum number of iterations for the Lamb's coefficient convergence.
 ******
 */

/****v* bluebottle/lamb_residual
 * NAME
 *  lamb_residual
 * TYPE
 */
extern real lamb_residual;
/*
 * PURPOSE
 *  The maximum desired residual for the Lamb's coefficient iteration process.
 ******
 */

/****v* bluebottle/lamb_relax
 * NAME
 *  lamb_relax
 * TYPE
 */
extern real lamb_relax;
/*
 * PURPOSE
 *  The underrelaxation factor for the Lamb's coefficient iteration process.
 *  Zero weights toward old value and unity weights toward new value.
 ******
 */

/****v* bluebottle/lamb_cut
 * NAME
 *  lamb_cut
 * TYPE
 */
extern real lamb_cut;
/*
 * PURPOSE
 *  The magnitude below which errors in Lamb's coefficients are ignored,
 *  compared to the coefficient with greatest magnitude. The lower this number,
 *  the more coefficients will be considered important when computing the error.
 *  To improve convergence rate, decrease this number. It should never be
 *  greater than 1e-2.
 ******
 */

/****v* bluebottle/out_plane
 * NAME
 *  out_plane
 * TYPE
 */
extern int out_plane;
/*
 * PURPOSE
 *  Define which plane is the outlet plane.
 ******
 */

/****v* bluebottle/stepnum
 * NAME
 *  stepnum
 * TYPE
 */
extern int stepnum;
/*
 * PURPOSE
 *  The current timestep number for the simulation.  The initial configuration
 *  is given by stepnum = 0.
 ******
 */

/****v* bluebottle/rec_flow_field_stepnum_out
 * NAME
 *  rec_flow_field_stepnum_out
 * TYPE
 */
extern int rec_flow_field_stepnum_out;
/*
 * PURPOSE
 *  The output timestep number for the simulation.  The initial configuration
 *  is given by stepnum = 0.
 ******
 */

/****v* bluebottle/rec_prec_flow_field_stepnum_out
 * NAME
 *  rec_prec_flow_field_stepnum_out
 * TYPE
 */
extern int rec_prec_flow_field_stepnum_out;
/*
 * PURPOSE
 *  The output timestep number for the simulation.  The initial configuration
 *  is given by stepnum = 0.
 ******
 */

/****v* bluebottle/rec_paraview_stepnum_out
 * NAME
 *  rec_paraview_stepnum_out
 * TYPE
 */
extern int rec_paraview_stepnum_out;
/*
 * PURPOSE
 *  The output timestep number for the simulation.  The initial configuration
 *  is given by stepnum = 0.
 ******
 */

/****v* bluebottle/rec_prec_stepnum_out
 * NAME
 *  rec_prec_stepnum_out
 * TYPE
 */
extern int rec_prec_stepnum_out;
/*
 * PURPOSE
 *  The output timestep number for the precursor simulation.  The initial
 *  configuration is given by stepnum = 0.
 ******
 */

/****v* bluebottle/rec_particle_stepnum_out
 * NAME
 *  rec_particle_stepnum_out
 * TYPE
 */
extern int rec_particle_stepnum_out;
/*
 * PURPOSE
 *  The output timestep number for the simulation.  The initial configuration
 *  is given by stepnum = 0.
 ******
 */

/****v* bluebottle/rec_flow_field_dt
 * NAME
 *  rec_flow_field_dt
 * TYPE
 */
extern real rec_flow_field_dt;
/*
 * PURPOSE
 *  Recorder flow field output timestep size.
 ******
 */

/****v* bluebottle/rec_flow_field_ttime_out
 * NAME
 *  rec_flow_field_ttime_out
 * TYPE
 */
extern real rec_flow_field_ttime_out;
/*
 * PURPOSE
 *  Recorder flow field output time since last output.
 ******
 */

/****v* bluebottle/rec_flow_field_vel
 * NAME
 *  rec_flow_field_vel
 * TYPE
 */
extern int rec_flow_field_vel;
/*
 * PURPOSE
 *  Recorder flow field output velocity field boolean.
 ******
 */

/****v* bluebottle/rec_flow_field_p
 * NAME
 *  rec_flow_field_p
 * TYPE
 */
extern int rec_flow_field_p;
/*
 * PURPOSE
 *  Recorder flow field output pressure field boolean.
 ******
 */

/****v* bluebottle/rec_flow_field_phase
 * NAME
 *  rec_flow_field_phase
 * TYPE
 */
extern int rec_flow_field_phase;
/*
 * PURPOSE
 *  Recorder flow field output phase field boolean.
 ******
 */

/****v* bluebottle/rec_prec_flow_field_dt
 * NAME
 *  rec_prec_flow_field_dt
 * TYPE
 */

extern real rec_prec_flow_field_dt;
/*
 * PURPOSE
 *  Recorder turbulence flow field output timestep size.
 ******
 */

/****v* bluebottle/rec_prec_flow_field_ttime_out
 * NAME
 *  rec_prec_flow_field_ttime_out
 * TYPE
 */
extern real rec_prec_flow_field_ttime_out;
/*
 * PURPOSE
 *  Recorder turbulence flow field output time since last output.
 ******
 */

/****v* bluebottle/rec_prec_flow_field_vel
 * NAME
 *  rec_prec_flow_field_vel
 * TYPE
 */
extern int rec_prec_flow_field_vel;
/*
 * PURPOSE
 *  Recorder turbulence flow field output velocity field boolean.
 ******
 */

/****v* bluebottle/rec_prec_flow_field_p
 * NAME
 *  rec_prec_flow_field_p
 * TYPE
 */
extern int rec_prec_flow_field_p;
/*
 * PURPOSE
 *  Recorder turbulence flow field output pressure field boolean.
 ******
 */

/****v* bluebottle/rec_prec_flow_field_phase
 * NAME
 *  rec_prec_flow_field_phase
 * TYPE
 */
extern int rec_prec_flow_field_phase;
/*
 * PURPOSE
 *  Recorder turbulence flow field output phase field boolean.
 ******
 */

/****v* bluebottle/rec_paraview_dt
 * NAME
 *  rec_paraview_dt
 * TYPE
 */

/****v* bluebottle/rec_paraview_dt
 * NAME
 *  rec_paraview_dt
 * TYPE
 */
extern real rec_paraview_dt;
/*
 * PURPOSE
 *  Recorder ParaView output timestep size.
 ******
 */

/****v* bluebottle/rec_paraview_ttime_out
 * NAME
 *  rec_paraview_ttime_out
 * TYPE
 */
extern real rec_paraview_ttime_out;
/*
 * PURPOSE
 *  Recorder paraview output time since last output.
 ******
 */

/****v* bluebottle/rec_prec_dt
 * NAME
 *  rec_prec_dt
 * TYPE
 */
extern real rec_prec_dt;
/*
 * PURPOSE
 *  Recorder ParaView precursor output timestep size.
 ******
 */

/****v* bluebottle/rec_prec_ttime_out
 * NAME
 *  rec_prec_ttime_out
 * TYPE
 */
extern real rec_prec_ttime_out;
/*
 * PURPOSE
 *  Recorder precursor output time since last output.
 ******
 */

/****v* bluebottle/rec_restart_dt
 * NAME
 *  rec_restart_dt
 * TYPE
 */
extern real rec_restart_dt;
/*
 * PURPOSE
 *  Recorder restart output timestep size.
 ******
 */

/****v* bluebottle/rec_restart_stop
 * NAME
 *  rec_restart_stop
 * TYPE
 */
extern int rec_restart_stop;
/*
 * PURPOSE
 *  Recorder restart and stop output boolean: 0 = continue; 1 = stop
 ******
 */

/****v* bluebottle/rec_restart_ttime_out
 * NAME
 *  rec_restart_ttime_out
 * TYPE
 */
extern real rec_restart_ttime_out;
/*
 * PURPOSE
 *  Recorder restart output time since last output.
 ******
 */

/****v* bluebottle/rec_particle_dt
 * NAME
 *  rec_particle_dt
 * TYPE
 */
extern real rec_particle_dt;
/*
 * PURPOSE
 *  Recorder particle output timestep size.
 ******
 */

/****v* bluebottle/rec_particle_ttime_out
 * NAME
 *  rec_particle_ttime_out
 * TYPE
 */
extern real rec_particle_ttime_out;
/*
 * PURPOSE
 *  Recorder particle output time since last output.
 ******
 */

/****v* bluebottle/rec_particle_pos
 * NAME
 *  rec_particle_pos
 * TYPE
 */
extern int rec_particle_pos;
/*
 * PURPOSE
 *  Recorder particle output position.
 ******
 */

/****v* bluebottle/rec_particle_a
 * NAME
 *  rec_particle_a
 * TYPE
 */
extern int rec_particle_a;
/*
 * PURPOSE
 *  Recorder particle output radius.
 ******
 */

/****v* bluebottle/rec_particle_vel
 * NAME
 *  rec_particle_vel
 * TYPE
 */
extern int rec_particle_vel;
/*
 * PURPOSE
 *  Recorder particle output velocity.
 ******
 */

/****v* bluebottle/rec_particle_omega
 * NAME
 *  rec_particle_omega
 * TYPE
 */
extern int rec_particle_omega;
/*
 * PURPOSE
 *  Recorder particle output angular velocity.
 ******
 */

/****v* bluebottle/rec_particle_force
 * NAME
 *  rec_particle_force
 * TYPE
 */
extern int rec_particle_force;
/*
 * PURPOSE
 *  Recorder particle output hydrodynamic force.
 ******
 */

/****v* bluebottle/rec_particle_moment
 * NAME
 *  rec_particle_moment
 * TYPE
 */
extern int rec_particle_moment;
/*
 * PURPOSE
 *  Recorder particle output hydrodynamic moment.
 ******
 */

/****v* bluebottle/rho_f
 * NAME
 *  rho_f
 * TYPE
 */
extern real rho_f;
/*
 * PURPOSE
 *  The fluid density.
 ******
 */

/****v* bluebottle/mu
 * NAME
 *  mu
 * TYPE
 */
extern real mu;
/*
 * PURPOSE
 *  The fluid dynamic viscosity.
 ******
 */

/****v* bluebottle/nu
 * NAME
 *  nu
 * TYPE
 */
extern real nu;
/*
 * PURPOSE
 *  The fluid kinematic viscosity.
 ******
 */

/****s* bluebottle/g_struct
 * NAME
 *  g_struct
 * TYPE
 */
typedef struct g_struct {
  real x;
  real xm;
  real xa;
  real y;
  real ym;
  real ya;
  real z;
  real zm;
  real za;
} g_struct;
/*
 * PURPOSE
 *  Body force on flow.
 * MEMBERS
 *  * x -- x-direction
 *  * y -- y-direction
 *  * z -- z-direction
 ******
 */

/****v* bluebottle/g
 * NAME
 *  g
 * TYPE
 */
extern g_struct g;
/*
 * PURPOSE
 *  Create an instance of the struct g_struct to carry body forces.
 ******
 */

/****s* bluebottle/BC
 * NAME
 *  BC
 * TYPE
 */
typedef struct BC {
  int pW;
  real pWD;
  int pE;
  real pED;
  int pS;
  real pSD;
  int pN;
  real pND;
  int pB;
  real pBD;
  int pT;
  real pTD;
  int uW;
  real uWDm;
  real uWD;
  real uWDa;
  int uE;
  real uEDm;
  real uED;
  real uEDa;
  int uS;
  real uSDm;
  real uSD;
  real uSDa;
  int uN;
  real uNDm;
  real uND;
  real uNDa;
  int uB;
  real uBDm;
  real uBD;
  real uBDa;
  int uT;
  real uTDm;
  real uTD;
  real uTDa;
  int vW;
  real vWDm;
  real vWD;
  real vWDa;
  int vE;
  real vEDm;
  real vED;
  real vEDa;
  int vS;
  real vSDm;
  real vSD;
  real vSDa;
  int vN;
  real vNDm;
  real vND;
  real vNDa;
  int vB;
  real vBDm;
  real vBD;
  real vBDa;
  int vT;
  real vTDm;
  real vTD;
  real vTDa;
  int wW;
  real wWDm;
  real wWD;
  real wWDa;
  int wE;
  real wEDm;
  real wED;
  real wEDa;
  int wS;
  real wSDm;
  real wSD;
  real wSDa;
  int wN;
  real wNDm;
  real wND;
  real wNDa;
  int wB;
  real wBDm;
  real wBD;
  real wBDa;
  int wT;
  real wTDm;
  real wTD;
  real wTDa;
  real dsW;
  real dsE;
  real dsS;
  real dsN;
  real dsB;
  real dsT;
} BC;
/*
 * PURPOSE
 *  Carry the type of boundary condition on each side of the domain.  Possible
 *  types include:
 *  * PERIODIC
 *  * DIRICHLET
 *  * NEUMANN
 *  * PRECURSOR
 *  If the boundary type is DIRICHLET or PRECURSOR, the value of the field
 *  variable on the boundary must be defined.
 * MEMBERS
 *  * pW -- the boundary condition type
 *  * pWD -- the DIRICHLET boundary conditon value
 *  * pE -- the boundary condition type
 *  * pED -- the DIRICHLET boundary conditon value
 *  * pS -- the boundary condition type
 *  * pSD -- the DIRICHLET boundary conditon value
 *  * pN -- the boundary condition type
 *  * pND -- the DIRICHLET boundary conditon value
 *  * pB -- the boundary condition type
 *  * pBD -- the DIRICHLET boundary conditon value
 *  * pT -- the boundary condition type
 *  * pTD -- the DIRICHLET boundary conditon value
 *  * uW -- the boundary condition type
 *  * uWDm -- the maximum DIRICHLET boundary condition value
 *  * uWD -- the current DIRICHLET boundary conditon value
 *  * uWDa -- the DIRICHLET boundary condition value acceleration
 *  * uE -- the boundary condition type
 *  * uEDm -- the maximum DIRICHLET boundary condition value
 *  * uED -- the current DIRICHLET boundary conditon value
 *  * uEDa -- the DIRICHLET boundary condition value acceleration
 *  * uS -- the boundary condition type
 *  * uSDm -- the maximum DIRICHLET boundary condition value
 *  * uSD -- the current DIRICHLET boundary conditon value
 *  * uSDa -- the DIRICHLET boundary condition value acceleration
 *  * uN -- the boundary condition type
 *  * uND -- the maximum DIRICHLET boundary condition valuem
 *  * uND -- the current DIRICHLET boundary conditon value
 *  * uNDa -- the DIRICHLET boundary condition value acceleration
 *  * uB -- the boundary condition type
 *  * uBDm -- the maximum DIRICHLET boundary condition value
 *  * uBD -- the current DIRICHLET boundary conditon value
 *  * uBDa -- the DIRICHLET boundary condition value acceleration
 *  * uT -- the boundary condition type
 *  * uTDm -- the maximum DIRICHLET boundary condition value
 *  * uTD -- the current DIRICHLET boundary conditon value
 *  * uTDa -- the DIRICHLET boundary condition value acceleration
 *  * vW -- the boundary condition type
 *  * vWDm -- the maximum DIRICHLET boundary condition value
 *  * vWD -- the current DIRICHLET boundary conditon value
 *  * vWDa -- the DIRICHLET boundary condition value acceleration
 *  * vE -- the boundary condition type
 *  * vEDm -- the maximum DIRICHLET boundary condition value
 *  * vED -- the current DIRICHLET boundary conditon value
 *  * vEDa -- the DIRICHLET boundary condition value acceleration
 *  * vS -- the boundary condition type
 *  * vSDm -- the maximum DIRICHLET boundary condition value
 *  * vSD -- the current DIRICHLET boundary conditon value
 *  * vSDa -- the DIRICHLET boundary condition value acceleration
 *  * vN -- the boundary condition type
 *  * vNDm -- the maximum DIRICHLET boundary condition value
 *  * vND -- the current DIRICHLET boundary conditon value
 *  * vNDa -- the DIRICHLET boundary condition value acceleration
 *  * vB -- the boundary condition type
 *  * vBDm -- the maximum DIRICHLET boundary condition value
 *  * vBD -- the current DIRICHLET boundary conditon value
 *  * vBDa -- the DIRICHLET boundary condition value acceleration
 *  * vT -- the boundary condition type
 *  * vTDm -- the maximum DIRICHLET boundary condition value
 *  * vTD -- the current DIRICHLET boundary conditon value
 *  * vTDa -- the DIRICHLET boundary condition value acceleration
 *  * wW -- the boundary condition type
 *  * wWDm -- the maximum DIRICHLET boundary condition value
 *  * wWD -- the current DIRICHLET boundary conditon value
 *  * wWDa -- the DIRICHLET boundary condition value acceleration
 *  * wE -- the boundary condition type
 *  * wEDm -- the maximum DIRICHLET boundary condition value
 *  * wED -- the current DIRICHLET boundary conditon value
 *  * wEDa -- the DIRICHLET boundary condition value acceleration
 *  * wS -- the boundary condition type
 *  * wSDm -- the maximum DIRICHLET boundary condition value
 *  * wSD -- the current DIRICHLET boundary conditon value
 *  * wSDa -- the DIRICHLET boundary condition value acceleration
 *  * wN -- the boundary condition type
 *  * wNDm -- the maximum DIRICHLET boundary condition value
 *  * wND -- the current DIRICHLET boundary conditon value
 *  * wNDa -- the DIRICHLET boundary condition value acceleration
 *  * wB -- the boundary condition type
 *  * wBDm -- the maximum DIRICHLET boundary condition value
 *  * wBD -- the current DIRICHLET boundary conditon value
 *  * wBDa -- the DIRICHLET boundary condition value acceleration
 *  * wT -- the boundary condition type
 *  * wTDm -- the maximum DIRICHLET boundary condition value
 *  * wTD -- the current DIRICHLET boundary conditon value
 *  * wTDa -- the DIRICHLET boundary condition value acceleration
 *  * dsW -- the SCREEN boundary condition offset value
 *  * dsE -- the SCREEN boundary condition offset value
 *  * dsS -- the SCREEN boundary condition offset value
 *  * dsN -- the SCREEN boundary condition offset value
 *  * dsB -- the SCREEN boundary condition offset value
 *  * dsT -- the SCREEN boundary condition offset value
 ******
 */

/****v* bluebottle/bc
 * NAME
 *  bc
 * TYPE
 */
extern BC bc;
/*
 * PURPOSE
 *  Create an instance of the struct BC to carry boundary condition types.
 ******
 */

/****s* bluebottle/gradP_struct
 * NAME
 *  gradP_struct
 * TYPE
 */
typedef struct gradP_struct {
  real x;
  real xm;
  real xa;
  real y;
  real ym;
  real ya;
  real z;
  real zm;
  real za;
} gradP_struct;
/* 
 * PURPOSE
 *  Carry imposed pressure gradient values.
 ******
 */

/****v* bluebottle/gradP
 * NAME
 *  gradP
 * TYPE
 */
extern gradP_struct gradP;
/*
 * PURPOSE
 *  Create an instance of the struct gradP_struct to carry imposed pressure
 *  gradient values.
 ******
 */

/****v* bluebottle/turbA
 * NAME
 *  turbA
 * TYPE
 */
extern real turbA;
/*
 * PURPOSE
 *  Turbulence precursor linear forcing magnitude (see Lundgren 2003, Rosales
 *  and Meneveau 2005, Carroll and Blanquart 2013).
 ******
 */

/****v* bluebottle/turbl
 * NAME
 *  turbl
 * TYPE
 */
extern real turbl;
/*
 * PURPOSE
 *  Turbulence precursor linear forcing integral scale (see Lundgren 2003,
 *  Rosales and Meneveau 2005, Carroll and Blanquart 2013).
 ******
 */

/****f* bluebottle/cuda_dom_malloc()
 * NAME
 *  cuda_dom_malloc()
 * USAGE
 */
void cuda_dom_malloc(void);
/*
 * FUNCTION
 *  Allocate device memory reference pointers on host and device memory on
 *  device for the flow domain.
 ******
 */

/****f* bluebottle/cuda_part_malloc()
 * NAME
 *  cuda_part_malloc()
 * USAGE
 */
void cuda_part_malloc(void);
/*
 * FUNCTION
 *  Allocate device memory reference pointers on host and device memory on
 *  device for the particles.
 ******
 */

/****f* bluebottle/cuda_dom_push()
 * NAME
 *  cuda_dom_push()
 * USAGE
 */
void cuda_dom_push(void);
/*
 * FUNCTION
 *  Copy p, u, v, and w from host to device.
 ******
 */

/****f* bluebottle/cuda_dom_turb_planes_push()
 * NAME
 *  cuda_dom_turb_planes_push()
 * USAGE
 */
void cuda_dom_turb_planes_push(int *bc_configs);
/*
 * FUNCTION
 *  Copy u, v, and w turbulence precursor planes from host to device.
 ******
 */

/****f* bluebottle/cuda_part_push()
 * NAME
 *  cuda_part_push()
 * USAGE
 */
void cuda_part_push(void);
/*
 * FUNCTION
 *  Copy particle data from host to device.
 ******
 */

/****f* bluebottle/cuda_dom_pull()
 * NAME
 *  cuda_dom_pull()
 * USAGE
 */
void cuda_dom_pull(void);
/*
 * FUNCTION
 *  Copy p, u, v, and w turbulence precursor planes from device to host.
 ******
 */

/****f* bluebottle/cuda_dom_turb_planes_pull()
 * NAME
 *  cuda_dom_turb_planes_pull()
 * USAGE
 */
void cuda_dom_turb_planes_pull(int *bc_configs);
/*
 * FUNCTION
 *  Copy u, v, and w from device to host.
 ******
 */

/****f* bluebottle/cuda_part_pull()
 * NAME
 *  cuda_part_pull()
 * USAGE
 */
void cuda_part_pull(void);
/*
 * FUNCTION
 *  Copy particle data from device to host.
 ******
 */

/****f* bluebottle/cuda_dom_free()
 * NAME
 *  cuda_dom_free()
 * USAGE
 */
void cuda_dom_free(void);
/*
 * FUNCTION
 *  Free device memory for the domain on device and device memory reference
 *  pointers on host.
 ******
 */

/****f* bluebottle/cuda_part_free()
 * NAME
 *  cuda_part_free()
 * USAGE
 */
void cuda_part_free(void);
/*
 * FUNCTION
 *  Free device memory for the particles on device and device memory reference
 *  pointers on host.
 ******
 */

/****f* bluebottle/cuda_dom_BC()
 * NAME
 *  cuda_dom_BC()
 * USAGE
 */
void cuda_dom_BC(void);
/*
 * FUNCTION
 *  Enforce boundary conditions in velocity and pressure fields.  *NOTE:*
 *  cuda_BC() implements PERIODIC boundary conditions for single-GPU only.
 *  Dirichlet and Neumann boundary conditions are supported on multi-GPU
 *  domain decompositions.
 ******
 */

/****f* bluebottle/cuda_dom_BC_star()
 * NAME
 *  cuda_dom_BC_star()
 * USAGE
 */
void cuda_dom_BC_star(void);
/*
 * FUNCTION
 *  Enforce boundary conditions in intermediate velocity fields.  *See note
 *  in cuda_BC().
 ******
 */

/****f* bluebottle/cuda_part_BC_star()
 * NAME
 *  cuda_part_BC_star()
 * USAGE
 */
void cuda_part_BC_star(void);
/*
 * FUNCTION
 *  Enforce boundary conditions in intermediate velocity fields.  *See note
 *  in cuda_BC().
 ******
 */

/****f* bluebottle/cuda_part_BC_p()
 * NAME
 *  cuda_part_BC_p()
 * USAGE
 */
void cuda_part_BC_p(int dev);
/*
 * FUNCTION
 *  Enforce boundary conditions in pressure-Poisson problem.
 * ARGUMENTS
 *  * dev -- the device number
 ******
 */

/****f* bluebottle/cuda_part_p_fill()
 * NAME
 *  cuda_part_p_fill()
 * USAGE
 */
void cuda_part_p_fill();
/*
 * FUNCTION
 *  Enforce boundary conditions in pressure-Poisson problem.
 * ARGUMENTS
 *  * dev -- the device number
 ******
 */

/****f* bluebottle/cuda_solvability()
 * NAME
 *  cuda_solvability()
 * USAGE
 */
void cuda_solvability(void);
/*
 * FUNCTION
 *  Enforce the solvability condition on the Poisson problem.
 ******
 */

/****f* bluebottle/cuda_dom_BC_p()
 * NAME
 *  cuda_dom_BC_p()
 * USAGE
 */
void cuda_dom_BC_p(void);
void cuda_dom_BC_phi(void);
/*
 * FUNCTION
 *  Enforce zero normal gradient boundary conditions on projected pressure.
 ******
 */

/****f* bluebottle/cuda_U_star_2()
 * NAME
 *  cuda_U_star_2()
 * USAGE
 */
void cuda_U_star_2(void);
/*
 * FUNCTION
 *  Compute u_star, v_star, w_star to second-order time accuracy.
 ******
 */

/****f* bluebottle/cuda_project()
 * NAME
 *  cuda_project()
 * USAGE
 */
void cuda_project(void);
/*
 * FUNCTION
 *  Project the intermediate velocity U* onto a divergence-free space
 *  via the projected pressure.
 ******
 */

/****f* bluebottle/cuda_U_star_test_exp()
 * NAME
 *  cuda_U_star_test_exp()
 * USAGE
 */
void cuda_U_star_test_exp(void);
/*
 * FUNCTION
 *  cuda_U_star_1 validation testbed using u = exp(x), v = exp(y), w = exp(z).
 ******
 */

/****f* bluebottle/cuda_U_star_test_cos()
 * NAME
 *  cuda_U_star_test_cos()
 * USAGE
 */
void cuda_U_star_test_cos(void);
/*
 * FUNCTION
 *  cuda_U_star_1 validation testbed using u = cos(x), v = cos(y), w = cos(z).
 ******
 */

/****f* bluebottle/cuda_ustar_helmholtz
 * NAME
 *  cuda_ustar_helmholtz()
 * USAGE
 */
void cuda_ustar_helmholtz(int rank);
/*
 * FUNCTION
 *  Entry point for the CUSP BiCGSTAB and required computations for the
 *  u_star Helmholtz problem.
 * ARGUMENTS
 *  * rank -- The MPI rank of the process to differentiate between flow and
 *    precursor domains.
 ******
 */

/****f* bluebottle/cuda_vstar_helmholtz
 * NAME
 *  cuda_vstar_helmholtz()
 * USAGE
 */
void cuda_vstar_helmholtz(int rank);
/*
 * FUNCTION
 *  Entry point for the CUSP BiCGSTAB and required computations for the
 *  v_star Helmholtz problem.
 * ARGUMENTS
 *  * rank -- The MPI rank of the process to differentiate between flow and
 *    precursor domains.
 ******
 */

/****f* bluebottle/cuda_wstar_helmholtz
 * NAME
 *  cuda_wstar_helmholtz()
 * USAGE
 */
void cuda_wstar_helmholtz(int rank);
/*
 * FUNCTION
 *  Entry point for the CUSP BiCGSTAB and required computations for the
 *  w_star Helmholtz problem.
 * ARGUMENTS
 *  * rank -- The MPI rank of the process to differentiate between flow and
 *    precursor domains.
 ******
 */

/****f* bluebottle/cuda_PP_bicgstab()
 * NAME
 *  cuda_PP_bicgstab()
 * USAGE
 */
void cuda_PP_bicgstab(int rank);
/*
 * FUNCTION
 *  Entry point for the CUSP BiCGSTAB and required computations for the
 *  pressure-Poisson problem.
 * ARGUMENTS
 *  * rank -- The MPI rank of the process to differentiate between flow and
 *    precursor domains.
 ******
 */

/****f* bluebottle/cuda_PP_bicgstab_jacobi()
 * NAME
 *  cuda_PP_bicgstab_jacobi()
 * USAGE
 */
void cuda_PP_bicgstab_jacobi(void);
/*
 * FUNCTION
 *  Entry point for the Jacobi iterative solver and required computations for
 *  the pressure-Poisson problem.
 ******
 */

/****f* bluebottle/cuda_div_U()
 * NAME
 *  cuda_div_U(void);
 * USAGE
 */
void cuda_div_U(void);
/*
 * FUNCTION
 *  Calculate the velocity divergence field.
 ******
 */

/****f* bluebottle/cuda_build_cages()
 * NAME
 *  cuda_build_cages()
 * USAGE
 */
void cuda_build_cages(void);
/*
 * FUNCTION
 *  Build the particle cages on the devices and generate phase array.
 ******
 */

/****f* bluebottle/cuda_part_BC()
 * NAME
 *  cuda_part_BC()
 * USAGE
 */
void cuda_part_BC(void);
/*
 * FUNCTION
 *  Apply the particle velocity boundary conditions to the flow domain.
 ******
 */

/****f* bluebottle/cuda_project_test()
 * NAME
 *  cuda_project_test()
 * USAGE
 */
void cuda_project_test(void);
/*
 * FUNCTION
 *  A validation testbed for cuda_project()using u = exp(x), v = exp(y),
 *  w = exp(z), p = exp(x) + exp(y) + exp(z).
 ******
 */

/****f* bluebottle/cuda_BC_test()
 * NAME
 *  cuda_BC_test()
 * USAGE
 */
void cuda_BC_test(void);
/*
 * FUNCTION
 *  Validation testbed for cuda_BC.  Domain must be -1 <= x <= 1,
 *  -1 <= y <= 1, -1 <= z <= 1.  In order to test cuda_BC, ensure that there
 *  are no pressure gradients applied to the system.
 ******
 */

/****f* bluebottle/cuda_quad_interp_test()
 * NAME
 *  cuda_quad_interp_test()
 * USAGE
 */
void cuda_quad_interp_test(void);
/*
 * FUNCTION
 *  Validation testbed for Lebedev quadrature interpolation function.
 ******
 */

/****f* bluebottle/cuda_lamb_test()
 * NAME
 *  cuda_lamb_test()
 * USAGE
 */
void cuda_lamb_test(void);
/*
 * FUNCTION
 *  Validation testbed for Lamb's coefficient calculation algorithm.
 ******
 */

/****f* bluebottle/cuda_find_dt()
 * NAME
 *  cuda_find_dt()
 * USAGE
 */
real cuda_find_dt(void);
/*
 * FUNCTION
 *  Determine the timestep to use based on the current flow fields and
 *  the Courant-Friedrichs-Lewy (CFL) condition.
 * RESULT
 *  The timestep to be used for the next step forward.
 ******
 */

/****f* bluebottle/cuda_store_u()
 * NAME
 *  cuda_store_u()
 * USAGE
 */
void cuda_store_u(void);
/*
 * FUNCTION
 *  Store the previous u, v, w components for use in the next timestep.
 ******
 */

/****f* bluebottle/cuda_update_p()
 * NAME
 *  cuda_update_p()
 * USAGE
 */
void cuda_update_p(void);
/*
 * FUNCTION
 *  Update the pressure.
 ******
 */

/****f* bluebottle/cuda_store_coeffs()
 * NAME
 *  cuda_store_coeffs()
 * USAGE
 */
void cuda_store_coeffs(void);
/*
 * FUNCTION
 *  Store the previous Lamb's coefficietns.
 ******
 */

/****f* bluebottle/cuda_compute_forcing()
 * NAME
 *  cuda_compute_forcing()
 * USAGE
 */
void cuda_compute_forcing(real *pid_int, real *pid_back, real Kp, real Ki,
  real Kd);
/*
 * FUNCTION
 *  Set up the forcing array for this time step. It also partially incorporates
 *  a PID controller to push the net force on particles to zero by adjusting
 *  the pressure gradient. It currently only works in the z-direction.
 * MEMBERS
 *  * pid_int - the integral up to this point in time
 *  * pid_back - the previous value for the net force, to be used by the
 *    derivative term
 *  * Kp - the proportional gain
 *  * Ki - the integral gain
 *  * Kd - the derivative gain
 ******
 */

/****f* bluebottle/cuda_compute_turb_forcing()
 * NAME
 *  cuda_compute_turb_forcing()
 * USAGE
 */
void cuda_compute_turb_forcing(void);
/*
 * FUNCTION
 *  Add in the turbulence precursor linear turbulence forcing.
 ******
 */

/****f* bluebottle/cuda_move_parts_sub()
 * NAME
 *  cuda_move_parts_sub()
 * USAGE
 */
void cuda_move_parts_sub(void);
/*
 * FUNCTION
 *  Calculate new particle velocities and positions for sub-timestep implicit
 *  iteration, but do not update the new positions.
 ******
 */

/****f* bluebottle/cuda_move_parts()
 * NAME
 *  cuda_move_parts()
 * USAGE
 */
void cuda_move_parts(void);
/*
 * FUNCTION
 *  Calculate new particle velocities and positions.
 ******
 */

/****f* bluebottle/cuda_parts_internal()
 * NAME
 *  cuda_parts_internal()
 * USAGE
 */
void cuda_parts_internal(void);
/*
 * FUNCTION
 *  Apply particle solid-body motion to internal velocity nodes.
 ******
 */

/****f* bluebottle/cuda_yank_turb_planes()
 * NAME
 *  cuda_yank_turb_planes()
 * USAGE
 */
void cuda_yank_turb_planes(int *bc_flow_configs, real *bc_plane_pos, real *vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor flow field planes to be used
 *  for inflow boundary conditions.
 * ARGUMENTS
 *  * bc_plane_pos -- a 9-element array containing the (x,y,z) positions of the
 *    centers of each turbulent interpolation plane (WE, SN, BT) in this order:
 *    WEx, SNx, BTx, WEy, SNy, BTy, WEz, SNz, BTz
 ******
 */

/****f* bluebottle/cuda_collisions()
 * NAME
 *  cuda_collisions()
 * USAGE
 */
void cuda_collisions(void);
/*
 * FUNCTION
 *  Calculate extra collision forces.
 ******
 */

/****f* bluebottle/seeder_read_input()
 * NAME
 *  seeder_read_input()
 * USAGE
 */
void seeder_read_input(int Nx, int Ny, int Nz, double ddz, double bias, 
  int nperturb);
/*
  * FUNCTION
  *   Read parts.config for seeder initialization
  * ARGUMENTS
  *  * Nx -- number of particles in the x direction
  *  * Ny -- number of particles in the y direction
  *  * Nz -- number of particles in the z direction
  *  * ddz -- distance b/t top of one layer and middle of next
  *  * bias -- perturbation amount
  *  * nperturb -- number of times to perturb
  ******
  */

/****f* bluebottle/seeder()
 * NAME
 *  seeder()
 * USAGE
 */
void seeder(int nparts, real loa, real a, real aFx, real aFy, real aFz, 
  real aLx, real aLy, real aLz, real rho, real E, real sigma, real e_dry,
  int o, real rs_r, real spring_k, real spring_x, real spring_y,
  real spring_z, real spring_l, int t, int r);
/*
 * FUNCTION
 *  Randomly seed nparts particles in the domain. To use this function, run
 *  'bluebottle -s', using the arguments as defined below. A new
 *  file called part_seeder.config will be created and the main bluebottle
 *  simulation code will not be run. To run a simulation using this input file,
 *  change its name to part.config and run bluebottle normally.
 * ARGUMENTS
 *  * nparts -- the number of particles
 *  * loa -- particle interaction compact support length
 *  * a -- the radius of the particles
 *  * aFx -- particle x forcing
 *  * aFy -- particle y forcing
 *  * aFz -- particle z forcing
 *  * aLx -- particle x angular forcing
 *  * aLy -- particle y angular forcing
 *  * aLz -- particle z angular forcing
 *  * d -- the density of the particles (rho)
 *  * E -- the particle Young's modulus
 *  * s -- the particle Poisson ratio (-1 < s <= 0.5)
 *  * e_dry -- dry coefficient of restitution
 *  * o -- the order of truncation of the Lamb's solution
 *  * rs_r -- cage extent ratio
 *  * spring_k -- particle spring constant
 *  * spring_x -- particle spring x location
 *  * spring_y -- particle spring y location
 *  * spring_z -- particle spring z location
 *  * spring_l -- particle spring length
 *  * t -- one if translating, zero if not
 *  * r -- one if rotating, zero if not
 ******
 */

 /****f* bluebottle/seeder_array()
 * NAME
 *  seeder_array()
 * USAGE
 */
void seeder_array(int Nx, int Ny, int Nz, real loa, real a, real aFx, real aFy, 
  real aFz, real aLx, real aLy, real aLz, real rho, real E, real sigma, 
  real e_dry, int o, real rs_r, real spring_k, real spring_x, 
  real spring_y, real spring_z, real spring_l, int t, int r);
/* FUNCTION
 *  Seed Nx*Ny*Nz particles in the domain as a regular array shape. To use this 
 *  function, run 'bluebottle -s. A new file called part_seeder_array.config 
 *  will be created and the main bluebottle simulation code will not be run. To 
 *  run a simulation using this input file, change its name to part.config and 
 *  run bluebottle normally.
 * ARGUMENTS
 *  * Nx -- particle number in x direction
 *  * Ny -- particle number in y direction
 *  * Nz -- particle number is z direction
 *  * loa -- particle interaction compact support length
 *  * a -- the radius of the particles
 *  * aFx -- particle x forcing
 *  * aFy -- particle y forcing
 *  * aFz -- particle z forcing
 *  * aLx -- particle x angular forcing
 *  * aLy -- particle y angular forcing
 *  * aLz -- particle z angular forcing
 *  * d -- the density of the particles (rho)
 *  * E -- the particle Young's modulus
 *  * s -- the particle Poisson ratio (-1 < s <= 0.5)
 *  * e_dry -- dry coefficient of restitution
 *  * o -- the order of truncation of the Lamb's solution
 *  * rs_r -- cage extent ratio
 *  * spring_k -- particle spring constant
 *  * spring_x -- particle spring x location
 *  * spring_y -- particle spring y location
 *  * spring_z -- particle spring z location
 *  * spring_l -- particle spring length
 *  * t -- one if translating, zero if not
 *  * r -- one if rotating, zero if not
******
*/

/****f* bluebottle/seeder_hex()
 * NAME
 *  seeder_hex()
 * USAGE
 */
void seeder_hex(int Nx, int Ny, int Nz, double ddz, real loa, real a, real aFx, 
  real aFy, real aFz, real aLx, real aLy, real aLz, real rho, real E, 
  real sigma, real e_dry, int o, real rs_r, real spring_k, 
  real spring_x, real spring_y, real spring_z, real spring_l, int t, int r);
/* FUNCTION
 *  Seed Nx*Ny*Nz particles in the domain in a hex shape. To use this function, 
 *  run 'bluebottle -s. A new file called part_seeder_hex.config will be created
 *  and the main bluebottle simulation code will not be run. To run a simulation
 *  using this input file, change its name to part.config and run bluebottle 
 *  normally.
 * ARGUMENTS
 * 	Nx -- particle number in x direction
 * 	Ny -- particle number in y direction
 *  Nz -- particle number is z direction
 *	ddz -- distance from top of one layer to center of next layer
 *  * loa -- particle interaction compact support length
 *  * a -- the radius of the particles
 *  * aFx -- particle x forcing
 *  * aFy -- particle y forcing
 *  * aFz -- particle z forcing
 *  * aLx -- particle x angular forcing
 *  * aLy -- particle y angular forcing
 *  * aLz -- particle z angular forcing
 *  * d -- the density of the particles (rho)
 *  * E -- the particle Young's modulus
 *  * s -- the particle Poisson ratio (-1 < s <= 0.5)
 *  * e_dry -- dry coefficient of restitution
 *  * o -- the order of truncation of the Lamb's solution
 *  * rs_r -- cage extent ratio
 *  * spring_k -- particle spring constant
 *  * spring_x -- particle spring x location
 *  * spring_y -- particle spring y location
 *  * spring_z -- particle spring z location
 *  * spring_l -- particle spring length
 *  * t -- one if translating, zero if not
 *  * r -- one if rotating, zero if not
 ******
*/

/****f* bluebottle/seeder_high_vol_random()
 * NAME
 *  seeder_high_vol_random()
 * USAGE
 */
void seeder_high_vol_random(int Nx, int Ny, int Nz, double bias, int nperturb, 
  real loa, real a, real aFx, real aFy, real aFz, real aLx, real aLy, real aLz, 
  real rho, real E, real sigma, real e_dry, int o, real rs, 
  real spring_k, real spring_x, real spring_y, real spring_z, real spring_l, 
  int t, int r);
/* FUNCTION
 *  Seed Nx*Ny*Nz particles in the domain randomly by perturbing the regular 
 *  array nperturb times. To use this function, run 'bluebottle -s'. A new file 
 *  called part_seeder_hex.config will be created and the main bluebottle
 *  simulation code will not be run. To run a simulation using this input file,
 *  change its name to part.config and run bluebottle normally.
 * ARGUMENTS
 * 	Nx - particle number in x direction
 * 	Ny - particle number in y direction
 *  Nz - particle number is z direction
 *	bias -- magnitude of pertubation
 *	nperturb -- number of pertubation times
 *  * loa -- particle interaction compact support length
 *  * a -- the radius of the particles
 *  * aFx -- particle x forcing
 *  * aFy -- particle y forcing
 *  * aFz -- particle z forcing
 *  * aLx -- particle x angular forcing
 *  * aLy -- particle y angular forcing
 *  * aLz -- particle z angular forcing
 *  * d -- the density of the particles (rho)
 *  * E -- the particle Young's modulus
 *  * s -- the particle Poisson ratio (-1 < s <= 0.5)
 *  * e_dry -- dry coefficient of restitution
 *  * o -- the order of truncation of the Lamb's solution
 *  * rs_r -- cage extent ratio
 *  * spring_k -- particle spring constant
 *  * spring_x -- particle spring x location
 *  * spring_y -- particle spring y location
 *  * spring_z -- particle spring z location
 *  * spring_l -- particle spring length
 *  * t -- one if translating, zero if not
 *  * r -- one if rotating, zero if not
 ******
*/

/****f* bluebottle/out_restart()
 * NAME
 *  out_restart()
 * USAGE
 */
void out_restart(void);
/*
 * FUNCTION
 *  Write the data required for restarting the simulation to file.
 ******
 */

/****f* bluebottle/in_restart()
 * NAME
 *  in_restart()
 * USAGE
 */
void in_restart(void);
/*
 * FUNCTION
 *  Read the data required for restarting the simulation to file.
 ******
 */

/****f* bluebottle/out_restart_turb()
 * NAME
 *  out_restart_turb()
 * USAGE
 */
void out_restart_turb(void);
/*
 * FUNCTION
 *  Write the data required for restarting the precursor simulation to file.
 ******
 */

/****f* bluebottle/in_restart_turb()
 * NAME
 *  in_restart_turb()
 * USAGE
 */
void in_restart_turb(void);
/*
 * FUNCTION
 *  Read the data required for restarting the precursor simulation to file.
 ******
 */

/****v* bluebottle/cpumem
 * NAME
 *  cpumem
 * TYPE
 */
extern long int cpumem;
/*
 * PURPOSE
 *  Counts total cpu usage 
 ******
 */

/****v* bluebottle/gpumem
 * NAME
 *  gpumem
 * TYPE
 */
extern long int gpumem;
/*
 * PURPOSE
 *  Counts total gpu usage 
 ******
 */

/****v* bluebottle/init_cond
 * NAME
 *  init_cond
 * TYPE
 */
extern int init_cond;
/*
 * PURPOSE
 *  Carries the initial condition type. For now, the only options are QUIESCENT
 *  and SHEAR.
 ******
 */

/****v* bluebottle/bc_flow_configs
 * NAME
 *  bc_flow_configs
 * TYPE
 */
extern int bc_flow_configs[18];
/*
 * PURPOSE
 *  Pertinent boundary condition configuration information from flow domain
 *  used for turbulent inflow communication. It contains the boundary condition
 *  configuration in this order:
 *  [uW uE uS uN uB uT vW vE vS vN vB vT wW wE wS wN wB wT]
 ******
 */

/****v* bluebottle/bc_flow_vels
 * NAME
 *  bc_flow_vels
 * TYPE
 */
extern real bc_flow_vels[18];
/*
 * PURPOSE
 *  Inflow domain boundary velocity information from flow domain
 *  used for turbulent inflow communication. It contains the boundary velocities
 *  in this order:
 *  [uWD uED uSD uND uBD uTD vWD vED vSD vND vBD vTD wWD wED wSD wND wBD wTD]
 ******
 */

/****v* bluebottle/bc_plane_pos
 * NAME
 *  bc_plane_pos
 * TYPE
 */
extern real bc_plane_pos[9];
/*
 * PURPOSE
 *  Inflow domain plane position in precursor domain used for turbulent
 *  inflow communication. It contains the boundary positions in this order:
 *  [WEx SNx BTx WEy SNy BTy WEz SNz BTz]
 ******
 */

/****v* bluebottle/pid_int
 * NAME
 *  pid_int
 * TYPE
 */
extern real pid_int;
/*
 * PURPOSE
 *  Store the integral of the PID controller target.
 ******
 */

/****v* bluebottle/pid_back
 * NAME
 *  pid_back
 * TYPE
 */
extern real pid_back;
/*
 * PURPOSE
 *  Store the previous value of the PID controller target for derivative
 *  term calculation.
 ******
 */

/****v* bluebottle/Kp
 * NAME
 *  Kp
 * TYPE
 */
extern real Kp;
/*
 * PURPOSE
 *  PID controller proportional gain.
 ******
 */

/****v* bluebottle/Ki
 * NAME
 *  Ki
 * TYPE
 */
extern real Ki;
/*
 * PURPOSE
 *  PID controller integral gain.
 ******
 */

/****v* bluebottle/Kd
 * NAME
 *  Kd
 * TYPE
 */
extern real Kd;
/*
 * PURPOSE
 *  PID controller derivative gain.
 ******
 */

/****f* bluebottle/cuda_compute_energy()
 * NAME
 *  cuda_compute_energy()
 * TYPE
 */
real cuda_compute_energy(void);
/*
 * PURPOSE
 *  Compute the kinetic energy k = 1/2 * mean(u*u).
 * RESULT
 *  The kinetic energy in the domain.
 ******
 */

/****f* bluebottle/cuda_colocate_Gfx()
 * NAME
 *  cuda_colocate_Gfx()
 * TYPE
 */
void cuda_colocate_Gfx(real *_u, real *_u_co, dom_struct *_dom);
/*
 * PURPOSE
 *  Interpolate a Gfx grid (e.g. u-velocity) to a Gcc grid. This routine
 *  interpolates from a Gfx grid with ghost cell buffers to a Gcc grid with
 *  no ghost cell buffers.
 * ARGUMENTS 
 *  * _u -- device Gfx field with ghost cell buffers
 *  * _u_co -- device Gcc field with no ghost cell buffers
 *  * _dom -- device domain data
 ******
 */

/****f* bluebottle/cuda_colocate_Gfy()
 * NAME
 *  cuda_colocate_Gfy()
 * TYPE
 */
void cuda_colocate_Gfy(real *_v, real *_v_co, dom_struct *_dom);
/*
 * PURPOSE
 *  Interpolate a Gfy grid (e.g. v-velocity) to a Gcc grid. This routine
 *  interpolates from a Gfy grid with ghost cell buffers to a Gcc grid with
 *  no ghost cell buffers.
 * ARGUMENTS 
 *  * _v -- device Gfy field with ghost cell buffers
 *  * _v_co -- device Gcc field with no ghost cell buffers
 *  * _dom -- device domain data
 ******
 */

/****f* bluebottle/cuda_colocate_Gfz()
 * NAME
 *  cuda_colocate_Gfz()
 * TYPE
 */
void cuda_colocate_Gfz(real *_w, real *_w_co, dom_struct *_dom);
/*
 * PURPOSE
 *  Interpolate a Gfz grid (e.g. w-velocity) to a Gcc grid. This routine
 *  interpolates from a Gfz grid with ghost cell buffers to a Gcc grid with
 *  no ghost cell buffers.
 * ARGUMENTS 
 *  * _w -- device Gfz field with ghost cell buffers
 *  * _w_co -- device Gcc field with no ghost cell buffers
 *  * _dom -- device domain data
 ******
 */

#endif
