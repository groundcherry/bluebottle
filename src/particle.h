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

/****h* Bluebottle/particle
 * NAME
 *  particle
 * FUNCTION
 *  Bluebottle particle functions.
 ******
 */

#ifndef _PARTICLE_H
#define _PARTICLE_H

#include "bluebottle.h"

/****d* particle/NNODES
 * NAME
 *  NNODES
 * TYPE
 */
#define NNODES 26
/*
 * PURPOSE
 *  Define the number of nodes used for the Lebedev quadrature scheme.
 ******
 */

/****v* particle/phase
 * NAME
 *  phase
 * TYPE
 */
extern int *phase;
/*
 * PURPOSE
 *  The phase of a discretized cell (Gcc-type grid).  If cell C is not inside
 *  a particle, then phase[C] = -1.  Otherwise, phase[C] is equal to the index 
 *  assigned to the particle in which the cell resides.
 ******
 */

/****v* particle/_phase
 * NAME
 *  _phase
 * TYPE
 */
extern int **_phase;
/*
 * PURPOSE
 *  CUDA device analog for phase.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* particle/phase_shell
 * NAME
 *  phase_shell
 * TYPE
 */
extern int *phase_shell;
/*
 * PURPOSE
 *  The outermost shell of phase, denoting the positions for Dirichlet pressure
 *  boundary conditions
 ******
 */

/****v* particle/_phase_shell
 * NAME
 *  _phase_shell
 * TYPE
 */
extern int **_phase_shell;
/*
 * PURPOSE
 *  CUDA device analog for phase_shell.  It contains pointers to arrays
 *  containing the subdomain fields on which each device operates.
 ******
 */

/****v* particle/flag_u
 * NAME
 *  flag_u
 * TYPE
 */
extern int *flag_u;
/* 
 * PURPOSE
 *  Flag x-direction components of velocity field that are set as boundaries.
 ******
 */

/****v* particle/_flag_u
 * NAME
 *  _flag_u
 * TYPE
 */
extern int **_flag_u;
/*
 * PURPOSE
 *  CUDA device analog for flag_u.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* particle/flag_v
 * NAME
 *  flag_v
 * TYPE
 */
extern int *flag_v;
/* 
 * PURPOSE
 *  Flag y-direction components of velocity field that are set as boundaries.
 ******
 */

/****v* particle/_flag_v
 * NAME
 *  _flag_v
 * TYPE
 */
extern int **_flag_v;
/*
 * PURPOSE
 *  CUDA device analog for flag_v.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* particle/flag_w
 * NAME
 *  flag_w
 * TYPE
 */
extern int *flag_w;
/* 
 * PURPOSE
 *  Flag z-direction components of velocity field that are set as boundaries.
 ******
 */

/****v* particle/_flag_w
 * NAME
 *  _flag_w
 * TYPE
 */
extern int **_flag_w;
/*
 * PURPOSE
 *  CUDA device analog for flag_w.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****s* particle/cage_struct
 * NAME
 *  cage_struct
 * TYPE
 */
typedef struct cage_struct {
  int cx;
  int cy;
  int cz;
  int is;
  int ibs;
  int ie;
  int ibe;
  int in;
  int js;
  int jbs;
  int je;
  int jbe;
  int jn;
  int ks;
  int kbs;
  int ke;
  int kbe;
  int kn;
} cage_struct;
/*
 * PURPOSE
 *  Carry domain location index information for the cage of cells containing
 *  a particle.
 * MEMBERS
 *  * cx -- the index location of the cell containing the particle center
 *  * cy -- the index location of the cell containing the particle center
 *  * cz -- the index location of the cell containing the particle center
 *  * is -- i start index
 *  * ibs -- i start break index (for periodic boundary straddling)
 *  * ie -- i end index
 *  * ibe -- i end break index (for periodic boundary straddling)
 *  * in -- number of indices in the x-direction
 *  * js -- j start index
 *  * jbs -- j start break index (for periodic boundary straddling)
 *  * je -- j end index
 *  * jbe -- j end break index (for periodic boundary straddling)
 *  * jn -- number of indices in the y-direction
 *  * ks -- k start index
 *  * kbs -- k start break index (for periodic boundary straddling)
 *  * ke -- k end index
 *  * kbe -- k end break index (for periodic boundary straddling)
 *  * kn -- number of indices in the z-direction
 ******
 */

/****s* particle/part_struct
 * NAME
 *  part_struct
 * TYPE
 */
typedef struct part_struct {
  real r;
  real x;
  real y;
  real z;
  real x0;
  real y0;
  real z0;
  real u;
  real v;
  real w;
  real u0;
  real v0;
  real w0;
  real udot;
  real vdot;
  real wdot;
  real udot0;
  real vdot0;
  real wdot0;
  real axx;
  real axy;
  real axz;
  real ayx;
  real ayy;
  real ayz;
  real azx;
  real azy;
  real azz;
  real ox;
  real oy;
  real oz;
  real ox0;
  real oy0;
  real oz0;
  real oxdot;
  real oydot;
  real ozdot;
  real oxdot0;
  real oydot0;
  real ozdot0;
  real Fx;
  real Fy;
  real Fz;
  real Lx;
  real Ly;
  real Lz;
  real aFx;
  real aFy;
  real aFz;
  real aLx;
  real aLy;
  real aLz;
  real kFx;
  real kFy;
  real kFz;
  real iFx;
  real iFy;
  real iFz;
  real iLx;
  real iLy;
  real iLz;
  int nodes[NNODES];
  real rho;
  real E;
  real sigma;
  cage_struct cage;
  int order;
  real rs;
  int ncoeff;
  real spring_k;
  real spring_x;
  real spring_y;
  real spring_z;
  real spring_l;
  int translating;
  int rotating;
  int bin;
  real St;
  real e_dry;
  real l_rough;
} part_struct;
/*
 * PURPOSE
 *  Carry physical information regarding a particle.
 * MEMBERS
 *  * r -- the particle radius size
 *  * x -- the particle location component
 *  * y -- the particle location component
 *  * z -- the particle location component
 *  * x0 -- the particle location component (previous timestep)
 *  * y0 -- the particle location component (previous timestep)
 *  * z0 -- the particle location component (previous timestep)
 *  * u -- linear velocity in x-direction
 *  * v -- linear velocity in y-direction
 *  * w -- linear velocity in z-direction
 *  * u0 -- initial linear velocity in x-direction
 *  * v0 -- initial linear velocity in y-direction
 *  * w0 -- initial linear velocity in z-direction
 *  * udot -- linear acceleration in x-direction
 *  * vdot -- linear acceleration in y-direction
 *  * wdot -- linear acceleration in z-direction
 *  * axx -- x-component of basis vector initially coincident with x-axis
 *  * axy -- y-component of basis vector initially coincident with x-axis
 *  * axz -- z-component of basis vector initially coincident with x-axis
 *  * ayx -- x-component of basis vector initially coincident with y-axis
 *  * ayy -- y-component of basis vector initially coincident with y-axis
 *  * ayz -- z-component of basis vector initially coincident with y-axis
 *  * azx -- x-component of basis vector initially coincident with z-axis
 *  * azy -- y-component of basis vector initially coincident with z-axis
 *  * azz -- z-component of basis vector initially coincident with z-axis
 *  * ox -- angular velocity in x-direction
 *  * oy -- angular velocity in y-direction
 *  * oz -- angular velocity in z-direction
 *  * oxdot -- angular acceleration in x-direction
 *  * oydot -- angular acceleration in y-direction
 *  * ozdot -- angular acceleration in z-direction
 *  * Fx -- hydrodynamic force in the x-direction
 *  * Fy -- hydrodynamic force in the y-direction
 *  * Fz -- hydrodynamic force in the z-direction
 *  * Lx -- hydrodynamic moment in the x-direction
 *  * Ly -- hydrodynamic moment in the y-direction
 *  * Lz -- hydrodynamic moment in the z-direction
 *  * aFx -- applied force in the x-direction
 *  * aFy -- applied force in the y-direction
 *  * aFz -- applied force in the z-direction
 *  * aLx -- applied moment in the x-direction
 *  * aLy -- applied moment in the y-direction
 *  * aLz -- applied moment in the z-direction
 *  * kFx -- applied spring force in the x-direction
 *  * kFy -- applied spring force in the y-direction
 *  * kFz -- applied spring force in the z-direction
 *  * iFx -- interaction force in the x-direction
 *  * iFy -- interaction force in the y-direction
 *  * iFz -- interaction force in the z-direction
 *  * iLx -- interaction moment in the x-direction
 *  * iLy -- interaction moment in the y-direction
 *  * iLz -- interaction moment in the z-direction
 *  * nodes -- the status of the nodes
 *  * rho -- particle density
 *  * E -- particle Young's modulus
 *  * sigma -- particle Poisson ratio (-1 < sigma <= 0.5)
 *  * cage -- the cage_struct that defines the particle within the domain
 *  * order -- the order above which to truncate the Lamb's series solution
 *  * rs -- the radius of integration for scalar products
 *  * ncoeff -- the number of Lamb's coefficients required order truncation
 *  * spring_k -- strength of spring pulling particle back to origin
 *  * spring_x -- x location of spring connection
 *  * spring_y -- y location of spring connection
 *  * spring_z -- z location of spring connection
 *  * spring_l -- the relaxed length of the spring
 *  * translating -- 1 if allowed to translate, 0 if not
 *  * rotating -- 1 if allowed to rotate, 0 if not
 *  * bin -- which bin the particle resides in
 *  * e_dry -- dry coefficient of restitution
 *  * l_rough -- particle surface roughness length
 *  * St -- particle-wall interaction Stokes number
 ******
 */

/****v* particle/nparts
 * NAME
 *  nparts
 * TYPE
 */
extern int nparts;
/*
 * PURPOSE
 *  Define the total number of particles in the domain.
 ******
 */

/****v* particle/interactionLength
 * NAME
 *  interactionLength
 * TYPE
 */
extern real interactionLength;
/*
 * PURPOSE
 *  Defines the particle-particle interaction length
 ******
 */

/****v* particle/binDom
 * NAME
 *  binDom
 * TYPE
 */
extern dom_struct binDom;
/*
 * PURPOSE
 *  A domain struct for the bin
 ******
 */

 /****v* particle/_binDom
 * NAME
 *  _binDom
 * TYPE
 */
extern dom_struct *_binDom;
/*
 * PURPOSE
 *  A domain struct for the bin (device)
 ******
 */

/****v* particle/parts
 * NAME
 *  parts
 * TYPE
 */
extern part_struct *parts;
/*
 * PURPOSE
 *  A list of all particles.
 ******
 */

/****v* particle/coeff_stride
 * NAME
 *  coeff_stride
 * TYPE
 */
extern int coeff_stride;
/*
 * PURPOSE
 *  Hold the stride length for pnm_re, pnm_im, phinm_re, phinm_im, chinm_re,
 *  chinm_im.
 ******
 */

/****v* particle/pnm_re
 * NAME
 *  pnm_re
 * TYPE
 */
extern real *pnm_re;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient p_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/pnm_im
 * NAME
 *  pnm_im
 * TYPE
 */
extern real *pnm_im;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient p_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/phinm_re
 * NAME
 *  phinm_re
 * TYPE
 */
extern real *phinm_re;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient phi_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/phinm_im
 * NAME
 *  phinm_im
 * TYPE
 */
extern real *phinm_im;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient phi_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/chinm_re
 * NAME
 *  chinm_re
 * TYPE
 */
extern real *chinm_re;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient chi_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/chinm_im
 * NAME
 *  chinm_im
 * TYPE
 */
extern real *chinm_im;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient chi_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/pnm_re0
 * NAME
 *  pnm_re0
 * TYPE
 */
extern real *pnm_re0;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient p_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/pnm_im0
 * NAME
 *  pnm_im0
 * TYPE
 */
extern real *pnm_im0;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient p_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/phinm_re0
 * NAME
 *  phinm_re0
 * TYPE
 */
extern real *phinm_re0;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient phi_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/phinm_im0
 * NAME
 *  phinm_im0
 * TYPE
 */
extern real *phinm_im0;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient phi_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/chinm_re0
 * NAME
 *  chinm_re0
 * TYPE
 */
extern real *chinm_re0;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient chi_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/chinm_im0
 * NAME
 *  chinm_im0
 * TYPE
 */
extern real *chinm_im0;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient chi_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/pnm_re00
 * NAME
 *  pnm_re00
 * TYPE
 */
extern real *pnm_re00;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient p_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/pnm_im00
 * NAME
 *  pnm_im00
 * TYPE
 */
extern real *pnm_im00;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient p_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/phinm_re00
 * NAME
 *  phinm_re00
 * TYPE
 */
extern real *phinm_re00;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient phi_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/phinm_im00
 * NAME
 *  phinm_im00
 * TYPE
 */
extern real *phinm_im00;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient phi_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/chinm_re00
 * NAME
 *  chinm_re00
 * TYPE
 */
extern real *chinm_re00;
/*
 * PURPOSE
 *  Contains all real parts of Lamb's coefficient chi_nm for each particle.  It
 *  strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/chinm_im00
 * NAME
 *  chinm_im00
 * TYPE
 */
extern real *chinm_im00;
/*
 * PURPOSE
 *  Contains all imaginary parts of Lamb's coefficient chi_nm for each particle.
 *  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_pnm_re
 * NAME
 *  _pnm_re
 * TYPE
 */
extern real **_pnm_re;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient p_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/_pnm_im
 * NAME
 *  _pnm_im
 * TYPE
 */
extern real **_pnm_im;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient p_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/_phinm_re
 * NAME
 *  _phinm_re
 * TYPE
 */
extern real **_phinm_re;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient phi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/_phinm_im
 * NAME
 *  _phinm_im
 * TYPE
 */
extern real **_phinm_im;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient phi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/_chinm_re
 * NAME
 *  _chinm_re
 * TYPE
 */
extern real **_chinm_re;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient chi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/_chinm_im
 * NAME
 *  _chinm_im
 * TYPE
 */
extern real **_chinm_im;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient chi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 ******
 */

/****v* particle/_pnm_re0
 * NAME
 *  _pnm_re0
 * TYPE
 */
extern real **_pnm_re0;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient p_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/_pnm_im0
 * NAME
 *  _pnm_im0
 * TYPE
 */
extern real **_pnm_im0;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient p_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/_phinm_re0
 * NAME
 *  _phinm_re0
 * TYPE
 */
extern real **_phinm_re0;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient phi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/_phinm_im0
 * NAME
 *  _phinm_im0
 * TYPE
 */
extern real **_phinm_im0;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient phi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/_chinm_re0
 * NAME
 *  _chinm_re0
 * TYPE
 */
extern real **_chinm_re0;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient chi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/_chinm_im0
 * NAME
 *  _chinm_im0
 * TYPE
 */
extern real **_chinm_im0;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient chi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous iteration.
 ******
 */

/****v* particle/_pnm_re00
 * NAME
 *  _pnm_re00
 * TYPE
 */
extern real **_pnm_re00;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient p_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_pnm_im00
 * NAME
 *  _pnm_im00
 * TYPE
 */
extern real **_pnm_im00;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient p_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_phinm_re00
 * NAME
 *  _phinm_re00
 * TYPE
 */
extern real **_phinm_re00;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient phi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_phinm_im00
 * NAME
 *  _phinm_im00
 * TYPE
 */
extern real **_phinm_im00;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient phi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_chinm_re00
 * NAME
 *  _chinm_re00
 * TYPE
 */
extern real **_chinm_re00;
/*
 * PURPOSE
 *  CUDA analog for real parts of Lamb's coefficient chi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_chinm_im00
 * NAME
 *  _chinm_im00
 * TYPE
 */
extern real **_chinm_im00;
/*
 * PURPOSE
 *  CUDA analog for imaginary parts of Lamb's coefficient chi_nm for each
 *  particle.  It strides by the coefficient for each particle varying fastest.
 *  This is the stored value at the previous timestep.
 ******
 */

/****v* particle/_parts
 * NAME
 *  _parts
 * TYPE
 */
extern part_struct **_parts;
/*
 * PURPOSE
 *  CUDA device analog for parts.  It contains pointers to arrays containing
 *  the particles in domain on which each device operates.  NOTE: for now,
 *  particles are designed to function on only one device.
 ******
 */

/****f* particle/parts_read_input()
 * NAME
 *  parts_read_input()
 * USAGE
 */
void parts_read_input(int turb);
/*
 * FUNCTION
 *  Read particle specfications from part.input.
 * ARGUMENTS
 *  * turb -- boolean operator for removing particles from turbulence precursor
 ******
 */

/****f* particle/parts_show_config()
 * NAME
 *  parts_show_config()
 * USAGE
 */
void parts_show_config(void);
/*
 * FUNCTION
 *  Write particle specifications to screen.
 ******
 */

/****f* particle/bin_show_config()
 * NAME
 *  bin_show_config()
 * USAGE
 */
void bin_show_config(void);
/*
 * FUNCTION
 *  Write bin specifications to screen.
 ******
 */

/****f* particle/parts_init()
 * NAME
 *  parts_init()
 * USAGE
 *
 */
int parts_init(void);
/*
 * FUNCTION
 *  Initialize the particles (build initial cages) and phase.
 * RESULT
 *  EXIT_SUCCESS if successful, EXIT_FAILURE otherwise.
 ******
 */

/****f* particle/init_cage()
 * NAME
 *  init_cage()
 * USAGE
 */
void init_cage(int part);
/*
 * FUNCTION
 *  Initialize the cage for particle part. The cage will be constructed
 *  by cuda_build_cages<<<>>>().
 * ARGUMENTS
 *  * part -- the particle for which to build a cage
 ******
 */

/****f* particle/flags_reset()
 * NAME
 *  flags_reset()
 * USAGE
 */
void flags_reset(void);
/*
 * FUNCTION
 *  Reinitializes flag arrays to no boundaries (1).
 ******
 */

/****f* particle/binDom_init()
 * NAME
 *  binDom_init()
 * USAGE
 *
 */
int binDom_init(void);
/*
 * FUNCTION
 *  Initialize the binDom structure
 * RESULT
 *  EXIT_SUCCESS if successful, EXIT_FAILURE otherwise.
 ******
 */

/****f* particle/parts_clean()
 * NAME
 *  parts_clean()
 * USAGE
 */
void parts_clean(void);
/*
 * FUNCTION
 *  Clean up.  Free any allocated host memory.
 ******
 */

/****f* particle/cuda_quad_check_nodes()
 * NAME
 *  cuda_quad_check_nodes()
 * USAGE
 */
void cuda_quad_check_nodes(int dev, real *node_t, real *node_p, int nnodes);
/*
 * FUNCTION
 *  Check if Lebedev quadrature nodes are inside another particle.
 * ARGUMENTS
 *  * dev -- the device on which to operate
 *  * node_t -- Lebedev quadrature node list theta component
 *  * node_p -- Lebedev quadrature node list phi component
 *  * nnodes -- the number of quadrature nodes in the list
 ******
 */

/****f* particle/cuda_quad_interp()
 * NAME
 *  cuda_quad_interp()
 * USAGE
 */
void cuda_quad_interp(int dev,
  real *node_t, real *node_p, int nnodes,
  real *pp, real *ur, real *ut, real *up);
/*
 * FUNCTION
 *  Interpolate fields to Lebedev quadrature nodes.
 * ARGUMENTS
 *  * dev -- the device on which to operate
 *  * node_t -- Lebedev quadrature node list theta component
 *  * node_p -- Lebedev quadrature node list phi component
 *  * nnodes -- the number of quadrature nodes in the list
 *  * pp -- the interpolated pressure at Lebedev quadrature nodes
 *  * ur -- the interpolated radial velocity at Lebedev quadrature nodes
 *  * ut -- the interpolated theta velocity at Lebedev quadrature nodes
 *  * up -- the interpolated phi velocity at Lebedev quadrature nodes
 ******
 */

/****f* particle/cuda_Lamb()
 * NAME
 *  cuda_Lamb()
 * USAGE
 */
void cuda_Lamb(void);
/*
 * FUNCTION
 *  Compute the Lamb's coefficients.
 ******
 */

/****f* particle/cuda_lamb_err()
 * NAME
 *  cuda_lamb_err()
 * USAGE
 */
real cuda_lamb_err(void);
/*
 * FUNCTION
 *  Compute the error between the current and previous sets of Lamb's
 *  coefficients.  It calculates the error then copies the current set
 *  of coefficients to the storage array to be saved for the next iteration
 *  error calculation.
 ******
 */

#endif
