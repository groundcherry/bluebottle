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

/****h* Bluebottle/domain
 * NAME
 *  domain
 * FUNCTION
 *  Low-level domain functions.
 ******
 */

#ifndef _DOMAIN_H
#define _DOMAIN_H

/****s* domain/grid_info
 * NAME
 *  grid_info
 * TYPE
 */
typedef struct grid_info {
  int is;
  int ie;
  int in;
  int isb;
  int ieb;
  int inb;
  int js;
  int je;
  int jn;
  int jsb;
  int jeb;
  int jnb;
  int ks;
  int ke;
  int kn;
  int ksb;
  int keb;
  int knb;
  int s1;
  int s1b;
  int s2;
  int s2b;
  int s3;
  int s3b;
  int _is;
  int _ie;
  int _in;
  int _isb;
  int _ieb;
  int _inb;
  int _js;
  int _je;
  int _jn;
  int _jsb;
  int _jeb;
  int _jnb;
  int _ks;
  int _ke;
  int _kn;
  int _ksb;
  int _keb;
  int _knb;
  int _s1;
  int _s1b;
  int _s2;
  int _s2b;
  int _s3;
  int _s3b;
} grid_info;
/*
 * PURPOSE
 *  Carry information related to the different discretization grids.
 * MEMBERS
 *  * is -- the domain start index in the x-direction
 *  * ie -- the domain end index in the x-direction (one greater than the
 *    index associated with the end index so a < can be used in loops instead
 *    of a <=)
 *  * in -- the number of elements in the domain in the x-direction
 *  * isb -- the domain start index in the x-direction plus boundary ghost
 *    elements
 *  * ieb -- the domain end index in the x-direction plus boundary ghost
 *    elements (one greater than the index associated with the end index so a
 *    < can be used in loops instead of a <=)
 *  * inb -- the number of elements in the domain in the x-direction plus
 *    the boundary ghost elements
 *  * js -- the domain start index in the y-direction
 *  * je -- the domain end index in the y-direction (one greater than the
 *    index associated with the end index so a < can be used in loops instead
 *    of a <=)
 *  * jn -- the number of elements in the domain in the y-direction
 *  * jsb -- the domain start index in the y-direction plus boundary ghost
 *    elements
 *  * jeb -- the domain end index in the y-direction plus boundary ghost
 *    elements (one greater than the index associated with the end index so a
 *    < can be used in loops instead of a <=)
 *  * jnb -- the number of elements in the domain in the y-direction plus
 *    the boundary ghost elements
 *  * ks -- the domain start index in the z-direction
 *  * ke -- the domain end index in the z-direction (one greater than the
 *    index associated with the end index so a < can be used in loops instead
 *    of a <=)
 *  * kn -- the number of elements in the domain in the z-direction
 *  * ksb -- the domain start index in the z-direction plus boundary ghost
 *    elements
 *  * keb -- the domain end index in the z-direction plus boundary ghost
 *    elements (one greater than the index associated with the end index so a
 *    < can be used in loops instead of a <=)
 *  * knb -- the number of elements in the domain in the z-direction plus
 *    the boundary ghost elements
 *  * s1 -- the looping stride length for the fastest-changing variable (x)
 *  * s1b -- the looping stride length for the fastest-changing variable (x)
 *    plus the boundary ghost elements
 *  * s2 -- the looping stride length for the second-changing variable (y)
 *  * s2b -- the looping stride length for the second-changing variable (y)
 *    plus the boundary ghost elements
 *  * s3 -- the looping stride length for the slowest-changing variable (z)
 *  * s3b -- the looping stride length for the slowest-changing variable (z)
 *    plus the boundary ghost elements
 *  * _is -- the domain start index in the x-direction on the device
 *  * _ie -- the domain end index in the x-direction (one greater than the
 *    index associated with the end index so a < can be used in loops instead
 *    of a <=) on the device
 *  * _in -- the number of elements in the domain in the x-direction on the
 *    device
 *  * _isb -- the domain start index in the x-direction plus boundary ghost
 *    elements on the device
 *  * _ieb -- the domain end index in the x-direction plus boundary ghost
 *    elements (one greater than the index associated with the end index so a
 *    < can be used in loops instead of a <=) on the device
 *  * _inb -- the number of elements in the domain in the x-direction plus
 *    the boundary ghost elements on the device
 *  * _js -- the domain start index in the y-direction on the device
 *  * _je -- the domain end index in the y-direction (one greater than the
 *    index associated with the end index so a < can be used in loops instead
 *    of a <=) on the device
 *  * _jn -- the number of elements in the domain in the y-direction on the
 *    device
 *  * _jsb -- the domain start index in the y-direction plus boundary ghost
 *    elements on the device
 *  * _jeb -- the domain end index in the y-direction plus boundary ghost
 *    elements (one greater than the index associated with the end index so a
 *    < can be used in loops instead of a <=) on the device
 *  * _jnb -- the number of elements in the domain in the y-direction plus
 *    the boundary ghost elements on the device
 *  * _ks -- the domain start index in the z-direction on the device
 *  * _ke -- the domain end index in the z-direction (one greater than the
 *    index associated with the end index so a < can be used in loops instead
 *    of a <=) on the device
 *  * _kn -- the number of elements in the domain in the z-direction
 *    on the device
 *  * _ksb -- the domain start index in the z-direction plus boundary ghost
 *    elements on the device
 *  * _keb -- the domain end index in the z-direction plus boundary ghost
 *    elements (one greater than the index associated with the end index so a
 *    < can be used in loops instead of a <=) on the device
 *  * _knb -- the number of elements in the domain in the z-direction plus
 *    the boundary ghost elements on the device
 *  * _s1 -- the looping stride length for the fastest-changing variable (x)
 *    on the device
 *  * _s1b -- the looping stride length for the fastest-changing variable (x)
 *    plus the boundary ghost elements on the device
 *  * _s2 -- the looping stride length for the second-changing variable (y)
 *    on the device
 *  * _s2b -- the looping stride length for the second-changing variable (y)
 *    plus the boundary ghost elements on the device
 *  * _s3 -- the looping stride length for the slowest-changing variable (z)
 *    on the device
 *  * _s3b -- the looping stride length for the slowest-changing variable (z)
 *    plus the boundary ghost elements on the device
 ******
 */

/****v* domain/nsubdom
 * NAME
 *  nsubdom
 * TYPE
 */
extern int nsubdom;
/*
 * PURPOSE
 *  Number of subdomains into which the domain should be decomposed for GPU
 *  subdivision of the domain.
 ******
 */

/****s* domain/dom_struct
 * NAME
 *  dom_struct
 * TYPE
 */
typedef struct dom_struct {
  grid_info Gcc;
  grid_info Gfx;
  grid_info Gfy;
  grid_info Gfz;
  real xs;
  real xe;
  real xl;
  int xn;
  real dx;
  real ys;
  real ye;
  real yl;
  int yn;
  real dy;
  real zs;
  real ze;
  real zl;
  int zn;
  real dz;
  int E;
  int W;
  int N;
  int S;
  int T;
  int B;
} dom_struct;
/*
 * PURPOSE
 *  Carry information related to a subdomain.
 * MEMBERS
 *  * Gcc -- cell-centered grid information
 *  * Gfx -- FACE_X face-centered grid information
 *  * Gfy -- FACE_Y face-centered grid information
 *  * Gfz -- FACE_Z face-centered grid information
 *  * xs -- physical start position in the x-direction
 *  * xe -- physical end position in the x-direction
 *  * xl -- physical length of the subdomain in the x-direction
 *  * xn -- number of discrete cells in the x-direction
 *  * dx -- cell size in the x-direction
 *  * ys -- physical start position in the y-direction
 *  * ye -- physical end position in the y-direction
 *  * yl -- physical length of the subdomain in the y-direction
 *  * yn -- number of discrete cells in the y-direction
 *  * dy -- cell size in the y-direction
 *  * zs -- physical start position in the z-direction
 *  * ze -- physical end position in the z-direction
 *  * zl -- physical length of the subdomain in the z-direction
 *  * zn -- number of discrete cells in the z-direction
 *  * dz -- cell size in the z-direction
 *  * E -- the subdomain adjacent to the east face of the cell
 *    (E = -1 if the face is a domain boundary)
 *  * W -- the subdomain adjacent to the west face of the cell
 *    (W = -1 if the face is a domain boundary)
 *  * N -- the subdomain adjacent to the north face of the cell
 *    (N = -1 if the face is a domain boundary)
 *  * S -- the subdomain adjacent to the south face of the cell
 *    (S = -1 if the face is a domain boundary)
 *  * T -- the subdomain adjacent to the top face of the cell
 *    (T = -1 if the face is a domain boundary)
 *  * B -- the subdomain adjacent to the bottom face of the cell
 *    (B = -1 if the face is a domain boundary)
 ******
 */

/****f* domain/domain_read_input()
 * NAME
 *  domain_read_input()
 * USAGE
 */
void domain_read_input(void);
/*
 * FUNCTION
 *  Read domain specifications and simulation parameters from flow.config.
 ******
 */

/****f* domain/turb_read_input()
 * NAME
 *  turb_read_input()
 * USAGE
 */
void turb_read_input(void);
/*
 * FUNCTION
 *  Read turbulent precursor domain specifications and simulation parameters
 *  from turb.config.
 ******
 */

/****f* domain/domain_show_config()
 * NAME
 *  domain_show_config()
 * USAGE
 */
void domain_show_config(void);
/*
 * FUNCTION
 *  Write domain specifications and simulation parameters to screen.
 ******
 */

/****f* domain/domain_init()
 * NAME
 *  domain_init()
 * USAGE
 */
int domain_init(void);
/*
 * FUNCTION
 *  Initialize the domain on the host.
 * RESULT
 *  EXIT_SUCCESS if successful, EXIT_FAILURE otherwise.
 ******
 */

/****f* domain/domain_init_turb()
 * NAME
 *  domain_init_turb()
 * USAGE
 */
int domain_init_turb(void);
/*
 * FUNCTION
 *  Initialize the turbulent domain on the host.
 * RESULT
 *  EXIT_SUCCESS if successful, EXIT_FAILURE otherwise.
 ******
 */

/****f* domain/domain_clean()
 * NAME
 *  domain_clean()
 * USAGE
 */
void domain_clean(void);
/*
 * FUNCTION
 *   Clean up.  Free any allocated host memory.
 ******
 */

/****f* domain/compute_vel_BC()
 * NAME
 *  compute_vel_BC()
 * USAGE
 */
void compute_vel_BC(void);
/*
 * FUNCTION
 *  Set up the Dirichlet velocity boundary conditions for this time step
 *  given user-inputted maximum velocity and acceleration.
 ******
 */

/****f* domain/count_mem()
 * NAME
 *  count_mem()
 * USAGE
 */
void count_mem(void);
/*
 * FUNCTION
 *  Counts the total memory usage on the CPU and GPU
 ******
 */ 
#endif
