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

#include "cuda_quadrature.h"

__device__ void rtp2xyz(real r, real theta, real phi, real *x, real *y, real *z)
{
  *x = r * sin(theta) * cos(phi);
  *y = r * sin(theta) * sin(phi);
  *z = r * cos(theta);
}

__device__ void cart2sphere(real u, real v, real w, real theta, real phi,
  real *ur, real *ut, real *up)
{
  real st = sin(theta);
  real ct = cos(theta);
  real sp = sin(phi);
  real cp = cos(phi);

  *ur = st * (u * cp + v * sp) + w * ct;
  *ut = ct * (u * cp + v * sp) - w * st;
  *up = -u * sp + v * cp;
}

__global__ void check_nodes(int nparts, part_struct *parts, dom_struct *dom,
  real *theta, real *phi, int nnodes, BC bc)
{
  int node = threadIdx.x;
  int part = blockIdx.x;

  // convert node (r, theta, phi) to (x, y, z)
  real xp, yp, zp;  // Cartesian radial vector
  real x, y, z;   // Cartesian location of node
  rtp2xyz(parts[part].rs, theta[node], phi[node], &xp, &yp, &zp);

  // shift from particle center
  x = xp + parts[part].x;
  y = yp + parts[part].y;
  z = zp + parts[part].z;

  // start off with all -1's
  parts[part].nodes[node] = -1;

  // check if the node is interfered with by a wall
  // compute distance between node and walls
  // set equal to some number to identify which wall is interefering
  if(x - dom->xs < 0) {
    if(bc.uW == DIRICHLET || bc.vW == DIRICHLET || bc.wW == DIRICHLET)
      if(parts[part].nodes[node] == -1)
        parts[part].nodes[node] = -10;
  } if(x - dom->xe > 0) {
    if(bc.uE == DIRICHLET || bc.vE == DIRICHLET || bc.wE == DIRICHLET)
      if(parts[part].nodes[node] == -1)
        parts[part].nodes[node] = -11;
  } if(y - dom->ys < 0) {
    if(bc.uS == DIRICHLET || bc.vS == DIRICHLET || bc.wS == DIRICHLET)
      if(parts[part].nodes[node] == -1)
        parts[part].nodes[node] = -12;
  } if(y - dom->ye > 0) {
    if(bc.uN == DIRICHLET || bc.vN == DIRICHLET || bc.wN == DIRICHLET)
      if(parts[part].nodes[node] == -1)
        parts[part].nodes[node] = -13;
  } if(z - dom->zs < 0) {
    if(bc.uB == DIRICHLET || bc.vB == DIRICHLET || bc.wB == DIRICHLET)
      if(parts[part].nodes[node] == -1)
        parts[part].nodes[node] = -14;
  } if(z - dom->ze > 0) {
    if(bc.uT == DIRICHLET || bc.vT == DIRICHLET || bc.wT == DIRICHLET)
      if(parts[part].nodes[node] == -1)
        parts[part].nodes[node] = -15;
  }

}

__global__ void interpolate_nodes(real *p0, real *p, real *u, real *v, real *w,
  real rho_f, real nu, gradP_struct gradP,
  part_struct *parts, dom_struct *dom, real *theta, real *phi, int nnodes,
  real *pp, real *ur, real *ut, real *up, real dt0, real dt, BC bc)
{
  int node = threadIdx.x;
  int part = blockIdx.x;

  // the node number of the intersecting node
  int intnode = parts[part].nodes[node];
  if(intnode < 0) intnode = part;

  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;

  real ox = parts[part].ox;
  real oy = parts[part].oy;
  real oz = parts[part].oz;
  real oxdot = parts[part].oxdot;
  real oydot = parts[part].oydot;
  real ozdot = parts[part].ozdot;
  real udot = parts[part].udot;
  real vdot = parts[part].vdot;
  real wdot = parts[part].wdot;

  real uu, vv, ww;  // temporary nodes for Cartesian result of interpolation
  real uuwall, vvwall, wwwall;

  // convert node (r, theta, phi) to (x, y, z)
  real xp, yp, zp;  // Cartesian radial vector
  real x, y, z;   // Cartesian location of node
  rtp2xyz(parts[part].rs, theta[node], phi[node], &xp, &yp, &zp);

  // shift from particle center
  x = xp + parts[part].x;
  y = yp + parts[part].y;
  z = zp + parts[part].z;

  if(x < dom->xs && bc.uW == PERIODIC) x = x + dom->xl;
  else if(x > dom->xe && bc.uE == PERIODIC) x = x - dom->xl;
  if(y < dom->ys && bc.vS == PERIODIC) y = y + dom->yl;
  else if(y > dom->ye && bc.vN == PERIODIC) y = y - dom->yl;
  if(z < dom->zs && bc.wB == PERIODIC) z = z + dom->zl;
  else if(z > dom->ze && bc.wT == PERIODIC) z = z - dom->zl;

  __syncthreads();

  // find index of cell containing node
  int i = floor((x - dom->xs) * ddx) + DOM_BUF;
  int j = floor((y - dom->ys) * ddy) + DOM_BUF;
  int k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gcc.is) i = dom->Gcc.is;
  if(j < dom->Gcc.js) j = dom->Gcc.js;
  if(k < dom->Gcc.ks) k = dom->Gcc.ks;
  if(i > dom->Gcc.ie-1) i = dom->Gcc.ie-1;
  if(j > dom->Gcc.je-1) j = dom->Gcc.je-1;
  if(k > dom->Gcc.ke-1) k = dom->Gcc.ke-1;
  int C = i + j*dom->Gcc.s1b + k*dom->Gcc.s2b;
  // Cartesian location of center of cell
  real xx = (i-0.5) * dom->dx + dom->xs;
  real yy = (j-0.5) * dom->dy + dom->ys;
  real zz = (k-0.5) * dom->dz + dom->zs;

  // interpolate pressure
  real pc = p[C];
  real pw = p[C-1];
  real pe = p[C+1];
  real ps = p[C-dom->Gcc.s1b];
  real pn = p[C+dom->Gcc.s1b];
  real pb = p[C-dom->Gcc.s2b];
  real pt = p[C+dom->Gcc.s2b];
  real dpdx = 0.5*(pe - pw) * ddx;
  real dpdy = 0.5*(pn - ps) * ddy;
  real dpdz = 0.5*(pt - pb) * ddz;
  pp[node+nnodes*part] = pc + dpdx*(x-xx) + dpdy*(y-yy) + dpdz*(z-zz);
  // switch to particle rest frame
  real ocrossr2 = (oy*zp - oz*yp) * (oy*zp - oz*yp);
  ocrossr2 += (ox*zp - oz*xp) * (ox*zp - oz*xp);
  ocrossr2 += (ox*yp - oy*xp) * (ox*yp - oy*xp);
  real rhoV = rho_f;
  real accdotr = (-gradP.x/rhoV - udot)*xp + (-gradP.y/rhoV - vdot)*yp
    + (-gradP.z/rhoV - wdot)*zp;
  pp[node+nnodes*part] -= 0.5 * rho_f * ocrossr2 + rho_f * accdotr;
  // zero if this node intersects wall
  pp[node+nnodes*part] = (parts[part].nodes[node]==-1)*pp[node+part*nnodes];

  // interpolate velocities
  // don't work with cell-center anymore;
  // find closest cell face in x-direction

  // interpolate u-velocity
  i = round((x - dom->xs) * ddx) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gfx.is) i = dom->Gfx.is;
  if(j < dom->Gfx.js) j = dom->Gfx.js;
  if(k < dom->Gfx.ks) k = dom->Gfx.ks;
  if(i > dom->Gfx.ie-1) i = dom->Gfx.ie-1;
  if(j > dom->Gfx.je-1) j = dom->Gfx.je-1;
  if(k > dom->Gfx.ke-1) k = dom->Gfx.ke-1;
  xx = (i-DOM_BUF) * dom->dx + dom->xs;
  yy = (j-0.5) * dom->dy + dom->ys;
  zz = (k-0.5) * dom->dz + dom->zs;
  C = i + j*dom->Gfx.s1b + k*dom->Gfx.s2b;
  real dudx = 0.5*(u[C+1] - u[C-1]) * ddx;
  real dudy = 0.5*(u[C+dom->Gfx.s1b] - u[C-dom->Gfx.s1b]) * ddy;
  real dudz = 0.5*(u[C+dom->Gfx.s2b] - u[C-dom->Gfx.s2b]) * ddz;
  uu = u[C] + dudx * (x - xx) + dudy * (y - yy) + dudz * (z - zz);
  // set uuwall equal to interfering wall u-velocity
  uuwall = (parts[part].nodes[node] == -10)*bc.uWD
            + (parts[part].nodes[node] == -11)*bc.uED
            + (parts[part].nodes[node] == -12)*bc.uSD
            + (parts[part].nodes[node] == -13)*bc.uND
            + (parts[part].nodes[node] == -14)*bc.uBD
            + (parts[part].nodes[node] == -15)*bc.uTD;
  // switch to particle rest frame
  real rs3 = parts[part].rs*parts[part].rs*parts[part].rs;
  real rs5 = rs3*parts[part].rs*parts[part].rs;
  real a5 = parts[part].r*parts[part].r*parts[part].r*parts[part].r*parts[part].r;
  real ocrossr_x = oy*zp - oz*yp;
  real odotcrossr_x = oydot*zp - ozdot*yp;
  uu -= parts[part].u + ocrossr_x;
  uu -= 0.1/nu *(rs5-a5)/rs3 * odotcrossr_x;
  uuwall -= parts[part].u + ocrossr_x;
  uuwall -= 0.1/nu *(rs5-a5)/rs3 * odotcrossr_x;
  // set actual node value based on whether it is interfered with
  uu = (parts[part].nodes[node]==-1)*uu
    + (parts[part].nodes[node]<-1)*uuwall;
//printf("uu = %f uuwall = %f\n", uu + parts[part].u + ocrossr_x + 0.1 / nu / rs3 * (rs5 - r5) * odotcrossr_x, uuwall + parts[part].u + ocrossr_x + 0.1 / nu / rs3 * (rs5 - r5) * odotcrossr_x);

  // interpolate v-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF;
  j = round((y - dom->ys) * ddy) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gfy.is) i = dom->Gfy.is;
  if(j < dom->Gfy.js) j = dom->Gfy.js;
  if(k < dom->Gfy.ks) k = dom->Gfy.ks;
  if(i > dom->Gfy.ie-1) i = dom->Gfy.ie-1;
  if(j > dom->Gfy.je-1) j = dom->Gfy.je-1;
  if(k > dom->Gfy.ke-1) k = dom->Gfy.ke-1;
  xx = (i-0.5) * dom->dx + dom->xs;
  yy = (j-DOM_BUF) * dom->dy + dom->ys;
  zz = (k-0.5) * dom->dz + dom->zs;
  C = i + j*dom->Gfy.s1b + k*dom->Gfy.s2b;
  real dvdx = 0.5*(v[C+1] - v[C-1]) * ddx;
  real dvdy = 0.5*(v[C+dom->Gfy.s1b] - v[C-dom->Gfy.s1b]) * ddy;
  real dvdz = 0.5*(v[C+dom->Gfy.s2b] - v[C-dom->Gfy.s2b]) * ddz;
  vv = v[C] + dvdx * (x - xx) + dvdy * (y - yy) + dvdz * (z - zz);
  // set vvwall equal to interfering wall v-velocity
  vvwall = (parts[part].nodes[node] == -10)*bc.vWD
            + (parts[part].nodes[node] == -11)*bc.vED
            + (parts[part].nodes[node] == -12)*bc.vSD
            + (parts[part].nodes[node] == -13)*bc.vND
            + (parts[part].nodes[node] == -14)*bc.vBD
            + (parts[part].nodes[node] == -15)*bc.vTD;
  // switch to particle rest frame
  real ocrossr_y = -(ox*zp - oz*xp);
  real odotcrossr_y = -(oxdot*zp - ozdot*xp);
  vv -= parts[part].v + ocrossr_y;
  vv -= 0.1/nu *(rs5-a5)/rs3 * odotcrossr_y;
  vvwall -= parts[part].v + ocrossr_y;
  vvwall -= 0.1/nu *(rs5-a5)/rs3 * odotcrossr_y;
  // set actual node value based on whether it is interfered with
  vv = (parts[part].nodes[node]==-1)*vv
    + (parts[part].nodes[node]<-1)*vvwall;

  // interpolate w-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = round((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gfz.is) i = dom->Gfz.is;
  if(j < dom->Gfz.js) j = dom->Gfz.js;
  if(k < dom->Gfz.ks) k = dom->Gfz.ks;
  if(i > dom->Gfz.ie-1) i = dom->Gfz.ie-1;
  if(j > dom->Gfz.je-1) j = dom->Gfz.je-1;
  if(k > dom->Gfz.ke-1) k = dom->Gfz.ke-1;
  xx = (i-0.5) * dom->dx + dom->xs;
  yy = (j-0.5) * dom->dy + dom->ys;
  zz = (k-DOM_BUF) * dom->dz + dom->zs;
  C = i + j*dom->Gfz.s1b + k*dom->Gfz.s2b;
  real dwdx = 0.5*(w[C+1] - w[C-1]) * ddx;
  real dwdy = 0.5*(w[C+dom->Gfz.s1b] - w[C-dom->Gfz.s1b]) * ddy;
  real dwdz = 0.5*(w[C+dom->Gfz.s2b] - w[C-dom->Gfz.s2b]) * ddz;
  ww = w[C] + dwdx * (x - xx) + dwdy * (y - yy) + dwdz * (z - zz);
  // set uuwall equal to interfering wall u-velocity
  wwwall = (parts[part].nodes[node] == -10)*bc.wWD
            + (parts[part].nodes[node] == -11)*bc.wED
            + (parts[part].nodes[node] == -12)*bc.wSD
            + (parts[part].nodes[node] == -13)*bc.wND
            + (parts[part].nodes[node] == -14)*bc.wBD
            + (parts[part].nodes[node] == -15)*bc.wTD;
  // switch to particle rest frame
  real ocrossr_z = ox*yp - oy*xp;
  real odotcrossr_z = oxdot*yp - oydot*xp;
  ww -= parts[part].w + ocrossr_z;
  ww -= 0.1/nu *(rs5-a5)/rs3 * odotcrossr_z;
  wwwall -= parts[part].w + ocrossr_z;
  wwwall -= 0.1/nu *(rs5-a5)/rs3 * odotcrossr_z;
  // set actual node value based on whether it is interfered with
  ww = (parts[part].nodes[node]==-1)*ww
    + (parts[part].nodes[node]<-1)*wwwall;

  // convert (uu, vv, ww) to (u_r, u_theta, u_phi) and write to node arrays
  cart2sphere(uu, vv, ww, theta[node], phi[node],
    &ur[node+part*nnodes], &ut[node+part*nnodes], &up[node+part*nnodes]);

//printf("%e %e u = %e v = %e w = %e\n", theta[node], phi[node], uu,vv,ww);
}

__device__ real nnm(int n, int m)
{
  real fact_top = 1;
  real fact_bot = 1;

  for(int i = 1; i <= (n-m); i++) fact_top *= (real)i;
  for(int i = 1; i <= (n+m); i++) fact_bot *= (real)i;

  return sqrt((2.*n+1.) / 4. / PI * fact_top / fact_bot);
}

__device__ real pnm(int n, int m, real theta)
{
  real x = cos(theta);
  real y = sin(theta);

  switch(n) {
    case 0: return 1;
    case 1:
      switch(m) {
        //case -1: return -0.5*y;
        case 0: return x;
        case 1: return -y;
      }
    case 2:
      switch(m) {
        //case -2: return 0.125*y*y;
        //case -1: return -0.5*x*y;
        case 0: return 0.5*(3.*x*x - 1.);
        case 1: return -3.*x*y;
        case 2: return 3.*y*y;
      }
    case 3:
      switch(m) {
        //case -3: return -0.02083333333333*y*y*y;
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
        //case -3: return -0.02083333333333*x*y*y*y*y;
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
        //case -5: return -0.000260416666667*y*y*y*y*y;
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
  real lambrelax)
{
  int node = threadIdx.x;
  int part = blockIdx.x;
  int coeff = blockIdx.y;
  real ars = parts[part].r / parts[part].rs;
  real rsa = parts[part].rs / parts[part].r;

  int i;  // iterator

  if(coeff < parts[part].ncoeff) {
    // calculate integrand at each node
    int j = part*nnodes*ncoeffs + coeff*nnodes + node;

    int n = nn[coeff];
    int m = mm[coeff];
    real theta = node_t[node];
    real phi = node_p[node];
    real N_nm = nnm(n,m);
    real P_nm = pnm(n,m,theta);
    real P_n1m = pnm(n+1.,m,theta);
    real dPdt = (n-m+1.)*P_n1m-(n+1.)*cos(theta)*P_nm;
    real dPdp = m*P_nm;

    int_Yp_re[j] = N_nm*P_nm*pp[node+part*nnodes]*cos(m*phi);
    int_Yp_im[j] = -N_nm*P_nm*pp[node+part*nnodes]*sin(m*phi);

    int_rDYu_re[j] = N_nm/sin(theta)*(dPdt*ut[node+part*nnodes]*cos(m*phi)
      - dPdp*up[node+part*nnodes]*sin(m*phi));
    int_rDYu_im[j] = N_nm/sin(theta)*(-dPdt*ut[node+part*nnodes]*sin(m*phi)
      - dPdp*up[node+part*nnodes]*cos(m*phi));

    int_xXDYu_re[j] = N_nm/sin(theta)*(dPdp*ut[node+part*nnodes]*sin(m*phi)
      + dPdt*up[node+part*nnodes]*cos(m*phi));
    int_xXDYu_im[j] = N_nm/sin(theta)*(dPdp*ut[node+part*nnodes]*cos(m*phi)
      - dPdt*up[node+part*nnodes]*sin(m*phi));

//if(n == 1 && m == 1) {
//  printf("int_Yp_re(%f, %f) = %e\n", theta, phi, int_Yp_re[j]);
//  printf("int_rDYu_re(%f, %f) = %e\n", theta, phi, int_rDYu_re[j]);
//  printf("int_xXDYu_re(%f, %f) = %e ut = %e up = %e\n", theta, phi, int_xXDYu_re[j], ut[node+part*nnodes], up[node+part*nnodes]);
//}

    __syncthreads();

    // compute scalar products
    // put sum into first node position for each coeff for each particle
    if(node == 0) {
      int_Yp_re[j] *= A1;
      int_Yp_im[j] *= A1;
      int_rDYu_re[j] *= A1;
      int_rDYu_im[j] *= A1;
      int_xXDYu_re[j] *= A1;
      int_xXDYu_im[j] *= A1;
      for(i = 1; i < 6; i++) {
        int_Yp_re[j] += A1 * int_Yp_re[j+i];
        int_Yp_im[j] += A1 * int_Yp_im[j+i];
        int_rDYu_re[j] += A1 * int_rDYu_re[j+i];
        int_rDYu_im[j] += A1 * int_rDYu_im[j+i];
        int_xXDYu_re[j] += A1 * int_xXDYu_re[j+i];
        int_xXDYu_im[j] += A1 * int_xXDYu_im[j+i];
      }
      for(i = 6; i < 18; i++) {
        int_Yp_re[j] += A2 * int_Yp_re[j+i];
        int_Yp_im[j] += A2 * int_Yp_im[j+i];
        int_rDYu_re[j] += A2 * int_rDYu_re[j+i];
        int_rDYu_im[j] += A2 * int_rDYu_im[j+i];
        int_xXDYu_re[j] += A2 * int_xXDYu_re[j+i];
        int_xXDYu_im[j] += A2 * int_xXDYu_im[j+i];
      }
      for(i = 18; i < 26; i++) {
        int_Yp_re[j] += A3 * int_Yp_re[j+i];
        int_Yp_im[j] += A3 * int_Yp_im[j+i];
        int_rDYu_re[j] += A3 * int_rDYu_re[j+i];
        int_rDYu_im[j] += A3 * int_rDYu_im[j+i];
        int_xXDYu_re[j] += A3 * int_xXDYu_re[j+i];
        int_xXDYu_im[j] += A3 * int_xXDYu_im[j+i];
      }

      /*for(i = 26; i < 50; i++) {
        int_Yp_re[j] += B * int_Yp_re[j+i];
        int_Yp_im[j] += B * int_Yp_im[j+i];
        int_rDYu_re[j] += B * int_rDYu_re[j+i];
        int_rDYu_im[j] += B * int_rDYu_im[j+i];
        int_xXDYu_re[j] += B * int_xXDYu_re[j+i];
        int_xXDYu_im[j] += B * int_xXDYu_im[j+i];
      }
*/

//if(n == 1 && m == 1) {
//  printf("int_Yp_re = %e\n", int_Yp_re[j]);
//  printf("int_rDYu_re = %e\n", int_rDYu_re[j]);
//  printf("int_xXDYu_re = %e\n", int_xXDYu_re[j]);
//}
      
#ifdef TEST
      real relax = 1.0;
#else
      real relax = lambrelax;
#endif
      if(n == 0) {
        pnm_re[stride*part+coeff] = pnm_re0[stride*part+coeff]
          + relax*(parts[part].r*parts[part].r/mu/nu*int_Yp_re[j]*pow(ars,n)
          - pnm_re0[stride*part+coeff]);
        pnm_im[stride*part+coeff] = pnm_im0[stride*part+coeff]
          + relax*(parts[part].r*parts[part].r/mu/nu*int_Yp_im[j]*pow(ars,n)
          - pnm_im0[stride*part+coeff]);
        phinm_re[stride*part+coeff] = 0.;
        phinm_im[stride*part+coeff] = 0.;
        chinm_re[stride*part+coeff] = 0.;
        chinm_im[stride*part+coeff] = 0.;
      } else {
        // calculate p_nm and phi_nm
        real A = (1.-0.5*n*(2.*n-1.)/(n+1.)*pow(ars,2.*n+1.))*pow(rsa,n);
        real B = n*(2.*n-1.)*(2.*n+1.)/(n+1.)*pow(ars,n+1.);
        real C = 0.25*n*(2.*(n+3.)/(2.*n+3.)
          + (n-2.-n*(2.*n+1.)/(2.*n+3.)*ars*ars)*pow(ars,2.*n+1.))*pow(rsa,n+1.);
        real D = n*(n+1.+0.5*((n-2.)*(2.*n+1.)*rsa*rsa
          - n*(2.*n-1.))*pow(ars,2.*n+1.))*pow(rsa,n-1.);

        pnm_re[stride*part+coeff] = (parts[part].r*parts[part].r/mu/nu
          *int_Yp_re[j]*D + parts[part].r/nu*int_rDYu_re[j]*B) / (A*D+B*C);
        pnm_im[stride*part+coeff] = (parts[part].r*parts[part].r/mu/nu
          *int_Yp_im[j]*D + parts[part].r/nu*int_rDYu_im[j]*B) / (A*D+B*C);

        phinm_re[stride*part+coeff] = (parts[part].r/nu*int_rDYu_re[j]*A
          - parts[part].r*parts[part].r/mu/nu*int_Yp_re[j]*C) / (A*D+B*C);
        phinm_im[stride*part+coeff] = (parts[part].r/nu*int_rDYu_im[j]*A
          - parts[part].r*parts[part].r/mu/nu*int_Yp_im[j]*C) / (A*D+B*C);

        // calculate chi_nm
        real E = n*(n+1.)*(pow(ars,2.*n+1.)-1.)*pow(rsa, n);
        chinm_re[stride*part+coeff] = parts[part].r/nu*int_xXDYu_re[j] / E;
        chinm_im[stride*part+coeff] = parts[part].r/nu*int_xXDYu_im[j] / E;

        // apply underrelaxation
        pnm_re[stride*part+coeff] = pnm_re0[stride*part+coeff]*(1.-relax)
          + relax*pnm_re[stride*part+coeff];
        pnm_im[stride*part+coeff] = pnm_im0[stride*part+coeff]*(1.-relax)
          + relax*pnm_im[stride*part+coeff];
        phinm_re[stride*part+coeff] = phinm_re0[stride*part+coeff]*(1.-relax)
          + relax*phinm_re[stride*part+coeff];
        phinm_im[stride*part+coeff] = phinm_im0[stride*part+coeff]*(1.-relax)
          + relax*phinm_im[stride*part+coeff];
        chinm_re[stride*part+coeff] = chinm_re0[stride*part+coeff]*(1.-relax)
          + relax*chinm_re[stride*part+coeff];
        chinm_im[stride*part+coeff] = chinm_im0[stride*part+coeff]*(1.-relax)
          + relax*chinm_im[stride*part+coeff];

//printf("pnm_re(%d,%d) = %e\n", n,m, pnm_re[stride*part+coeff]);
      }
    }
  }
}

__global__ void cuda_calc_forces(dom_struct *dom, part_struct *parts,
  int nparts, gradP_struct gradP,
  real rho_f, real mu, real nu, int stride,
  real *pnm_re, real *pnm_im,
  real *phinm_re, real *phinm_im,
  real *chinm_re, real *chinm_im)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x; // particle number

  if(pp < nparts) {
    real vol = 4./3. * PI *  parts[pp].r*parts[pp].r*parts[pp].r;
    real N10 = sqrt(3./4./PI);
    real N11 = sqrt(3./8./PI);

    parts[pp].Fx = rho_f * vol * (parts[pp].udot + gradP.x/rho_f)
      - PI * mu * nu * 2.*N11 * (pnm_re[stride*pp + 2]
      + 6.*phinm_re[stride*pp + 2]);
    parts[pp].Fy = rho_f * vol * (parts[pp].vdot + gradP.y/rho_f)
      + PI * mu * nu * 2.*N11 * (pnm_im[stride*pp + 2]
      + 6.*phinm_im[stride*pp + 2]);
    parts[pp].Fz = rho_f * vol * (parts[pp].wdot + gradP.z/rho_f)
      + PI * mu * nu * N10 * (pnm_re[stride*pp + 1]
      + 6.*phinm_re[stride*pp + 1]);

    parts[pp].Lx = rho_f * vol * parts[pp].r*parts[pp].r * parts[pp].oxdot
      - 8. * PI * mu * nu * 2.*N11 * parts[pp].r * chinm_re[stride*pp + 2];
    parts[pp].Ly = rho_f * vol * parts[pp].r*parts[pp].r * parts[pp].oydot
      + 8. * PI * mu * nu * 2.*N11 * parts[pp].r * chinm_im[stride*pp + 2];
    parts[pp].Lz = rho_f * vol * parts[pp].r*parts[pp].r * parts[pp].ozdot
      + 8. * PI * mu * nu * N10 * parts[pp].r * chinm_re[stride*pp + 1];
  }
}

__global__ void compute_error(real lamb_cut, int stride, int nparts,
  real *pnm_re, real *pnm_re0, real *pnm_im, real *pnm_im0,
  real *phinm_re, real *phinm_re0, real *phinm_im, real *phinm_im0,
  real *chinm_re, real *chinm_re0, real *chinm_im, real *chinm_im0,
  real *coeffs, real *errors, real *part_errors, dom_struct *dom, real nu)
{
  int part = blockIdx.x;
  int i,j;
  real tmp = FLT_MIN;
  int loc = 0;
  real avg = 0;
  real div = 0;

  // create shared memory space
  __shared__ real s_coeffs[6*21];  // ** have to hard-code this length **
  __shared__ real s_coeffs0[6*21]; // ** have to hard-code this length **
                            // using 6 coefficient sets, each holding
                            // a maximum of 21 coefficients (5th-order
                            // truncation)

  // copy coeffs for this particle into shared memory
  for(i = 0; i < stride; i++) {
    s_coeffs[i] = pnm_re[part*stride+i];
    s_coeffs[i+1*stride] = pnm_im[part*stride+i];
    s_coeffs[i+2*stride] = phinm_re[part*stride+i];
    s_coeffs[i+3*stride] = phinm_im[part*stride+i];
    s_coeffs[i+4*stride] = chinm_re[part*stride+i];
    s_coeffs[i+5*stride] = chinm_im[part*stride+i];
    s_coeffs0[i] = pnm_re0[part*stride+i];
    s_coeffs0[i+1*stride] = pnm_im0[part*stride+i];
    s_coeffs0[i+2*stride] = phinm_re0[part*stride+i];
    s_coeffs0[i+3*stride] = phinm_im0[part*stride+i];
    s_coeffs0[i+4*stride] = chinm_re0[part*stride+i];
    s_coeffs0[i+5*stride] = chinm_im0[part*stride+i];
  }

  // compute the average of the coefficients
  for(i = 0; i < stride*6; i++) {
    avg += s_coeffs[i]*s_coeffs[i];
  }
  avg = avg / (stride*6.);

  // sort the coefficients in shared memory and calculate errors along the way
  for(i = 0; i < 6*stride; i++) {
    // search for the largest magnitude value in shared and store its location
    tmp = FLT_MIN;
    for(j = 0; j < 6*stride; j++) {
      if(s_coeffs[j]*s_coeffs[j] > tmp) {
        tmp = s_coeffs[j]*s_coeffs[j];
        loc = j;
      }
    }

    // move the largest value into sorted list
    coeffs[part*stride+i] = s_coeffs[loc];

    // if its corresponding coefficient has large enough magnitude,
    // compute error for this coefficient
    if(fabs(s_coeffs[loc]) > lamb_cut*fabs(coeffs[part*stride+0])) {
      div = fabs(s_coeffs[loc]);// + fabs(avg)*1e-4;
      if(div < 1e-16) div = 1e-16;
      errors[part*stride+i] = fabs((s_coeffs[loc] - s_coeffs0[loc]) / div);
    } else errors[part*stride+i] = 0.;

    // discard this value since we've used it once
    s_coeffs[loc] = 0.;
  }

  // find the largest error for each particle
  tmp = FLT_MIN;
  for(i = 0; i < 6*stride; i++) {
    if(errors[part*stride+i] > tmp) tmp = errors[part*stride+i];
  }

  // write error to return for each particle
  part_errors[part] = tmp;
}
