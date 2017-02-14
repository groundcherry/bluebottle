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

#include "particle.h"

int *phase;
int **_phase;
int *phase_shell;
int **_phase_shell;
int *flag_u;
int **_flag_u;
int *flag_v;
int **_flag_v;
int *flag_w;
int **_flag_w;
int nparts;
real interactionLengthRatio;
part_struct *parts;
part_struct **_parts;
dom_struct binDom;
dom_struct *_binDom;
int coeff_stride;
real *pnm_re;
real *pnm_im;
real *phinm_re;
real *phinm_im;
real *chinm_re;
real *chinm_im;
real *pnm_re0;
real *pnm_im0;
real *phinm_re0;
real *phinm_im0;
real *chinm_re0;
real *chinm_im0;
real *pnm_re00;
real *pnm_im00;
real *phinm_re00;
real *phinm_im00;
real *chinm_re00;
real *chinm_im00;
real **_pnm_re;
real **_pnm_im;
real **_phinm_re;
real **_phinm_im;
real **_chinm_re;
real **_chinm_im;
real **_pnm_re0;
real **_pnm_im0;
real **_phinm_re0;
real **_phinm_im0;
real **_chinm_re0;
real **_chinm_im0;
real **_pnm_re00;
real **_pnm_im00;
real **_phinm_re00;
real **_phinm_im00;
real **_chinm_re00;
real **_chinm_im00;

void parts_read_input(int turb)
{
  int i;  // iterator

  int fret = 0;
  fret = fret; // prevent compiler warning

  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/part.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // read particle list
  fret = fscanf(infile, "n %d\n", &nparts);
  if(turb) nparts = 0; // remove particles from turbulence precursor simulation
  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  // allocate bin domain
  // read nbody parameters
  #ifdef DOUBLE  
    fret = fscanf(infile, "(l/a) %lf\n", &interactionLengthRatio);
  #else
    fret = fscanf(infile, "(l/a) %f\n", &interactionLengthRatio);
  #endif

  // read nparts particles
  for(i = 0; i < nparts; i++) {
    fret = fscanf(infile, "\n");
#ifdef DOUBLE
    fret = fscanf(infile, "r %lf\n", &parts[i].r);
    fret = fscanf(infile, "(x, y, z) %lf %lf %lf\n",
      &parts[i].x, &parts[i].y, &parts[i].z);
    fret = fscanf(infile, "(aFx, aFy, aFz) %lf %lf %lf\n",
      &parts[i].aFx, &parts[i].aFy, &parts[i].aFz);
    fret = fscanf(infile, "(aLx, aLy, aLz) %lf %lf %lf\n",
      &parts[i].aLx, &parts[i].aLy, &parts[i].aLz);
    fret = fscanf(infile, "rho %lf\n", &parts[i].rho);
    fret = fscanf(infile, "E %lf\n", &parts[i].E);
    fret = fscanf(infile, "sigma %lf\n", &parts[i].sigma);
    fret = fscanf(infile, "e_dry %lf\n", &parts[i].e_dry);
    fret = fscanf(infile, "coeff_fric %lf\n", &parts[i].coeff_fric);
#else // single precision
    fret = fscanf(infile, "r %f\n", &parts[i].r);
    fret = fscanf(infile, "(x, y, z) %f %f %f\n",
      &parts[i].x, &parts[i].y, &parts[i].z);
    fret = fscanf(infile, "(aFx, aFy, aFz) %f %f %f\n",
      &parts[i].aFx, &parts[i].aFy, &parts[i].aFz);
    fret = fscanf(infile, "(aLx, aLy, aLz) %f %f %f\n",
      &parts[i].aLx, &parts[i].aLy, &parts[i].aLz);
    fret = fscanf(infile, "rho %f\n", &parts[i].rho);
    fret = fscanf(infile, "E %f\n", &parts[i].E);
    fret = fscanf(infile, "sigma %f\n", &parts[i].sigma);
    fret = fscanf(infile, "e_dry %f\n", &parts[i].e_dry);
    fret = fscanf(infile, "coeff_fric %f\n", &parts[i].coeff_fric);
#endif
    fret = fscanf(infile, "order %d\n", &parts[i].order);
#ifdef DOUBLE
    fret = fscanf(infile, "rs/r %lf\n", &parts[i].rs);
    fret = fscanf(infile, "spring_k %lf\n", &parts[i].spring_k);
    fret = fscanf(infile, "spring (x, y, z) %lf %lf %lf\n",
      &parts[i].spring_x, &parts[i].spring_y, &parts[i].spring_z);
    fret = fscanf(infile, "spring_l %lf\n", &parts[i].spring_l);
#else // single precision
    fret = fscanf(infile, "rs/r %f\n", &parts[i].rs);
    fret = fscanf(infile, "spring_k %f\n", &parts[i].spring_k);
    fret = fscanf(infile, "spring (x, y, z) %f %f %f\n",
      &parts[i].spring_x, &parts[i].spring_y, &parts[i].spring_z);
    fret = fscanf(infile, "spring_l %f\n", &parts[i].spring_l);
#endif
    fret = fscanf(infile, "translating %d\n", &parts[i].translating);
    fret = fscanf(infile, "rotating %d\n", &parts[i].rotating);
  }

  fclose(infile);
}

void parts_show_config(void)
{
  int i;  // iterator

  printf("Particles:\n");
  #ifdef DOUBLE  
    printf("  Interaction Support Length Ratio (l/a) = %lf\n", 
      interactionLengthRatio);
  #else
    printf("  Interaction Support Length Ratio (l/a) = %f\n", 
      interactionLengthRatio);
  #endif
  for(i = 0; i < nparts; i++) {
    printf("  Particle %d:\n", i);
    printf("    r = %e\n", parts[i].r);
    printf("    (x, y, z) = (%e, %e, %e)\n",
      parts[i].x, parts[i].y, parts[i].z);
    printf("    (u, v, w) = (%e, %e, %e)\n",
      parts[i].u, parts[i].v, parts[i].w);
    printf("    (udot, vdot, wdot) = (%e, %e, %e)\n",
      parts[i].udot, parts[i].vdot, parts[i].wdot);
    printf("    (axx, axy, axz) = (%e, %e, %e)\n",
      parts[i].axx, parts[i].axy, parts[i].axz);
    printf("    (ayx, ayy, ayz) = (%e, %e, %e)\n",
      parts[i].ayx, parts[i].ayy, parts[i].ayz);
    printf("    (azx, azy, azz) = (%e, %e, %e)\n",
      parts[i].azx, parts[i].azy, parts[i].azz);
    printf("    (ox, oy, oz) = (%e, %e, %e)\n",
      parts[i].ox, parts[i].oy, parts[i].oz);
    printf("    (oxdot, oydot, ozdot) = (%e, %e, %e)\n",
      parts[i].oxdot, parts[i].oydot, parts[i].ozdot);
    printf("    (Fx, Fy, Fz) = (%e %e %e)\n",
      parts[i].Fx, parts[i].Fy, parts[i].Fz);
    printf("    (Lx, Ly, Lz) = (%e %e %e)\n",
      parts[i].Lx, parts[i].Ly, parts[i].Lz);
    printf("    (aFx, aFy, aFz) = (%e %e %e)\n",
      parts[i].aFx, parts[i].aFy, parts[i].aFz);
    printf("    (aLx, aLy, aLz) = (%e %e %e)\n",
      parts[i].aLx, parts[i].aLy, parts[i].aLz);
    printf("    Cage:\n");
    printf("      cx = %d, cy = %d, cz = %d\n",
      parts[i].cage.cx, parts[i].cage.cy, parts[i].cage.cz);
    printf("      is = %d, ibs = %d, ibe = %d, ie = %d, in = %d\n",
      parts[i].cage.is, parts[i].cage.ibs, parts[i].cage.ibe,
      parts[i].cage.ie, parts[i].cage.in);
    printf("      js = %d, jbs = %d, jbe = %d, je = %d, jn = %d\n",
      parts[i].cage.js, parts[i].cage.jbs, parts[i].cage.jbe,
      parts[i].cage.je, parts[i].cage.jn);
    printf("      ks = %d, kbs = %d, kbe = %d, ke = %d, kn = %d\n",
      parts[i].cage.ks, parts[i].cage.kbs, parts[i].cage.kbe,
      parts[i].cage.ke, parts[i].cage.kn);
    printf("    rho = %f\n", parts[i].rho);
    printf("    E = %e\n", parts[i].E);
    printf("    sigma = %e\n", parts[i].sigma);
    printf("    e_dry = %e\n", parts[i].e_dry);
    printf("    coeff_fric = %e\n", parts[i].coeff_fric);
    printf("    order = %d\n", parts[i].order);
    printf("    rs = %e\n", parts[i].rs);
    printf("    spring_k = %f\n", parts[i].spring_k);
    printf("    spring (x, y, z) = (%e %e %e)\n", parts[i].spring_x,
      parts[i].spring_y, parts[i].spring_z);
    printf("    spring_l = %f\n", parts[i].spring_l);
    printf("    ncoeff = %d\n", parts[i].ncoeff);
    printf("    translating = %d\n", parts[i].translating);
    printf("    rotating = %d\n", parts[i].rotating);
  }
}

void bin_show_config(void) {

  printf(" Bin Domain:\n");
  printf("  X: (%f, %f), dX = %f\n", binDom.xs, binDom.xe, binDom.dx);
  printf("  Y: (%f, %f), dY = %f\n", binDom.ys, binDom.ye, binDom.dy);
  printf("  Z: (%f, %f), dZ = %f\n", binDom.zs, binDom.ze, binDom.dz);
  printf("  Xn = %d, Yn = %d, Zn = %d\n", binDom.xn, binDom.yn, binDom.zn);
  printf("binDomain Grids:\n");
  printf("  binDom.Gcc:\n");
  printf("    is = %d, ie = %d, in = %d\n", binDom.Gcc.is, binDom.Gcc.ie, binDom.Gcc.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", binDom.Gcc.isb, binDom.Gcc.ieb,
    binDom.Gcc.inb);
  printf("    js = %d, je = %d, jn = %d\n", binDom.Gcc.js, binDom.Gcc.je, binDom.Gcc.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", binDom.Gcc.jsb, binDom.Gcc.jeb,
    binDom.Gcc.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", binDom.Gcc.ks, binDom.Gcc.ke, binDom.Gcc.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", binDom.Gcc.ksb, binDom.Gcc.keb,
    binDom.Gcc.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", binDom.Gcc.s1, binDom.Gcc.s2,
    binDom.Gcc.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", binDom.Gcc.s1b, binDom.Gcc.s2b,
    binDom.Gcc.s3b);
  printf("  binDom.Gfx:\n");
  printf("    is = %d, ie = %d, in = %d\n", binDom.Gfx.is, binDom.Gfx.ie, binDom.Gfx.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", binDom.Gfx.isb, binDom.Gfx.ieb,
    binDom.Gfx.inb);
  printf("    js = %d, je = %d, jn = %d\n", binDom.Gfx.js, binDom.Gfx.je, binDom.Gfx.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", binDom.Gfx.jsb, binDom.Gfx.jeb,
    binDom.Gfx.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", binDom.Gfx.ks, binDom.Gfx.ke, binDom.Gfx.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", binDom.Gfx.ksb, binDom.Gfx.keb,
    binDom.Gfx.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", binDom.Gfx.s1, binDom.Gfx.s2,
    binDom.Gfx.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", binDom.Gfx.s1b, binDom.Gfx.s2b,
    binDom.Gfx.s3b);
  printf("  binDom.Gfy:\n");
  printf("    is = %d, ie = %d, in = %d\n", binDom.Gfy.is, binDom.Gfy.ie, binDom.Gfy.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", binDom.Gfy.isb, binDom.Gfy.ieb,
    binDom.Gfy.inb);
  printf("    js = %d, je = %d, jn = %d\n", binDom.Gfy.js, binDom.Gfy.je, binDom.Gfy.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", binDom.Gfy.jsb, binDom.Gfy.jeb,
    binDom.Gfy.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", binDom.Gfy.ks, binDom.Gfy.ke, binDom.Gfy.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", binDom.Gfy.ksb, binDom.Gfy.keb,
    binDom.Gfy.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", binDom.Gfy.s1, binDom.Gfy.s2,
    binDom.Gfy.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", binDom.Gfy.s1b, binDom.Gfy.s2b,
    binDom.Gfy.s3b);
  printf("  binDom.Gfz:\n");
  printf("    is = %d, ie = %d, in = %d\n", binDom.Gfz.is, binDom.Gfz.ie, binDom.Gfz.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", binDom.Gfz.isb, binDom.Gfz.ieb,
    binDom.Gfz.inb);
  printf("    js = %d, je = %d, jn = %d\n", binDom.Gfz.js, binDom.Gfz.je, binDom.Gfz.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", binDom.Gfz.jsb, binDom.Gfz.jeb,
    binDom.Gfz.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", binDom.Gfz.ks, binDom.Gfz.ke, binDom.Gfz.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", binDom.Gfz.ksb, binDom.Gfz.keb,
    binDom.Gfz.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", binDom.Gfz.s1, binDom.Gfz.s2,
    binDom.Gfz.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", binDom.Gfz.s1b, binDom.Gfz.s2b,
    binDom.Gfz.s3b);
}

int parts_init(void)
{
  int i, j;  // iterators

  // allocate and initialize phase array
  phase = (int*) malloc(Dom.Gcc.s3b * sizeof(int));
  cpumem += Dom.Gcc.s3b * sizeof(int);
  phase_shell = (int*) malloc(Dom.Gcc.s3b * sizeof(int));
  cpumem += Dom.Gcc.s3b * sizeof(int);
  flag_u = (int*) malloc(Dom.Gfx.s3b * sizeof(int));
  cpumem += Dom.Gfx.s3b * sizeof(int);
  flag_v = (int*) malloc(Dom.Gfy.s3b * sizeof(int));
  cpumem += Dom.Gfy.s3b * sizeof(int);
  flag_w = (int*) malloc(Dom.Gfz.s3b * sizeof(int));
  cpumem += Dom.Gfz.s3b * sizeof(int);
  flags_reset();

  // build the particle cages
  for(i = 0; i < nparts; i++) init_cage(i);

  // TODO: check that the particles are in valid locations (i.e. exist
  // completely inside the domain and do not overlap)

  coeff_stride = 0;

  for(i = 0; i < nparts; i++) {
    parts[i].rs = parts[i].rs * parts[i].r;
    // set rs as one cell away from surface of particle
    //parts[i].rs = parts[i].r + 2.*(Dom.dx + Dom.dy + Dom.dz)/3.;

    // calculate the number of coefficients needed
    parts[i].ncoeff = 0;
    for(j = 0; j <= parts[i].order; j++) {
      parts[i].ncoeff += j + 1;
    }

    // initialize previous position
    parts[i].x0 = parts[i].x;
    parts[i].y0 = parts[i].y;
    parts[i].z0 = parts[i].z;

    // initialize velocity and acceleration to zero (default: QUIESCENT)
    parts[i].u = 0.;
    parts[i].v = 0.;
    parts[i].w = 0.;
    parts[i].u0 = 0.;
    parts[i].v0 = 0.;
    parts[i].w0 = 0.;
    parts[i].udot = 0.;
    parts[i].vdot = 0.;
    parts[i].wdot = 0.;
    parts[i].udot0 = 0.;
    parts[i].vdot0 = 0.;
    parts[i].wdot0 = 0.;
    /* set initial position of particle reference basis to match the global
     * domain basis */
    parts[i].axx = 1.;
    parts[i].axy = 0.;
    parts[i].axz = 0.;
    parts[i].ayx = 0.;
    parts[i].ayy = 1.;
    parts[i].ayz = 0.;
    parts[i].azx = 0.;
    parts[i].azy = 0.;
    parts[i].azz = 1.;
    parts[i].ox = 0.;
    parts[i].oy = 0.;
    parts[i].oz = 0.;
    parts[i].ox0 = 0.;
    parts[i].oy0 = 0.;
    parts[i].oz0 = 0.;
    parts[i].oxdot = 0.;
    parts[i].oydot = 0.;
    parts[i].ozdot = 0.;
    parts[i].oxdot0 = 0.;
    parts[i].oydot0 = 0.;
    parts[i].ozdot0 = 0.;

    if(init_cond == SHEAR) {
      // initialize SHEAR flow
      if(parts[i].translating) { // if translating
        // set linear velocity according to position
        parts[i].u = (bc.uNDm-bc.uSDm)*(parts[i].y-Dom.ys)/Dom.yl + bc.uSDm;
        parts[i].u += (bc.uTDm-bc.uBDm)*(parts[i].z-Dom.zs)/Dom.zl + bc.uBDm;
        parts[i].u0 = parts[i].u;
        parts[i].v = (bc.vEDm-bc.vWDm)*(parts[i].x-Dom.xs)/Dom.xl + bc.vWDm;
        parts[i].v += (bc.vTDm-bc.vBDm)*(parts[i].z-Dom.zs)/Dom.zl + bc.vBDm;
        parts[i].v0 = parts[i].v;
        parts[i].w = (bc.wEDm-bc.wWDm)*(parts[i].x-Dom.xs)/Dom.xl + bc.wWDm;
        parts[i].w += (bc.wNDm-bc.wSDm)*(parts[i].y-Dom.ys)/Dom.yl + bc.wSDm;
        parts[i].w0 = parts[i].w;

        // initialize previous position
        parts[i].x0 = parts[i].x - dt * parts[i].u;
        parts[i].y0 = parts[i].y - dt * parts[i].v;
        parts[i].z0 = parts[i].z - dt * parts[i].w;
      }
      if(parts[i].rotating) { // if rotating
        // set angular velocity according to (one-half of) the shear rate
        parts[i].ox = -0.5*(bc.vTDm-bc.vBDm)/Dom.zl;
        parts[i].ox += 0.5*(bc.wNDm-bc.wSDm)/Dom.yl;
        parts[i].oy = 0.5*(bc.uTDm-bc.uBDm)/Dom.zl;
        parts[i].oy += -0.5*(bc.wEDm-bc.wWDm)/Dom.xl;
        parts[i].oz = -0.5*(bc.uNDm-bc.uSDm)/Dom.yl;
        parts[i].oz += 0.5*(bc.vEDm-bc.vWDm)/Dom.xl;
        parts[i].ox0 = parts[i].ox;
        parts[i].oy0 = parts[i].oy;
        parts[i].oz0 = parts[i].oz;
      }
    } else if(init_cond == CHANNEL) {
      // initialize CHANNEL flow
      if(parts[i].translating) { // if translating
        // set linear velocity according to position
        real x = parts[i].x;
        real y = parts[i].y;
        real z = parts[i].z;
        parts[i].u = 0.5/mu*gradP.xm*(y*y-(Dom.ys+Dom.ye)*y+Dom.ys*Dom.ye)
          * (bc.uS == DIRICHLET)
          + 0.5/mu*gradP.xm*(z*z-(Dom.zs+Dom.ze)*z+Dom.zs*Dom.ze)
          * (bc.uB == DIRICHLET);
        parts[i].u0 = parts[i].u;
        parts[i].v = 0.5/mu*gradP.ym*(x*x-(Dom.xs+Dom.xe)*x+Dom.xs*Dom.xe)
          * (bc.vW == DIRICHLET)
          + 0.5/mu*gradP.ym*(z*z-(Dom.zs+Dom.ze)*z+Dom.zs*Dom.ze)
          * (bc.vB == DIRICHLET);
        parts[i].v0 = parts[i].v;
        parts[i].w = 0.5/mu*gradP.zm*(x*x-(Dom.xs+Dom.xe)*x+Dom.xs*Dom.xe)
          * (bc.wW == DIRICHLET)
          + 0.5/mu*gradP.zm*(y*y-(Dom.ys+Dom.ye)*y+Dom.ys*Dom.ye)
          * (bc.wS == DIRICHLET);
        parts[i].w0 = parts[i].w;

        // initialize previous position
        parts[i].x0 = parts[i].x - dt * parts[i].u;
        parts[i].y0 = parts[i].y - dt * parts[i].v;
        parts[i].z0 = parts[i].z - dt * parts[i].w;
      }


      // initialize no rotation component
    }

    // initialize the hydrodynamic forces and moments to zero
    parts[i].Fx = 0.;
    parts[i].Fy = 0.;
    parts[i].Fz = 0.;
    parts[i].Lx = 0.;
    parts[i].Ly = 0.;
    parts[i].Lz = 0.;

    // initialize the particle spring force to zero
    parts[i].kFx = 0.;
    parts[i].kFy = 0.;
    parts[i].kFz = 0.;

    // initialize the particle interaction force to zero
    parts[i].iFx = 0.;
    parts[i].iFy = 0.;
    parts[i].iFz = 0.;
    parts[i].iLx = 0.;
    parts[i].iLy = 0.;
    parts[i].iLz = 0.;

    // set the Lamb's coefficient matrix stride length equal to the
    // largest number of coefficients any particle needs to hold
    if(parts[i].ncoeff > coeff_stride) coeff_stride = parts[i].ncoeff;

    // initialize nodes array
    for(j = 0; j < NNODES; j++) {
      parts[i].nodes[j] = -1;
    }

    // initialize Stokes number lists to -1
    for(j = 0; j < MAX_NEIGHBORS; j++) {
      parts[i].St[j] = 0.;
      parts[i].iSt[j] = -1;
    }
  }

  // allocate Lamb's coefficients
  pnm_re = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  pnm_im = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  phinm_re = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  phinm_im = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  chinm_re = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  chinm_im = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  pnm_re0 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  pnm_im0 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  phinm_re0 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  phinm_im0 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  chinm_re0 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  chinm_im0 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  pnm_re00 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  pnm_im00 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  phinm_re00 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  phinm_im00 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  chinm_re00 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);
  chinm_im00 = (real*) malloc(coeff_stride * nparts * sizeof(real));
  cpumem += coeff_stride * nparts * sizeof(real);

  // initialize Lamb's coefficients
  for(i = 0; i < coeff_stride * nparts; i++) {
    pnm_re[i] = 0.;
    pnm_im[i] = 0.;
    phinm_re[i] = 0.;
    phinm_im[i] = 0.;
    chinm_re[i] = 0.;
    chinm_im[i] = 0.;
    pnm_re0[i] = 0.;
    pnm_im0[i] = 0.;
    phinm_re0[i] = 0.;
    phinm_im0[i] = 0.;
    chinm_re0[i] = 0.;
    chinm_im0[i] = 0.;
    pnm_re00[i] = 0.;
    pnm_im00[i] = 0.;
    phinm_re00[i] = 0.;
    phinm_im00[i] = 0.;
    chinm_re00[i] = 0.;
    chinm_im00[i] = 0.;
  }

  return EXIT_SUCCESS;
}

void init_cage(int part)
{
  parts[part].cage.cx = 0;
  parts[part].cage.cy = 0;
  parts[part].cage.cz = 0;
  parts[part].cage.is = 0;
  parts[part].cage.ibs = 0;
  parts[part].cage.ie = 0;
  parts[part].cage.ibe = 0;
  parts[part].cage.in = 0;
  parts[part].cage.js = 0;
  parts[part].cage.jbs = 0;
  parts[part].cage.je = 0;
  parts[part].cage.jbe = 0;
  parts[part].cage.jn = 0;
  parts[part].cage.ks = 0;
  parts[part].cage.kbs = 0;
  parts[part].cage.ke = 0;
  parts[part].cage.kbe = 0;
  parts[part].cage.kn = 0;
}

void flags_reset(void)
{
  int i, j, k;  // iterators
  int C;        // cell location
  for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        phase[C] = 1;  // all fluid
        phase_shell[C] = 1;  // all fluid
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        flag_u[C] = 1;  // not a boundary
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        flag_v[C] = 1;  // not a boundary
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        flag_w[C] = 1;  // not a boundary
      }
    }
  }
}

int binDom_init(void)
{
  // find max of radii to set up bins
  real rmax = 0;
  for (int i = 0; i < nparts; i++) {
    rmax = rmax + (parts[i].r > rmax)*(parts[i].r - rmax);
  }
  binDom.xs = Dom.xs;
  binDom.xe = Dom.xe;
  binDom.xl = Dom.xl;
  binDom.xn = floor(Dom.xl/(2.*rmax + interactionLengthRatio*rmax));
  if (binDom.xn == 0) { // to avoid dividing by zero and having infinite bin
    binDom.xn = 1;
    binDom.dx = Dom.xl;
  } else {
    binDom.dx = Dom.xl / binDom.xn;
  }

  binDom.ys = Dom.ys;
  binDom.ye = Dom.ye;
  binDom.yl = Dom.yl;
  binDom.yn = floor(Dom.yl/(2.*rmax + interactionLengthRatio*rmax));
  if (binDom.yn == 0) {
    binDom.yn = 1;
    binDom.dy = Dom.yl;
  } else {
    binDom.dy = Dom.yl / binDom.yn;
  }

  binDom.zs = Dom.zs;
  binDom.ze = Dom.ze;
  binDom.zl = Dom.zl;
  binDom.zn = floor(Dom.zl/(2.*rmax + interactionLengthRatio*rmax));
  if (binDom.zn == 0) {
    binDom.zn = 1;
    binDom.dz = Dom.zl;
  } else {
    binDom.dz = Dom.zl / binDom.zn;
  }

  binDom.E = 0;
  binDom.W = 0;
  binDom.N = 0;
  binDom.S = 0;
  binDom.T = 0;
  binDom.B = 0;

  //GCC
  binDom.Gcc.is = DOM_BUF;
  binDom.Gcc.isb = 0;
  binDom.Gcc.in = binDom.xn;
  binDom.Gcc.inb = 0;
  binDom.Gcc.ie = binDom.Gcc.is + binDom.Gcc.in;
  binDom.Gcc.ieb = 0;

  binDom.Gcc.js = DOM_BUF;
  binDom.Gcc.jsb = 0;
  binDom.Gcc.jn = binDom.yn;
  binDom.Gcc.jnb = 0;
  binDom.Gcc.je = binDom.Gcc.js + binDom.Gcc.jn;
  binDom.Gcc.jeb = 0;

  binDom.Gcc.ks = DOM_BUF;
  binDom.Gcc.ksb = 0;
  binDom.Gcc.kn = binDom.zn;
  binDom.Gcc.knb = 0;
  binDom.Gcc.ke = DOM_BUF + binDom.Gcc.kn;
  binDom.Gcc.keb = 0;

  binDom.Gcc.s1 = binDom.Gcc.in;
  binDom.Gcc.s2 = binDom.Gcc.s1 * binDom.Gcc.jn;
  binDom.Gcc.s3 = binDom.Gcc.s2 * binDom.Gcc.kn;
  binDom.Gcc.s1b = 0;
  binDom.Gcc.s2b = 0;
  binDom.Gcc.s3b = 0;

  //GFX
  binDom.Gfx.is = 0;
  binDom.Gfx.in = 0;
  binDom.Gfx.ie = 0;
  binDom.Gfx.isb = 0;
  binDom.Gfx.inb = 0;
  binDom.Gfx.ieb = 0;

  binDom.Gfx.js = 0;
  binDom.Gfx.jn = 0;
  binDom.Gfx.je = 0;
  binDom.Gfx.jsb = 0;
  binDom.Gfx.jnb = 0;
  binDom.Gfx.jeb = 0;
  
  binDom.Gfx.ks = 0;
  binDom.Gfx.kn = 0;
  binDom.Gfx.ke = 0;
  binDom.Gfx.ksb = 0;
  binDom.Gfx.knb = 0;
  binDom.Gfx.keb = 0;

  binDom.Gfx.s1 = 0;
  binDom.Gfx.s1b = 0;
  binDom.Gfx.s2 = 0;
  binDom.Gfx.s2b = 0;
  binDom.Gfx.s3 = 0;
  binDom.Gfx.s3b = 0;

  //GFY
  binDom.Gfy.is = 0;
  binDom.Gfy.in = 0;
  binDom.Gfy.ie = 0;
  binDom.Gfy.isb = 0;
  binDom.Gfy.inb = 0;
  binDom.Gfy.ieb = 0;

  binDom.Gfy.js = 0;
  binDom.Gfy.jn = 0;
  binDom.Gfy.je = 0;
  binDom.Gfy.jsb = 0;
  binDom.Gfy.jnb = 0;
  binDom.Gfy.jeb = 0;
  
  binDom.Gfy.ks = 0;
  binDom.Gfy.kn = 0;
  binDom.Gfy.ke = 0;
  binDom.Gfy.ksb = 0;
  binDom.Gfy.knb = 0;
  binDom.Gfy.keb = 0;

  binDom.Gfy.s1 = 0;
  binDom.Gfy.s1b = 0;
  binDom.Gfy.s2 = 0;
  binDom.Gfy.s2b = 0;
  binDom.Gfy.s3 = 0;
  binDom.Gfy.s3b = 0;

  //GFZ
  binDom.Gfz.is = 0;
  binDom.Gfz.in = 0;
  binDom.Gfz.ie = 0;
  binDom.Gfz.isb = 0;
  binDom.Gfz.inb = 0;
  binDom.Gfz.ieb = 0;

  binDom.Gfz.js = 0;
  binDom.Gfz.jn = 0;
  binDom.Gfz.je = 0;
  binDom.Gfz.jsb = 0;
  binDom.Gfz.jnb = 0;
  binDom.Gfz.jeb = 0;
  
  binDom.Gfz.ks = 0;
  binDom.Gfz.kn = 0;
  binDom.Gfz.ke = 0;
  binDom.Gfz.ksb = 0;
  binDom.Gfz.knb = 0;
  binDom.Gfz.keb = 0;

  binDom.Gfz.s1 = 0;
  binDom.Gfz.s1b = 0;
  binDom.Gfz.s2 = 0;
  binDom.Gfz.s2b = 0;
  binDom.Gfz.s3 = 0;
  binDom.Gfz.s3b = 0;

  return EXIT_SUCCESS;
}

void parts_clean(void)
{
  // free Lamb's coefficients
  free(pnm_re);
  free(pnm_im);
  free(phinm_re);
  free(phinm_im);
  free(chinm_re);
  free(chinm_im);
  free(pnm_re0);
  free(pnm_im0);
  free(phinm_re0);
  free(phinm_im0);
  free(chinm_re0);
  free(chinm_im0);
  free(pnm_re00);
  free(pnm_im00);
  free(phinm_re00);
  free(phinm_im00);
  free(chinm_re00);
  free(chinm_im00);

  free(parts);
  free(phase);
  free(phase_shell);
  free(flag_u);
  free(flag_v);
  free(flag_w);
}
