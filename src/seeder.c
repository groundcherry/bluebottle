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

#include <time.h>

#include "bluebottle.h"
#include "domain.h"
#include "particle.h"

void seeder_read_input(int Nx, int Ny, int Nz, double ddz, double bias, 
  int nperturb)
{

  int N;               // number of parts
  real loa;            // interaction length

  real a;              // particle radius
  real x,y,z;          // particle positions
  real aFx, aFy, aFz;  // particle linear forcing
  real aLx, aLy, aLz;  // particle angular forcing
  real rho;            // density
  real E;              // youngs modulus
  real sigma;          // poisson ratio
  real e_dry;          // dry coefficient of restitution
  int order;           // lamb truncation order
  real rs_r;           // cage ratio extents
  real spring_k;       // particle spring constant
  real spring_x;       // spring attachment locations
  real spring_y;
  real spring_z;
  real spring_l;       // spring length
  int trans;           // particle is allowed to translate
  int rot;             // particle is allowed to rotate

  int fret = 0;
  fret = fret;    // prevent compiler warning

  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/part.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }


  // read particle list
  fret = fscanf(infile, "n %d\n", &N);

#ifdef DOUBLE
  fret = fscanf(infile, "(l/a) %lf\n", &loa);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "r %lf\n", &a);
  fret = fscanf(infile, "(x, y, z) %lf %lf %lf\n", &x, &y, &z);
  fret = fscanf(infile, "(aFx, aFy, aFz) %lf %lf %lf\n",
    &aFx, &aFy, &aFz);
  fret = fscanf(infile, "(aLx, aLy, aLz) %lf %lf %lf\n",
    &aLx, &aLy, &aLz);
  fret = fscanf(infile, "rho %lf\n", &rho);
  fret = fscanf(infile, "E %lf\n", &E);
  fret = fscanf(infile, "sigma %lf\n", &sigma);
  fret = fscanf(infile, "e_dry %lf\n", &e_dry);
#else
  fret = fscanf(infile, "loa %f\n", &loa);
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "r %f\n", &a);
  fret = fscanf(infile, "(x, y, z) %f %f %f\n", &x, &y, &z);
  fret = fscanf(infile, "(aFx, aFy, aFz) %f %f %f\n",
    &aFx, &aFy, &aFz);
  fret = fscanf(infile, "(aLx, aLy, aLz) %f %f %f\n",
    &aLx, &aLy, &aLz);
  fret = fscanf(infile, "rho %f\n", &rho);
  fret = fscanf(infile, "E %f\n", &E);
  fret = fscanf(infile, "sigma %f\n", &sigma);
  fret = fscanf(infile, "e_dry %f\n", &e_dry);
#endif
  fret = fscanf(infile, "order %d\n", &order);
#ifdef DOUBLE
  fret = fscanf(infile, "rs/r %lf\n", &rs_r);
  fret = fscanf(infile, "spring_k %lf\n", &spring_k);
  fret = fscanf(infile, "spring (x, y, z) %lf %lf %lf\n",
    &spring_x, &spring_y, &spring_z);
  fret = fscanf(infile, "spring_l %lf\n", &spring_l);
#else // single precision
  fret = fscanf(infile, "rs/r %f\n", &rs_r);
  fret = fscanf(infile, "spring_k %f\n", &spring_k);
  fret = fscanf(infile, "spring (x, y, z) %f %f %f\n",
    &spring_x, &spring_y, &spring_z);
  fret = fscanf(infile, "spring_l %f\n", &spring_l);
#endif
  fret = fscanf(infile, "translating %d\n", &trans);
  fret = fscanf(infile, "rotating %d\n", &rot);

  // check parameters
  if (N < 1) {
    printf("Error: N must be > 1\n");
    exit(EXIT_FAILURE);
  } else if (a < 0) {
    printf("Error: a must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (rho < 0) {
    printf("Error: rho must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (E < 0) {
    printf("Error: E must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (sigma > 0.5 || sigma <= -1) {
    printf("Error: sigma must be between -1 < sigma <= 0.5\n");
    exit(EXIT_FAILURE);
  } else if (order < 0) {
    printf("Error: order must be > 0\n");
    exit(EXIT_FAILURE);
  } else if (trans > 1) {
    printf("Error: translating must be 0 or 1\n");
    exit(EXIT_FAILURE);
  } else if (rot > 1) {
    printf("Error: rotating must be 0 or 1\n");
    exit(EXIT_FAILURE);
  }

  // SEEDER

  if (Nx*Ny*Nz != 0 ) {
    N = Nx*Ny*Nz;
  }
  printf("Requested Parameters:\n");
  printf("       N = %d\n", N);
#ifdef DOUBLE
  printf("       (l/a) = %lf\n", loa);
  printf("       r = %lf\n", a);
  printf("       (aFx, aFy, aFz) = (%lf, %lf, %lf)\n", aFx, aFy, aFz);
  printf("       (aLx, aLy, aLz) = (%lf, %lf, %lf)\n", aLx, aLy, aLz);
  printf("       rho = %lf\n", rho);
  printf("       E = %lf\n", E);
  printf("       sigma = %lf\n", sigma);
  printf("       e_dry = %lf\n", e_dry);
  printf("       order = %d\n", order);
#else
  printf("       (l/a) = %f\n", loa);
  printf("       r = %f\n", a);
  printf("       (aFx, aFy, aFz) = (%f, %f, %f)\n", aFx, aFy, aFz);
  printf("       (aLx, aLy, aLz) = (%f, %f, %f)\n", aLx, aLy, aLz);
  printf("       rho = %f\n", rho);
  printf("       E = %f\n", E);
  printf("       sigma = %f\n", sigma);
  printf("       e_dry = %f\n", e_dry);
  printf("       order = %d\n", order);
#endif
#ifdef DOUBLE
  printf("       rs_r = %lf\n", rs_r);
  printf("       spring_k = %lf\n", spring_k);
  printf("       spring (x, y, z) = (%lf, %lf, %lf)\n", 
    spring_x, spring_y, spring_z);
  printf("       spring_l = %lf\n", spring_l);
#else
  printf("       rs_r = %f\n", rs_r);
  printf("       spring_k = %f\n", spring_k);
  printf("       spring (x, y, z) = (%f, %f, %f)\n", 
    spring_x, spring_y, spring_z);
  printf("       spring_l = %f\n", spring_l);
#endif
  printf("       translating = %d\n", trans);
  printf("       rotating = %d\n", rot);

  fflush(stdout);

  if(Nx == 0 && Ny == 0 && Nz == 0){ // random case
    seeder(N, loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, sigma, e_dry, 
      order, rs_r, spring_k, spring_x, spring_y, spring_z, spring_l,
      trans, rot);   
  }
  else if(ddz == 0.0 && bias == 0.0 && nperturb == 0){ // array case
    seeder_array(Nx, Ny, Nz, loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, 
    sigma, e_dry, order, rs_r, spring_k, spring_x, spring_y, spring_z, 
    spring_l,trans, rot);      
  }
  else if(ddz != 0.0 && bias == 0.0 && nperturb == 0){ // hex case
    seeder_hex(Nx, Ny, Nz, ddz, loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, 
    sigma, e_dry, order, rs_r, spring_k, spring_x, spring_y, spring_z, 
    spring_l, trans, rot);
  }
  else if(ddz == 0.0 && bias != 0 && nperturb != 0){ // perturb case
    seeder_high_vol_random(Nx, Ny, Nz, bias, nperturb, loa, a, aFx, aFy, aFz, 
    aLx, aLy, aLz, rho, E, sigma, e_dry, order, rs_r, spring_k, 
    spring_x, spring_y, spring_z, spring_l, trans, rot);    
  }
  else{
    printf("The input parameters are not correct! Please check!\n");
    fflush(stdout);
  }
}

void seeder(int nparts, real loa, real a, real aFx, real aFy, real aFz, 
  real aLx, real aLy, real aLz, real rho, real E, real sigma, real e_dry,
  int o, real rs, real spring_k, real spring_x, real spring_y,
  real spring_z, real spring_l, int t, int r) {

  printf("Running bluebottle seeder for %d particles...\n\n", nparts);
  fflush(stdout);
  real xx, yy, zz;
  int fits = 1;
  int attempts = 1;
  int fail = 0;
  int redo = 1;

  // seed the random number generator
  srand(time(NULL));

  // read domain input
  domain_read_input();
  domain_init();

  // domain size accounting for screen
  real xs = Dom.xs + bc.dsW;
  real xe = Dom.xe - bc.dsE;
  real xl = Dom.xl - bc.dsE - bc.dsW;
  real ys = Dom.ys + bc.dsS;
  real ye = Dom.ye - bc.dsN;
  real yl = Dom.yl - bc.dsN - bc.dsS;
  real zs = Dom.zs + bc.dsB;
  real ze = Dom.ze - bc.dsT;
  real zl = Dom.zl - bc.dsT - bc.dsB;
  
  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  real gap = 1.00;

  // place the first particle
  parts[0].r = a;
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].x = rand() / (real)RAND_MAX;
    parts[0].x *= xl;
    parts[0].x += xs;
    if((bc.uW != PERIODIC) && (parts[0].x < (xs + gap*parts[0].r)))
      redo = 1;
    if((bc.uE != PERIODIC) && (parts[0].x > (xe - gap*parts[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].y = rand() / (real)RAND_MAX;
    parts[0].y *= yl;
    parts[0].y += ys;
    if((bc.vS != PERIODIC) && (parts[0].y < (ys + gap*parts[0].r)))
      redo = 1;
    if((bc.vN != PERIODIC) && (parts[0].y > (ye - gap*parts[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].z = rand() / (real)RAND_MAX;
    //parts[0].z = acos(2.*(parts[0].z-0.5))/PI;
    parts[0].z *= zl;
    parts[0].z += zs;
    if((bc.wB != PERIODIC) && (parts[0].z < (zs + gap*parts[0].r)))
      redo = 1;
    if((bc.wT != PERIODIC) && (parts[0].z > (ze - gap*parts[0].r)))
      redo = 1;
  }

  parts[0].u = 0;
  parts[0].v = 0;
  parts[0].w = 0;
  parts[0].aFx = aFx;
  parts[0].aFy = aFy;
  parts[0].aFz = aFz;
  parts[0].aLx = aLx;
  parts[0].aLy = aLy;
  parts[0].aLz = aLz;
  parts[0].rho = rho;
  parts[0].E = E;
  parts[0].sigma = sigma;
  parts[0].e_dry = e_dry;
  parts[0].order = o;
  parts[0].rs = rs;
  parts[0].spring_k = spring_k;
  parts[0].spring_x = spring_x;
  parts[0].spring_y = spring_y;
  parts[0].spring_z = spring_z;
  parts[0].spring_l = spring_l;
  parts[0].ncoeff = 0;
  parts[0].translating = t;
  parts[0].rotating = r;

  // place the rest of the particles
  int i = 0;
  for(i = 1; i < nparts; i++) {
    fits = !fits;
    if(fail) break;
    while(!fits) {
      attempts++;
      // place particle
      parts[i].r = a;
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].x = rand() / (real)RAND_MAX;
        parts[i].x *= xl;
        parts[i].x += xs;
        if((bc.uW != PERIODIC) && (parts[i].x < (xs + gap*parts[i].r)))
          redo = 1;
        if((bc.uE != PERIODIC) && (parts[i].x > (xe - gap*parts[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].y = rand() / (real)RAND_MAX;
        parts[i].y *= yl;
        parts[i].y += ys;
        if((bc.vS != PERIODIC) && (parts[i].y < (ys + gap*parts[i].r)))
          redo = 1;
        if((bc.vN != PERIODIC) && (parts[i].y > (ye - gap*parts[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].z = rand() / (real)RAND_MAX;
        parts[i].z *= zl;
        parts[i].z += zs;
        if((bc.wB != PERIODIC) && (parts[i].z < (zs + gap*parts[i].r)))
          redo = 1;
        if((bc.wT != PERIODIC) && (parts[i].z > (ze - gap*parts[i].r)))
          redo = 1;
      }

      parts[i].u = 0;
      parts[i].v = 0;
      parts[i].w = 0;
      parts[i].aFx = aFx;
      parts[i].aFy = aFy;
      parts[i].aFz = aFz;
      parts[i].aLx = aLx;
      parts[i].aLy = aLy;
      parts[i].aLz = aLz;
      parts[i].rho = rho;
      parts[i].E = E;
      parts[i].sigma = sigma;
      parts[i].e_dry = e_dry;
      parts[i].order = o;
      parts[i].rs = rs;
      parts[i].spring_k = spring_k;
      parts[i].spring_x = spring_x;
      parts[i].spring_y = spring_y;
      parts[i].spring_z = spring_z;
      parts[i].spring_l = spring_l;
      parts[i].ncoeff = 0;
      parts[i].translating = t;
      parts[i].rotating = r;

      // check that this particle does not intersect any other particle
      fits = !fits;
      for(int j = 0; j < i; j++) {
        xx = parts[i].x - parts[j].x;
        xx = xx * xx;
        yy = parts[i].y - parts[j].y;
        yy = yy * yy;
        zz = parts[i].z - parts[j].z;
        zz = zz * zz;
        if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
          fits = !fits;
          break;
        }

        // also use virtual particle to check if a particle is too close in
        // a periodic direction
        // x only
        if(bc.uW == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // y only
        if(bc.vS == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // z only
        if(bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // x and y
        if(bc.uW == PERIODIC && bc.vS == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // y and z
        if(bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // z and x
        if(bc.uW == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // x, y, and z
        if(bc.uW == PERIODIC && bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[i].x - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          zz = parts[i].z - parts[j].z;
          zz = zz * zz;
          if(parts[i].x < (xs + parts[i].r))
            xx = parts[i].x + xl - parts[j].x;
          if(parts[i].x > (xe - parts[i].r))
            xx = parts[i].x - xl - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (ys + parts[i].r))
            yy = parts[i].y + yl - parts[j].y;
          if(parts[i].y > (ye - parts[i].r))
            yy = parts[i].y - yl - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (zs + parts[i].r))
            zz = parts[i].z + zl - parts[j].z;
          if(parts[i].z > (ze - parts[i].r))
            zz = parts[i].z - zl - parts[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[i].r + gap*parts[j].r)) {
            fits = !fits;
            break;
          }
        }

        // check both ways
        // x only
        if(bc.uW == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // y only
        if(bc.vS == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // z only
        if(bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // x and y
        if(bc.uW == PERIODIC && bc.vS == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // y and z
        if(bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // z and x
        if(bc.uW == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

        // x, y, and z
        if(bc.uW == PERIODIC && bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = parts[j].x - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          zz = parts[j].z - parts[i].z;
          zz = zz * zz;
          if(parts[j].x < (xs + parts[j].r))
            xx = parts[j].x + xl - parts[i].x;
          if(parts[j].x > (xe - parts[j].r))
            xx = parts[j].x - xl - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (ys + parts[j].r))
            yy = parts[j].y + yl - parts[i].y;
          if(parts[j].y > (ye - parts[j].r))
            yy = parts[j].y - yl - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (zs + parts[j].r))
            zz = parts[j].z + zl - parts[i].z;
          if(parts[j].z > (ze - parts[j].r))
            zz = parts[j].z - zl - parts[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (gap*parts[j].r + gap*parts[i].r)) {
            fits = !fits;
            break;
          }
        }

      }
      if(attempts == 1e5*nparts) {
        fail = !fail;
        break;
      }
    }
  }

  if(fail) {
    printf("After %d attempts, the seeder has placed", attempts);
    printf(" %d of %d particles (a = %f).\n\n", i-1, nparts, a);
    printf("...bluebottle seeder done.\n\n");
    exit(EXIT_FAILURE);
  }

  printf("It took %d attempts to place %d", attempts, nparts);
  printf(" particles (a = %f) with no intersections.\n\n", a);
  fflush(stdout);

  printf("Writing part_seeder.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean();
}

void seeder_array(int Nx, int Ny, int Nz, real loa, real a, real aFx, real aFy, 
  real aFz, real aLx, real aLy, real aLz, real rho, real E, real sigma, 
  real e_dry, int o, real rs, real spring_k, real spring_x, 
  real spring_y, real spring_z, real spring_l, int t, int r)
{
  printf("Running bluebottle seeder for %d particles...\n\n", Nx*Ny*Nz);
  fflush(stdout);
  int fail = 0;

  // read domain input
  domain_read_input();
  domain_init();

  nparts = Nx*Ny*Nz;
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  //dx is the distance between centers of two nearby particles in x direction
  real dx = Dom.xl/Nx;
  real dy = Dom.yl/Ny;
  real dz = Dom.zl/Nz;

  if(dx < 2.*a){
    printf("Too many particles in x direction!\n");
    fail = 1;
  }
  if(dy < 2.*a){
    printf("Too many particles in y direction!\n");
    fail = 1; 
  }
  if(dz < 2.*a){
    printf("Too many particles in z direction!\n");
    fail = 1; 
  } 
  if(fail == 1) {
    printf("...bluebottle seeder failed.\n\n");
    exit(EXIT_FAILURE);
  }

  for(int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        parts[i + j*Nx + k*(Nx*Ny)].x = Dom.xs + 0.5*dx + i*dx;
        parts[i + j*Nx + k*(Nx*Ny)].y = Dom.ys + 0.5*dy + j*dy;
        parts[i + j*Nx + k*(Nx*Ny)].z = Dom.zs + 0.5*dz + k*dz;
        parts[i + j*Nx + k*(Nx*Ny)].r = a;
        parts[i + j*Nx + k*(Nx*Ny)].u = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].v = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].w = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].aFx = aFx;
        parts[i + j*Nx + k*(Nx*Ny)].aFy = aFy;
        parts[i + j*Nx + k*(Nx*Ny)].aFz = aFz;
        parts[i + j*Nx + k*(Nx*Ny)].aLx = aLx;
        parts[i + j*Nx + k*(Nx*Ny)].aLy = aLy;
        parts[i + j*Nx + k*(Nx*Ny)].aLz = aLz;
        parts[i + j*Nx + k*(Nx*Ny)].rho = rho;
        parts[i + j*Nx + k*(Nx*Ny)].E = E;
        parts[i + j*Nx + k*(Nx*Ny)].sigma = sigma;
        parts[i + j*Nx + k*(Nx*Ny)].e_dry = e_dry;
        parts[i + j*Nx + k*(Nx*Ny)].order = o;
        parts[i + j*Nx + k*(Nx*Ny)].rs = rs;
        parts[i + j*Nx + k*(Nx*Ny)].spring_k = spring_k;
        parts[i + j*Nx + k*(Nx*Ny)].spring_x = spring_x;
        parts[i + j*Nx + k*(Nx*Ny)].spring_y = spring_y;
        parts[i + j*Nx + k*(Nx*Ny)].spring_z = spring_z;
        parts[i + j*Nx + k*(Nx*Ny)].spring_l = spring_l;
        parts[i + j*Nx + k*(Nx*Ny)].ncoeff = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].translating = t;
        parts[i + j*Nx + k*(Nx*Ny)].rotating = r;       
      }
    }
  }
  printf("Writing part_seeder_array.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder_array.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean(); 
}
void seeder_hex(int Nx, int Ny, int Nz, double ddz, real loa, real a, real aFx, 
  real aFy, real aFz, real aLx, real aLy, real aLz, real rho, real E, 
  real sigma, real e_dry, int o, real rs, real spring_k, 
  real spring_x, real spring_y, real spring_z, real spring_l, int t, int r)
{

  // ddz is the vertical distance from the top of one layer to the middle of the
  // next
  //real ddz = 0.87; //the miminum value of this one is "0.732"


  // read domain input
  domain_read_input();
  domain_init();

  // total number of particles
  // half are large, half are small. if odd number of layers, top layer is large
  int layers = floor(0.5*((real) Nz));
  int large = Nx*Ny;
  int small = (Nx - 1)*(Ny - 1);
  nparts = layers*(large + small) + (Nz % 2)*Nx*Ny;

  printf("The total number of particle is %d \n\n", nparts);
  printf("Running bluebottle seeder for %d particles...\n\n", nparts);
  fflush(stdout);
  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  // dx in the distance between centers of two nearby particles in x direction
  real dx = Dom.xl/Nx;
  real dy = Dom.yl/Ny;

  if(Nz*ddz + 2.*a > Dom.zl){
    printf("Too many layers in z direction");
    exit(EXIT_FAILURE);
  }
  if(Nx*2.*a > Dom.xl){
    printf("Too many layers in x direction");
    exit(EXIT_FAILURE);
  }
  if(Ny*2.*a > Dom.yl){
    printf("Too many layers in y direction");
    exit(EXIT_FAILURE);
  }  
  if(ddz < 0.732){
    printf("Too small distance in neighbouring layer");
    exit(EXIT_FAILURE);
  }
  
  int nx, ny;
  int point;
  int index = 0;
    
  for(int k = 0; k < Nz; k++)
  { 
    point = k % 2;
    if (point == 1) {
      nx = Nx - 1;
      ny = Ny - 1;
    } else {
      nx = Nx;
      ny = Ny;
    }   
    for (int j = 0; j < ny; j++)
    {
      for (int i = 0; i < nx; i++)
      {   
        parts[index].x = Dom.xs + 0.5*dx + i*dx + 0.5*point*dx;
        parts[index].y = Dom.ys + 0.5*dy + j*dy + 0.5*point*dy;
        parts[index].z = Dom.zs + a + k*ddz*a;        
        parts[index].r = a;
        parts[index].u = 0.;
        parts[index].v = 0.;
        parts[index].w = 0.;
        parts[index].aFx = aFx;
        parts[index].aFy = aFy;
        parts[index].aFz = aFz;
        parts[index].aLx = aLx;
        parts[index].aLy = aLy;
        parts[index].aLz = aLz;
        parts[index].rho = rho;
        parts[index].E = E;
        parts[index].sigma = sigma;
        parts[index].e_dry = e_dry;
        parts[index].order = o;
        parts[index].rs = rs;
        parts[index].spring_k = spring_k;
        parts[index].spring_x = spring_x;
        parts[index].spring_y = spring_y;
        parts[index].spring_z = spring_z;
        parts[index].spring_l = spring_l;        
        parts[index].ncoeff = 0.;
        parts[index].translating = t;
        parts[index].rotating = r;
        index = index + 1;       
      }
    }
  }
  printf("Writing part_seeder_hex.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder_hex.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean();    
}
void seeder_high_vol_random(int Nx, int Ny, int Nz, double bias, int nperturb, 
  real loa, real a, real aFx, real aFy, real aFz, real aLx, real aLy, real aLz,
  real rho, real E, real sigma, real e_dry, int o, real rs, 
  real spring_k, real spring_x, real spring_y, real spring_z, real spring_l, 
  int t, int r)
{

  // Generate a random field by perturbing a regular arrary npertrub times and
  // checkinng for interactions

  printf("Running bluebottle seeder for %d particles...\n\n", Nx*Ny*Nz);
  fflush(stdout);
  int fail = 0;
  
  // bias is the ratio of how much it can move 
  // read domain input
  domain_read_input();
  domain_init();
  nparts = Nx*Ny*Nz;

  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  // dx in the distance between centers of two nearby particles in x direction
  real dx = Dom.xl/Nx;
  real dy = Dom.yl/Ny;
  real dz = Dom.zl/Nz;

  if(dx < 2.*a){
    printf(" Too many particles in x direction\n");
    fail = 1;
  }
  if(dy < 2.*a){
    printf("Too many particles in y direction\n");
    fail = 1; 
  }
  if(dz < 2.*a){
    printf("Too many particles in z direction\n");
    fail = 1; 
  } 
  if(fail == 1) {
    printf("...bluebottle seeder done.\n\n");
    exit(EXIT_FAILURE);
  }
  // Set the initial regular domain
  for(int k = 0; k < Nz; k++)
  {
    for (int j = 0; j < Ny; j++)
    {
      for (int i = 0; i < Nx; i++)
      {
        parts[i + j*Nx + k*(Nx*Ny)].x = Dom.xs + (2.*i + 1)*dx*0.5; 
        parts[i + j*Nx + k*(Nx*Ny)].y = Dom.ys + (2.*j + 1)*dy*0.5;
        parts[i + j*Nx + k*(Nx*Ny)].z = Dom.zs + (2.*k + 1)*dz*0.5;  
      }
    }
  }

  real x_new = 0.;
  real y_new = 0.;
  real z_new = 0.;
  real d_min = 100.*a;
  real d_pair= 0.;

  real x_pert = 0.;
  real y_pert = 0.;
  real z_pert = 0.;

  for (int t = 0; t < nperturb; t++) {
    for(int k = 0; k < Nz; k++) {
      for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
          d_min = 100.*a;
          
          x_pert = -1. + 2.*rand() / (real)RAND_MAX;
          y_pert = -1. + 2.*rand() / (real)RAND_MAX;
          z_pert = -1. + 2.*rand() / (real)RAND_MAX;

          x_new = parts[i + j*Nx + k*(Nx*Ny)].x + bias*a*x_pert;
          y_new = parts[i + j*Nx + k*(Nx*Ny)].y + bias*a*y_pert;
          z_new = parts[i + j*Nx + k*(Nx*Ny)].z + bias*a*z_pert;
          if (x_new > Dom.xs && x_new < Dom.xe && 
              y_new > Dom.ys && y_new < Dom.ye && 
              z_new > Dom.zs && z_new < Dom.ze) {
            for (int n = 0; n < Nz; n++) {
              for (int m = 0; m < Ny; m++) {
                for (int l = 0; l < Nx; l++) {
                  // if it calculates the distance to itself
                  if(i == l && j == m && k == n) { 
                    d_pair = 100.*a;
                  } else {
                    d_pair = (x_new - parts[l + m*Nx + n*(Nx*Ny)].x)*
                             (x_new - parts[l + m*Nx + n*(Nx*Ny)].x) +
                             (y_new - parts[l + m*Nx + n*(Nx*Ny)].y)*
                             (y_new - parts[l + m*Nx + n*(Nx*Ny)].y) +
                             (z_new - parts[l + m*Nx + n*(Nx*Ny)].z)*
                             (z_new - parts[l + m*Nx + n*(Nx*Ny)].z);
                    d_pair = sqrt(d_pair);
                  } 
                  if (d_pair < d_min) { 
                    //find the minimum distance between particle pairs
                    d_min = d_pair;
                  } 
                }
              }
            }
            if (d_min > 2.*a) {
            // accept the perturbation if no particle interactions
              parts[i + j*Nx + k*(Nx*Ny)].x = x_new;
              parts[i + j*Nx + k*(Nx*Ny)].y = y_new;
              parts[i + j*Nx + k*(Nx*Ny)].z = z_new;
            } 
          }
        }
      }
    }
  }

 for(int k = 0; k < Nz; k++) {
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        parts[i + j*Nx + k*(Nx*Ny)].r = a;
        parts[i + j*Nx + k*(Nx*Ny)].u = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].v = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].w = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].aFx = aFx;
        parts[i + j*Nx + k*(Nx*Ny)].aFy = aFy;
        parts[i + j*Nx + k*(Nx*Ny)].aFz = aFz;
        parts[i + j*Nx + k*(Nx*Ny)].aLx = aLx;
        parts[i + j*Nx + k*(Nx*Ny)].aLy = aLy;
        parts[i + j*Nx + k*(Nx*Ny)].aLz = aLz;
        parts[i + j*Nx + k*(Nx*Ny)].rho = rho;
        parts[i + j*Nx + k*(Nx*Ny)].E = E;
        parts[i + j*Nx + k*(Nx*Ny)].sigma = sigma;
        parts[i + j*Nx + k*(Nx*Ny)].e_dry = e_dry;
        parts[i + j*Nx + k*(Nx*Ny)].order = o;
        parts[i + j*Nx + k*(Nx*Ny)].rs = rs;
        parts[i + j*Nx + k*(Nx*Ny)].spring_k = spring_k;
        parts[i + j*Nx + k*(Nx*Ny)].spring_x = spring_x;
        parts[i + j*Nx + k*(Nx*Ny)].spring_y = spring_y;
        parts[i + j*Nx + k*(Nx*Ny)].spring_z = spring_z;
        parts[i + j*Nx + k*(Nx*Ny)].spring_l = spring_l;
        parts[i + j*Nx + k*(Nx*Ny)].ncoeff = 0.;
        parts[i + j*Nx + k*(Nx*Ny)].translating = t;
        parts[i + j*Nx + k*(Nx*Ny)].rotating = r;       
      }
    }
  }
  printf("Writing part_seeder_perturbed.config...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder_perturbed.config", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles and compact support length
  fprintf(ofile, "n %d\n", nparts);
  fprintf(ofile, "(l/a) %f\n", loa);

  // write each particle configuration
  for(int i = 0; i < nparts; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", parts[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", parts[i].x, parts[i].y, parts[i].z);
    fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", parts[i].aFx, parts[i].aFy,
      parts[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", parts[i].aLx, parts[i].aLy,
      parts[i].aLz);
    fprintf(ofile, "rho %f\n", parts[i].rho);
    fprintf(ofile, "E %f\n", parts[i].E);
    fprintf(ofile, "sigma %f\n", parts[i].sigma);
    fprintf(ofile, "e_dry %f\n", parts[i].e_dry);
    fprintf(ofile, "order %d\n", parts[i].order);
    fprintf(ofile, "rs/r %f\n", parts[i].rs);
    fprintf(ofile, "spring_k %f\n", parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n",
      parts[i].spring_x, parts[i].spring_y, parts[i].spring_z);
    fprintf(ofile, "spring_l %f\n", parts[i].spring_l);
    fprintf(ofile, "translating %d\n", parts[i].translating);
    fprintf(ofile, "rotating %d\n", parts[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  parts_clean(); 
}  
