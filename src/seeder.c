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

void seeder_read_input()
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
  real l_rough;        // particle surface roughness
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
  printf("Running bluebottle seeder for %d particles...\n\n", N);
  fflush(stdout);

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
  fret = fscanf(infile, "l_rough %lf\n", &l_rough);
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
  fret = fscanf(infile, "l_rough %f\n", &l_rough);
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
  printf("       l_rough = %lf\n", l_rough);
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
  printf("       l_rough = %f\n", l_rough);
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

  seeder(N, loa, a, aFx, aFy, aFz, aLx, aLy, aLz, rho, E, sigma, e_dry, l_rough,
    order, rs_r, spring_k, spring_x, spring_y, spring_z, spring_l,
    trans, rot);
}

void seeder(int N, real loa, real a, real aFx, real aFy, real aFz, 
  real aLx, real aLy, real aLz, real rho, real E, real sigma, real e_dry,
  real l_rough, int o, real rs, real spring_k, real spring_x, real spring_y,
  real spring_z, real spring_l, int t, int r) {
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
  
  nparts = N;

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
  parts[0].l_rough = l_rough;
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
      parts[i].l_rough = l_rough;
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
    fprintf(ofile, "l_rough %f\n", parts[i].l_rough);
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
