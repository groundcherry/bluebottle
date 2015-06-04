/*******************************************************************************
 ******************************* BLUEBOTTLE-1.0 ********************************
 *******************************************************************************
 *
 *  Copyright 2012 - 2014 Adam Sierakowski, The Johns Hopkins University
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

void seeder(int N, real a, real rho, real E, real sigma, int o, int t, int r) {
  printf("Running bluebottle seeder for %d particles...\n\n", N);
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

  nparts = N;

  // allocate particle list
  parts = (part_struct*) malloc(nparts * sizeof(part_struct));
  cpumem += nparts * sizeof(part_struct);

  real e_dry = 1.;
  real l_rough = 1e-4;
  real gap = 1.00;
  real rs = 1.2;

  // place the first particle
  parts[0].r = a;
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].x = rand() / (real)RAND_MAX;
    parts[0].x *= Dom.xl;
    parts[0].x += Dom.xs;
    if((bc.uW != PERIODIC) && (parts[0].x < (Dom.xs + gap*parts[0].r)))
      redo = 1;
    if((bc.uE != PERIODIC) && (parts[0].x > (Dom.xe - gap*parts[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].y = rand() / (real)RAND_MAX;
    parts[0].y *= Dom.yl;
    parts[0].y += Dom.ys;
    if((bc.vS != PERIODIC) && (parts[0].y < (Dom.ys + gap*parts[0].r)))
      redo = 1;
    if((bc.vN != PERIODIC) && (parts[0].y > (Dom.ye - gap*parts[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    parts[0].z = rand() / (real)RAND_MAX;
    //parts[0].z = acos(2.*(parts[0].z-0.5))/PI;
    parts[0].z *= Dom.zl;
    parts[0].z += Dom.zs;
    if((bc.wB != PERIODIC) && (parts[0].z < (Dom.zs + gap*parts[0].r)))
      redo = 1;
    if((bc.wT != PERIODIC) && (parts[0].z > (Dom.ze - gap*parts[0].r)))
      redo = 1;
  }

  parts[0].u = 0;
  parts[0].v = 0;
  parts[0].w = 0;
  parts[0].aFx = 0;
  parts[0].aFy = 0;
  parts[0].aFz = 0;
  parts[0].aLx = 0;
  parts[0].aLy = 0;
  parts[0].aLz = 0;
  parts[0].rho = rho;
  parts[0].E = E;
  parts[0].sigma = sigma;
  parts[0].e_dry = e_dry;
  parts[0].l_rough = l_rough;
  parts[0].order = o;
  parts[0].rs = rs;
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
        parts[i].x *= Dom.xl;
        parts[i].x += Dom.xs;
        if((bc.uW != PERIODIC) && (parts[i].x < (Dom.xs + gap*parts[i].r)))
          redo = 1;
        if((bc.uE != PERIODIC) && (parts[i].x > (Dom.xe - gap*parts[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].y = rand() / (real)RAND_MAX;
        parts[i].y *= Dom.yl;
        parts[i].y += Dom.ys;
        if((bc.vS != PERIODIC) && (parts[i].y < (Dom.ys + gap*parts[i].r)))
          redo = 1;
        if((bc.vN != PERIODIC) && (parts[i].y > (Dom.ye - gap*parts[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        parts[i].z = rand() / (real)RAND_MAX;
        parts[i].z *= Dom.zl;
        parts[i].z += Dom.zs;
        if((bc.wB != PERIODIC) && (parts[i].z < (Dom.zs + gap*parts[i].r)))
          redo = 1;
        if((bc.wT != PERIODIC) && (parts[i].z > (Dom.ze - gap*parts[i].r)))
          redo = 1;
      }

      parts[i].u = 0;
      parts[i].v = 0;
      parts[i].w = 0;
      parts[i].aFx = 0;
      parts[i].aFy = 0;
      parts[i].aFz = 0;
      parts[i].aLx = 0;
      parts[i].aLy = 0;
      parts[i].aLz = 0;
      parts[i].rho = rho;
      parts[i].E = E;
      parts[i].sigma = sigma;
      parts[i].e_dry = e_dry;
      parts[i].l_rough = l_rough;
      parts[i].order = o;
      parts[i].rs = rs;
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
          if(parts[i].x < (Dom.xs + parts[i].r))
            xx = parts[i].x + Dom.xl - parts[j].x;
          if(parts[i].x > (Dom.xe - parts[i].r))
            xx = parts[i].x - Dom.xl - parts[j].x;
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
          if(parts[i].y < (Dom.ys + parts[i].r))
            yy = parts[i].y + Dom.yl - parts[j].y;
          if(parts[i].y > (Dom.ye - parts[i].r))
            yy = parts[i].y - Dom.yl - parts[j].y;
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
          if(parts[i].z < (Dom.zs + parts[i].r))
            zz = parts[i].z + Dom.zl - parts[j].z;
          if(parts[i].z > (Dom.ze - parts[i].r))
            zz = parts[i].z - Dom.zl - parts[j].z;
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
          if(parts[i].x < (Dom.xs + parts[i].r))
            xx = parts[i].x + Dom.xl - parts[j].x;
          if(parts[i].x > (Dom.xe - parts[i].r))
            xx = parts[i].x - Dom.xl - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (Dom.ys + parts[i].r))
            yy = parts[i].y + Dom.yl - parts[j].y;
          if(parts[i].y > (Dom.ye - parts[i].r))
            yy = parts[i].y - Dom.yl - parts[j].y;
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
          if(parts[i].y < (Dom.ys + parts[i].r))
            yy = parts[i].y + Dom.yl - parts[j].y;
          if(parts[i].y > (Dom.ye - parts[i].r))
            yy = parts[i].y - Dom.yl - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (Dom.zs + parts[i].r))
            zz = parts[i].z + Dom.zl - parts[j].z;
          if(parts[i].z > (Dom.ze - parts[i].r))
            zz = parts[i].z - Dom.zl - parts[j].z;
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
          if(parts[i].x < (Dom.xs + parts[i].r))
            xx = parts[i].x + Dom.xl - parts[j].x;
          if(parts[i].x > (Dom.xe - parts[i].r))
            xx = parts[i].x - Dom.xl - parts[j].x;
          xx = xx * xx;
          yy = parts[i].y - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (Dom.zs + parts[i].r))
            zz = parts[i].z + Dom.zl - parts[j].z;
          if(parts[i].z > (Dom.ze - parts[i].r))
            zz = parts[i].z - Dom.zl - parts[j].z;
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
          if(parts[i].x < (Dom.xs + parts[i].r))
            xx = parts[i].x + Dom.xl - parts[j].x;
          if(parts[i].x > (Dom.xe - parts[i].r))
            xx = parts[i].x - Dom.xl - parts[j].x;
          xx = xx * xx;
          if(parts[i].y < (Dom.ys + parts[i].r))
            yy = parts[i].y + Dom.yl - parts[j].y;
          if(parts[i].y > (Dom.ye - parts[i].r))
            yy = parts[i].y - Dom.yl - parts[j].y;
          yy = yy * yy;
          if(parts[i].z < (Dom.zs + parts[i].r))
            zz = parts[i].z + Dom.zl - parts[j].z;
          if(parts[i].z > (Dom.ze - parts[i].r))
            zz = parts[i].z - Dom.zl - parts[j].z;
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
          if(parts[j].x < (Dom.xs + parts[j].r))
            xx = parts[j].x + Dom.xl - parts[i].x;
          if(parts[j].x > (Dom.xe - parts[j].r))
            xx = parts[j].x - Dom.xl - parts[i].x;
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
          if(parts[j].y < (Dom.ys + parts[j].r))
            yy = parts[j].y + Dom.yl - parts[i].y;
          if(parts[j].y > (Dom.ye - parts[j].r))
            yy = parts[j].y - Dom.yl - parts[i].y;
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
          if(parts[j].z < (Dom.zs + parts[j].r))
            zz = parts[j].z + Dom.zl - parts[i].z;
          if(parts[j].z > (Dom.ze - parts[j].r))
            zz = parts[j].z - Dom.zl - parts[i].z;
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
          if(parts[j].x < (Dom.xs + parts[j].r))
            xx = parts[j].x + Dom.xl - parts[i].x;
          if(parts[j].x > (Dom.xe - parts[j].r))
            xx = parts[j].x - Dom.xl - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (Dom.ys + parts[j].r))
            yy = parts[j].y + Dom.yl - parts[i].y;
          if(parts[j].y > (Dom.ye - parts[j].r))
            yy = parts[j].y - Dom.yl - parts[i].y;
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
          if(parts[j].y < (Dom.ys + parts[j].r))
            yy = parts[j].y + Dom.yl - parts[i].y;
          if(parts[j].y > (Dom.ye - parts[j].r))
            yy = parts[j].y - Dom.yl - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (Dom.zs + parts[j].r))
            zz = parts[j].z + Dom.zl - parts[i].z;
          if(parts[j].z > (Dom.ze - parts[j].r))
            zz = parts[j].z - Dom.zl - parts[i].z;
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
          if(parts[j].x < (Dom.xs + parts[j].r))
            xx = parts[j].x + Dom.xl - parts[i].x;
          if(parts[j].x > (Dom.xe - parts[j].r))
            xx = parts[j].x - Dom.xl - parts[i].x;
          xx = xx * xx;
          yy = parts[j].y - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (Dom.zs + parts[j].r))
            zz = parts[j].z + Dom.zl - parts[i].z;
          if(parts[j].z > (Dom.ze - parts[j].r))
            zz = parts[j].z - Dom.zl - parts[i].z;
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
          if(parts[j].x < (Dom.xs + parts[j].r))
            xx = parts[j].x + Dom.xl - parts[i].x;
          if(parts[j].x > (Dom.xe - parts[j].r))
            xx = parts[j].x - Dom.xl - parts[i].x;
          xx = xx * xx;
          if(parts[j].y < (Dom.ys + parts[j].r))
            yy = parts[j].y + Dom.yl - parts[i].y;
          if(parts[j].y > (Dom.ye - parts[j].r))
            yy = parts[j].y - Dom.yl - parts[i].y;
          yy = yy * yy;
          if(parts[j].z < (Dom.zs + parts[j].r))
            zz = parts[j].z + Dom.zl - parts[i].z;
          if(parts[j].z > (Dom.ze - parts[j].r))
            zz = parts[j].z - Dom.zl - parts[i].z;
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

  printf("Writing part_seeder.input...");
  fflush(stdout);
  // write particle configuration to file
  char fname[FILE_NAME_SIZE] = "";
  // open file for writing
  sprintf(fname, "%spart_seeder.input", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of particles
  fprintf(ofile, "n %d\n", nparts);

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
    fprintf(ofile, "spring_k %f\n", 0.);//parts[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n", 0., 0., 0.);
    fprintf(ofile, "spring_l %f\n", 0.);//parts[i].spring_l);
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
