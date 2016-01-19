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

#include "recorder.h"
#include <unistd.h>

void recorder_read_config(void)
{
  int i, j;  // iterator

  char *fret = malloc(sizeof(char) * CHAR_BUF_SIZE);
  fret = fret;  // prevent compiler warning

  // initialize recorder values
  // if rec_X_dt < 0, the output is off
  // otherwise, the fields are all boolean values with 0 = false and 1 = true
  rec_flow_field_dt = -1;
  rec_flow_field_vel = 0;
  rec_flow_field_p = 0;
  rec_flow_field_phase = 0;

  rec_paraview_dt = -1;

  rec_particle_dt = -1;
  rec_particle_pos = 0;
  rec_particle_a = 0;
  rec_particle_vel = 0;
  rec_particle_omega = 0;
  rec_particle_force = 0;
  rec_particle_moment = 0;

  rec_restart_dt = -1;
  rec_prec_dt = -1;
  rec_restart_stop = 0;

  rec_prec_flow_field_dt = -1;
  rec_prec_flow_field_vel = 0;
  rec_prec_flow_field_p = 0;
  rec_prec_flow_field_phase = 0;

  // read the config file
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/record.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  char buf[CHAR_BUF_SIZE] = "";  // character read buffer

  /** list of recognized output configurations **/
  char ***configs;
  int nconfigs = 7;
  int *nconfigsn = (int*) malloc(nconfigs * sizeof(int));
  // cpumem += nconfigs * sizeof(int);
  configs = (char***) malloc(nconfigs * sizeof(char**));
  // cpumum += nconfigs * sizeof(char**);
  int N;
  // flow field info
  N = 0;
  nconfigsn[N] = 3;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigsn[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "FLOW_FIELD");
  sprintf(configs[N][1], "velocity");
  sprintf(configs[N][2], "pressure");
  // ParaView info
  N = 1;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigsn[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "PARAVIEW");
  // particle info
  N = 2;
  nconfigsn[N] = 2;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "PARTICLE");
  sprintf(configs[N][1], "position");
  // restart info
  N = 3;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "RESTART");
  // precursor info
  N = 4;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
  }
  sprintf(configs[N][0], "PRECURSOR_VTK");
  // restart_stop info
  N = 5;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "RESTART_STOP");
  // turbulence flow field info
  N = 6;
  nconfigsn[N] = 3;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigsn[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "PRECURSOR");
  sprintf(configs[N][1], "velocity");
  sprintf(configs[N][2], "pressure");

  // read config file
  char bufconfig[CHAR_BUF_SIZE];
  while(fgets(buf, CHAR_BUF_SIZE, infile) != NULL) {
    // compare first line with predefined configurations
    int foundi = 0;
    for(i = 0; i < nconfigs; i++) {
      // remove newline character
      real dtout = 0;
      sscanf(buf, "%s %lf\n", bufconfig, &dtout);
      // check config list
      if(strcmp(configs[i][0], bufconfig) == 0) {
        // found: continue reading config
        foundi = 1;
        switch(i) {
          case 0: // FLOW_FIELD
            rec_flow_field_dt = dtout;
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            while(strcmp("\n", buf) != 0) {
              int foundj = 0;
              // remove trailing newline
              int ln = strlen(buf) - 1;
              if(buf[ln] == '\n') buf[ln] = '\0';
              // check options
              for(j = 1; j < nconfigsn[i]; j++) {
                if(strcmp(configs[i][j], buf) == 0) {
                  foundj = 1;
                  switch(j) {
                    case 1: // velocity
                      rec_flow_field_vel = 1;
                      break;
                    case 2: // pressure
                      rec_flow_field_p = 1;
                      break;
                    case 3: // phase
                      rec_flow_field_phase = 1;
                      break;
                    default:
                      printf("UNRECOGNIZED OPTION\n");
                  }
                }
              }
              if(!foundj) {
                fprintf(stderr, "Unrecognized record.config ");
                fprintf(stderr, "%s option %s\n", bufconfig, buf);

                // clean up
                fclose(infile);
                for(i = 0; i < nconfigs; i++) {
                  for(j = 0; j < nconfigsn[i]; j++) {
                    free(configs[i][j]);
                  }
                  free(configs[i]);
                }
                free(configs);
                free(nconfigsn);

                exit(EXIT_FAILURE);
              }
              // read next line
              fret = fgets(buf, CHAR_BUF_SIZE, infile);
            }
            break;
          case 1: // PARAVIEW
            rec_paraview_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 2: // PARTICLE
            rec_particle_dt = dtout;
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            while(strcmp("\n", buf) != 0) {
              int foundj = 0;
              // remove trailing newline
              int ln = strlen(buf) - 1;
              if(buf[ln] == '\n') buf[ln] = '\0';
              // check options
              for(j = 1; j < nconfigsn[i]; j++) {
                if(strcmp(configs[i][j], buf) == 0) {
                  foundj = 1;
                  switch(j) {
                    case 1: // position
                      rec_particle_pos = 1;
                      break;
                    case 2: // radius
                      rec_particle_a = 1;
                      break;
                    case 3: // velocity
                      rec_particle_vel = 1;
                      break;
                    case 4: // angular velocity
                      rec_particle_omega = 1;
                      break;
                    case 5: // hydrodynamic force
                      rec_particle_force = 1;
                      break;
                    case 6: // hydrodynamic moment
                      rec_particle_moment = 1;
                      break;
                    default:
                      printf("UNRECOGNIZED OPTION\n");
                  }
                }
              }
              if(!foundj) {
                fprintf(stderr, "Unrecognized record.config ");
                fprintf(stderr, "%s option %s\n", bufconfig, buf);

                // clean up
                fclose(infile);
                for(i = 0; i < nconfigs; i++) {
                  for(j = 0; j < nconfigsn[i]; j++) {
                    free(configs[i][j]);
                  }
                  free(configs[i]);
                }
                free(configs);
                free(nconfigsn);

                exit(EXIT_FAILURE);
              }
              // read next line
              fret = fgets(buf, CHAR_BUF_SIZE, infile);
            }
            break;
          case 3: // RESTART
            rec_restart_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 4: // PRECURSOR_VTK
            rec_prec_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 5: // RESTART_STOP
            rec_restart_stop = (int)dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 6: // PRECURSOR
            rec_prec_flow_field_dt = dtout;
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            /* while(strcmp("\n", buf) != 0) {
              int foundj = 0;
              // remove trailing newline
              int ln = strlen(buf) - 1;
              if(buf[ln] == '\n') buf[ln] = '\0';
              // check options
              for(j = 1; j < nconfigsn[i]; j++) {
                if(strcmp(configs[i][j], buf) == 0) {
                  foundj = 1;
                  switch(j) {
                    case 1: // velocity
                      rec_prec_flow_field_vel = 1;
                      break;
                    case 2: // pressure
                      rec_prec_flow_field_p = 1;
                      break;
                    case 3: // phase
                      rec_prec_flow_field_phase = 1;
                      break;
                    default:
                      printf("UNRECOGNIZED OPTION\n");
                  }
                }
              }
              if(!foundj) {
                fprintf(stderr, "Unrecognized record.config ");
                fprintf(stderr, "%s option %s\n", bufconfig, buf);

                // clean up
                fclose(infile);
                for(i = 0; i < nconfigs; i++) {
                  for(j = 0; j < nconfigsn[i]; j++) {
                    free(configs[i][j]);
                  }
                  free(configs[i]);
                }
                free(configs);
                free(nconfigsn);

                exit(EXIT_FAILURE);
              }*/
            break;
          default:
            printf("UNRECOGNIZED TYPE\n");
        }
      }
    }
    // if we've parsed the entire list but don't find anything: error
    if(!foundi) {
      fprintf(stderr, "Unrecognized record.config type %s\n", bufconfig);

      // clean up
      fclose(infile);
      for(i = 0; i < nconfigs; i++) {
        for(j = 0; j < nconfigsn[i]; j++) {
          free(configs[i][j]);
        }
        free(configs[i]);
      }
      free(configs);
      free(nconfigsn);

      exit(EXIT_FAILURE);
    }
  }

  // clean up
  fclose(infile);
  for(i = 0; i < nconfigs; i++) {
    for(j = 0; j < nconfigsn[i]; j++) {
      free(configs[i][j]);
    }
    free(configs[i]);
  }
  free(configs);
  free(nconfigsn);
}

void cgns_grid(void)
{
  // create the file
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/output/%s", ROOT_DIR, "grid.cgns");
  int fn;
  int bn;
  int zn;
  int gn;
  int cn;
  if (access(fname, F_OK) != -1) {
    // file exists
  } else {
    // file does not exist
    cg_open(fname, CG_MODE_WRITE, &fn);
    cg_base_write(fn, "Base", 3, 3, &bn);
    cgsize_t size[9];
    size[0] = Dom.xn+1; // cells -> vertices
    size[1] = Dom.yn+1;
    size[2] = Dom.zn+1;
    size[3] = Dom.xn;
    size[4] = Dom.yn;
    size[5] = Dom.zn;
    size[6] = 0;
    size[7] = 0;
    size[8] = 0;
    cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
    cg_grid_write(fn, bn, zn, "GridCoordinates", &gn);

    real *x = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
    // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
    real *y = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
    // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
    real *z = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
    // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
    for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke+1; k++) {
      for(int j = Dom.Gcc.js; j < Dom.Gcc.je+1; j++) {
        for(int i = Dom.Gcc.is; i < Dom.Gcc.ie+1; i++) {
          int C = (i-1) + (j-1)*(Dom.xn+1) + (k-1)*(Dom.xn+1)*(Dom.yn+1);
          x[C] = Dom.xs + (i-1)*Dom.dx;
          y[C] = Dom.ys + (j-1)*Dom.dy;
          z[C] = Dom.zs + (k-1)*Dom.dz;
        }
      }
    }
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &cn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &cn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &cn);

    free(x);
    free(y);
    free(z);

    cg_close(fn);
  }
}

void cgns_turb_grid(void)
{
  // create the file
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/output/%s", ROOT_DIR, "prec-grid.cgns");
  int fn;
  int bn;
  int zn;
  int gn;
  int cn;
  if(access(fname, F_OK) != -1) {
    // file exists
  } else {
    cg_open(fname, CG_MODE_WRITE, &fn);
    cg_base_write(fn, "Base", 3, 3, &bn);
    cgsize_t size[9];
    size[0] = Dom.xn+1; // cells -> vertices
    size[1] = Dom.yn+1;
    size[2] = Dom.zn+1;
    size[3] = Dom.xn;
    size[4] = Dom.yn;
    size[5] = Dom.zn;
    size[6] = 0;
    size[7] = 0;
    size[8] = 0;
    cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
    cg_grid_write(fn, bn, zn, "GridCoordinates", &gn);

    real *x = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
    // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
    real *y = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
    // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
    real *z = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
    // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
    for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke+1; k++) {
        for(int j = Dom.Gcc.js; j < Dom.Gcc.je+1; j++) {
            for(int i = Dom.Gcc.is; i < Dom.Gcc.ie+1; i++) {
                int C = (i-1) + (j-1)*(Dom.xn+1) + (k-1)*(Dom.xn+1)*(Dom.yn+1);
                x[C] = Dom.xs + (i-1)*Dom.dx;
                y[C] = Dom.ys + (j-1)*Dom.dy;
                z[C] = Dom.zs + (k-1)*Dom.dz;
            }
        }
    }

    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &cn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &cn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &cn);

    free(x);
    free(y);
    free(z);

    cg_close(fn);
  }
}

void cgns_flow_field(real dtout)
{
  // create the solution file
  char fname[FILE_NAME_SIZE] = "";
  char fname2[FILE_NAME_SIZE] = "";
  char fnameall[FILE_NAME_SIZE] = "";
  char fnameall2[FILE_NAME_SIZE] = "";
  char gname[FILE_NAME_SIZE] = "";
  char gnameall[FILE_NAME_SIZE] = "";
  real tout = ttime; //  = rec_flow_field_stepnum_out * dtout;
  char format[CHAR_BUF_SIZE] = "";
  char snodename[CHAR_BUF_SIZE] = "";
  char snodenameall[CHAR_BUF_SIZE] = "";
  int sigfigs = ceil(log10(1. / dtout));
  if(sigfigs < 1) sigfigs = 1;
  sprintf(format, "%%.%df", sigfigs);
  sprintf(fname2, "flow-%s.cgns", format);
  sprintf(fnameall2, "%s/output/flow-%s.cgns", ROOT_DIR, format);
  sprintf(snodename, "Solution-");
  sprintf(snodenameall, "/Base/Zone0/Solution-");
  sprintf(snodename, "%s%s", snodename, format);
  sprintf(snodenameall, "%s%s", snodenameall, format);
  sprintf(fname, fname2, tout);
  sprintf(fnameall, fnameall2, tout);
  sprintf(snodename, snodename, tout);
  sprintf(snodenameall, snodenameall, tout);
  sprintf(gname, "grid.cgns");
  sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
  int fn;
  int bn;
  int zn;
  int sn;
  int fnpress;
  int fnphase;
  int fnu;
  int fnv;
  int fnw;
  cg_open(fnameall, CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 3, 3, &bn);
  cgsize_t size[9];
  size[0] = Dom.xn+1; // cells -> vertices
  size[1] = Dom.yn+1;
  size[2] = Dom.zn+1;
  size[3] = Dom.xn;
  size[4] = Dom.yn;
  size[5] = Dom.zn;
  size[6] = 0;
  size[7] = 0;
  size[8] = 0;
  cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
  cg_goto(fn, bn, "Zone_t", zn, "end");
  // check that grid.cgns exists
  /*int fng;
  if(cg_open(gnameall, CG_MODE_READ, &fng) != 0) {
    fprintf(stderr, "CGNS flow field write failure: no grid.cgns\n");
    exit(EXIT_FAILURE);
  } else {
    cg_close(fng);
  }
    cg_close(fng);
*/
  
  cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");

  cg_sol_write(fn, bn, zn, "Solution", CellCenter, &sn);
  real *pout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        pout[C] = p[CC];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "Pressure", pout, &fnpress);

  real *uout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(int i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = (i-1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC1 = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC2 = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC3 = (i+2) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        uout[C] = -0.0625*u[CC0] + 0.5625*u[CC1] + 0.5625*u[CC2] - 0.0625*u[CC3];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", uout, &fnu);

  real *vout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(int j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
      for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC1 = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC2 = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC3 = i + (j+2)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        vout[C] = -0.0625*v[CC0] + 0.5625*v[CC1] + 0.5625*v[CC2] - 0.0625*v[CC3];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", vout, &fnv);

  real *wout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
    for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
        int CC1 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        int CC2 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
        int CC3 = i + j*Dom.Gfz.s1b + (k+2)*Dom.Gfz.s2b;
        wout[C] = -0.0625*w[CC0] + 0.5625*w[CC1] + 0.5625*w[CC2] - 0.0625*w[CC3];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", wout, &fnw);

  int *phaseout = malloc(Dom.Gcc.s3 * sizeof(int));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        phaseout[C] = phase[CC];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, Integer, "Phase", phaseout, &fnphase);

  cg_user_data_write("Etc");
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cgsize_t *N = malloc(sizeof(cgsize_t));
  N[0] = 1;
  cg_array_write("Time", RealDouble, 1, N, &ttime);
  cg_array_write("Density", RealDouble, 1, N, &rho_f);
  cg_array_write("KinematicViscosity", RealDouble, 1, N, &nu);
  free(N);

  cg_close(fn);
  free(pout);
  free(uout);
  free(vout);
  free(wout);
  free(phaseout);

//  // now link this timestep into time series
//  // create the time series file if it doesn't exist
//  char sname[FILE_NAME_SIZE];
//  sprintf(sname, "%s/output/flow.cgns", ROOT_DIR);
//  int fns;
//  int bns;
//  int depth;
//  char *label = malloc(CHAR_BUF_SIZE * sizeof(real));
//  cpumem += CHAR_BUF_SIZE * sizeof(real);
//  char *bitername = malloc(CHAR_BUF_SIZE * sizeof(real));
//  cpumem += CHAR_BUF_SIZE * sizeof(real);
//  int index;
//  int zns;
//  real *tseries = malloc(sizeof(real));
//  cpumem += sizeof(real);
//  tseries[0] = 0;
//  cgsize_t tsize[1];
//  tsize[0] = 1;
//
//  if(cg_is_cgns(sname, &fns) != 0) {
//    // if it does not exist, create it
//    cg_open(sname, CG_MODE_WRITE, &fns);
//    cg_base_write(fn, "Base", 3, 3, &bns);
//    //cgsize_t size[9];
//    size[0] = Dom.xn+1; // cells -> vertices
//    size[1] = Dom.yn+1;
//    size[2] = Dom.zn+1;
//    size[3] = Dom.xn;
//    size[4] = Dom.yn;
//    size[5] = Dom.zn;
//    size[6] = 0;
//    size[7] = 0;
//    size[8] = 0;
//    cg_zone_write(fns, bns, "Zone0", size, Structured, &zns);
//    cg_goto(fns, bns, "Zone_t", zns, "end");
//    cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");
//    char solname2[CHAR_BUF_SIZE];
//    sprintf(solname2, "Solution-%s", format);
//    sprintf(solname2, solname2, tout);
//    cg_link_write(solname2, fname, "Base/Zone0/Solution");
//    cg_biter_write(fns, bns, "BaseIterativeData", 1);
//    cg_goto(fns, bns, "BaseIterativeData_t", 1, "end");
//    int N = 1;
//    cg_array_write("TimeValues", RealDouble, 1, &N, &ttime);
//    cg_ziter_write(fns, bns, zns, "ZoneIterativeData");
//    cg_goto(fns, bns, "Zone_t", zns, "ZoneIterativeData_t", 1, "end");
//    cg_simulation_type_write(fns, bns, TimeAccurate);
//    cg_close(fns);
//  } else {
//    // modify the timeseries file
//    cg_open(sname, CG_MODE_MODIFY, &fns);
//    cg_gopath(fns, "/Base");
//    cg_where(&fns, &bns, &depth, &label, &index);
//    int N = 0;
//    cg_biter_read(fns, bns, bitername, &N);
//    free(tseries);
//    tseries = malloc((N+1) * sizeof(real));
//    cpumem = (N+1) * sizeof(real);
//    cg_goto(fns, bns, "BaseIterativeData_t", 1, "end");
//    cg_array_read(1, tseries);
//    cg_biter_write(fns, bns, "BaseIterativeData", N+1);
//    cg_goto(fns, bns, "BaseIterativeData_t", 1, "end");
//    tseries[N] = ttime;
//    tsize[0] = N+1;
//    cg_array_write("TimeValues", RealDouble, 1, tsize, tseries);
//    cg_goto(fns, bns, "Zone_t", 1, "end");
//    char solname2[CHAR_BUF_SIZE];
//    sprintf(solname2, "Solution-%s", format);
//    sprintf(solname2, solname2, tout);
//    cg_link_write(solname2, fname, "Base/Zone0/Solution");
//    cg_close(fns);
//   
//    free(tseries);
//    free(label);
//    free(bitername);
//  }
}

void cgns_turb_flow_field(real dtout)
{
  // create the solution file
  char fname[FILE_NAME_SIZE] = "";
  char fname2[FILE_NAME_SIZE] = "";
  char fnameall[FILE_NAME_SIZE] = "";
  char fnameall2[FILE_NAME_SIZE] = "";
  char gname[FILE_NAME_SIZE] = "";
  char gnameall[FILE_NAME_SIZE] = "";
  real tout = ttime; //  = rec_flow_field_stepnum_out * dtout;
  char format[CHAR_BUF_SIZE] = "";
  char snodename[CHAR_BUF_SIZE] = "";
  char snodenameall[CHAR_BUF_SIZE] = "";
  int sigfigs = ceil(log10(1. / dtout));
  if(sigfigs < 1) sigfigs = 1;
  sprintf(format, "%%.%df", sigfigs);
  sprintf(fname2, "prec-flow-%s.cgns", format);
  sprintf(fnameall2, "%s/output/prec-flow-%s.cgns", ROOT_DIR, format);
  sprintf(snodename, "Solution-");
  sprintf(snodenameall, "/Base/Zone0/Solution-");
  sprintf(snodename, "%s%s", snodename, format);
  sprintf(snodenameall, "%s%s", snodenameall, format);
  sprintf(fname, fname2, tout);
  sprintf(fnameall, fnameall2, tout);
  sprintf(snodename, snodename, tout);
  sprintf(snodenameall, snodenameall, tout);
  sprintf(gname, "prec-grid.cgns");
  sprintf(gnameall, "%s/output/%s", ROOT_DIR, "prec-grid.cgns");
  int fn;
  int bn;
  int zn;
  int sn;
  int fnpress;
  int fnphase;
  int fnu;
  int fnv;
  int fnw;
  cg_open(fnameall, CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 3, 3, &bn);
  cgsize_t size[9];
  size[0] = Dom.xn+1; // cells -> vertices
  size[1] = Dom.yn+1;
  size[2] = Dom.zn+1;
  size[3] = Dom.xn;
  size[4] = Dom.yn;
  size[5] = Dom.zn;
  size[6] = 0;
  size[7] = 0;
  size[8] = 0;
  cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
  cg_goto(fn, bn, "Zone_t", zn, "end");
  // check that grid.cgns exists
  /*int fng;
  if(cg_open(gnameall, CG_MODE_READ, &fng) != 0) {
    fprintf(stderr, "CGNS flow field write failure: no grid.cgns\n");
    exit(EXIT_FAILURE);
  } else {
    cg_close(fng);
  }
    cg_close(fng);
*/

  cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");

  cg_sol_write(fn, bn, zn, "Solution", CellCenter, &sn);
  real *pout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        pout[C] = p[CC];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "Pressure", pout, &fnpress);

  real *uout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(int i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = (i-1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC1 = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC2 = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC3 = (i+2) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        uout[C] = -0.0625*u[CC0] + 0.5625*u[CC1] + 0.5625*u[CC2] - 0.0625*u[CC3];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", uout, &fnu);

  real *vout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(int j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
      for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC1 = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC2 = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC3 = i + (j+2)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        vout[C] = -0.0625*v[CC0] + 0.5625*v[CC1] + 0.5625*v[CC2] - 0.0625*v[CC3];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", vout, &fnv);

  real *wout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
    for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
        int CC1 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        int CC2 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
        int CC3 = i + j*Dom.Gfz.s1b + (k+2)*Dom.Gfz.s2b;
        wout[C] = -0.0625*w[CC0] + 0.5625*w[CC1] + 0.5625*w[CC2] - 0.0625*w[CC3];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", wout, &fnw);

  int *phaseout = malloc(Dom.Gcc.s3 * sizeof(int));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        phaseout[C] = phase[CC];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, Integer, "Phase", phaseout, &fnphase);

  cg_user_data_write("Etc");
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cgsize_t *N = malloc(sizeof(cgsize_t));
  N[0] = 1;
  cg_array_write("Time", RealDouble, 1, N, &ttime);
  free(N);

  cg_close(fn);
  free(pout);
  free(uout);
  free(vout);
  free(wout);
  free(phaseout);
}

void cgns_particles(real dtout)
{
  if(nparts > 0) {
    // create the solution file
    char fname[FILE_NAME_SIZE] = "";
    char fname2[FILE_NAME_SIZE] = "";
    char fnameall[FILE_NAME_SIZE] = "";
    char fnameall2[FILE_NAME_SIZE] = "";
    real tout = ttime; // = rec_particle_stepnum_out * dtout;
    char format[CHAR_BUF_SIZE] = "";
    int sigfigs = ceil(log10(1. / dtout));
    if(sigfigs < 1) sigfigs = 1;
    sprintf(format, "%%.%df", sigfigs);
    sprintf(fname2, "part-%s.cgns", format);
    sprintf(fnameall2, "%s/output/part-%s.cgns", ROOT_DIR, format);
    sprintf(fname, fname2, tout);
    sprintf(fnameall, fnameall2, tout);
    int fn;
    int bn;
    int zn;
    int en;
    int sn;
    int Xn;
    int Yn;
    int Zn;
    int fnr;
    cg_open(fnameall, CG_MODE_WRITE, &fn);
    cg_base_write(fn, "Base", 3, 3, &bn);
    cgsize_t size[3][1];
    size[0][0] = nparts;
    size[1][0] = 0;
    size[2][0] = 0;
    cg_zone_write(fn, bn, "Zone0", size[0], Unstructured, &zn);

    // write particle locations
    real *x = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *y = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *z = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    cgsize_t *conn = malloc(nparts * sizeof(cgsize_t));
    // cpumem += nparts * sizeof(int);
    real *a = malloc(nparts * sizeof(real));
    real *rho = malloc(nparts * sizeof(real));
    real *E = malloc(nparts * sizeof(real));
    real *sigma = malloc(nparts * sizeof(real));
    real *e_dry = malloc(nparts * sizeof(real));
    real *coeff_fric = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(int);
    int *order = malloc(nparts * sizeof(int));
    // cpumem += nparts * sizeof(real);
    real *u = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *v = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *w = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *udot = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *vdot = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *wdot = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *ox = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *oy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *oz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *Fx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *Fy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *Fz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *Lx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *Ly = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *Lz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real); 
    real *iFx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *iFy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *iFz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *iLx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *iLy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *iLz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real); 
    real *kFx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *kFy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *kFz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *hFx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *hFy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *hFz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *hLx = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *hLy = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real);
    real *hLz = malloc(nparts * sizeof(real));
    // cpumem += nparts * sizeof(real); 
    real *axx = malloc(nparts * sizeof(real));
    real *axy = malloc(nparts * sizeof(real));
    real *axz = malloc(nparts * sizeof(real));
    real *ayx = malloc(nparts * sizeof(real));
    real *ayy = malloc(nparts * sizeof(real));
    real *ayz = malloc(nparts * sizeof(real));
    real *azx = malloc(nparts * sizeof(real));
    real *azy = malloc(nparts * sizeof(real));
    real *azz = malloc(nparts * sizeof(real));

    real   *p00_r = malloc(nparts * sizeof(real));
    real   *p00_i = malloc(nparts * sizeof(real));
    real *phi00_r = malloc(nparts * sizeof(real));
    real *phi00_i = malloc(nparts * sizeof(real));
    real *chi00_r = malloc(nparts * sizeof(real));
    real *chi00_i = malloc(nparts * sizeof(real));
    real   *p10_r = malloc(nparts * sizeof(real));
    real   *p10_i = malloc(nparts * sizeof(real));
    real *phi10_r = malloc(nparts * sizeof(real));
    real *phi10_i = malloc(nparts * sizeof(real));
    real *chi10_r = malloc(nparts * sizeof(real));
    real *chi10_i = malloc(nparts * sizeof(real));
    real   *p11_r = malloc(nparts * sizeof(real));
    real   *p11_i = malloc(nparts * sizeof(real));
    real *phi11_r = malloc(nparts * sizeof(real));
    real *phi11_i = malloc(nparts * sizeof(real));
    real *chi11_r = malloc(nparts * sizeof(real));
    real *chi11_i = malloc(nparts * sizeof(real));
    real   *p20_r = malloc(nparts * sizeof(real));
    real   *p20_i = malloc(nparts * sizeof(real));
    real *phi20_r = malloc(nparts * sizeof(real));
    real *phi20_i = malloc(nparts * sizeof(real));
    real *chi20_r = malloc(nparts * sizeof(real));
    real *chi20_i = malloc(nparts * sizeof(real));
    real   *p21_r = malloc(nparts * sizeof(real));
    real   *p21_i = malloc(nparts * sizeof(real));
    real *phi21_r = malloc(nparts * sizeof(real));
    real *phi21_i = malloc(nparts * sizeof(real));
    real *chi21_r = malloc(nparts * sizeof(real));
    real *chi21_i = malloc(nparts * sizeof(real));
    real   *p22_r = malloc(nparts * sizeof(real));
    real   *p22_i = malloc(nparts * sizeof(real));
    real *phi22_r = malloc(nparts * sizeof(real));
    real *phi22_i = malloc(nparts * sizeof(real));
    real *chi22_r = malloc(nparts * sizeof(real));
    real *chi22_i = malloc(nparts * sizeof(real));
    real   *p30_r = malloc(nparts * sizeof(real));
    real   *p30_i = malloc(nparts * sizeof(real));
    real *phi30_r = malloc(nparts * sizeof(real));
    real *phi30_i = malloc(nparts * sizeof(real));
    real *chi30_r = malloc(nparts * sizeof(real));
    real *chi30_i = malloc(nparts * sizeof(real));
    real   *p31_r = malloc(nparts * sizeof(real));
    real   *p31_i = malloc(nparts * sizeof(real));
    real *phi31_r = malloc(nparts * sizeof(real));
    real *phi31_i = malloc(nparts * sizeof(real));
    real *chi31_r = malloc(nparts * sizeof(real));
    real *chi31_i = malloc(nparts * sizeof(real));
    real   *p32_r = malloc(nparts * sizeof(real));
    real   *p32_i = malloc(nparts * sizeof(real));
    real *phi32_r = malloc(nparts * sizeof(real));
    real *phi32_i = malloc(nparts * sizeof(real));
    real *chi32_r = malloc(nparts * sizeof(real));
    real *chi32_i = malloc(nparts * sizeof(real));
    real   *p33_r = malloc(nparts * sizeof(real));
    real   *p33_i = malloc(nparts * sizeof(real));
    real *phi33_r = malloc(nparts * sizeof(real));
    real *phi33_i = malloc(nparts * sizeof(real));
    real *chi33_r = malloc(nparts * sizeof(real));
    real *chi33_i = malloc(nparts * sizeof(real));
    real   *p40_r = malloc(nparts * sizeof(real));
    real   *p40_i = malloc(nparts * sizeof(real));
    real *phi40_r = malloc(nparts * sizeof(real));
    real *phi40_i = malloc(nparts * sizeof(real));
    real *chi40_r = malloc(nparts * sizeof(real));
    real *chi40_i = malloc(nparts * sizeof(real));
    real   *p41_r = malloc(nparts * sizeof(real));
    real   *p41_i = malloc(nparts * sizeof(real));
    real *phi41_r = malloc(nparts * sizeof(real));
    real *phi41_i = malloc(nparts * sizeof(real));
    real *chi41_r = malloc(nparts * sizeof(real));
    real *chi41_i = malloc(nparts * sizeof(real));
    real   *p42_r = malloc(nparts * sizeof(real));
    real   *p42_i = malloc(nparts * sizeof(real));
    real *phi42_r = malloc(nparts * sizeof(real));
    real *phi42_i = malloc(nparts * sizeof(real));
    real *chi42_r = malloc(nparts * sizeof(real));
    real *chi42_i = malloc(nparts * sizeof(real));
    real   *p43_r = malloc(nparts * sizeof(real));
    real   *p43_i = malloc(nparts * sizeof(real));
    real *phi43_r = malloc(nparts * sizeof(real));
    real *phi43_i = malloc(nparts * sizeof(real));
    real *chi43_r = malloc(nparts * sizeof(real));
    real *chi43_i = malloc(nparts * sizeof(real));
    real   *p44_r = malloc(nparts * sizeof(real));
    real   *p44_i = malloc(nparts * sizeof(real));
    real *phi44_r = malloc(nparts * sizeof(real));
    real *phi44_i = malloc(nparts * sizeof(real));
    real *chi44_r = malloc(nparts * sizeof(real));
    real *chi44_i = malloc(nparts * sizeof(real));
    real   *p50_r = malloc(nparts * sizeof(real));
    real   *p50_i = malloc(nparts * sizeof(real));
    real *phi50_r = malloc(nparts * sizeof(real));
    real *phi50_i = malloc(nparts * sizeof(real));
    real *chi50_r = malloc(nparts * sizeof(real));
    real *chi50_i = malloc(nparts * sizeof(real));
    real   *p51_r = malloc(nparts * sizeof(real));
    real   *p51_i = malloc(nparts * sizeof(real));
    real *phi51_r = malloc(nparts * sizeof(real));
    real *phi51_i = malloc(nparts * sizeof(real));
    real *chi51_r = malloc(nparts * sizeof(real));
    real *chi51_i = malloc(nparts * sizeof(real));
    real   *p52_r = malloc(nparts * sizeof(real));
    real   *p52_i = malloc(nparts * sizeof(real));
    real *phi52_r = malloc(nparts * sizeof(real));
    real *phi52_i = malloc(nparts * sizeof(real));
    real *chi52_r = malloc(nparts * sizeof(real));
    real *chi52_i = malloc(nparts * sizeof(real));
    real   *p53_r = malloc(nparts * sizeof(real));
    real   *p53_i = malloc(nparts * sizeof(real));
    real *phi53_r = malloc(nparts * sizeof(real));
    real *phi53_i = malloc(nparts * sizeof(real));
    real *chi53_r = malloc(nparts * sizeof(real));
    real *chi53_i = malloc(nparts * sizeof(real));
    real   *p54_r = malloc(nparts * sizeof(real));
    real   *p54_i = malloc(nparts * sizeof(real));
    real *phi54_r = malloc(nparts * sizeof(real));
    real *phi54_i = malloc(nparts * sizeof(real));
    real *chi54_r = malloc(nparts * sizeof(real));
    real *chi54_i = malloc(nparts * sizeof(real));
    real   *p55_r = malloc(nparts * sizeof(real));
    real   *p55_i = malloc(nparts * sizeof(real));
    real *phi55_r = malloc(nparts * sizeof(real));
    real *phi55_i = malloc(nparts * sizeof(real));
    real *chi55_r = malloc(nparts * sizeof(real));
    real *chi55_i = malloc(nparts * sizeof(real));

    for(int i = 0; i < nparts; i++) {
      real mass = 4./3.*PI*(parts[i].rho-rho_f)
        *parts[i].r*parts[i].r*parts[i].r;
      x[i] = parts[i].x;
      y[i] = parts[i].y;
      z[i] = parts[i].z;
      conn[i] = nparts-i;
      a[i] = parts[i].r;
      rho[i] = parts[i].rho;
      E[i] = parts[i].E;
      sigma[i] = parts[i].sigma;
      e_dry[i] = parts[i].e_dry;
      coeff_fric[i] = parts[i].coeff_fric;
      order[i] = parts[i].order;
      u[i] = parts[i].u;
      v[i] = parts[i].v;
      w[i] = parts[i].w;
      udot[i] = parts[i].udot;
      vdot[i] = parts[i].vdot;
      wdot[i] = parts[i].wdot;
      kFx[i] = parts[i].kFx;
      kFy[i] = parts[i].kFy;
      kFz[i] = parts[i].kFz;
      iFx[i] = parts[i].iFx;
      iFy[i] = parts[i].iFy;
      iFz[i] = parts[i].iFz;
      iLx[i] = parts[i].iLx;
      iLy[i] = parts[i].iLy;
      iLz[i] = parts[i].iLz;
      hFx[i] = parts[i].Fx;
      hFy[i] = parts[i].Fy;
      hFz[i] = parts[i].Fz;
      hLx[i] = parts[i].Lx;
      hLy[i] = parts[i].Ly;
      hLz[i] = parts[i].Lz;
      Fx[i] = kFx[i] + iFx[i] + hFx[i] + parts[i].aFx + mass*g.x;
      Fy[i] = kFy[i] + iFy[i] + hFy[i] + parts[i].aFy + mass*g.y;
      Fz[i] = kFz[i] + iFz[i] + hFz[i] + parts[i].aFz + mass*g.z;
      Lx[i] = iLx[i] + hLx[i];
      Ly[i] = iLy[i] + hLy[i];
      Lz[i] = iLz[i] + hLz[i];
      ox[i] = parts[i].ox;
      oy[i] = parts[i].oy;
      oz[i] = parts[i].oz;
      axx[i] = parts[i].axx;
      axy[i] = parts[i].axy;
      axz[i] = parts[i].axz;
      ayx[i] = parts[i].ayx;
      ayy[i] = parts[i].ayy;
      ayz[i] = parts[i].ayz;
      azx[i] = parts[i].azx;
      azy[i] = parts[i].azy;
      azz[i] = parts[i].azz;

        p00_r[i] = 0;
        p00_i[i] = 0;
      phi00_r[i] = 0;
      phi00_i[i] = 0;
      chi00_r[i] = 0;
      chi00_i[i] = 0;
        p10_r[i] = 0;
        p10_i[i] = 0;
      phi10_r[i] = 0;
      phi10_i[i] = 0;
      chi10_r[i] = 0;
      chi10_i[i] = 0;
        p11_r[i] = 0;
        p11_i[i] = 0;
      phi11_r[i] = 0;
      phi11_i[i] = 0;
      chi11_r[i] = 0;
      chi11_i[i] = 0;
        p20_r[i] = 0;
        p20_i[i] = 0;
      phi20_r[i] = 0;
      phi20_i[i] = 0;
      chi20_r[i] = 0;
      chi20_i[i] = 0;
        p21_r[i] = 0;
        p21_i[i] = 0;
      phi21_r[i] = 0;
      phi21_i[i] = 0;
      chi21_r[i] = 0;
      chi21_i[i] = 0;
        p22_r[i] = 0;
        p22_i[i] = 0;
      phi22_r[i] = 0;
      phi22_i[i] = 0;
      chi22_r[i] = 0;
      chi22_i[i] = 0;
        p30_r[i] = 0;
        p30_i[i] = 0;
      phi30_r[i] = 0;
      phi30_i[i] = 0;
      chi30_r[i] = 0;
      chi30_i[i] = 0;
        p31_r[i] = 0;
        p31_i[i] = 0;
      phi31_r[i] = 0;
      phi31_i[i] = 0;
      chi31_r[i] = 0;
      chi31_i[i] = 0;
        p32_r[i] = 0;
        p32_i[i] = 0;
      phi32_r[i] = 0;
      phi32_i[i] = 0;
      chi32_r[i] = 0;
      chi32_i[i] = 0;
        p33_r[i] = 0;
        p33_i[i] = 0;
      phi33_r[i] = 0;
      phi33_i[i] = 0;
      chi33_r[i] = 0;
      chi33_i[i] = 0;
        p40_r[i] = 0;
        p40_i[i] = 0;
      phi40_r[i] = 0;
      phi40_i[i] = 0;
      chi40_r[i] = 0;
      chi40_i[i] = 0;
        p41_r[i] = 0;
        p41_i[i] = 0;
      phi41_r[i] = 0;
      phi41_i[i] = 0;
      chi41_r[i] = 0;
      chi41_i[i] = 0;
        p42_r[i] = 0;
        p42_i[i] = 0;
      phi42_r[i] = 0;
      phi42_i[i] = 0;
      chi42_r[i] = 0;
      chi42_i[i] = 0;
        p43_r[i] = 0;
        p43_i[i] = 0;
      phi43_r[i] = 0;
      phi43_i[i] = 0;
      chi43_r[i] = 0;
      chi43_i[i] = 0;
        p44_r[i] = 0;
        p44_i[i] = 0;
      phi44_r[i] = 0;
      phi44_i[i] = 0;
      chi44_r[i] = 0;
      chi44_i[i] = 0;
        p50_r[i] = 0;
        p50_i[i] = 0;
      phi50_r[i] = 0;
      phi50_i[i] = 0;
      chi50_r[i] = 0;
      chi50_i[i] = 0;
        p51_r[i] = 0;
        p51_i[i] = 0;
      phi51_r[i] = 0;
      phi51_i[i] = 0;
      chi51_r[i] = 0;
      chi51_i[i] = 0;
        p52_r[i] = 0;
        p52_i[i] = 0;
      phi52_r[i] = 0;
      phi52_i[i] = 0;
      chi52_r[i] = 0;
      chi52_i[i] = 0;
        p53_r[i] = 0;
        p53_i[i] = 0;
      phi53_r[i] = 0;
      phi53_i[i] = 0;
      chi53_r[i] = 0;
      chi53_i[i] = 0;
        p54_r[i] = 0;
        p54_i[i] = 0;
      phi54_r[i] = 0;
      phi54_i[i] = 0;
      chi54_r[i] = 0;
      chi54_i[i] = 0;
        p55_r[i] = 0;
        p55_i[i] = 0;
      phi55_r[i] = 0;
      phi55_i[i] = 0;
      chi55_r[i] = 0;
      chi55_i[i] = 0;

      switch(parts[i].ncoeff) {
        case(21):
            p55_r[i] =   pnm_re[coeff_stride*i + 20];
            p55_i[i] =   pnm_im[coeff_stride*i + 20];
          phi55_r[i] = phinm_re[coeff_stride*i + 20];
          phi55_i[i] = phinm_im[coeff_stride*i + 20];
          chi55_r[i] = chinm_re[coeff_stride*i + 20];
          chi55_i[i] = chinm_im[coeff_stride*i + 20];
            p54_r[i] =   pnm_re[coeff_stride*i + 19];
            p54_i[i] =   pnm_im[coeff_stride*i + 19];
          phi54_r[i] = phinm_re[coeff_stride*i + 19];
          phi54_i[i] = phinm_im[coeff_stride*i + 19];
          chi54_r[i] = chinm_re[coeff_stride*i + 19];
          chi54_i[i] = chinm_im[coeff_stride*i + 19];
            p53_r[i] =   pnm_re[coeff_stride*i + 18];
            p53_i[i] =   pnm_im[coeff_stride*i + 18];
          phi53_r[i] = phinm_re[coeff_stride*i + 18];
          phi53_i[i] = phinm_im[coeff_stride*i + 18];
          chi53_r[i] = chinm_re[coeff_stride*i + 18];
          chi53_i[i] = chinm_im[coeff_stride*i + 18];
            p52_r[i] =   pnm_re[coeff_stride*i + 17];
            p52_i[i] =   pnm_im[coeff_stride*i + 17];
          phi52_r[i] = phinm_re[coeff_stride*i + 17];
          phi52_i[i] = phinm_im[coeff_stride*i + 17];
          chi52_r[i] = chinm_re[coeff_stride*i + 17];
          chi52_i[i] = chinm_im[coeff_stride*i + 17];
            p51_r[i] =   pnm_re[coeff_stride*i + 16];
            p51_i[i] =   pnm_im[coeff_stride*i + 16];
          phi51_r[i] = phinm_re[coeff_stride*i + 16];
          phi51_i[i] = phinm_im[coeff_stride*i + 16];
          chi51_r[i] = chinm_re[coeff_stride*i + 16];
          chi51_i[i] = chinm_im[coeff_stride*i + 16];
            p50_r[i] =   pnm_re[coeff_stride*i + 15];
            p50_i[i] =   pnm_im[coeff_stride*i + 15];
          phi50_r[i] = phinm_re[coeff_stride*i + 15];
          phi50_i[i] = phinm_im[coeff_stride*i + 15];
          chi50_r[i] = chinm_re[coeff_stride*i + 15];
          chi50_i[i] = chinm_im[coeff_stride*i + 15];
        case(15):
            p44_r[i] =   pnm_re[coeff_stride*i + 14];
            p44_i[i] =   pnm_im[coeff_stride*i + 14];
          phi44_r[i] = phinm_re[coeff_stride*i + 14];
          phi44_i[i] = phinm_im[coeff_stride*i + 14];
          chi44_r[i] = chinm_re[coeff_stride*i + 14];
          chi44_i[i] = chinm_im[coeff_stride*i + 14];
            p43_r[i] =   pnm_re[coeff_stride*i + 13];
            p43_i[i] =   pnm_im[coeff_stride*i + 13];
          phi43_r[i] = phinm_re[coeff_stride*i + 13];
          phi43_i[i] = phinm_im[coeff_stride*i + 13];
          chi43_r[i] = chinm_re[coeff_stride*i + 13];
          chi43_i[i] = chinm_im[coeff_stride*i + 13];
            p42_r[i] =   pnm_re[coeff_stride*i + 12];
            p42_i[i] =   pnm_im[coeff_stride*i + 12];
          phi42_r[i] = phinm_re[coeff_stride*i + 12];
          phi42_i[i] = phinm_im[coeff_stride*i + 12];
          chi42_r[i] = chinm_re[coeff_stride*i + 12];
          chi42_i[i] = chinm_im[coeff_stride*i + 12];
            p41_r[i] =   pnm_re[coeff_stride*i + 11];
            p41_i[i] =   pnm_im[coeff_stride*i + 11];
          phi41_r[i] = phinm_re[coeff_stride*i + 11];
          phi41_i[i] = phinm_im[coeff_stride*i + 11];
          chi41_r[i] = chinm_re[coeff_stride*i + 11];
          chi41_i[i] = chinm_im[coeff_stride*i + 11];
            p40_r[i] =   pnm_re[coeff_stride*i + 10];
            p40_i[i] =   pnm_im[coeff_stride*i + 10];
          phi40_r[i] = phinm_re[coeff_stride*i + 10];
          phi40_i[i] = phinm_im[coeff_stride*i + 10];
          chi40_r[i] = chinm_re[coeff_stride*i + 10];
          chi40_i[i] = chinm_im[coeff_stride*i + 10];
        case(10):
            p33_r[i] =   pnm_re[coeff_stride*i +  9];
            p33_i[i] =   pnm_im[coeff_stride*i +  9];
          phi33_r[i] = phinm_re[coeff_stride*i +  9];
          phi33_i[i] = phinm_im[coeff_stride*i +  9];
          chi33_r[i] = chinm_re[coeff_stride*i +  9];
          chi33_i[i] = chinm_im[coeff_stride*i +  9];
            p32_r[i] =   pnm_re[coeff_stride*i +  8];
            p32_i[i] =   pnm_im[coeff_stride*i +  8];
          phi32_r[i] = phinm_re[coeff_stride*i +  8];
          phi32_i[i] = phinm_im[coeff_stride*i +  8];
          chi32_r[i] = chinm_re[coeff_stride*i +  8];
          chi32_i[i] = chinm_im[coeff_stride*i +  8];
            p31_r[i] =   pnm_re[coeff_stride*i +  7];
            p31_i[i] =   pnm_im[coeff_stride*i +  7];
          phi31_r[i] = phinm_re[coeff_stride*i +  7];
          phi31_i[i] = phinm_im[coeff_stride*i +  7];
          chi31_r[i] = chinm_re[coeff_stride*i +  7];
          chi31_i[i] = chinm_im[coeff_stride*i +  7];
            p30_r[i] =   pnm_re[coeff_stride*i +  6];
            p30_i[i] =   pnm_im[coeff_stride*i +  6];
          phi30_r[i] = phinm_re[coeff_stride*i +  6];
          phi30_i[i] = phinm_im[coeff_stride*i +  6];
          chi30_r[i] = chinm_re[coeff_stride*i +  6];
          chi30_i[i] = chinm_im[coeff_stride*i +  6];
        case( 6):
            p22_r[i] =   pnm_re[coeff_stride*i +  5];
            p22_i[i] =   pnm_im[coeff_stride*i +  5];
          phi22_r[i] = phinm_re[coeff_stride*i +  5];
          phi22_i[i] = phinm_im[coeff_stride*i +  5];
          chi22_r[i] = chinm_re[coeff_stride*i +  5];
          chi22_i[i] = chinm_im[coeff_stride*i +  5];
            p21_r[i] =   pnm_re[coeff_stride*i +  4];
            p21_i[i] =   pnm_im[coeff_stride*i +  4];
          phi21_r[i] = phinm_re[coeff_stride*i +  4];
          phi21_i[i] = phinm_im[coeff_stride*i +  4];
          chi21_r[i] = chinm_re[coeff_stride*i +  4];
          chi21_i[i] = chinm_im[coeff_stride*i +  4];
            p20_r[i] =   pnm_re[coeff_stride*i +  3];
            p20_i[i] =   pnm_im[coeff_stride*i +  3];
          phi20_r[i] = phinm_re[coeff_stride*i +  3];
          phi20_i[i] = phinm_im[coeff_stride*i +  3];
          chi20_r[i] = chinm_re[coeff_stride*i +  3];
          chi20_i[i] = chinm_im[coeff_stride*i +  3];
        case( 3):
            p11_r[i] =   pnm_re[coeff_stride*i +  2];
            p11_i[i] =   pnm_im[coeff_stride*i +  2];
          phi11_r[i] = phinm_re[coeff_stride*i +  2];
          phi11_i[i] = phinm_im[coeff_stride*i +  2];
          chi11_r[i] = chinm_re[coeff_stride*i +  2];
          chi11_i[i] = chinm_im[coeff_stride*i +  2];
            p10_r[i] =   pnm_re[coeff_stride*i +  1];
            p10_i[i] =   pnm_im[coeff_stride*i +  1];
          phi10_r[i] = phinm_re[coeff_stride*i +  1];
          phi10_i[i] = phinm_im[coeff_stride*i +  1];
          chi10_r[i] = chinm_re[coeff_stride*i +  1];
          chi10_i[i] = chinm_im[coeff_stride*i +  1];
        case( 1):
            p00_r[i] =   pnm_re[coeff_stride*i +  0];
            p00_i[i] =   pnm_im[coeff_stride*i +  0];
          phi00_r[i] = phinm_re[coeff_stride*i +  0];
          phi00_i[i] = phinm_im[coeff_stride*i +  0];
          chi00_r[i] = chinm_re[coeff_stride*i +  0];
          chi00_i[i] = chinm_im[coeff_stride*i +  0];
      }
    }

    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &Xn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &Yn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &Zn);

    cg_section_write(fn, bn, zn, "Elements", NODE, 0, nparts-1, 0, conn, &en);

    cg_sol_write(fn, bn, zn, "Solution", Vertex, &sn);
    cg_field_write(fn, bn, zn, sn, RealDouble, "Radius", a, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "Density", rho, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "YoungsModulus", E, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "PoissonsRatio", sigma, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "DryCoeffRest", e_dry, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "FricCoeff", coeff_fric, &fnr);
    cg_field_write(fn, bn, zn, sn, Integer, "LambOrder", order, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", u, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", v, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", w, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AccelerationX", udot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AccelerationY", vdot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AccelerationZ", wdot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularVelocityX", ox, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularVelocityY", oy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularVelocityZ", oz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "HydroForceX", hFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "HydroForceY", hFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "HydroForceZ", hFz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SpringForceX", kFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SpringForceY", kFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SpringForceZ", kFz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceX", iFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceY", iFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceZ", iFz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceX", Fx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceY", Fy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceZ", Fz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentX", Lx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentY", Ly, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentZ", Lz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosXx", axx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosXy", axy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosXz", axz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosYx", ayx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosYy", ayy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosYz", ayz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosZx", azx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosZy", azy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularPosZz", azz, &fnr);
    

    switch(coeff_stride) {
      case(21):
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p55_re",   p55_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p55_im",   p55_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi55_re", phi55_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi55_im", phi55_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi55_re", chi55_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi55_im", chi55_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p54_re",   p54_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p54_im",   p54_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi54_re", phi54_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi54_im", phi54_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi54_re", chi54_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi54_im", chi54_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p53_re",   p53_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p53_im",   p53_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi53_re", phi53_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi53_im", phi53_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi53_re", chi53_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi53_im", chi53_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p52_re",   p52_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p52_im",   p52_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi52_re", phi52_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi52_im", phi52_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi52_re", chi52_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi52_im", chi52_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p51_re",   p51_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p51_im",   p51_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi51_re", phi51_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi51_im", phi51_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi51_re", chi51_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi51_im", chi51_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p50_re",   p50_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p50_im",   p50_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi50_re", phi50_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi50_im", phi50_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi50_re", chi50_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi50_im", chi50_i, &fnr);
      case(15):
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p44_re",   p44_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p44_im",   p44_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi44_re", phi44_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi44_im", phi44_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi44_re", chi44_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi44_im", chi44_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p43_re",   p43_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p43_im",   p43_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi43_re", phi43_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi43_im", phi43_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi43_re", chi43_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi43_im", chi43_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p42_re",   p42_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p42_im",   p42_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi42_re", phi42_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi42_im", phi42_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi42_re", chi42_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi42_im", chi42_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p41_re",   p41_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p41_im",   p41_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi41_re", phi41_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi41_im", phi41_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi41_re", chi41_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi41_im", chi41_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p40_re",   p40_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p40_im",   p40_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi40_re", phi40_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi40_im", phi40_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi40_re", chi40_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi40_im", chi40_i, &fnr);
      case(10):
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p33_re",   p33_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p33_im",   p33_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi33_re", phi33_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi33_im", phi33_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi33_re", chi33_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi33_im", chi33_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p32_re",   p32_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p32_im",   p32_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi32_re", phi32_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi32_im", phi32_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi32_re", chi32_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi32_im", chi32_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p31_re",   p31_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p31_im",   p31_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi31_re", phi31_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi31_im", phi31_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi31_re", chi31_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi31_im", chi31_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p30_re",   p30_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p30_im",   p30_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi30_re", phi30_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi30_im", phi30_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi30_re", chi30_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi30_im", chi30_i, &fnr);
      case( 6):
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p22_re",   p22_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p22_im",   p22_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi22_re", phi22_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi22_im", phi22_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi22_re", chi22_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi22_im", chi22_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p21_re",   p21_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p21_im",   p21_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi21_re", phi21_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi21_im", phi21_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi21_re", chi21_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi21_im", chi21_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p20_re",   p20_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p20_im",   p20_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi20_re", phi20_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi20_im", phi20_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi20_re", chi20_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi20_im", chi20_i, &fnr);
      case( 3):
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p11_re",   p11_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p11_im",   p11_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi11_re", phi11_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi11_im", phi11_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi11_re", chi11_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi11_im", chi11_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p10_re",   p10_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p10_im",   p10_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi10_re", phi10_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi10_im", phi10_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi10_re", chi10_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi10_im", chi10_i, &fnr);
      case( 1):
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p00_re",   p00_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble,   "p00_im",   p00_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi00_re", phi00_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "phi00_im", phi00_i, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi00_re", chi00_r, &fnr);
        cg_field_write(fn, bn, zn, sn, RealDouble, "chi00_im", chi00_i, &fnr);
    }

    cg_goto(fn, bn, "Zone_t", zn, "end");
    cg_user_data_write("Etc");
    cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
    cgsize_t *N = malloc(sizeof(cgsize_t));
    N[0] = 1;
    cg_array_write("Time", RealDouble, 1, N, &ttime);
    free(N);

    cg_close(fn);
    free(x);
    free(y);
    free(z);
    free(conn);
    free(a);
    free(rho);
    free(E);
    free(sigma);
    free(e_dry);
    free(coeff_fric);

    free(order);
    free(u);
    free(v);
    free(w);
    free(udot);
    free(vdot);
    free(wdot);
    free(kFx);
    free(kFy);
    free(kFz);
    free(iFx);
    free(iFy);
    free(iFz);
    free(iLx);
    free(iLy);
    free(iLz);
    free(hFx);
    free(hFy);
    free(hFz);
    free(hLx);
    free(hLy);
    free(hLz);
    free(ox);
    free(oy);
    free(oz);
    free(Fx);
    free(Fy);
    free(Fz);
    free(Lx);
    free(Ly);
    free(Lz);
    free(axx);
    free(axy);
    free(axz);
    free(ayx);
    free(ayy);
    free(ayz);
    free(azx);
    free(azy);
    free(azz);

    free(  p00_r);
    free(  p00_i);
    free(phi00_r);
    free(phi00_i);
    free(chi00_r);
    free(chi00_i);
    free(  p10_r);
    free(  p10_i);
    free(phi10_r);
    free(phi10_i);
    free(chi10_r);
    free(chi10_i);
    free(  p11_r);
    free(  p11_i);
    free(phi11_r);
    free(phi11_i);
    free(chi11_r);
    free(chi11_i);
    free(  p20_r);
    free(  p20_i);
    free(phi20_r);
    free(phi20_i);
    free(chi20_r);
    free(chi20_i);
    free(  p21_r);
    free(  p21_i);
    free(phi21_r);
    free(phi21_i);
    free(chi21_r);
    free(chi21_i);
    free(  p22_r);
    free(  p22_i);
    free(phi22_r);
    free(phi22_i);
    free(chi22_r);
    free(chi22_i);
    free(  p30_r);
    free(  p30_i);
    free(phi30_r);
    free(phi30_i);
    free(chi30_r);
    free(chi30_i);
    free(  p31_r);
    free(  p31_i);
    free(phi31_r);
    free(phi31_i);
    free(chi31_r);
    free(chi31_i);
    free(  p32_r);
    free(  p32_i);
    free(phi32_r);
    free(phi32_i);
    free(chi32_r);
    free(chi32_i);
    free(  p33_r);
    free(  p33_i);
    free(phi33_r);
    free(phi33_i);
    free(chi33_r);
    free(chi33_i);
    free(  p40_r);
    free(  p40_i);
    free(phi40_r);
    free(phi40_i);
    free(chi40_r);
    free(chi40_i);
    free(  p41_r);
    free(  p41_i);
    free(phi41_r);
    free(phi41_i);
    free(chi41_r);
    free(chi41_i);
    free(  p42_r);
    free(  p42_i);
    free(phi42_r);
    free(phi42_i);
    free(chi42_r);
    free(chi42_i);
    free(  p43_r);
    free(  p43_i);
    free(phi43_r);
    free(phi43_i);
    free(chi43_r);
    free(chi43_i);
    free(  p44_r);
    free(  p44_i);
    free(phi44_r);
    free(phi44_i);
    free(chi44_r);
    free(chi44_i);
    free(  p50_r);
    free(  p50_i);
    free(phi50_r);
    free(phi50_i);
    free(chi50_r);
    free(chi50_i);
    free(  p51_r);
    free(  p51_i);
    free(phi51_r);
    free(phi51_i);
    free(chi51_r);
    free(chi51_i);
    free(  p52_r);
    free(  p52_i);
    free(phi52_r);
    free(phi52_i);
    free(chi52_r);
    free(chi52_i);
    free(  p53_r);
    free(  p53_i);
    free(phi53_r);
    free(phi53_i);
    free(chi53_r);
    free(chi53_i);
    free(  p54_r);
    free(  p54_i);
    free(phi54_r);
    free(phi54_i);
    free(chi54_r);
    free(chi54_i);
    free(  p55_r);
    free(  p55_i);
    free(phi55_r);
    free(phi55_i);
    free(chi55_r);
    free(chi55_i);

  //  // now link this timestep into time series
  //  // create the time series file if it doesn't exist
  //  char sname[FILE_NAME_SIZE];
  //  sprintf(sname, "%s/output/part.cgns", ROOT_DIR);
  //  int fns;
  //  int bns;
  //  int depth;
  //  char *label = malloc(CHAR_BUF_SIZE * sizeof(real));
  //  cpumem += CHAR_BUF_SIZE * sizeof(real);
  //  char *bitername = malloc(CHAR_BUF_SIZE * sizeof(real));
  //  cpumem += CHAR_BUF_SIZE * sizeof(real);
  //  int index;
  //  int zns;
  //  real *tseries = malloc(sizeof(real));
  //  cpumem += sizeof(real);
  //  tseries[0] = 0;
  //  cgsize_t tsize[1];
  //  tsize[0] = 1;
  //
  //  if(cg_is_cgns(sname, &fns) != 0) {
  //    // if it does not exist, create it
  //    cg_open(sname, CG_MODE_WRITE, &fns);
  //    cg_base_write(fn, "Base", 3, 3, &bns);
  //    size[0][0] = nparts;
  //    size[1][0] = 0;
  //    size[2][0] = 0;
  //    cg_zone_write(fns, bns, "Zone0", size[0], Unstructured, &zns);
  //    cg_goto(fns, bns, "Zone_t", zns, "end");
  //    char gridname2[CHAR_BUF_SIZE];
  //    sprintf(gridname2, "GridCoordinates-%s", format);
  //    sprintf(gridname2, gridname2, tout);
  //    cg_link_write(gridname2, fname, "Base/Zone0/GridCoordinates");
  //    char solname2[CHAR_BUF_SIZE];
  //    sprintf(solname2, "Solution-%s", format);
  //    sprintf(solname2, solname2, tout);
  //    cg_link_write(solname2, fname, "Base/Zone0/Solution");
  //    cg_biter_write(fns, bns, "BaseIterativeData", 1);
  //    cg_goto(fns, bns, "BaseIterativeData_t", 1, "end");
  //    int N = 1;
  //    cg_array_write("TimeValues", RealDouble, 1, &N, &ttime);
  //    cg_ziter_write(fns, bns, zns, "ZoneIterativeData");
  //    cg_goto(fns, bns, "Zone_t", zns, "ZoneIterativeData_t", 1, "end");
  //    cg_simulation_type_write(fns, bns, TimeAccurate);
  //    cg_close(fns);
  //  } else {
  //    // modify the timeseries file
  //    cg_open(sname, CG_MODE_MODIFY, &fns);
  //    cg_gopath(fns, "/Base");
  //    cg_where(&fns, &bns, &depth, &label, &index);
  //    int N = 0;
  //    cg_biter_read(fns, bns, bitername, &N);
  //    free(tseries);
  //    tseries = malloc((N+1) * sizeof(real));
  //    cpumem += (N+1) * sizeof(real);
  //    cg_goto(fns, bns, "BaseIterativeData_t", 1, "end");
  //    cg_array_read(1, tseries);
  //    cg_biter_write(fns, bns, "BaseIterativeData", N+1);
  //    cg_goto(fns, bns, "BaseIterativeData_t", 1, "end");
  //    tseries[N] = ttime;
  //    tsize[0] = N+1;
  //    cg_array_write("TimeValues", RealDouble, 1, tsize, tseries);
  //    cg_goto(fns, bns, "Zone_t", 1, "end");
  //    char gridname2[CHAR_BUF_SIZE];
  //    sprintf(gridname2, "GridCoordinates-%s", format);
  //    sprintf(gridname2, gridname2, tout);
  //    cg_link_write(gridname2, fname, "Base/Zone0/GridCoordinates");
  //    char solname2[CHAR_BUF_SIZE];
  //    sprintf(solname2, "Solution-%s", format);
  //    sprintf(solname2, solname2, tout);
  //    cg_link_write(solname2, fname, "Base/Zone0/Solution");
  //    cg_close(fns);
  //   
  //    free(tseries);
  //    free(label);
  //    free(bitername);
  //  }
  }
}

void recorder_bicgstab_init(char *name)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  fprintf(rec, "%-12s", "stepnum");
  fprintf(rec, "%-15s", "ttime");
  fprintf(rec, "%-15s", "dt");
  fprintf(rec, "%-8s", "niter");
  fprintf(rec, "%-15s", "resid");

  // close the file
  fclose(rec);
}

void recorder_lamb_init(char *name)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  // close the file
  fclose(rec);
}

void recorder_bicgstab(char *name, int niter, real resid)
{
  // open the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");
  if(rec == NULL) {
    recorder_bicgstab_init(name);
    rec = fopen(path, "r+");
    #ifdef IMPLICIT
      recorder_bicgstab_init(name);
    #endif
  }

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "\n");
  fprintf(rec, "%-12d", stepnum);
  fprintf(rec, "%-15e", ttime);
  fprintf(rec, "%-15e", dt);
  fprintf(rec, "%-8d", niter);
  fprintf(rec, "%-15e", resid);

  // close the file
  fclose(rec);
}

void recorder_lamb(char *name, int iter)
{
  // open the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");
  if(rec == NULL) {
    recorder_lamb_init(name);
    rec = fopen(path, "r+");
  }

  int n = 0;
  int m = 0;
  int c = 0;

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "ttime = %e, iter = %d\n", ttime, iter);
  for(int i = 0; i < nparts; i++) {
    fprintf(rec, "particle %d\n", i);
    fprintf(rec, "%4s%4s%12s%12s", "n", "m", "pnm_re", "pnm_im");
    fprintf(rec, "%12s%12s", "phinm_re", "phinm_im");
    fprintf(rec, "%12s%12s\n", "chinm_re", "chinm_im");
    n = 0;
    m = 0;
    c = 0;
    while(c < parts[i].ncoeff) {
      fprintf(rec, "%4d%4d", n, m);
      fprintf(rec, "%12.3e%12.3e",
        pnm_re[coeff_stride*i+c], pnm_im[coeff_stride*i+c]);
      fprintf(rec, "%12.3e%12.3e",
        phinm_re[coeff_stride*i+c], phinm_im[coeff_stride*i+c]);
      fprintf(rec, "%12.3e%12.3e\n",
        chinm_re[coeff_stride*i+c], chinm_im[coeff_stride*i+c]);
      m++;
      c++;
      if(m > n) {
        n++;
        m = 0;
      }
    }
  }

  // close the file
  fclose(rec);
}
