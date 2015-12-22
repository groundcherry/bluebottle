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

int nsubdom;

void domain_read_input(void)
{
  int i;  // iterator

  int fret = 0;
  fret = fret; // prevent compiler warning

  cpumem = 0;
  gpumem = 0;

  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/flow.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  char buf[CHAR_BUF_SIZE] = "";  // character read buffer

  // read domain
  fret = fscanf(infile, "DOMAIN\n");
#ifdef DOUBLE
  fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &Dom.xs, &Dom.xe, &Dom.xn);
  fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &Dom.ys, &Dom.ye, &Dom.yn);
  fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &Dom.zs, &Dom.ze, &Dom.zn);
  fret = fscanf(infile, "\n"); // \n
#else // single precision
  fret = fscanf(infile, "(Xs, Xe, Xn) %f %f %d\n", &Dom.xs, &Dom.xe, &Dom.xn);
  fret = fscanf(infile, "(Ys, Ye, Yn) %f %f %d\n", &Dom.ys, &Dom.ye, &Dom.yn);
  fret = fscanf(infile, "(Zs, Ze, Zn) %f %f %d\n", &Dom.zs, &Dom.ze, &Dom.zn);
  fret = fscanf(infile, "\n"); // \n
#endif
  fret = fscanf(infile, "GPU DOMAIN DECOMPOSITION\n");
  fret = fscanf(infile, "DEV RANGE %d %d\n", &dev_start, &dev_end);
  fret = fscanf(infile, "n %d\n", &nsubdom);
  // allocate subdomain data structure
  dom = (dom_struct*) malloc(nsubdom * sizeof(dom_struct));
  cpumem += nsubdom * sizeof(dom_struct);
  for(i = 0; i < nsubdom; i++) {  // read subdomains
#ifdef DOUBLE
    fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &dom[i].xs, &dom[i].xe,
      &dom[i].xn);
    fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &dom[i].ys, &dom[i].ye,
      &dom[i].yn);
    fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &dom[i].zs, &dom[i].ze,
      &dom[i].zn);
#else // single
    fret = fscanf(infile, "(Xs, Xe, Xn) %f %f %d\n", &dom[i].xs, &dom[i].xe,
      &dom[i].xn);
    fret = fscanf(infile, "(Ys, Ye, Yn) %f %f %d\n", &dom[i].ys, &dom[i].ye,
      &dom[i].yn);
    fret = fscanf(infile, "(Zs, Ze, Zn) %f %f %d\n", &dom[i].zs, &dom[i].ze,
      &dom[i].zn);
#endif
    fret = fscanf(infile, "E %d W %d N %d S %d T %d B %d\n", &dom[i].E, &dom[i].W,
      &dom[i].N, &dom[i].S, &dom[i].T, &dom[i].B);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "PHYSICAL PARAMETERS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "rho_f %lf\n", &rho_f);
  fret = fscanf(infile, "nu %lf\n", &nu);
#else // single
  fret = fscanf(infile, "rho_f %f\n", &rho_f);
  fret = fscanf(infile, "nu %f\n", &nu);
#endif
  mu = nu * rho_f;  // set dynamic viscosity
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SIMULATION PARAMETERS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "duration %lf\n", &duration);
  fret = fscanf(infile, "CFL %lf\n", &CFL);
  fret = fscanf(infile, "pp_max_iter %d\n", &pp_max_iter);
  fret = fscanf(infile, "pp_residual %lf\n", &pp_residual);
  fret = fscanf(infile, "lamb_max_iter %d\n", &lamb_max_iter);
  fret = fscanf(infile, "lamb_residual %lf\n", &lamb_residual);
  fret = fscanf(infile, "lamb_relax %lf\n", &lamb_relax);
  fret = fscanf(infile, "lamb_cut %lf\n", &lamb_cut);
#else
  fret = fscanf(infile, "duration %f\n", &duration);
  fret = fscanf(infile, "CFL %f\n", &CFL);
  fret = fscanf(infile, "pp_max_iter %d\n", &pp_max_iter);
  fret = fscanf(infile, "pp_residual %f\n", &pp_residual);
  fret = fscanf(infile, "lamb_max_iter %d\n", &lamb_max_iter);
  fret = fscanf(infile, "lamb_residual %f\n", &lamb_residual);
  fret = fscanf(infile, "lamb_relax %f\n", &lamb_relax);
  fret = fscanf(infile, "lamb_cut %f\n", &lamb_cut);
#endif
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "BOUNDARY CONDITIONS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "vel_tDelay %lf\n", &vel_tDelay);
#else  
  fret = fscanf(infile, "vel_tDelay %f\n", &vel_tDelay);
#endif
  fret = fscanf(infile, "PRESSURE\n");
  fret = fscanf(infile, "bc.pW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pW = PERIODIC;
    bc.dsW = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pW = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsW);
  } else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pE = PERIODIC;
    bc.dsE = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pE = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsE);
  } else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pS = PERIODIC;
    bc.dsS = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pS = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsS);
  } else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pN = PERIODIC;
    bc.dsN = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pN = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsN);
  } else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pB = PERIODIC;
    bc.dsB = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pB = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsB);
  } else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pT = PERIODIC;
    bc.dsT = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pT = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsT);
  } else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  if (bc.dsE < 0 || bc.dsW < 0 || 
      bc.dsN < 0 || bc.dsS < 0 || 
      bc.dsT < 0 || bc.dsB < 0) {
    fprintf(stderr, "flow.config read error -- screen offsets must be >= 0\n");
    exit(EXIT_FAILURE);
  }

  fret = fscanf(infile, "X-VELOCITY\n");

  fret = fscanf(infile, "bc.uW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uW = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uWDm, &bc.uWDa);
    bc.uWD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uW = NEUMANN;
    bc.uWDm = 0;
    bc.uWDa = 0;
    bc.uWD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.uW = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.uWDm, &bc.uWDa);
    bc.uWD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uE = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uEDm, &bc.uEDa);
    bc.uED = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uE = NEUMANN;
    bc.uEDm = 0;
    bc.uEDa = 0;
    bc.uED = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.uE = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.uEDm, &bc.uEDa);
    bc.uED = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uS = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uSDm, &bc.uSDa);
    bc.uSD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uS = NEUMANN;
    bc.uSDm = 0;
    bc.uSDa = 0;
    bc.uSD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.uS = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.uSDm, &bc.uSDa);
    bc.uSD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uN = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uNDm, &bc.uNDa);
    bc.uND = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uN = NEUMANN;
    bc.uNDm = 0;
    bc.uNDa = 0;
    bc.uND = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.uN = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.uNDm, &bc.uNDa);
    bc.uND = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uB = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uBDm, &bc.uBDa);
    bc.uBD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uB = NEUMANN;
    bc.uBDm = 0;
    bc.uBDa = 0;
    bc.uBD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.uB = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.uBDm, &bc.uBDa);
    bc.uBD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uT = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uTDm, &bc.uTDa);
    bc.uTD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uT = NEUMANN;
    bc.uTDm = 0;
    bc.uTDa = 0;
    bc.uTD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.uT = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.uTDm, &bc.uTDa);
    bc.uTD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "Y-VELOCITY\n");

  fret = fscanf(infile, "bc.vW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vW = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vWDm, &bc.vWDa);
    bc.vWD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vW = NEUMANN;
    bc.vWDm = 0;
    bc.vWDa = 0;
    bc.vWD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.vW = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.vWDm, &bc.vWDa);
    bc.vWD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vE = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vEDm, &bc.vEDa);
    bc.vED = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vE = NEUMANN;
    bc.vEDm = 0;
    bc.vEDa = 0;
    bc.vED = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.vE = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.vEDm, &bc.vEDa);
    bc.vED = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vS = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vSDm, &bc.vSDa);
    bc.vSD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vS = NEUMANN;
    bc.vSDm = 0;
    bc.vSDa = 0;
    bc.vSD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.vS = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.vSDm, &bc.vSDa);
    bc.vSD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vN = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vNDm, &bc.vNDa);
    bc.vND = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vN = NEUMANN;
    bc.vNDm = 0;
    bc.vNDa = 0;
    bc.vND = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.vN = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.vNDm, &bc.vNDa);
    bc.vND = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vB = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vBDm, &bc.vBDa);
    bc.vBD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vB = NEUMANN;
    bc.vBDm = 0;
    bc.vBDa = 0;
    bc.vBD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.vB = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.vBDm, &bc.vBDa);
    bc.vBD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vT = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vTDm, &bc.vTDa);
    bc.vTD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vT = NEUMANN;
    bc.vTDm = 0;
    bc.vTDa = 0;
    bc.vTD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.vT = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.vTDm, &bc.vTDa);
    bc.vTD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "Z-VELOCITY\n");

  fret = fscanf(infile, "bc.wW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wW = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wWDm, &bc.wWDa);
    bc.wWD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wW = NEUMANN;
    bc.wWDm = 0;
    bc.wWDa = 0;
    bc.wWD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.wW = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.wWDm, &bc.wWDa);
    bc.wWD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wE = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wEDm, & bc.wEDa);
    bc.wED = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wE = NEUMANN;
    bc.wEDm = 0;
    bc.wEDa = 0;
    bc.wED = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.wE = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.wEDm, &bc.wEDa);
    bc.wED = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wS = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wSDm, &bc.wSDa);
    bc.wSD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wS = NEUMANN;
    bc.wSDm = 0;
    bc.wSDa = 0;
    bc.wSD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.wS = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.wSDm, &bc.wSDa);
    bc.wSD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wN = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wNDm, &bc.wNDa);
    bc.wND = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wN = NEUMANN;
    bc.wNDm = 0;
    bc.wNDa = 0;
    bc.wND = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.wN = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.wNDm, &bc.wNDa);
    bc.wND = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wB = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wBDm, &bc.wBDa);
    bc.wBD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wB = NEUMANN;
    bc.wBDm = 0;
    bc.wBDa = 0;
    bc.wBD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.wB = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.wBDm, &bc.wBDa);
    bc.wBD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wT = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wTDm, &bc.wTDa);
    bc.wTD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wT = NEUMANN;
    bc.wTDm = 0;
    bc.wTDa = 0;
    bc.wTD = 0;
  } else if(strcmp(buf, "PRECURSOR") == 0) {
    bc.wT = PRECURSOR;
    fret = fscanf(infile, "%lf %lf", &bc.wTDm, &bc.wTDa);
    bc.wTD = 0;
  }
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "INITIAL CONDITION\n");
  fret = fscanf(infile, "init_cond %s", buf);
  if(strcmp(buf, "QUIESCENT") == 0)
    init_cond = QUIESCENT;
  else if(strcmp(buf, "SHEAR") == 0)
    init_cond = SHEAR;
  else if(strcmp(buf, "CHANNEL") == 0)
    init_cond = CHANNEL;
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SOLVABILITY ENFORCEMENT PLANE\n");
  fret = fscanf(infile, "out_plane %s", buf);
  if(strcmp(buf, "WEST") == 0)
    out_plane = WEST;
  else if(strcmp(buf, "EAST") == 0)
    out_plane = EAST;
  else if(strcmp(buf, "SOUTH") == 0)
    out_plane = SOUTH;
  else if(strcmp(buf, "NORTH") == 0)
    out_plane = NORTH;
  else if(strcmp(buf, "BOTTOM") == 0)
    out_plane = BOTTOM;
  else if(strcmp(buf, "TOP") == 0)
    out_plane = TOP;
  else if(strcmp(buf, "HOMOGENEOUS") == 0)
    out_plane = HOMOGENEOUS;
  else {
    fprintf(stderr, "flow.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SIMULATION DRIVING CONDITIONS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "p_tDelay %lf\n", &p_tDelay);
#else  
  fret = fscanf(infile, "p_tDelay %f\n", &p_tDelay);
#endif
  turbA = 0.; // prevent turbulence linear forcing from being used
#ifdef DOUBLE
  fret = fscanf(infile, "gradP.x %lf %lf\n", &gradP.xm, &gradP.xa);
  gradP.x = 0;
  fret = fscanf(infile, "gradP.y %lf %lf\n", &gradP.ym, &gradP.ya);
  gradP.y = 0;
  fret = fscanf(infile, "gradP.z %lf %lf\n", &gradP.zm, &gradP.za);
  gradP.z = 0;
#ifdef DOUBLE
  fret = fscanf(infile, "g_tDelay %lf\n", &g_tDelay);
#else  
  fret = fscanf(infile, "g_tDelay %f\n", &g_tDelay);
#endif
  fret = fscanf(infile, "g.x %lf %lf\n", &g.xm, &g.xa);
  g.x = 0;
  fret = fscanf(infile, "g.y %lf %lf\n", &g.ym, &g.ya);
  g.y = 0;
  fret = fscanf(infile, "g.z %lf %lf\n", &g.zm, &g.za);
  g.z = 0;
#else
  fret = fscanf(infile, "gradP.x %f %f\n", &gradP.xm, &gradP.xa);
  gradP.x = 0;
  fret = fscanf(infile, "gradP.y %f %f\n", &gradP.ym, &gradP.ya);
  gradP.y = 0;
  fret = fscanf(infile, "gradP.z %f %f\n", &gradP.zm, &gradP.za);
  gradP.z = 0;
  fret = fscanf(infile, "g.x %f %f\n", &g.xm, &g.xa);
  g.x = 0;
  fret = fscanf(infile, "g.y %f %f\n", &g.ym, &g.ya);
  g.y = 0;
  fret = fscanf(infile, "g.z %f %f\n", &g.zm, &g.za);
  g.z = 0;
#endif
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "PID CONTROLLER GAINS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "Kp %lf\n", &Kp);
  fret = fscanf(infile, "Ki %lf\n", &Ki);
  fret = fscanf(infile, "Kd %lf\n", &Kd);
#else
  fret = fscanf(infile, "Kp %f\n", &Kp);
  fret = fscanf(infile, "Ki %f\n", &Ki);
  fret = fscanf(infile, "Kd %f\n", &Kd);
#endif

  fclose(infile);
}

void turb_read_input(void)
{
  int i;  // iterator

  int fret = 0;
  fret = fret; // prevent compiler warning

  cpumem = 0;
  gpumem = 0;

  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/turb.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  char buf[CHAR_BUF_SIZE] = "";  // character read buffer

  // read domain
  fret = fscanf(infile, "DOMAIN\n");
#ifdef DOUBLE
  fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &Dom.xs, &Dom.xe, &Dom.xn);
  fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &Dom.ys, &Dom.ye, &Dom.yn);
  fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &Dom.zs, &Dom.ze, &Dom.zn);
  fret = fscanf(infile, "\n"); // \n
#else // single precision
  fret = fscanf(infile, "(Xs, Xe, Xn) %f %f %d\n", &Dom.xs, &Dom.xe, &Dom.xn);
  fret = fscanf(infile, "(Ys, Ye, Yn) %f %f %d\n", &Dom.ys, &Dom.ye, &Dom.yn);
  fret = fscanf(infile, "(Zs, Ze, Zn) %f %f %d\n", &Dom.zs, &Dom.ze, &Dom.zn);
  fret = fscanf(infile, "\n"); // \n
#endif
  fret = fscanf(infile, "GPU DOMAIN DECOMPOSITION\n");
  fret = fscanf(infile, "DEV RANGE %d %d\n", &dev_start, &dev_end);
  fret = fscanf(infile, "n %d\n", &nsubdom);
  // allocate subdomain data structure
  dom = (dom_struct*) malloc(nsubdom * sizeof(dom_struct));
  cpumem += nsubdom * sizeof(dom_struct);
  for(i = 0; i < nsubdom; i++) {  // read subdomains
#ifdef DOUBLE
    fret = fscanf(infile, "(Xs, Xe, Xn) %lf %lf %d\n", &dom[i].xs, &dom[i].xe,
      &dom[i].xn);
    fret = fscanf(infile, "(Ys, Ye, Yn) %lf %lf %d\n", &dom[i].ys, &dom[i].ye,
      &dom[i].yn);
    fret = fscanf(infile, "(Zs, Ze, Zn) %lf %lf %d\n", &dom[i].zs, &dom[i].ze,
      &dom[i].zn);
#else // single
    fret = fscanf(infile, "(Xs, Xe, Xn) %f %f %d\n", &dom[i].xs, &dom[i].xe,
      &dom[i].xn);
    fret = fscanf(infile, "(Ys, Ye, Yn) %f %f %d\n", &dom[i].ys, &dom[i].ye,
      &dom[i].yn);
    fret = fscanf(infile, "(Zs, Ze, Zn) %f %f %d\n", &dom[i].zs, &dom[i].ze,
      &dom[i].zn);
#endif
    fret = fscanf(infile, "E %d W %d N %d S %d T %d B %d\n", &dom[i].E, &dom[i].W,
      &dom[i].N, &dom[i].S, &dom[i].T, &dom[i].B);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "PHYSICAL PARAMETERS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "rho_f %lf\n", &rho_f);
  fret = fscanf(infile, "nu %lf\n", &nu);
#else // single
  fret = fscanf(infile, "rho_f %f\n", &rho_f);
  fret = fscanf(infile, "nu %f\n", &nu);
#endif
  mu = nu * rho_f;  // set dynamic viscosity
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SIMULATION PARAMETERS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "duration %lf\n", &duration);
  fret = fscanf(infile, "CFL %lf\n", &CFL);
  fret = fscanf(infile, "pp_max_iter %d\n", &pp_max_iter);
  fret = fscanf(infile, "pp_residual %lf\n", &pp_residual);
#else
  fret = fscanf(infile, "duration %f\n", &duration);
  fret = fscanf(infile, "CFL %f\n", &CFL);
  fret = fscanf(infile, "pp_max_iter %d\n", &pp_max_iter);
  fret = fscanf(infile, "pp_residual %f\n", &pp_residual);
#endif
  // these are unnecessary in the precursor
  lamb_max_iter = 0.;
  lamb_residual = 0.;
  lamb_relax = 0.;
  lamb_cut = 0.;
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "BOUNDARY CONDITIONS\n");
  vel_tDelay = 0;
  fret = fscanf(infile, "PRESSURE\n");
  fret = fscanf(infile, "bc.pW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pW = PERIODIC;
    bc.dsW = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pW = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsW);
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pE = PERIODIC;
    bc.dsE = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pE = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsE);
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pS = PERIODIC;
    bc.dsS = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pS = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsS);
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pN = PERIODIC;
    bc.dsN = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pN = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsN);
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pB = PERIODIC;
    bc.dsB = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pB = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsB);
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.pT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0) {
    bc.pT = PERIODIC;
    bc.dsT = 0;
  } else if(strcmp(buf, "NEUMANN") == 0) {
    bc.pT = NEUMANN;
    fret = fscanf(infile, "%lf", &bc.dsT);
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "X-VELOCITY\n");

  fret = fscanf(infile, "bc.uW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uW = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uWDm, &bc.uWDa);
    bc.uWD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uW = NEUMANN;
    bc.uWDm = 0;
    bc.uWDa = 0;
    bc.uWD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uE = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uEDm, &bc.uEDa);
    bc.uED = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uE = NEUMANN;
    bc.uEDm = 0;
    bc.uEDa = 0;
    bc.uED = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uS = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uSDm, &bc.uSDa);
    bc.uSD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uS = NEUMANN;
    bc.uSDm = 0;
    bc.uSDa = 0;
    bc.uSD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uN = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uNDm, &bc.uNDa);
    bc.uND = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uN = NEUMANN;
    bc.uNDm = 0;
    bc.uNDa = 0;
    bc.uND = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uB = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uBDm, &bc.uBDa);
    bc.uBD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uB = NEUMANN;
    bc.uBDm = 0;
    bc.uBDa = 0;
    bc.uBD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.uT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.uT = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.uT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.uTDm, &bc.uTDa);
    bc.uTD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.uT = NEUMANN;
    bc.uTDm = 0;
    bc.uTDa = 0;
    bc.uTD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "Y-VELOCITY\n");

  fret = fscanf(infile, "bc.vW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vW = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vWDm, &bc.vWDa);
    bc.vWD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vW = NEUMANN;
    bc.vWDm = 0;
    bc.vWDa = 0;
    bc.vWD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vE = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vEDm, &bc.vEDa);
    bc.vED = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vE = NEUMANN;
    bc.vEDm = 0;
    bc.vEDa = 0;
    bc.vED = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vS = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vSDm, &bc.vSDa);
    bc.vSD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vS = NEUMANN;
    bc.vSDm = 0;
    bc.vSDa = 0;
    bc.vSD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vN = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vNDm, &bc.vNDa);
    bc.vND = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vN = NEUMANN;
    bc.vNDm = 0;
    bc.vNDa = 0;
    bc.vND = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vB = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vBDm, &bc.vBDa);
    bc.vBD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vB = NEUMANN;
    bc.vBDm = 0;
    bc.vBDa = 0;
    bc.vBD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.vT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.vT = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.vT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.vTDm, &bc.vTDa);
    bc.vTD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.vT = NEUMANN;
    bc.vTDm = 0;
    bc.vTDa = 0;
    bc.vTD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "Z-VELOCITY\n");

  fret = fscanf(infile, "bc.wW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wW = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wWDm, &bc.wWDa);
    bc.wWD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wW = NEUMANN;
    bc.wWDm = 0;
    bc.wWDa = 0;
    bc.wWD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wE = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wEDm, & bc.wEDa);
    bc.wED = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wE = NEUMANN;
    bc.wEDm = 0;
    bc.wEDa = 0;
    bc.wED = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wS = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wSDm, &bc.wSDa);
    bc.wSD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wS = NEUMANN;
    bc.wSDm = 0;
    bc.wSDa = 0;
    bc.wSD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wN = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wNDm, &bc.wNDa);
    bc.wND = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wN = NEUMANN;
    bc.wNDm = 0;
    bc.wNDa = 0;
    bc.wND = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wB = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wBDm, &bc.wBDa);
    bc.wBD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wB = NEUMANN;
    bc.wBDm = 0;
    bc.wBDa = 0;
    bc.wBD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "bc.wT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    bc.wT = PERIODIC;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    bc.wT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &bc.wTDm, &bc.wTDa);
    bc.wTD = 0;
  }
  else if(strcmp(buf, "NEUMANN") == 0) {
    bc.wT = NEUMANN;
    bc.wTDm = 0;
    bc.wTDa = 0;
    bc.wTD = 0;
  } else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SOLVABILITY ENFORCEMENT PLANE\n");
  fret = fscanf(infile, "out_plane %s", buf);
  if(strcmp(buf, "WEST") == 0)
    out_plane = WEST;
  else if(strcmp(buf, "EAST") == 0)
    out_plane = EAST;
  else if(strcmp(buf, "SOUTH") == 0)
    out_plane = SOUTH;
  else if(strcmp(buf, "NORTH") == 0)
    out_plane = NORTH;
  else if(strcmp(buf, "BOTTOM") == 0)
    out_plane = BOTTOM;
  else if(strcmp(buf, "TOP") == 0)
    out_plane = TOP;
  else if(strcmp(buf, "HOMOGENEOUS") == 0)
    out_plane = HOMOGENEOUS;
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "SIMULATION DRIVING CONDITIONS\n");
#ifdef DOUBLE
  fret = fscanf(infile, "turbA %lf\n", &turbA);
#else
  fret = fscanf(infile, "turbA %f\n", &turbA);
#endif
  // these are unnecessary in the precursor
  p_tDelay = 0;
  gradP.x = 0.;
  gradP.xm = 0.;
  gradP.xa = 0.;
  gradP.y = 0.;
  gradP.ym = 0.;
  gradP.ya = 0.;
  gradP.z = 0.;
  gradP.zm = 0.;
  gradP.za = 0.;
  g_tDelay = 0; 
  g.x = 0.;
  g.xm = 0.;
  g.xa = 0.;
  g.y = 0.;
  g.ym = 0.;
  g.ya = 0.;
  g.z = 0.;
  g.zm = 0.;
  g.za = 0.;

  fclose(infile);
}

void domain_show_config(void)
{
  int i;  // iterator

  printf("Domain:\n");
  printf("  X: (%f, %f), dX = %f\n", Dom.xs, Dom.xe, Dom.dx);
  printf("  Y: (%f, %f), dY = %f\n", Dom.ys, Dom.ye, Dom.dy);
  printf("  Z: (%f, %f), dZ = %f\n", Dom.zs, Dom.ze, Dom.dz);
  printf("  Xn = %d, Yn = %d, Zn = %d\n", Dom.xn, Dom.yn, Dom.zn);
  printf("Domain Grids:\n");
  printf("  Dom.Gcc:\n");
  printf("    is = %d, ie = %d, in = %d\n", Dom.Gcc.is, Dom.Gcc.ie, Dom.Gcc.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", Dom.Gcc.isb, Dom.Gcc.ieb,
    Dom.Gcc.inb);
  printf("    js = %d, je = %d, jn = %d\n", Dom.Gcc.js, Dom.Gcc.je, Dom.Gcc.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", Dom.Gcc.jsb, Dom.Gcc.jeb,
    Dom.Gcc.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", Dom.Gcc.ks, Dom.Gcc.ke, Dom.Gcc.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", Dom.Gcc.ksb, Dom.Gcc.keb,
    Dom.Gcc.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", Dom.Gcc.s1, Dom.Gcc.s2,
    Dom.Gcc.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", Dom.Gcc.s1b, Dom.Gcc.s2b,
    Dom.Gcc.s3b);
  printf("  Dom.Gfx:\n");
  printf("    is = %d, ie = %d, in = %d\n", Dom.Gfx.is, Dom.Gfx.ie, Dom.Gfx.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", Dom.Gfx.isb, Dom.Gfx.ieb,
    Dom.Gfx.inb);
  printf("    js = %d, je = %d, jn = %d\n", Dom.Gfx.js, Dom.Gfx.je, Dom.Gfx.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", Dom.Gfx.jsb, Dom.Gfx.jeb,
    Dom.Gfx.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", Dom.Gfx.ks, Dom.Gfx.ke, Dom.Gfx.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", Dom.Gfx.ksb, Dom.Gfx.keb,
    Dom.Gfx.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", Dom.Gfx.s1, Dom.Gfx.s2,
    Dom.Gfx.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", Dom.Gfx.s1b, Dom.Gfx.s2b,
    Dom.Gfx.s3b);
  printf("  Dom.Gfy:\n");
  printf("    is = %d, ie = %d, in = %d\n", Dom.Gfy.is, Dom.Gfy.ie, Dom.Gfy.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", Dom.Gfy.isb, Dom.Gfy.ieb,
    Dom.Gfy.inb);
  printf("    js = %d, je = %d, jn = %d\n", Dom.Gfy.js, Dom.Gfy.je, Dom.Gfy.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", Dom.Gfy.jsb, Dom.Gfy.jeb,
    Dom.Gfy.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", Dom.Gfy.ks, Dom.Gfy.ke, Dom.Gfy.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", Dom.Gfy.ksb, Dom.Gfy.keb,
    Dom.Gfy.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", Dom.Gfy.s1, Dom.Gfy.s2,
    Dom.Gfy.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", Dom.Gfy.s1b, Dom.Gfy.s2b,
    Dom.Gfy.s3b);
  printf("  Dom.Gfz:\n");
  printf("    is = %d, ie = %d, in = %d\n", Dom.Gfz.is, Dom.Gfz.ie, Dom.Gfz.in);
  printf("    isb = %d, ieb = %d, inb = %d\n", Dom.Gfz.isb, Dom.Gfz.ieb,
    Dom.Gfz.inb);
  printf("    js = %d, je = %d, jn = %d\n", Dom.Gfz.js, Dom.Gfz.je, Dom.Gfz.jn);
  printf("    jsb = %d, jeb = %d, jnb = %d\n", Dom.Gfz.jsb, Dom.Gfz.jeb,
    Dom.Gfz.jnb);
  printf("    ks = %d, ke = %d, kn = %d\n", Dom.Gfz.ks, Dom.Gfz.ke, Dom.Gfz.kn);
  printf("    ksb = %d, keb = %d, knb = %d\n", Dom.Gfz.ksb, Dom.Gfz.keb,
    Dom.Gfz.knb);
  printf("    s1 = %d, s2 = %d, s3 = %d\n", Dom.Gfz.s1, Dom.Gfz.s2,
    Dom.Gfz.s3);
  printf("    s1b = %d, s2b = %d, s3b = %d\n", Dom.Gfz.s1b, Dom.Gfz.s2b,
    Dom.Gfz.s3b);

  printf("GPU Domain Decomposition:\n");
  printf("  nsubdom = %d\n", nsubdom);
  for(i = 0; i < nsubdom; i++) {
    printf("Subdomain %d:\n", i);
    printf("  X: (%lf, %lf), dX = %f\n", dom[i].xs, dom[i].xe, dom[i].dx);
    printf("  Y: (%lf, %lf), dY = %f\n", dom[i].ys, dom[i].ye, dom[i].dy);
    printf("  Z: (%lf, %lf), dZ = %f\n", dom[i].zs, dom[i].ze, dom[i].dz);
    printf("  Xn = %d, Yn = %d, Zn = %d\n", dom[i].xn, dom[i].yn, dom[i].zn);
    printf("Subdomain connectivity:\n");
    printf("  E: %d, W: %d, N: %d, S: %d, T: %d, B: %d\n", dom[i].E, dom[i].W,
      dom[i].N, dom[i].S, dom[i].T, dom[i].B);
    printf("Subdomain %d Grids:\n", i);
    printf("  dom[%d].Gcc:\n", i);
    printf("    is = %d, ie = %d, in = %d\n", dom[i].Gcc.is, dom[i].Gcc.ie, dom[i].Gcc.in);
    printf("    isb = %d, ieb = %d, inb = %d\n", dom[i].Gcc.isb, dom[i].Gcc.ieb,
      dom[i].Gcc.inb);
    printf("    js = %d, je = %d, jn = %d\n", dom[i].Gcc.js, dom[i].Gcc.je, dom[i].Gcc.jn);
    printf("    jsb = %d, jeb = %d, jnb = %d\n", dom[i].Gcc.jsb, dom[i].Gcc.jeb,
      dom[i].Gcc.jnb);
    printf("    ks = %d, ke = %d, kn = %d\n", dom[i].Gcc.ks, dom[i].Gcc.ke, dom[i].Gcc.kn);
    printf("    ksb = %d, keb = %d, knb = %d\n", dom[i].Gcc.ksb, dom[i].Gcc.keb,
      dom[i].Gcc.knb);
    printf("    s1 = %d, s2 = %d, s3 = %d\n", dom[i].Gcc.s1, dom[i].Gcc.s2,
      dom[i].Gcc.s3);
    printf("    s1b = %d, s2b = %d, s3b = %d\n", dom[i].Gcc.s1b, dom[i].Gcc.s2b,
      dom[i].Gcc.s3b);
    printf("    _is = %d, _ie = %d, _in = %d\n", dom[i].Gcc._is, dom[i].Gcc._ie,
      dom[i].Gcc._in);
    printf("    _isb = %d, _ieb = %d, _inb = %d\n", dom[i].Gcc._isb, dom[i].Gcc._ieb,
      dom[i].Gcc._inb);
    printf("    _js = %d, _je = %d, _jn = %d\n", dom[i].Gcc._js, dom[i].Gcc._je,
      dom[i].Gcc._jn);
    printf("    _jsb = %d, _jeb = %d, _jnb = %d\n", dom[i].Gcc._jsb, dom[i].Gcc._jeb,
      dom[i].Gcc._jnb);
    printf("    _ks = %d, _ke = %d, _kn = %d\n", dom[i].Gcc._ks, dom[i].Gcc._ke,
      dom[i].Gcc._kn);
    printf("    _ksb = %d, _keb = %d, _knb = %d\n", dom[i].Gcc._ksb, dom[i].Gcc._keb,
      dom[i].Gcc._knb);
    printf("    _s1 = %d, _s2 = %d, _s3 = %d\n", dom[i].Gcc._s1, dom[i].Gcc._s2,
      dom[i].Gcc._s3);
    printf("    _s1b = %d, _s2b = %d, _s3b = %d\n", dom[i].Gcc._s1b, dom[i].Gcc._s2b,
      dom[i].Gcc._s3b);
    printf("  dom[%d].Gfx:\n", i);
    printf("    is = %d, ie = %d, in = %d\n", dom[i].Gfx.is, dom[i].Gfx.ie, dom[i].Gfx.in);
    printf("    isb = %d, ieb = %d, inb = %d\n", dom[i].Gfx.isb, dom[i].Gfx.ieb,
      dom[i].Gfx.inb);
    printf("    js = %d, je = %d, jn = %d\n", dom[i].Gfx.js, dom[i].Gfx.je, dom[i].Gfx.jn);
    printf("    jsb = %d, jeb = %d, jnb = %d\n", dom[i].Gfx.jsb, dom[i].Gfx.jeb,
      dom[i].Gfx.jnb);
    printf("    ks = %d, ke = %d, kn = %d\n", dom[i].Gfx.ks, dom[i].Gfx.ke, dom[i].Gfx.kn);
    printf("    ksb = %d, keb = %d, knb = %d\n", dom[i].Gfx.ksb, dom[i].Gfx.keb,
      dom[i].Gfx.knb);
    printf("    s1 = %d, s2 = %d, s3 = %d\n", dom[i].Gfx.s1, dom[i].Gfx.s2,
      dom[i].Gfx.s3);
    printf("    s1b = %d, s2b = %d, s3b = %d\n", dom[i].Gfx.s1b, dom[i].Gfx.s2b,
      dom[i].Gfx.s3b);
    printf("    _is = %d, _ie = %d, _in = %d\n", dom[i].Gfx._is, dom[i].Gfx._ie,
      dom[i].Gfx._in);
    printf("    _isb = %d, _ieb = %d, _inb = %d\n", dom[i].Gfx._isb, dom[i].Gfx._ieb,
      dom[i].Gfx._inb);
    printf("    _js = %d, _je = %d, _jn = %d\n", dom[i].Gfx._js, dom[i].Gfx._je,
      dom[i].Gfx._jn);
    printf("    _jsb = %d, _jeb = %d, _jnb = %d\n", dom[i].Gfx._jsb, dom[i].Gfx._jeb,
      dom[i].Gfx._jnb);
    printf("    _ks = %d, _ke = %d, _kn = %d\n", dom[i].Gfx._ks, dom[i].Gfx._ke,
      dom[i].Gfx._kn);
    printf("    _ksb = %d, _keb = %d, _knb = %d\n", dom[i].Gfx._ksb, dom[i].Gfx._keb,
      dom[i].Gfx._knb);
    printf("    _s1 = %d, _s2 = %d, _s3 = %d\n", dom[i].Gfx._s1, dom[i].Gfx._s2,
      dom[i].Gfx._s3);
    printf("    _s1b = %d, _s2b = %d, _s3b = %d\n", dom[i].Gfx._s1b, dom[i].Gfx._s2b,
      dom[i].Gfx._s3b);
    printf("  dom[%d].Gfy:\n", i);
    printf("    is = %d, ie = %d, in = %d\n", dom[i].Gfy.is, dom[i].Gfy.ie, dom[i].Gfy.in);
    printf("    isb = %d, ieb = %d, inb = %d\n", dom[i].Gfy.isb, dom[i].Gfy.ieb,
      dom[i].Gfy.inb);
    printf("    js = %d, je = %d, jn = %d\n", dom[i].Gfy.js, dom[i].Gfy.je, dom[i].Gfy.jn);
    printf("    jsb = %d, jeb = %d, jnb = %d\n", dom[i].Gfy.jsb, dom[i].Gfy.jeb,
      dom[i].Gfy.jnb);
    printf("    ks = %d, ke = %d, kn = %d\n", dom[i].Gfy.ks, dom[i].Gfy.ke, dom[i].Gfy.kn);
    printf("    ksb = %d, keb = %d, knb = %d\n", dom[i].Gfy.ksb, dom[i].Gfy.keb,
      dom[i].Gfy.knb);
    printf("    s1 = %d, s2 = %d, s3 = %d\n", dom[i].Gfy.s1, dom[i].Gfy.s2,
      dom[i].Gfy.s3);
    printf("    s1b = %d, s2b = %d, s3b = %d\n", dom[i].Gfy.s1b, dom[i].Gfy.s2b,
      dom[i].Gfy.s3b);
    printf("    _is = %d, _ie = %d, _in = %d\n", dom[i].Gfy._is, dom[i].Gfy._ie,
      dom[i].Gfy._in);
    printf("    _isb = %d, _ieb = %d, _inb = %d\n", dom[i].Gfy._isb, dom[i].Gfy._ieb,
      dom[i].Gfy._inb);
    printf("    _js = %d, _je = %d, _jn = %d\n", dom[i].Gfy._js, dom[i].Gfy._je,
      dom[i].Gfy._jn);
    printf("    _jsb = %d, _jeb = %d, _jnb = %d\n", dom[i].Gfy._jsb, dom[i].Gfy._jeb,
      dom[i].Gfy._jnb);
    printf("    _ks = %d, _ke = %d, _kn = %d\n", dom[i].Gfy._ks, dom[i].Gfy._ke,
      dom[i].Gfy._kn);
    printf("    _ksb = %d, _keb = %d, _knb = %d\n", dom[i].Gfy._ksb, dom[i].Gfy._keb,
      dom[i].Gfy._knb);
    printf("    _s1 = %d, _s2 = %d, _s3 = %d\n", dom[i].Gfy._s1, dom[i].Gfy._s2,
      dom[i].Gfy._s3);
    printf("    _s1b = %d, _s2b = %d, _s3b = %d\n", dom[i].Gfy._s1b, dom[i].Gfy._s2b,
      dom[i].Gfy._s3b);
    printf("  dom[%d].Gfz:\n", i);
    printf("    is = %d, ie = %d, in = %d\n", dom[i].Gfz.is, dom[i].Gfz.ie, dom[i].Gfz.in);
    printf("    isb = %d, ieb = %d, inb = %d\n", dom[i].Gfz.isb, dom[i].Gfz.ieb,
      dom[i].Gfz.inb);
    printf("    js = %d, je = %d, jn = %d\n", dom[i].Gfz.js, dom[i].Gfz.je, dom[i].Gfz.jn);
    printf("    jsb = %d, jeb = %d, jnb = %d\n", dom[i].Gfz.jsb, dom[i].Gfz.jeb,
      dom[i].Gfz.jnb);
    printf("    ks = %d, ke = %d, kn = %d\n", dom[i].Gfz.ks, dom[i].Gfz.ke, dom[i].Gfz.kn);
    printf("    ksb = %d, keb = %d, knb = %d\n", dom[i].Gfz.ksb, dom[i].Gfz.keb,
      dom[i].Gfz.knb);
    printf("    s1 = %d, s2 = %d, s3 = %d\n", dom[i].Gfz.s1, dom[i].Gfz.s2,
      dom[i].Gfz.s3);
    printf("    s1b = %d, s2b = %d, s3b = %d\n", dom[i].Gfz.s1b, dom[i].Gfz.s2b,
      dom[i].Gfz.s3b);
    printf("    _is = %d, _ie = %d, _in = %d\n", dom[i].Gfz._is, dom[i].Gfz._ie,
      dom[i].Gfz._in);
    printf("    _isb = %d, _ieb = %d, _inb = %d\n", dom[i].Gfz._isb, dom[i].Gfz._ieb,
      dom[i].Gfz._inb);
    printf("    _js = %d, _je = %d, _jn = %d\n", dom[i].Gfz._js, dom[i].Gfz._je,
      dom[i].Gfz._jn);
    printf("    _jsb = %d, _jeb = %d, _jnb = %d\n", dom[i].Gfz._jsb, dom[i].Gfz._jeb,
      dom[i].Gfz._jnb);
    printf("    _ks = %d, _ke = %d, _kn = %d\n", dom[i].Gfz._ks, dom[i].Gfz._ke,
      dom[i].Gfz._kn);
    printf("    _ksb = %d, _keb = %d, _knb = %d\n", dom[i].Gfz._ksb, dom[i].Gfz._keb,
      dom[i].Gfz._knb);
    printf("    _s1 = %d, _s2 = %d, _s3 = %d\n", dom[i].Gfz._s1, dom[i].Gfz._s2,
      dom[i].Gfz._s3);
    printf("    _s1b = %d, _s2b = %d, _s3b = %d\n", dom[i].Gfz._s1b, dom[i].Gfz._s2b,
      dom[i].Gfz._s3b);
  }
  printf("Physical Parameters:\n");
  printf("  rho_f = %e\n", rho_f);
  printf("  nu = %e\n", nu);
  printf("  mu = %e\n", mu);
  printf("Simulation Parameters:\n");
  printf("  duration = %e\n", duration);
  printf("  CFL = %e\n", CFL);
  printf("  pp_max_iter = %d\n", pp_max_iter);
  printf("  pp_residual = %e\n", pp_residual);
  printf("Boundary Conditions: (0 = PERIODIC, 1 = DIRICHLET, 2 = NEUMANN)\n");
  printf("  bc.pW = %d", bc.pW);
  if(bc.pW == DIRICHLET) printf(" %f", bc.pWD);
  printf(", bc.pE = %d", bc.pE);
  if(bc.pE == DIRICHLET) printf(" %f", bc.pED);
  printf("\n");
  printf("  bc.pS = %d", bc.pS);
  if(bc.pS == DIRICHLET) printf(" %f", bc.pSD);
  printf(", bc.pN = %d", bc.pN);
  if(bc.pN == DIRICHLET) printf(" %f", bc.pND);
  printf("\n");
  printf("  bc.pB = %d", bc.pB);
  if(bc.pB == DIRICHLET) printf(" %f", bc.pBD);
  printf(", bc.pT = %d", bc.pT);
  if(bc.pT == DIRICHLET) printf(" %f", bc.pTD);
  printf("\n");
  printf("  bc.uW = %d", bc.uW);
  if(bc.uW == DIRICHLET) printf(" %f", bc.uWDm);
  printf(", bc.uE = %d", bc.uE);
  if(bc.uE == DIRICHLET) printf(" %f", bc.uEDm);
  printf("\n");
  printf("  bc.uS = %d", bc.uS);
  if(bc.uS == DIRICHLET) printf(" %f", bc.uSDm);
  printf(", bc.uN = %d", bc.uN);
  if(bc.uN == DIRICHLET) printf(" %f", bc.uNDm);
  printf("\n");
  printf("  bc.uB = %d", bc.uB);
  if(bc.uB == DIRICHLET) printf(" %f", bc.uBDm);
  printf(", bc.uT = %d", bc.uT);
  if(bc.uT == DIRICHLET) printf(" %f", bc.uTDm);
  printf("\n");
  printf("  bc.vW = %d", bc.vW);
  if(bc.vW == DIRICHLET) printf(" %f", bc.vWDm);
  printf(", bc.vE = %d", bc.vE);
  if(bc.vE == DIRICHLET) printf(" %f", bc.vEDm);
  printf("\n");
  printf("  bc.vS = %d", bc.vS);
  if(bc.vS == DIRICHLET) printf(" %f", bc.vSDm);
  printf(", bc.vN = %d", bc.vN);
  if(bc.vN == DIRICHLET) printf(" %f", bc.vNDm);
  printf("\n");
  printf("  bc.vB = %d", bc.vB);
  if(bc.vB == DIRICHLET) printf(" %f", bc.vBDm);
  printf(", bc.vT = %d", bc.vT);
  if(bc.vT == DIRICHLET) printf(" %f", bc.vTDm);
  printf("\n");
  printf("  bc.wW = %d", bc.wW);
  if(bc.wW == DIRICHLET) printf(" %f", bc.wWDm);
  printf(", bc.wE = %d", bc.wE);
  if(bc.wE == DIRICHLET) printf(" %f", bc.wEDm);
  printf("\n");
  printf("  bc.wS = %d", bc.wS);
  if(bc.wS == DIRICHLET) printf(" %f", bc.wSDm);
  printf(", bc.wN = %d", bc.wN);
  if(bc.wN == DIRICHLET) printf(" %f", bc.wNDm);
  printf("\n");
  printf("  bc.wB = %d", bc.wB);
  if(bc.wB == DIRICHLET) printf(" %f", bc.wBDm);
  printf(", bc.wT = %d", bc.wT);
  if(bc.wT == DIRICHLET) printf(" %f", bc.wTDm);
  printf("\n");

  printf("Screen offsets:\n");
  printf("  bc.dsW = %f\n", bc.dsW);
  printf("  bc.dsE = %f\n", bc.dsE);
  printf("  bc.dsS = %f\n", bc.dsS);
  printf("  bc.dsN = %f\n", bc.dsN);
  printf("  bc.dsB = %f\n", bc.dsB);
  printf("  bc.dsT = %f\n", bc.dsT);
  printf("\n");

  printf("Applied Pressure Gradient:\n");
  printf("  gradP.x = %f\n", gradP.x);
  printf("  gradP.xm = %f\n", gradP.xm);
  printf("  gradP.xa = %f\n", gradP.xa);
  printf("  gradP.y = %f\n", gradP.y);
  printf("  gradP.ym = %f\n", gradP.ym);
  printf("  gradP.ya = %f\n", gradP.ya);
  printf("  gradP.z = %f\n", gradP.z);
  printf("  gradP.zm = %f\n", gradP.zm);
  printf("  gradP.za = %f\n", gradP.za);
  printf("Applied Body Forces:\n");
  printf("  g.x = %f\n", g.x);
  printf("  g.xm = %f\n", g.xm);
  printf("  g.xa = %f\n", g.xa);
  printf("  g.y = %f\n", g.y);
  printf("  g.ym = %f\n", g.ym);
  printf("  g.ya = %f\n", g.ya);
  printf("  g.z = %f\n", g.z);
  printf("  g.zm = %f\n", g.zm);
  printf("  g.za = %f\n", g.za);
}

int domain_init(void)
{
  int i,j,k;  // iterators

  // make sure there are enough GPU devices in the given range
  if(nsubdom > dev_end - dev_start + 1) {
    return EXIT_FAILURE;
  }

  // calculate domain sizes
  Dom.xl = Dom.xe - Dom.xs;
  Dom.yl = Dom.ye - Dom.ys;
  Dom.zl = Dom.ze - Dom.zs;

  // calculate cell sizes
  Dom.dx = Dom.xl / Dom.xn;
  Dom.dy = Dom.yl / Dom.yn;
  Dom.dz = Dom.zl / Dom.zn;

  // set up grids
  // Gcc
  Dom.Gcc.is = DOM_BUF;
  Dom.Gcc.isb = Dom.Gcc.is - DOM_BUF;
  Dom.Gcc.in = Dom.xn;
  Dom.Gcc.inb = Dom.Gcc.in + 2 * DOM_BUF;
  Dom.Gcc.ie = Dom.Gcc.is + Dom.Gcc.in;
  Dom.Gcc.ieb = Dom.Gcc.ie + DOM_BUF;

  Dom.Gcc.js = DOM_BUF;
  Dom.Gcc.jsb = Dom.Gcc.js - DOM_BUF;
  Dom.Gcc.jn = Dom.yn;
  Dom.Gcc.jnb = Dom.Gcc.jn + 2 * DOM_BUF;
  Dom.Gcc.je = Dom.Gcc.js + Dom.Gcc.jn;
  Dom.Gcc.jeb = Dom.Gcc.je + DOM_BUF;

  Dom.Gcc.ks = DOM_BUF;
  Dom.Gcc.ksb = Dom.Gcc.ks - DOM_BUF;
  Dom.Gcc.kn = Dom.zn;
  Dom.Gcc.knb = Dom.Gcc.kn + 2 * DOM_BUF;
  Dom.Gcc.ke = DOM_BUF + Dom.Gcc.kn;
  Dom.Gcc.keb = Dom.Gcc.ke + DOM_BUF;

  Dom.Gcc.s1 = Dom.Gcc.in;
  Dom.Gcc.s2 = Dom.Gcc.s1 * Dom.Gcc.jn;
  Dom.Gcc.s3 = Dom.Gcc.s2 * Dom.Gcc.kn;
  Dom.Gcc.s1b = Dom.Gcc.inb;
  Dom.Gcc.s2b = Dom.Gcc.s1b * Dom.Gcc.jnb;
  Dom.Gcc.s3b = Dom.Gcc.s2b * Dom.Gcc.knb;

  // Gfx
  Dom.Gfx.is = DOM_BUF;
  Dom.Gfx.isb = Dom.Gfx.is - DOM_BUF;
  Dom.Gfx.in = Dom.xn + 1;
  Dom.Gfx.inb = Dom.Gfx.in + 2 * DOM_BUF;
  Dom.Gfx.ie = Dom.Gfx.is + Dom.Gfx.in;
  Dom.Gfx.ieb = Dom.Gfx.ie + DOM_BUF;

  Dom.Gfx.js = DOM_BUF;
  Dom.Gfx.jsb = Dom.Gfx.js - DOM_BUF;
  Dom.Gfx.jn = Dom.yn;
  Dom.Gfx.jnb = Dom.Gfx.jn + 2 * DOM_BUF;
  Dom.Gfx.je = Dom.Gfx.js + Dom.Gfx.jn;
  Dom.Gfx.jeb = Dom.Gfx.je + DOM_BUF;

  Dom.Gfx.ks = DOM_BUF;
  Dom.Gfx.ksb = Dom.Gfx.ks - DOM_BUF;
  Dom.Gfx.kn = Dom.zn;
  Dom.Gfx.knb = Dom.Gfx.kn + 2 * DOM_BUF;
  Dom.Gfx.ke = Dom.Gfx.ks + Dom.Gfx.kn;
  Dom.Gfx.keb = Dom.Gfx.ke + DOM_BUF;

  Dom.Gfx.s1 = Dom.Gfx.in;
  Dom.Gfx.s2 = Dom.Gfx.s1 * Dom.Gfx.jn;
  Dom.Gfx.s3 = Dom.Gfx.s2 * Dom.Gfx.kn;
  Dom.Gfx.s1b = Dom.Gfx.inb;
  Dom.Gfx.s2b = Dom.Gfx.s1b * Dom.Gfx.jnb;
  Dom.Gfx.s3b = Dom.Gfx.s2b * Dom.Gfx.knb;

  // Gfy
  Dom.Gfy.is = DOM_BUF;
  Dom.Gfy.isb = Dom.Gfy.is - DOM_BUF;
  Dom.Gfy.in = Dom.xn;
  Dom.Gfy.inb = Dom.Gfy.in + 2 * DOM_BUF;
  Dom.Gfy.ie = Dom.Gfy.is + Dom.Gfy.in;
  Dom.Gfy.ieb = Dom.Gfy.ie + DOM_BUF;

  Dom.Gfy.js = DOM_BUF;
  Dom.Gfy.jsb = Dom.Gfy.js - DOM_BUF;
  Dom.Gfy.jn = Dom.yn + 1;
  Dom.Gfy.jnb = Dom.Gfy.jn + 2 * DOM_BUF;
  Dom.Gfy.je = Dom.Gfy.js + Dom.Gfy.jn;
  Dom.Gfy.jeb = Dom.Gfy.je + DOM_BUF;

  Dom.Gfy.ks = DOM_BUF;
  Dom.Gfy.ksb = Dom.Gfy.ks - DOM_BUF;
  Dom.Gfy.kn = Dom.zn;
  Dom.Gfy.knb = Dom.Gfy.kn + 2 * DOM_BUF;
  Dom.Gfy.ke = Dom.Gfy.ks + Dom.Gfy.kn;
  Dom.Gfy.keb = Dom.Gfy.ke + DOM_BUF;

  Dom.Gfy.s1 = Dom.Gfy.in;
  Dom.Gfy.s2 = Dom.Gfy.s1 * Dom.Gfy.jn;
  Dom.Gfy.s3 = Dom.Gfy.s2 * Dom.Gfy.kn;
  Dom.Gfy.s1b = Dom.Gfy.inb;
  Dom.Gfy.s2b = Dom.Gfy.s1b * Dom.Gfy.jnb;
  Dom.Gfy.s3b = Dom.Gfy.s2b * Dom.Gfy.knb;

  // Gfz
  Dom.Gfz.is = DOM_BUF;
  Dom.Gfz.isb = Dom.Gfz.is - DOM_BUF;
  Dom.Gfz.in = Dom.xn;
  Dom.Gfz.inb = Dom.Gfz.in + 2 * DOM_BUF;
  Dom.Gfz.ie = Dom.Gfz.is + Dom.Gfz.in;
  Dom.Gfz.ieb = Dom.Gfz.ie + DOM_BUF;

  Dom.Gfz.js = DOM_BUF;
  Dom.Gfz.jsb = Dom.Gfz.js - DOM_BUF;
  Dom.Gfz.jn = Dom.yn;
  Dom.Gfz.jnb = Dom.Gfz.jn + 2 * DOM_BUF;
  Dom.Gfz.je = Dom.Gfz.js + Dom.Gfz.jn;
  Dom.Gfz.jeb = Dom.Gfz.je + DOM_BUF;

  Dom.Gfz.ks = DOM_BUF;
  Dom.Gfz.ksb = Dom.Gfz.ks - DOM_BUF;
  Dom.Gfz.kn = Dom.zn + 1;
  Dom.Gfz.knb = Dom.Gfz.kn + 2 * DOM_BUF;
  Dom.Gfz.ke = Dom.Gfz.ks + Dom.Gfz.kn;
  Dom.Gfz.keb = Dom.Gfz.ke + DOM_BUF;

  Dom.Gfz.s1 = Dom.Gfz.in;
  Dom.Gfz.s2 = Dom.Gfz.s1 * Dom.Gfz.jn;
  Dom.Gfz.s3 = Dom.Gfz.s2 * Dom.Gfz.kn;
  Dom.Gfz.s1b = Dom.Gfz.inb;
  Dom.Gfz.s2b = Dom.Gfz.s1b * Dom.Gfz.jnb;
  Dom.Gfz.s3b = Dom.Gfz.s2b * Dom.Gfz.knb;

  // initialize subdomains
  for(i = 0; i < nsubdom; i++) {
    dom[i].xl = dom[i].xe - dom[i].xs;
    dom[i].yl = dom[i].ye - dom[i].ys;
    dom[i].zl = dom[i].ze - dom[i].zs;
    dom[i].dx = dom[i].xl / dom[i].xn;
    dom[i].dy = dom[i].yl / dom[i].yn;
    dom[i].dz = dom[i].zl / dom[i].zn;

    // TODO: this algorithm will fail if subdomains are not numbered in
    // increasing order.  This will need to be fixed before going to a machine
    // with more than two GPU's.
    // Gcc
    if(dom[i].W > -1)
      dom[i].Gcc.is = dom[dom[i].W].Gcc.ie;
    else
      dom[i].Gcc.is = DOM_BUF;
    dom[i].Gcc.isb = dom[i].Gcc.is - DOM_BUF;
    dom[i].Gcc.in = dom[i].xn;
    dom[i].Gcc.inb = dom[i].Gcc.in + 2 * DOM_BUF;
    dom[i].Gcc.ie = dom[i].Gcc.is + dom[i].Gcc.in;
    dom[i].Gcc.ieb = dom[i].Gcc.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gcc.js = dom[dom[i].S].Gcc.je;
    else
      dom[i].Gcc.js = DOM_BUF;
    dom[i].Gcc.jsb = dom[i].Gcc.js - DOM_BUF;
    dom[i].Gcc.jn = dom[i].yn;
    dom[i].Gcc.jnb = dom[i].Gcc.jn + 2 * DOM_BUF;
    dom[i].Gcc.je = dom[i].Gcc.js + dom[i].Gcc.jn;
    dom[i].Gcc.jeb = dom[i].Gcc.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gcc.ks = dom[dom[i].B].Gcc.ke;
    else
      dom[i].Gcc.ks = DOM_BUF;
    dom[i].Gcc.ksb = dom[i].Gcc.ks - DOM_BUF;
    dom[i].Gcc.kn = dom[i].zn;
    dom[i].Gcc.knb = dom[i].Gcc.kn + 2 * DOM_BUF;
    dom[i].Gcc.ke = DOM_BUF + dom[i].Gcc.kn;
    dom[i].Gcc.keb = dom[i].Gcc.ke + DOM_BUF;

    dom[i].Gcc.s1 = dom[i].Gcc.in;
    dom[i].Gcc.s2 = dom[i].Gcc.s1 * dom[i].Gcc.jn;
    dom[i].Gcc.s3 = dom[i].Gcc.s2 * dom[i].Gcc.kn;
    dom[i].Gcc.s1b = dom[i].Gcc.inb;
    dom[i].Gcc.s2b = dom[i].Gcc.s1b * dom[i].Gcc.jnb;
    dom[i].Gcc.s3b = dom[i].Gcc.s2b * dom[i].Gcc.knb;

    dom[i].Gcc._is = DOM_BUF;
    dom[i].Gcc._isb = dom[i].Gcc._is - DOM_BUF;
    dom[i].Gcc._in = dom[i].xn;
    dom[i].Gcc._inb = dom[i].Gcc._in + 2 * DOM_BUF;
    dom[i].Gcc._ie = dom[i].Gcc._is + dom[i].Gcc._in;
    dom[i].Gcc._ieb = dom[i].Gcc._ie + DOM_BUF;

    dom[i].Gcc._js = DOM_BUF;
    dom[i].Gcc._jsb = dom[i].Gcc._js - DOM_BUF;
    dom[i].Gcc._jn = dom[i].yn;
    dom[i].Gcc._jnb = dom[i].Gcc._jn + 2 * DOM_BUF;
    dom[i].Gcc._je = dom[i].Gcc._js + dom[i].Gcc._jn;
    dom[i].Gcc._jeb = dom[i].Gcc._je + DOM_BUF;

    dom[i].Gcc._ks = DOM_BUF;
    dom[i].Gcc._ksb = dom[i].Gcc._ks - DOM_BUF;
    dom[i].Gcc._kn = dom[i].zn;
    dom[i].Gcc._knb = dom[i].Gcc._kn + 2 * DOM_BUF;
    dom[i].Gcc._ke = dom[i].Gcc._ks + dom[i].Gcc._kn;
    dom[i].Gcc._keb = dom[i].Gcc._ke + DOM_BUF;

    dom[i].Gcc._s1 = dom[i].Gcc._in;
    dom[i].Gcc._s2 = dom[i].Gcc._s1 * dom[i].Gcc._jn;
    dom[i].Gcc._s3 = dom[i].Gcc._s2 * dom[i].Gcc._kn;
    dom[i].Gcc._s1b = dom[i].Gcc._inb;
    dom[i].Gcc._s2b = dom[i].Gcc._s1b * dom[i].Gcc._jnb;
    dom[i].Gcc._s3b = dom[i].Gcc._s2b * dom[i].Gcc._knb;

    // Gfx
    if(dom[i].W > -1)
      dom[i].Gfx.is = dom[dom[i].W].Gfx.ie - 1;
    else
      dom[i].Gfx.is = DOM_BUF;
    dom[i].Gfx.isb = dom[i].Gfx.is - DOM_BUF;
    dom[i].Gfx.in = dom[i].xn + 1;
    dom[i].Gfx.inb = dom[i].Gfx.in + 2 * DOM_BUF;
    dom[i].Gfx.ie = dom[i].Gfx.is + dom[i].Gfx.in;
    dom[i].Gfx.ieb = dom[i].Gfx.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfx.js = dom[dom[i].S].Gfx.je - 1;
    else
      dom[i].Gfx.js = DOM_BUF;
    dom[i].Gfx.jsb = dom[i].Gfx.js - DOM_BUF;
    dom[i].Gfx.jn = dom[i].yn;
    dom[i].Gfx.jnb = dom[i].Gfx.jn + 2 * DOM_BUF;
    dom[i].Gfx.je = dom[i].Gfx.js + dom[i].Gfx.jn;
    dom[i].Gfx.jeb = dom[i].Gfx.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfx.ks = dom[dom[i].B].Gfx.ke - 1;
    else
      dom[i].Gfx.ks = DOM_BUF;
    dom[i].Gfx.ksb = dom[i].Gfx.ks - DOM_BUF;
    dom[i].Gfx.kn = dom[i].zn;
    dom[i].Gfx.knb = dom[i].Gfx.kn + 2 * DOM_BUF;
    dom[i].Gfx.ke = dom[i].Gfx.ks + dom[i].Gfx.kn;
    dom[i].Gfx.keb = dom[i].Gfx.ke + DOM_BUF;

    dom[i].Gfx.s1 = dom[i].Gfx.in;
    dom[i].Gfx.s2 = dom[i].Gfx.s1 * dom[i].Gfx.jn;
    dom[i].Gfx.s3 = dom[i].Gfx.s2 * dom[i].Gfx.kn;
    dom[i].Gfx.s1b = dom[i].Gfx.inb;
    dom[i].Gfx.s2b = dom[i].Gfx.s1b * dom[i].Gfx.jnb;
    dom[i].Gfx.s3b = dom[i].Gfx.s2b * dom[i].Gfx.knb;

    dom[i].Gfx._is = DOM_BUF;
    dom[i].Gfx._isb = dom[i].Gfx._is - DOM_BUF;
    dom[i].Gfx._in = dom[i].xn + 1;
    dom[i].Gfx._inb = dom[i].Gfx._in + 2 * DOM_BUF;
    dom[i].Gfx._ie = dom[i].Gfx._is + dom[i].Gfx._in;
    dom[i].Gfx._ieb = dom[i].Gfx._ie + DOM_BUF;

    dom[i].Gfx._js = DOM_BUF;
    dom[i].Gfx._jsb = dom[i].Gfx._js - DOM_BUF;
    dom[i].Gfx._jn = dom[i].yn;
    dom[i].Gfx._jnb = dom[i].Gfx._jn + 2 * DOM_BUF;
    dom[i].Gfx._je = dom[i].Gfx._js + dom[i].Gfx._jn;
    dom[i].Gfx._jeb = dom[i].Gfx._je + DOM_BUF;

    dom[i].Gfx._ks = DOM_BUF;
    dom[i].Gfx._ksb = dom[i].Gfx._ks - DOM_BUF;
    dom[i].Gfx._kn = dom[i].zn;
    dom[i].Gfx._knb = dom[i].Gfx._kn + 2 * DOM_BUF;
    dom[i].Gfx._ke = dom[i].Gfx._ks + dom[i].Gfx._kn;
    dom[i].Gfx._keb = dom[i].Gfx._ke + DOM_BUF;

    dom[i].Gfx._s1 = dom[i].Gfx._in;
    dom[i].Gfx._s2 = dom[i].Gfx._s1 * dom[i].Gfx._jn;
    dom[i].Gfx._s3 = dom[i].Gfx._s2 * dom[i].Gfx._kn;
    dom[i].Gfx._s1b = dom[i].Gfx._inb;
    dom[i].Gfx._s2b = dom[i].Gfx._s1b * dom[i].Gfx._jnb;
    dom[i].Gfx._s3b = dom[i].Gfx._s2b * dom[i].Gfx._knb;

    // Gfy
    if(dom[i].W > -1)
      dom[i].Gfy.is = dom[dom[i].W].Gfy.ie;
    else
      dom[i].Gfy.is = DOM_BUF;
    dom[i].Gfy.isb = dom[i].Gfy.is - DOM_BUF;
    dom[i].Gfy.in = dom[i].xn;
    dom[i].Gfy.inb = dom[i].Gfy.in + 2 * DOM_BUF;
    dom[i].Gfy.ie = dom[i].Gfy.is + dom[i].Gfy.in;
    dom[i].Gfy.ieb = dom[i].Gfy.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfy.js = dom[dom[i].S].Gfy.je;
    else
      dom[i].Gfy.js = DOM_BUF;
    dom[i].Gfy.jsb = dom[i].Gfy.js - DOM_BUF;
    dom[i].Gfy.jn = dom[i].yn + 1;
    dom[i].Gfy.jnb = dom[i].Gfy.jn + 2 * DOM_BUF;
    dom[i].Gfy.je = dom[i].Gfy.js + dom[i].Gfy.jn;
    dom[i].Gfy.jeb = dom[i].Gfy.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfy.ks = dom[dom[i].B].Gfy.ke;
    else
      dom[i].Gfy.ks = DOM_BUF;
    dom[i].Gfy.ksb = dom[i].Gfy.ks - DOM_BUF;
    dom[i].Gfy.kn = dom[i].zn;
    dom[i].Gfy.knb = dom[i].Gfy.kn + 2 * DOM_BUF;
    dom[i].Gfy.ke = dom[i].Gfy.ks + dom[i].Gfy.kn;
    dom[i].Gfy.keb = dom[i].Gfy.ke + DOM_BUF;

    dom[i].Gfy.s1 = dom[i].Gfy.in;
    dom[i].Gfy.s2 = dom[i].Gfy.s1 * dom[i].Gfy.jn;
    dom[i].Gfy.s3 = dom[i].Gfy.s2 * dom[i].Gfy.kn;
    dom[i].Gfy.s1b = dom[i].Gfy.inb;
    dom[i].Gfy.s2b = dom[i].Gfy.s1b * dom[i].Gfy.jnb;
    dom[i].Gfy.s3b = dom[i].Gfy.s2b * dom[i].Gfy.knb;

    dom[i].Gfy._is = DOM_BUF;
    dom[i].Gfy._isb = dom[i].Gfy._is - DOM_BUF;
    dom[i].Gfy._in = dom[i].xn;
    dom[i].Gfy._inb = dom[i].Gfy._in + 2 * DOM_BUF;
    dom[i].Gfy._ie = dom[i].Gfy._is + dom[i].Gfy._in;
    dom[i].Gfy._ieb = dom[i].Gfy._ie + DOM_BUF;

    dom[i].Gfy._js = DOM_BUF;
    dom[i].Gfy._jsb = dom[i].Gfy._js - DOM_BUF;
    dom[i].Gfy._jn = dom[i].yn + 1;
    dom[i].Gfy._jnb = dom[i].Gfy._jn + 2 * DOM_BUF;
    dom[i].Gfy._je = dom[i].Gfy._js + dom[i].Gfy._jn;
    dom[i].Gfy._jeb = dom[i].Gfy._je + DOM_BUF;

    dom[i].Gfy._ks = DOM_BUF;
    dom[i].Gfy._ksb = dom[i].Gfy._ks - DOM_BUF;
    dom[i].Gfy._kn = dom[i].zn;
    dom[i].Gfy._knb = dom[i].Gfy._kn + 2 * DOM_BUF;
    dom[i].Gfy._ke = dom[i].Gfy._ks + dom[i].Gfy._kn;
    dom[i].Gfy._keb = dom[i].Gfy._ke + DOM_BUF;

    dom[i].Gfy._s1 = dom[i].Gfy._in;
    dom[i].Gfy._s2 = dom[i].Gfy._s1 * dom[i].Gfy._jn;
    dom[i].Gfy._s3 = dom[i].Gfy._s2 * dom[i].Gfy._kn;
    dom[i].Gfy._s1b = dom[i].Gfy._inb;
    dom[i].Gfy._s2b = dom[i].Gfy._s1b * dom[i].Gfy._jnb;
    dom[i].Gfy._s3b = dom[i].Gfy._s2b * dom[i].Gfy._knb;

    // Gfz
    if(dom[i].W > -1)
      dom[i].Gfz.is = dom[dom[i].W].Gfz.ie;
    else
      dom[i].Gfz.is = DOM_BUF;
    dom[i].Gfz.isb = dom[i].Gfz.is - DOM_BUF;
    dom[i].Gfz.in = dom[i].xn;
    dom[i].Gfz.inb = dom[i].Gfz.in + 2 * DOM_BUF;
    dom[i].Gfz.ie = dom[i].Gfz.is + dom[i].Gfz.in;
    dom[i].Gfz.ieb = dom[i].Gfz.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfz.js = dom[dom[i].S].Gfz.je;
    else
      dom[i].Gfz.js = DOM_BUF;
    dom[i].Gfz.jsb = dom[i].Gfz.js - DOM_BUF;
    dom[i].Gfz.jn = dom[i].yn;
    dom[i].Gfz.jnb = dom[i].Gfz.jn + 2 * DOM_BUF;
    dom[i].Gfz.je = dom[i].Gfz.js + dom[i].Gfz.jn;
    dom[i].Gfz.jeb = dom[i].Gfz.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfz.ks = dom[dom[i].B].Gfz.ke;
    else
      dom[i].Gfz.ks = DOM_BUF;
    dom[i].Gfz.ksb = dom[i].Gfz.ks - DOM_BUF;
    dom[i].Gfz.kn = dom[i].zn + 1;
    dom[i].Gfz.knb = dom[i].Gfz.kn + 2 * DOM_BUF;
    dom[i].Gfz.ke = dom[i].Gfz.ks + dom[i].Gfz.kn;
    dom[i].Gfz.keb = dom[i].Gfz.ke + DOM_BUF;

    dom[i].Gfz.s1 = dom[i].Gfz.in;
    dom[i].Gfz.s2 = dom[i].Gfz.s1 * dom[i].Gfz.jn;
    dom[i].Gfz.s3 = dom[i].Gfz.s2 * dom[i].Gfz.kn;
    dom[i].Gfz.s1b = dom[i].Gfz.inb;
    dom[i].Gfz.s2b = dom[i].Gfz.s1b * dom[i].Gfz.jnb;
    dom[i].Gfz.s3b = dom[i].Gfz.s2b * dom[i].Gfz.knb;

    dom[i].Gfz._is = DOM_BUF;
    dom[i].Gfz._isb = dom[i].Gfz._is - DOM_BUF;
    dom[i].Gfz._in = dom[i].xn;
    dom[i].Gfz._inb = dom[i].Gfz._in + 2 * DOM_BUF;
    dom[i].Gfz._ie = dom[i].Gfz._is + dom[i].Gfz._in;
    dom[i].Gfz._ieb = dom[i].Gfz._ie + DOM_BUF;

    dom[i].Gfz._js = DOM_BUF;
    dom[i].Gfz._jsb = dom[i].Gfz._js - DOM_BUF;
    dom[i].Gfz._jn = dom[i].yn;
    dom[i].Gfz._jnb = dom[i].Gfz._jn + 2 * DOM_BUF;
    dom[i].Gfz._je = dom[i].Gfz._js + dom[i].Gfz._jn;
    dom[i].Gfz._jeb = dom[i].Gfz._je + DOM_BUF;

    dom[i].Gfz._ks = DOM_BUF;
    dom[i].Gfz._ksb = dom[i].Gfz._ks - DOM_BUF;
    dom[i].Gfz._kn = dom[i].zn + 1;
    dom[i].Gfz._knb = dom[i].Gfz._kn + 2 * DOM_BUF;
    dom[i].Gfz._ke = dom[i].Gfz._ks + dom[i].Gfz._kn;
    dom[i].Gfz._keb = dom[i].Gfz._ke + DOM_BUF;

    dom[i].Gfz._s1 = dom[i].Gfz._in;
    dom[i].Gfz._s2 = dom[i].Gfz._s1 * dom[i].Gfz._jn;
    dom[i].Gfz._s3 = dom[i].Gfz._s2 * dom[i].Gfz._kn;
    dom[i].Gfz._s1b = dom[i].Gfz._inb;
    dom[i].Gfz._s2b = dom[i].Gfz._s1b * dom[i].Gfz._jnb;
    dom[i].Gfz._s3b = dom[i].Gfz._s2b * dom[i].Gfz._knb;
  }

  // set up grid index structs
  // allocate and initialize pressure and velocity vectors
  p0 = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  p = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  phi = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3 * sizeof(real);
  //divU = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u0 = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v0 = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w0 = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
#ifndef IMPLICIT
  diff0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
#endif
  conv0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  diff_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  conv_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_star = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v_star = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w_star = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_WE = (real*) malloc(Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real);
  u_SN_S = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real);
  u_SN_N = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real);
  u_BT_B = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real);
  u_BT_T = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real);
  v_WE_W = (real*) malloc(Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real);
  v_WE_E = (real*) malloc(Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real);
  v_SN = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real);
  v_BT_B = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real);
  v_BT_T = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real);
  w_WE_W = (real*) malloc(Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real);
  w_WE_E = (real*) malloc(Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real);
  w_SN_S = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real);
  w_SN_N = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real);
  w_BT = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real);
  f_x = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  f_y = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  f_z = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);

  // initialize QUIESCENT flow (default)
  for(i = 0; i < Dom.Gcc.s3b; i++) {
    p0[i] = 0.;
    p[i] = 0.;
    phi[i] = 0.;
    //divU[i] = 0.;
  }

  for(i = 0; i < Dom.Gfx.s3b; i++) {
    u[i] = 0.;
#ifndef IMPLICIT
    diff0_u[i] = 0.;
#endif
    diff_u[i] = 0.;
    u0[i] = 0.;
    conv0_u[i] = 0.;
    conv_u[i] = 0.;
    f_x[i] = 0.;
    u_star[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.s3b; i++) {
    v[i] = 0.;
#ifndef IMPLICIT
    diff0_v[i] = 0.;
#endif
    diff_v[i] = 0.;
    v0[i] = 0.;
    conv0_v[i] = 0.;
    conv_v[i] = 0.;
    f_y[i] = 0.;
    v_star[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.s3b; i++) {
    w[i] = 0.;
#ifndef IMPLICIT
    diff0_w[i] = 0.;
#endif
    diff_w[i] = 0.;
    w0[i] = 0.;
    conv0_w[i] = 0.;
    conv_w[i] = 0.;
    f_z[i] = 0.;
    w_star[i] = 0.;
  }

/*
  // initialize TAYLOR_GREEN vortex
  //dt = 0.015199416342412e-4;
  //dt = 0.060562015503876e-4;
  dt = 0.240384615384615e-4/2.;
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        real x = (i-1.0) * Dom.dx + Dom.xs;
        real y = (j-0.5) * Dom.dy + Dom.ys;
        real z = (k-0.5) * Dom.dz + Dom.zs;
        int C = i+j*Dom.Gfx.s1b+k*Dom.Gfx.s2b;
        u[C] = cos(2.*PI*x)*sin(2.*PI*y);
        u0[C] = u[C] * exp(16.*PI*PI*dt);
        conv0_u[C] = -4.*PI*sin(2.*PI*x)*cos(2.*PI*x)*sin(2.*PI*y)*sin(2.*PI*y)
          + 2.*PI*sin(2.*PI*x)*cos(2.*PI*x)
          *(sin(2.*PI*y)*sin(2.*PI*y)-cos(2.*PI*y)*cos(2.*PI*y));
        conv0_u[C] *= exp(16.*PI*PI*dt);
        diff0_u[C] = -8.*PI*PI*nu*cos(2.*PI*x)*sin(2.*PI*y);
        diff0_u[C] *= exp(16.*PI*PI*dt);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        real x = (i-0.5) * Dom.dx + Dom.xs;
        real y = (j-1.0) * Dom.dy + Dom.ys;
        real z = (k-0.5) * Dom.dz + Dom.zs;
        int C = i+j*Dom.Gfy.s1b+k*Dom.Gfy.s2b;
        v[C] = -sin(2.*PI*x)*cos(2.*PI*y);
        v0[C] = v[C] * exp(16.*PI*PI*dt);
        conv0_v[C] = -4.*PI*sin(2.*PI*x)*sin(2.*PI*x)*sin(2.*PI*y)*cos(2.*PI*y)
          + 2.*PI*sin(2.*PI*y)*cos(2.*PI*y)
          *(sin(2.*PI*x)*sin(2.*PI*x)-cos(2.*PI*x)*cos(2.*PI*x));
        conv0_v[C] *= exp(16.*PI*PI*dt);
        diff0_v[C] = 8.*PI*PI*nu*sin(2.*PI*x)*cos(2.*PI*y);
        diff0_v[C] *= exp(16.*PI*PI*dt);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        int C = i+j*Dom.Gfz.s1b+k*Dom.Gfz.s2b;
        w[C] = 0;
        w0[C] = w[C];
      }
    }
  }
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        real x = (i-0.5) * Dom.dx + Dom.xs;
        real y = (j-0.5) * Dom.dy + Dom.ys;
        real z = (k-0.5) * Dom.dz + Dom.zs;
        int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
        p0[C] = -0.25*rho_f*(cos(2.*2.*PI*x)+cos(2.*2.*PI*y))*exp(16.*PI*PI*dt/2.);
        p[C] = -0.25*rho_f*(cos(2.*2.*PI*x)+cos(2.*2.*PI*y))*exp(16.*PI*PI*dt/2.);
      }
    }
  }
*/

  // initialize SHEAR flow
  if(init_cond == SHEAR) {
    for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
      for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
        for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
          real y = ((j-0.5) * Dom.dy) + Dom.ys;
          real z = ((k-0.5) * Dom.dz) + Dom.zs;
          int C = i+j*Dom.Gfx.s1b+k*Dom.Gfx.s2b;
          u[C] = (bc.uNDm-bc.uSDm)*(y-Dom.ys)/Dom.yl + bc.uSDm;
          u[C] += (bc.uTDm-bc.uBDm)*(z-Dom.zs)/Dom.zl + bc.uBDm;
          u0[C] = u[C];
        }
      }
    }
    for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
      for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
        for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
          real x = ((i-0.5) * Dom.dx) + Dom.xs;
          real z = ((k-0.5) * Dom.dz) + Dom.zs;
          int C = i+j*Dom.Gfy.s1b+k*Dom.Gfy.s2b;
          v[C] = (bc.vEDm-bc.vWDm)*(x-Dom.xs)/Dom.xl + bc.vWDm;
          v[C] += (bc.vTDm-bc.vBDm)*(z-Dom.zs)/Dom.ze + bc.vBDm;
          v0[C] = v[C];
        }
      }
    }
    for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
      for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
        for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
          real x = ((i-0.5) * Dom.dx) + Dom.xs;
          real y = ((j-0.5) * Dom.dy) + Dom.ys;
          int C = i+j*Dom.Gfz.s1b+k*Dom.Gfz.s2b;
          w[C] = (bc.wEDm-bc.wWDm)*(x-Dom.xs)/Dom.xl + bc.wWDm;
          w[C] += (bc.wNDm-bc.wSDm)*(y-Dom.ys)/Dom.yl + bc.wSDm;
          w0[C] = w[C];
        }
      }
    }
  } else if(init_cond == CHANNEL) {
    // initialize CHANNEL flow
    for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
      for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
        for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
          real y = ((j-0.5) * Dom.dy) + Dom.ys;
          real z = ((k-0.5) * Dom.dz) + Dom.zs;
          int C = i+j*Dom.Gfx.s1b+k*Dom.Gfx.s2b;
          u[C]  = 0.5/mu*gradP.xm*(y*y-(Dom.ys+Dom.ye)*y+Dom.ys*Dom.ye)
            * (bc.uS == DIRICHLET)
            + 0.5/mu*gradP.xm*(z*z-(Dom.zs+Dom.ze)*z+Dom.zs*Dom.ze)
            * (bc.uB == DIRICHLET);
          u0[C] = u[C];
        }
      }
    }
    for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
      for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
        for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
          real x = ((i-0.5) * Dom.dx) + Dom.xs;
          real z = ((k-0.5) * Dom.dz) + Dom.zs;
          int C = i+j*Dom.Gfy.s1b+k*Dom.Gfy.s2b;
          v[C]  = 0.5/mu*gradP.ym*(x*x-(Dom.xs+Dom.xe)*x+Dom.xs*Dom.xe)
            * (bc.vW == DIRICHLET)
            + 0.5/mu*gradP.ym*(z*z-(Dom.zs+Dom.ze)*z+Dom.zs*Dom.ze)
            * (bc.vB == DIRICHLET);
          v0[C] = v[C];
        }
      }
    }
    for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
      for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
        for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
          real x = ((i-0.5) * Dom.dx) + Dom.xs;
          real y = ((j-0.5) * Dom.dy) + Dom.ys;
          int C = i+j*Dom.Gfz.s1b+k*Dom.Gfz.s2b;
          w[C]  = 0.5/mu*gradP.zm*(x*x-(Dom.xs+Dom.xe)*x+Dom.xs*Dom.xe)
            * (bc.wW == DIRICHLET)
            + 0.5/mu*gradP.zm*(y*y-(Dom.ys+Dom.ye)*y+Dom.ys*Dom.ye)
            * (bc.wS == DIRICHLET);
          w0[C] = w[C];
        }
      }
    }
  }

  for(i = 0; i < Dom.Gfx.jnb*Dom.Gfx.knb; i++) {
    u_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.knb; i++) {
    u_SN_S[i] = 0.;
    u_SN_N[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.jnb; i++) {
    u_BT_B[i] = 0.;
    u_BT_T[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.jnb*Dom.Gfy.knb; i++) {
    v_WE_W[i] = 0.;
    v_WE_E[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.knb; i++) {
    v_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.jnb; i++) {
    v_BT_B[i] = 0.;
    v_BT_T[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.jnb*Dom.Gfz.knb; i++) {
    w_WE_W[i] = 0.;
    w_WE_E[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.knb; i++) {
    w_SN_S[i] = 0.;
    w_SN_N[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.jnb; i++) {
    w_BT[i] = 0.;
  }

  // initialize some variables
  real dx_min = Dom.dx;
  if(Dom.dy < dx_min) dx_min = Dom.dy;
  if(Dom.dz < dx_min) dx_min = Dom.dz;

  dt = 2. * nu / (dx_min * dx_min);
  dt += 2. * nu / (Dom.dy * Dom.dy);
  dt += 2. * nu / (Dom.dz * Dom.dz);
  dt = CFL / dt;
  dt0 = -1.;
  stepnum = 0;
  rec_flow_field_stepnum_out = 0;
  rec_paraview_stepnum_out = 0;
  rec_particle_stepnum_out = 0;
  rec_prec_stepnum_out = 0;
  rec_flow_field_ttime_out = 0;
  rec_paraview_ttime_out = 0;
  rec_particle_ttime_out = 0;
  rec_restart_ttime_out = 0;
  rec_prec_ttime_out = 0;

  turbl = 0.; // integral scale
              // (prevent turbulence linear forcing from being used)

  return EXIT_SUCCESS;
}

int domain_init_turb(void)
{
  int i, j, k;    // iterator
  int C, W, E, S, N, B, T;
  real tmp;

  // make sure there are enough GPU devices in the given range
  if(nsubdom > dev_end - dev_start + 1) {
    return EXIT_FAILURE;
  }

  // calculate domain sizes
  Dom.xl = Dom.xe - Dom.xs;
  Dom.yl = Dom.ye - Dom.ys;
  Dom.zl = Dom.ze - Dom.zs;

  // calculate cell sizes
  Dom.dx = Dom.xl / Dom.xn;
  Dom.dy = Dom.yl / Dom.yn;
  Dom.dz = Dom.zl / Dom.zn;

  // set up grids
  // Gcc
  Dom.Gcc.is = DOM_BUF;
  Dom.Gcc.isb = Dom.Gcc.is - DOM_BUF;
  Dom.Gcc.in = Dom.xn;
  Dom.Gcc.inb = Dom.Gcc.in + 2 * DOM_BUF;
  Dom.Gcc.ie = Dom.Gcc.is + Dom.Gcc.in;
  Dom.Gcc.ieb = Dom.Gcc.ie + DOM_BUF;

  Dom.Gcc.js = DOM_BUF;
  Dom.Gcc.jsb = Dom.Gcc.js - DOM_BUF;
  Dom.Gcc.jn = Dom.yn;
  Dom.Gcc.jnb = Dom.Gcc.jn + 2 * DOM_BUF;
  Dom.Gcc.je = Dom.Gcc.js + Dom.Gcc.jn;
  Dom.Gcc.jeb = Dom.Gcc.je + DOM_BUF;

  Dom.Gcc.ks = DOM_BUF;
  Dom.Gcc.ksb = Dom.Gcc.ks - DOM_BUF;
  Dom.Gcc.kn = Dom.zn;
  Dom.Gcc.knb = Dom.Gcc.kn + 2 * DOM_BUF;
  Dom.Gcc.ke = DOM_BUF + Dom.Gcc.kn;
  Dom.Gcc.keb = Dom.Gcc.ke + DOM_BUF;

  Dom.Gcc.s1 = Dom.Gcc.in;
  Dom.Gcc.s2 = Dom.Gcc.s1 * Dom.Gcc.jn;
  Dom.Gcc.s3 = Dom.Gcc.s2 * Dom.Gcc.kn;
  Dom.Gcc.s1b = Dom.Gcc.inb;
  Dom.Gcc.s2b = Dom.Gcc.s1b * Dom.Gcc.jnb;
  Dom.Gcc.s3b = Dom.Gcc.s2b * Dom.Gcc.knb;

  // Gfx
  Dom.Gfx.is = DOM_BUF;
  Dom.Gfx.isb = Dom.Gfx.is - DOM_BUF;
  Dom.Gfx.in = Dom.xn + 1;
  Dom.Gfx.inb = Dom.Gfx.in + 2 * DOM_BUF;
  Dom.Gfx.ie = Dom.Gfx.is + Dom.Gfx.in;
  Dom.Gfx.ieb = Dom.Gfx.ie + DOM_BUF;

  Dom.Gfx.js = DOM_BUF;
  Dom.Gfx.jsb = Dom.Gfx.js - DOM_BUF;
  Dom.Gfx.jn = Dom.yn;
  Dom.Gfx.jnb = Dom.Gfx.jn + 2 * DOM_BUF;
  Dom.Gfx.je = Dom.Gfx.js + Dom.Gfx.jn;
  Dom.Gfx.jeb = Dom.Gfx.je + DOM_BUF;

  Dom.Gfx.ks = DOM_BUF;
  Dom.Gfx.ksb = Dom.Gfx.ks - DOM_BUF;
  Dom.Gfx.kn = Dom.zn;
  Dom.Gfx.knb = Dom.Gfx.kn + 2 * DOM_BUF;
  Dom.Gfx.ke = Dom.Gfx.ks + Dom.Gfx.kn;
  Dom.Gfx.keb = Dom.Gfx.ke + DOM_BUF;

  Dom.Gfx.s1 = Dom.Gfx.in;
  Dom.Gfx.s2 = Dom.Gfx.s1 * Dom.Gfx.jn;
  Dom.Gfx.s3 = Dom.Gfx.s2 * Dom.Gfx.kn;
  Dom.Gfx.s1b = Dom.Gfx.inb;
  Dom.Gfx.s2b = Dom.Gfx.s1b * Dom.Gfx.jnb;
  Dom.Gfx.s3b = Dom.Gfx.s2b * Dom.Gfx.knb;

  // Gfy
  Dom.Gfy.is = DOM_BUF;
  Dom.Gfy.isb = Dom.Gfy.is - DOM_BUF;
  Dom.Gfy.in = Dom.xn;
  Dom.Gfy.inb = Dom.Gfy.in + 2 * DOM_BUF;
  Dom.Gfy.ie = Dom.Gfy.is + Dom.Gfy.in;
  Dom.Gfy.ieb = Dom.Gfy.ie + DOM_BUF;

  Dom.Gfy.js = DOM_BUF;
  Dom.Gfy.jsb = Dom.Gfy.js - DOM_BUF;
  Dom.Gfy.jn = Dom.yn + 1;
  Dom.Gfy.jnb = Dom.Gfy.jn + 2 * DOM_BUF;
  Dom.Gfy.je = Dom.Gfy.js + Dom.Gfy.jn;
  Dom.Gfy.jeb = Dom.Gfy.je + DOM_BUF;

  Dom.Gfy.ks = DOM_BUF;
  Dom.Gfy.ksb = Dom.Gfy.ks - DOM_BUF;
  Dom.Gfy.kn = Dom.zn;
  Dom.Gfy.knb = Dom.Gfy.kn + 2 * DOM_BUF;
  Dom.Gfy.ke = Dom.Gfy.ks + Dom.Gfy.kn;
  Dom.Gfy.keb = Dom.Gfy.ke + DOM_BUF;

  Dom.Gfy.s1 = Dom.Gfy.in;
  Dom.Gfy.s2 = Dom.Gfy.s1 * Dom.Gfy.jn;
  Dom.Gfy.s3 = Dom.Gfy.s2 * Dom.Gfy.kn;
  Dom.Gfy.s1b = Dom.Gfy.inb;
  Dom.Gfy.s2b = Dom.Gfy.s1b * Dom.Gfy.jnb;
  Dom.Gfy.s3b = Dom.Gfy.s2b * Dom.Gfy.knb;

  // Gfz
  Dom.Gfz.is = DOM_BUF;
  Dom.Gfz.isb = Dom.Gfz.is - DOM_BUF;
  Dom.Gfz.in = Dom.xn;
  Dom.Gfz.inb = Dom.Gfz.in + 2 * DOM_BUF;
  Dom.Gfz.ie = Dom.Gfz.is + Dom.Gfz.in;
  Dom.Gfz.ieb = Dom.Gfz.ie + DOM_BUF;

  Dom.Gfz.js = DOM_BUF;
  Dom.Gfz.jsb = Dom.Gfz.js - DOM_BUF;
  Dom.Gfz.jn = Dom.yn;
  Dom.Gfz.jnb = Dom.Gfz.jn + 2 * DOM_BUF;
  Dom.Gfz.je = Dom.Gfz.js + Dom.Gfz.jn;
  Dom.Gfz.jeb = Dom.Gfz.je + DOM_BUF;

  Dom.Gfz.ks = DOM_BUF;
  Dom.Gfz.ksb = Dom.Gfz.ks - DOM_BUF;
  Dom.Gfz.kn = Dom.zn + 1;
  Dom.Gfz.knb = Dom.Gfz.kn + 2 * DOM_BUF;
  Dom.Gfz.ke = Dom.Gfz.ks + Dom.Gfz.kn;
  Dom.Gfz.keb = Dom.Gfz.ke + DOM_BUF;

  Dom.Gfz.s1 = Dom.Gfz.in;
  Dom.Gfz.s2 = Dom.Gfz.s1 * Dom.Gfz.jn;
  Dom.Gfz.s3 = Dom.Gfz.s2 * Dom.Gfz.kn;
  Dom.Gfz.s1b = Dom.Gfz.inb;
  Dom.Gfz.s2b = Dom.Gfz.s1b * Dom.Gfz.jnb;
  Dom.Gfz.s3b = Dom.Gfz.s2b * Dom.Gfz.knb;

  // initialize subdomains
  for(i = 0; i < nsubdom; i++) {
    dom[i].xl = dom[i].xe - dom[i].xs;
    dom[i].yl = dom[i].ye - dom[i].ys;
    dom[i].zl = dom[i].ze - dom[i].zs;
    dom[i].dx = dom[i].xl / dom[i].xn;
    dom[i].dy = dom[i].yl / dom[i].yn;
    dom[i].dz = dom[i].zl / dom[i].zn;

    // TODO: this algorithm will fail if subdomains are not numbered in
    // increasing order.  This will need to be fixed before going to a machine
    // with more than two GPU's.
    // Gcc
    if(dom[i].W > -1)
      dom[i].Gcc.is = dom[dom[i].W].Gcc.ie;
    else
      dom[i].Gcc.is = DOM_BUF;
    dom[i].Gcc.isb = dom[i].Gcc.is - DOM_BUF;
    dom[i].Gcc.in = dom[i].xn;
    dom[i].Gcc.inb = dom[i].Gcc.in + 2 * DOM_BUF;
    dom[i].Gcc.ie = dom[i].Gcc.is + dom[i].Gcc.in;
    dom[i].Gcc.ieb = dom[i].Gcc.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gcc.js = dom[dom[i].S].Gcc.je;
    else
      dom[i].Gcc.js = DOM_BUF;
    dom[i].Gcc.jsb = dom[i].Gcc.js - DOM_BUF;
    dom[i].Gcc.jn = dom[i].yn;
    dom[i].Gcc.jnb = dom[i].Gcc.jn + 2 * DOM_BUF;
    dom[i].Gcc.je = dom[i].Gcc.js + dom[i].Gcc.jn;
    dom[i].Gcc.jeb = dom[i].Gcc.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gcc.ks = dom[dom[i].B].Gcc.ke;
    else
      dom[i].Gcc.ks = DOM_BUF;
    dom[i].Gcc.ksb = dom[i].Gcc.ks - DOM_BUF;
    dom[i].Gcc.kn = dom[i].zn;
    dom[i].Gcc.knb = dom[i].Gcc.kn + 2 * DOM_BUF;
    dom[i].Gcc.ke = DOM_BUF + dom[i].Gcc.kn;
    dom[i].Gcc.keb = dom[i].Gcc.ke + DOM_BUF;

    dom[i].Gcc.s1 = dom[i].Gcc.in;
    dom[i].Gcc.s2 = dom[i].Gcc.s1 * dom[i].Gcc.jn;
    dom[i].Gcc.s3 = dom[i].Gcc.s2 * dom[i].Gcc.kn;
    dom[i].Gcc.s1b = dom[i].Gcc.inb;
    dom[i].Gcc.s2b = dom[i].Gcc.s1b * dom[i].Gcc.jnb;
    dom[i].Gcc.s3b = dom[i].Gcc.s2b * dom[i].Gcc.knb;

    dom[i].Gcc._is = DOM_BUF;
    dom[i].Gcc._isb = dom[i].Gcc._is - DOM_BUF;
    dom[i].Gcc._in = dom[i].xn;
    dom[i].Gcc._inb = dom[i].Gcc._in + 2 * DOM_BUF;
    dom[i].Gcc._ie = dom[i].Gcc._is + dom[i].Gcc._in;
    dom[i].Gcc._ieb = dom[i].Gcc._ie + DOM_BUF;

    dom[i].Gcc._js = DOM_BUF;
    dom[i].Gcc._jsb = dom[i].Gcc._js - DOM_BUF;
    dom[i].Gcc._jn = dom[i].yn;
    dom[i].Gcc._jnb = dom[i].Gcc._jn + 2 * DOM_BUF;
    dom[i].Gcc._je = dom[i].Gcc._js + dom[i].Gcc._jn;
    dom[i].Gcc._jeb = dom[i].Gcc._je + DOM_BUF;

    dom[i].Gcc._ks = DOM_BUF;
    dom[i].Gcc._ksb = dom[i].Gcc._ks - DOM_BUF;
    dom[i].Gcc._kn = dom[i].zn;
    dom[i].Gcc._knb = dom[i].Gcc._kn + 2 * DOM_BUF;
    dom[i].Gcc._ke = dom[i].Gcc._ks + dom[i].Gcc._kn;
    dom[i].Gcc._keb = dom[i].Gcc._ke + DOM_BUF;

    dom[i].Gcc._s1 = dom[i].Gcc._in;
    dom[i].Gcc._s2 = dom[i].Gcc._s1 * dom[i].Gcc._jn;
    dom[i].Gcc._s3 = dom[i].Gcc._s2 * dom[i].Gcc._kn;
    dom[i].Gcc._s1b = dom[i].Gcc._inb;
    dom[i].Gcc._s2b = dom[i].Gcc._s1b * dom[i].Gcc._jnb;
    dom[i].Gcc._s3b = dom[i].Gcc._s2b * dom[i].Gcc._knb;

    // Gfx
    if(dom[i].W > -1)
      dom[i].Gfx.is = dom[dom[i].W].Gfx.ie - 1;
    else
      dom[i].Gfx.is = DOM_BUF;
    dom[i].Gfx.isb = dom[i].Gfx.is - DOM_BUF;
    dom[i].Gfx.in = dom[i].xn + 1;
    dom[i].Gfx.inb = dom[i].Gfx.in + 2 * DOM_BUF;
    dom[i].Gfx.ie = dom[i].Gfx.is + dom[i].Gfx.in;
    dom[i].Gfx.ieb = dom[i].Gfx.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfx.js = dom[dom[i].S].Gfx.je - 1;
    else
      dom[i].Gfx.js = DOM_BUF;
    dom[i].Gfx.jsb = dom[i].Gfx.js - DOM_BUF;
    dom[i].Gfx.jn = dom[i].yn;
    dom[i].Gfx.jnb = dom[i].Gfx.jn + 2 * DOM_BUF;
    dom[i].Gfx.je = dom[i].Gfx.js + dom[i].Gfx.jn;
    dom[i].Gfx.jeb = dom[i].Gfx.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfx.ks = dom[dom[i].B].Gfx.ke - 1;
    else
      dom[i].Gfx.ks = DOM_BUF;
    dom[i].Gfx.ksb = dom[i].Gfx.ks - DOM_BUF;
    dom[i].Gfx.kn = dom[i].zn;
    dom[i].Gfx.knb = dom[i].Gfx.kn + 2 * DOM_BUF;
    dom[i].Gfx.ke = dom[i].Gfx.ks + dom[i].Gfx.kn;
    dom[i].Gfx.keb = dom[i].Gfx.ke + DOM_BUF;

    dom[i].Gfx.s1 = dom[i].Gfx.in;
    dom[i].Gfx.s2 = dom[i].Gfx.s1 * dom[i].Gfx.jn;
    dom[i].Gfx.s3 = dom[i].Gfx.s2 * dom[i].Gfx.kn;
    dom[i].Gfx.s1b = dom[i].Gfx.inb;
    dom[i].Gfx.s2b = dom[i].Gfx.s1b * dom[i].Gfx.jnb;
    dom[i].Gfx.s3b = dom[i].Gfx.s2b * dom[i].Gfx.knb;

    dom[i].Gfx._is = DOM_BUF;
    dom[i].Gfx._isb = dom[i].Gfx._is - DOM_BUF;
    dom[i].Gfx._in = dom[i].xn + 1;
    dom[i].Gfx._inb = dom[i].Gfx._in + 2 * DOM_BUF;
    dom[i].Gfx._ie = dom[i].Gfx._is + dom[i].Gfx._in;
    dom[i].Gfx._ieb = dom[i].Gfx._ie + DOM_BUF;

    dom[i].Gfx._js = DOM_BUF;
    dom[i].Gfx._jsb = dom[i].Gfx._js - DOM_BUF;
    dom[i].Gfx._jn = dom[i].yn;
    dom[i].Gfx._jnb = dom[i].Gfx._jn + 2 * DOM_BUF;
    dom[i].Gfx._je = dom[i].Gfx._js + dom[i].Gfx._jn;
    dom[i].Gfx._jeb = dom[i].Gfx._je + DOM_BUF;

    dom[i].Gfx._ks = DOM_BUF;
    dom[i].Gfx._ksb = dom[i].Gfx._ks - DOM_BUF;
    dom[i].Gfx._kn = dom[i].zn;
    dom[i].Gfx._knb = dom[i].Gfx._kn + 2 * DOM_BUF;
    dom[i].Gfx._ke = dom[i].Gfx._ks + dom[i].Gfx._kn;
    dom[i].Gfx._keb = dom[i].Gfx._ke + DOM_BUF;

    dom[i].Gfx._s1 = dom[i].Gfx._in;
    dom[i].Gfx._s2 = dom[i].Gfx._s1 * dom[i].Gfx._jn;
    dom[i].Gfx._s3 = dom[i].Gfx._s2 * dom[i].Gfx._kn;
    dom[i].Gfx._s1b = dom[i].Gfx._inb;
    dom[i].Gfx._s2b = dom[i].Gfx._s1b * dom[i].Gfx._jnb;
    dom[i].Gfx._s3b = dom[i].Gfx._s2b * dom[i].Gfx._knb;

    // Gfy
    if(dom[i].W > -1)
      dom[i].Gfy.is = dom[dom[i].W].Gfy.ie;
    else
      dom[i].Gfy.is = DOM_BUF;
    dom[i].Gfy.isb = dom[i].Gfy.is - DOM_BUF;
    dom[i].Gfy.in = dom[i].xn;
    dom[i].Gfy.inb = dom[i].Gfy.in + 2 * DOM_BUF;
    dom[i].Gfy.ie = dom[i].Gfy.is + dom[i].Gfy.in;
    dom[i].Gfy.ieb = dom[i].Gfy.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfy.js = dom[dom[i].S].Gfy.je;
    else
      dom[i].Gfy.js = DOM_BUF;
    dom[i].Gfy.jsb = dom[i].Gfy.js - DOM_BUF;
    dom[i].Gfy.jn = dom[i].yn + 1;
    dom[i].Gfy.jnb = dom[i].Gfy.jn + 2 * DOM_BUF;
    dom[i].Gfy.je = dom[i].Gfy.js + dom[i].Gfy.jn;
    dom[i].Gfy.jeb = dom[i].Gfy.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfy.ks = dom[dom[i].B].Gfy.ke;
    else
      dom[i].Gfy.ks = DOM_BUF;
    dom[i].Gfy.ksb = dom[i].Gfy.ks - DOM_BUF;
    dom[i].Gfy.kn = dom[i].zn;
    dom[i].Gfy.knb = dom[i].Gfy.kn + 2 * DOM_BUF;
    dom[i].Gfy.ke = dom[i].Gfy.ks + dom[i].Gfy.kn;
    dom[i].Gfy.keb = dom[i].Gfy.ke + DOM_BUF;

    dom[i].Gfy.s1 = dom[i].Gfy.in;
    dom[i].Gfy.s2 = dom[i].Gfy.s1 * dom[i].Gfy.jn;
    dom[i].Gfy.s3 = dom[i].Gfy.s2 * dom[i].Gfy.kn;
    dom[i].Gfy.s1b = dom[i].Gfy.inb;
    dom[i].Gfy.s2b = dom[i].Gfy.s1b * dom[i].Gfy.jnb;
    dom[i].Gfy.s3b = dom[i].Gfy.s2b * dom[i].Gfy.knb;

    dom[i].Gfy._is = DOM_BUF;
    dom[i].Gfy._isb = dom[i].Gfy._is - DOM_BUF;
    dom[i].Gfy._in = dom[i].xn;
    dom[i].Gfy._inb = dom[i].Gfy._in + 2 * DOM_BUF;
    dom[i].Gfy._ie = dom[i].Gfy._is + dom[i].Gfy._in;
    dom[i].Gfy._ieb = dom[i].Gfy._ie + DOM_BUF;

    dom[i].Gfy._js = DOM_BUF;
    dom[i].Gfy._jsb = dom[i].Gfy._js - DOM_BUF;
    dom[i].Gfy._jn = dom[i].yn + 1;
    dom[i].Gfy._jnb = dom[i].Gfy._jn + 2 * DOM_BUF;
    dom[i].Gfy._je = dom[i].Gfy._js + dom[i].Gfy._jn;
    dom[i].Gfy._jeb = dom[i].Gfy._je + DOM_BUF;

    dom[i].Gfy._ks = DOM_BUF;
    dom[i].Gfy._ksb = dom[i].Gfy._ks - DOM_BUF;
    dom[i].Gfy._kn = dom[i].zn;
    dom[i].Gfy._knb = dom[i].Gfy._kn + 2 * DOM_BUF;
    dom[i].Gfy._ke = dom[i].Gfy._ks + dom[i].Gfy._kn;
    dom[i].Gfy._keb = dom[i].Gfy._ke + DOM_BUF;

    dom[i].Gfy._s1 = dom[i].Gfy._in;
    dom[i].Gfy._s2 = dom[i].Gfy._s1 * dom[i].Gfy._jn;
    dom[i].Gfy._s3 = dom[i].Gfy._s2 * dom[i].Gfy._kn;
    dom[i].Gfy._s1b = dom[i].Gfy._inb;
    dom[i].Gfy._s2b = dom[i].Gfy._s1b * dom[i].Gfy._jnb;
    dom[i].Gfy._s3b = dom[i].Gfy._s2b * dom[i].Gfy._knb;

    // Gfz
    if(dom[i].W > -1)
      dom[i].Gfz.is = dom[dom[i].W].Gfz.ie;
    else
      dom[i].Gfz.is = DOM_BUF;
    dom[i].Gfz.isb = dom[i].Gfz.is - DOM_BUF;
    dom[i].Gfz.in = dom[i].xn;
    dom[i].Gfz.inb = dom[i].Gfz.in + 2 * DOM_BUF;
    dom[i].Gfz.ie = dom[i].Gfz.is + dom[i].Gfz.in;
    dom[i].Gfz.ieb = dom[i].Gfz.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfz.js = dom[dom[i].S].Gfz.je;
    else
      dom[i].Gfz.js = DOM_BUF;
    dom[i].Gfz.jsb = dom[i].Gfz.js - DOM_BUF;
    dom[i].Gfz.jn = dom[i].yn;
    dom[i].Gfz.jnb = dom[i].Gfz.jn + 2 * DOM_BUF;
    dom[i].Gfz.je = dom[i].Gfz.js + dom[i].Gfz.jn;
    dom[i].Gfz.jeb = dom[i].Gfz.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfz.ks = dom[dom[i].B].Gfz.ke;
    else
      dom[i].Gfz.ks = DOM_BUF;
    dom[i].Gfz.ksb = dom[i].Gfz.ks - DOM_BUF;
    dom[i].Gfz.kn = dom[i].zn + 1;
    dom[i].Gfz.knb = dom[i].Gfz.kn + 2 * DOM_BUF;
    dom[i].Gfz.ke = dom[i].Gfz.ks + dom[i].Gfz.kn;
    dom[i].Gfz.keb = dom[i].Gfz.ke + DOM_BUF;

    dom[i].Gfz.s1 = dom[i].Gfz.in;
    dom[i].Gfz.s2 = dom[i].Gfz.s1 * dom[i].Gfz.jn;
    dom[i].Gfz.s3 = dom[i].Gfz.s2 * dom[i].Gfz.kn;
    dom[i].Gfz.s1b = dom[i].Gfz.inb;
    dom[i].Gfz.s2b = dom[i].Gfz.s1b * dom[i].Gfz.jnb;
    dom[i].Gfz.s3b = dom[i].Gfz.s2b * dom[i].Gfz.knb;

    dom[i].Gfz._is = DOM_BUF;
    dom[i].Gfz._isb = dom[i].Gfz._is - DOM_BUF;
    dom[i].Gfz._in = dom[i].xn;
    dom[i].Gfz._inb = dom[i].Gfz._in + 2 * DOM_BUF;
    dom[i].Gfz._ie = dom[i].Gfz._is + dom[i].Gfz._in;
    dom[i].Gfz._ieb = dom[i].Gfz._ie + DOM_BUF;

    dom[i].Gfz._js = DOM_BUF;
    dom[i].Gfz._jsb = dom[i].Gfz._js - DOM_BUF;
    dom[i].Gfz._jn = dom[i].yn;
    dom[i].Gfz._jnb = dom[i].Gfz._jn + 2 * DOM_BUF;
    dom[i].Gfz._je = dom[i].Gfz._js + dom[i].Gfz._jn;
    dom[i].Gfz._jeb = dom[i].Gfz._je + DOM_BUF;

    dom[i].Gfz._ks = DOM_BUF;
    dom[i].Gfz._ksb = dom[i].Gfz._ks - DOM_BUF;
    dom[i].Gfz._kn = dom[i].zn + 1;
    dom[i].Gfz._knb = dom[i].Gfz._kn + 2 * DOM_BUF;
    dom[i].Gfz._ke = dom[i].Gfz._ks + dom[i].Gfz._kn;
    dom[i].Gfz._keb = dom[i].Gfz._ke + DOM_BUF;

    dom[i].Gfz._s1 = dom[i].Gfz._in;
    dom[i].Gfz._s2 = dom[i].Gfz._s1 * dom[i].Gfz._jn;
    dom[i].Gfz._s3 = dom[i].Gfz._s2 * dom[i].Gfz._kn;
    dom[i].Gfz._s1b = dom[i].Gfz._inb;
    dom[i].Gfz._s2b = dom[i].Gfz._s1b * dom[i].Gfz._jnb;
    dom[i].Gfz._s3b = dom[i].Gfz._s2b * dom[i].Gfz._knb;
  }

  // set up grid index structs
  // allocate and initialize pressure and velocity vectors
  p0 = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  p = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  phi = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3 * sizeof(real);
  //divU = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u0 = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v0 = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w0 = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
#ifndef IMPLICIT
  diff0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
#endif
  conv0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  diff_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  conv_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_star = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v_star = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w_star = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_WE = (real*) malloc(Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real);
  u_SN_S = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real);
  u_SN_N = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real);
  u_BT_B = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real);
  u_BT_T = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real);
  v_WE_W = (real*) malloc(Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real);
  v_WE_E = (real*) malloc(Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real);
  v_SN = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real);
  v_BT_B = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real);
  v_BT_T = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real);
  w_WE_W = (real*) malloc(Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real);
  w_WE_E = (real*) malloc(Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real);
  w_SN_S = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real);
  w_SN_N = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real);
  w_BT = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real);
  f_x = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  f_y = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  f_z = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);

  // set up the random number generator
  srand(time(NULL));

  for(i = 0; i < Dom.Gcc.s3b; i++) {
    p0[i] = 0.;
    p[i] = 0.;
    phi[i] = 0.;
    //divU[i] = 0.;
  }

  for(i = 0; i < Dom.Gfx.s3b; i++) {
    u[i] = 0.;
#ifndef IMPLICIT
    diff0_u[i] = 0.;
#endif
    diff_u[i] = 0.;
    u0[i] = 0.;
    conv0_u[i] = 0.;
    conv_u[i] = 0.;
    f_x[i] = 0.;
    u_star[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.s3b; i++) {
    v[i] = 0.;
#ifndef IMPLICIT
    diff0_v[i] = 0.;
#endif
    diff_v[i] = 0.;
    v0[i] = 0.;
    conv0_v[i] = 0.;
    conv_v[i] = 0.;
    f_y[i] = 0.;
    v_star[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.s3b; i++) {
    w[i] = 0.;
#ifndef IMPLICIT
    diff0_w[i] = 0.;
#endif
    diff_w[i] = 0.;
    w0[i] = 0.;
    conv0_w[i] = 0.;
    conv_w[i] = 0.;
    f_z[i] = 0.;
    w_star[i] = 0.;
  }

  // integral scale
  turbl = (Dom.xl + Dom.yl + Dom.zl) / 3.;
  real urms = 3*turbA*turbl;

  // randomly initialize velocity components
  for(i = 0; i < Dom.Gfx.s3b; i++) {
    tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
    u[i] = tmp;
    u0[i] = tmp;
  }
  for(i = 0; i < Dom.Gfy.s3b; i++) {
    tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
    v[i] = tmp;
    v0[i] = tmp;
  }
  for(i = 0; i < Dom.Gfz.s3b; i++) {
    tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
    w[i] = tmp;
    w0[i] = tmp;
  }

  // calculate the divergence of U
  real vol = Dom.xn*Dom.yn*Dom.zn;
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        W = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        E = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        S = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        N = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        B = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        T = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p[C] = (u[E]-u[W])/Dom.dx + (v[N]-v[S])/Dom.dy + (w[T]-w[B])/Dom.dz;
        p[C] = p[C] / vol;
      }
    }
  }

  real umean = 0.;
  real vmean = 0.;
  real wmean = 0.;
  // subtract off the divergence of U to make U solenoidal
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        W = (i-1) + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        E = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        u[C] = u[C] - 0.5*(p[W] + p[E]);
        u0[C] = u0[C] - 0.5*(p[W] + p[E]);
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        S = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        N = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = v[C] - 0.5*(p[S] + p[N]);
        v0[C] = v0[C] - 0.5*(p[S] + p[N]);
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        B = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
        T = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = w[C] - 0.5*(p[B] + p[T]);
        w0[C] = w0[C] - 0.5*(p[B] + p[T]);
      }
    }
  }

  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        umean += u[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        vmean += v[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        wmean += w[C];
      }
    }
  } 

  umean /= vol;
  vmean /= vol;
  wmean /= vol;

  // re-scale to give zero mean velocity in each direction
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = u[C] - umean;
        u0[C] = u0[C] - umean;
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = v[C] - vmean;
        v0[C] = v0[C] - vmean;
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = w[C] - wmean;
        w0[C] = w0[C] - wmean;
      }
    }
  }

  for(i = 0; i < Dom.Gfx.jnb*Dom.Gfx.knb; i++) {
    u_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.knb; i++) {
    u_SN_S[i] = 0.;
    u_SN_N[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.jnb; i++) {
    u_BT_B[i] = 0.;
    u_BT_T[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.jnb*Dom.Gfy.knb; i++) {
    v_WE_W[i] = 0.;
    v_WE_E[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.knb; i++) {
    v_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.jnb; i++) {
    v_BT_B[i] = 0.;
    v_BT_T[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.jnb*Dom.Gfz.knb; i++) {
    w_WE_W[i] = 0.;
    w_WE_E[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.knb; i++) {
    w_SN_S[i] = 0.;
    w_SN_N[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.jnb; i++) {
    w_BT[i] = 0.;
  }

  // initialize some variables
  real dx_min = Dom.dx;
  if(Dom.dy < dx_min) dx_min = Dom.dy;
  if(Dom.dz < dx_min) dx_min = Dom.dz;

  dt = 2. * nu / (dx_min * dx_min);
  //dt += 2. * nu / (Dom.dy * Dom.dy);
  //dt += 2. * nu / (Dom.dz * Dom.dz);
  dt = CFL / dt;
  dt0 = -1.;
  stepnum = 0;
  rec_flow_field_stepnum_out = 0;
  rec_paraview_stepnum_out = 0;
  rec_particle_stepnum_out = 0;
  rec_prec_stepnum_out = 0;
  rec_flow_field_ttime_out = 0;
  rec_paraview_ttime_out = 0;
  rec_particle_ttime_out = 0;
  rec_restart_ttime_out = 0;
  rec_prec_ttime_out = 0;

  return EXIT_SUCCESS;
}

void domain_clean(void)
{
  free(dom);
  free(p0);
  free(p);
  free(phi);
  //free(divU);
  free(u);
  free(v);
  free(w);
  free(u0);
  free(v0);
  free(w0);
#ifndef IMPLICIT
  free(diff0_u);
  free(diff0_v);
  free(diff0_w);
#endif
  free(conv0_u);
  free(conv0_v);
  free(conv0_w);
  free(diff_u);
  free(diff_v);
  free(diff_w);
  free(conv_u);
  free(conv_v);
  free(conv_w);
  free(u_star);
  free(v_star);
  free(w_star);
  free(u_WE);
  free(u_SN_S);
  free(u_SN_N);
  free(u_BT_B);
  free(u_BT_T);
  free(v_WE_W);
  free(v_WE_E);
  free(v_SN);
  free(v_BT_B);
  free(v_BT_T);
  free(w_WE_W);
  free(w_WE_E);
  free(w_SN_S);
  free(w_SN_N);
  free(w_BT);
  free(f_x);
  free(f_y);
  free(f_z);
}

void compute_vel_BC(void)
{
  real delta = ttime - vel_tDelay;
  if (delta >= 0) { // this works since bc.u*D is initialized to zero
    // uWD
    if(bc.uWDa == 0) bc.uWD = bc.uWDm;
    else if(fabs(delta*bc.uWDa) > fabs(bc.uWDm)) bc.uWD = bc.uWDm;
    else bc.uWD = delta*bc.uWDa;
    // uED
    if(bc.uEDa == 0) bc.uED = bc.uEDm;
    else if(fabs(delta*bc.uEDa) > fabs(bc.uEDm)) bc.uED = bc.uEDm;
    else bc.uED = delta*bc.uEDa;
    // uSD
    if(bc.uSDa == 0) bc.uSD = bc.uSDm;
    else if(fabs(delta*bc.uSDa) > fabs(bc.uSDm)) bc.uSD = bc.uSDm;
    else bc.uSD = delta*bc.uSDa;
    // uND
    if(bc.uNDa == 0) bc.uND = bc.uNDm;
    else if(fabs(delta*bc.uNDa) > fabs(bc.uNDm)) bc.uND = bc.uNDm;
    else bc.uND = delta*bc.uNDa;
    // uBD
    if(bc.uBDa == 0) bc.uBD = bc.uBDm;
    else if(fabs(delta*bc.uBDa) > fabs(bc.uBDm)) bc.uBD = bc.uBDm;
    else bc.uBD = delta*bc.uBDa;
    // uTD
    if(bc.uTDa == 0) bc.uTD = bc.uTDm;
    else if(fabs(delta*bc.uTDa) > fabs(bc.uTDm)) bc.uTD = bc.uTDm;
    else bc.uTD = delta*bc.uTDa;
    // vWD
    if(bc.vWDa == 0) bc.vWD = bc.vWDm;
    else if(fabs(delta*bc.vWDa) > fabs(bc.vWDm)) bc.vWD = bc.vWDm;
    else bc.vWD = delta*bc.vWDa;
    // vED
    if(bc.vEDa == 0) bc.vED = bc.vEDm;
    else if(fabs(delta*bc.vEDa) > fabs(bc.vEDm)) bc.vED = bc.vEDm;
    else bc.vED = delta*bc.vEDa;
    // vSD
    if(bc.vSDa == 0) bc.vSD = bc.vSDm;
    else if(fabs(delta*bc.vSDa) > fabs(bc.vSDm)) bc.vSD = bc.vSDm;
    else bc.vSD = delta*bc.vSDa;
    // vND
    if(bc.vNDa == 0) bc.vND = bc.vNDm;
    else if(fabs(delta*bc.vNDa) > fabs(bc.vNDm)) bc.vND = bc.vNDm;
    else bc.vND = delta*bc.vNDa;
    // vBD
    if(bc.vBDa == 0) bc.vBD = bc.vBDm;
    else if(fabs(delta*bc.vBDa) > fabs(bc.vBDm)) bc.vBD = bc.vBDm;
    else bc.vBD = delta*bc.vBDa;
    // vTD
    if(bc.vTDa == 0) bc.vTD = bc.vTDm;
    else if(fabs(delta*bc.vTDa) > fabs(bc.vTDm)) bc.vTD = bc.vTDm;
    else bc.vTD = delta*bc.vTDa;
    // wWD
    if(bc.wWDa == 0) bc.wWD = bc.wWDm;
    else if(fabs(delta*bc.wWDa) > fabs(bc.wWDm)) bc.wWD = bc.wWDm;
    else bc.wWD = delta*bc.wWDa;
    // wED
    if(bc.wEDa == 0) bc.wED = bc.wEDm;
    else if(fabs(delta*bc.wEDa) > fabs(bc.wEDm)) bc.wED = bc.wEDm;
    else bc.wED = delta*bc.wEDa;
    // wSD
    if(bc.wSDa == 0) bc.wSD = bc.wSDm;
    else if(fabs(delta*bc.wSDa) > fabs(bc.wSDm)) bc.wSD = bc.wSDm;
    else bc.wSD = delta*bc.wSDa;
    // wND
    if(bc.wNDa == 0) bc.wND = bc.wNDm;
    else if(fabs(delta*bc.wNDa) > fabs(bc.wNDm)) bc.wND = bc.wNDm;
    else bc.wND = delta*bc.wNDa;
    // wBD
    if(bc.wBDa == 0) bc.wBD = bc.wBDm;
    else if(fabs(delta*bc.wBDa) > fabs(bc.wBDm)) bc.wBD = bc.wBDm;
    else bc.wBD = delta*bc.wBDa;
    // wTD
    if(bc.wTDa == 0) bc.wTD = bc.wTDm;
    else if(fabs(delta*bc.wTDa) > fabs(bc.wTDm)) bc.wTD = bc.wTDm;
    else bc.wTD = delta*bc.wTDa;
  }
}

void out_restart(void)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/input/restart.config", ROOT_DIR);
  FILE *rest = fopen(path, "w");
  if(rest == NULL) {
    fprintf(stderr, "Could not open file restart.input.\n");
    exit(EXIT_FAILURE);
  }

  fwrite(&ttime, sizeof(real), 1, rest);
  fwrite(&dt0, sizeof(real), 1, rest);
  fwrite(&dt, sizeof(real), 1, rest);
  fwrite(&stepnum, sizeof(int), 1, rest);
  fwrite(&rec_paraview_stepnum_out, sizeof(int), 1, rest);

  fwrite(u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u0, sizeof(real), Dom.Gfx.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_u, sizeof(real), Dom.Gfx.s3b, rest);
#endif
  fwrite(conv0_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(diff_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(conv_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u_star, sizeof(real), Dom.Gfx.s3b, rest);

  fwrite(v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v0, sizeof(real), Dom.Gfy.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_v, sizeof(real), Dom.Gfy.s3b, rest);
#endif
  fwrite(conv0_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(diff_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(conv_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v_star, sizeof(real), Dom.Gfy.s3b, rest);

  fwrite(w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w0, sizeof(real), Dom.Gfz.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_w, sizeof(real), Dom.Gfz.s3b, rest);
#endif
  fwrite(conv0_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(diff_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(conv_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w_star, sizeof(real), Dom.Gfz.s3b, rest);

  fwrite(p, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phi, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(p0, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phase, sizeof(int), Dom.Gcc.s3b, rest);
  fwrite(phase_shell, sizeof(int), Dom.Gcc.s3b, rest);

  fwrite(flag_u, sizeof(int), Dom.Gfx.s3b, rest);
  fwrite(flag_v, sizeof(int), Dom.Gfy.s3b, rest);
  fwrite(flag_w, sizeof(int), Dom.Gfz.s3b, rest);

  fwrite(bc_plane_pos, sizeof(real), 9, rest);

  fwrite(parts, sizeof(part_struct), nparts, rest);

  fwrite(&coeff_stride, sizeof(int), 1, rest);

  fwrite(pnm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(pnm_im00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(phinm_im00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im0, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_re00, sizeof(real), nparts*coeff_stride, rest);
  fwrite(chinm_im00, sizeof(real), nparts*coeff_stride, rest);

  fwrite(&rec_flow_field_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_paraview_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_particle_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_restart_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_prec_ttime_out, sizeof(real), 1, rest);
  fwrite(&pid_int, sizeof(real), 1, rest);
  fwrite(&pid_back, sizeof(real), 1, rest);
  fwrite(&gradP.z, sizeof(real), 2, rest);

  // close the file
  fclose(rest);
}

void in_restart(void)
{
  int fret = 0;
  fret = fret; // prevent compiler warning
  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/restart.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  fret = fread(&ttime, sizeof(real), 1, infile);
  fret = fread(&dt0, sizeof(real), 1, infile);
  fret = fread(&dt, sizeof(real), 1, infile);
  fret = fread(&stepnum, sizeof(int), 1, infile);
  fret = fread(&rec_paraview_stepnum_out, sizeof(int), 1, infile);

  fret = fread(u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(u0, sizeof(real), Dom.Gfx.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_u, sizeof(real), Dom.Gfx.s3b, infile);
#endif
  fret = fread(conv0_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(diff_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(conv_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(u_star, sizeof(real), Dom.Gfx.s3b, infile);

  fret = fread(v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(v0, sizeof(real), Dom.Gfy.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_v, sizeof(real), Dom.Gfy.s3b, infile);
#endif
  fret = fread(conv0_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(diff_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(conv_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(v_star, sizeof(real), Dom.Gfy.s3b, infile);

  fret = fread(w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(w0, sizeof(real), Dom.Gfz.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_w, sizeof(real), Dom.Gfz.s3b, infile);
#endif
  fret = fread(conv0_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(diff_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(conv_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(w_star, sizeof(real), Dom.Gfz.s3b, infile);

  fret = fread(p, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(phi, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(p0, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(phase, sizeof(int), Dom.Gcc.s3b, infile);
  fret = fread(phase_shell, sizeof(int), Dom.Gcc.s3b, infile);

  fret = fread(flag_u, sizeof(int), Dom.Gfx.s3b, infile);
  fret = fread(flag_v, sizeof(int), Dom.Gfy.s3b, infile);
  fret = fread(flag_w, sizeof(int), Dom.Gfz.s3b, infile);

  fret = fread(bc_plane_pos, sizeof(real), 9, infile);

  fret = fread(parts, sizeof(part_struct), nparts, infile);

  fret = fread(&coeff_stride, sizeof(int), 1, infile);

  fret = fread(pnm_re, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_im, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_re0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_im0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_re00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(pnm_im00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_re, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_im, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_re0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_im0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_re00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(phinm_im00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_re, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_im, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_re0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_im0, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_re00, sizeof(real), nparts*coeff_stride, infile);
  fret = fread(chinm_im00, sizeof(real), nparts*coeff_stride, infile);

  fret = fread(&rec_flow_field_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_paraview_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_particle_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_restart_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_prec_ttime_out, sizeof(real), 1, infile);
  fret = fread(&pid_int, sizeof(real), 1, infile);
  fret = fread(&pid_back, sizeof(real), 1, infile);
  fret = fread(&gradP.z, sizeof(real), 1, infile);

  // close file
  fclose(infile);
}

void out_restart_turb(void)
{
  // create the file
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/input/restart_turb.config", ROOT_DIR);
  FILE *rest = fopen(path, "w");
  if(rest == NULL) {
    fprintf(stderr, "Could not open file restart_turb.input.\n");
    exit(EXIT_FAILURE);
  }

  fwrite(&ttime, sizeof(real), 1, rest);
  fwrite(&dt0, sizeof(real), 1, rest);
  fwrite(&dt, sizeof(real), 1, rest);
  fwrite(&stepnum, sizeof(int), 1, rest);
  fwrite(&rec_paraview_stepnum_out, sizeof(int), 1, rest);

  fwrite(u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u0, sizeof(real), Dom.Gfx.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_u, sizeof(real), Dom.Gfx.s3b, rest);
#endif
  fwrite(conv0_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(diff_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(conv_u, sizeof(real), Dom.Gfx.s3b, rest);
  fwrite(u_star, sizeof(real), Dom.Gfx.s3b, rest);

  fwrite(v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v0, sizeof(real), Dom.Gfy.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_v, sizeof(real), Dom.Gfy.s3b, rest);
#endif
  fwrite(conv0_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(diff_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(conv_v, sizeof(real), Dom.Gfy.s3b, rest);
  fwrite(v_star, sizeof(real), Dom.Gfy.s3b, rest);

  fwrite(w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w0, sizeof(real), Dom.Gfz.s3b, rest);
#ifndef IMPLICIT
  fwrite(diff0_w, sizeof(real), Dom.Gfz.s3b, rest);
#endif
  fwrite(conv0_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(diff_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(conv_w, sizeof(real), Dom.Gfz.s3b, rest);
  fwrite(w_star, sizeof(real), Dom.Gfz.s3b, rest);

  fwrite(p, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phi, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(p0, sizeof(real), Dom.Gcc.s3b, rest);
  fwrite(phase, sizeof(int), Dom.Gcc.s3b, rest);
  fwrite(phase_shell, sizeof(int), Dom.Gcc.s3b, rest);

  fwrite(flag_u, sizeof(int), Dom.Gfx.s3b, rest);
  fwrite(flag_v, sizeof(int), Dom.Gfy.s3b, rest);
  fwrite(flag_w, sizeof(int), Dom.Gfz.s3b, rest);

  fwrite(bc_plane_pos, sizeof(real), 9, rest);

  fwrite(&rec_prec_ttime_out, sizeof(real), 1, rest);
  fwrite(&rec_prec_flow_field_ttime_out, sizeof(real), 1, rest);
  // close the file
  fclose(rest);
}

void in_restart_turb(void)
{
  int fret = 0;
  fret = fret; // prevent compiler warning
  // open configuration file for reading
  char fname[FILE_NAME_SIZE] = "";
  sprintf(fname, "%s/input/restart_turb.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  fret = fread(&ttime, sizeof(real), 1, infile);
  fret = fread(&dt0, sizeof(real), 1, infile);
  fret = fread(&dt, sizeof(real), 1, infile);
  fret = fread(&stepnum, sizeof(int), 1, infile);
  fret = fread(&rec_paraview_stepnum_out, sizeof(int), 1, infile);

  fret = fread(u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(u0, sizeof(real), Dom.Gfx.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_u, sizeof(real), Dom.Gfx.s3b, infile);
#endif
  fret = fread(conv0_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(diff_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(conv_u, sizeof(real), Dom.Gfx.s3b, infile);
  fret = fread(u_star, sizeof(real), Dom.Gfx.s3b, infile);

  fret = fread(v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(v0, sizeof(real), Dom.Gfy.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_v, sizeof(real), Dom.Gfy.s3b, infile);
#endif
  fret = fread(conv0_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(diff_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(conv_v, sizeof(real), Dom.Gfy.s3b, infile);
  fret = fread(v_star, sizeof(real), Dom.Gfy.s3b, infile);

  fret = fread(w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(w0, sizeof(real), Dom.Gfz.s3b, infile);
#ifndef IMPLICIT
  fret = fread(diff0_w, sizeof(real), Dom.Gfz.s3b, infile);
#endif
  fret = fread(conv0_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(diff_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(conv_w, sizeof(real), Dom.Gfz.s3b, infile);
  fret = fread(w_star, sizeof(real), Dom.Gfz.s3b, infile);

  fret = fread(p, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(phi, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(p0, sizeof(real), Dom.Gcc.s3b, infile);
  fret = fread(phase, sizeof(int), Dom.Gcc.s3b, infile);
  fret = fread(phase_shell, sizeof(int), Dom.Gcc.s3b, infile);

  fret = fread(flag_u, sizeof(int), Dom.Gfx.s3b, infile);
  fret = fread(flag_v, sizeof(int), Dom.Gfy.s3b, infile);
  fret = fread(flag_w, sizeof(int), Dom.Gfz.s3b, infile);

  fret = fread(bc_plane_pos, sizeof(real), 9, infile);

  fret = fread(&rec_prec_ttime_out, sizeof(real), 1, infile);
  fret = fread(&rec_prec_flow_field_ttime_out, sizeof(real), 1, infile);
  // close file
  fclose(infile);
}

void count_mem(void)
{
  cpumem = cpumem * 1e-6;
  gpumem = gpumem * 1e-6;
  real gpureal = 1.9 * gpumem + 80;
  printf("\nTotal CPU memory usage... %ldMB\n",cpumem);
  printf("Maximum estimated GPU memory usage... %1.0fMB\n",gpureal);
}

