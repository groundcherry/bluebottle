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

extern "C"
{
#include "bluebottle.h"
}

#include <cuda.h>

#include "cuda_quadrature.h"

extern "C"
void cuda_U_star_test_exp(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("\nIntermediate velocity calculation validation:\n\n");
  printf("  u = exp(x), v = exp(y), w = exp(z)\n\n");

  // set up expected solution
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_a[C] = nu;
        u_a[C] -= 2 * exp((i-1.5)*Dom.dx);
        u_a[C] -= exp((j-1.0)*Dom.dy);
        u_a[C] -= exp((k-1.0)*Dom.dz);
        u_a[C] *= dt * exp((i-1.5)*Dom.dx);
        u_a[C] += exp((i-1.5)*Dom.dx);
        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_a[C] = nu;
        v_a[C] -= 2 * exp((j-1.5)*Dom.dy);
        v_a[C] -= exp((k-1.0)*Dom.dz);
        v_a[C] -= exp((i-1.0)*Dom.dx);
        v_a[C] *= dt * exp((j-1.5)*Dom.dy);
        v_a[C] += exp((j-1.5)*Dom.dy);
        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_a[C] = nu;
        w_a[C] -= 2 * exp((k-1.5)*Dom.dz);
        w_a[C] -= exp((i-1.0)*Dom.dx);
        w_a[C] -= exp((j-1.0)*Dom.dy);
        w_a[C] *= dt * exp((k-1.5)*Dom.dz);
        w_a[C] += exp((k-1.5)*Dom.dz);
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // initialize input velocity fields
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = exp((i-1.5)*Dom.dx);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = exp((j-1.5)*Dom.dy);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = exp((k-1.5)*Dom.dz);
      }
    }
  }

  // write initial fields
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
#ifdef EXPLICIT
  printf("  Running cuda_U_star_2()...");
  cuda_U_star_2();
  printf("done.\n");
#endif

  // pull fields back to host
  printf("  Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u[C];
        u_e[C] = (u_c[C] - u_a[C]) / u_a[C];
        if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
        u[C] = u_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v[C];
        v_e[C] = (v_c[C] - v_a[C]) / v_a[C];
        if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w[C];
        w_e[C] = (w_c[C] - w_a[C]) / w_a[C];
        if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Velocity component:     minimum error:     maximum error:\n");
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}

extern "C"
void cuda_U_star_test_cos(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("Intermediate velocity calculation validation:\n\n");
  printf("  u = cos(y), v = cos(z), w = cos(x)\n\n");

dt = 1;
//dt0 = dt;
printf("dt = %f, dt0 = %f\n", dt, dt0);
  // set up expected solution
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
/*
        u_a[C] = nu * sin((j-1.0)*Dom.dy);
        u_a[C] += cos((j-1.0)*Dom.dy) * sin((k-1.0)*Dom.dz);
        u_a[C] *= -dt;
        u_a[C] += sin((j-1.0)*Dom.dy);
*/
        real x = 2.*PI*(i-1.0)*Dom.dx;
        real y = 2.*PI*(j-0.5)*Dom.dy;
        //real z = 2.*PI*(k-0.5)*Dom.dz;
        //u_a[C] = -4.*PI*sin(x)*cos(x)*sin(y)*sin(y);
        //u_a[C] += 2.*PI*sin(x)*cos(x)*(sin(y)*sin(y)-cos(y)*cos(y));
        u_a[C] = 8.*PI*PI*nu*cos(x)*sin(y);
        //u_a[C] += PI*sin(2.*x);
        u_a[C] *= -dt;
        //u_a[C] = cos(x)*sin(y) * exp(-16.*PI*PI*1.0*dt);
        u[C] = u_a[C];
        conv0_u[C] = -4.*PI*sin(x)*cos(x)*sin(y)*sin(y)
          + 2.*PI*sin(x)*cos(x)
          *(sin(y)*sin(y)-cos(y)*cos(y));
        conv0_u[C] *= exp(16.*PI*PI*1.0*dt);
#ifdef EXPLICIT
        diff0_u[C] = -8.*PI*PI*nu*cos(x)*sin(y);
        diff0_u[C] *= exp(16.*PI*PI*1.0*dt);
#endif
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
/*
        v_a[C] = nu * sin((k-1.0)*Dom.dz);
        v_a[C] += cos((k-1.0)*Dom.dz) * sin((i-1.0)*Dom.dx);
        v_a[C] *= -dt;
        v_a[C] += sin((k-1.0)*Dom.dz);
*/
        real x = 2.*PI*(i-0.5)*Dom.dx;
        real y = 2.*PI*(j-1.0)*Dom.dy;
        //real z = 2.*PI*(k-0.5)*Dom.dz;
        //v_a[C] = 2.*PI*cos(y)*sin(y)*(sin(x)*sin(x)-cos(x)*cos(x));
        //v_a[C] += -4.*PI*sin(x)*sin(x)*cos(y)*sin(y);
        v_a[C] = -8.*PI*PI*nu*sin(x)*cos(y);
//        v_a[C] += PI*sin(2.*y);
        v_a[C] *= -dt;
        //v_a[C] = -sin(x)*cos(y) * exp(-16.*PI*PI*1.0*dt);
        v[C] = v_a[C];
        conv0_v[C] = -4.*PI*sin(x)*sin(x)*sin(y)*cos(y)
          + 2.*PI*sin(y)*cos(y)
          *(sin(x)*sin(x)-cos(x)*cos(x));
        conv0_v[C] *= exp(16.*PI*PI*1.0*dt);
#ifdef EXPLICIT
        diff0_v[C] = 8.*PI*PI*nu*sin(x)*cos(y);
        diff0_v[C] *= exp(16.*PI*PI*1.0*dt);
#endif
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_a[C] = nu * sin((i-1.0)*Dom.dx);
        w_a[C] += cos((i-1.0)*Dom.dx) * sin((j-1.0)*Dom.dy);
        w_a[C] *= -dt;
        w_a[C] += sin((i-1.0)*Dom.dx);
        w_a[C] = 0;
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // initialize input pressure and velocity fields
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        real x = 2.*PI*(i-1.0)*Dom.dx;
        real y = 2.*PI*(j-0.5)*Dom.dy;
        //real z = 2.*PI*(k-0.5)*Dom.dz;
        //u[C] = sin((j-1.0)*Dom.dy);
        u[C] = cos(x)*sin(y);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        real x = 2.*PI*(i-0.5)*Dom.dx;
        real y = 2.*PI*(j-1.0)*Dom.dy;
        //real z = 2.*PI*(k-0.5)*Dom.dz;
        //v[C] = sin((k-1.0)*Dom.dz);
        v[C] = -sin(x)*cos(y);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        //w[C] = sin((i-1.0)*Dom.dx);
        w[C] = 0;
      }
    }
  }
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        real x = 2.*PI*(i-0.5)*Dom.dx;
        real y = 2.*PI*(j-0.5)*Dom.dy;
        //real z = 2.*PI*(k-0.5)*Dom.dz;
        p[C] = -0.25*rho_f*(cos(2.*x)+cos(2.*y));
        p0[C] = -0.25*rho_f*(cos(2.*x)+cos(2.*y));
      }
    }
  }

  // write initial fields
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
#ifdef EXPLICIT
  printf("  Running cuda_U_star_2()...");
  cuda_U_star_2();
  printf("done.\n");
#endif

  // pull fields back to host
  printf("  Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
u_err_max = 0;
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u_star[C];
        u_e[C] = (u_c[C] - u_a[C]);// / u_a[C];
        u_err_max += fabs(u_e[C]);
        //if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
        u[C] = u_e[C];
      }
    }
  }
v_err_max = 0;
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v_star[C];
        v_e[C] = (v_c[C] - v_a[C]);// / v_a[C];
        v_err_max += fabs(v_e[C]);
        //if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w_star[C];
        w_e[C] = (w_c[C] - w_a[C]);// / w_a[C];
        if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Velocity component:     minimum error:     maximum error:\n");
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}

__global__ void memcpy_u_star_test(real *dst, real *src, dom_struct *dom)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
    dst[i + j * dom->Gfx._s1b + k * dom->Gfx._s2b] = src[i + j * dom->Gfx._s1b
      + k * dom->Gfx._s2b];
  }
}

__global__ void memcpy_v_star_test(real *dst, real *src, dom_struct *dom)
{
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
    dst[i + j * dom->Gfy._s1b + k * dom->Gfy._s2b] = src[i + j * dom->Gfy._s1b
      + k * dom->Gfy._s2b];
  }
}

__global__ void memcpy_w_star_test(real *dst, real *src, dom_struct *dom)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
    dst[i + j * dom->Gfz._s1b + k * dom->Gfz._s2b] = src[i + j * dom->Gfz._s1b
      + k * dom->Gfz._s2b];
  }
}

__global__ void PP_memcpy_p_test(real *dst, real *src, dom_struct *dom)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  for(int i = dom->Gcc._is-DOM_BUF; i < dom->Gcc._ie-DOM_BUF; i++) {
    dst[(i+DOM_BUF) + (j+DOM_BUF) * dom->Gfx._s1b
      + (k+DOM_BUF) * dom->Gcc._s2b] = src[i + j * dom->Gcc._s1b
      + k * dom->Gcc._s2b];
  }
}

extern "C"
void cuda_BC_test(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *p_p_i = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_p_o = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_p_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_d_i = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_d_o = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_d_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_n_i = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_n_o = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_n_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gcc.s3b * sizeof(real);

  real *u_p_i = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_p_o = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_p_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_d_i = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_d_o = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_d_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_n_i = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_n_o = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_n_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gfx.s3b * sizeof(real);

  real *v_p_i = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_p_o = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_p_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_d_i = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_d_o = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_d_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_n_i = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_n_o = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_n_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gfy.s3b * sizeof(real);

  real *w_p_i = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_p_o = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_p_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_d_i = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_d_o = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_d_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_n_i = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_n_o = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_n_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real p_p_err_min = FLT_MAX;
  real p_p_err_max = FLT_MIN;
  real p_d_err_min = FLT_MAX;
  real p_d_err_max = FLT_MIN;
  real p_n_err_min = FLT_MAX;
  real p_n_err_max = FLT_MIN;
  real u_p_err_min = FLT_MAX;
  real u_p_err_max = FLT_MIN;
  real u_d_err_min = FLT_MAX;
  real u_d_err_max = FLT_MIN;
  real u_n_err_min = FLT_MAX;
  real u_n_err_max = FLT_MIN;
  real v_p_err_min = FLT_MAX;
  real v_p_err_max = FLT_MIN;
  real v_d_err_min = FLT_MAX;
  real v_d_err_max = FLT_MIN;
  real v_n_err_min = FLT_MAX;
  real v_n_err_max = FLT_MIN;
  real w_p_err_min = FLT_MAX;
  real w_p_err_max = FLT_MIN;
  real w_d_err_min = FLT_MAX;
  real w_d_err_max = FLT_MIN;
  real w_n_err_min = FLT_MAX;
  real w_n_err_max = FLT_MIN;

  printf("\nBoundary condition application validation:\n");

  // periodic field (on -1 <= x <= 1, -1 <= y <= 1, -1 <= z <= 1)
  printf("\n  Periodic boundary conditions:\n");
  printf("    p = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    u = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    v = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    w = cos(pi*x) + cos(pi*y) + cos(pi*z)\n\n");

  // write input fields
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_p_i[C] = cos(PI * (i-0.5)*Dom.dx);
        p_p_i[C] += cos(PI * (j-0.5)*Dom.dy);
        p_p_i[C] += cos(PI * (k-0.5)*Dom.dz);
        p[C] = p_p_i[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
      C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
      u_p_i[C] = cos(PI * (i-1.0)*Dom.dx);
      u_p_i[C] += cos(PI * (j-0.5)*Dom.dy);
      u_p_i[C] += cos(PI * (k-0.5)*Dom.dz);
      u[C] = u_p_i[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
      C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
      v_p_i[C] = cos(PI * (i-0.5)*Dom.dx);
      v_p_i[C] += cos(PI * (j-1.0)*Dom.dy);
      v_p_i[C] += cos(PI * (k-0.5)*Dom.dz);
      v[C] = v_p_i[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
      C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
      w_p_i[C] = cos(PI * (i-0.5)*Dom.dx);
      w_p_i[C] += cos(PI * (j-0.5)*Dom.dy);
      w_p_i[C] += cos(PI * (k-1.0)*Dom.dz);
      w[C] = w_p_i[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("    Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push to device
  printf("\n    Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // set BC to PERIODIC (overwrite input file for testing)
  printf("    Overwriting input file boundary conditions to ensure ");
  printf("PERIODIC...");
  bc.pW = PERIODIC;
  bc.pE = PERIODIC;
  bc.pS = PERIODIC;
  bc.pN = PERIODIC;
  bc.pB = PERIODIC;
  bc.pT = PERIODIC;
  bc.uW = PERIODIC;
  bc.uE = PERIODIC;
  bc.uS = PERIODIC;
  bc.uN = PERIODIC;
  bc.uB = PERIODIC;
  bc.uT = PERIODIC;
  bc.vW = PERIODIC;
  bc.vE = PERIODIC;
  bc.vS = PERIODIC;
  bc.vN = PERIODIC;
  bc.vB = PERIODIC;
  bc.vT = PERIODIC;
  bc.wW = PERIODIC;
  bc.wE = PERIODIC;
  bc.wS = PERIODIC;
  bc.wN = PERIODIC;
  bc.wB = PERIODIC;
  bc.wT = PERIODIC;
  printf("done.\n");

  // apply BC
  printf("    Running cuda_BC()...");
  cuda_dom_BC();
  printf("done.\n");

  // pull to host
  printf("    Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n    Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_p_o[C] = p[C];
        p_p_e[C] = p_p_o[C] - p_p_i[C];
        if(fabs(p_p_e[C]) > p_p_err_max) p_p_err_max = fabs(p_p_e[C]);
        if(fabs(p_p_e[C]) < p_p_err_min) p_p_err_min = fabs(p_p_e[C]);
        p[C] = p_p_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_p_o[C] = u[C];
        u_p_e[C] = u_p_o[C] - u_p_i[C];
        if(fabs(u_p_e[C]) > u_p_err_max) u_p_err_max = fabs(u_p_e[C]);
        if(fabs(u_p_e[C]) < u_p_err_min) u_p_err_min = fabs(u_p_e[C]);
        u[C] = u_p_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_p_o[C] = v[C];
        v_p_e[C] = v_p_o[C] - v_p_i[C];
        if(fabs(v_p_e[C]) > v_p_err_max) v_p_err_max = fabs(v_p_e[C]);
        if(fabs(v_p_e[C]) < v_p_err_min) v_p_err_min = fabs(v_p_e[C]);
        v[C] = v_p_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_p_o[C] = w[C];
        w_p_e[C] = w_p_o[C] - w_p_i[C];
        if(fabs(w_p_e[C]) > w_p_err_max) w_p_err_max = fabs(w_p_e[C]);
        if(fabs(w_p_e[C]) < w_p_err_min) w_p_err_min = fabs(w_p_e[C]);
        w[C] = w_p_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("    Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n    Error summary:\n");
  printf("    Field component:      minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_p_err_min, p_p_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_p_err_min, u_p_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_p_err_min, v_p_err_max);
  printf("          w              %12.3e       %12.3e\n",
    w_p_err_min, w_p_err_max);

  // Dirichlet field (on -1 <= x <= 1, -1 <= y <= 1, -1 <= z <= 1)
  printf("\n  Dirichlet boundary conditions:\n");
  printf("    p = sin(pi*x) * sin(pi*y) * sin(pi*z)\n");
  printf("    u = sin(pi*x) * sin(pi*y) * sin(pi*z)\n");
  printf("    v = sin(pi*x) * sin(pi*y) * sin(pi*z)\n");
  printf("    w = sin(pi*x) * sin(pi*y) * sin(pi*z)\n\n");

  // write input field
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_d_i[C] = sin(PI * (i-0.5)*Dom.dx);
        p_d_i[C] *= sin(PI * (j-0.5)*Dom.dy);
        p_d_i[C] *= sin(PI * (k-0.5)*Dom.dz);
        p[C] = p_d_i[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_d_i[C] = sin(PI * (i-1.0)*Dom.dx);
        u_d_i[C] *= sin(PI * (j-0.5)*Dom.dy);
        u_d_i[C] *= sin(PI * (k-0.5)*Dom.dz);
        u[C] = u_d_i[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_d_i[C] = sin(PI * (i-0.5)*Dom.dx);
        v_d_i[C] *= sin(PI * (j-1.0)*Dom.dy);
        v_d_i[C] *= sin(PI * (k-0.5)*Dom.dz);
        v[C] = v_d_i[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_d_i[C] = sin(PI * (i-0.5)*Dom.dx);
        w_d_i[C] *= sin(PI * (j-0.5)*Dom.dy);
        w_d_i[C] *= sin(PI * (k-1.0)*Dom.dz);
        w[C] = w_d_i[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("    Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push to device
  printf("\n    Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // set BC to DIRICHLET (overwrite input file for testing)
  printf("    Overwriting input file boundary conditions to ensure ");
  printf("DIRICHLET...");
  bc.pW = DIRICHLET;
  bc.pWD = 0.;
  bc.pE = DIRICHLET;
  bc.pED = 0.;
  bc.pS = DIRICHLET;
  bc.pSD = 0.;
  bc.pN = DIRICHLET;
  bc.pND = 0.;
  bc.pB = DIRICHLET;
  bc.pBD = 0.;
  bc.pT = DIRICHLET;
  bc.pTD = 0.;

  bc.uW = DIRICHLET;
  bc.uWD = 0.;
  bc.uE = DIRICHLET;
  bc.uED = 0.;
  bc.uS = DIRICHLET;
  bc.uSD = 0.;
  bc.uN = DIRICHLET;
  bc.uND = 0.;
  bc.uB = DIRICHLET;
  bc.uBD = 0.;
  bc.uT = DIRICHLET;
  bc.uTD = 0.;

  bc.vW = DIRICHLET;
  bc.vWD = 0.;
  bc.vE = DIRICHLET;
  bc.vED = 0.;
  bc.vS = DIRICHLET;
  bc.vSD = 0.;
  bc.vN = DIRICHLET;
  bc.vND = 0.;
  bc.vB = DIRICHLET;
  bc.vBD = 0.;
  bc.vT = DIRICHLET;
  bc.vTD = 0.;

  bc.wW = DIRICHLET;
  bc.wWD = 0.;
  bc.wE = DIRICHLET;
  bc.wED = 0.;
  bc.wS = DIRICHLET;
  bc.wSD = 0.;
  bc.wN = DIRICHLET;
  bc.wND = 0.;
  bc.wB = DIRICHLET;
  bc.wBD = 0.;
  bc.wT = DIRICHLET;
  bc.wTD = 0.;
  printf("done.\n");

  // apply BC
  printf("    Running cuda_BC()...");
  cuda_dom_BC();
  printf("done.\n");

  // pull to host
  printf("    Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n    Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_d_o[C] = p[C];
        p_d_e[C] = p_d_o[C] - p_d_i[C];
        if(fabs(p_d_e[C]) > p_d_err_max) p_d_err_max = fabs(p_d_e[C]);
        if(fabs(p_d_e[C]) < p_d_err_min) p_d_err_min = fabs(p_d_e[C]);
        p[C] = p_d_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_d_o[C] = u[C];
        u_d_e[C] = u_d_o[C] - u_d_i[C];
        if(fabs(u_d_e[C]) > u_d_err_max) u_d_err_max = fabs(u_d_e[C]);
        if(fabs(u_d_e[C]) < u_d_err_min) u_d_err_min = fabs(u_d_e[C]);
        u[C] = u_d_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_d_o[C] = v[C];
        v_d_e[C] = v_d_o[C] - v_d_i[C];
        if(fabs(v_d_e[C]) > v_d_err_max) v_d_err_max = fabs(v_d_e[C]);
        if(fabs(v_d_e[C]) < v_d_err_min) v_d_err_min = fabs(v_d_e[C]);
        v[C] = v_d_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_d_o[C] = w[C];
        w_d_e[C] = w_d_o[C] - w_d_i[C];
        if(fabs(w_d_e[C]) > w_d_err_max) w_d_err_max = fabs(w_d_e[C]);
        if(fabs(w_d_e[C]) < w_d_err_min) w_d_err_min = fabs(w_d_e[C]);
        w[C] = w_d_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("    Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n    Error summary:\n");
  printf("    Field component:      minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_d_err_min, p_d_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_d_err_min, u_d_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_d_err_min, v_d_err_max);
  printf("          w              %12.3e       %12.3e\n",
    w_d_err_min, w_d_err_max);

  // Neumann field (on -1 <= x <= 1, -1 <= y <= 1, -1 <= z <= 1)
  printf("\n  Neumann boundary conditions:\n");
  printf("    p = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    u = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    v = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    w = cos(pi*x) + cos(pi*y) + cos(pi*z)\n\n");

  // write input field
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_n_i[C] = cos(PI * (i-0.5)*Dom.dx);
        p_n_i[C] += cos(PI * (j-0.5)*Dom.dy);
        p_n_i[C] += cos(PI * (k-0.5)*Dom.dz);
        p[C] = p_n_i[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_n_i[C] = cos(PI * (i-1.0)*Dom.dx);
        u_n_i[C] += cos(PI * (j-0.5)*Dom.dy);
        u_n_i[C] += cos(PI * (k-0.5)*Dom.dz);
        u[C] = u_n_i[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_n_i[C] = cos(PI * (i-0.5)*Dom.dx);
        v_n_i[C] += cos(PI * (j-1.0)*Dom.dy);
        v_n_i[C] += cos(PI * (k-0.5)*Dom.dz);
        v[C] = v_n_i[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_n_i[C] = cos(PI * (i-0.5)*Dom.dx);
        w_n_i[C] += cos(PI * (j-0.5)*Dom.dy);
        w_n_i[C] += cos(PI * (k-1.0)*Dom.dz);
        w[C] = w_n_i[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("    Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push to device
  printf("\n    Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // set BC to NEUMANN (overwrite input file for testing)
  printf("    Overwriting input file boundary conditions to ensure ");
  printf("NEUMANN...");
  bc.pW = NEUMANN;
  bc.pE = NEUMANN;
  bc.pS = NEUMANN;
  bc.pN = NEUMANN;
  bc.pB = NEUMANN;
  bc.pT = NEUMANN;
  bc.uW = NEUMANN;
  bc.uE = NEUMANN;
  bc.uS = NEUMANN;
  bc.uN = NEUMANN;
  bc.uB = NEUMANN;
  bc.uT = NEUMANN;
  bc.vW = NEUMANN;
  bc.vE = NEUMANN;
  bc.vS = NEUMANN;
  bc.vN = NEUMANN;
  bc.vB = NEUMANN;
  bc.vT = NEUMANN;
  bc.wW = NEUMANN;
  bc.wE = NEUMANN;
  bc.wS = NEUMANN;
  bc.wN = NEUMANN;
  bc.wB = NEUMANN;
  bc.wT = NEUMANN;
  printf("done.\n");

  // apply BC
  printf("    Running cuda_BC()...");
  cuda_dom_BC();
  printf("done.\n");

  // pull to host
  printf("    Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n    Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_n_o[C] = p[C];
        p_n_e[C] = p_n_o[C] - p_n_i[C];
        if(fabs(p_n_e[C]) > p_n_err_max) p_n_err_max = fabs(p_n_e[C]);
        if(fabs(p_n_e[C]) < p_n_err_min) p_n_err_min = fabs(p_n_e[C]);
        p[C] = p_n_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_n_o[C] = u[C];
        u_n_e[C] = u_n_o[C] - u_n_i[C];
        if(fabs(u_n_e[C]) > u_n_err_max) u_n_err_max = fabs(u_n_e[C]);
        if(fabs(u_n_e[C]) < u_n_err_min) u_n_err_min = fabs(u_n_e[C]);
        u[C] = u_n_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_n_o[C] = v[C];
        v_n_e[C] = v_n_o[C] - v_n_i[C];
        if(fabs(v_n_e[C]) > v_n_err_max) v_n_err_max = fabs(v_n_e[C]);
        if(fabs(v_n_e[C]) < v_n_err_min) v_n_err_min = fabs(v_n_e[C]);
        v[C] = v_n_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_n_o[C] = w[C];
        w_n_e[C] = w_n_o[C] - w_n_i[C];
        if(fabs(w_n_e[C]) > w_n_err_max) w_n_err_max = fabs(w_n_e[C]);
        if(fabs(w_n_e[C]) < w_n_err_min) w_n_err_min = fabs(w_n_e[C]);
        w[C] = w_n_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("    Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n    Error summary:\n");
  printf("    Field component:      minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_n_err_min, p_n_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_n_err_min, u_n_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_n_err_min, v_n_err_max);
  printf("          w              %12.3e       %12.3e\n",
    w_n_err_min, w_n_err_max);

  // clean up
  free(p_p_i);
  free(p_p_o);
  free(p_p_e);
  free(p_d_i);
  free(p_d_o);
  free(p_d_e);
  free(p_n_i);
  free(p_n_o);
  free(p_n_e);
  free(u_p_i);
  free(u_p_o);
  free(u_p_e);
  free(u_d_i);
  free(u_d_o);
  free(u_d_e);
  free(u_n_i);
  free(u_n_o);
  free(u_n_e);
  free(v_p_i);
  free(v_p_o);
  free(v_p_e);
  free(v_d_i);
  free(v_d_o);
  free(v_d_e);
  free(v_n_i);
  free(v_n_o);
  free(v_n_e);
  free(w_p_i);
  free(w_p_o);
  free(w_p_e);
  free(w_d_i);
  free(w_d_o);
  free(w_d_e);
  free(w_n_i);
  free(w_n_o);
  free(w_n_e);
}

void cuda_project_test(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gyz.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("\nVelocity projection calculation validation:\n\n");
  printf("  u = exp(x), v = exp(y), w = exp(z), ");
  printf("p = exp(x) + exp(y) + exp(z)\n\n");

  // set up expected solution
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_a[C] = nu;
        u_a[C] -= 1. / rho_f;
        u_a[C] -= 2 * exp((i-1.5)*Dom.dx);
        u_a[C] -= exp((j-1.0)*Dom.dy);
        u_a[C] -= exp((k-1.0)*Dom.dz);
        u_a[C] *= dt * exp((i-1.5)*Dom.dx);
        u_a[C] += exp((i-1.5)*Dom.dx);
        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_a[C] = nu;
        v_a[C] -= 1. / rho_f;
        v_a[C] -= 2 * exp((j-1.5)*Dom.dy);
        v_a[C] -= exp((k-1.0)*Dom.dz);
        v_a[C] -= exp((i-1.0)*Dom.dx);
        v_a[C] *= dt * exp((j-1.5)*Dom.dy);
        v_a[C] += exp((j-1.5)*Dom.dy);
        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_a[C] = nu;
        w_a[C] -= 1. / rho_f;
        w_a[C] -= 2 * exp((k-1.5)*Dom.dz);
        w_a[C] -= exp((i-1.0)*Dom.dx);
        w_a[C] -= exp((j-1.0)*Dom.dy);
        w_a[C] *= dt * exp((k-1.5)*Dom.dz);
        w_a[C] += exp((k-1.5)*Dom.dz);
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // initialize input pressure and velocity fields
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p[C] = exp((i-1.0)*Dom.dx) + exp((j-1.0)*Dom.dy) + exp((k-1.0)*Dom.dz);
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = exp((i-1.5)*Dom.dx);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = exp((j-1.5)*Dom.dy);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = exp((k-1.5)*Dom.dz);
      }
    }
  }

  // write initial fields
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
#ifdef EXPLICIT
  printf("  Running cuda_U_star_2()...");
  cuda_U_star_2();
  cuda_project();
  printf("done.\n");
#endif

  // pull fields back to host
  printf("  Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u[C];
        u_e[C] = (u_c[C] - u_a[C]) / u_a[C];
        if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
        u[C] = u_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v[C];
        v_e[C] = (v_c[C] - v_a[C]) / v_a[C];
        if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w[C];
        w_e[C] = (w_c[C] - w_a[C]) / w_a[C];
        if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Velocity component:     minimum error:     maximum error:\n");
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}

extern "C"
void cuda_quad_interp_test(void)
{
  int i, j, k;  // iterators
  int C;
  real *p_a = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real x, y, z;

  // min and max error search
  real *p_err_min = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *p_err_max = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *u_err_min = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *u_err_max = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *v_err_min = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *v_err_max = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *w_err_min = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);
  real *w_err_max = (real*) malloc(nparts * sizeof(real));
  // cpumem += nparts * sizeof(real);

  (cudaSetDevice(dev_start));

  printf("\nLebedev quadrature interpolation validation:\n\n");
  printf("  p = u = v = w = exp(x) + exp(y) + exp(z)\n\n");

  // create analytic result and push to device
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.ksb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;
        y = (j-0.5)*Dom.dy + Dom.ys;
        z = (k-0.5)*Dom.dz + Dom.zs;

        p_a[C] = exp(x) + exp(y) + exp(z);
        p[C] = p_a[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.ksb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;

        x = (i-1.0)*Dom.dx + Dom.xs;
        y = (j-0.5)*Dom.dy + Dom.ys;
        z = (k-0.5)*Dom.dz + Dom.zs;

        u_a[C] = exp(x) + exp(y) + exp(z);
        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.ksb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;
        y = (j-1.0)*Dom.dy + Dom.ys;
        z = (k-0.5)*Dom.dz + Dom.zs;

        v_a[C] = exp(x) + exp(y) + exp(z);
        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.ksb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;
        y = (j-0.5)*Dom.dy + Dom.ys;
        z = (k-1.0)*Dom.dz + Dom.zs;

        w_a[C] = exp(x) + exp(y) + exp(z);
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
  printf("  Running cuda_quad_interp()...");
  // set up quadrature nodes for 7th-order Lebedev quadrature
  /*real PI14 = 0.25 * PI;
  real PI12 = 0.5 * PI;
  real PI34 = 0.75 * PI;
  real PI54 = 1.25 * PI;
  real PI32 = 1.5 * PI;
  real PI74 = 1.75 * PI;
  real alph1 = 0.955316618124509;
  real alph2 = 2.186276035465284;

  // nodes
  real a1_t[6] = {PI12, PI12, PI12, PI12, 0., 0.};
  real a1_p[6] = {0., PI12, PI, PI32, 0., PI};
  real a2_t[12] = {PI12, PI12, PI12, PI12,
                   PI14, PI14, PI14, PI14,
                   PI34, PI34, PI34, PI34};
  real a2_p[12] = {PI14, PI34, PI54, PI74,
                   0., PI12, PI, PI32,
                   0., PI12, PI, PI32};
  real a3_t[8] = {alph1, alph1, alph1, alph1,
                  alph2, alph2, alph2, alph2};
  real a3_p[8] = {PI14, PI34, PI54, PI74,
                  PI14, PI34, PI54, PI74};
  int nnodes = 26;

  // put all quadrature nodes together for interpolation
  real node_t[26];
  real node_p[26];
  for(i = 0; i < 6; i++) {
    node_t[i] = a1_t[i];
    node_p[i] = a1_p[i];
  }
  for(i = 0; i < 12; i++) {
    node_t[6+i] = a2_t[i];
    node_p[6+i] = a2_p[i];
  }
  for(i = 0; i < 8; i++) {
    node_t[18+i] = a3_t[i];
    node_p[18+i] = a3_p[i];
  }
*/


  real PI14 = 0.25 * PI;
  real PI12 = 0.5 * PI;
  real PI34 = 0.75 * PI;
  real PI54 = 1.25 * PI;
  real PI32 = 1.5 * PI;
  real PI74 = 1.75 * PI;
  real alph1 = 0.955316618124509; //54.736
  real alph2 = 2.186276035465284; //125.264
  /*real alph3 = 0.440510663004698; //25.239
  real alph4 = 2.701081990585095; //154.761
  real alph5 = 1.264518957625227; //72.452
  real alph6 = 1.877073695964566; //107.548
  real alph7 = 1.249045772398254; //71.565
  real alph8 = 1.892546881191539; //108.435
  real alph9 = 0.321750554396642; //18.435
  real alph10 = 2.819842099193151; //161.565
*/

  // nodes TODO: find a more elegant way of fixing the divide by sin(0)
  real a1_t[6] = {PI12, PI12, PI12, PI12, 0.+DIV_ST, PI-DIV_ST};
  real a1_p[6] = {0., PI12, PI, PI32, 0., 0.};
  real a2_t[12] = {PI12, PI12, PI12, PI12,
                   PI14, PI14, PI14, PI14,
                   PI34, PI34, PI34, PI34};
  real a2_p[12] = {PI14, PI34, PI54, PI74,
                   0., PI12, PI, PI32,
                   0., PI12, PI, PI32};
  real a3_t[8] = {alph1, alph1, alph1, alph1,
                  alph2, alph2, alph2, alph2};
  real a3_p[8] = {PI14, PI34, PI54, PI74,
                  PI14, PI34, PI54, PI74};
  /*real b_t[24] = {alph3, alph4, alph3, alph4,
                  alph3, alph4, alph3, alph4,
                  alph5, alph5, alph6, alph6,
                  alph5, alph5, alph6, alph6,
                  alph5, alph5, alph6, alph6,
                  alph5, alph5, alph6, alph6};
  real b_p[24] = {PI14, PI74, PI14, PI74,
                  PI34, PI54, PI34, PI54,
                  alph7, -alph7, alph7, -alph7,
                  alph8, -alph8, alph8, -alph8,
                  alph9, alph10, alph9, alph10,
                  -alph9, -alph10, -alph9, -alph10};
*/

  int nnodes = 26;

  // put all quadrature nodes together for interpolation
  real node_t[nnodes];
  real node_p[nnodes];
  for(i = 0; i < 6; i++) {
    node_t[i] = a1_t[i];
    node_p[i] = a1_p[i];
  }
  for(i = 0; i < 12; i++) {
    node_t[6+i] = a2_t[i];
    node_p[6+i] = a2_p[i];
  }
  for(i = 0; i < 8; i++) {
    node_t[18+i] = a3_t[i];
    node_p[18+i] = a3_p[i];
  }
  /*for(i = 0; i < 24; i++) {
    node_t[26+i] = b_t[i];
    node_p[26+i] = b_p[i];
  }
*/

  // create a place to temporarily store field variables at quadrature nodes
  real *_node_t;
  real *_node_p;
  (cudaMalloc((void**) &_node_t, nnodes * sizeof(real)));
  gpumem += nnodes * sizeof(real);
  (cudaMalloc((void**) &_node_p, nnodes * sizeof(real)));
  gpumem += nnodes * sizeof(real);
  (cudaMemcpy(_node_t, node_t, nnodes * sizeof(real),
    cudaMemcpyHostToDevice));
  (cudaMemcpy(_node_p, node_p, nnodes * sizeof(real),
    cudaMemcpyHostToDevice));
  real *_pp;
  real *_ur;
  real *_ut;
  real *_up;
  real *pp = (real*) malloc(nnodes * nparts * sizeof(real));
  // cpumem += nnodes * nparts * sizeof(real);
  real *ur = (real*) malloc(nnodes * nparts * sizeof(real));
  // cpumem += nnodes * nparts * sizeof(real);
  real *ut = (real*) malloc(nnodes * nparts * sizeof(real));
  // cpumem += nnodes * nparts * sizeof(real);
  real *up = (real*) malloc(nnodes * nparts * sizeof(real));
  // cpumem += nnodes * nparts * sizeof(real);
  (cudaMalloc((void**) &_pp, nnodes * nparts * sizeof(real)));
  gpumem += nnodes * nparts * sizeof(real);
  (cudaMalloc((void**) &_ur, nnodes * nparts * sizeof(real)));
  gpumem += nnodes * nparts * sizeof(real);
  (cudaMalloc((void**) &_ut, nnodes * nparts * sizeof(real)));
  gpumem += nnodes * nparts * sizeof(real);
  (cudaMalloc((void**) &_up, nnodes * nparts * sizeof(real)));
  gpumem += nnodes * nparts * sizeof(real);

  cuda_quad_interp(dev_start, _node_t, _node_p, nnodes, _pp, _ur, _ut, _up);
  printf("done.\n");

  // pull fields back to host
  printf("  Pulling fields back to host...");
  (cudaMemcpy(pp, _pp, nnodes * nparts * sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaMemcpy(ur, _ur, nnodes * nparts * sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaMemcpy(ut, _ut, nnodes * nparts * sizeof(real),
    cudaMemcpyDeviceToHost));
  (cudaMemcpy(up, _up, nnodes * nparts * sizeof(real),
    cudaMemcpyDeviceToHost));
  printf("done.\n");

  for(i = 0; i < nnodes; i++) {
    printf("xx[%d] = %f yy[%d] = %f zz[%d] = %f\n", i, ur[i], i, ut[i], i, up[i]);
  }

  // write computed solution
  printf("\n  Writing summarized solution to: out_%d.interp...", rec_paraview_stepnum_out);
  char path[FILE_NAME_SIZE] = "";
  sprintf(path, "%s/output/out_%d.interp", ROOT_DIR, rec_paraview_stepnum_out);
  FILE *file = fopen(path, "w");
  if(file == NULL) {
    fprintf(stderr, "Could not open file out_%d.interp", rec_paraview_stepnum_out);
    exit(EXIT_FAILURE);
  }
  for(int part = 0; part < nparts; part++) {
    p_err_min[part] = FLT_MAX;
    p_err_max[part] = FLT_MIN;
    u_err_min[part] = FLT_MAX;
    u_err_max[part] = FLT_MIN;
    v_err_min[part] = FLT_MAX;
    v_err_max[part] = FLT_MIN;
    w_err_min[part] = FLT_MAX;
    w_err_max[part] = FLT_MIN;
    fprintf(file, "parts[%d].rs = %f\n", part, parts[part].rs);
    fprintf(file, "%11s%11s%11s%11s%11s%11s%11s\n",
      "theta", "phi", "expected", "p_err", "u_err", "v_err", "w_err");
    for(int n = 0; n < nnodes; n++) {
      real x_tmp = parts[part].rs*sin(node_t[n])*cos(node_p[n])
        + parts[part].x;
      real y_tmp = parts[part].rs*sin(node_t[n])*sin(node_p[n])
        + parts[part].y;
      real z_tmp = parts[part].rs*cos(node_t[n])
        + parts[part].z;
      real pa_tmp = exp(x_tmp) + exp(y_tmp) + exp(z_tmp);
      real u_tmp = ur[n+part*nnodes]*sin(node_t[n])*cos(node_p[n]);
      u_tmp += ut[n+part*nnodes]*cos(node_t[n])*cos(node_p[n]);
      u_tmp -= up[n+part*nnodes]*sin(node_p[n]);
      real v_tmp = ur[n+part*nnodes]*sin(node_t[n])*sin(node_p[n]);
      v_tmp += ut[n+part*nnodes]*cos(node_t[n])*sin(node_p[n]);
      v_tmp += up[n+part*nnodes]*cos(node_p[n]);
      real w_tmp = ur[n+part*nnodes]*cos(node_t[n]);
      w_tmp -= ut[n+part*nnodes]*sin(node_t[n]);

      real p_out = (pa_tmp-pp[n+part*nnodes]) / pa_tmp;
      real u_out = (pa_tmp-u_tmp) / pa_tmp;
      real v_out = (pa_tmp-v_tmp) / pa_tmp;
      real w_out = (pa_tmp-w_tmp) / pa_tmp;

      if(fabs(p_out) < p_err_min[part]) p_err_min[part] = fabs(p_out);
      if(fabs(p_out) > p_err_max[part]) p_err_max[part] = fabs(p_out);
      if(fabs(u_out) < u_err_min[part]) u_err_min[part] = fabs(u_out);
      if(fabs(u_out) > u_err_max[part]) u_err_max[part] = fabs(u_out);
      if(fabs(v_out) < v_err_min[part]) v_err_min[part] = fabs(v_out);
      if(fabs(v_out) > v_err_max[part]) v_err_max[part] = fabs(v_out);
      if(fabs(w_out) < w_err_min[part]) w_err_min[part] = fabs(w_out);
      if(fabs(w_out) > w_err_max[part]) w_err_max[part] = fabs(w_out);

      fprintf(file, "%11.7f%11.7f%11.7f%11.3e%11.3e%11.3e%11.3e\n",
        node_t[n], node_p[n], pa_tmp, p_out, u_out, v_out, w_out);
    }
  }
  fclose(file);
  printf("done.\n");

  printf("\n  Error summary:\n");
  for(int a = 0; a < nparts; a++) {
    printf("  Particle %d\n", a);
    printf("    Field component:     minimum error:     maximum error:\n");
    printf("          p              %12.3e       %12.3e\n",
      p_err_min[a], p_err_max[a]);
    printf("          u              %12.3e       %12.3e\n",
      u_err_min[a], u_err_max[a]);
    printf("          v              %12.3e       %12.3e\n",
      v_err_min[a], v_err_max[a]);
    printf("          w              %12.3e       %12.3e\n\n",
      w_err_min[a], w_err_max[a]);
  }
  free(p_a);
  free(u_a);
  free(v_a);
  free(w_a);
  free(p_err_min);
  free(p_err_max);
  free(u_err_min);
  free(u_err_max);
  free(v_err_min);
  free(v_err_max);
  free(w_err_min);
  free(w_err_max);
  free(pp);
  free(ur);
  free(ut);
  free(up);
  cudaFree(_node_t);
  cudaFree(_node_p);
  cudaFree(_pp);
  cudaFree(_ur);
  cudaFree(_ut);
  cudaFree(_up);
}

extern "C"
void cuda_lamb_test(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations
  real x, y, z;
  real r, theta, phi;
  real a = parts[0].r;

  real *p_a = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_c = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gcc.s3b * sizeof(real)
  real *p_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gcc.s3b * sizeof(real)
  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real)
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real)
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real)
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real)
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real)
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfy.s3b * sizeof(real)
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real)
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real)
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real)

  // min and max error search
  real p_err_min = FLT_MAX;
  real p_err_max = FLT_MIN;
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("\nLamb's coefficient calculation validation:\n\n");
  printf("  u = exp(x), v = exp(y), w = exp(z), ");
  printf("p = exp(x) + exp(y) + exp(z)\n\n");

  real U = 0.;
  real V = 0.;
  real W = 0.01;

  // set up expected solution
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-0.5)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-0.5)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        p_a[C] = -1.5*(U*st*cp + V*st*sp + W*ct)*a/r/r;

        p[C] = p_a[C];
      }
    }
  }
 
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;

        x = (i-1.0)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-0.5)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-0.5)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        u_a[C] = -0.75*a/r*(U + (U*st*cp + V*st*sp + W*ct)*st*cp)
          - 0.25*a*a*a/r/r/r*(U - 3.*(U*st*cp + V*st*sp + W*ct)*st*cp)
          + U;

        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-1.0)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-0.5)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        v_a[C] = -0.75*a/r*(V + (U*st*cp + V*st*sp + W*ct)*st*sp)
          - 0.25*a*a*a/r/r/r*(V - 3.*(U*st*cp + V*st*sp + W*ct)*st*sp)
          + V;

        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-0.5)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-1.0)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        w_a[C] = -0.75*a/r*(W + (U*st*cp + V*st*sp + W*ct)*ct)
          - 0.25*a*a*a/r/r/r*(W - 3.*(U*st*cp + V*st*sp + W*ct)*ct)
          + W;

        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK();
  printf("done.\n");

  // set up expected solution
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-0.5)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-0.5)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        p_a[C] = -1.5*(U*st*cp + V*st*sp + W*ct)*a/r/r;

        p[C] = p_a[C];
      }
    }
  }
 
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;

        x = (i-1.0)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-0.5)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-0.5)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        u_a[C] = -0.75*a/r*(U + (U*st*cp + V*st*sp + W*ct)*st*cp)
          - 0.25*a*a*a/r/r/r*(U - 3.*(U*st*cp + V*st*sp + W*ct)*st*cp)
          + U;

        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-1.0)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-0.5)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        v_a[C] = -0.75*a/r*(V + (U*st*cp + V*st*sp + W*ct)*st*sp)
          - 0.25*a*a*a/r/r/r*(V - 3.*(U*st*cp + V*st*sp + W*ct)*st*sp)
          + V;

        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;

        x = (i-0.5)*Dom.dx + Dom.xs;// - parts[0].x;
        y = (j-0.5)*Dom.dy + Dom.ys;// - parts[0].y;
        z = (k-1.0)*Dom.dz + Dom.zs;// - parts[0].z;

        r = sqrt(x*x+y*y+z*z);
        theta = acos(z/r);
        phi = acos(x/sqrt(x*x+y*y));
        if(y<0.) phi = 2.*PI-phi;

        real st = sin(theta);
        real ct = cos(theta);
        real sp = sin(phi);
        real cp = cos(phi);

        w_a[C] = -0.75*a/r*(W + (U*st*cp + V*st*sp + W*ct)*ct)
          - 0.25*a*a*a/r/r/r*(W - 3.*(U*st*cp + V*st*sp + W*ct)*ct)
          + W;

        w[C] = w_a[C];
      }
    }
  }

  // write initial fields (same as expected solution)
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
  printf("  Running cuda_part_BC()...");
  cuda_Lamb();
  cuda_part_BC();
  cuda_part_pull();
  char nam[FILE_NAME_SIZE] = "lamb.rec";
  recorder_lamb(nam,0);
  printf("done.\n");

  // pull fields back to host
  printf("  Pulling fields back to host...");
  //cuda_div_U();
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_c[C] = p[C];
        if(p_c[C] != 0)
          p_e[C] = (p_c[C] - p_a[C]) / p_c[C];
        else
          p_e[C] = (p_c[C] - p_a[C]);
        if(fabs(p_e[C]) > p_err_max) p_err_max = fabs(p_e[C]);
        if(fabs(p_e[C]) < p_err_min) p_err_min = fabs(p_e[C]);
        p[C] = p_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u[C];
        /*if(u_c[C] != 0)
          u_e[C] = (u_c[C] - u_a[C]) / u_c[C];
        else
*/
          u_e[C] = (u_c[C] - u_a[C]);
        /*if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
*/
        u[C] = u_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v[C];
        /*if(v_c[C] != 0)
          v_e[C] = (v_c[C] - v_a[C]) / v_c[C];
        else
*/
          v_e[C] = (v_c[C] - v_a[C]);
        /*if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
*/
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w[C];
        /*if(w_c[C] != 0)
          w_e[C] = (w_c[C] - w_a[C]) / w_c[C];
        else
*/
          w_e[C] = (w_c[C] - w_a[C]);
        /*if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
*/
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Field variable:     minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_err_min, p_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(p_a);
  free(p_c);
  free(p_e);
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}
