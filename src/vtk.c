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

#include "bluebottle.h"
#include "particle.h"

void init_VTK(void)
{
  char fname[FILE_NAME_SIZE] = "";

  // open PVD file for writing
  sprintf(fname, "%sout.pvd", OUTPUT_DIR);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write PVD file header and footer
  fprintf(outfile, "<VTKFile type=\"Collection\">\n");
  fprintf(outfile, "<Collection>\n");
  fprintf(outfile, "</Collection>\n");
  fprintf(outfile, "</VTKFile>");

  // close the file
  fclose(outfile);
}

void out_VTK(void)
{
  char fname_pvd[FILE_NAME_SIZE] = ""; // pvd filename
  char fname_pvtr[FILE_NAME_SIZE] = ""; // pvtr filename
  char fname_vtp[FILE_NAME_SIZE] = ""; // vtp filename
  char fnamenodes_vtp[FILE_NAME_SIZE] = ""; // vtp filename

  sprintf(fname_pvd, "%sout.pvd", OUTPUT_DIR);
  sprintf(fname_pvtr, "out_%d.pvtr", rec_paraview_stepnum_out);
  sprintf(fname_vtp, "out_%d.vtp", rec_paraview_stepnum_out);
  sprintf(fnamenodes_vtp, "out_nodes_%d.vtp", rec_paraview_stepnum_out);

  FILE *pvdfile = fopen(fname_pvd, "r+");
  if(pvdfile == NULL) {
    init_VTK();
    pvdfile = fopen(fname_pvd, "r+");
  }
  // moves back 2 lines from the end of the file (above the footer)
  fseek(pvdfile, -24, SEEK_END);

  fprintf(pvdfile, "<DataSet timestep=\"%e\" part=\"0\" file=\"%s\"/>\n",
    ttime, fname_pvtr);
  fprintf(pvdfile, "<DataSet timestep=\"%e\" part=\"1\" file=\"%s\"/>\n",
    ttime, fname_vtp);
  fprintf(pvdfile, "<DataSet timestep=\"%e\" part=\"2\" file=\"%s\"/>\n",
    ttime, fnamenodes_vtp);
  fprintf(pvdfile, "</Collection>\n");
  fprintf(pvdfile, "</VTKFile>");
  fclose(pvdfile);

  dom_out_VTK();
  part_out_VTK();
  quadnodes_out_VTK();
}

void dom_out_VTK(void)
{
  int i, j, k, l; // iterators
  char fname[FILE_NAME_SIZE] = ""; // output filename
  char fname_dom[FILE_NAME_SIZE] = ""; // subdomain filename
  int C;  // cell center index
  int Cx;  // cell center index for interpolation
  int Cy;  // cell center index for interpolation
  int Cz;  // cell center index for interpolation

  // number of cells in a subdomain (length, start, end)
  int ncx_l, ncx_s, ncx_e;
  int ncy_l, ncy_s, ncy_e;
  int ncz_l, ncz_s, ncz_e;

  sprintf(fname, "%sout_%d.pvtr", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write Paraview pvtr file
  fprintf(outfile, "<VTKFile type=\"PRectilinearGrid\">\n");
  fprintf(outfile, "<PRectilinearGrid WholeExtent=");
  fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
    Dom.xn, Dom.yn, Dom.zn);
  //fprintf(outfile, "<PCellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
  fprintf(outfile, "<PCellData Scalars=\"p phase\" Vectors=\"vel\">\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"p\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"divU\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase_shell\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"flag\"");
  //fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "</PCellData>\n");
  fprintf(outfile, "<PCoordinates>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"x\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"y\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"z\"/>\n");
  fprintf(outfile, "</PCoordinates>\n");
  for(l = 0; l < 6 * nsubdom; l += 6) {
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\" ", ncz_s, ncz_e);
    sprintf(fname_dom, "out_%d_%d.vtr", rec_paraview_stepnum_out, l/6);
    fprintf(outfile, "Source=\"%s\"/>\n", fname_dom);
  }
  fprintf(outfile, "</PRectilinearGrid>\n");
  fprintf(outfile, "</VTKFile>\n");
  fclose(outfile);

  // interpolate velocities to cell centers
  // cell-center working arrays
  real *uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        // interpolate velocity
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        uu[C] = 0.5 * (u[Cx] + u[Cx+1]);
        vv[C] = 0.5 * (v[Cy] + v[Cy+Dom.Gfy.s1b]);
        ww[C] = 0.5 * (w[Cz] + w[Cz+Dom.Gfz.s2b]);
        // interpolate flags
        //flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        //flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        //flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }

  // write each subdomain file
  for(l = 0; l < 6 * nsubdom; l += 6) {
    // number of cells in the subdomain
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;

    // open file for writing
    sprintf(fname, "%s/out_%d_%d.vtr", OUTPUT_DIR, rec_paraview_stepnum_out, l/6);
    FILE *outfile = fopen(fname, "w");
    if(outfile == NULL) {
      fprintf(stderr, "Could not open file %s\n", fname);
      exit(EXIT_FAILURE);
    }

    fprintf(outfile, "<VTKFile type=\"RectilinearGrid\">\n");
    fprintf(outfile, "<RectilinearGrid WholeExtent=");
    fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
      Dom.xn, Dom.yn, Dom.zn);
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\">\n", ncz_s, ncz_e);
    //fprintf(outfile, "<CellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
    fprintf(outfile, "<CellData Scalars=\"p phase\" Vectors=\"vel\">\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"p\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", p[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    /*fprintf(outfile, "<DataArray type=\"Float32\" Name=\"divU\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", divU[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    /*fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase_shell\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase_shell[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/

    // write velocity vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu[C], vv[C], ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    // write flag vector
    /*fprintf(outfile, "<DataArray type=\"Float32\" Name=\"flag\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
          fprintf(outfile, "%lf %lf %lf ", flag_uu[C], flag_vv[C], flag_ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
*/
    fprintf(outfile, "</CellData>\n");

    fprintf(outfile, "<Coordinates>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"x\">\n");
    for(i = 0; i <= ncx_l; i++) {
      fprintf(outfile, "%lf ", dom[l].xs + Dom.dx * i);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"y\">\n");
    for(j = 0; j <= ncy_l; j++) {
      fprintf(outfile, "%lf ", dom[l].ys + Dom.dy * j);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"z\">\n");
    for(k = 0; k <= ncz_l; k++) {
      fprintf(outfile, "%lf ", dom[l].zs + Dom.dz * k);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "</Coordinates>\n");
    fprintf(outfile, "</Piece>\n");
    fprintf(outfile, "</RectilinearGrid>\n");
    fprintf(outfile, "</VTKFile>\n");
    fclose(outfile);
  }

  // clean up interpolated fields
  free(uu);
  free(vv);
  free(ww);
  //free(flag_uu);
  //free(flag_vv);
  //free(flag_ww);
}

void init_VTK_turb(void)
{
  char fname[FILE_NAME_SIZE] = "";

  // open PVD file for writing
  sprintf(fname, "%sout_turb.pvd", OUTPUT_DIR);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write PVD file header and footer
  fprintf(outfile, "<VTKFile type=\"Collection\">\n");
  fprintf(outfile, "<Collection>\n");
  fprintf(outfile, "</Collection>\n");
  fprintf(outfile, "</VTKFile>");

  // close the file
  fclose(outfile);
}

void out_VTK_turb(void)
{
  char fname_pvd[FILE_NAME_SIZE] = ""; // pvd filename
  char fname_pvtr[FILE_NAME_SIZE] = ""; // pvtr filename

  sprintf(fname_pvd, "%sout_turb.pvd", OUTPUT_DIR);
  sprintf(fname_pvtr, "out_turb_%d.pvtr", rec_prec_stepnum_out);

  FILE *pvdfile= fopen(fname_pvd, "r+");
  if (pvdfile == NULL) {
    init_VTK_turb();
    pvdfile= fopen(fname_pvd, "r+");
  }
  // moves back 2 lines from the end of the file (above the footer)
  fseek(pvdfile, -24, SEEK_END);

  fprintf(pvdfile, "<DataSet timestep=\"%e\" part=\"0\" file=\"%s\"/>\n",
    ttime, fname_pvtr);
  fprintf(pvdfile, "</Collection>\n");
  fprintf(pvdfile, "</VTKFile>");
  fclose(pvdfile);

  dom_out_VTK_turb();
}

void dom_out_VTK_turb(void)
{
  int i, j, k, l; // iterators
  char fname[FILE_NAME_SIZE] = ""; // output filename
  char fname_dom[FILE_NAME_SIZE] = ""; // subdomain filename
  int C;  // cell center index
  int Cx;  // cell center index for interpolation
  int Cy;  // cell center index for interpolation
  int Cz;  // cell center index for interpolation

  // number of cells in a subdomain (length, start, end)
  int ncx_l, ncx_s, ncx_e;  
  int ncy_l, ncy_s, ncy_e;  
  int ncz_l, ncz_s, ncz_e;  

  sprintf(fname, "%sout_turb_%d.pvtr", OUTPUT_DIR, rec_prec_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write Paraview pvtr file
  fprintf(outfile, "<VTKFile type=\"PRectilinearGrid\">\n");
  fprintf(outfile, "<PRectilinearGrid WholeExtent=");
  fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
    Dom.xn, Dom.yn, Dom.zn);
  //fprintf(outfile, "<PCellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
  fprintf(outfile, "<PCellData Scalars=\"p\" Vectors=\"vel\">\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"p\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"divU\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase_shell\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"flag\"");
  //fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "</PCellData>\n");
  fprintf(outfile, "<PCoordinates>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"x\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"y\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"z\"/>\n");
  fprintf(outfile, "</PCoordinates>\n");
  for(l = 0; l < 6 * nsubdom; l += 6) {
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\" ", ncz_s, ncz_e);
    sprintf(fname_dom, "out_turb_%d_%d.vtr", rec_prec_stepnum_out, l/6);
    fprintf(outfile, "Source=\"%s\"/>\n", fname_dom);
  }
  fprintf(outfile, "</PRectilinearGrid>\n");
  fprintf(outfile, "</VTKFile>\n");
  fclose(outfile);

  // interpolate velocities to cell centers
  // cell-center working arrays
  real *uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        // interpolate velocity
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        uu[C] = 0.5 * (u[Cx] + u[Cx+1]);
        vv[C] = 0.5 * (v[Cy] + v[Cy+Dom.Gfy.s1b]);
        ww[C] = 0.5 * (w[Cz] + w[Cz+Dom.Gfz.s2b]);
        // interpolate flags
        //flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        //flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        //flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }

  // write each subdomain file
  for(l = 0; l < 6 * nsubdom; l += 6) {
    // number of cells in the subdomain
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;

    // open file for writing
    sprintf(fname, "%s/out_turb_%d_%d.vtr", OUTPUT_DIR, rec_prec_stepnum_out, l/6);
    FILE *outfile = fopen(fname, "w");
    if(outfile == NULL) {
      fprintf(stderr, "Could not open file %s\n", fname);
      exit(EXIT_FAILURE);
    }

    fprintf(outfile, "<VTKFile type=\"RectilinearGrid\">\n");
    fprintf(outfile, "<RectilinearGrid WholeExtent=");
    fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
      Dom.xn, Dom.yn, Dom.zn);
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\">\n", ncz_s, ncz_e);
    //fprintf(outfile, "<CellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
    fprintf(outfile, "<CellData Scalars=\"p\" Vectors=\"vel\">\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"p\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", p[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    //fprintf(outfile, "<DataArray type=\"Float32\" Name=\"divU\">\n");
    // write pressure for this subdomain
    /*for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", divU[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/
    /*fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase_shell\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase_shell[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/

    // write velocity vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu[C], vv[C], ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    // write flag vector
    /*fprintf(outfile, "<DataArray type=\"Float32\" Name=\"flag\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
          fprintf(outfile, "%lf %lf %lf ", flag_uu[C], flag_vv[C], flag_ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
*/
    fprintf(outfile, "</CellData>\n");

    fprintf(outfile, "<Coordinates>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"x\">\n");
    for(i = 0; i <= ncx_l; i++) {
      fprintf(outfile, "%lf ", dom[l].xs + Dom.dx * i);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"y\">\n");
    for(j = 0; j <= ncy_l; j++) {
      fprintf(outfile, "%lf ", dom[l].ys + Dom.dy * j);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"z\">\n");
    for(k = 0; k <= ncz_l; k++) {
      fprintf(outfile, "%lf ", dom[l].zs + Dom.dz * k);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "</Coordinates>\n");
    fprintf(outfile, "</Piece>\n");
    fprintf(outfile, "</RectilinearGrid>\n");
    fprintf(outfile, "</VTKFile>\n");
    fclose(outfile);
  }

  // clean up interpolated fields
  free(uu);
  free(vv);
  free(ww);
  //free(flag_uu);
  //free(flag_vv);
  //free(flag_ww);
}

void init_VTK_ghost(void)
{
  char fname[FILE_NAME_SIZE] = "";

  // open PVD file for writing
  sprintf(fname, "%sout_ghost.pvd", OUTPUT_DIR);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write PVD file header and footer
  fprintf(outfile, "<VTKFile type=\"Collection\">\n");
  fprintf(outfile, "<Collection>\n");
  fprintf(outfile, "</Collection>\n");
  fprintf(outfile, "</VTKFile>");

  // close the file
  fclose(outfile);
}

void out_VTK_ghost(void)
{
  char fname_pvd[FILE_NAME_SIZE] = ""; // pvd filename
  char fname_pvtr[FILE_NAME_SIZE] = ""; // pvtr filename
  char fname_vtp[FILE_NAME_SIZE] = ""; // vtp filename

  sprintf(fname_pvd, "%sout_ghost.pvd", OUTPUT_DIR);
  sprintf(fname_pvtr, "out_ghost_%d.pvtr", rec_paraview_stepnum_out);
  sprintf(fname_vtp, "out_%d.vtp", rec_paraview_stepnum_out);

  FILE *pvdfile= fopen(fname_pvd, "r+");
  if (pvdfile == NULL) {
    init_VTK_ghost();
    pvdfile= fopen(fname_pvd, "r+");
  }
  // moves back 2 lines from the end of the file (above the footer)
  fseek(pvdfile, -24, SEEK_END);

  fprintf(pvdfile, "<DataSet timestep=\"%e\" part=\"0\" file=\"%s\"/>\n",
    ttime, fname_pvtr);
  fprintf(pvdfile, "<DataSet timestep=\"%e\" part=\"1\" file=\"%s\"/>\n",
    ttime, fname_vtp);
  fprintf(pvdfile, "</Collection>\n");
  fprintf(pvdfile, "</VTKFile>");
  fclose(pvdfile);

  dom_out_VTK_ghost();
  part_out_VTK();
}

void dom_out_VTK_ghost(void)
{
  int i, j, k, l; // iterators
  char fname[FILE_NAME_SIZE] = ""; // output filename
  char fname_dom[FILE_NAME_SIZE] = ""; // subdomain filename
  int C;  // cell center index
  int Cx;  // cell center index for interpolation
  int Cy;  // cell center index for interpolation
  int Cz;  // cell center index for interpolation

  // number of cells in a subdomain (length, start, end)
  int ncx_l, ncx_s, ncx_e;  
  int ncy_l, ncy_s, ncy_e;  
  int ncz_l, ncz_s, ncz_e;  

  sprintf(fname, "%sout_ghost_%d.pvtr", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write Paraview pvtr file
  fprintf(outfile, "<VTKFile type=\"PRectilinearGrid\">\n");
  fprintf(outfile, "<PRectilinearGrid WholeExtent=");
  fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
    Dom.xn+2*DOM_BUF, Dom.yn+2*DOM_BUF, Dom.zn+2*DOM_BUF);
  //fprintf(outfile, "<PCellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
  fprintf(outfile, "<PCellData Scalars=\"p phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"p\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"divU\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase_shell\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel_star\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"flag\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "</PCellData>\n");
  fprintf(outfile, "<PCoordinates>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"x\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"y\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"z\"/>\n");
  fprintf(outfile, "</PCoordinates>\n");
  for(l = 0; l < 6 * nsubdom; l += 6) {
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx + 2*DOM_BUF;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy + 2*DOM_BUF;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz + 2*DOM_BUF;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\" ", ncz_s, ncz_e);
    sprintf(fname_dom, "out_ghost_%d_%d.vtr", rec_paraview_stepnum_out, l/6);
    fprintf(outfile, "Source=\"%s\"/>\n", fname_dom);
  }
  fprintf(outfile, "</PRectilinearGrid>\n");
  fprintf(outfile, "</VTKFile>\n");
  fclose(outfile);

  // interpolate velocities to cell centers
  // cell-center working arrays
  real *uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *uu_star = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv_star = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww_star = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *flag_uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *flag_vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *flag_ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        // interpolate velocity
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        uu[C] = 0.5 * (u[Cx] + u[Cx+1]);
        vv[C] = 0.5 * (v[Cy] + v[Cy+Dom.Gfy.s1b]);
        ww[C] = 0.5 * (w[Cz] + w[Cz+Dom.Gfz.s2b]);
        uu_star[C] = 0.5 * (u_star[Cx] + u_star[Cx+1]);
        vv_star[C] = 0.5 * (v_star[Cy] + v_star[Cy+Dom.Gfy.s1b]);
        ww_star[C] = 0.5 * (w_star[Cz] + w_star[Cz+Dom.Gfz.s2b]);
        // interpolate flags
        flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }

  // write each subdomain file
  for(l = 0; l < 6 * nsubdom; l += 6) {
    // number of cells in the subdomain
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx + 2*DOM_BUF;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy + 2*DOM_BUF;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz + 2*DOM_BUF;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;

    // open file for writing
    sprintf(fname, "%s/out_ghost_%d_%d.vtr", OUTPUT_DIR, rec_paraview_stepnum_out, l/6);
    FILE *outfile = fopen(fname, "w");
    if(outfile == NULL) {
      fprintf(stderr, "Could not open file %s\n", fname);
      exit(EXIT_FAILURE);
    }

    fprintf(outfile, "<VTKFile type=\"RectilinearGrid\">\n");
    fprintf(outfile, "<RectilinearGrid WholeExtent=");
    fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
      Dom.xn+2*DOM_BUF, Dom.yn+2*DOM_BUF, Dom.zn+2*DOM_BUF);
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\">\n", ncz_s, ncz_e);
    //fprintf(outfile, "<CellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
    fprintf(outfile, "<CellData Scalars=\"p phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"p\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", p[C]);
        }
      }
    }
    //fprintf(outfile, "\n");
    //fprintf(outfile, "</DataArray>\n");
    //fprintf(outfile, "<DataArray type=\"Float32\" Name=\"divU\">\n");
    // write pressure for this subdomain
    /*for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          if(i >= Dom.Gcc.is && i < Dom.Gcc.ie
            && j >= Dom.Gcc.js && j < Dom.Gcc.je
            && k >= Dom.Gcc.ks && k < Dom.Gcc.ke) {
            fprintf(outfile, "%1.16e ", divU[C]);
          } else {
            fprintf(outfile, "%1.16e ", 0.);
          }
        }
      }
    }
*/
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase_shell\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase_shell[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    // write velocity vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu[C], vv[C], ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    // write velocity star vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel_star\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu_star[C], vv_star[C], ww_star[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");

    // write flag vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"flag\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%lf %lf %lf ", flag_uu[C], flag_vv[C], flag_ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    fprintf(outfile, "</CellData>\n");

    fprintf(outfile, "<Coordinates>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"x\">\n");
    for(i = -1; i < ncx_l; i++) {
      fprintf(outfile, "%lf ", dom[l].xs + Dom.dx * i);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"y\">\n");
    for(j = -1; j < ncy_l; j++) {
      fprintf(outfile, "%lf ", dom[l].ys + Dom.dy * j);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"z\">\n");
    for(k = -1; k < ncz_l; k++) {
      fprintf(outfile, "%lf ", dom[l].zs + Dom.dz * k);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "</Coordinates>\n");
    fprintf(outfile, "</Piece>\n");
    fprintf(outfile, "</RectilinearGrid>\n");
    fprintf(outfile, "</VTKFile>\n");
    fclose(outfile);
  }

  // clean up interpolated fields
  free(uu);
  free(vv);
  free(ww);
  free(uu_star);
  free(vv_star);
  free(ww_star);
  free(flag_uu);
  free(flag_vv);
  free(flag_ww);
}

void part_out_VTK(void)
{
  int i; // iterator
  char fname[FILE_NAME_SIZE] = ""; // output filename
  int nparts_plot = 0;    // the number of particles to plot; may be greater
                          // may be greater than nparts if particles straddle

  // open file for writing
  sprintf(fname, "%s/out_%d.vtp", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // go through all particles and check if any straddle boundaries
  nparts_plot = nparts;
  for(i = 0; i < nparts; i++) {
    int straddle_x = 0;
    int straddle_y = 0;
    int straddle_z = 0;
    // figure out which boudaries this particle straddles
    if(parts[i].x < (Dom.xs + parts[i].r)) straddle_x = 1;
    if(parts[i].x > (Dom.xe - parts[i].r)) straddle_x = 1;
    if(parts[i].y < (Dom.ys + parts[i].r)) straddle_y = 1;
    if(parts[i].y > (Dom.ye - parts[i].r)) straddle_y = 1;
    if(parts[i].z < (Dom.zs + parts[i].r)) straddle_z = 1;
    if(parts[i].z > (Dom.ze - parts[i].r)) straddle_z = 1;

    // now add the appropriate number of particles to plot
    if(straddle_x) nparts_plot += 1;
    if(straddle_y) nparts_plot += 1;
    if(straddle_z) nparts_plot += 1;
    if(straddle_x && straddle_y) nparts_plot += 1;
    if(straddle_x && straddle_z) nparts_plot += 1;
    if(straddle_y && straddle_z) nparts_plot += 1;
    if(straddle_x && straddle_y && straddle_z) nparts_plot += 1;
  }

  // create a temporary parts list containing virtual particles
  part_struct *parts_virt = (part_struct*) malloc(nparts_plot
    * sizeof(part_struct));
  // cpumem += nparts_plot * sizeof(part_struct);
  int *ind = (int*) malloc(nparts_plot * sizeof(int));
  // cpumem += nparts_plot * sizeof(int);

  int j = 0;  // virtual particle counter
  for(i = 0; i < nparts; i++) {
    // these take care of all of the particles stradding boundaries
    if(parts[i].x < (Dom.xs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    }
    // if a particle straddles two boundaries, a fourth particle is needed
    if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].z < (Dom.zs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].z < (Dom.zs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].z > (Dom.ze - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].z > (Dom.ze - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    }

    // if a particle straddles all three boundaries, an eighth is needed
    if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
      j++;
    }

    ind[j] = i;
    parts_virt[j].r = parts[i].r;
      parts_virt[j].rho = parts[i].rho;
    parts_virt[j].x = parts[i].x;
    parts_virt[j].y = parts[i].y;
    parts_virt[j].z = parts[i].z;
    parts_virt[j].u = parts[i].u;
    parts_virt[j].v = parts[i].v;
    parts_virt[j].w = parts[i].w;
    parts_virt[j].udot = parts[i].udot;
    parts_virt[j].vdot = parts[i].vdot;
    parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
    parts_virt[j].Fx = parts[i].Fx;
    parts_virt[j].Fy = parts[i].Fy;
    parts_virt[j].Fz = parts[i].Fz;
    parts_virt[j].Lx = parts[i].Lx;
    parts_virt[j].Ly = parts[i].Ly;
    parts_virt[j].Lz = parts[i].Lz;
      parts_virt[j].iFx = parts[i].iFx;
      parts_virt[j].iFy = parts[i].iFy;
      parts_virt[j].iFz = parts[i].iFz;
      parts_virt[j].iLx = parts[i].iLx;
      parts_virt[j].iLy = parts[i].iLy;
      parts_virt[j].iLz = parts[i].iLz;
      parts_virt[j].kFx = parts[i].kFx;
      parts_virt[j].kFy = parts[i].kFy;
      parts_virt[j].kFz = parts[i].kFz;
    j++;
  }

  fprintf(outfile, "<VTKFile type=\"PolyData\">\n");
  fprintf(outfile, "<PolyData>\n");
  fprintf(outfile, "<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" ",
    nparts_plot);
  fprintf(outfile, "NumberOfLines=\"0\" NumberOfStrips=\"0\" ");
  fprintf(outfile, "NumberOfPolys=\"0\">\n");
  fprintf(outfile, "<Points>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write the locations of the particle centers
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].x, parts_virt[i].y,
      parts_virt[i].z);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "</Points>\n");
  fprintf(outfile, "<PointData Scalars=\"n r\" Vectors=\"pvel pacc");
  fprintf(outfile, " pax pay paz pomega pomegadot F L Fh Fk Fi\">\n");
  fprintf(outfile, "<DataArray type=\"Int32\" Name=\"n\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write index of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%d ", ind[i]);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" Name=\"r\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write radius of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f ", parts_virt[i].r);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pvel\" format=\"ascii\">\n");

  // write velocity of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].u, parts_virt[i].v,
      parts_virt[i].w);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pacc\" format=\"ascii\">\n");

  // write acceleration of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].udot, parts_virt[i].vdot,
      parts_virt[i].wdot);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pax\" format=\"ascii\">\n");

  // write angular position of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].axx, parts_virt[i].axy,
      parts_virt[i].axz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pay\" format=\"ascii\">\n");

  // write angular position of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].ayx, parts_virt[i].ayy,
      parts_virt[i].ayz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"paz\" format=\"ascii\">\n");

  // write angular position of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].azx, parts_virt[i].azy,
      parts_virt[i].azz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pomega\" format=\"ascii\">\n");

  // write angular velocity of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].ox, parts_virt[i].oy,
      parts_virt[i].oz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pomegadot\" format=\"ascii\">\n");

  // write angular acceleration of each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ", parts_virt[i].oxdot, parts_virt[i].oydot,
      parts_virt[i].ozdot);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"F\" format=\"ascii\">\n");

  // write linear force on each particle
  for(i = 0; i < nparts_plot; i++) {
    real mass = 4./3.*PI*(parts_virt[i].rho-rho_f)
      *parts_virt[i].r*parts_virt[i].r*parts_virt[i].r;
    fprintf(outfile, "%f %f %f ",
      parts_virt[i].Fx + parts_virt[i].iFx + parts_virt[i].kFx + mass*g.x,
      parts_virt[i].Fy + parts_virt[i].iFy + parts_virt[i].kFy + mass*g.y,
      parts_virt[i].Fz + parts_virt[i].iFz + parts_virt[i].kFz + mass*g.z);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"L\" format=\"ascii\">\n");

  // write angular moment on each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      parts_virt[i].Lx + parts_virt[i].iLx,
      parts_virt[i].Ly + parts_virt[i].iLy,
      parts_virt[i].Lz + parts_virt[i].iLz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"Fh\" format=\"ascii\">\n");

  // write linear force on each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      parts_virt[i].Fx,
      parts_virt[i].Fy,
      parts_virt[i].Fz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"Fk\" format=\"ascii\">\n");

  // write linear force on each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      parts_virt[i].kFx,
      parts_virt[i].kFy,
      parts_virt[i].kFz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"Fi\" format=\"ascii\">\n");

  // write linear force on each particle
  for(i = 0; i < nparts_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      parts_virt[i].iFx,
      parts_virt[i].iFy,
      parts_virt[i].iFz);
  }

  fprintf(outfile, "\n</DataArray>\n");

  fprintf(outfile, "</PointData>\n");
  fprintf(outfile, "</Piece>\n");
  fprintf(outfile, "</PolyData>\n");
  fprintf(outfile, "</VTKFile>\n");

  // clean up
  free(parts_virt);
  free(ind);

  fclose(outfile);
}

void quadnodes_out_VTK(void)
{
  int i, j; // iterator
  char fname[FILE_NAME_SIZE] = ""; // output filename
  int nparts_plot = 0;    // the number of particles to plot; may be greater
                          // may be greater than nparts if particles straddle

  // create quadratue node locations
  real PI14 = 0.25 * PI;
  real PI12 = 0.5 * PI;
  real PI34 = 0.75 * PI;
  real PI54 = 1.25 * PI;
  real PI32 = 1.5 * PI;
  real PI74 = 1.75 * PI;
  real alph1 = 0.955316618124509; //54.736
  real alph2 = 2.186276035465284; //125.264

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

  // put all quadrature nodes together for interpolation
  real node_t[NNODES];
  real node_p[NNODES];
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


  // open file for writing
  sprintf(fname, "%s/out_nodes_%d.vtp", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // go through all particles and check if any straddle boundaries
  nparts_plot = nparts;
  for(i = 0; i < nparts; i++) {
    int straddle_x = 0;
    int straddle_y = 0;
    int straddle_z = 0;
    // figure out which boudaries this particle straddles
    if(parts[i].x < (Dom.xs + parts[i].r)) straddle_x = 1;
    if(parts[i].x > (Dom.xe - parts[i].r)) straddle_x = 1;
    if(parts[i].y < (Dom.ys + parts[i].r)) straddle_y = 1;
    if(parts[i].y > (Dom.ye - parts[i].r)) straddle_y = 1;
    if(parts[i].z < (Dom.zs + parts[i].r)) straddle_z = 1;
    if(parts[i].z > (Dom.ze - parts[i].r)) straddle_z = 1;

    // now add the appropriate number of particles to plot
    if(straddle_x) nparts_plot += 1;
    if(straddle_y) nparts_plot += 1;
    if(straddle_z) nparts_plot += 1;
    if(straddle_x && straddle_y) nparts_plot += 1;
    if(straddle_x && straddle_z) nparts_plot += 1;
    if(straddle_y && straddle_z) nparts_plot += 1;
    if(straddle_x && straddle_y && straddle_z) nparts_plot += 1;
  }

  // create a temporary parts list containing virtual particles
  part_struct *parts_virt = (part_struct*) malloc(nparts_plot
    * sizeof(part_struct));
  // cpumem += nparts_plot * sizeof(part_struct);
  int *ind = (int*) malloc(nparts_plot * sizeof(int));
  // cpumem += nparts_plot * sizeof(int);

  /***** TODO Need to include PERIODIC boolean checker for virtual particles */
  j = 0;  // virtual particle counter
  for(i = 0; i < nparts; i++) {
    // these take care of all of the particles stradding boundaries
    if(parts[i].x < (Dom.xs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    }
    // if a particle straddles two boundaries, a fourth particle is needed
    if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].z < (Dom.zs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].z < (Dom.zs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].z > (Dom.ze - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].z > (Dom.ze - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    }

    // if a particle straddles all three boundaries, an eighth is needed
    if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z < (Dom.zs + parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z + Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y < (Dom.ys + parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y + Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x < (Dom.xs + parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x + Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    } if(parts[i].x > (Dom.xe - parts[i].r)
      && parts[i].y > (Dom.ye - parts[i].r)
      && parts[i].z > (Dom.ze - parts[i].r)) {
      ind[j] = i;
      parts_virt[j].r = parts[i].r;
      parts_virt[j].x = parts[i].x - Dom.xl;
      parts_virt[j].y = parts[i].y - Dom.yl;
      parts_virt[j].z = parts[i].z - Dom.zl;
      parts_virt[j].u = parts[i].u;
      parts_virt[j].v = parts[i].v;
      parts_virt[j].w = parts[i].w;
      parts_virt[j].udot = parts[i].udot;
      parts_virt[j].vdot = parts[i].vdot;
      parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
      parts_virt[j].Fx = parts[i].Fx;
      parts_virt[j].Fy = parts[i].Fy;
      parts_virt[j].Fz = parts[i].Fz;
      parts_virt[j].Lx = parts[i].Lx;
      parts_virt[j].Ly = parts[i].Ly;
      parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
      j++;
    }

    ind[j] = i;
    parts_virt[j].r = parts[i].r;
    parts_virt[j].x = parts[i].x;
    parts_virt[j].y = parts[i].y;
    parts_virt[j].z = parts[i].z;
    parts_virt[j].u = parts[i].u;
    parts_virt[j].v = parts[i].v;
    parts_virt[j].w = parts[i].w;
    parts_virt[j].udot = parts[i].udot;
    parts_virt[j].vdot = parts[i].vdot;
    parts_virt[j].wdot = parts[i].wdot;
      parts_virt[j].axx = parts[i].axx;
      parts_virt[j].axy = parts[i].axy;
      parts_virt[j].axz = parts[i].axz;
      parts_virt[j].ayx = parts[i].ayx;
      parts_virt[j].ayy = parts[i].ayy;
      parts_virt[j].ayz = parts[i].ayz;
      parts_virt[j].azx = parts[i].azx;
      parts_virt[j].azy = parts[i].azy;
      parts_virt[j].azz = parts[i].azz;
      parts_virt[j].ox = parts[i].ox;
      parts_virt[j].oy = parts[i].oy;
      parts_virt[j].oz = parts[i].oz;
      parts_virt[j].oxdot = parts[i].oxdot;
      parts_virt[j].oydot = parts[i].oydot;
      parts_virt[j].ozdot = parts[i].ozdot;
    parts_virt[j].Fx = parts[i].Fx;
    parts_virt[j].Fy = parts[i].Fy;
    parts_virt[j].Fz = parts[i].Fz;
    parts_virt[j].Lx = parts[i].Lx;
    parts_virt[j].Ly = parts[i].Ly;
    parts_virt[j].Lz = parts[i].Lz;
      for(int k = 0; k < NNODES; k++) {
        parts_virt[j].nodes[k] = parts[i].nodes[k];
      }
    j++;
  }

  fprintf(outfile, "<VTKFile type=\"PolyData\">\n");
  fprintf(outfile, "<PolyData>\n");
  fprintf(outfile, "<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" ",
    NNODES*nparts_plot);
  fprintf(outfile, "NumberOfLines=\"0\" NumberOfStrips=\"0\" ");
  fprintf(outfile, "NumberOfPolys=\"0\">\n");
  fprintf(outfile, "<Points>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  real r, X, Y, Z;

  // write the locations of the nodes
  for(i = 0; i < nparts_plot; i++) {
    for(j = 0; j < NNODES; j++) {
      r = parts[i].rs;
      // convert node (r, theta, phi) to (x, y, z)
      X = r * sin(node_t[j]) * cos(node_p[j]);
      Y = r * sin(node_t[j]) * sin(node_p[j]);
      Z = r * cos(node_t[j]);

      // shift from particle center
      X = X + parts_virt[i].x;
      Y = Y + parts_virt[i].y;
      Z = Z + parts_virt[i].z;

      if(X < dom->xs) X = X + dom->xl;
      else if(X > dom->xe) X = X - dom->xl;
      if(Y < dom->ys) Y = Y + dom->yl;
      else if(Y > dom->ye) Y = Y - dom->yl;
      if(Z < dom->zs) Z = Z + dom->zl;
      else if(Z > dom->xe) Z = Z - dom->zl;
      
      // write to file
      fprintf(outfile, "%f %f %f ", X, Y, Z);
    }
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "</Points>\n");
  fprintf(outfile, "<PointData Scalars=\"n status\">\n");
  fprintf(outfile, "<DataArray type=\"Int32\" Name=\"n\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write index of each particle to each of the nodes
  for(i = 0; i < nparts_plot; i++) {
    for(j = 0; j < NNODES; j++) {
      fprintf(outfile, "%d ", ind[i]);
    }
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Int32\" Name=\"status\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write status of each node
  for(i = 0; i < nparts_plot; i++) {
    for(j = 0; j < NNODES; j++) {
      fprintf(outfile, "%d ", parts_virt[i].nodes[j]);
    }
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "</PointData>\n");
  fprintf(outfile, "</Piece>\n");
  fprintf(outfile, "</PolyData>\n");
  fprintf(outfile, "</VTKFile>\n");

  // clean up
  free(parts_virt);
  free(ind);

  fclose(outfile);
}
