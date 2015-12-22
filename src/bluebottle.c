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

#include <mpi.h>
#include "bluebottle.h"
#include "particle.h"
#include "precursor.h"

// define global variables that were declared in header file
int dev_start;
int dev_end;
dom_struct *dom;
dom_struct **_dom;
dom_struct Dom;
real *p0;
real *p;
real *phi;
//real *divU;
real **_p0;
real **_p;
real **_phi;
//real **_divU;
real *u;
real *u0;
real **_u;
real **_u0;
real *v;
real *v0;
real **_v;
real **_v0;
real *w;
real *w0;
real **_w;
real **_w0;
real *f_x;
real *f_y;
real *f_z;
real **_f_x;
real **_f_y;
real **_f_z;
#ifndef IMPLICIT
real *diff0_u;
real *diff0_v;
real *diff0_w;
real **_diff0_u;
real **_diff0_v;
real **_diff0_w;
#endif
real *conv0_u;
real *conv0_v;
real *conv0_w;
real **_conv0_u;
real **_conv0_v;
real **_conv0_w;
real *diff_u;
real *diff_v;
real *diff_w;
real **_diff_u;
real **_diff_v;
real **_diff_w;
real *conv_u;
real *conv_v;
real *conv_w;
real **_conv_u;
real **_conv_v;
real **_conv_w;
real **_u_star;
real **_v_star;
real **_w_star;
real *u_star;
real *v_star;
real *w_star;
real *u_WE;
real *u_SN_S;
real *u_SN_N;
real *u_BT_B;
real *u_BT_T;
real *v_WE_W;
real *v_WE_E;
real *v_SN;
real *v_BT_B;
real *v_BT_T;
real *w_WE_W;
real *w_WE_E;
real *w_SN_S;
real *w_SN_N;
real *w_BT;
real **_u_WE;
real **_u_SN_S;
real **_u_SN_N;
real **_u_BT_B;
real **_u_BT_T;
real **_v_WE_W;
real **_v_WE_E;
real **_v_SN;
real **_v_BT_B;
real **_v_BT_T;
real **_w_WE_W;
real **_w_WE_E;
real **_w_SN_S;
real **_w_SN_N;
real **_w_BT;
real **_rhs_p;
real duration;
real ttime;
real vel_tDelay;
real p_tDelay;
real g_tDelay;
real dt;
real dt0;
real CFL;
int pp_max_iter;
real pp_residual;
int lamb_max_iter;
real lamb_residual;
real lamb_relax;
real lamb_cut;
int out_plane;
int stepnum;
int rec_flow_field_stepnum_out;
int rec_paraview_stepnum_out;
int rec_particle_stepnum_out;
int rec_prec_stepnum_out;
int rec_prec_flow_field_stepnum_out;
real rec_flow_field_dt;
real rec_flow_field_ttime_out;
int rec_flow_field_vel;
int rec_flow_field_p;
int rec_flow_field_phase;
real rec_prec_flow_field_dt;
real rec_prec_flow_field_ttime_out;
int rec_prec_flow_field_vel;
int rec_prec_flow_field_p;
int rec_prec_flow_field_phase;
real rec_paraview_dt;
real rec_paraview_ttime_out;
real rec_particle_dt;
real rec_particle_ttime_out;
real rec_restart_dt;
int rec_restart_stop;
real rec_restart_ttime_out;
real rec_prec_dt;
real rec_prec_ttime_out;
int rec_particle_pos;
int rec_particle_a;
int rec_particle_vel;
int rec_particle_omega;
int rec_particle_force;
int rec_particle_moment;
g_struct g;
real rho_f;
real mu;
real nu;
BC bc;
int init_cond;
gradP_struct gradP;
real turbA;
real turbl;
long int cpumem;
long int gpumem;
int bc_flow_configs[18];
real bc_flow_vels[18];
real bc_plane_pos[9];
real pid_int;
real pid_back;
real Kp;
real Ki;
real Kd;

int main(int argc, char *argv[]) {

  int np = 0;     // number of MPI processes
  int rank = 0;   // number assigned to this MPI process
  int restart_stop = 0; // boolean to determine when to stop restart loop

  // set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  time_t startwalltime = time(NULL);
  time_t timestepwalltime = time(NULL);
  time_t diffwalltime = difftime(timestepwalltime, startwalltime);

  if(rank == MASTER) {
    int turb = 0;
    MPI_Status status; // for MPI communication

    // parse command-line arguments
    // if none given, run normally
    // if -s given, run seeder program only
    // if anything else given, exit
    // this code comes from K&R, p. 117
    int lambflag = 0;
    int argin;
    int runseeder = 0;
    int runrestart = 0;
    while(--argc > 0 && (*++argv)[0] == '-') {
      while((argin = *++argv[0])) {
        switch(argin) {
          case 's':
            runseeder = 1;
            break;
          case 'r':
            runrestart = 1;
            break;
          default:
            runseeder = 2;
            runrestart = 2;
            printf("bluebottle: illegal option %c\n", argin);
            argc = 0;
            break;
        }
      }
    }
    
    if(runseeder == 1) {
      int fret = 0;
      fret = fret;  // prevent compiler warning
      printf("Seed particles according to parameters specified in");
      printf(" parts.config? (y/N)\n");
      fflush(stdout);
      int c = getchar();
      int tmp = getchar();
      tmp = tmp;  // prevent compiler warning
      if (c == 'Y' || c == 'y') {
        printf("Seed particles for which kind of");
        printf(" array? (r)andom / (a)rray / (h)ex / (p)erturbed?\n");
        fflush(stdout);
        int type = getchar();
        tmp = getchar();
        int Nx = 0; int Ny = 0; int Nz = 0; 
        double ddz = 0.0; double bias = 0.0; int nperturb = 0;
        if(type == 'r'){
          seeder_read_input(Nx, Ny, Nz, ddz, bias, nperturb);
        }
        if(type == 'a'){
          printf("Please input the number of particles in the x direction\n");
          fflush(stdout);
          fret = scanf("%d",&Nx);
          printf("Please input the number of particles in the y direction\n");
          fflush(stdout);
          fret = scanf("%d",&Ny);
          printf("Please input the number of particles in the z direction\n");
          fflush(stdout);
          fret = scanf("%d",&Nz);
          printf("Nx Ny Nz is: %d %d %d\n",Nx, Ny, Nz);  
          seeder_read_input(Nx, Ny, Nz, ddz, bias, nperturb);               
        }
        if(type == 'h'){
          printf("Please input the number of particles in the x direction\n");
          fflush(stdout);
          fret = scanf("%d",&Nx);
          printf("Please input the number of particles in the y direction\n");
          fflush(stdout);
          fret = scanf("%d",&Ny);
          printf("Please input the number of particles in the z direction\n");
          fflush(stdout);
          fret = scanf("%d",&Nz);
          printf("Please input the layer distance for hex array\n");
          fflush(stdout);
          fret = scanf("%lf",&ddz);
          printf("Nx Ny Nz ddz is: %d %d %d %lf\n",Nx, Ny, Nz, ddz);
          seeder_read_input(Nx, Ny, Nz, ddz, bias, nperturb);  
        }
        if(type == 'p'){
          printf("Please input the particle number in x direction\n");
          fflush(stdout);
          fret = scanf("%d",&Nx);
          printf("Please input the particle number in y direction\n");
          fflush(stdout);
          fret = scanf("%d",&Ny);
          printf("Please input the particle number in z direction\n");
          fflush(stdout);
          fret = scanf("%d",&Nz);
          printf("Please input the perturbation magnitude");
          printf(" (between 0 and 1)\n");
          fflush(stdout);
          fret = scanf("%lf",&bias);
          printf("Please input the number of perturbations");
          printf(" (should be larger than 100,000)\n");
          fflush(stdout);
          fret = scanf("%d",&nperturb);
          printf("Nx Ny Nz bias nperturb is: %d %d %d %lf %d\n",Nx, Ny, Nz, bias,
            nperturb);
          seeder_read_input(Nx, Ny, Nz, ddz, bias, nperturb);
        }
        return EXIT_SUCCESS;
      } else {
        printf("Please specify the desired parameters in parts.config\n\n");
        fflush(stdout);
        return EXIT_FAILURE;
      }
    } else if(runrestart == 1 && argc > 0) {
      printf("Usage restart simulation: bluebottle -r\n");
      return EXIT_FAILURE;
    } else if(runseeder == 2) {
      return EXIT_FAILURE;
    } else if(runrestart == 2) {
      return EXIT_FAILURE;
    } else {
      // read recorder config file
      recorder_read_config();

      // read simulation input configuration file
      printf("\nRunning bluebottle_0.1...\n\n");
      printf("Reading the domain and particle input files...\n\n");
      domain_read_input();
      parts_read_input(turb);
      fflush(stdout);
      //printf("EXPD: Using devices %d through %d.\n\n", dev_start, dev_end);
      //fflush(stdout);

      /********* Messy hack for taking advantage of CUDA_VISIBLE_DEVICES
       ********* treatment in SLURM resource manager. */

      // read the environment variable
      char *cvdin;
      cvdin = getenv("CUDA_VISIBLE_DEVICES");
      if(cvdin != NULL) {
        // number of devices
        int n_CUDA_VISIBLE_DEVICES = 0.5*(strlen(cvdin)+1.);
        // list of devices
        int *CUDA_VISIBLE_DEVICES = malloc(n_CUDA_VISIBLE_DEVICES*sizeof(int));
        // fill list of devices assuming single-character separation
        int j = 0;
        for(int i = 0; i < 2*n_CUDA_VISIBLE_DEVICES-1; i+=2) {
          CUDA_VISIBLE_DEVICES[j] = cvdin[i] - '0';
          j++;
        }
        // use the first device available (devices are re-indexed by CUDA)
        if(n_CUDA_VISIBLE_DEVICES > 0) {
          dev_start = 0;
          dev_end = 0;
        } else { // exit if things aren't just right
          printf("Environment variable CUDA_VISIBLE_DEVICES is empty:\n");
          printf("  a. To use the config files to specify the device number,\n");
          printf("     type 'unset CUDA_VISIBLE_DEVICES'\n");
          printf("  b. To use CUDA_VISIBLE_DEVICES to specify the device number,\n");
          printf("     type 'export CUDA_VISIBLE_DEVICES=N1,N2,...',\n");
          printf("     where N1,N2 are comma-separated device numbers\n");
          exit(EXIT_FAILURE);
        }
      }

      /********* End messy CUDA_VISIBLE_DEVICES hack. */

      if(runrestart != 1) {
        // start BICGSTAB recorder
        recorder_bicgstab_init("solver_expd.rec");
        #ifdef IMPLICIT
          // start Helmholtz recorder
          recorder_bicgstab_init("solver_helmholtz_expd.rec");
        #endif
        // start Lamb's coefficient recorder
        // commented out because it should now automatically init itself
        // from recorder_lamb(...) if the file doesn't already exist
        //recorder_lamb_init("lamb.rec");
      }

      // initialize the domain
      printf("Initializing domain variables...");
      fflush(stdout);
      int domain_init_flag = domain_init();
      printf("done.\n");
      fflush(stdout);
      if(domain_init_flag == EXIT_FAILURE) {
        printf("\nThe number of devices in DEV RANGE is insufficient\n");
        printf("for the given domain decomposition.  Exiting now.\n");
        return EXIT_FAILURE;
      }

      // set up the boundary condition config info to send to precursor
      expd_init_BC(np);

      // initialize the particles
      printf("Initializing particle variables...");
      fflush(stdout);
      int parts_init_flag = parts_init();
      int binDom_init_flag = binDom_init();
      printf("done.\n");
      fflush(stdout);
      if(parts_init_flag == EXIT_FAILURE) {
        printf("\nThe initial particle configuration is not allowed.\n");
        return EXIT_FAILURE;
      } else if(binDom_init_flag == EXIT_FAILURE) {
        printf("\nThe bin configuration is not allowed.\n");
        return EXIT_FAILURE;
      }

      // allocate device memory
      printf("Allocating domain CUDA device memory...");
      fflush(stdout);
      cuda_dom_malloc();
      printf("...done.\n");
      fflush(stdout);
      printf("Allocating particle CUDA device memory...");
      fflush(stdout);
      cuda_part_malloc();
      printf("...done.\n");
      fflush(stdout);

      // copy host data to devices
      printf("Copying host domain data to devices...");
      fflush(stdout);
      cuda_dom_push();
      printf("done.\n");
      fflush(stdout);
      printf("Copying host particle data to devices...");
      fflush(stdout);
      cuda_part_push();
      printf("done.\n");
      fflush(stdout);

      count_mem();

      // initialize ParaView VTK output PVD file
      if(runrestart != 1) {
        #ifdef DDEBUG
          init_VTK_ghost();
        #else
          if(rec_paraview_dt > 0) {
            init_VTK();
          }
        #endif
      }

      // set up particles
      cuda_build_cages();
      cuda_part_pull();

      // run restart if requested
      if(runrestart == 1) {
        printf("\nRestart requested.\n\n");
        printf("Reading restart file...");
        fflush(stdout);
        in_restart();
        printf("done.\n");
        fflush(stdout);
        printf("Copying host domain data to devices...");
        fflush(stdout);
        cuda_dom_push();
        printf("done.\n");
        fflush(stdout);
        printf("Copying host particle data to devices...");
        fflush(stdout);
        cuda_part_push();
        printf("done.\n");
        fflush(stdout);
        cgns_grid();
        if(ttime >= duration) {
          printf("\n...simulation completed.\n");
          restart_stop = 1;
        }
      }

      #ifdef DDEBUG
        // write config to screen
        printf("\n=====DEBUG");
        printf("================================");
        printf("======================================\n");
        fflush(stdout);
        cuda_dom_pull();
        cuda_part_pull();
        domain_show_config();
        parts_show_config();
        bin_show_config();
        printf("========================================");
        printf("========================================\n\n");
        fflush(stdout);
      #endif

      #ifdef TEST // run test code
        // ** note that some of these work only for DEV RANGE 0 0 **
        // test CUDA functionality
        printf("\n=====TEST");
        printf("=================================");
        printf("======================================\n");
        fflush(stdout);
        dt = 1.;
        dt0 = -1.;
        cuda_compute_forcing(&pid_int, &pid_back, Kp, Ki, Kd);
        rec_flow_field_stepnum_out = -1;
        rec_paraview_stepnum_out = -1;
        rec_particle_stepnum_out = -1;
        //rec_restart_stepnum_out = -1;
        rec_prec_stepnum_out = -1;
        cuda_part_pull();
        //cuda_BC_test();
        //cuda_U_star_test_exp();
        //cuda_U_star_test_cos();
        //cuda_project_test();
        //cuda_quad_interp_test();
        cuda_lamb_test();
        printf("========================================");
        printf("========================================\n\n");
        fflush(stdout);

      #else // run simulation
        // begin simulation
        printf("\n=====BLUEBOTTLE");
        printf("===========================");
        printf("======================================\n");
        fflush(stdout);

        // get initial dt; this is an extra check for the SHEAR initialization
        dt = cuda_find_dt();
	printf("cuda_find_dt\n");
        // share this with the precursor domain
        expd_compare_dt(np, status);
	printf("expd_compare_dt\n");
        // update the boundary condition config info to share with precursor
        expd_update_BC(np, status);
        // apply boundary conditions to field variables
        if(nparts > 0) {
          cuda_part_BC();
        }
        cuda_dom_BC();

        // write particle internal flow equal to solid body velocity
        cuda_parts_internal();
        cuda_dom_BC();

        // write initial fields
        if(runrestart != 1) {
          cuda_dom_pull();
          cuda_part_pull();

        #ifdef DDEBUG
            printf("Writing ParaView file %d (t = %e)...",
            rec_paraview_stepnum_out, ttime);
            fflush(stdout);
            out_VTK_ghost();
            rec_paraview_stepnum_out++;
            printf("done.               \n");
            fflush(stdout);
        #else
          if(rec_flow_field_dt > 0) {
            printf("Writing flow field file t = %e...", ttime);
            fflush(stdout);
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
            printf("done.               \n");
            fflush(stdout);
          }
          if(rec_particle_dt > 0) {
            printf("Writing particle file t = %e...", ttime);
            fflush(stdout);
            cgns_particles(rec_particle_dt);
            recorder_lamb("lamb.rec", 0);
            rec_particle_stepnum_out++;
            printf("done.               \n");
            fflush(stdout);
          }
          if(rec_paraview_dt > 0) {
            printf("Writing ParaView file %d (t = %e)...",
              rec_paraview_stepnum_out, ttime);
            fflush(stdout);
            out_VTK();
            rec_paraview_stepnum_out++;
            printf("done.               \n");
            fflush(stdout);
          }

        #endif
        }

        /******************************************************************/
        /** Begin the main timestepping loop in the experimental domain. **/
        /******************************************************************/
        while(ttime <= duration) {
          ttime += dt;
          rec_flow_field_ttime_out += dt;
          rec_paraview_ttime_out += dt;
          rec_particle_ttime_out += dt;
          rec_restart_ttime_out += dt;
          stepnum++;
          printf("EXPD: Time = %e of %e (dt = %e).\n", ttime, duration, dt);
          fflush(stdout);

          cuda_compute_forcing(&pid_int, &pid_back, Kp, Ki, Kd);
          compute_vel_BC();
          // update the boundary condition config info and share with precursor
          expd_update_BC(np, status);
          // TODO: save work by rebuilding only the cages that need to be rebuilt
          cuda_build_cages();

          int iter = 0;
          real iter_err = FLT_MAX;

          while(iter_err > lamb_residual) {  // iterate for Lamb's coefficients
            #ifndef BATCHRUN
              printf("  Iteration %d: ", iter);
              fflush(stdout);
            #endif

            // solve for U_star
            #ifndef IMPLICIT
              cuda_U_star_2();
            #else
              cuda_ustar_helmholtz(rank);
              cuda_vstar_helmholtz(rank);
              cuda_wstar_helmholtz(rank);
            #endif

            // apply boundary conditions to U_star
            if(nparts > 0) {
              cuda_part_BC_star();
            }
            cuda_dom_BC_star();
            // enforce solvability condition
            cuda_solvability();
            if(nparts > 0) {
              cuda_part_BC_star();
            }
            cuda_dom_BC_star();
            // solve for pressure
            cuda_PP_bicgstab(rank);
            cuda_dom_BC_phi();
            // solve for U
            cuda_project();
            // apply boundary conditions to field variables
            if(nparts > 0) {
              cuda_part_BC();
            }
            cuda_dom_BC();
            // update pressure
            cuda_update_p();
            if(nparts > 0) {
              cuda_part_BC();
              cuda_part_p_fill();
            }
            cuda_dom_BC_p();

            // update Lamb's coefficients
            cuda_move_parts_sub();
            cuda_Lamb();

            #ifdef STEPS // force no sub-timestep iteration
              iter_err = -1;
            #else
              // check error between this set of coefficients and previous set
              // of coefficients
              iter_err = cuda_lamb_err();
              // TODO: write error to lamb.rec
            #endif
            #ifndef BATCHRUN
              printf("Error = %f\r", iter_err);
            #endif
            iter++;
            // check iteration limit
            if(iter == lamb_max_iter) {
              //lambflag = !lambflag;
              //printf("Reached the maximum number of Lamb's");
              //printf(" coefficient iterations.");
              //printf(" CONTINUING simulation.\n");
              break;
            }
          }

          printf("  The Lamb's coefficients converged in");
          printf(" %d iterations.\n", iter);

          if(!lambflag) {
            // update particle position
            cuda_move_parts();

            // write particle internal flow equal to solid body velocity
            cuda_parts_internal();
            cuda_dom_BC();

            // store u, conv, and coeffs for use in next timestep
            cuda_store_u();
            if(nparts > 0)
              cuda_store_coeffs();

            // compute div(U)
            //cuda_div_U();

            // compute next timestep size
            dt0 = dt;
            dt = cuda_find_dt();

            // compare this timestep size to that in the precursor and
            // and synchronize the result
            expd_compare_dt(np, status);

          } else {
            return EXIT_FAILURE;
          }

          if(rec_flow_field_dt > 0) {
            if(rec_flow_field_ttime_out >= rec_flow_field_dt) {
              // pull back data and write fields
              cuda_dom_pull();
              cuda_part_pull();
              #ifndef BATCHRUN
                printf("  Writing flow field file t = %e...                  \r",
                  ttime);
                fflush(stdout);
              #endif
              cgns_flow_field(rec_flow_field_dt);
              printf("  Writing flow field file t = %e...done.\n", ttime);
              fflush(stdout);
              rec_flow_field_ttime_out = rec_flow_field_ttime_out
                - rec_flow_field_dt;
              rec_flow_field_stepnum_out++;
            }
          }
          if(rec_paraview_dt > 0) {
            if(rec_paraview_ttime_out >= rec_paraview_dt) {
              // pull back data and write fields
              cuda_dom_pull();
              cuda_part_pull();
              #ifndef BATCHRUN
                printf("  Writing ParaView output file");
                printf(" %d (t = %e)...                  \r",
                  rec_paraview_stepnum_out, ttime);
                fflush(stdout);
              #endif
              #ifdef DDEBUG
                out_VTK_ghost();
              #else
                out_VTK();
              #endif
              printf("  Writing ParaView file %d (t = %e)...done.\n",
                rec_paraview_stepnum_out, ttime);
              rec_paraview_stepnum_out++;
              fflush(stdout);
              rec_paraview_ttime_out = rec_paraview_ttime_out - rec_paraview_dt;
            }
          }
          if(rec_particle_dt > 0) {
            if(rec_particle_ttime_out >= rec_particle_dt) {
              // pull back data and write fields
              cuda_part_pull();
              #ifndef BATCHRUN
                printf("  Writing particle file t = %e...                  \r",
                  ttime);
                fflush(stdout);
              #endif
              #ifdef DDEBUG
                recorder_lamb("lamb.rec", iter);
              #else
                cgns_particles(rec_particle_dt);
                recorder_lamb("lamb.rec", iter);
              #endif
              printf("  Writing particle file t = %e...done.\n", ttime);
              fflush(stdout);
              rec_particle_ttime_out = rec_particle_ttime_out - rec_particle_dt;
              rec_particle_stepnum_out++;
            }
          }

          // write a restart file and exit when the time is appropriate
          timestepwalltime = time(NULL);
          diffwalltime = difftime(timestepwalltime, startwalltime);
          int rest_com = (rec_restart_dt > 0)
            && ((real)diffwalltime/60. > rec_restart_dt);

          // communicate write restart with precursor domain
          expd_comm_restart_write(np, rest_com);

          if(rest_com) {
            printf("  Writing restart file (t = %e)...", ttime);
            fflush(stdout);
            cuda_dom_pull();
            cuda_part_pull();
            out_restart();
            printf("done.               \n");
            fflush(stdout);
            rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
            startwalltime = time(NULL);
            if(rec_restart_stop)
              break; // exit!
          }

          // check for blow-up condition
          if(dt < 1e-20) {
            printf("The solution has diverged.  Ending simulation.              \n");
            return EXIT_FAILURE;
          }
        }

        if(rec_restart_dt > 0 && ttime >= duration && !restart_stop) {
          printf("  Writing final restart file (t = %e)...", ttime);
          fflush(stdout);
          cuda_dom_pull();
          cuda_part_pull();
          out_restart();
          printf("done.               \n");
          fflush(stdout);
          rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
          startwalltime = time(NULL);
        }

        printf("========================================");
        printf("========================================\n\n");
        fflush(stdout);
      #endif

      // clean up devices
      printf("Cleaning up domain data on devices...");
      fflush(stdout);
      cuda_dom_free();
      printf("done.     \n");
      fflush(stdout);
      printf("Cleaning up particle data on devices...");
      fflush(stdout);
      cuda_part_free();
      printf("done.\n");
      fflush(stdout);

      // clean up host
      printf("Cleaning up particles...");
      fflush(stdout);
      parts_clean();
      printf("done.\n");
      fflush(stdout);
      printf("Cleaning up domain...");
      fflush(stdout);
      domain_clean();
      printf("done.\n");
      fflush(stdout);

      printf("\n...bluebottle_0.1 done.\n\n");
    }
  } else {
    int turb = 1;   // boolean
    MPI_Status status; // for MPI communication

    // parse command-line arguments
    // if none given, run normally
    // if -s given, run seeder program only
    // if anything else given, exit
    // this code comes from K&R, p. 117
    int argin;
    int runrestart = 0;
    while(--argc > 0 && (*++argv)[0] == '-') {
      while((argin = *++argv[0])) {
        switch(argin) {
          case 'r':
            runrestart = 1;
            break;
          default:
            runrestart = 2;
            printf("bluebottle: illegal option %c\n", argin);
            argc = 0;
            break;
        }
      }
    }

    // read simulation input configuration file
    recorder_read_config();
    turb_read_input();
    parts_read_input(turb);

    //printf("PREC: Using devices %d through %d.\n\n", dev_start, dev_end);
    //fflush(stdout);

    /********* Messy hack for taking advantage of CUDA_VISIBLE_DEVICES
     ********* treatment in SLURM resource manager. */

    // read th environment variable
    char *cvdin;
    cvdin = getenv("CUDA_VISIBLE_DEVICES");
    if(cvdin != NULL) {
      // number of devices
      int n_CUDA_VISIBLE_DEVICES = 0.5*(strlen(cvdin)+1.);
      // list of devices
      int *CUDA_VISIBLE_DEVICES = malloc(n_CUDA_VISIBLE_DEVICES*sizeof(int));
      // fill list of devices assuming single-character separation
      int j = 0;
      for(int i = 0; i < 2*n_CUDA_VISIBLE_DEVICES-1; i+=2) {
        CUDA_VISIBLE_DEVICES[j] = cvdin[i] - '0';
        j++;
      }
      // use the second device available (devices are re-indexed by CUDA)
      if(n_CUDA_VISIBLE_DEVICES > 1) {
        dev_start = 1;
        dev_end = 1;
      // if only one device is available, try to use it
      } else if(n_CUDA_VISIBLE_DEVICES > 0) {
        dev_start = 0;
        dev_end = 0;
      } else { // exit if things aren't just right
        printf("Environment variable CUDA_VISIBLE_DEVICES is empty:\n");
        printf("  a. To use the config files to specify the device number,\n");
        printf("     type 'unset CUDA_VISIBLE_DEVICES'\n");
        printf("  b. To use CUDA_VISIBLE_DEVICES to specify the device number,\n");
        printf("     type 'export CUDA_VISIBLE_DEVICES=N1,N2,...',\n");
        printf("     where N1,N2 are comma-separated device numbers\n");
        exit(EXIT_FAILURE);
      }
    }

    /********* End messy CUDA_VISIBLE_DEVICES hack. */

    if(runrestart != 1) {
      // start BICGSTAB recorder
      recorder_bicgstab_init("solver_prec.rec");
      #ifdef IMPLICIT
        // start Helmholtz recorder
        recorder_bicgstab_init("solver_helmholtz_prec.rec");
      #endif
    }

    // initialize the domain
    int domain_init_flag = domain_init_turb();
    if(domain_init_flag == EXIT_FAILURE) {
      printf("\nThe number of devices in DEV RANGE is insufficient\n");
      printf("for the given turbulence domain decomposition.  Exiting now.\n");
      return EXIT_FAILURE;
    }

    // receive the boundary condition config info from MASTER
    prec_init_BC(np, rank, status);

    // initialize the particles
    int parts_init_flag = parts_init();
    int binDom_init_flag = binDom_init();
    if(parts_init_flag == EXIT_FAILURE) {
      printf("\nThe initial particle configuration is not allowed.\n");
      return EXIT_FAILURE;
    } else if(binDom_init_flag == EXIT_FAILURE) {
      printf("\nThe bin configuration is not allowed.\n");
      return EXIT_FAILURE;
    }

    // allocate device memory
    cuda_dom_malloc();
    cuda_part_malloc();

    // copy host data to devices
    cuda_dom_push();
    cuda_dom_turb_planes_push(bc_flow_configs);
    cuda_part_push();

    //count_mem();

    // initialize ParaView VTK output PVD file
    if(rec_prec_dt > 0) {
      init_VTK_turb();
    }

    real rec_prec_ttime_out = 0.;
    real rec_restart_ttime_out = 0.;

    cuda_build_cages();
    cuda_part_pull();

    // run restart if requested
    if(runrestart == 1) {
      printf("\nRestart requested.\n\n");
      printf("Reading restart file...");
      fflush(stdout);
      cgns_turb_grid();
      in_restart_turb();
      printf("done.\n");
      fflush(stdout);
      printf("Copying host domain data to devices...");
      fflush(stdout);
      cuda_dom_push();
      printf("done.\n");
      fflush(stdout);
      printf("Copying host particle data to devices...");
      fflush(stdout);
      cuda_part_push();
      printf("done.\n");
      fflush(stdout);
      cgns_grid();
    }

    // initialize timestep size since turbulent velocity is nonzero
    dt = cuda_find_dt();

    // share this dt with the experimental domain
    prec_compare_dt(np, rank, status);

    // update the boundary conditions according to the experimental domain
    prec_update_BC(np, rank, status);

    // begin simulation
    // apply boundary conditions to field variables
    cuda_dom_BC();

    // write initial fields
    if(rec_prec_dt > 0 && runrestart != 1) {
      cuda_dom_pull();
      printf("Writing precursor file %d (t = %e)...",
        rec_prec_stepnum_out, ttime);
      fflush(stdout);
      out_VTK_turb();
      rec_prec_stepnum_out++;
      printf("done.               \n");
      fflush(stdout);
    }
    if(rec_prec_flow_field_dt > 0 && runrestart != 1) {
      cuda_dom_pull();
      printf("Writing turbulence flow field file t = %e...", ttime);
      fflush(stdout);
      cgns_turb_grid();
      cgns_turb_flow_field(rec_prec_flow_field_dt);
      rec_prec_flow_field_stepnum_out++;
      printf("done.               \n");
      fflush(stdout);
    }
    /***************************************************************/
    /** Begin the main timestepping loop in the precursor domain. **/
    /***************************************************************/
    while(ttime <= duration) {
      ttime += dt;
      rec_prec_flow_field_ttime_out += dt;
      rec_prec_ttime_out += dt;
      rec_restart_ttime_out += dt;
      stepnum++;
      printf("PREC: Time = %e of %e (dt = %e).\n", ttime, duration, dt);

      cuda_compute_forcing(&pid_int, &pid_back, Kp, Ki, Kd);
      cuda_compute_turb_forcing();
      compute_vel_BC();

      // solve for U_star
      #ifndef IMPLICIT
        cuda_U_star_2();
      #else
        cuda_ustar_helmholtz(rank);
        cuda_vstar_helmholtz(rank);
        cuda_wstar_helmholtz(rank);
      #endif
      // apply boundary conditions to U_star
      cuda_dom_BC_star();
      // force solvability condition
      cuda_solvability();
      cuda_dom_BC_star();
      // solve for pressure
      cuda_PP_bicgstab(rank);
      cuda_dom_BC_phi();
      // solve for U
      cuda_project();
      // apply boundary conditions to field variables
      cuda_dom_BC();
      // update pressure
      cuda_update_p();
      cuda_dom_BC_p();

      cuda_store_u();

      // compute next timestep size
      dt0 = dt;
      dt = cuda_find_dt();
     
      // check for blow-up condition
      if(dt < 1e-20) {
        printf("The solution has diverged.  Ending simulation.              \n");
        return EXIT_FAILURE;
      }

      // communicate the boundary condition with the experimental domain
      prec_send_BC(np, rank, status);

      if(rec_prec_dt > 0) {
        if(rec_prec_ttime_out >= rec_prec_dt) {
          // pull back data and write fields
          cuda_dom_pull();
          #ifndef BATCHRUN
            printf("  Writing precursor output file");
            printf(" %d (t = %e)...                  \r",
              rec_prec_stepnum_out, ttime);
            fflush(stdout);
          #endif
          #ifdef DDEBUG
            out_VTK_ghost();
          #else
            out_VTK_turb();
          #endif
          printf("  Writing precursor file %d (t = %e)...done.\n",
            rec_prec_stepnum_out, ttime);
          rec_prec_stepnum_out++;
          fflush(stdout);
          rec_prec_ttime_out = rec_prec_ttime_out - rec_prec_dt;
        }
      }
      if(rec_prec_flow_field_dt > 0) {
       if(rec_prec_flow_field_ttime_out >= rec_prec_flow_field_dt) {
          // pull back data and write fields
          cuda_dom_pull();
          cgns_turb_flow_field(rec_prec_flow_field_dt);
          printf("  Writing precursor flow field file t = %e...done.\n",
            ttime);
          fflush(stdout);
          rec_prec_flow_field_ttime_out = rec_prec_flow_field_ttime_out
            - rec_prec_flow_field_dt;
          rec_prec_flow_field_stepnum_out++;
        }
      }
      // write a restart file and exit when the time is appropriate
      int rest_com;
      prec_comm_restart_write(np, &rest_com, rank, status);
      if(rest_com) {
        printf("  Writing precursor restart file (t = %e)...", ttime);
        fflush(stdout);
        cuda_dom_pull();
        out_restart_turb();
        printf("done.               \n");
        fflush(stdout);
        rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
        startwalltime = time(NULL);
        if(rec_restart_stop)
          break; // exit!
      }
    }

    if(rec_restart_dt > 0 && ttime >= duration && !restart_stop) {
      printf("  Writing final precursor restart file (t = %e)...", ttime);
      fflush(stdout);
      cuda_dom_pull();
      out_restart_turb();
      printf("done.               \n");
      fflush(stdout);
      rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
      startwalltime = time(NULL);
    }

    // clean up devices
    cuda_dom_free();
    cuda_part_free();

    // clean up host
    parts_clean();
    domain_clean();
  }

  // finalize MPI
  MPI_Finalize();

  if(restart_stop) return EXIT_FAILURE;
  else return EXIT_SUCCESS;
}
