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

#include "precursor.h"

void expd_init_BC(int np)
{
  bc_flow_configs[ 0] = bc.uW;    bc_flow_vels[ 0] = bc.uWD;
  bc_flow_configs[ 1] = bc.uE;    bc_flow_vels[ 1] = bc.uED;
  bc_flow_configs[ 2] = bc.uS;    bc_flow_vels[ 2] = bc.uSD;
  bc_flow_configs[ 3] = bc.uN;    bc_flow_vels[ 3] = bc.uND;
  bc_flow_configs[ 4] = bc.uB;    bc_flow_vels[ 4] = bc.uBD;
  bc_flow_configs[ 5] = bc.uT;    bc_flow_vels[ 5] = bc.uTD;
  bc_flow_configs[ 6] = bc.vW;    bc_flow_vels[ 6] = bc.vWD;
  bc_flow_configs[ 7] = bc.vE;    bc_flow_vels[ 7] = bc.vED;
  bc_flow_configs[ 8] = bc.vS;    bc_flow_vels[ 8] = bc.vSD;
  bc_flow_configs[ 9] = bc.vN;    bc_flow_vels[ 9] = bc.vND;
  bc_flow_configs[10] = bc.vB;    bc_flow_vels[10] = bc.vBD;
  bc_flow_configs[11] = bc.vT;    bc_flow_vels[11] = bc.vTD;
  bc_flow_configs[12] = bc.wW;    bc_flow_vels[12] = bc.wWD;
  bc_flow_configs[13] = bc.wE;    bc_flow_vels[13] = bc.wED;
  bc_flow_configs[14] = bc.wS;    bc_flow_vels[14] = bc.wSD;
  bc_flow_configs[15] = bc.wN;    bc_flow_vels[15] = bc.wND;
  bc_flow_configs[16] = bc.wB;    bc_flow_vels[16] = bc.wBD;
  bc_flow_configs[17] = bc.wT;    bc_flow_vels[17] = bc.wTD;

  // share configuration with precursor domain
  // only one precursor domain at rank = 1 for now
  if(np > 1) {
    MPI_Send(&bc_flow_configs, 18, MPI_INT, 1, 1, MPI_COMM_WORLD);
    MPI_Send(&bc_flow_vels, 18, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
  }
}

void expd_update_BC(int np, MPI_Status status)
{
  if(np > 1) {
    // update the boundary condition config info to send to precursor
    bc_flow_vels[ 0] = bc.uWD;
    bc_flow_vels[ 1] = bc.uED;
    bc_flow_vels[ 2] = bc.uSD;
    bc_flow_vels[ 3] = bc.uND;
    bc_flow_vels[ 4] = bc.uBD;
    bc_flow_vels[ 5] = bc.uTD;
    bc_flow_vels[ 6] = bc.vWD;
    bc_flow_vels[ 7] = bc.vED;
    bc_flow_vels[ 8] = bc.vSD;
    bc_flow_vels[ 9] = bc.vND;
    bc_flow_vels[10] = bc.vBD;
    bc_flow_vels[11] = bc.vTD;
    bc_flow_vels[12] = bc.wWD;
    bc_flow_vels[13] = bc.wED;
    bc_flow_vels[14] = bc.wSD;
    bc_flow_vels[15] = bc.wND;
    bc_flow_vels[16] = bc.wBD;
    bc_flow_vels[17] = bc.wTD;

    // only one precursor domain at rank = 1 for now
    MPI_Send(&bc_flow_vels, 18, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);

    // receive the planes to be used as inflow velocities
    if(bc_flow_configs[ 0] == PRECURSOR || bc_flow_configs[ 1] == PRECURSOR)
      MPI_Recv(u_WE, Dom.Gfx.jnb*Dom.Gfx.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[ 2] == PRECURSOR || bc_flow_configs[ 3] == PRECURSOR)
      MPI_Recv(u_SN_S, Dom.Gfx.inb*Dom.Gfx.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
      MPI_Recv(u_SN_N, Dom.Gfx.inb*Dom.Gfx.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[ 4] == PRECURSOR || bc_flow_configs[ 5] == PRECURSOR)
      MPI_Recv(u_BT_B, Dom.Gfx.inb*Dom.Gfx.jnb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
      MPI_Recv(u_BT_T, Dom.Gfx.inb*Dom.Gfx.jnb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[ 6] == PRECURSOR || bc_flow_configs[ 7] == PRECURSOR)
      MPI_Recv(v_WE_W, Dom.Gfy.jnb*Dom.Gfy.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
      MPI_Recv(v_WE_E, Dom.Gfy.jnb*Dom.Gfy.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[ 8] == PRECURSOR || bc_flow_configs[ 9] == PRECURSOR)
      MPI_Recv(v_SN, Dom.Gfy.inb*Dom.Gfy.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[10] == PRECURSOR || bc_flow_configs[11] == PRECURSOR)
      MPI_Recv(v_BT_B, Dom.Gfy.inb*Dom.Gfy.jnb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
      MPI_Recv(v_BT_T, Dom.Gfy.inb*Dom.Gfy.jnb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[12] == PRECURSOR || bc_flow_configs[13] == PRECURSOR)
      MPI_Recv(w_WE_W, Dom.Gfz.jnb*Dom.Gfz.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
      MPI_Recv(w_WE_E, Dom.Gfz.jnb*Dom.Gfz.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[14] == PRECURSOR || bc_flow_configs[15] == PRECURSOR)
      MPI_Recv(w_SN_S, Dom.Gfz.inb*Dom.Gfz.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
      MPI_Recv(w_SN_N, Dom.Gfz.inb*Dom.Gfz.knb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);
    if(bc_flow_configs[16] == PRECURSOR || bc_flow_configs[17] == PRECURSOR)
      MPI_Recv(w_BT, Dom.Gfz.inb*Dom.Gfz.jnb, MPI_DOUBLE, 1, 1,
        MPI_COMM_WORLD, &status);

    // push planes to devices
    cuda_dom_turb_planes_push(bc_flow_configs);
  }
}

void expd_compare_dt(int np, MPI_Status status)
{
  if(np > 1) {
    // receive dt from other simulation for comparison
    // only have one other domain now, so this is easy
    real dt_turb = 0;
    //printf("Checking timestep requested for each simulation...\n");
    //fflush(stdout);
    MPI_Recv(&dt_turb, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &status);
    //printf("Comparing dt_turb = %e, dt = %e...\n", dt_turb, dt);
    //fflush(stdout);
    // do comparison
    if(dt_turb < dt) dt = dt_turb;
    //printf("Choosing dt = %e...\n", dt);
    //fflush(stdout);
    // send the result back
    // only one precursor domain at rank = 1 for now
    MPI_Send(&dt, 1, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD);
  }
}

void prec_init_BC(int np, int rank, MPI_Status status)
{
  if(np > 1) {
    MPI_Recv(&bc_flow_configs, 18, MPI_INT, MASTER, rank,
      MPI_COMM_WORLD, &status);
    MPI_Recv(&bc_flow_vels, 18, MPI_DOUBLE, MASTER, rank,
      MPI_COMM_WORLD, &status);
  }
  for(int i = 0; i < 9; i++) bc_plane_pos[i] = 0.;
}

void prec_update_BC(int np, int rank, MPI_Status status)
{
  if(np > 1) {

    MPI_Recv(&bc_flow_vels, 18, MPI_DOUBLE, MASTER, rank,
      MPI_COMM_WORLD, &status);

    // yank the appropriate planes for the inflow from the precursor domain
    // and put into their own plane arrays
    cuda_yank_turb_planes(bc_flow_configs, bc_plane_pos, bc_flow_vels);

    // update plane position
    if(bc_flow_configs[ 0] == PRECURSOR) {
      bc_plane_pos[ 0] = bc_plane_pos[ 0] - bc_flow_vels[ 0] * dt;
      if(bc_plane_pos[ 0] > Dom.xe)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] - Dom.xl;
      else if(bc_plane_pos[ 0] < Dom.xs)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] - Dom.xl;
    } else if(bc_flow_configs[ 1] == PRECURSOR) {
      bc_plane_pos[ 0] = bc_plane_pos[ 0] - bc_flow_vels[ 1] * dt;
      if(bc_plane_pos[ 0] > Dom.xe)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] - Dom.xl;
      else if(bc_plane_pos[ 0] < Dom.xs)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] - Dom.xl;
    }
    if(bc_flow_configs[ 2] == PRECURSOR) {
      bc_plane_pos[ 1] = bc_plane_pos[ 1] - bc_flow_vels[ 2] * dt;
      if(bc_plane_pos[ 1] > Dom.xe)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] - Dom.xl;
      else if(bc_plane_pos[ 1] < Dom.xs)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] - Dom.xl;
    } else if(bc_flow_configs[ 3] == PRECURSOR) {
      bc_plane_pos[ 1] = bc_plane_pos[ 1] - bc_flow_vels[ 3] * dt;
      if(bc_plane_pos[ 1] > Dom.xe)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] - Dom.xl;
      else if(bc_plane_pos[ 1] < Dom.xs)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] - Dom.xl;
    }
    if(bc_flow_configs[ 4] == PRECURSOR) {
      bc_plane_pos[ 2] = bc_plane_pos[ 2] - bc_flow_vels[ 4] * dt;
      if(bc_plane_pos[ 2] > Dom.xe)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] - Dom.xl;
      else if(bc_plane_pos[ 2] < Dom.xs)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] - Dom.xl;
    } else if(bc_flow_configs[ 5] == PRECURSOR) {
      bc_plane_pos[ 2] = bc_plane_pos[ 2] - bc_flow_vels[ 5] * dt;
      if(bc_plane_pos[ 2] > Dom.xe)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] - Dom.xl;
      else if(bc_plane_pos[ 2] < Dom.xs)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] - Dom.xl;
    }
    if(bc_flow_configs[ 6] == PRECURSOR) {
      bc_plane_pos[ 3] = bc_plane_pos[ 3] - bc_flow_vels[ 6] * dt;
      if(bc_plane_pos[ 3] > Dom.ye)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] - Dom.yl;
      else if(bc_plane_pos[ 3] < Dom.ys)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] - Dom.yl;
    } else if(bc_flow_configs[ 7] == PRECURSOR) {
      bc_plane_pos[ 3] = bc_plane_pos[ 3] - bc_flow_vels[ 7] * dt;
      if(bc_plane_pos[ 3] > Dom.ye)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] - Dom.yl;
      else if(bc_plane_pos[ 3] < Dom.ys)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] - Dom.yl;
    }
    if(bc_flow_configs[ 8] == PRECURSOR) {
      bc_plane_pos[ 4] = bc_plane_pos[ 4] - bc_flow_vels[ 8] * dt;
      if(bc_plane_pos[ 4] > Dom.ye)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] - Dom.yl;
      else if(bc_plane_pos[ 4] < Dom.ys)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] - Dom.yl;
    } else if(bc_flow_configs[ 9] == PRECURSOR) {
      bc_plane_pos[ 4] = bc_plane_pos[ 4] - bc_flow_vels[ 9] * dt;
      if(bc_plane_pos[ 4] > Dom.ye)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] - Dom.yl;
      else if(bc_plane_pos[ 4] < Dom.ys)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] - Dom.yl;
    }
    if(bc_flow_configs[10] == PRECURSOR) {
      bc_plane_pos[ 5] = bc_plane_pos[ 5] - bc_flow_vels[10] * dt;
      if(bc_plane_pos[ 5] > Dom.ye)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] - Dom.yl;
      else if(bc_plane_pos[ 5] < Dom.ys)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] - Dom.yl;
    } else if(bc_flow_configs[11] == PRECURSOR) {
      bc_plane_pos[ 5] = bc_plane_pos[ 5] - bc_flow_vels[11] * dt;
      if(bc_plane_pos[ 5] > Dom.ye)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] - Dom.yl;
      else if(bc_plane_pos[ 5] < Dom.ys)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] - Dom.yl;
    }
    if(bc_flow_configs[12] == PRECURSOR) {
      bc_plane_pos[ 6] = bc_plane_pos[ 6] - bc_flow_vels[12] * dt;
      if(bc_plane_pos[ 6] > Dom.ze)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] - Dom.zl;
      else if(bc_plane_pos[ 6] < Dom.zs)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] - Dom.zl;
    } else if(bc_flow_configs[13] == PRECURSOR) {
      bc_plane_pos[ 6] = bc_plane_pos[ 6] - bc_flow_vels[13] * dt;
      if(bc_plane_pos[ 6] > Dom.ze)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] - Dom.zl;
      else if(bc_plane_pos[ 6] < Dom.zs)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] - Dom.zl;
    }
    if(bc_flow_configs[14] == PRECURSOR) {
      bc_plane_pos[ 7] = bc_plane_pos[ 7] - bc_flow_vels[14] * dt;
      if(bc_plane_pos[ 7] > Dom.ze)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] - Dom.zl;
      else if(bc_plane_pos[ 7] < Dom.zs)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] - Dom.zl;
    } else if(bc_flow_configs[15] == PRECURSOR) {
      bc_plane_pos[ 7] = bc_plane_pos[ 7] - bc_flow_vels[15] * dt;
      if(bc_plane_pos[ 7] > Dom.ze)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] - Dom.zl;
      else if(bc_plane_pos[ 7] < Dom.zs)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] - Dom.zl;
    }
    if(bc_flow_configs[16] == PRECURSOR) {
      bc_plane_pos[ 8] = bc_plane_pos[ 8] - bc_flow_vels[16] * dt;
      if(bc_plane_pos[ 8] > Dom.ze)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] - Dom.zl;
      else if(bc_plane_pos[ 8] < Dom.zs)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] - Dom.zl;
    } else if(bc_flow_configs[17] == PRECURSOR) {
      bc_plane_pos[ 8] = bc_plane_pos[ 8] - bc_flow_vels[17] * dt;
      if(bc_plane_pos[ 8] > Dom.ze)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] - Dom.zl;
      else if(bc_plane_pos[ 8] < Dom.zs)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] - Dom.zl;
    }

    // pull the planes for the inflow
    cuda_dom_turb_planes_pull(bc_flow_configs);

    // send the planes to MASTER to be used as inflow velocities
    if(bc_flow_configs[ 0] == PRECURSOR || bc_flow_configs[ 1] == PRECURSOR)
      MPI_Send(u_WE, Dom.Gfx.jnb*Dom.Gfx.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 2] == PRECURSOR || bc_flow_configs[ 3] == PRECURSOR)
      MPI_Send(u_SN_S, Dom.Gfx.inb*Dom.Gfx.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
     MPI_Send(u_SN_N, Dom.Gfx.inb*Dom.Gfx.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 4] == PRECURSOR || bc_flow_configs[ 5] == PRECURSOR)
      MPI_Send(u_BT_B, Dom.Gfx.inb*Dom.Gfx.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(u_BT_T, Dom.Gfx.inb*Dom.Gfx.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 6] == PRECURSOR || bc_flow_configs[ 7] == PRECURSOR)
      MPI_Send(v_WE_W, Dom.Gfy.jnb*Dom.Gfy.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(v_WE_E, Dom.Gfy.jnb*Dom.Gfy.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 8] == PRECURSOR || bc_flow_configs[ 9] == PRECURSOR)
      MPI_Send(v_SN, Dom.Gfy.inb*Dom.Gfy.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[10] == PRECURSOR || bc_flow_configs[11] == PRECURSOR)
      MPI_Send(v_BT_B, Dom.Gfy.inb*Dom.Gfy.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(v_BT_T, Dom.Gfy.inb*Dom.Gfy.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[12] == PRECURSOR || bc_flow_configs[13] == PRECURSOR)
      MPI_Send(w_WE_W, Dom.Gfz.jnb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(w_WE_E, Dom.Gfz.jnb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[14] == PRECURSOR || bc_flow_configs[15] == PRECURSOR)
      MPI_Send(w_SN_S, Dom.Gfz.inb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(w_SN_N, Dom.Gfz.inb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[16] == PRECURSOR || bc_flow_configs[17] == PRECURSOR)
      MPI_Send(w_BT, Dom.Gfz.inb*Dom.Gfz.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
  }
}

void prec_send_BC(int np, int rank, MPI_Status status)
{
  if(np > 1) {
    // receive the boundary condition config info from MASTER
    MPI_Recv(&bc_flow_vels, 18, MPI_DOUBLE, MASTER, rank,
      MPI_COMM_WORLD, &status);

    // yank the appropriate planes for the inflow from the precursor domain
    // and put into their own plane arrays
    cuda_yank_turb_planes(bc_flow_configs, bc_plane_pos, bc_flow_vels);

    // pull the planes for the inflow
    cuda_dom_turb_planes_pull(bc_flow_configs);

    // send the planes to MASTER to be used as inflow velocities
    if(bc_flow_configs[ 0] == PRECURSOR || bc_flow_configs[ 1] == PRECURSOR)
      MPI_Send(u_WE, Dom.Gfx.jnb*Dom.Gfx.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 2] == PRECURSOR || bc_flow_configs[ 3] == PRECURSOR)
      MPI_Send(u_SN_S, Dom.Gfx.inb*Dom.Gfx.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(u_SN_N, Dom.Gfx.inb*Dom.Gfx.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 4] == PRECURSOR || bc_flow_configs[ 5] == PRECURSOR)
      MPI_Send(u_BT_B, Dom.Gfx.inb*Dom.Gfx.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(u_BT_T, Dom.Gfx.inb*Dom.Gfx.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 6] == PRECURSOR || bc_flow_configs[ 7] == PRECURSOR)
      MPI_Send(v_WE_W, Dom.Gfy.jnb*Dom.Gfy.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(v_WE_E, Dom.Gfy.jnb*Dom.Gfy.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[ 8] == PRECURSOR || bc_flow_configs[ 9] == PRECURSOR)
      MPI_Send(v_SN, Dom.Gfy.inb*Dom.Gfy.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[10] == PRECURSOR || bc_flow_configs[11] == PRECURSOR)
      MPI_Send(v_BT_B, Dom.Gfy.inb*Dom.Gfy.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(v_BT_T, Dom.Gfy.inb*Dom.Gfy.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[12] == PRECURSOR || bc_flow_configs[13] == PRECURSOR)
      MPI_Send(w_WE_W, Dom.Gfz.jnb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(w_WE_E, Dom.Gfz.jnb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[14] == PRECURSOR || bc_flow_configs[15] == PRECURSOR)
      MPI_Send(w_SN_S, Dom.Gfz.inb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
      MPI_Send(w_SN_N, Dom.Gfz.inb*Dom.Gfz.knb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);
    if(bc_flow_configs[16] == PRECURSOR || bc_flow_configs[17] == PRECURSOR)
      MPI_Send(w_BT, Dom.Gfz.inb*Dom.Gfz.jnb, MPI_DOUBLE, MASTER, rank,
        MPI_COMM_WORLD);

    // update plane position
    if(bc_flow_configs[ 0] == PRECURSOR) {
      bc_plane_pos[ 0] = bc_plane_pos[ 0] - bc_flow_vels[ 0] * dt;
      if(bc_plane_pos[ 0] > Dom.xe)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] - Dom.xl;
      else if(bc_plane_pos[ 0] < Dom.xs)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] + Dom.xl;
    } else if(bc_flow_configs[ 1] == PRECURSOR) {
      bc_plane_pos[ 0] = bc_plane_pos[ 0] - bc_flow_vels[ 1] * dt;
      if(bc_plane_pos[ 0] > Dom.xe)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] - Dom.xl;
      else if(bc_plane_pos[ 0] < Dom.xs)
        bc_plane_pos[ 0] = bc_plane_pos[ 0] + Dom.xl;
    }
    if(bc_flow_configs[ 2] == PRECURSOR) {
      bc_plane_pos[ 1] = bc_plane_pos[ 1] - bc_flow_vels[ 2] * dt;
      if(bc_plane_pos[ 1] > Dom.xe)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] - Dom.xl;
      else if(bc_plane_pos[ 1] < Dom.xs)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] + Dom.xl;
    } else if(bc_flow_configs[ 3] == PRECURSOR) {
      bc_plane_pos[ 1] = bc_plane_pos[ 1] - bc_flow_vels[ 3] * dt;
      if(bc_plane_pos[ 1] > Dom.xe)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] - Dom.xl;
      else if(bc_plane_pos[ 1] < Dom.xs)
        bc_plane_pos[ 1] = bc_plane_pos[ 1] + Dom.xl;
    }
    if(bc_flow_configs[ 4] == PRECURSOR) {
      bc_plane_pos[ 2] = bc_plane_pos[ 2] - bc_flow_vels[ 4] * dt;
      if(bc_plane_pos[ 2] > Dom.xe)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] - Dom.xl;
      else if(bc_plane_pos[ 2] < Dom.xs)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] + Dom.xl;
    } else if(bc_flow_configs[ 5] == PRECURSOR) {
      bc_plane_pos[ 2] = bc_plane_pos[ 2] - bc_flow_vels[ 5] * dt;
      if(bc_plane_pos[ 2] > Dom.xe)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] - Dom.xl;
      else if(bc_plane_pos[ 2] < Dom.xs)
        bc_plane_pos[ 2] = bc_plane_pos[ 2] + Dom.xl;
    }
    if(bc_flow_configs[ 6] == PRECURSOR) {
      bc_plane_pos[ 3] = bc_plane_pos[ 3] - bc_flow_vels[ 6] * dt;
      if(bc_plane_pos[ 3] > Dom.ye)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] - Dom.yl;
      else if(bc_plane_pos[ 3] < Dom.ys)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] + Dom.yl;
    } else if(bc_flow_configs[ 7] == PRECURSOR) {
      bc_plane_pos[ 3] = bc_plane_pos[ 3] - bc_flow_vels[ 7] * dt;
      if(bc_plane_pos[ 3] > Dom.ye)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] - Dom.yl;
      else if(bc_plane_pos[ 3] < Dom.ys)
        bc_plane_pos[ 3] = bc_plane_pos[ 3] + Dom.yl;
    }
    if(bc_flow_configs[ 8] == PRECURSOR) {
      bc_plane_pos[ 4] = bc_plane_pos[ 4] - bc_flow_vels[ 8] * dt;
      if(bc_plane_pos[ 4] > Dom.ye)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] - Dom.yl;
      else if(bc_plane_pos[ 4] < Dom.ys)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] + Dom.yl;
    } else if(bc_flow_configs[ 9] == PRECURSOR) {
      bc_plane_pos[ 4] = bc_plane_pos[ 4] - bc_flow_vels[ 9] * dt;
      if(bc_plane_pos[ 4] > Dom.ye)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] - Dom.yl;
      else if(bc_plane_pos[ 4] < Dom.ys)
        bc_plane_pos[ 4] = bc_plane_pos[ 4] + Dom.yl;
    }
    if(bc_flow_configs[10] == PRECURSOR) {
      bc_plane_pos[ 5] = bc_plane_pos[ 5] - bc_flow_vels[10] * dt;
      if(bc_plane_pos[ 5] > Dom.ye)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] - Dom.yl;
      else if(bc_plane_pos[ 5] < Dom.ys)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] + Dom.yl;
    } else if(bc_flow_configs[11] == PRECURSOR) {
      bc_plane_pos[ 5] = bc_plane_pos[ 5] - bc_flow_vels[11] * dt;
      if(bc_plane_pos[ 5] > Dom.ye)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] - Dom.yl;
      else if(bc_plane_pos[ 5] < Dom.ys)
        bc_plane_pos[ 5] = bc_plane_pos[ 5] + Dom.yl;
    }
    if(bc_flow_configs[12] == PRECURSOR) {
      bc_plane_pos[ 6] = bc_plane_pos[ 6] - bc_flow_vels[12] * dt;
      if(bc_plane_pos[ 6] > Dom.ze)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] - Dom.zl;
      else if(bc_plane_pos[ 6] < Dom.zs)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] + Dom.zl;
    } else if(bc_flow_configs[13] == PRECURSOR) {
      bc_plane_pos[ 6] = bc_plane_pos[ 6] - bc_flow_vels[13] * dt;
      if(bc_plane_pos[ 6] > Dom.ze)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] - Dom.zl;
      else if(bc_plane_pos[ 6] < Dom.zs)
        bc_plane_pos[ 6] = bc_plane_pos[ 6] + Dom.zl;
    }
    if(bc_flow_configs[14] == PRECURSOR) {
      bc_plane_pos[ 7] = bc_plane_pos[ 7] - bc_flow_vels[14] * dt;
      if(bc_plane_pos[ 7] > Dom.ze)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] - Dom.zl;
      else if(bc_plane_pos[ 7] < Dom.zs)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] + Dom.zl;
    } else if(bc_flow_configs[15] == PRECURSOR) {
      bc_plane_pos[ 7] = bc_plane_pos[ 7] - bc_flow_vels[15] * dt;
      if(bc_plane_pos[ 7] > Dom.ze)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] - Dom.zl;
      else if(bc_plane_pos[ 7] < Dom.zs)
        bc_plane_pos[ 7] = bc_plane_pos[ 7] + Dom.zl;
    }
    if(bc_flow_configs[16] == PRECURSOR) {
      bc_plane_pos[ 8] = bc_plane_pos[ 8] - bc_flow_vels[16] * dt;
      if(bc_plane_pos[ 8] > Dom.ze)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] - Dom.zl;
      else if(bc_plane_pos[ 8] < Dom.zs)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] + Dom.zl;
    } else if(bc_flow_configs[17] == PRECURSOR) {
      bc_plane_pos[ 8] = bc_plane_pos[ 8] - bc_flow_vels[17] * dt;
      if(bc_plane_pos[ 8] > Dom.ze)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] - Dom.zl;
      else if(bc_plane_pos[ 8] < Dom.zs)
        bc_plane_pos[ 8] = bc_plane_pos[ 8] + Dom.zl;
    }

    // send this timestep size to MASTER for comparison
    MPI_Send(&dt, 1, MPI_DOUBLE, MASTER, rank, MPI_COMM_WORLD);
    // wait for a response holding the smallest of the dts
    MPI_Recv(&dt, 1, MPI_DOUBLE, MASTER, rank, MPI_COMM_WORLD, &status);
  }
}

void prec_compare_dt(int np, int rank, MPI_Status status)
{
  if(np > 1) {
    // send this timestep size to MASTER for initialization
    MPI_Send(&dt, 1, MPI_DOUBLE, MASTER, rank, MPI_COMM_WORLD);
    // receive the boundary condition config info from MASTER
    MPI_Recv(&dt, 1, MPI_DOUBLE, MASTER, rank, MPI_COMM_WORLD, &status);
  }
}

void expd_comm_restart_write(int np, int rest_com) {
  if(np > 1)
    MPI_Send(&rest_com, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
}

void prec_comm_restart_write(int np, int *rest_com, int rank,
  MPI_Status status) {
  if(np > 1)
    MPI_Recv(rest_com, 1, MPI_INT, MASTER, rank, MPI_COMM_WORLD, &status);
}
