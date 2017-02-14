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

/****h* Bluebottle/vtk
 * NAME
 *  vtk
 * FUNCTION
 *  Write Paraview VTK files.
 ******
 */

#ifndef _VTK_H
#define _VTK_H

/****f* vtk/init_VTK()
 * NAME
 *  init_VTK()
 * USAGE
 */
void init_VTK(void);
/*
 * FUNCTION
 *  Set up the Paraview VTK output PVD file for writing timestep information
 *  about each output file.
 ******
 */

/****f* vtk/out_VTK()
 * NAME
 *  out_VTK()
 * USAGE
 */
void out_VTK(void);
/*
 * FUNCTION
 *  Output ParaView VTK files for the timestep (flow domain and particle data).
 ******
 */

/****f* vtk/dom_out_VTK()
 * NAME
 *  dom_out_VTK()
 * USAGE
 */
void dom_out_VTK(void);
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing velocity and pressure fields.  The
 *  velocity field is interpolated from the faces and assigned to a cell.
 ******
 */

/****f* vtk/init_VTK_ghost()
 * NAME
 *  init_VTK_ghost()
 * USAGE
 */
void init_VTK_ghost(void);
/*
 * FUNCTION
 *  Set up the Paraview VTK output PVD file for writing timestep information
 *  about each output file. This version outputs ghost cells.
 ******
 */

/****f* vtk/out_VTK_ghost()
 * NAME
 *  out_VTK_ghost()
 * USAGE
 */
void out_VTK_ghost(void);
/*
 * FUNCTION
 *  Output ParaView VTK files for the timestep (flow domain and particle data).
 *  This version outputs ghost cells.
 ******
 */

/****f* vtk/dom_out_VTK_ghost()
 * NAME
 *  dom_out_VTK_ghost()
 * USAGE
 */
void dom_out_VTK_ghost(void);
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing velocity and pressure fields.  The
 *  velocity field is interpolated from the faces and assigned to a cell. This
 *  version outputs ghost cells.
 ******
 */

/****f* vtk/init_VTK_turb()
 * NAME
 *  init_VTK_turb()
 * USAGE
 */
void init_VTK_turb(void);
/*
 * FUNCTION
 *  Set up the Paraview VTK output PVD file for writing timestep information
 *  about each output file.
 ******
 */

/****f* vtk/out_VTK_turb()
 * NAME
 *  out_VTK_turb()
 * USAGE
 */
void out_VTK_turb(void);
/*
 * FUNCTION
 *  Output ParaView VTK files for the timestep (flow domain and particle data).
 ******
 */

/****f* vtk/dom_out_VTK_turb()
 * NAME
 *  dom_out_VTK_turb()
 * USAGE
 */
void dom_out_VTK_turb(void);
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing velocity and pressure fields.  The
 *  velocity field is interpolated from the faces and assigned to a cell.
 ******
 */

/****f* vtk/part_out_VTK()
 * NAME
 *  part_out_VTK()
 * USAGE
 */
void part_out_VTK();
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing particle data.
 ******
 */

/****f* vtk/quadnodes_out_VTK()
 * NAME
 *  quadnodes_out_VTK()
 * USAGE
 */
void quadnodes_out_VTK();
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing particle quadrature node data.
 ******
 */

#endif
