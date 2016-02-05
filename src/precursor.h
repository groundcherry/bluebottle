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

/****h* Bluebottle/precursor
 * NAME
 *  precursor
 * FUNCTION
 *  A utility for performing the various data manipulations needed for
 *  transferring precursor flow data to the inflow plane of the experimental
 *  domain.
 ******
 */

#ifndef _PRECURSOR_H
#define _PRECURSOR_H

#include <mpi.h>
#include "bluebottle.h"

/****f* precursor/expd_init_BC()
 * NAME
 *  expd_init_BC()
 * TYPE
 */
void expd_init_BC(int np);
/*
 * FUNCTION
 *  Communiciate the boundary condition configuration information with the
 *  precursor.
 * ARGUMENTS
 *  * np -- number of processors
 ******
 */

/****f* precursor/expd_update_BC()
 * NAME
 *  expd_update_BC()
 * TYPE
 */
void expd_update_BC(int np, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition configuration information with the
 *  precursor.
 * ARGUMENTS
 *  * np -- number of processors
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/expd_compare_dt()
 * NAME
 *  expd_compare_dt()
 * TYPE
 */
void expd_compare_dt(int np, MPI_Status status);
/*
 * FUNCTION
 *  Negotiate with the precursor domain on the appropriate timestep size to
 *  take.
 * ARGUMENTS
 *  * np -- number of processors
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/prec_init_BC()
 * NAME
 *  prec_init_BC()
 * TYPE
 */
void prec_init_BC(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition configuration information with the
 *  experimental domain.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/prec_update_BC()
 * NAME
 *  prec_update_BC()
 * TYPE
 */
void prec_update_BC(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition configuration information with the
 *  experimental domain.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 ******
 */

/****f* precursor/prec_send_BC()
 * NAME
 *  prec_send_BC()
 * TYPE
 */
void prec_send_BC(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition to the experimental domain.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/prec_compare_dt()
 * NAME
 *  prec_compare_dt()
 * TYPE
 */
void prec_compare_dt(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Negotiate with the precursor domain on the appropriate timestep size to
 *  take.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/expd_comm_restart_write()
 * NAME
 *  expd_comm_restart_write()
 * TYPE
 */
void expd_comm_restart_write(int np, int rest_com);
/*
 * FUNCTION
 *  Send a flag to write a restart file to the turbulent precursor.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rest_com -- boolean-valued restart flag
 ******
 */

/****f* precursor/prec_comm_restart_write()
 * NAME
 *  prec_comm_restart_write()
 * TYPE
 */
void prec_comm_restart_write(int np, int *rest_com, int rank,
  MPI_Status status);
/*
 * FUNCTION
 *  Receive a flag to write a restart file to the turbulent precursor.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rest_com -- boolean-valued restart flag
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

#endif
