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

/****h* Bluebottle/recorder
 * NAME
 *  recorder
 * FUNCTION
 *  A utility for recording simulation metrics.
 ******
 */

#ifndef _RECORDER_H
#define _RECORDER_H

#include "bluebottle.h"
#include "particle.h"
#include <cgnslib.h>

/****f* recorder/recorder_read_config()
 * NAME
 *  recorder_read_config()
 * TYPE
 */
void recorder_read_config(void);
/*
 * FUNCTION
 *  Read the record.config file to determine what output to write.
 ******
 */

/****f* recorder/cgns_grid()
 * NAME
 *  cgns_grid()
 * TYPE
 */
void cgns_grid(void);
/*
 * FUNCTION
 *  Write the CGNS grid output file.
 ******
 */

/****f* recorder/cgns_turb_grid()
 * NAME
 *  cgns_turb_grid()
 * TYPE
 */
void cgns_turb_grid(void);
/*
 * FUNCTION
 *  Write the CGNS grid output file.
 ******
 */

/****f* recorder/cgns_flow_field()
 * NAME
 *  cgns_flow_field()
 * TYPE
 */
void cgns_flow_field(real dtout);
/*
 * FUNCTION
 *  Write the CGNS flow_field output file.
 * ARGUMENTS
 *  * dtout -- the output timestep size
 ******
 */

/****f* recorder/cgns_turb_flow_field()
 * NAME
 *  cgns_turb_flow_field()
 * TYPE
 */
void cgns_turb_flow_field(real dtout);
/*
 * FUNCTION
 *  Write the CGNS flow_field output file.
 * ARGUMENTS
 *  * dtout -- the output timestep size
 ******
 */

/****f* recorder/cgns_particles()
 * NAME
 *  cgns_particles()
 * TYPE
 */
void cgns_particles(real dtout);
/*
 * FUNCTION
 *  Write the CGNS particles output file.
 * ARGUMENTS
 *  * dtout -- the output timestep size
 ******
 */

/****f* recorder/recorder_bicgstab_init()
 * NAME
 *  recorder_bicgstab_init()
 * TYPE
 */
void recorder_bicgstab_init(char *name);
/*
 * FUNCTION
 *  Create the file name for writing and summarize fields to be written for the
 *  BICGSTAB solver.
 * ARGUMENTS
 *  * name -- the name of the file to be written
 ******
 */

/****f* recorder/recorder_lamb_init()
 * NAME
 *  recorder_lamb_init()
 * TYPE
 */
void recorder_lamb_init(char *name);
/*
 * FUNCTION
 *  Create the file name for writing Lamb's coefficients.
 * ARGUMENTS
 *  * name -- the name of the file to be written
 ******
 */

/****f* recorder/recorder_bicgstab()
 * NAME
 *  recorder_bicgstab()
 * TYPE
 */
void recorder_bicgstab(char *name, int niter, real resid);
/*
 * FUNCTION 
 *  Write out BICGSTAB solver information to file name.
 * ARGUMENTS
 *  * name -- the name of the file to which to write
 *  * niter -- the number of iterations to convergence
 *  * resid -- the residual at convergence
 ******
 */

/****f* recorder/recorder_lamb()
 * NAME
 *  recorder_lamb()
 * TYPE
 */
void recorder_lamb(char *name, int iter);
/*
 * FUNCTION 
 *  Write out Lamb's coefficients to file.
 * ARGUMENTS
 *  * name -- the name of the file to which to write
 *  * iter -- the number of the iteration
 ******
 */

#endif
