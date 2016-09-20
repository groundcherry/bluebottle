################################################################################
################################## BLUEBOTTLE ##################################
################################################################################
#
#   Copyright 2012 - 2016 Adam Sierakowski, The Johns Hopkins University
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#       http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
# 
#   Please contact the Johns Hopkins University to use Bluebottle for
#   commercial and/or for-profit applications.
################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EDIT: DEPENDENCIES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MPI_DIR = /usr/lib/openmpi
HDF5_DIR = /usr/local/hdf5
CGNS_DIR = /usr/local/cgns
CUSP_DIR = /usr/local/cusplibrary
CUDA_DIR = /usr/local/cuda
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ EDIT: COMPILERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
MPICC = mpicc.openmpi
NVCC = nvcc
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

PREC = DOUBLE

SRC_DIR = src
SIM_DIR = sim

COPT = -std=c99 -pedantic -Wall -Wextra -fopenmp -D$(PREC)
LDINCS = -I $(MPI_DIR)/include -I $(CGNS_DIR)/include -I $(CUSP_DIR)
LDLIBS = -lm -L $(HDF5_DIR)/lib -L $(CGNS_DIR)/lib -lcgns -lhdf5 

CUDAOPT = -arch=sm_30 -Xcompiler -fopenmp -m64 -D$(PREC)

CUDALIBS = -L $(CUDA_DIR)/lib64 -lcudart

SRCC =	bluebottle.c	\
	domain.c	\
	particle.c	\
	precursor.c	\
	recorder.c	\
	seeder.c	\
	vtk.c		\
	rng.c

SRCCUDA = cuda_bluebottle.cu	\
	cuda_bicgstab.cu	\
	cuda_particle.cu	\
	cuda_quadrature.cu	\
	cuda_testing.cu		\
	entrySearch.cu		\
	bluebottle_kernel.cu	\
	bicgstab_kernel.cu	\
	entrySearch_kernel.cu	\
	particle_kernel.cu	\
	quadrature_kernel.cu

EXTRA = Makefile		\
	bluebottle.h		\
	cuda_bluebottle.h	\
	cuda_bicgstab.h		\
	cuda_particle.h		\
	cuda_quadrature.h	\
	cuda_testing.h		\
	domain.h		\
	entrySearch.h		\
	particle.h		\
	precursor.h		\
	recorder.h		\
	vtk.h			\
	rng.h

# compile normally
all: COPT += -O2
all: CUDAOPT += -O2
all: bluebottle

# compile with explicit numerics
implicit: COPT += -O2 -DIMPLICIT
implicit: CUDAOPT += -O2 -DIMPLICIT
implicit: bluebottle

# compile with Stokes flow only
stokes: COPT += -O2 -DSTOKESFLOW
stokes: CUDAOPT += -O2 -DSTOKESFLOW
stokes: bluebottle

# compile without pressure correction for U*
nocorrection: COPT += -O2 -DNOPRESSURECORRECTION
nocorrection: CUDAOPT += -O2 -DNOPRESSURECORRECTION
nocorrection: bluebottle

# compile for batch job submission
batch: COPT += -O2 -DBATCHRUN
batch: CUDAOPT += -O2
batch: bluebottle

# compile with stair-stepped interior boundaries
steps: COPT += -DSTEPS -O2
steps: CUDAOPT += -DSTEPS -O2
steps: bluebottle

# compile with debug output
debug: COPT += -DDDEBUG -g
debug: CUDAOPT += -DDDEBUG -g -G
debug: bluebottle

# compile with testing code
test: COPT += -DDDEBUG -DTEST -g
test: CUDAOPT += -DDDEBUG -DTEST -g -G
test: bluebottle

# write robodoc documentation
doc:
	cd .. && robodoc --html --multidoc --doc doc/robodoc && robodoc --latex --singledoc --sections --doc doc/LaTeX/Bluebottle_0.1_robodoc && cd doc/LaTeX && pdflatex Bluebottle_0.1_robodoc.tex && pdflatex Bluebottle_0.1_robodoc.tex && pdflatex Bluebottle_0.1_robodoc.tex && echo '\nmake doc: Complete.'

OBJS = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCC))))
OBJSCUDA = $(addprefix $(SRC_DIR)/, $(addsuffix .o, $(basename $(SRCCUDA))))

$(OBJSCUDA):$(SRC_DIR)/%.o:$(SRC_DIR)/%.cu
	$(NVCC) $(CUDAOPT) -dc $< $(LDINCS) -o $@

$(OBJS):$(SRC_DIR)/%.o:$(SRC_DIR)/%.c
	$(MPICC) $(COPT) -c $< $(LDINCS) -o $@

$(SRC_DIR)/bblib.o:$(OBJSCUDA)
	$(NVCC) $(CUDAOPT) -dlink $+ $(CUDALIBS) -o $(SRC_DIR)/bblib.o

bluebottle: $(OBJSCUDA) $(SRC_DIR)/bblib.o $(OBJS)
	$(MPICC) $(COPT) $+ $(LDLIBS) -lcudadevrt -lstdc++ $(CUDALIBS) -o $(SIM_DIR)/bluebottle

clean:
	rm -f $(SRC_DIR)/*.o $(SIM_DIR)/bluebottle
