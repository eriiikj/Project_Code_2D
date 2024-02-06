#
# Only compile the modulus
#

CF77=ifort
CF90=ifort
F77FLAGS = -03
#F77FLAGS= -O3  -ftrapuv -check all -check bound -inline_debug_info -warn all -traceback -heap-arrays  -heap-arrays -fstack-protector-all -fstack-security-check
#F77FLAGS= -O3   -warn all -traceback -fstack-security-check -gen-interfaces -warn interfaces -check arg_temp_created
#F77FLAGS= -O3 
F90FLAGS = ${F77FLAGS}

# modules
# note order is important if module used in other module
SolidLib=some_constants.o matrix_util.o fem_util.o memory_util.o \
	sparse_util.o fem_system.o \
	mater_small.o mater_large.o mater_hyperel.o  \
	mater_magnlin.o \
	abaqus_util.o readind_util.o\
	elem_small_cont_2d.o elem_large_cont_2d.o \
	elem_large_cont_3d.o  \
	elem_magn_2d.o \
	mater_J2iso.o  \
	wrt2vtk.o wrt2silo.o elem_flow_2d.o mater_mises.o\
	mater_magnjav.o mater_magnnonl1.o mater_magnkho.o mater_magn.o matlab_util.o


.SUFFIXES:
.SUFFIXES: .f .f90 .mod .o .x

# making objects
.f90.o:
	${CF90} ${F90FLAGS} -c $*.f90

.f.o:
	${CF77} ${F77FLAGS} -c $*.f

# making module
.f90.mod:
	${CF90} ${F90FLAGS} -c $*.f90

default: solid_mod.o

solid_mod.o: ${SolidLib}

clean:
	/bin/rm -f *.x *.o *.mod
