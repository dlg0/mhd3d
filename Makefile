
ifeq (${MACHINE},dlghp)
	HOME = /home/dg6/code
	FC	= ${HOME}/intel/ifort/bin/ifort 
	IMSL_DIR = ${HOME}/imsl/imsl_lib/imsl/fnl600/lnxin100i32
	NETCDF = -I ${HOME}/netcdf/netcdf_intel/include -L ${HOME}/netcdf/netcdf_intel/lib -lnetcdf -lnetcdf_c++
	DISLIN = -I${HOME}/dislin/ifc -L${HOME}/dislin -ldislin
endif

ifeq (${MACHINE},fitacf)
	HOME = /home/superdarn
	FC	= ${HOME}/intel/Compiler/11.0/074/bin/intel64/ifort 
	IMSL_DIR = ${HOME}/imsl/imsl/fnl600/lnxin100e64
	NETCDF = -I ${HOME}/netcdf/netcdf_intel/include -L ${HOME}/netcdf/netcdf_intel/lib -lnetcdf -lnetcdf_c++
	DISLIN = -I${HOME}/dislin/ifc -L${HOME}/dislin -ldislin
endif
	
MOD = mod
OBJ = obj
SRC = src

FFLAGS	= -openmp -p -g -warn all -check bounds -module ${MOD}
IMSL = -L${IMSL_DIR}/lib -Bdynamic -limsl -limslsuperlu -limslscalar -limslblas -limslmpistub -lm -Xlinker -rpath -Xlinker ${IMSL_DIR}/lib -I ${IMSL_DIR}/include
EXEC = xMhd3d


OBJECTS = ${OBJ}/newtonsRule.o \
					${OBJ}/mhd_grid.o \
					${OBJ}/timer.o \
					${OBJ}/dlg.o \
					${OBJ}/vA_profile.o \
					${OBJ}/spherHarmFns.o

mhd3d : ${OBJECTS} ${SRC}/mhd3d.f90 
	${FC} ${FFLAGS} ${IMSL} ${DISLIN} ${SRC}/mhd3d.f90 ${OBJECTS} -o ${EXEC} ${NETCDF}

${OBJ}/spherHarmFns.o : ${SRC}/spherHarmFns.f90 ${OBJ}/mhd_grid.o ${MOD}/mhd_grid.mod
	${FC} ${FFLAGS} ${IMSL} ${DISLIN} -c ${SRC}/spherHarmFns.f90 -o ${OBJ}/spherHarmFns.o

${OBJ}/vA_profile.o : ${SRC}/vA_profile.f90 ${OBJ}/mhd_grid.o ${MOD}/mhd_grid.mod
	${FC} ${FFLAGS} -c ${SRC}/vA_profile.f90 -o ${OBJ}/vA_profile.o

${OBJ}/mhd_grid.o : ${SRC}/mhd_grid.f90
	${FC} ${FFLAGS} -c ${SRC}/mhd_grid.f90 -o ${OBJ}/mhd_grid.o

${OBJ}/newtonsRule.o : ${SRC}/newtonsRule.f90
	${FC} ${FFLAGS} -c ${SRC}/newtonsRule.f90 -o ${OBJ}/newtonsRule.o

${OBJ}/timer.o : ${SRC}/timer.f90
	${FC} ${FFLAGS} -c ${SRC}/timer.f90 -o ${OBJ}/timer.o

${OBJ}/dlg.o : ${SRC}/dlg.f90
	${FC} ${FFLAGS} -c ${SRC}/dlg.f90 ${NETCDF} -o ${OBJ}/dlg.o

clean:
	rm ${EXEC} ${MOD}/*.mod ${OBJ}/*.o
