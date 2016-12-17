CC=mpicc
CXX=mpic++

DIR_MUMPS_SPARSE=LinearSolvers/MumpsSparse
DIR_EIGEN_SPARSE_LU=LinearSolvers/EigenSparseLU

SRCS_BASIS=${wildcard *.cpp *.c}
SRCS_MUMPS_SPARSE=${DIR_MUMPS_SPARSE}/MumpsSparse.cpp
SRCS_EIGEN_SPARSE_LU=${DIR_EIGEN_SPARSE_LU}/EigenSparseLU.cpp
SRCS_SHLIB=${SRCS_BASIS} ${SRCS_MUMPS_SPARSE} ${SRCS_EIGEN_SPARSE_LU}

all: library child

library:
	MAKEFLAGS="CC=${CC} CXX=${CXX}" &&\
	R CMD SHLIB ${SRCS_SHLIB} -o fdaPDE.so

child:
	$(CXX) -std=c++11 -ldmumps ${DIR_MUMPS_SPARSE}/child.cpp -o child

clean:
	rm -f *.o ${DIR_MUMPS_SPARSE}/*.o ${DIR_EIGEN_SPARSE_LU}/*.o
