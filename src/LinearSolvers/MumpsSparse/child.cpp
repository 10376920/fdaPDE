#include <mpi.h>
#include "dmumps_c.h"
#include <iostream>
int main (int argc,char** argv){
	int rank;
	int size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm parent;
	MPI_Comm_get_parent(&parent);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm intracomm;
	MPI_Intercomm_merge(parent, 1, &intracomm);
	DMUMPS_STRUC_C id;
	id.job = -1;
	id.par=1;
	id.comm_fortran= (MUMPS_INT) MPI_Comm_c2f(intracomm);
	dmumps_c(&id);
	do{
		MPI_Bcast(&id.job, 1, MPI_INT, 0, parent);
		id.sym=0;
		dmumps_c(&id);
	} while(id.job != -2);
	
	MPI_Finalize ();
}
