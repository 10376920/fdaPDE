#include "MumpsSparse.h"
#include <limits>
#include <cmath>

namespace LinearSolvers {

MumpsSparse::MumpsSparse() {
	setDefault();
}

MumpsSparse::MumpsSparse(const ParameterList &list) {
	setDefault();
	setParameters(list);
}

void MumpsSparse::setDefault() {
	_parallel_flag = false;
	_children_is_empty = true;
	_id.sym = 0;
	_id.par = 1;
	_id.comm_fortran = _use_comm_world;
	_nproc = 1;
	_id.job = _job_init;
	dmumps_c(&_id);
	//qui
	//_id.icntl[13]=100;
}

MumpsSparse::~MumpsSparse(){
	//andare out of scope perchè funzioni
	// ricordarsi di allineare i broadcast con i broadcast di child
	_id.job = _job_end;
	if (_parallel_flag == true) {
		MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, _children);
	}
	dmumps_c(&_id);
}

void MumpsSparse::factorize(const Eigen::SparseMatrix<double> &A) {
	std::vector<int> irn;
	std::vector<int> jcn;
	std::vector<double> a;
	_id.n = A.cols();
	bool is_symmetric = (_id.sym == 1 || _id.sym == 2);
	for (int j=0; j<A.outerSize(); ++j){
		for (Eigen::SparseMatrix<double>::InnerIterator it(A,j); it; ++it){
		// If the matrix is symmetric we store only the lower triangular part.
			if (!is_symmetric || (is_symmetric && it.row() <= it.col())) {
				irn.push_back(it.row()+1);
				jcn.push_back(it.col()+1);
				a.push_back(it.value());
			}
		}
	}
	_id.nz = irn.size();//A.nonZeros();
	_id.irn=irn.data();
	_id.jcn=jcn.data();
	_id.a=a.data();
	_id.job = _job_analyze_factorize;
	//fare broadcast per coere_nz tra i processi
	//MPI_Bcast(_id.icntl, 40, MPI_INT, MPI_ROOT, _children);
	if (_parallel_flag == true) {
		MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, _children);
	}
	dmumps_c(&_id);
}

void MumpsSparse::solve(const Eigen::MatrixXd &b) {
	int n = b.cols()*b.rows();
	double *rhs = new double[n];
	for(int i=0; i<n; ++i) {
		rhs[i] = b(i);
	}
	_id.rhs = rhs;
	_id.nrhs = b.cols();
	_id.lrhs = b.rows();
	_id.job = _job_solve;
	//fare broadcast per coere_nz tra i processi
	//MPI_Bcast(_id.icntl, 40, MPI_INT, MPI_ROOT, _children);
	if (_parallel_flag == true) {
		MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, _children);
	}
	dmumps_c(&_id);

	_sol.resize(b.rows(), b.cols());
	for(int i=0; i<n; ++i) {
		_sol(i) = rhs[i];
	}
	delete[] rhs;
}

void MumpsSparse::setParameters(const ParameterList &list) {
	//occhio che incntl devono esssere sfasati tra umano e pc
	bool changed_sym   = list.getValue<int>("sym", 0)   != _id.sym;
	bool changed_par   = list.getValue<int>("par", 1)   != _id.par;
	bool changed_nproc = list.getValue<int>("nproc", 1) != _nproc;
	// se sono cambiati rispetto al giro prima 
	//quindi se non ci sono nella lista vogliamo che ritornino al valore di default  
	if (changed_sym || changed_par || changed_nproc) {
		_id.job = _job_end;
		if (! _children_is_empty) {
			MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, _children);
		}
		dmumps_c(&_id);
		int par = list.getValue<int>("par", 1);
		_nproc  = list.getValue<int>("nproc", 1);
	
		if (_nproc > 1) {
			_parallel_flag = true;
			MPI_Comm_spawn(CHILD_PATH, NULL, _nproc - par, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &_children, err);
			MPI_Comm intracomm;
			MPI_Intercomm_merge(_children, 0, &intracomm);
			_children_is_empty = false;
			_id.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(intracomm);
		}
		else {
			_parallel_flag = false;
			_id.par = 1;
			_id.comm_fortran = _use_comm_world;
			if (! _children_is_empty) {
				MPI_Comm_free(&_children);
				_children_is_empty = true;
			}
		}
		_id.par = list.getValue<int>("par", 1);
		_id.sym = list.getValue<int>("sym", 0);
		_id.job = _job_init;
		dmumps_c(&_id);
	}
	
	_id.icntl[0]  = list.getValue<int>("icntl[1]", 6);
	_id.icntl[1]  = list.getValue<int>("icntl[2]", 0);
	_id.icntl[2]  = list.getValue<int>("icntl[3]", 6);
	_id.icntl[3]  = list.getValue<int>("icntl[4]", 2);
	_id.icntl[6]  = list.getValue<int>("icntl[7]", 7);
	_id.icntl[9]  = list.getValue<int>("icntl[10]", 0);
	_id.icntl[10] = list.getValue<int>("icntl[11]", 0);
	_id.icntl[12] = list.getValue<int>("icntl[13]", 0);
	_id.icntl[13] = list.getValue<int>("icntl[14]", 20);
	_id.icntl[22] = list.getValue<int>("icntl[23]", 0);
	_id.icntl[26] = list.getValue<int>("icntl[27]", -8);
	_id.icntl[27] = list.getValue<int>("icntl[28]", 0);
	_id.icntl[30] = list.getValue<int>("icntl[31]", 0);
	_id.icntl[31] = list.getValue<int>("icntl[32]", 0);

	
	// if (_id.sym == 1 || _id.sym == 2) {
	// 	_id.cntl[0] = list.getValue<double>("cntl[1]", 0.01);
	// }
	// else {
	// 	_id.cntl[0] = list.getValue<double>("cntl[1]", 0);
	// }
	// double cntl1_default = sqrt(std::numeric_limits<double>::epsilon());
	// _id.cntl[1] = list.getValue<double>("cntl[2]", cntl1_default);
	// _id.cntl[3] = list.getValue<double>("cntl[4]", -1.0);
}

} // End Namespace LinearSolvers

