#include "MumpsSparse.h"

namespace LinearSolvers {

	MumpsSparse::MumpsSparse() {
		int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &_myid);
		_id.sym = 0;
		_id.par = 1;
		_id.comm_fortran = _use_comm_world;
		_id.job = _job_init;
		_id.icntl[0]=-1;
		_id.icntl[1]=-1;
		_id.icntl[2]=-1;
		_id.icntl[3]=4;
		dmumps_c(&_id);
	}

	void MumpsSparse::factorize(const Eigen::SparseMatrix<double> &A) {
		std::vector<int> irn;
		std::vector<int> jcn;
		std::vector<double> a;
		if(_myid == 0){
			_id.n = A.cols();
			_id.nz = A.nonZeros();
			for (int j=0; j<A.outerSize(); ++j){
				for (Eigen::SparseMatrix<double>::InnerIterator it(A,j); it; ++it){
					irn.push_back(it.row()+1);
					jcn.push_back(it.col()+1);
					a.push_back(it.value());
				}
			}
			_id.irn=irn.data();
			_id.jcn=jcn.data();
			_id.a=a.data();
		}
		_id.job = _job_analyze_factorize;
		dmumps_c(&_id);
	}

	void MumpsSparse::solve(const Eigen::MatrixXd &b) {
		int n = b.cols()*b.rows();
		double *rhs = new double[n];
		if (_myid == 0) {
			for(int i=0; i<n; ++i) {
				rhs[i] = b(i);
			}
			_id.rhs = rhs;
			_id.nrhs = b.cols();
			_id.lrhs = b.rows();
		}
		_id.job = _job_solve;
		dmumps_c(&_id);
		_sol.resize(b.rows(), b.cols());
		for(int i=0; i<n; ++i) {
			_sol(i) = rhs[i];
		}
		delete[] rhs;
	}

	void MumpsSparse::setParameters(const ParameterList &list) {
		std::string name;
		for (ParameterList::iterator it = list.begin(); it != list.end(); ++it) {
			name = list.getName(it);
			if(name == "sym") {
				_id.sym = list.getValue<int>(it);
			}
			else if(name == "comm") {
				_id.comm_fortran = list.getValue<int>(it);
			}
			else if(name == "par") {
				_id.par = list.getValue<int>(it);
			}
			else if(name == "icntl[0]") {
				_id.icntl[0] = list.getValue<int>(it);
			}
			else if(name == "icntl[1]") {
				_id.icntl[1] = list.getValue<int>(it);
			}
			else if(name == "icntl[2]") {
				_id.icntl[2] = list.getValue<int>(it);
			}
			else if(name == "icntl[3]") {
				_id.icntl[3] = list.getValue<int>(it);
			}
		}
	}

}

