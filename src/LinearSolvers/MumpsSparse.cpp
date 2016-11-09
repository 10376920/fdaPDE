#include "MumpsSparse.h"

namespace LinearSolvers {

	MumpsSparse::MumpsSparse() {
		int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &_myid);
		// TO DO : dare all'utente la possibilità di decidere quanti processori spannare
		//eventualmente fare una versione seriale n=0
		MPI_Comm_spawn("/home/pacs_student/Desktop/fdaPDE/src/LinearSolvers/child", NULL, 1, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &children, err);
		MPI_Comm intracomm;
 		MPI_Intercomm_merge(children, 0, &intracomm);
		_id.comm_fortran = (MUMPS_INT) MPI_Comm_c2f(intracomm);
		_id.sym = 0;
		_id.par = 1;
		_id.job = _job_init;
		dmumps_c(&_id);
		// _id.icntl[0]=6;
		// _id.icntl[1]=6; //
		// _id.icntl[2]=6;
		// _id.icntl[3]=2;
		_id.icntl[13]=200;
	}

	MumpsSparse::~MumpsSparse(){
		//andare out of scope perchè funzioni
		// ricordarsi di allineare i broadcast con i broadcast di child
		_id.job = _job_end;
		MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, children);
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
		//fare broadcast per coere_nz tra i processi
		//MPI_Bcast(_id.icntl, 40, MPI_INT, MPI_ROOT, children);
		MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, children);
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
		//fare broadcast per coere_nz tra i processi
		//MPI_Bcast(_id.icntl, 40, MPI_INT, MPI_ROOT, children);
		MPI_Bcast(&_id.job, 1, MPI_INT, MPI_ROOT, children);
		dmumps_c(&_id);

		_sol.resize(b.rows(), b.cols());
		for(int i=0; i<n; ++i) {
			_sol(i) = rhs[i];
		}
		delete[] rhs;
	}

	void MumpsSparse::setParameters(const ParameterList &list) {
		//occhio che incntl devono esssere sfasati tra umano e pc
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
			else if(name == "icntl[1]") {
				_id.icntl[0] = list.getValue<int>(it);
			}
			else if(name == "icntl[2]") {
				_id.icntl[1] = list.getValue<int>(it);
			}
			else if(name == "icntl[3]") {
				_id.icntl[2] = list.getValue<int>(it);
			}
			else if(name == "icntl[4]") {
				_id.icntl[3] = list.getValue<int>(it);
			}
		}
	}

}

