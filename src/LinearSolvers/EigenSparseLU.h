#ifndef LINEAR_SOLVERS_EIGEN_SPARSE_LU_H
#define LINEAR_SOLVERS_EIGEN_SPARSE_LU_H

#include "SpLinearSolver.h"

namespace LinearSolvers {

class EigenSparseLU: public SpLinearSolver {
	private:
	Eigen::SparseLU<Eigen::SparseMatrix<double>> _solver;
	public:
	virtual void factorize(const Eigen::SparseMatrix<double> &);
	virtual void solve(const Eigen::MatrixXd &);
};

}

#endif
