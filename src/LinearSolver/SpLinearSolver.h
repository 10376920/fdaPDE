#ifndef LINEAR_SOLVERS_SPLINEAR_SOLVER_H
#define LINEAR_SOLVERS_SPLINEAR_SOLVER_H

#include <Eigen/Sparse>

namespace LinearSolvers {

class SpLinearSolver {
	protected:
	Eigen::MatrixXd _sol;
	public:
	virtual ~SpLinearSolver() = default;
	virtual void factorize(const Eigen::SparseMatrix<double> &)=0;
	virtual void solve(const Eigen::MatrixXd &)=0;
	virtual const Eigen::MatrixXd &getSolution() {return _sol;}
};

}

#endif
