#include "EigenSparseLU.h"

namespace LinearSolvers {

	void EigenSparseLU::factorize(const Eigen::SparseMatrix<double> &mat) {
		_solver.compute(mat);
	}

	void EigenSparseLU::solve(const Eigen::MatrixXd &b) {
		_sol = _solver.solve(b);
	}

}
