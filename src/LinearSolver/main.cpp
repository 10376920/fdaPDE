#include <iostream>
//#include "SpLinearSolver.h"
#include "MumpsSparse.h"
//#include "EigenSparseLU"
#include <Eigen/Sparse>
typedef Eigen::Triplet<double> coeff;
int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int n = 10;
    std::vector<Eigen::Triplet<double>> tripletAll;
    for (int i=1; i<n-1; ++i) {
        tripletAll.push_back(coeff(i, i, 2.0));
        tripletAll.push_back(coeff(i, i+1, -1.0));
        tripletAll.push_back(coeff(i, i-1, -1.0));
    }
    tripletAll.push_back(coeff(0,0,2.0));
    tripletAll.push_back(coeff(0,1,-1.0));
    tripletAll.push_back(coeff(n-1,n-1,2.0));
    tripletAll.push_back(coeff(n-1,n-2,-1.0));
    Eigen::SparseMatrix<double> A(n,n);
    A.setFromTriplets(tripletAll.begin(), tripletAll.end());
    A.makeCompressed();
	Eigen::VectorXd b(VectorXd::Zero(n));
	b(0)=1;
//	EigenSparseLU solver1;
	MumpsSparse solver2;
//	solver1.factorize(A);
//	auto x1 = solver1.solve(b);
//	std::cout << "x1 = " << x1;
	int myid;
	int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	solver2.factorize(A);
	solver2.solve(b);
	if (myid == 0) {
		auto x2 = solver2.getsolution();
		std::cout << "x2 = " << x2;
	}
	MPI_Finalize();
	return 0;
}
