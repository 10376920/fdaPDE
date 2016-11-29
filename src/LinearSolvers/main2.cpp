#include <iostream>
#include <mpi.h>
#include <Eigen/Sparse>
#include "EigenSparseLU.h"
#include "MumpsSparse.h"
#include "MarketIO.h"
#include "timing.h"
#include "SpLinearSolver.h"


typedef Eigen::Triplet<double> coeff;
int main(int argc, char **argv) {
{
	MPI_Init(&argc,&argv);
	int me;
	int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &me);
	using namespace LinearSolvers;
	int n = 10;
	Eigen::SparseMatrix<double> A(n,n);
	Eigen::VectorXd b(Eigen::VectorXd::Zero(n));
//	loadMarket(A, "");
	b(0)=1;
	if (me == 0) {
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
		A.setFromTriplets(tripletAll.begin(), tripletAll.end());
		A.makeCompressed();

		EigenSparseLU solver1;
		solver1.factorize(A);
		solver1.solve(b);
		auto x1 = solver1.getSolution();
		std::cout << "x1 = " << x1 << std::endl;
	}
	
	
	ParameterList list;
	char one[7] = "1";
	list.set("sym", 1);
	if (!strcmp(argv[1], one)) {
		std::cout << "Parallel" << std::endl;
		std::cout << "argv[1]: "  << argv[1]<< std::endl;
		list.set("nproc", 2);
	}
	else {
		std::cout << "Serial" << std::endl;
		std::cout << "argv[1]: " << argv[1]<< std::endl;
		list.set("nproc", 1);
	}
	MumpsSparse solver2(list);
	timer clock;
	clock.start();
	solver2.factorize(A);
	clock.stop();
	solver2.solve(b);
	if (me == 0) {
		auto x2 = solver2.getSolution();
		std::cout << "x2 = " << x2 << std::endl;
	}
	list.set("nproc", 7);
	list.set("par", 0);
	list.set("sym", 1);
	solver2.setParameters(list);
	clock.start();
	solver2.factorize(A);
	clock.stop();
	solver2.solve(b);
	if (me == 0) {
		auto x2 = solver2.getSolution();
		std::cout << "x2 = " << x2 << std::endl;
	}
	}
	
	MPI_Finalize();

	return 0;
}
