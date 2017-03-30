#ifndef __MIXEDFEREGRESSION_HPP__
#define __MIXEDFEREGRESSION_HPP__

#include "fdaPDE.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "mesh.h"
#include "param_functors.h"
#include "regressionData.h"
#include "solver.h"
#include <memory>
#include "EigenSparseLU.h"
#include "MumpsSparse.h"
#include "Proxy.hpp"
#include "Factory.hpp"

//! A LinearSystem class: A class for the linear system construction and resolution.

template<typename InputHandler, typename Integrator, UInt ORDER>
class MixedFERegression{
	private:
		using LSFactory=GenericFactory::Factory<LinearSolvers::SpLinearSolver,std::string>;
		template<class C>
		using LSProxy=GenericFactory::Proxy<LSFactory,C>;
		const MeshHandler<ORDER> &mesh_;
		const InputHandler& regressionData_;
		std::vector<coeff> tripletsData_;

		SpMat A_;		// System matrix with psi^T*psi in north-west block
		SpMat AMat_;	// North-east block of system matrix A_
		SpMat MMat_;	// South-east block of system matrix A_
		SpMat psi_;
		MatrixXr U_;	// psi^T*W padded with zeros
		
		
		std::unique_ptr<LinearSolvers::SpLinearSolver> Adec_; // Stores the factorization of A_
		Eigen::PartialPivLU<MatrixXr> Gdec_;	// Stores factorization of G =  C + [V * A^-1 * U]
		Eigen::PartialPivLU<MatrixXr> WTWinv_;	// Stores the factorization of W^T * W
		bool isWTWfactorized_;

		VectorXr _b;                     //!A Eigen::VectorXr: Stores the system right hand side.
		std::vector<VectorXr> _solution; //!A Eigen::VectorXr : Stores the system solution
		std::vector<Real> _dof;
		std::vector<Real> _var;
		std::string _finalRNGstate;

		void setPsi();
		void buildA(const SpMat& Psi,  const SpMat& AMat,  const SpMat& MMat);
		MatrixXr LeftMultiplybyQ(const MatrixXr& u);

	public:
		//!A Constructor.
		MixedFERegression(const MeshHandler<ORDER>& mesh, const InputHandler& regressionData):
			mesh_(mesh),
			regressionData_(regressionData),
			isWTWfactorized_(false)
		{
			std::string solver_name = regressionData_.getSolver();
			LSProxy<LinearSolvers::EigenSparseLU> dummy1("EigenSparseLU");
			LSProxy<LinearSolvers::MumpsSparse> dummy2("MumpsSparse");
			LSFactory & LSfactory=LSFactory::Instance();
			Adec_ = LSfactory.create(solver_name); 
			// Definition of a list of parameters for the solver
			LinearSolvers::ParameterList list;
			if(solver_name == "MumpsSparse") {
				list.set("icntl[1]", -1);
				list.set("icntl[2]", -1);
				list.set("icntl[3]", -1);
				list.set("icntl[4]", 0);
				list.set("icntl[14]", 200);
				list.set("sym", 2);
				list.set("nproc", regressionData_.getnprocessors());
				list.set("hosts", regressionData_.getHosts());
			}
			Adec_->setParameters(list);
		};
		
		void smoothLaplace();
		void smoothEllipticPDE();
		void smoothEllipticPDESpaceVarying();

		//  //! A template member for the system resolution.
  //        ! the template P is the solutor, possible choices are: SpLu, SpQR, SpCholesky,SpConjGrad.
  //         *  the solution is stored in _solution		
		// template<typename P>
		// void solve(UInt output_index);

		inline std::vector<VectorXr> const & getSolution() const{return _solution;};
		inline std::vector<Real> const & getDOF() const{return _dof;};
		inline std::vector<Real> const & getVar() const{return _var;};
		inline std::string const & getFinalRNGstate() const{return _finalRNGstate;}

		void addDirichletBC();
		void getRightHandData(VectorXr& rightHandData);
		void computeDegreesOfFreedom(UInt output_index, Real lambda);
		void computeDegreesOfFreedomExact(UInt output_index, Real lambda);
		void computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);

		void system_factorize();
		template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);
};

#include "mixedFERegression_imp.h"

#endif
