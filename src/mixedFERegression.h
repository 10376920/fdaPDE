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
//#include "LinearSolvers/SpLinearSolver.h"
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

		
		SpMat psi_;
		SpMat AMat_;
		SpMat MMat_;
		SpMat A_;		// system matrix with psi^T*psi in north-west block
		MatrixXr U_;	// psi^T*W padded with zeros
		
		
		std::unique_ptr<LinearSolvers::SpLinearSolver> Adec_;
		Eigen::PartialPivLU<MatrixXr> Gdec_;
		
		VectorXr _b;                     //!A Eigen::VectorXr: Stores the system right hand side.
		std::vector<VectorXr> _solution; //!A Eigen::VectorXr : Stores the system solution
		std::vector<Real> _dof;
		std::vector<Real> _var;
		std::string _finalRNGstate;

		Eigen::PartialPivLU<MatrixXr> WTWinv_;
		bool isWTWfactorized_;

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
			LSProxy<LinearSolvers::EigenSparseLU> dummy1("EigenSparseLU");
			LSProxy<LinearSolvers::MumpsSparse> dummy2("MumpsSparse");
			LSFactory & LSfactory=LSFactory::Instance();
			Adec_ = LSfactory.create("EigenSparseLU"); //TODO get name from input
		};
		
		//!A Destructor
		//~Model(){};
		
		void smoothLaplace();
		void smoothEllipticPDE();
		void smoothEllipticPDESpaceVarying();

		 //! A template member for the system resolution.
         /*! the template P is the solutor, possible choices are: SpLu, SpQR, SpCholesky,SpConjGrad.
          *  the solution is stored in _solution		
		*/

		template<typename P>
		void solve(UInt output_index);
		void system_factorize();
		template<typename Derived>
		MatrixXr system_solve(const Eigen::MatrixBase<Derived>&);
		//! A inline member that returns a VectorXr, returns the whole _solution. 
		inline std::vector<VectorXr> const & getSolution() const{return _solution;};
		inline std::vector<Real> const & getDOF() const{return _dof;};
		inline std::vector<Real> const & getVar() const{return _var;};
		inline std::string const & getFinalRNGstate() const{return _finalRNGstate;}

		//Real eval_sol(MeshTria const & mesh,VectorXr Point p);
		//! A member for printing the solution.
		//void printSolution(std::ostream & out) {out<<_solution; out<<"dim"<<_solution.rows()<<"x"<<_solution.cols();};
		//apply dirichlet boundary condition (first attempt)
		//void applyDir(int* bcindex, double* bcvalues, int n_b, UInt order);
		 //! A normal member taking two arguments: Dirichlet Boundary condition
		 /*!
		  * This member applies Dirichlet boundary conditions on the linear system with the penalization method.
		\param bcindex is a const reference to vector<int> : the global indexes of the nodes to which the boundary condition has to be applied.
		\param bcvalues is a const reference to vector<double> : the values of the boundary conditions relative to bcindex.
		\the method modifies _coeffmatrix and _b
		*/		
		void addDirichletBC();
		void getRightHandData(VectorXr& rightHandData);
		void computeDegreesOfFreedom(UInt output_index, Real lambda);
		void computeDegreesOfFreedomExact(UInt output_index, Real lambda);
		void computeDegreesOfFreedomStochastic(UInt output_index, Real lambda);
};

#include "mixedFERegression_imp_woodbury_whole.h"

#endif
