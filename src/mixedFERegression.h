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
#include "LinearSolvers/EigenSparseLU.h"
#include "LinearSolvers/MumpsSparse.h"
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

		//SpMat NWblock_;
		SpMat psi_;
		SpMat DMat_;
		SpMat AMat_;
		SpMat MMat_;
		SpMat A_;		// _coeffmatrix with psi^T*psi in north-west block
		MatrixXr U_;	// psi^T*W padded with zeros
		
		
		std::unique_ptr<LinearSolvers::SpLinearSolver> Adec_;
		Eigen::PartialPivLU<MatrixXr> Gdec_;
		

		MatrixXr Q_;
		MatrixXr H_;

		SpMat _coeffmatrix;
		SpMat _coeffmatrix2;
		VectorXr _b;                     //!A Eigen::VectorXr: Stores the system right hand side.
		std::vector<VectorXr> _solution; //!A Eigen::VectorXr : Stores the system solution
		std::vector<Real> _dof;
		std::vector<Real> _var;

		Eigen::PartialPivLU<MatrixXr> WTWinv_;
		bool isWTWfactorized_;

		void setPsi();
		void setQ();
		void setH();

		void buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat);
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
		
		 //! A normal member taking six arguments. Build the linear system.
		 /*!
		 * This member builds the system matrix and system right hand side, storing them in _coeffmatr and _b
		\param L is a const reference Eigen::SparseMatrix<> : NW block of the system matrix .
		\param opMat is a const reference Eigen::SparseMatrix<> : NE block of the system matrix .
		\param opMat2 is a const reference Eigen::SparseMatrix<> : SW block of the system matrix .
		\param mass is a const reference Eigen::SparseMatrix<> : SE block of the system matrix .
		\the method modifies _coeffmatrix and _b.
		*/
//		void build( SpMat& L,  SpMat& opMat,  SpMat& opMat2,  SpMat& mass,
//						  const VectorXr& rightside, const VectorXr& forcing_term );

		//
		//		|DMat | AMat^T  |
		//		|AMat | MMat	|



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
		void addDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues);
		void getDataMatrix(SpMat& DMat);
		void getDataMatrixByIndices(SpMat& DMat);
		void getRightHandData(VectorXr& rightHandData);
#ifdef ORIGINAL_VERSION
		void computeDegreesOfFreedom(UInt output_index);
#else
        void computeDegreesOfFreedom(UInt output_index, Real lambda);
#endif
};

#ifdef ORIGINAL_VERSION
#include "mixedFERegression_imp_original.h"
#endif

#ifdef SCHUR_VERSION
#include "mixedFERegression_imp_SCHUR.h"
#endif

//#ifdef MUMPS_VERSION
//#include "mixedFERegression_imp_MUMPS.h"
//#endif

#ifdef MUMPS_WHOLE_VERSION
#include "mixedFERegression_imp_mumps_whole.h"
#endif

#ifdef STOCHASTIC_VERSION
#include "mixedFERegression_imp_stochastic.h"
#endif

#ifdef TEMP_VERSION
#include "mixedFERegression_imp_temp.h"
#endif

#ifdef WOODBURY_WHOLE_VERSION
#include "mixedFERegression_imp_woodbury_whole.h"
#endif

#ifdef WOODBURY_DECOMPOSEQ_VERSION
#include "mixedFERegression_imp_woodbury_decomposeQ.h"
#endif

#endif
