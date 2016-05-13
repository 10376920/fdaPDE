#ifndef __MIXEDFEREGRESSION_IMP_HPP__
#define __MIXEDFEREGRESSION_IMP_HPP__

#include <iostream>
#include<random>
#include "timing.h"
#include "woodbury.hpp"

////build system matrix in sparse format SWest non serve??
//void MixedFE::build(SpMat & L, SpMat& opMat, SpMat& opMat2, SpMat& mass, const VectorXr& righthand, const VectorXr& forcing_term )
//{
//	if(opMat2.rows()!=mass.rows() || L.rows()!=opMat.rows() )
//	 std::cerr<<"incompatible blocks, dimension mismatch"<<std::endl;
//
//	if(righthand.size()+forcing_term.size()!=2*opMat.rows() )
//	 std::cerr<<"incompatible right hand side forcing term, dimension mismatch"<<std::endl;
//
//
//	UInt nnodes=opMat.rows();
//
//
//	/*SpMat tempStiff(_femhandl.getStiff()),
//			tempMass(_femhandl.getMass());*/
//
//	//I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
//	_coeffmatrix.resize(2*nnodes,2*nnodes);
//	std::vector<int> entries(2*nnodes);
//	UInt number=0,number2=0;
//
//	for(auto k=0; k<nnodes; k++)
//	{
//
//	 number=L.col(k).nonZeros()+opMat.col(k).nonZeros();
//
//	 entries[k]=number;
//
//	 number2=opMat2.col(k).nonZeros()+mass.col(k).nonZeros();
//
//	 entries[nnodes+k]=number2;
//	}
//
//	_coeffmatrix.reserve(entries);
//	////building system matrix
//	for(auto i=0; i<nnodes; i++)
//	{
//		// north-west block from matrix L, cycling over non-zero elements
//		for(SpMat::InnerIterator it(L,i); it; ++it)
//		{
//			_coeffmatrix.insert(it.index(),i)=it.value();
//		}
//		// north-east cycling over non-zero elements
//		for(SpMat::InnerIterator it(opMat,i); it; ++it)
//		{
//			_coeffmatrix.insert(it.index(),i+nnodes)=it.value(); //north-east block
//			//_coeffmatrix.insert(it.index()+nnodes,i)=it.value(); //south-west block
//		}
//		// south-west block cycling over non-zero elements
//		for(SpMat::InnerIterator it(opMat2,i); it; ++it)
//		{
//		    _coeffmatrix.insert(it.index()+nnodes,i)=it.value(); //south-west block
//	    }
//
//		//south-east block from Mass matrix, cycling over non-zero elements
//		for(SpMat::InnerIterator it(mass,i); it; ++it)
//		{
//			_coeffmatrix.insert(it.index()+nnodes,i+nnodes)=it.value();
//		}
//	}
//	_coeffmatrix.makeCompressed();
//
//	_b.resize(nnodes*2);
//
//	_b.topRows(righthand.rows())=righthand;
//
//	_b.bottomRows(forcing_term.rows())=forcing_term;
//
//}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::buildCoeffMatrix(const SpMat& DMat,  const SpMat& AMat,  const SpMat& MMat)
{
	//I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
	//_coeffmatrix.setFromTriplets(tripletA.begin(),tripletA.end());

	UInt nnodes = mesh_.num_nodes();

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(DMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
	  }
	for (int k=0; k<MMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(MMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	_coeffmatrix.setZero();
	_coeffmatrix.resize(2*nnodes,2*nnodes);
	_coeffmatrix.setFromTriplets(tripletAll.begin(),tripletAll.end());
	_coeffmatrix.makeCompressed();
	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::addDirichletBC(const vector<int>& bcindex, const vector<Real>& bcvalues)
{


	UInt id1,id3;

	UInt nnodes = mesh_.num_nodes();

	const std::vector<UInt>& bc_indices = regressionData_.getDirichletIndices();
	const std::vector<Real>& bc_values = regressionData_.getDirichletValues();
	UInt nbc_indices = bc_indices.size();

	Real pen=10e20;

	for( auto i=0; i<nbc_indices; i++)
	 {
			id1=bcindex[i];
			id3=id1+nnodes;

			//_coeffmatrix.prune([id1,id3](int i, int j, Real) { return (i!=id1 && i!=id3); });

			_coeffmatrix.coeffRef(id1,id1)=pen;
			_coeffmatrix.coeffRef(id3,id3)=pen;


			_b(id1)+=bc_values[i]*pen;
			_b(id3)=0;
	 }

	_coeffmatrix.makeCompressed();
}

//construct NW block of the system matrix in Ramsay with covariates format
//!! Depends on setPsi and setQ
template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::getDataMatrix(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		//UInt nlocations = regressionData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		if (regressionData_.getCovariates().rows() == 0)
			DMat = psi_.transpose()*psi_;
		else
		{
			DMat = (SpMat(psi_.transpose())*Q_*psi_).sparseView();
		}
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::getDataMatrixByIndices(SpMat& DMat)
{
		UInt nnodes = mesh_.num_nodes();
		UInt nlocations = regressionData_.getNumberofObservations();

		DMat.resize(nnodes,nnodes);

		if (regressionData_.getCovariates().rows() == 0)
		{
			DMat.reserve(1);
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index = regressionData_.getObservationsIndices()[i];
				DMat.insert(index,index) = 1;
			}
		}
		else
		{
			//May be inefficient
			for (auto i = 0; i<nlocations; ++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				for (auto j = 0; j<nlocations; ++j)
				{
					auto index_j = regressionData_.getObservationsIndices()[j];
					DMat.insert(index_i,index_j) = Q_(i,j);
				}
			}
		}
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::setPsi(){

	//std::cout<<"Data Matrix Computation by Basis Evaluation.."<<std::endl;
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	//cout<<"Nodes number "<<nnodes<<"Locations number "<<nlocations<<endl;

	//std::vector<coeff> entries;
	//entries.resize((ORDER * 3)*nlocations);


	psi_.resize(nlocations, nnodes);
	//psi_.reserve(Eigen::VectorXi::Constant(nlocations,ORDER*3));
	if (regressionData_.isLocationsByNodes()){
		std::vector<coeff> tripletAll;
		auto k = regressionData_.getObservationsIndices();
		tripletAll.reserve(k.size());
		for (int i = 0; i< k.size(); ++i){
			tripletAll.push_back(coeff(i,k[i],1.0));
		}
    	psi_.setFromTriplets(tripletAll.begin(),tripletAll.end());
    	psi_.makeCompressed();
    }
    else {
		Triangle<ORDER*3> tri_activated;
		Eigen::Matrix<Real,ORDER * 3,1> coefficients;

		Real evaluator;

		for(UInt i=0; i<nlocations;i++)
		{
			tri_activated = mesh_.findLocationNaive(regressionData_.getLocations()[i]);
			if(tri_activated.getId() == Identifier::NVAL)
			{
				#ifdef R_VERSION_
				Rprintf("ERROR: Point %d is not in the domain, remove point and re-perform smoothing\n", i+1);
				#else
				std::cout << "ERROR: Point " << i+1 <<" is not in the domain\n";
				#endif
			}else
			{
				for(UInt node = 0; node < ORDER*3 ; ++node)
				{
					coefficients = Eigen::Matrix<Real,ORDER * 3,1>::Zero();
					coefficients(node) = 1; //Activates only current base
					evaluator = evaluate_point<ORDER>(tri_activated, regressionData_.getLocations()[i], coefficients);
					psi_.insert(i, tri_activated[node].getId()) = evaluator;
				}
			}
		}

		psi_.makeCompressed();
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::setQ()
{
	//std::cout<<"Computing Orthogonal Space Projection Matrix"<<std::endl;
	Q_.resize(H_.rows(),H_.cols());
	Q_ = -H_;
	for (int i=0; i<H_.rows();++i)
	{
		Q_(i,i) += 1;
	}
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::setH()
{
	//std::cout<<"Computing Projection Matrix"<<std::endl;
	//UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	//regressionData_.printCovariates(std::cout);
	MatrixXr W(this->regressionData_.getCovariates());
	//std::cout<<"W "<< W <<std::endl;
	//total number of mesh nodes
	//UInt nnodes = mesh_.num_nodes();
	if(regressionData_.isLocationsByNodes())
	{
		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
		for (auto i=0; i<nlocations;++i)
		{
			auto index_i = regressionData_.getObservationsIndices()[i];
			for (auto j=0; j<W.cols();++j)
			{
				W_reduced(i,j) = W(index_i,j);
			}
		}
		W = W_reduced;
	}


	MatrixXr WTW(W.transpose()*W);

	H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::getRightHandData(VectorXr& rightHandData)
{
	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();
	//rightHandData.resize(nnodes);
	rightHandData = VectorXr::Zero(nnodes);

	if(regressionData_.getCovariates().rows() == 0)
	{
		if(regressionData_.isLocationsByNodes())
		{

			for (auto i=0; i<nlocations;++i)
			{
				auto index_i = regressionData_.getObservationsIndices()[i];
				rightHandData(index_i) = regressionData_.getObservations()[i];
			}
		}
		else
		{
			rightHandData=psi_.transpose()*regressionData_.getObservations();
		}
	}
	else
	{
		rightHandData=psi_.transpose()*LeftMultiplybyQ(regressionData_.getObservations());
	}
}


template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::computeDegreesOfFreedom(UInt output_index, Real lambda)
{
	std::cout << "Starting GCV computation" << std::endl;
	timer clock1;
	clock1.start();

	UInt nnodes = mesh_.num_nodes();
	UInt nlocations = regressionData_.getNumberofObservations();

	// genero matrice aleatoria
	std::default_random_engine generator;
	std::bernoulli_distribution distribution(0.5);
	int nrealizations=100;
	MatrixXr u(nlocations, nrealizations);
	for (int j=0; j<nrealizations; ++j) {
		for (int i=0; i<nlocations; ++i) {
			if (distribution(generator)) {
				u(i,j) = 1.0;
			}
			else {
				u(i,j) = -1.0;
			}
		}
	}

	MatrixXr b = MatrixXr::Zero(2*nnodes,u.cols());
	b.topRows(nnodes) = psi_.transpose()* LeftMultiplybyQ(u);

	MatrixXr x = system_solve(b);
	MatrixXr uTpsi = u.transpose()*psi_;
	Real edf = 0;
	for (int i=0; i<nrealizations; ++i) {
		edf += uTpsi.row(i).dot(x.col(i).head(nnodes));
	}
	edf /= nrealizations;
	if (regressionData_.getCovariates().rows() != 0) {
		edf += regressionData_.getCovariates().cols();
	}
	_dof[output_index] = edf;
	clock1.stop();
}



//Implementation kept from Sangalli et al
template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::smoothLaplace()
{
	//std::cout<<"Laplace Penalization - Order: "<<ORDER<<std::endl;

	//UInt ndata=regressionData_.getObservations().size();
	UInt nnodes=mesh_.num_nodes();

	FiniteElement<Integrator, ORDER> fe;

	typedef EOExpr<Mass> ETMass;
	typedef EOExpr<Stiff> ETStiff;

	Mass EMass;
	Stiff EStiff;

	ETMass mass(EMass);
	ETStiff stiff(EStiff);

	setPsi();


//	if(!regressionData_.getCovariates().rows() == 0)
//	{
//		setH();
//		setQ();
//	}

//    if(!regressionData_.isLocationsByNodes())
//    {
//    	getDataMatrix(DMat_);
//    }
//    else
//    {
//    	getDataMatrixByIndices(DMat_);
//    }
    //std::cout<<"Block Data"<<DMat_<<std::endl;


    Assembler::operKernel(stiff, mesh_, fe, AMat_);
    Assembler::operKernel(mass, mesh_, fe, MMat_);

    VectorXr rightHandData;
    getRightHandData(rightHandData);
    _b = VectorXr::Zero(2*nnodes);
    _b.topRows(nnodes)=rightHandData;
    //std::cout<<"b vector"<<_b;

    _solution.resize(regressionData_.getLambda().size());
    _dof.resize(regressionData_.getLambda().size());

    for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
    	//build(tripletsData_,(-regressionData_.getLambda())*stiff, (-regressionData_.getLambda())*mass, righthand, forcing);

    	Real lambda = regressionData_.getLambda()[i];
    	SpMat AMat_lambda = (-lambda)*AMat_;
    	SpMat MMat_lambda = (-lambda)*MMat_;
    	this->buildA(psi_, AMat_lambda, MMat_lambda);

    	//Appling border conditions if necessary
    	if(regressionData_.getDirichletIndices().size() != 0)
    		addDirichletBC(regressionData_.getDirichletIndices(), regressionData_.getDirichletValues());

    	//prova.solveSystem<SpConjGrad>();
    	system_factorize();
    	_solution[i] = this->template system_solve(this->_b);
    	if(regressionData_.computeDOF())
    		computeDegreesOfFreedom(i, lambda);
    	else
    		_dof[i] = -1;

	}

}

//Implementation kept from Sangalli et al
template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::smoothEllipticPDE()
{
	//std::cout<<"Elliptic PDE Penalization - Order: "<<ORDER<<std::endl;

	//UInt ndata=regressionData_.getObservations().size();
	UInt nnodes=mesh_.num_nodes();

	FiniteElement<Integrator, ORDER> fe;

	typedef EOExpr<Mass> ETMass;
	typedef EOExpr<Stiff> ETStiff;
	typedef EOExpr<Grad> ETGrad;

	Mass EMass;
	Stiff EStiff;
	Grad EGrad;

	ETMass mass(EMass);
	ETStiff stiff(EStiff);
	ETGrad grad(EGrad);

	setPsi();

	if(!regressionData_.getCovariates().rows() == 0)
	{
		setH();
		setQ();
	}

    if(!regressionData_.isLocationsByNodes())
    {
    	getDataMatrix(DMat_);
    }
    else
    {
    	getDataMatrixByIndices(DMat_);
    }
    //std::cout<<"Block Data"<<DMat<<std::endl;



    const Real& c = regressionData_.getC();
    const Eigen::Matrix<Real,2,2>& K = regressionData_.getK();
    const Eigen::Matrix<Real,2,1>& beta = regressionData_.getBeta();
    Assembler::operKernel(c*mass+stiff[K]+dot(beta,grad), mesh_, fe, AMat_);
    Assembler::operKernel(mass, mesh_, fe, MMat_);

    VectorXr rightHandData;
    getRightHandData(rightHandData);
    _b = VectorXr::Zero(2*nnodes);
    _b.topRows(nnodes)=rightHandData;
    //std::cout<<"b vector"<<_b;

    _solution.resize(regressionData_.getLambda().size());
    _dof.resize(regressionData_.getLambda().size());

    for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
    	//build(tripletsData_,(-regressionData_.getLambda())*stiff, (-regressionData_.getLambda())*mass, righthand, forcing);

    	Real lambda = regressionData_.getLambda()[i];
    	SpMat AMat_lambda = (-lambda)*AMat_;
    	SpMat MMat_lambda = (-lambda)*MMat_;
    	this->buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);
    	this->buildA(psi_, AMat_lambda, MMat_lambda);

    	//std::cout<<"AMat"<<std::endl<<_coeffmatrix;


    	//Appling border conditions if necessary
    	if(regressionData_.getDirichletIndices().size() != 0)
    		addDirichletBC(regressionData_.getDirichletIndices(), regressionData_.getDirichletValues());

    	//prova.solveSystem<SpConjGrad>();
    	this-> template solve<SpLU>(i);
    	if(regressionData_.computeDOF())
    		computeDegreesOfFreedom(i, lambda);
    	else
    		_dof[i] = -1;

	}

}


//Implementation kept from Sangalli et al
template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::smoothEllipticPDESpaceVarying()
{
	//std::cout<<"Space-varying Coefficient - Elliptic PDE Penalization - Order: "<<ORDER<<std::endl;
		//UInt ndata=regressionData_.getObservations().size();
		UInt nnodes=mesh_.num_nodes();

		FiniteElement<Integrator, ORDER> fe;

		typedef EOExpr<Mass> ETMass;
		typedef EOExpr<Stiff> ETStiff;
		typedef EOExpr<Grad> ETGrad;

		Mass EMass;
		Stiff EStiff;
		Grad EGrad;

		ETMass mass(EMass);
		ETStiff stiff(EStiff);
		ETGrad grad(EGrad);

	setPsi();

	if(!regressionData_.getCovariates().rows() == 0)
	{
		setH();
		setQ();
	}

	if(!regressionData_.isLocationsByNodes())
	{
		getDataMatrix(DMat_);
	}
	else
	{
		getDataMatrixByIndices(DMat_);
	}
	//std::cout<<"Block Data"<<DMat_<<std::endl;

	const Reaction& c = regressionData_.getC();
	const Diffusivity& K = regressionData_.getK();
	const Advection& beta = regressionData_.getBeta();
	Assembler::operKernel(c*mass+stiff[K]+dot(beta,grad), mesh_, fe, AMat_);
	Assembler::operKernel(mass, mesh_, fe, MMat_);

	const ForcingTerm& u = regressionData_.getU();
	//for(auto i=0;i<18;i++) std::cout<<u(i)<<std::endl;
	VectorXr forcingTerm;
	Assembler::forcingTerm(mesh_,fe, u, forcingTerm);

	VectorXr rightHandData;
	getRightHandData(rightHandData);
	//_b.resize(2*nnodes);
	_b = VectorXr::Zero(2*nnodes);
	_b.topRows(nnodes)=rightHandData;
	//std::cout<<"Forcing Term "<<std::cout<<forcingTerm<<"END";
	_b.bottomRows(nnodes)=forcingTerm;
	//std::cout<<"b vector"<<_b;

	_solution.resize(regressionData_.getLambda().size());
	_dof.resize(regressionData_.getLambda().size());

	for(UInt i = 0; i<regressionData_.getLambda().size(); ++i)
	{
		//build(tripletsData_,(-regressionData_.getLambda())*stiff, (-regressionData_.getLambda())*mass, righthand, forcing);

		Real lambda = regressionData_.getLambda()[i];
		SpMat AMat_lambda = (-lambda)*AMat_;
		SpMat MMat_lambda = (-lambda)*MMat_;
		this->buildCoeffMatrix(DMat_, AMat_lambda, MMat_lambda);
		this->buildA(psi_, AMat_lambda, MMat_lambda);

		//std::cout<<"AMat"<<std::endl<<_coeffmatrix;


		//Appling border conditions if necessary
		if(regressionData_.getDirichletIndices().size() != 0)
			addDirichletBC(regressionData_.getDirichletIndices(), regressionData_.getDirichletValues());

		//prova.solveSystem<SpConjGrad>();

		//std::cout<< _coeffmatrix;
		this-> template solve<SpLU>(i);
		if(regressionData_.computeDOF())
			computeDegreesOfFreedom(i, lambda);
		else
			_dof[i] = -1;

	}

}


//solve sparse system with P method

template<typename InputHandler, typename Integrator, UInt ORDER>
template <typename P>
void MixedFERegression<InputHandler,Integrator,ORDER>::solve(UInt output_index)
{
	//std::cout<<this->_coeffmatrix;
	this->_solution[output_index].resize(this->_coeffmatrix.rows());
	P::solve(this->_coeffmatrix,this->_b,this->_solution[output_index]);
}

template<typename InputHandler, typename Integrator, UInt ORDER>
MatrixXr MixedFERegression<InputHandler,Integrator,ORDER>::LeftMultiplybyQ(const MatrixXr& u)
{	
	if (regressionData_.getCovariates().rows() == 0){
		return u;
	}
	else{
		MatrixXr W(this->regressionData_.getCovariates());
		if (isWTWfactorized_ == false ){
			WTWinv_.compute(W.transpose()*W);
			isWTWfactorized_=true;
		}
		MatrixXr Pu= W*WTWinv_.solve(W.transpose()*u);
		return u-Pu;
	}

}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::buildA(const SpMat& Psi,  const SpMat& AMat,  const SpMat& MMat) {
		//I reserve the exact memory for the nonzero entries of each row of the coeffmatrix for boosting performance
	//_coeffmatrix.setFromTriplets(tripletA.begin(),tripletA.end());

	UInt nnodes = mesh_.num_nodes();
	
	SpMat DMat = Psi.transpose()*Psi;

	std::vector<coeff> tripletAll;
	tripletAll.reserve(DMat.nonZeros() + 2*AMat.nonZeros() + MMat.nonZeros());

	for (int k=0; k<DMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(DMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row(), it.col(),it.value()));
	  }
	for (int k=0; k<MMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(MMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.col(), it.row()+nnodes,it.value()));
	  }
	for (int k=0; k<AMat.outerSize(); ++k)
	  for (SpMat::InnerIterator it(AMat,k); it; ++it)
	  {
		  tripletAll.push_back(coeff(it.row()+nnodes, it.col(), it.value()));
	  }

	A_.setZero();
	A_.resize(2*nnodes,2*nnodes);
	A_.setFromTriplets(tripletAll.begin(),tripletAll.end());
	A_.makeCompressed();
	//std::cout<<"Coefficients' Matrix Set Correctly"<<std::endl;
}

template<typename InputHandler, typename Integrator, UInt ORDER>
void MixedFERegression<InputHandler,Integrator,ORDER>::system_factorize() {
	timer clock1, clock2;
	std::cout << "Decomposing A" << std::endl;
	clock1.start();
	Adec_.compute(A_);
	clock1.stop();
	std::cout << "Decomposing D" << std::endl;
	clock2.start();
	UInt nnodes = mesh_.num_nodes();
	if (regressionData_.getCovariates().rows() != 0) {
		MatrixXr W(this->regressionData_.getCovariates());
		U_ = MatrixXr::Zero(2*nnodes, W.cols());
		U_.topRows(nnodes) = psi_.transpose()*W;
		MatrixXr G = -W.transpose()*W + U_.transpose()*Adec_.solve(U_);
		Gdec_.compute(G);
	}
	clock2.stop();
}

template<typename InputHandler, typename Integrator, UInt ORDER>
template<typename Derived>
MatrixXr MixedFERegression<InputHandler,Integrator,ORDER>::system_solve(const Eigen::MatrixBase<Derived> &b) {
	timer clock1, clock2;
	std::cout << "Solving FEM: 1" << std::endl;
	clock1.start();
	MatrixXr x1 = Adec_.solve(b);
	clock1.stop();
	clock2.start();
	if (regressionData_.getCovariates().rows() != 0) {
		std::cout << "Solving FEM: 2" << std::endl;
		clock2.start();
		MatrixXr x2 = Gdec_.solve(U_.transpose()*x1);
		x1 -= Adec_.solve(U_*x2);
		clock2.stop();
	}
	return x1;
}

#endif
