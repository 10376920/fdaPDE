#ifndef LINEAR_SOLVERS_SPLINEAR_SOLVER_H
#define LINEAR_SOLVERS_SPLINEAR_SOLVER_H

#include <Eigen/Sparse>
#include <string>
#include <list>
#include <memory>
#include <iostream>

namespace LinearSolvers {

class AbstractParameter {
	private:
	std::string _name;
	public:
//	AbstractParameter() = default;
	AbstractParameter(const std::string &name): _name(name) {};
	virtual ~AbstractParameter() = default;
	std::string getName() {return _name;}
	void printInfo() {std::cout << "_name = " << _name;}
};

template <typename T>
class Parameter: public AbstractParameter {
	private:
	T _value;
	public:
//	Parameter() = default;
	Parameter(const std::string &name, const T &value):
		AbstractParameter(name),
		_value(value) {};
	void set(const T &value) {_value = value;}
	T getValue() const {return _value;}
	void printInfo() {
		AbstractParameter::printInfo();
		std::cout << "_value = " << _value;
	}
};

class ParameterList {
	private:
	std::list<std::unique_ptr<AbstractParameter>> _list;
	public:
	typedef std::list<std::unique_ptr<AbstractParameter>>::const_iterator iterator;
	ParameterList() = default;
	template <typename T>
	void push_back(const std::string &s, T value);
	iterator begin() const;
	iterator end() const;
	std::string getName(const iterator &) const;
	template <typename T>
	T getValue(const iterator &) const;
};

class SpLinearSolver {
	protected:
	Eigen::MatrixXd _sol;
	public:
	virtual ~SpLinearSolver() = default;
	virtual void factorize(const Eigen::SparseMatrix<double> &)=0;
	virtual void solve(const Eigen::MatrixXd &)=0;
	virtual const Eigen::MatrixXd &getSolution() {return _sol;}
	virtual void setParameters(const ParameterList &list) {};
};

template <typename T>
void ParameterList::push_back(const std::string &s, T value) {
	_list.push_back(std::unique_ptr<Parameter<T>>(new Parameter<T>(s, value)));
	_list.front()->printInfo();
}

template <typename T>
T ParameterList::getValue(const ParameterList::iterator & it) const {
	dynamic_cast<Parameter<T>*>(it->get())->getValue();
}

}

#endif
