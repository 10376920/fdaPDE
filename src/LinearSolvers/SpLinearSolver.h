#ifndef LINEAR_SOLVERS_SPLINEAR_SOLVER_H
#define LINEAR_SOLVERS_SPLINEAR_SOLVER_H

#include <Eigen/Sparse>
#include <string>
#include <unordered_map>
#include <utility>
#include <memory>
#include <iostream>

namespace LinearSolvers {

class AbstractParameter {
	public:
	virtual ~AbstractParameter() = default;
};

template <typename T>
class Parameter: public AbstractParameter {
	private:
	T _value;
	public:
	Parameter(const T &value): _value(value) {};
	void set(const T &value) {_value = value;}
	T getValue() const {return _value;}
	void printInfo() const{
		std::cout << "_value = " << _value << std::endl;;
	}
};

class ParameterList {
	private:
	std::unordered_map<std::string, std::unique_ptr<AbstractParameter>> _map;
	public:
	ParameterList() = default;
	template <typename T>
	void set(const std::string &s, T value);
	template <typename T>
	T getValue(const std::string & name, T defaultValue) const;
//	bool there_are(const std::string &s) const;
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
void ParameterList::set(const std::string &name, T value) {
	std::unique_ptr<Parameter<T>> value_ptr(new Parameter<T>(value));
	_map[name] = std::move(value_ptr);
}

template <typename T>
T ParameterList::getValue(const std::string & name, T defaultValue) const {
	auto iter = _map.find(name);
	if (iter != _map.end()) {
		auto abstract_ptr = iter->second.get();
		auto concrete_ptr = dynamic_cast<Parameter<T>*>(abstract_ptr);
		return concrete_ptr->getValue();
	}
	else {
		return defaultValue;
	}
}

//bool ParameterList::there_are(const std::string &s) const {
//	return _map.count(s);
//}

} // End namespace LinearSolvers

#endif
