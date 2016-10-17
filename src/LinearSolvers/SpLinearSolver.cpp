#include "SpLinearSolver.h"

namespace LinearSolvers {

ParameterList::iterator ParameterList::begin() const{
	return _list.begin();
}

ParameterList::iterator ParameterList::end() const{
	return _list.end();
}

std::string ParameterList::getName(const ParameterList::iterator &it) const{
	return (*it)->getName();
}

}
