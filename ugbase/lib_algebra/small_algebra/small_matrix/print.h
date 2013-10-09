/*
 * print.h
 *
 *  Created on: 18.09.2013
 *      Author: mrupp
 */

#ifndef SMALL_ALGEBRA_PRINT_H_
#define SMALL_ALGEBRA_PRINT_H_

#include "densematrix.h"
#include <string>
#include <sstream>
#include "common/log.h"

namespace ug{
template<typename TStorage>
void MaplePrint(const DenseMatrix<TStorage> &A, const char *name)
{
	UG_LOG(name << " = " << MapleString(A));
}

template<typename TStorage>
std::string NonzeroString(const DenseMatrix<TStorage> &A, const char *brackets, const char *seperator, int base, const char *name)
{
	std::stringstream ss;
	ss << std::setprecision(9);
	for(size_t r=0; r<A.num_rows(); ++r)
	{
		for(size_t c=0; c<A.num_cols(); ++c)
			if(A(r, c) != 0.0)
				ss << name << brackets[0] << r+base << ", " << c+base << brackets[1] << " = " << A(r, c) << seperator;
		ss << "\n";
	}
	return ss.str();
}

template<typename TStorage>
std::string
MatlabString(const DenseMatrix<TStorage> &A, const char *name)
{
	std::stringstream ss;
	ss << name << " = zeros(" << A.num_rows() << ", " << A.num_cols() << ")\n";
	ss << NonzeroString(A, "()", "\n", 1, name);
	return ss.str();
}


template<typename TStorage>
std::string
JuliaString(const DenseMatrix<TStorage> &A, const char *name)
{
	std::stringstream ss;
	ss << name << " = zeros(" << A.num_rows() << ", " << A.num_cols() << ")\n";
	ss << NonzeroString(A, "[]", "\n", 1, name);
	return ss.str();
}

template<typename TStorage>
std::string
CPPString(const DenseMatrix<TStorage> &A, const char *name)
{
	std::stringstream ss;
	ss << name << " = Matrix(" << A.num_rows() << ", " << A.num_cols() << ")\n";
	ss << NonzeroString(A, "()", ";\n", 0, name);
	return ss.str();
}

}
#endif /* SMALL_ALGEBRA_PRINT_H_ */
