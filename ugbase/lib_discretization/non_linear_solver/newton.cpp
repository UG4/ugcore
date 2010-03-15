/*
 * newton.c
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#include "newton.h"

namespace ug{

NewtonSolver::NewtonSolver(ILinearSolver& LinearSolver, int MaxIterations, number MinDefect, bool reallocate)
{
	m_LinearSolver = &LinearSolver;
	m_MaxIterations = MaxIterations;
	m_MinDefect = MinDefect;
	_reallocate = reallocate;
	_allocated = false;
}

NewtonSolver::~NewtonSolver()
{
	if(_allocated)
	{
		delete c; delete d; delete J; _allocated = false;
	}
}


}



