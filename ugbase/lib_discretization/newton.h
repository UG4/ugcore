/*
 * newton.h
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__NEWTON__
#define __H__LIBDISCRETIZATION__NEWTON__

#include "../lib_algebra/lib_algebra.h"
#include "numericalsolution.h"

namespace ug {

class NonLinearAssemble{

	public:
		virtual bool assemble_jacobian(Matrix& mat, NumericalSolution& u) = 0;
		virtual bool assemble_defect(Vector& vec, NumericalSolution& u) = 0;
		virtual bool assemble_solution(NumericalSolution& u) = 0;
};


class NewtonSolver {
	protected:
		typedef bool (*JacobianFct)(Matrix&, NumericalSolution&);
		typedef bool (*DefectFct)(Vector&, NumericalSolution&);

	public:
		NewtonSolver(LinearSolver& LinearSolver, int MaxIterations, number MinDefect);

		bool solve(NumericalSolution& u, NonLinearAssemble* NLAssemble);

	private:
		LinearSolver* m_LinearSolver;
		int m_MaxIterations;
		number m_MinDefect;

};

}

#endif /* __H__LIBDISCRETIZATION__NEWTON__ */
