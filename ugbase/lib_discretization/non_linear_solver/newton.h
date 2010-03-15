/*
 * newton.h
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__NEWTON__
#define __H__LIBDISCRETIZATION__NEWTON__

// other ug libraries
#include "../lib_algebra/lib_algebra.h"

// modul intern headers
#include "numericalsolution.h"
#include "assemble.h"

namespace ug {

class NewtonSolver {

	public:
		NewtonSolver(ILinearSolver& LinearSolver, int MaxIterations, number MinDefect, bool reallocate);

		template <int dim>
		bool solve(NumericalSolution<dim>& u, IAssemble<dim>& ass, uint level = 0);

		~NewtonSolver();

	private:
		ILinearSolver* m_LinearSolver;
		int m_MaxIterations;
		number m_MinDefect;
		Vector* d;
		Vector* c;
		Matrix* J;

		bool _reallocate;
		bool _allocated;
};


template <int dim>
bool NewtonSolver::solve(NumericalSolution<dim>& u, IAssemble<dim>& ass, uint level)
{
	number r0, r = 0.0;

	if(_reallocate || !(_allocated))
	{
		d = new Vector();  // defect
		c = new Vector();  // correction
		J = new Matrix();  // Jacobian

		int nDoF = u.GridVector(level)->length();

		c->create_vector(nDoF);
		d->create_vector(nDoF);
		J->create_matrix(nDoF,nDoF);
		_allocated = true;
	}

	//verbose
	std::cout << "###### NEWTON - METHOD ######" << std::endl;

	// Set Dirichlet - Nodes to exact values
	if(ass.assemble_solution(u, level) != IAssemble_OK) return false;

	// Compute first Defect
	d->set(0.0);
	if(ass.assemble_defect(*d, u, level) != IAssemble_OK) return false;

	const char RHSFile[] = "Defect";
	d->printToFile(RHSFile);

	// Compute first Residuum
	r0 = d->norm2();

	// verbose first Defect
	std::cout << "Start Defect: " << r0 << std::endl;
	if(r0 < m_MinDefect)
	{
		std::cout << "Initial defect smaller than desired defect. Finish newton method." << std::endl;
		return true;
	}

	//loop iteration
	int i;
	for(i = 0; i < m_MaxIterations; ++i)
	{
		// set c = 0
		c->set(0.0);

		// reset Jacobian
		J->set(0.0);

		// Compute Jacobian
		if(ass.assemble_jacobian(*J, u, level) != IAssemble_OK) return false;

		if(i == 0)
		{
			const char MatrixFile[] = "Jacobian";
			J->printToFile(MatrixFile);
		}

		// Solve Linearized System
		if(m_LinearSolver->solve(*J, *c, *d) == false) return false;

		if(i == 0)
		{
			const char CorrectionFile[] = "Correction";
			c->printToFile(CorrectionFile);
		}

		if(i == 0)
		{
			const char UNachherFile[] = "UVorher";
			u.GridVector(level)->printToFile(UNachherFile);
		}

		// update Solution
		*u.GridVector(level) -= *c;

		if(i == 0)
		{
			const char UNachherFile[] = "UNachher";
			u.GridVector(level)->printToFile(UNachherFile);
		}

		// compute new Defect
		d->set(0.0);
		if(ass.assemble_defect(*d, u, level) != IAssemble_OK) return false;

		//compute new Residuum
		r = d->norm2();

		// verbose new defect
		std::cout << i+1 << ". Iteration: defect = " << r << std::endl;

		// finish loop if MinDefect is reached
		if(r < m_MinDefect) break;

		// if MaxIterations reached without success -> abort
		if(i == m_MaxIterations - 1)
		{
			std::cout << "NewtonSolver: Maximum number of iterations reached. Aborting. " << std::endl;
			if(_reallocate)
			{
				delete c; delete d; delete J; _allocated = false;
			}
			return false;
		}
	}

	std::cout << "NewtonSolver: System solved after " << i + 1 << " Iteration(s) with Residuum r = " << r << std::endl;
	if(_reallocate)
	{
		delete c; delete d; delete J; _allocated = false;
	}
	return true;
}

}

#endif /* __H__LIBDISCRETIZATION__NEWTON__ */
