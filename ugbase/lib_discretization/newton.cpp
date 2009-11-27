/*
 * newton.c
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */
 
#include "newton.h"

namespace ug{

NewtonSolver::NewtonSolver(LinearSolver& LinearSolver, int MaxIterations, number MinDefect)
{
	m_LinearSolver = &LinearSolver;
	m_MaxIterations = MaxIterations;
	m_MinDefect = MinDefect;
}

bool NewtonSolver::solve(NumericalSolution& u, NonLinearAssemble* NLAssemble)
{
	number r0, r;
	Vector d;  // defect
	Vector c;  // correction
	//Matrix J;  // Jacobian

	int nDoF = u.GridVector()->length();

	c.create_vector(nDoF);
	d.create_vector(nDoF);
	//J.create_matrix(nDoF,nDoF);

	//verbose
	std::cout << "###### NEWTON - METHOD ######" << std::endl;

	// Set Dirichlet - Nodes to exact values
	if(NLAssemble->assemble_solution(u) == false) return false;

	// Compute first Defect
	d.set(0.0);
	if(NLAssemble->assemble_defect(d, u) == false) return false;

	// Compute first Residuum
	r0 = d.norm2();

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
		c.set(0.0);

		// HYPRE can not change matrix pattern otherwise. Therefore this ugly reallocating. Remove this!
		Matrix* J = new Matrix();
		J->create_matrix(nDoF,nDoF);

		// Compute Jacobian
		if(NLAssemble->assemble_jacobian(*J, u) == false) return false;

		if(i == 0)
		{
			const char MatrixFile[] = "Jacobian";
			const char RHSFile[] = "Defect";
			J->printToFile(MatrixFile);
			d.printToFile(RHSFile);
		}

		// Solve Linearized System
		if(m_LinearSolver->solve(*J, c, d) == false) return false;

		delete J;

		if(i == 0)
		{
			const char CorrectionFile[] = "Correction";
			c.printToFile(CorrectionFile);
		}

		// update Solution
		*u.GridVector() -= c;

		if(i == 0)
		{
			const char UNachherFile[] = "UNachher";
			u.GridVector()->printToFile(UNachherFile);
		}

		// compute new Defect
		d.set(0.0);
		if(NLAssemble->assemble_defect(d, u) == false) return false;

		//compute new Residuum
		r = d.norm2();

		// verbose new defect
		std::cout << i+1 << ". Iteration: defect = " << r << std::endl;

		// finish loop if MinDefect is reached
		if(r < m_MinDefect) break;

		// if MaxIterations reached without success -> abort
		if(i == m_MaxIterations - 1)
		{
			std::cout << "NewtonSolver: Maximum number of iterations reached. Aborting. " << std::endl;
			return false;
		}
	}

	std::cout << "NewtonSolver: System solved after " << i + 1 << " Iteration(s) with Residuum r = " << r << std::endl;
	return true;
}

}