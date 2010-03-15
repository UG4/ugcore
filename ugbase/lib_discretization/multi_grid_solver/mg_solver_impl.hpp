/*
 * mg_solver_impl.h
 *
 *  Created on: 04.01.2010
 *      Author: andreasvogel
 */

#ifndef MG_SOLVER_IMPL_H_
#define MG_SOLVER_IMPL_H_

namespace ug{

template <int dim>
MultiGridSolver<dim>::
MultiGridSolver(IAssemble<dim>& ass, Domain<dim>& domain, uint surfaceLevel, uint baseLevel, NumericalSolution<dim>& u,
		int maxIter, number tol, int cycle, LinearSolver& Smoother, number damp, int nu1, int nu2, LinearSolver& BaseSolver)
{
	_ass = &ass;
	_domain = &domain;
	_u = &u;

	_maxIter = maxIter;
	_tol = tol;

	_Smoother = &Smoother;
	_damp = damp;
	_nu1 = nu1;
	_nu2 = nu2;
	_num_cycle = cycle;

	_BaseSolver = &BaseSolver;
	_baseLevel = baseLevel;
	_surfaceLevel = surfaceLevel;

}


template <int dim>
bool
MultiGridSolver<dim>::
lmgc(uint l)
{
	if(l > _baseLevel)
	{
		// Presmooth
		for(int i = 0; i < _nu1; ++i)
		{
			_t[l].set(0.0);
			_Smoother->step(*_A[l], _t[l], _d[l], _damp);
			_c[l] += _t[l];
		}

		// Restrict Defect
		_I[l-1]->applyTransposed(_d[l-1], _d[l]);

		// reset correction
		_c[l-1].set(0.0);

		// apply lmgc on coarser grid
		for(int i = 0; i < _num_cycle; ++i)
		{
			if(lmgc(l-1) == false)
			{
				std::cout << "Error in lmgc on level " << l-1 << "." << std::endl;
				return false;
			}
		}

		//interpolate correction
		_I[l-1]->apply(_t[l], _c[l-1]);

		// update correction
		_c[l] += _t[l];

		//update defect
		_A[l]->matmul_minus(_d[l], _t[l]);

		// Postsmooth
		for(int i = 0; i < _nu2; ++i)
		{
			_t[l].set(0.0);
			_Smoother->step(*_A[l], _t[l], _d[l], _damp);
			_c[l] += _t[l];
		}

		return true;
	}
	else if(l == _baseLevel)
	{
		// exact solve
		_t[l].set(0.0);
		_BaseSolver->solve(*_A[l], _t[l], _d[l]);
		_c[l] += _t[l];
		return true;
	}
	else
	{
		std::cout << "Level index below 'baseLevel' in lmgc. ERROR." << std::endl;
		return false;
	}
}


template <int dim>
bool
MultiGridSolver<dim>::
solve(Matrix& A, Vector& x, Vector& b)
{
	if(_domain->get_grid_type() == GT_GRID)
	{
		std::cout << "Can not solve on a grid. Please use MultiGrid structure.\n";
		return false;
	}

	MultiGrid& mg = *_domain->get_multigrid();

	if(_surfaceLevel >= mg.num_levels())
	{
		std::cout << "SurfaceLevel " << _surfaceLevel << " does not exist" << std::endl;
		return false;
	}
	if(_baseLevel > _surfaceLevel)
	{
		std::cout << "Baselevel must be smaller than to equal to SurfaceLevel." << std::endl;
		return false;
	}

	// allocate Matrices and Vectors (currently to much, should be baseLevel to surfaceLevel+1 only)
	_A = new Matrix*[_surfaceLevel+1];
	_I = new Matrix*[_surfaceLevel];

	_c = new Vector[_surfaceLevel+1];
	_d = new Vector[_surfaceLevel+1];
	_t = new Vector[_surfaceLevel+1];

	// create matrices and vectors for coarser grids
	for(uint i = _baseLevel; i < _surfaceLevel; ++i)
	{
		_A[i] = new Matrix;
		_A[i]->create_matrix(_u->num_dofs(i), _u->num_dofs(i));
		_c[i].create_vector(_u->num_dofs(i));
		_d[i].create_vector(_u->num_dofs(i));
		_t[i].create_vector(_u->num_dofs(i));

		_d[i].set(0.0); //TODO: Needed ???
	}

	// set Matrix for surface Level
	_A[_surfaceLevel] = &A;

	// create defect, correction and help vector on surface Level
	_c[_surfaceLevel].create_vector(_u->num_dofs(_surfaceLevel));
	_d[_surfaceLevel].create_vector(_u->num_dofs(_surfaceLevel));
	_t[_surfaceLevel].create_vector(_u->num_dofs(_surfaceLevel));

	// assemble matrices on coarser grids
	for(uint i = _baseLevel; i < _surfaceLevel; ++i)
	{
		if(_ass->assemble_jacobian(*_A[i], *_u, i) != IAssemble_OK)
		{
			std::cout << " Error in assemble_stiffness matrix, aborting." << std::endl;
			for(uint j = _baseLevel; j < _surfaceLevel; ++j) {delete _A[j];} delete[] _A;delete[] _d;delete[] _c;delete[] _t;
			return false;
		}
	}

	// create interpolation matrices
	for(uint i = _baseLevel; i < _surfaceLevel; ++i)
	{
		_I[i] = new Matrix;
		_I[i]->create_matrix(_u->num_dofs(i+1), _u->num_dofs(i));

		if(assemble_interpolation(*_I[i], *_u->get_pattern(),
					mg.begin<Vertex>(i+1), mg.end<Vertex>(i+1),
					i, i+1) == false)
		{
			std::cout << "Error in assemble_interpolation, aborting." << std::endl;
			for(uint j = _baseLevel; j < _surfaceLevel; ++j) {delete _A[j]; delete _I[j];} delete[] _A;delete[] _d;delete[] _c;delete[] _t; delete[] _I;
			return false;
		}
	}

	// copy rhs
	_d[_surfaceLevel].set(0.0);
	_d[_surfaceLevel] += b;

	// build defect:  d := b - A*x
	_A[_surfaceLevel]->matmul_minus(_d[_surfaceLevel], x);

	// compute start norm ||d||_2
	number norm, norm_old;
	norm = norm_old = _d[_surfaceLevel].norm2();

	// Print Start information
	std::cout << "\n   ######### MGS #########\n";
	std::cout << "  Iter     Defect         Rate \n";
	std::cout << std::setw(4) << 0 << ":  " << std::scientific << norm_old <<  "      -------" << std::endl;

	// Iteration loop
	for(int i = 1; i <= _maxIter; ++i)
	{
		// check for convergence
		if(norm < _tol)
		{
			std::cout << "\n ##### Tolerance " << _tol << " reached. MultiGrid Solver converged.\n" << std::endl;
			for(uint j = _baseLevel; j < _surfaceLevel; ++j) {delete _A[j]; delete _I[j];} delete[] _A;delete[] _d;delete[] _c;delete[] _t; delete[] _I;
			return true;
		}

		// reset correction
		_c[_surfaceLevel].set(0.0);

		// perform one multigrid cycle
		if(lmgc(_surfaceLevel)==false)
		{
			std::cout << "Error in lmgc. Aborting." << std::endl;
			for(uint j = _baseLevel; j < _surfaceLevel; ++j) {delete _A[j]; delete _I[j];} delete[] _A;delete[] _d;delete[] _c;delete[] _t; delete[] _I;
			return false;
		}

		// add correction to solution
		x += _c[_surfaceLevel];

		// compute new defect norm
		norm = _d[_surfaceLevel].norm2();

		// print convergence rate
		std::cout << std::setw(4) << i << ":  " << std::scientific << norm << "    " <<norm/norm_old << std::endl;

		// remember current norm
		norm_old = norm;
	}

	std::cout << "\n ##### Tolerance " << _tol << " NOT reached after " << _maxIter << " Iterations. MultiGrid Solver did NOT CONVERGE.\n" << std::endl;
	for(uint j = _baseLevel; j < _surfaceLevel; ++j) {delete _A[j]; delete _I[j];} delete[] _A;delete[] _d;delete[] _c;delete[] _t; delete[] _I;
	return false;
}

template <int dim>
bool
MultiGridSolver<dim>::
print_matrices()
{
	if(_domain->get_grid_type() == GT_GRID)
	{
		std::cout << "Can not solve on a grid. Please use MultiGrid structure.\n";
		return false;
	}

	MultiGrid& mg = *_domain->get_multigrid();

	if(_surfaceLevel >= mg.num_levels())
	{
		std::cout << "SurfaceLevel " << _surfaceLevel << " does not exist" << std::endl;
		return false;
	}
	if(_baseLevel > _surfaceLevel)
	{
		std::cout << "Baselevel must be smaller than to equal to SurfaceLevel." << std::endl;
		return false;
	}

	_A = new Matrix*[_surfaceLevel];
	_I = new Matrix*[_surfaceLevel];

	// create matrices and vectors for coarser grids
	for(uint i = _baseLevel; i < _surfaceLevel; ++i)
	{
		_A[i] = new Matrix;
		_A[i]->create_matrix(_u->num_dofs(i), _u->num_dofs(i));
	}

	// assemble matrices on coarser grids
	for(uint i = _baseLevel; i < _surfaceLevel; ++i)
	{
		if(_ass->assemble_jacobian(*_A[i], *_u, i) != IAssemble_OK)
		{
			std::cout << " Error in assemble_stiffness matrix, aborting." << std::endl;
			for(uint j = _baseLevel; j < _surfaceLevel; ++j) delete _A[j]; delete[] _A;delete[] _I;
			return false;
		}

		char name[30];
        sprintf(name, "Matrix_Level_%i.txt", i);
		_A[i]->printToFile(name);
	}

	// create interpolation matrices
	for(uint i = _baseLevel; i < _surfaceLevel; ++i)
	{
		_I[i] = new Matrix;
		_I[i]->create_matrix(_u->num_dofs(i+1), _u->num_dofs(i));

		if(assemble_interpolation(*_I[i], *_u->get_pattern(),
					mg.begin<Vertex>(i+1), mg.end<Vertex>(i+1),
					i, i+1) == false)
		{
			std::cout << "Error in assemble_interpolation, aborting." << std::endl;
			for(uint j = _baseLevel; j < _surfaceLevel; ++j) delete _A[j]; delete[] _A;delete[] _I;
			return false;
		}

		char name[30];
		sprintf(name, "Interpolation_Matrix_Level_%i.txt", i);
		_I[i]->printToFile(name);
	}

	for(uint j = _baseLevel; j < _surfaceLevel; ++j) delete _A[j]; delete[] _A;delete[] _I;

	return true;

}
}


#endif /* MG_SOLVER_IMPL_H_ */
