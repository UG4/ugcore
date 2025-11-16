/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "common/util/histogramm.h"
namespace ug{
void checksub(const CPUAlgebra::matrix_type &A)
{
	if(A.num_rows() == 0)
	{
		UG_LOG("EMPTY MATRIX!\n");
		return;
	}
	A.print("A");
	UG_LOG(reset_floats);
	size_t N = A.num_rows();
	// check isolated
	size_t iIsolated = 0;
	for(size_t r=0; r<A.num_rows(); r++)
		if(A.is_isolated(r)) iIsolated++;

	using row_it = CPUAlgebra::matrix_type::const_row_iterator;
	// using value_type = CPUAlgebra::matrix_type::value_type;

	UG_LOG("Nr of dirichlet nodes: " << iIsolated << " (" << iIsolated*100.0/N << "% )\n");

	// check symmetric

	std::vector<double> alpha(N, 0);

	size_t iUnsymmetric=0;
	const double unsymmetricEps = 1e-8;
	for(size_t r=0; r<A.num_rows(); r++)
	{
		double dUnsymmetric = 0;
		if(A.is_isolated(r) ) continue;

		for(row_it it = A.begin_row(r); it != A.end_row(r); ++it)
		{
			size_t c = it.index();
			if(A.is_isolated(c) ) continue;
			double T = A(c, r);
			if(A(c, r) != MatrixTranspose(it.value()))
			{
				dUnsymmetric += dabs(T-it.value());
			}

		}
		double diag = A(r, r);
		if(diag == 0.0) continue;

		dUnsymmetric /= diag;

		alpha[r] = dUnsymmetric;
		if(dUnsymmetric > unsymmetricEps)
			iUnsymmetric++;

	}
	std::sort(alpha.begin(), alpha.end());

	// check sign condition

	if(iUnsymmetric==0)
	{	UG_LOG("Matrix is symmetric! (maximum alpha = " << alpha[N-1] << ")\n"); }
	else
	{
		UG_LOG("Matrix is unsymmetric in " << (iUnsymmetric*100)/(N-iIsolated) << "% of the non-dirchlet rows (" << iUnsymmetric << " total)\n");
		UG_LOG(" row i alpha-unsymmetric means: sum_{A_{ij} != 0} |A_{ij}-A{ji}| / A_{ii} >= alpha\n")

		UG_LOG("> alpha distribution:\n");
		UG_LOG(DistributionPercentage(alpha));
	}
	size_t signConditionMet=0;
	size_t zeroDiagonal=0;

	double minEW = 1e20;
	double maxEW = 1e-20;

	for(size_t r=0; r<A.num_rows(); r++)
	{
		if(A(r, r)==0.0)
		{
			zeroDiagonal++;
			continue;
		}
		bool bPos = A(r, r) > 0;
		double s=0.0;
		bool bSignCondMet = true;
		for(row_it it = A.begin_row(r); it != A.end_row(r); ++it)
		{
			if(it.index() == r) continue;
			if(bPos)
			{
				if(it.value() > 0)
					bSignCondMet = false;
			}
			else
			{
				if(it.value() < 0)
					bSignCondMet = false;
			}
			s += dabs(it.value());
		}
		if(bSignCondMet && A.is_isolated(r) == false)
			signConditionMet++;

		minEW = std::min(minEW, A(r,r)-s);
		maxEW = std::max(maxEW, A(r,r)+s);
	}

	if(signConditionMet == N-iIsolated)
	{	UG_LOG("Sign condition met in all nodes\n"); }
	else
	{
		UG_LOG("Sign condition met in " << (signConditionMet*100.0)/(N-iIsolated) << "% of the non-dirichlet rows (" << signConditionMet << " total)\n");
	}
	UG_LOG("Gershgorin Eigenvalues are within [" << minEW << ", " << maxEW << "]\n");

}
}
