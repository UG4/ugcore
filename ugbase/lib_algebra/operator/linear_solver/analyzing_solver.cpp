#include "lib_algebra/lib_algebra.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "common/util/histogramm.h"
namespace ug{
void checksub(const CPUAlgebra::matrix_type &A)
{
	UG_LOG(reset_floats);
	size_t N = A.num_rows();
	// check isolated
	size_t iIsolated = 0;
	for(size_t r=0; r<A.num_rows(); r++)
		if(A.is_isolated(r)) iIsolated++;

	typedef CPUAlgebra::matrix_type::const_row_iterator row_it;
	typedef CPUAlgebra::matrix_type::value_type value_type;

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
