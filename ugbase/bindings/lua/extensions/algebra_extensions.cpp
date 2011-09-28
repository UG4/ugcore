/*
 * algebra_extensions.cpp
 *
 *  Created on: 04.03.2011
 *      Author: kosta
 */

/**	This file is only temporary. Methods which are declared here
 * should be reworked and moved into the main ug-section.
 * Please take care to not introduce unwanted dependencies.
 */
#include "ug.h"
#include "algebra_extensions.h"
#include "ug_bridge/ug_bridge.h"
#include "lib_algebra/algebra_types.h"

#include "bindings/lua/user_data/user_data.h"

#include "lib_disc/dof_manager/p1conform/p1conform.h"
#include "lib_disc/spatial_discretization/disc_util/finite_volume_geometry.h"
#include "lib_disc/common/groups_util.h"

#include <iostream>
#include <sstream>

namespace ug
{

namespace bridge
{

template <typename TVector>
void KostaUpdate(TVector& xOut, const TVector& xOld, const TVector& v_m, const number dt,
			const LuaUserNumberNumberFunction& alpha, const LuaUserNumberNumberFunction& beta)
{
	UG_ASSERT(xOut.size() == xOld.size(), "Vector size does not match with first vector.");

	UG_ASSERT(xOut.size() == v_m.size(), "Vector size does not match with second vector.");


	for(size_t i = 0; i < xOut.size(); ++i)
	{
		const number V_m = v_m[i];
		const number x_old = xOld[i];

		xOut[i] = (x_old + dt*alpha(1, V_m))/(1 + dt*(alpha(1, V_m) + beta(1, V_m)));

	}
}


template <	typename TVector,
			typename TDoFDistribution,
			typename TDomain>
void transmembrane_current_as_vector(TVector& transmembrane_current, const TVector& membrane_potential_tn,
									const TVector& membrane_potential_tn1,
									const TVector& n_gate, const TVector& m_gate, const TVector& h_gate,
									const number dt,
									const IDoFDistribution<TDoFDistribution>& dofDistr,
									const TDomain& domain,
									const char* bndSubsetName)
{
	UG_ASSERT(transmembrane_current.size() == membrane_potential_tn.size(), "Vector size does not match with first vector.");
	UG_ASSERT(transmembrane_current.size() == membrane_potential_tn1.size(), "Vector size does not match with first vector.");

	UG_ASSERT(transmembrane_current.size() == n_gate.size(), "Vector size does not match with first vector.");
	UG_ASSERT(transmembrane_current.size() == m_gate.size(), "Vector size does not match with first vector.");
	UG_ASSERT(transmembrane_current.size() == h_gate.size(), "Vector size does not match with first vector.");

//	Reset current
	transmembrane_current.set(0.0);

	const static int dim = TDomain::dim;

//	get number of subset requested as boundary
	SubsetGroup bndSSGrp;
	if(!ConvertStringToSubsetGroup(bndSSGrp, dofDistr.get_function_pattern(),
								   bndSubsetName))
	{
		UG_LOG("ERROR in 'transmembrane_current_as_vector':"
				" Subsets '"<<bndSubsetName<<"' not"
				" all contained in ApproximationSpace.\n");
		throw(UGFatalError("Wrong Subset name."));
	}



	DimFV1Geometry<dim> geo;
	std::vector<MultiIndex<2> > ind;

	SubsetGroup unionSubsets;
	std::vector<SubsetGroup> vSSGrp;

//	get the subset handler
	const ISubsetHandler& sh = *dofDistr.get_function_pattern().get_subset_handler();

//	request the creation of bf of this bndSubset
	for(size_t i = 0; i < bndSSGrp.num_subsets(); ++i)
		geo.add_boundary_subset(bndSSGrp[i]);

//	component of solution vector
	const size_t _C_ = 0;

	typedef typename domain_traits<dim>::geometric_base_object TElem;

//	loop subsets
	for(int si = 0; si < dofDistr.num_subsets(); ++si)
	{
		if(si != 1) continue;
		// run over all elements of the grid restricted to the subset
		typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
		iterBegin = dofDistr.template begin<TElem>(si);
		iterEnd = dofDistr.template end<TElem>(si);
	// 	Loop over all elements
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		// 	get Element
			TElem* elem = *iter;

		// 	get global indices
			dofDistr.multi_indices(elem, _C_, ind);

		//	vector of corner coordinates
			std::vector<MathVector<dim> > vCornerCoords;

		//	get corner coords
			CollectCornerCoordinates(vCornerCoords, *elem, domain, true);

		//	update geom
			geo.update(elem, &vCornerCoords[0], &sh);


			for(size_t i = 0; i < bndSSGrp.num_subsets(); ++i)
				for(size_t k= 0; k<geo.num_bf(bndSSGrp[i]); ++k)
				{
					const typename DimFV1Geometry<dim>::BF& bf = geo.bf(bndSSGrp[i], k);

					const size_t locIndex = bf.node_id(); // loc = 0, ..., 3
					const size_t globIndex = ind[locIndex][0]; // i = 0, ...,  numDoFs
//					const size_t globAlpha = ind[locIndex][1]; // always 0, when only 1 component

					transmembrane_current[globIndex] += bf.volume()*(
								-(membrane_potential_tn1[globIndex] - membrane_potential_tn[globIndex])/dt +
									-120.0*pow(n_gate[globIndex],4)*(membrane_potential_tn1[globIndex]- 50) +
									-36.0*pow(m_gate[globIndex],3)*h_gate[globIndex]*(membrane_potential_tn1[globIndex] + 77) +
									-0.3*(membrane_potential_tn1[globIndex] + 54.4));
				}
		}
	}
}



template <typename TVector>
void HhFlux(TVector& hhFlux, const TVector& nGate, const TVector& mGate, const TVector& hGate,
			const TVector& vm, const TVector& injection)
{
	UG_ASSERT(hhFlux.size() == nGate.size(), "Vector size does not match with first vector.");

	UG_ASSERT(hhFlux.size() == mGate.size(), "Vector size does not match with first vector.");
	UG_ASSERT(hhFlux.size() == hGate.size(), "Vector size does not match with first vector.");
	UG_ASSERT(hhFlux.size() == vm.size(), "Vector size does not match with first vector.");
	UG_ASSERT(hhFlux.size() == injection.size(), "Vector size does not match with first vector.");


	for(size_t i = 0; i < hhFlux.size(); i++)
	{
		const number V_m = vm[i];
		const number n = nGate[i];
		const number m = mGate[i];
		const number h = hGate[i];
		const number inj = injection[i];

		hhFlux[i] = (120*m*m*m*h*(V_m - 50) + 36*n*n*n*n*(V_m + 77) + 0.3*(V_m + 54.4) - inj);
	//	UG_LOG("hhFlux[i]="<<hhFlux[i]);
	}
}


bool RegisterAlgebraExtensions(Registry& reg, const char* parentGroup)
{
//	get group string
	std::stringstream groupString; groupString << parentGroup << "/AlgebraExtensions";
	std::string grp = groupString.str();
	reg.add_function("transmembrane_current_as_vector",&transmembrane_current_as_vector<CPUAlgebra::vector_type,P1DoFDistribution, Domain2d>, grp.c_str() );
	reg.add_function("KostaUpdate", &KostaUpdate<CPUAlgebra::vector_type>, grp.c_str());
	reg.add_function("HhFlux", &HhFlux<CPUAlgebra::vector_type>, grp.c_str());
	return true;
}



}
}
