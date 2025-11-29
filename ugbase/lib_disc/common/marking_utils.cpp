/*
 * Copyright (c) 2011-2018:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#include "marking_utils.h"

#include "lib_disc/domain.h"
#include "lib_disc/domain_traits.h"
#include "lib_grid/algorithms/geom_obj_util/anisotropy_util.h"
#include "lib_grid/tools/grid_level.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_grid/tools/surface_view.h"


namespace ug {



template <typename TDomain>
void MarkGlobal(SmartPtr<IRefiner> refiner, SmartPtr<TDomain> domain)
{
	using elem_type = typename domain_traits<TDomain::dim>::element_type;
	using const_iterator = typename SurfaceView::traits<elem_type>::const_iterator;

	// get surface view
	SurfaceView sv(domain->subset_handler());

	// loop elements for marking
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; iter != iterEnd; ++iter)
		refiner->mark(*iter, RM_FULL);
}


template <typename TDomain>
void MarkSubsets
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	const std::vector<std::string>& vSubset
)
{
	using elem_type = typename domain_traits<TDomain::dim>::element_type;
	using const_iterator = typename SurfaceView::traits<elem_type>::const_iterator;

	// get subset handler
	SmartPtr<MGSubsetHandler> sh = domain->subset_handler();

	// transform subset names to indices
	std::vector contained(sh->num_subsets(), false);
	{
		SubsetGroup ssg(sh);
		try {ssg.add(vSubset);}
		UG_CATCH_THROW("MarkSubsets failed to add subsets.");

		const size_t nSs = ssg.size();
		for (size_t i = 0; i < nSs; ++i)
			contained[ssg[i]] = true;
	}

	// get surface view
	SurfaceView sv(sh);

	// loop elements for marking
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::ALL_BUT_SHADOW_COPY);
	for (; iter != iterEnd; ++iter)
		if (contained[sh->get_subset_index(*iter)])
			refiner->mark(*iter, RM_FULL);
}


template <typename TDomain>
void MarkAlongSurface
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	const std::vector<std::string>& surfaceSubsets,
	const std::vector<std::string>& volumeSubsets
)
{
	using elem_type = typename domain_traits<TDomain::dim>::element_type;
	using elem_list_type = typename MultiGrid::traits<elem_type>::secure_container;
	using iter_type = SurfaceView::traits<Vertex>::const_iterator;

	const size_t nSurf = surfaceSubsets.size();
	UG_COND_THROW(volumeSubsets.size() != nSurf, "Same number of surface and volume subsets required.");

	Grid* grid = refiner->grid();
	SmartPtr<MGSubsetHandler> spSH = domain->subset_handler();
	const SurfaceView sv(spSH);

	for (size_t i = 0; i < nSurf; ++i)
	{
		int surf_si;
		int vol_si;
		try
		{
			surf_si = spSH->get_subset_index(surfaceSubsets[i].c_str());
			vol_si = spSH->get_subset_index(volumeSubsets[i].c_str());
		}
		UG_CATCH_THROW("Cannot convert subset names " << surfaceSubsets[i] << " and "
			<< volumeSubsets[i] << " to indices.");

		//	loop elements for marking
		iter_type iter = sv.begin<Vertex>(surf_si, GridLevel(), SurfaceView::SurfaceConstants::MG_ALL);
		iter_type iterEnd = sv.end<Vertex>(surf_si, GridLevel(), SurfaceView::SurfaceConstants::MG_ALL);
		for (; iter != iterEnd; ++iter)
		{
			Vertex* v = *iter;
			elem_list_type el;
			grid->associated_elements(el, v);
			size_t el_sz = el.size();
			for (size_t j = 0; j < el_sz; ++j)
			{
				int elem_si = spSH->get_subset_index(el[j]);
				if (elem_si == vol_si && !sv.surface_state(el[j]).contains(SurfaceView::SurfaceConstants::MG_SHADOW_PURE))
				{
					// mark the element for anisotropic refinement
					refiner->mark(el[j], RefinementMark::RM_REFINE);
				}
			}
		}
	}
}


template <typename TDomain>
void MarkAnisotropic
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
)
{
	using elem_type = typename domain_traits<TDomain::dim>::element_type;
	using side_type = typename domain_traits<TDomain::dim>::side_type;
	using const_iterator = typename SurfaceView::traits<elem_type>::const_iterator;

	Grid& grid = *refiner->grid();
	typename TDomain::position_accessor_type aaPos = domain->position_accessor();

	// get surface view and prepare loop over all surface elements
	SurfaceView sv(domain->subset_handler());
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::SurfaceConstants::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::SurfaceConstants::ALL_BUT_SHADOW_COPY);

	// loop elements for marking
	std::vector<Edge*> longEdges;
	for (; iter != iterEnd; ++iter)
	{
		longEdges.clear();
		AnisotropyState state = long_edges_of_anisotropic_elem(*iter, grid, aaPos, thresholdRatio, longEdges);
		if (state != AnisotropyState::ISOTROPIC)
		{
			// mark elem
			refiner->mark(*iter, RefinementMark::RM_CLOSURE);

			// mark long edges
			const size_t nEdges = longEdges.size();
			for (size_t e = 0; e < nEdges; ++e)
				refiner->mark(longEdges[e], RefinementMark::RM_FULL);

			// mark all sides
			typename Grid::traits<side_type>::secure_container sl;
			grid.associated_elements(sl, *iter);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
				if (refiner->get_mark(sl[s]) != RefinementMark::RM_FULL)
					refiner->mark(sl[s], RefinementMark::RM_CLOSURE);
		}
	}
}


template <typename TDomain>
void MarkAnisotropicOnlyX
(
	SmartPtr<IRefiner> refiner,
	SmartPtr<TDomain> domain,
	number thresholdRatio
)
{
	using elem_type = typename domain_traits<TDomain::dim>::element_type;
	using side_type = typename domain_traits<TDomain::dim>::side_type;
	using const_iterator = typename SurfaceView::traits<elem_type>::const_iterator;

	Grid& grid = *refiner->grid();
	typename TDomain::position_accessor_type aaPos = domain->position_accessor();

	// get surface view and prepare loop over all surface elements
	SurfaceView sv(domain->subset_handler());
	const_iterator iter = sv.begin<elem_type>(GridLevel(), SurfaceView::SurfaceConstants::ALL_BUT_SHADOW_COPY);
	const_iterator iterEnd = sv.end<elem_type>(GridLevel(), SurfaceView::SurfaceConstants::ALL_BUT_SHADOW_COPY);

	// loop elements for marking
	std::vector<Edge*> longEdges;
	for (; iter != iterEnd; ++iter)
	{
		longEdges.clear();
		AnisotropyState state = long_edges_of_anisotropic_elem(*iter, grid, aaPos, thresholdRatio, longEdges);
		if (state == AnisotropyState::ISOTROPIC)
			continue;

		UG_COND_THROW(!longEdges.size(), "Element is anisotropic, but no long edges present.");

		// check whether edges point in x direction
		Edge* longEdge = longEdges[0];
		MathVector<TDomain::dim> dir;
		VecSubtract(dir, aaPos[longEdge->vertex(1)], aaPos[longEdge->vertex(0)]);
		VecNormalize(dir, dir);
		if (fabs(dir[0]) > 0.9)
		{
			// mark elem
			refiner->mark(*iter, RefinementMark::RM_CLOSURE);

			// mark long edges
			const size_t nEdges = longEdges.size();
			for (size_t e = 0; e < nEdges; ++e)
				refiner->mark(longEdges[e], RefinementMark::RM_FULL);

			// mark all sides
			typename Grid::traits<side_type>::secure_container sl;
			grid.associated_elements(sl, *iter);
			const size_t slSz = sl.size();
			for (size_t s = 0; s < slSz; ++s)
				if (refiner->get_mark(sl[s]) != RefinementMark::RM_FULL)
					refiner->mark(sl[s], RefinementMark::RM_CLOSURE);
		}
	}
}



// explicit template specializations
#ifdef UG_DIM_1
	template void MarkGlobal<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>);
	template void MarkSubsets<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>, const std::vector<std::string>&);
	template void MarkAlongSurface(SmartPtr<IRefiner>, SmartPtr<Domain1d>, const std::vector<std::string>&,
		const std::vector<std::string>&);
	template void MarkAnisotropic<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>, number);
	template void MarkAnisotropicOnlyX<Domain1d>(SmartPtr<IRefiner>, SmartPtr<Domain1d>, number);
#endif
#ifdef UG_DIM_2
	template void MarkGlobal<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>);
	template void MarkSubsets<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>, const std::vector<std::string>&);
	template void MarkAlongSurface(SmartPtr<IRefiner>, SmartPtr<Domain2d>, const std::vector<std::string>&,
		const std::vector<std::string>&);
	template void MarkAnisotropic<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>, number);
	template void MarkAnisotropicOnlyX<Domain2d>(SmartPtr<IRefiner>, SmartPtr<Domain2d>, number);
#endif
#ifdef UG_DIM_3
	template void MarkGlobal<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>);
	template void MarkSubsets<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>, const std::vector<std::string>&);
	template void MarkAlongSurface(SmartPtr<IRefiner>, SmartPtr<Domain3d>, const std::vector<std::string>&,
		const std::vector<std::string>&);
	template void MarkAnisotropic<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>, number);
	template void MarkAnisotropicOnlyX<Domain3d>(SmartPtr<IRefiner>, SmartPtr<Domain3d>, number);
#endif


} // end namespace ug
