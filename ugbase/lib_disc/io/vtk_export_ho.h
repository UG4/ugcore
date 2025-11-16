/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__IO__VTK_EXPORT_HO__
#define __H__UG__LIB_DISC__IO__VTK_EXPORT_HO__

#include <cmath>                                                                   // for ceil, log2
#include <cstddef>                                                                 // for size_t
#include <limits>                                                                  // for numeric_limits
#include <ostream>                                                                 // for string, operator<<, basic_ostream, endl
#include <string>                                                                  // for char_traits, allocator, operator+, basic_string
#include <vector>                                                                  // for vector

#include "common/assert.h"                                                         // for UG_ASSERT
#include "common/error.h"                                                          // for UG_CATCH_THROW, UG_COND_THROW
#include "common/types.h"                                                          // for number
#include "common/math/math_vector_matrix/math_vector.h"                            // for MathVector
#include "common/util/smart_pointer.h"                                             // for SmartPtr, ConstSmartPtr, make_sp
#include "lib_disc/common/function_group.h"                                        // for FunctionGroup
#include "lib_disc/common/multi_index.h"                                           // for DoFIndex, DoFRef
#include "lib_disc/dof_manager/dof_distribution.h"                                 // for DoFDistribution, DoFDistribution::traits
#include "lib_disc/function_spaces/dof_position_util.h"                            // for InnerDoFPosition
#include "lib_disc/function_spaces/grid_function_global_user_data.h"               // for GlobalGridFunctionNumberData
#include "lib_disc/io/vtkoutput.h"                                                 // for VTKOutput
#include "lib_disc/local_finite_element/local_finite_element_id.h"                 // for LFEID
#include "lib_grid/algorithms/debug_util.h"                                        // for ElementDebugInfo
#include "lib_grid/algorithms/selection_util.h"                                     // for SelectAssociatedGridObjects
#include "lib_grid/attachments/attachment_pipe.h"                                  // for AttachmentAccessor
#include "lib_grid/common_attachments.h"                                           // for AVertex
#include "lib_grid/grid/grid.h"                                                    // for Grid::VertexAttachmentAccessor::VertexAttachmentAccessor<TAt...
#include "lib_grid/grid/grid_base_object_traits.h"                                 // for VertexIterator
#include "lib_grid/grid/grid_base_objects.h"                                       // for CustomVertexGroup, Vertex, Edge (ptr only), Face (ptr only)
#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/parallel_refinement/parallel_refinement.h"  // for ParallelGlobalRefiner_MultiGrid
#else
	#include "lib_grid/refinement/global_multi_grid_refiner.h"                     // for GlobalMultiGridRefiner
#endif
#include "lib_grid/multi_grid.h"                                                   // for MultiGrid
#include "lib_grid/tools/selector_grid.h"                                          // for Selector
#include "lib_grid/tools/subset_group.h"                                           // for SubsetGroup


namespace ug {


template <typename TDomain, class TElem>
inline void CopySelectedElements
(
	SmartPtr<TDomain> destDom,
	SmartPtr<TDomain> srcDom,
	Selector& sel,
	AVertex& aNewVrt
)
{
	using SH_type = typename TDomain::subset_handler_type;

	MultiGrid& srcGrid = *srcDom->grid();
	MultiGrid& destGrid = *destDom->grid();

	ConstSmartPtr<SH_type> srcSH = srcDom->subset_handler();
	SmartPtr<SH_type> destSH = destDom->subset_handler();

	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	CustomVertexGroup vrts;
	using iter_t = typename Grid::traits<TElem>::iterator;

	for (iter_t eiter = sel.begin<TElem>();	eiter != sel.end<TElem>(); ++eiter)
	{
		TElem* e = *eiter;
		vrts.resize(e->num_vertices());
		for (size_t iv = 0; iv < e->num_vertices(); ++iv)
			vrts.set_vertex(iv, aaNewVrt[e->vertex(iv)]);

		TElem* ne;
		try {ne = *destGrid.create_by_cloning(e, vrts);}
		UG_CATCH_THROW("New element could not be created.");
		destSH->assign_subset(ne, srcSH->get_subset_index(e));
	}
}

template <typename TDomain>
inline void CopySelected
(
	SmartPtr<TDomain> destDom,
	SmartPtr<TDomain> srcDom,
	Selector& sel
)
{
	using APos_type = typename TDomain::position_attachment_type;
	using SH_type = typename TDomain::subset_handler_type;

	MultiGrid& srcGrid = *srcDom->grid();
	MultiGrid& destGrid = *destDom->grid();

	ConstSmartPtr<SH_type> srcSH = srcDom->subset_handler();
	SmartPtr<SH_type> destSH = destDom->subset_handler();

	APos_type& aPos = srcDom->position_attachment();
	UG_COND_THROW(!srcGrid.has_vertex_attachment(aPos), "Position attachment required.");
	Grid::VertexAttachmentAccessor<APos_type> aaPosSrc(srcGrid, aPos);
	if (!destGrid.has_vertex_attachment(aPos))
		destGrid.attach_to_vertices(aPos);
	Grid::VertexAttachmentAccessor<APos_type> aaPosDest(destGrid, aPos);

	AVertex aNewVrt;
	srcGrid.attach_to_vertices(aNewVrt);
	Grid::VertexAttachmentAccessor<AVertex> aaNewVrt(srcGrid, aNewVrt);

	for (int si = destSH->num_subsets(); si < srcSH->num_subsets(); ++si)
		destSH->subset_info(si) = srcSH->subset_info(si);

	SelectAssociatedGridObjects(sel);

	// create new vertices in destGrid
	for (VertexIterator viter = sel.begin<Vertex>(); viter != sel.end<Vertex>(); ++viter)
	{
		Vertex* v = *viter;
		Vertex* nv;
		try {nv = *destGrid.create_by_cloning(v);}
		UG_CATCH_THROW("New vertex could not be created.");
		aaNewVrt[v] = nv;
		aaPosDest[nv] = aaPosSrc[v];
		destSH->assign_subset(nv, srcSH->get_subset_index(v));
	}

	CopySelectedElements<TDomain,Edge>(destDom, srcDom, sel, aNewVrt);
	CopySelectedElements<TDomain,Face>(destDom, srcDom, sel, aNewVrt);
	CopySelectedElements<TDomain,Volume>(destDom, srcDom, sel, aNewVrt);

	srcGrid.detach_from_vertices(aNewVrt);
}


template <typename TElem, typename TGridFunction, typename TGGFND>
inline void interpolate_from_original_fct
(
	SmartPtr<TGridFunction> u_new,
	const TGGFND& u_orig,
	size_t fct,
	const LFEID& lfeid
)
{
	using dom_type = typename TGridFunction::domain_type;
	static constexpr int dim = dom_type::dim;
	using const_iter_type = typename DoFDistribution::traits<TElem>::const_iterator;

	// interpolate subset-wise
	size_t numSubsets = u_new->num_subsets();
	for (size_t si = 0; si < numSubsets; ++si)
	{
		if (!u_new->is_def_in_subset(fct, si)) continue;

		const_iter_type elem_iter = u_new->template begin<TElem>(si);
		const_iter_type iterEnd = u_new->template end<TElem>(si);

		std::vector<DoFIndex> ind;
		for (; elem_iter != iterEnd; ++elem_iter)
		{
			u_new->inner_dof_indices(*elem_iter, fct, ind);

			// get dof positions
			std::vector<MathVector<dim> > globPos;
			InnerDoFPosition<dom_type>(globPos, *elem_iter, *u_new->domain(), lfeid);

			UG_ASSERT(globPos.size() == ind.size(),
				"#DoF mismatch: InnerDoFPosition() found " << globPos.size()
				<< ", but grid function has " << ind.size() << std::endl
				<< "on " << ElementDebugInfo(*u_new->domain()->grid(), *elem_iter) << ".");

			// write values in new grid function
			for (size_t dof = 0; dof < ind.size(); ++dof)
			{
				if (!u_orig.evaluate(DoFRef(*u_new, ind[dof]), globPos[dof]))
				{
					DoFRef(*u_new, ind[dof]) = std::numeric_limits<number>::quiet_NaN();
					//UG_THROW("Interpolation onto new grid did not succeed.\n"
					//		 "DoF with coords " << globPos[dof] << " is out of range.");
				}
			}
		}
	}
}


template <typename TGridFunction, int trueDim>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename,
	size_t step,
	number time,
	const SubsetGroup& ssg
)
{
	// for order 1, use the given grid function as is
	if (order == 1)
	{
		// print out new grid function on desired subsets
		bool printAll = ssg.size() == (size_t) u->subset_handler()->num_subsets();
		if (printAll)
			vtkOutput->print(filename, *u, step, time);
		else
		{
			for (size_t s = 0; s < ssg.size(); ++s)
				vtkOutput->print_subset(filename, *u, ssg[s], step, time);
		}

		return;
	}

	using dom_type = typename TGridFunction::domain_type;
	using algebra_type = typename TGridFunction::algebra_type;
	using position_attachment_type = typename dom_type::position_attachment_type;
	using approx_space_type = typename TGridFunction::approximation_space_type;
	using elem_type = typename TGridFunction::template dim_traits<trueDim>::grid_base_object;

	SmartPtr<approx_space_type> srcApproxSpace = u->approx_space();
	SmartPtr<dom_type> srcDom = srcApproxSpace->domain();

	// select surface elements in old grid
	MultiGrid& srcGrid = *srcDom->grid();
	Selector srcSel(srcGrid);
	srcSel.select(u->template begin<elem_type>(), u->template end<elem_type>());

	// create new domain
	dom_type* dom_ptr;
	try	{dom_ptr = new dom_type();}
	UG_CATCH_THROW("Temporary domain could not be created.");
	SmartPtr<dom_type> destDom = make_sp(dom_ptr);

	// copy grid from old domain to new domain (into a flat grid!)
	try	{CopySelected(destDom, srcDom, srcSel);}
	UG_CATCH_THROW("Temporary grid could not be created.");

	// refine
#ifdef UG_PARALLEL
	ParallelGlobalRefiner_MultiGrid refiner(*destDom->distributed_grid_manager());
#else
	MultiGrid& destGrid = *destDom->grid();
	GlobalMultiGridRefiner refiner(destGrid);
#endif
	size_t numRefs = (size_t) ceil(log2(order));
	for (size_t iref = 0; iref < numRefs; ++iref)
	{
		try	{refiner.refine();}
		UG_CATCH_THROW("Refinement step " << iref << " could not be carried out.");
	}

	// retain function group for functions being exported
	FunctionGroup fg(srcApproxSpace->dof_distribution_info(), vFct);

	// create approx space and add functions
	approx_space_type* approx_ptr;
	try {approx_ptr = new approx_space_type(destDom, algebra_type::get_type());}
	UG_CATCH_THROW("Temporary approximation space could not be created.");
	SmartPtr<approx_space_type> destApproxSpace = make_sp(approx_ptr);

	for (size_t fct = 0; fct < fg.size(); ++fct)
	{
		if (fg.function_pattern()->is_def_everywhere(fg.unique_id(fct)))
			destApproxSpace->add(fg.name(fct), "Lagrange", 1);
		else
		{
			int num_subsets = fg.function_pattern()->num_subsets();
			std::string subsets;
			for (int si = 0; si < num_subsets; ++si)
			{
				if (fg.function_pattern()->is_def_in_subset(fct, si))
					subsets += std::string(",") + fg.function_pattern()->subset_name(si);
			}
			if (!subsets.empty())
				subsets.erase(0,1);	// delete leading ","

			destApproxSpace->add(fg.name(fct), "Lagrange", 1, subsets.c_str());
		}
	}
	destApproxSpace->init_top_surface();

	TGridFunction* gridFct_ptr;
	try {gridFct_ptr = new TGridFunction(destApproxSpace);}
	UG_CATCH_THROW("Temporary grid function could not be created.");
	SmartPtr<TGridFunction> u_new = make_sp(gridFct_ptr);

	// interpolate onto new grid
	for (size_t fct = 0; fct < fg.size(); ++fct)
	{
		const LFEID lfeid = u_new->dof_distribution()->lfeid(fct);

		GlobalGridFunctionNumberData<TGridFunction, trueDim> ggfnd =
			GlobalGridFunctionNumberData<TGridFunction, trueDim>(u, fg.name(fct));

		// iterate over DoFs in new function and evaluate
		// should be vertices only for Lagrange-1
		if (u_new->max_dofs(VERTEX))
			interpolate_from_original_fct<Vertex, TGridFunction, GlobalGridFunctionNumberData<TGridFunction, trueDim> >
				(u_new, ggfnd, fct, lfeid);
		/*
		if (u_new->max_dofs(EDGE))
			interpolate_from_original_fct<Edge, TGridFunction, GlobalGridFunctionNumberData<TGridFunction, trueDim> >
				(u_new, ggfnd, fct, lfeid);
		if (u_new->max_dofs(FACE))
			interpolate_from_original_fct<Face, TGridFunction, GlobalGridFunctionNumberData<TGridFunction, trueDim> >
				(u_new, ggfnd, fct, lfeid);
		if (u_new->max_dofs(VOLUME))
			interpolate_from_original_fct<Volume, TGridFunction, GlobalGridFunctionNumberData<TGridFunction, trueDim> >
				(u_new, ggfnd, fct, lfeid);
		*/
	}

#ifdef UG_PARALLEL
	// copy storage type and layouts
	u_new->set_storage_type(u->get_storage_mask());
	u_new->set_layouts(u->layouts());
#endif

	// print out new grid function on desired subsets
	bool printAll = ssg.size() == (size_t) u->subset_handler()->num_subsets();
	if (printAll)
		vtkOutput->print(filename, *u_new, step, time);
	else
	{
		for (size_t s = 0; s < ssg.size(); ++s)
			vtkOutput->print_subset(filename, *u_new, ssg[s], step, time);
	}
}


/// export a solution from a high-order ansatz space to vtk file(s)
/**
 *  This function will create a temporary domain, copy all elements from the domain
 *  which the grid function u is defined on to the temporary domain and then refine
 *  the resulting grid until it has at least as many vertices as the original grid
 *  functions has unknowns (e.g. a grid for a function of order 2 would be refined
 *  once, a grid for a function of order 4 would be refined twice, and so on).
 *  After refinement, a temporary grid function of order 1 (Lagrange) is defined on
 *  the refined grid and its values interpolated from the original function u.
 *  The temporary first-order grid function is then exported to vtk using the usual
 *  mechanisms.
 *
 * @param u			original high-order grid function to be exported
 * @param vFct		vector of function names (contained in grid function) to be exported
 * @param order		order to be used
 * @param vtkOutput	VTKOutput object to use for export of linearized function
 * @param filename	file name to be used in export
 *
 * @todo The order parameter might be left out and determined automatically from the
 * 		 grid function.
 */
template <typename TGridFunction>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename,
	size_t step,
	number time,
	const std::vector<std::string>& vSubset
)
{
	// construct subset group from given subsets
	SubsetGroup ssg(u->subset_handler());
	try
	{
		if (vSubset.empty())
			ssg.add_all();
		else
			ssg.add(vSubset);
	}
	UG_CATCH_THROW("Subsets are faulty.");

	// find highest dim that contains any elements
	MultiGrid& srcGrid = *u->approx_space()->domain()->grid();
	if (srcGrid.num_volumes())
		vtk_export_ho<TGridFunction, TGridFunction::dim >= 3 ? 3 : TGridFunction::dim>
			(u, vFct, order, vtkOutput, filename, step, time, ssg);
	else if (srcGrid.num_faces())
		vtk_export_ho<TGridFunction, TGridFunction::dim >= 2 ? 2 : TGridFunction::dim>
			(u, vFct, order, vtkOutput, filename, step, time, ssg);
	else if (srcGrid.num_edges())
		vtk_export_ho<TGridFunction, TGridFunction::dim >= 1 ? 1 : TGridFunction::dim>
			(u, vFct, order, vtkOutput, filename, step, time, ssg);
	else if (srcGrid.num_vertices())
		vtk_export_ho<TGridFunction, TGridFunction::dim >= 0 ? 0 : TGridFunction::dim>
			(u, vFct, order, vtkOutput, filename, step, time, ssg);
	else
		vtk_export_ho<TGridFunction, TGridFunction::dim>(u, vFct, order, vtkOutput, filename, step, time, ssg);
}


template <typename TGridFunction>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename,
	const std::vector<std::string>& vSubset
	)
{
	vtk_export_ho(u, vFct, order, vtkOutput, filename, 0, 0.0, vSubset);
}


template <typename TGridFunction>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename,
	size_t step,
	number time
)
{
	vtk_export_ho(u, vFct, order, vtkOutput, filename, step, time, std::vector<std::string>());
}


template <typename TGridFunction>
void vtk_export_ho
(
	SmartPtr<TGridFunction> u,
	const std::vector<std::string>& vFct,
	size_t order,
	SmartPtr<VTKOutput<TGridFunction::domain_type::dim> > vtkOutput,
	const char* filename
	)
{
	vtk_export_ho(u, vFct, order, vtkOutput, filename, 0, 0.0, std::vector<std::string>());
}


} // end namespace ug


#endif // __H__UG__LIB_DISC__IO__VTK_EXPORT_HO__
