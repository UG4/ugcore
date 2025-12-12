/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER_IMPL__

#include "std_transfer.h"

// #include "lib_disc/reference_element/reference_mapping_provider.h"
// #include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "lib_grid/algorithms/debug_util.h"								// ElementDebugInfo

namespace ug {


template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
assemble_prolongation_p1(matrix_type& P,
                         const DoFDistribution& fineDD,
                         const DoFDistribution& coarseDD)
{
	PROFILE_FUNC_GROUP("gmg");
// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.lfeid(fct).type() != LFEID::LAGRANGE ||
			fineDD.lfeid(fct).order() != 1)
			UG_THROW("AssembleStdProlongationForP1Lagrange:"
				"Interpolation only implemented for Lagrange P1 functions.");

//  resize matrix
	P.resize_and_clear(fineDD.num_indices(), coarseDD.num_indices());

//  iterators
	const MultiGrid& mg = *coarseDD.multi_grid();
	using const_iterator = DoFDistribution::traits<Vertex>::const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	std::vector<size_t> vParentIndex, vChildIndex;
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.begin<Vertex>(si);
		iterEnd = fineDD.end<Vertex>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Vertex* child = *iter;

		//  get father
			GridObject* parent = mg.get_parent(child);

		//	check if child contained in coarseDD. This should always be false
		//	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
		//	and an adaptive grid-part used by both dds. In such a case we can
		//	simply set identity.
			if(coarseDD.is_contained(child)){
			//	get indices
				coarseDD.inner_algebra_indices(child, vParentIndex);
				fineDD.inner_algebra_indices(child, vChildIndex);
				UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");

			//	set identity
				for(size_t i = 0; i < vParentIndex.size(); ++i)
					P(vChildIndex[i], vParentIndex[i]) = 1.0;

			//	this child is perfectly handled
				continue;
			}
			else{
			//	check if parent exists (this should always be the case, except in
			//	the case that 'child' is a v-slave)
				if(!parent) continue;

				if(!coarseDD.is_contained(parent)){
					UG_THROW("StdTransfer: Parent element \n"
							<< ElementDebugInfo(mg, parent) <<
							"is not contained in coarse-dd nor the child element\n"
							<< ElementDebugInfo(mg, child) <<
							" in the coarse-dd. This should not happen.")
				}
			}

		//	type of father
			const ReferenceObjectID roid = parent->reference_object_id();

		//	loop all components
			for(size_t fct = 0; fct < fineDD.num_fct(); fct++)
			{
			//	check that fct defined on subset
				if(!fineDD.is_def_in_subset(fct, si)) continue;

			//  get global indices
				fineDD.inner_dof_indices(child, fct, vChildDoF);

			//	detect type of father
				switch(roid)
				{
					case ROID_VERTEX:
					{
						auto vrt = dynamic_cast<Vertex*>(parent);
						coarseDD.inner_dof_indices(vrt, fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 1.0;
					}
					break;
					case ROID_EDGE:
					for(int i = 0; i < 2; ++i)
					{
						auto edge = dynamic_cast<Edge*>(parent);
						coarseDD.inner_dof_indices(edge->vertex(i), fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 0.5;
					}
					break;
					case ROID_QUADRILATERAL:
					for(int i = 0; i < 4; ++i)
					{
						auto face = dynamic_cast<Face*>(parent);
						coarseDD.inner_dof_indices(face->vertex(i), fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 0.25;
					}
					break;
					case ROID_HEXAHEDRON:
					for(int i = 0; i < 8; ++i)
					{
						auto hexaeder = dynamic_cast<Volume*>(parent);
						coarseDD.inner_dof_indices(hexaeder->vertex(i), fct, vParentDoF);
						DoFRef(P, vChildDoF[0], vParentDoF[0]) = 0.125;
					}
					break;
					default: UG_THROW("AssembleStdProlongationForP1Lagrange: Element father"
									 " is of unsupported type "<< roid << " for "
									 << ElementDebugInfo(mg, child) << ".");
				}
			}
		}
	}
}
/*
template <typename TDomain>
void ProjectGlobalPositionToElem(std::vector<MathVector<TDomain::dim> >& vGlobPos,
                                 GridObject* parent, const TDomain& domain)
{
	const int parentDim = parent->base_object_id();

	// vertex and full dim parent must match
	if(parentDim == 0 || parentDim == TDomain::dim)
		return;

//	get the vertices
	std::vector<MathVector<TDomain::dim> > vCornerCoord;
	switch(parentDim)
	{
		case EDGE:
		{
			CollectCornerCoordinates(vCornerCoord, *static_cast<Edge*>(parent), domain, true);
			MathVector<TDomain::dim> dir;
			VecSubtract(dir, vCornerCoord[1], vCornerCoord[0]);
			for(size_t p = 0; p < vGlobPos.size(); ++p){
				ProjectPointToRay(vGlobPos[p], vGlobPos[p], vCornerCoord[0], dir);
			}
		}
		break;
		case FACE:
		{
			CollectCornerCoordinates(vCornerCoord, *static_cast<Face*>(parent), domain, true);
			MathVector<TDomain::dim> normal;
			MathVector<TDomain::dim> a, b;
			VecSubtract(a, vCornerCoord[1], vCornerCoord[0]);
			VecSubtract(b, vCornerCoord[2], vCornerCoord[0]);
			VecCross(normal, a,b);

			for(size_t p = 0; p < vGlobPos.size(); ++p){
				ProjectPointToPlane(vGlobPos[p], vGlobPos[p], vCornerCoord[0], normal);
			}
		}
		break;
		default: UG_THROW( "Base Object type not found.");
	}
}
*/


template <typename TDomain, typename TAlgebra>
template <typename TChild>
void StdTransfer<TDomain, TAlgebra>::
assemble_prolongation(matrix_type& P,
                      const DoFDistribution& fineDD,
                      const DoFDistribution& coarseDD,
                      ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");

//  iterators
	MultiGrid& mg = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());
	using const_iterator = typename DoFDistribution::traits<TChild>::const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	std::vector<size_t> vParentIndex, vChildIndex;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.begin<TChild>(si);
		iterEnd = fineDD.end<TChild>(si);

	//	check, which cmps to consider on this subset
		std::vector<LFEID> vLFEID;
		std::vector<size_t> vFct;
		for(size_t fct = 0; fct < fineDD.num_fct(); ++fct){
			if(fineDD.max_fct_dofs(fct, TChild::dim, si) == 0) continue;
			vFct.push_back(fct);
			vLFEID.push_back(fineDD.lfeid(fct));
		}
		if(vFct.empty()) continue;

	//  loop elems on coarse level for subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get child
			TChild* child = *iter;

		//	get parent
			GridObject* parent = mg.get_parent(child);

		//	check if child contained in coarseDD. This should always be false
		//	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
		//	and an adaptive grid-part used by both dds. In such a case we can
		//	simply set identity.
			if(coarseDD.is_contained(child)){
			//	get indices
				coarseDD.inner_algebra_indices(child, vParentIndex);
				fineDD.inner_algebra_indices(child, vChildIndex);
				UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");

			//	set identity
				for(size_t i = 0; i < vParentIndex.size(); ++i)
					P(vChildIndex[i], vParentIndex[i]) = 1.0;

			//	this child is perfectly handled
				continue;
			}
			else{

			//	check if parent exists (this should always be the case, except in
			//	the case that 'child' is a v-slave)
				if(!parent) continue;

				if(!coarseDD.is_contained(parent)){
					UG_THROW("StdTransfer: A parent element is not contained in "
							" coarse-dd nor the child element in the coarse-dd. "
							"This should not happen.")
				}
			}

		//	loop all components
			for(size_t f = 0; f < vFct.size(); f++)
			{
			//	get comp and lfeid
				const size_t fct = vFct[f];
				const LFEID& lfeID = vLFEID[f];

			//  get global indices
				fineDD.inner_dof_indices(child, fct, vChildDoF);

			//	switch space type
				switch(lfeID.type())
				{
					case LFEID::PIECEWISE_CONSTANT:
					{
						coarseDD.dof_indices(parent, fct, vParentDoF);
						UG_ASSERT(vChildDoF.size() == 1, "Must be one.");
						UG_ASSERT(vParentDoF.size() == 1, "Must be one.");

						DoFRef(P, vChildDoF[0], vParentDoF[0]) =  1.0;
					}
					break;

					case LFEID::CROUZEIX_RAVIART:
					{
					//	get dimension of parent
						const int parentDim = parent->base_object_id();
						std::vector<GridObject*> vParent;

					//	check if to interpolate from neighbor elems
						if(parentDim == lfeID.dim()){
							// case: Side inner to parent. --> Parent fine.
							vParent.push_back(parent);
						} else if(parentDim == lfeID.dim() - 1){
							// case: parent is Side. --> Get neighbor elems
							using TElem = typename TChild::sideof;
							std::vector<TElem*> vElem;
							coarseDD.collect_associated(vElem, parent);
							for(size_t p = 0; p < vElem.size(); ++p)
								vParent.push_back(vElem[p]);

						} else {
							UG_THROW("StdTransfer: For CR parent must be full-dim "
									"elem or a side (dim-1). But has dim: "<<parentDim);
						}


					//	global positions of fine dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	loop contributions from parents
						for(size_t i = 0; i < vParent.size(); ++i)
						{
						//	get coarse indices
							coarseDD.dof_indices(vParent[i], fct, vParentDoF);

						//	get shapes at global positions
							std::vector<std::vector<number> > vvShape;
							ShapesAtGlobalPosition(vvShape, vDoFPos, vParent[i], *spDomain, lfeID);

						//	add restriction
							for(size_t ip = 0; ip < vvShape.size(); ++ip)
								for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
									DoFRef(P, vChildDoF[ip], vParentDoF[sh]) +=
											(1./vParent.size()) * vvShape[ip][sh];
						}
					}
					break;

					case LFEID::LAGRANGE:
					{
					//	get coarse indices
						coarseDD.dof_indices(parent, fct, vParentDoF);

					//	global positions of child dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	project
					//	ProjectGlobalPositionToElem(vDoFPos, parent, *spDomain);

					//	get shapes at global positions
						std::vector<std::vector<number> > vvShape;
						ShapesAtGlobalPosition(vvShape, vDoFPos, parent, *spDomain, lfeID);

					//	set restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								DoFRef(P, vChildDoF[ip], vParentDoF[sh]) = vvShape[ip][sh];
					}
					break;

					default:
						UG_THROW("StdTransfer: Local-Finite-Element: "<<lfeID<<
						         " is not supported by this Transfer.")

				} // end LFEID-switch
			} // end fct - cmps
		} // end fine - elements
	} // end subset
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
assemble_prolongation(matrix_type& P,
                      const DoFDistribution& fineDD,
                      const DoFDistribution& coarseDD,
                      ConstSmartPtr<TDomain> spDomain)
{
	//  resize matrix
	P.resize_and_clear(fineDD.num_indices(), coarseDD.num_indices());

	// loop all base types carrying indices on fine elems
	if(fineDD.max_dofs(VERTEX)) assemble_prolongation<Vertex>(P, fineDD, coarseDD, spDomain);
	if(fineDD.max_dofs(EDGE)) assemble_prolongation<Edge>(P, fineDD, coarseDD, spDomain);
	if(fineDD.max_dofs(FACE)) assemble_prolongation<Face>(P, fineDD, coarseDD, spDomain);
	if(fineDD.max_dofs(VOLUME)) assemble_prolongation<Volume>(P, fineDD, coarseDD, spDomain);
}


template <typename TDomain, typename TAlgebra>
template <typename TChild>
void StdTransfer<TDomain, TAlgebra>::
assemble_restriction(matrix_type& R,
                     const DoFDistribution& coarseDD,
                     const DoFDistribution& fineDD,
                     ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");

//  iterators
	MultiGrid& mg = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());
	using const_iterator = typename DoFDistribution::traits<TChild>::const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	std::vector<size_t> vParentIndex, vChildIndex;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.begin<TChild>(si);
		iterEnd = fineDD.end<TChild>(si);

	//	check, which cmps to consider on this subset
		std::vector<LFEID> vLFEID;
		std::vector<size_t> vFct;
		for(size_t fct = 0; fct < fineDD.num_fct(); ++fct){
			if(fineDD.max_fct_dofs(fct, TChild::dim, si) == 0) continue;
			vFct.push_back(fct);
			vLFEID.push_back(fineDD.lfeid(fct));
		}
		if(vFct.empty()) continue;

	//  loop elems on coarse level for subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get child
			TChild* child = *iter;

		//	get parent
			GridObject* parent = mg.get_parent(child);

		//	check if child contained in coarseDD. This should always be false
		//	for a GridLevel::LEVEL, but might be the case for GridLevel::SURFACE
		//	and an adaptive grid-part used by both dds. In such a case we can
		//	simply set identity.
			if(coarseDD.is_contained(child)){
			//	get indices
				coarseDD.inner_algebra_indices(child, vParentIndex);
				fineDD.inner_algebra_indices(child, vChildIndex);
				UG_ASSERT(vParentIndex.size() == vChildIndex.size(), "Size mismatch");

			//	set identity
				for(size_t i = 0; i < vParentIndex.size(); ++i)
					R(vParentIndex[i], vChildIndex[i]) = 1.0;

			//	this child is perfectly handled
				continue;
			}
			else{

			//	check if parent exists (this should always be the case, except in
			//	the case that 'child' is a v-slave)
				if(!parent) continue;

				if(!coarseDD.is_contained(parent)){
					UG_THROW("StdTransfer: A parent element is not contained in "
							" coarse-dd nor the child element in the coarse-dd. "
							"This should not happen.")
				}
			}

		//	loop all components
			for(size_t f = 0; f < vFct.size(); f++)
			{
			//	get comp and lfeid
				const size_t fct = vFct[f];
				const LFEID& lfeID = vLFEID[f];

			//  get global indices
				fineDD.inner_dof_indices(child, fct, vChildDoF);

			//	switch space type
				switch(lfeID.type())
				{
					case LFEID::PIECEWISE_CONSTANT:
					{
						coarseDD.dof_indices(parent, fct, vParentDoF);
						UG_ASSERT(vChildDoF.size() == 1, "Must be one.");
						UG_ASSERT(vParentDoF.size() == 1, "Must be one.");

						DoFRef(R, vParentDoF[0], vChildDoF[0]) =  1.0;
					}
					break;

					case LFEID::CROUZEIX_RAVIART:
					{
					//	get dimension of parent
						const int parentDim = parent->base_object_id();
						std::vector<GridObject*> vParent;

					//	check if to interpolate from neighbor elems
						if(parentDim == lfeID.dim()){
							// case: Side inner to parent. --> Parent fine.
							vParent.push_back(parent);
						} else if(parentDim == lfeID.dim() - 1){
							// case: parent is Side. --> Get neighbor elems
							using TElem = typename TChild::sideof;
							std::vector<TElem*> vElem;
							coarseDD.collect_associated(vElem, parent);
							for(size_t p = 0; p < vElem.size(); ++p){
							//	NOTE: This is not the transposed of the prolongation
							//		  in adaptive case, since we only restrict to
							//		  covered parts.
								if(mg.num_children<TElem>(vElem[p]) > 0)
									vParent.push_back(vElem[p]);
							}

						} else {
							UG_THROW("StdTransfer: For CR parent must be full-dim "
									"elem or a side (dim-1). But has dim: "<<parentDim);
						}


					//	global positions of fine dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	loop contributions from parents
						for(size_t i = 0; i < vParent.size(); ++i)
						{
						//	get coarse indices
							coarseDD.dof_indices(vParent[i], fct, vParentDoF);

						//	get shapes at global positions
							std::vector<std::vector<number> > vvShape;
							ShapesAtGlobalPosition(vvShape, vDoFPos, vParent[i], *spDomain, lfeID);

						//	add restriction
							for(size_t ip = 0; ip < vvShape.size(); ++ip)
								for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
									DoFRef(R, vParentDoF[sh], vChildDoF[ip]) +=
											(1./vParent.size()) * vvShape[ip][sh];
						}
					}
					break;

					case LFEID::LAGRANGE:
					{
					//	get coarse indices
						coarseDD.dof_indices(parent, fct, vParentDoF);

					//	global positions of child dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						InnerDoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	get shapes at global positions
						std::vector<std::vector<number> > vvShape;
						ShapesAtGlobalPosition(vvShape, vDoFPos, parent, *spDomain, lfeID);

					//	set restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								DoFRef(R, vParentDoF[sh], vChildDoF[ip]) = vvShape[ip][sh];
					}
					break;

					default:
						UG_THROW("StdTransfer: Local-Finite-Element: "<<lfeID<<
						         " is not supported by this Transfer.")

				} // end LFEID-switch
			} // end fct - cmps
		} // end fine - elements
	} // end subset
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
assemble_restriction(matrix_type& R,
                     const DoFDistribution& coarseDD,
                     const DoFDistribution& fineDD,
                     ConstSmartPtr<TDomain> spDomain)
{
	//  resize matrix
	R.resize_and_clear(coarseDD.num_indices(), fineDD.num_indices());

	// loop all base types carrying indices on fine elems
	if(fineDD.max_dofs(VERTEX)) assemble_restriction<Vertex>(R, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(EDGE)) assemble_restriction<Edge>(R, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(FACE)) assemble_restriction<Face>(R, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(VOLUME)) assemble_restriction<Volume>(R, coarseDD, fineDD, spDomain);
}


template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
StdTransfer<TDomain, TAlgebra>::
prolongation(const GridLevel& fineGL, const GridLevel& coarseGL,
             ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
	if(fineGL.level() - coarseGL.level() != 1)
		UG_THROW("StdTransfer: Can only project between successive level, "
				"but fine = "<<fineGL<<", coarse = "<<coarseGL);

	if(fineGL.type() != coarseGL.type())
		UG_THROW("StdTransfer: Can only project between dof distributions of "
				"same type, but fine = "<<fineGL<<", coarse = "<<coarseGL);

	// remove old revisions
	remove_outdated(m_mProlongation, spApproxSpace->revision());

	// key of this restriction
	TransferKey key(coarseGL, fineGL, spApproxSpace->revision());

	// check if must be created
	if(m_mProlongation.find(key) == m_mProlongation.end())
	{
		SmartPtr<matrix_type> P =
				m_mProlongation[key] = SmartPtr<matrix_type>(new matrix_type);

		ConstSmartPtr<DoFDistribution> spCoarseDD = spApproxSpace->dof_distribution(coarseGL);
		ConstSmartPtr<DoFDistribution> spFineDD = spApproxSpace->dof_distribution(fineGL);

		bool P1LagrangeOnly = false;
		if(m_p1LagrangeOptimizationEnabled){
			P1LagrangeOnly = true;
			for(size_t fct = 0; fct < spApproxSpace->num_fct(); ++fct)
				if(spApproxSpace->lfeid(fct).type() != LFEID::LAGRANGE ||
					spApproxSpace->lfeid(fct).order() != 1)
					P1LagrangeOnly = false;
		}

		if(P1LagrangeOnly){
			assemble_prolongation_p1(*P, *spFineDD, *spCoarseDD);
		} else{
			assemble_prolongation(*P, *spFineDD, *spCoarseDD, spApproxSpace->domain());
		}

		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_prolongation(*P, spFineDD, spCoarseDD, type);
			}
		}

		#ifdef UG_PARALLEL
		P->set_storage_type(PST_CONSISTENT);
		#endif

		write_debug(*P, "P", fineGL, coarseGL);
	}

	return m_mProlongation[key];
}

template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
StdTransfer<TDomain, TAlgebra>::
restriction(const GridLevel& coarseGL, const GridLevel& fineGL,
            ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
	if(fineGL.level() - coarseGL.level() != 1)
		UG_THROW("StdTransfer: Can only project between successive level, "
				"but fine = "<<fineGL<<", coarse = "<<coarseGL);

	if(fineGL.type() != coarseGL.type())
		UG_THROW("StdTransfer: Can only project between dof distributions of "
				"same type, but fine = "<<fineGL<<", coarse = "<<coarseGL);

	// remove old revisions
	remove_outdated(m_mRestriction, spApproxSpace->revision());

	// key of this restriction
	TransferKey key(coarseGL, fineGL, spApproxSpace->revision());

	// check if must be created
	if(m_mRestriction.find(key) == m_mRestriction.end())
	{
		SmartPtr<matrix_type> R =
				m_mRestriction[key] = SmartPtr<matrix_type>(new matrix_type);

		ConstSmartPtr<DoFDistribution> spCoarseDD = spApproxSpace->dof_distribution(coarseGL);
		ConstSmartPtr<DoFDistribution> spFineDD = spApproxSpace->dof_distribution(fineGL);

		if(m_bUseTransposed)
			R->set_as_transpose_of(*prolongation(fineGL, coarseGL, spApproxSpace));
		else
			assemble_restriction(*R, *spCoarseDD, *spFineDD, spApproxSpace->domain());


		#ifdef UG_PARALLEL
		R->set_storage_type(PST_CONSISTENT);
		#endif

		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_restriction(*R, spCoarseDD, spFineDD, type);
			}
		}

		write_debug(*R, "R", coarseGL, fineGL);
	}

	return m_mRestriction[key];
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
prolongate(GF& uFine, const GF& uCoarse)
{
	PROFILE_FUNC_GROUP("gmg");

	if(!bCached)
		UG_THROW("StdTransfer: currently only cached implemented.");

	const GridLevel& coarseGL = uCoarse.grid_level();
	const GridLevel& fineGL = uFine.grid_level();
	ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace = uFine.approx_space();
	if(uCoarse.approx_space() != spApproxSpace)
		UG_THROW("StdTransfer: cannot prolongate between grid functions from "
				"different approximation spaces.");

	try{
		//prolongation(fineGL, coarseGL, spApproxSpace)->apply(uFine, uCoarse);
#ifdef UG_PARALLEL
		MatMultDirect(uFine, m_dampProl, *prolongation(fineGL, coarseGL, spApproxSpace), uCoarse);
#else
		prolongation(fineGL, coarseGL, spApproxSpace)->axpy(uFine, 0.0, uFine, m_dampProl, uCoarse);
#endif

	// 	adjust using constraints
		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_prolongation(uFine, fineGL, uCoarse, coarseGL, type);
			}
		}

	}
	UG_CATCH_THROW("StdTransfer:prolongation: Failed for fine = "<<fineGL<<" and "
	               " coarse = "<<coarseGL);

// 	check CR functions
#ifdef UG_PARALLEL
	bool bCROnly = true;
	for(size_t fct = 0; fct < spApproxSpace->num_fct(); ++fct)
		if(spApproxSpace->lfeid(fct).type() != LFEID::CROUZEIX_RAVIART &&
				spApproxSpace->lfeid(fct).type() != LFEID::PIECEWISE_CONSTANT)
			bCROnly = false;

	if(bCROnly){
		ScaleLayoutValues(&uFine, uFine.layouts()->master(), 0.5);
		ScaleLayoutValues(&uFine, uFine.layouts()->slave(), 0.5);
		AdditiveToConsistent(&uFine, uFine.layouts()->master(), uFine.layouts()->slave(),
		                     &uFine.layouts()->comm());
	}
#endif
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
do_restrict(GF& uCoarse, const GF& uFine)
{
	PROFILE_FUNC_GROUP("gmg");

	if(!bCached)
		UG_THROW("StdTransfer: currently only cached implemented.");

	const GridLevel& coarseGL = uCoarse.grid_level();
	const GridLevel& fineGL = uFine.grid_level();
	ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace = uFine.approx_space();
	if(uCoarse.approx_space() != spApproxSpace)
		UG_THROW("StdTransfer: cannot prolongate between grid functions from "
				"different approximation spaces.");
	try{

		restriction(coarseGL, fineGL, spApproxSpace)->
				apply_ignore_zero_rows(uCoarse, m_dampRes, uFine);

	// 	adjust using constraints
		for (int type = 1; type < CT_ALL; type = type << 1)
		{
			for (size_t i = 0; i < m_vConstraint.size(); ++i)
			{
				if (m_vConstraint[i]->type() & type)
					m_vConstraint[i]->adjust_restriction(uCoarse, coarseGL, uFine, fineGL, type);
			}
		}

	} UG_CATCH_THROW("StdTransfer:do_restrict: Failed for fine = "<<fineGL<<" and "
	                 " coarse = "<<coarseGL);
}

template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferOperator<TDomain, TAlgebra> >
StdTransfer<TDomain, TAlgebra>::clone()
{
	SmartPtr<StdTransfer> op(new StdTransfer);
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		op->add_constraint(m_vConstraint[i]);
	op->set_restriction_damping(m_dampRes);
	op->set_prolongation_damping(m_dampProl);
	op->set_debug(m_spDebugWriter);
	op->enable_p1_lagrange_optimization(p1_lagrange_optimization_enabled());
	op->set_use_transposed(m_bUseTransposed);
	return op;
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
write_debug(const matrix_type& mat, std::string name,
            const GridLevel& glTo, const GridLevel& glFrom)
{
	PROFILE_FUNC_GROUP("debug");
//	if no debug writer set, we're done
	if(m_spDebugWriter.invalid()) return;

//	cast dbg writer
	SmartPtr<GridFunctionDebugWriter<TDomain, TAlgebra> > dbgWriter =
			m_spDebugWriter.template cast_dynamic<GridFunctionDebugWriter<TDomain, TAlgebra> >();

//	check success
	if(dbgWriter.invalid()) return;

//	add iter count to name
	name.append("_").append(ToString(glTo.level()));
	name.append("_").append(ToString(glFrom.level()));
	name.append(".mat");

//	write
	GridLevel gridLev = dbgWriter->grid_level();
	dbgWriter->set_grid_levels(glFrom, glTo);
	dbgWriter->write_matrix(mat, name.c_str());
	dbgWriter->set_grid_level(gridLev);
}

} // end namespace ug

#endif