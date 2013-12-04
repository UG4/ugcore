/*
 * std_transfer_impl.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER_IMPL__

#include "std_transfer.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/function_spaces/grid_function_util.h"

namespace ug{


template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
assemble_restriction_p1(matrix_type& mat,
                         const DoFDistribution& coarseDD, const DoFDistribution& fineDD)
{
	PROFILE_FUNC_GROUP("gmg");
// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.lfeid(fct).type() != LFEID::LAGRANGE ||
			fineDD.lfeid(fct).order() != 1)
			UG_THROW("AssembleStdProlongationForP1Lagrange:"
				"Interpolation only implemented for Lagrange P1 functions.");

//  resize matrix
	mat.resize_and_clear(coarseDD.num_indices(), fineDD.num_indices());

//  iterators
	const MultiGrid& mg = *coarseDD.multi_grid();
	typedef DoFDistribution::traits<VertexBase>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	std::vector<size_t> vParentIndex, vChildIndex;
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<VertexBase>(si);
		iterEnd = fineDD.template end<VertexBase>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			VertexBase* child = *iter;

		//  get father
			GeometricObject* parent = mg.get_parent(child);

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
					mat(vParentIndex[i], vChildIndex[i]) = 1.0;

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
						VertexBase* vrt = dynamic_cast<VertexBase*>(parent);
						coarseDD.inner_dof_indices(vrt, fct, vParentDoF);
						DoFRef(mat, vParentDoF[0], vChildDoF[0]) = 1.0;
					}
					break;
					case ROID_EDGE:
					for(int i = 0; i < 2; ++i)
					{
						EdgeBase* edge = dynamic_cast<EdgeBase*>(parent);
						coarseDD.inner_dof_indices(edge->vertex(i), fct, vParentDoF);
						DoFRef(mat, vParentDoF[0], vChildDoF[0]) = 0.5;
					}
					break;
					case ROID_QUADRILATERAL:
					for(int i = 0; i < 4; ++i)
					{
						Face* face = dynamic_cast<Face*>(parent);
						coarseDD.inner_dof_indices(face->vertex(i), fct, vParentDoF);
						DoFRef(mat, vParentDoF[0], vChildDoF[0]) = 0.25;
					}
					break;
					case ROID_HEXAHEDRON:
					for(int i = 0; i < 8; ++i)
					{
						Volume* hexaeder = dynamic_cast<Volume*>(parent);
						coarseDD.inner_dof_indices(hexaeder->vertex(i), fct, vParentDoF);
						DoFRef(mat, vParentDoF[0], vChildDoF[0]) = 0.125;
					}
					break;
					default: UG_THROW("AssembleStdProlongationForP1Lagrange: Element Father"
									 "is of unsupported type "<<roid);
				}
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
template <typename TChild>
void StdTransfer<TDomain, TAlgebra>::
assemble_restriction_elemwise(matrix_type& mat,
                              const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
                              ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");

//  resize matrix
	mat.resize_and_clear(coarseDD.num_indices(), fineDD.num_indices());

//  iterators
	MultiGrid& mg = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());
	typedef typename DoFDistribution::traits<TChild>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	std::vector<DoFIndex> vParentDoF, vChildDoF;
	std::vector<size_t> vParentIndex, vChildIndex;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<TChild>(si);
		iterEnd = fineDD.template end<TChild>(si);

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
			GeometricObject* parent = mg.get_parent(child);

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
					mat(vParentIndex[i], vChildIndex[i]) = 1.0;

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

						DoFRef(mat, vParentDoF[0], vChildDoF[0]) =  1.0;
					}
					break;

					case LFEID::CROUZEIX_RAVIART:
					{
					//	get dimension of parent
						const int parentDim = parent->base_object_id();
						std::vector<GeometricObject*> vParent;

					//	check if to interpolate from neighbor elems
						if(parentDim == lfeID.dim()){
							// case: Side inner to parent. --> Parent fine.
							vParent.push_back(parent);
						} else if(parentDim == lfeID.dim() - 1){
							// case: parent is Side. --> Get neighbor elems
							typedef typename TChild::sideof TElem;
							typename Grid::traits<TElem>::secure_container vElem;
							mg.associated_elements(vElem, parent);
							for(size_t p = 0; p < vElem.size(); ++p)
								vParent.push_back(vElem[p]);

						} else {
							UG_THROW("StdTransfer: For CR parent must be full-dim "
									"elem or a side (dim-1). But has dim: "<<parentDim);
						}


					//	global positions of fine dofs
						std::vector<MathVector<TDomain::dim> > vDoFPos;
						DoFPosition(vDoFPos, child, *spDomain, lfeID);

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
									DoFRef(mat, vParentDoF[sh], vChildDoF[ip]) +=
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
						DoFPosition(vDoFPos, child, *spDomain, lfeID);

					//	get shapes at global positions
						std::vector<std::vector<number> > vvShape;
						ShapesAtGlobalPosition(vvShape, vDoFPos, parent, *spDomain, lfeID);

					//	set restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								DoFRef(mat, vParentDoF[sh], vChildDoF[ip]) = vvShape[ip][sh];
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
assemble_restriction_elemwise(matrix_type& mat,
                              const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
                              ConstSmartPtr<TDomain> spDomain)
{
	// loop all base types carrying indices on fine elems
	if(fineDD.max_dofs(VERTEX)) assemble_restriction_elemwise<VertexBase>(mat, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(EDGE)) assemble_restriction_elemwise<EdgeBase>(mat, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(FACE)) assemble_restriction_elemwise<Face>(mat, coarseDD, fineDD, spDomain);
	if(fineDD.max_dofs(VOLUME)) assemble_restriction_elemwise<Volume>(mat, coarseDD, fineDD, spDomain);
}


template <typename TDomain, typename TAlgebra>
template <typename TElem>
void StdTransfer<TDomain, TAlgebra>::
set_identity_on_pure_surface(matrix_type& mat,
                             const DoFDistribution& coarseDD, const DoFDistribution& fineDD)
{
	PROFILE_FUNC_GROUP("gmg");

	std::vector<size_t> vCoarseIndex, vFineIndex;
	const MultiGrid& mg = *coarseDD.multi_grid();

//  iterators
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	for(int si = 0; si < coarseDD.num_subsets(); ++si)
	{
		iterBegin = coarseDD.template begin<TElem>(si);
		iterEnd = coarseDD.template end<TElem>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			TElem* coarseElem = *iter;

			const size_t numChild = mg.num_children<TElem, TElem>(coarseElem);
			if(numChild != 0) continue;

		//	get indices
			coarseDD.inner_algebra_indices(coarseElem, vCoarseIndex);
			fineDD.inner_algebra_indices(coarseElem, vFineIndex);
			UG_ASSERT(vCoarseIndex.size() == vFineIndex.size(), "Size mismatch");

		//	set identity
			for(size_t i = 0; i < vCoarseIndex.size(); ++i)
				mat(vCoarseIndex[i], vFineIndex[i]) = 1.0;
		}
	}
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
set_identity_on_pure_surface(matrix_type& mat,
                             const DoFDistribution& coarseDD, const DoFDistribution& fineDD)
{
	if(coarseDD.max_dofs(VERTEX)) set_identity_on_pure_surface<VertexBase>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(EDGE)) set_identity_on_pure_surface<EdgeBase>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(FACE)) set_identity_on_pure_surface<Face>(mat, coarseDD, fineDD);
	if(coarseDD.max_dofs(VOLUME)) set_identity_on_pure_surface<Volume>(mat, coarseDD, fineDD);
}

template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
StdTransfer<TDomain, TAlgebra>::
prolongation(const GridLevel& fineGL, const GridLevel& coarseGL,
             ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
	// remove old revisions
	remove_outdated(m_mProlongation, spApproxSpace->revision());

	// key of this prolongation
	TransferKey key(coarseGL, fineGL, spApproxSpace->revision());

	// check if must be created
	if(m_mProlongation.find(key) == m_mProlongation.end()){

		SmartPtr<matrix_type> P =
				m_mProlongation[key] = SmartPtr<matrix_type>(new matrix_type);

		P->set_as_transpose_of(*restriction(coarseGL, fineGL, spApproxSpace));
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

		bool P1LagrangeOnly = false;
		if(m_p1LagrangeOptimizationEnabled){
			P1LagrangeOnly = true;
			for(size_t fct = 0; fct < spApproxSpace->num_fct(); ++fct)
				if(spApproxSpace->lfeid(fct).type() != LFEID::LAGRANGE ||
					spApproxSpace->lfeid(fct).order() != 1)
					P1LagrangeOnly = false;
		}

		if(P1LagrangeOnly){
			assemble_restriction_p1(*R, *spCoarseDD, *spFineDD);
		} else{
			assemble_restriction_elemwise(*R, *spCoarseDD, *spFineDD, spApproxSpace->domain());
		}

		if(coarseGL.is_surface()){
//			set_identity_on_pure_surface(*R, *spCoarseDD, *spFineDD);
		}

		#ifdef UG_PARALLEL
		R->set_storage_type(PST_CONSISTENT);
		#endif

		for(size_t i = 0; i < m_vConstraint.size(); ++i){
			if (m_vConstraint[i]->type() & CT_DIRICHLET){
				m_vConstraint[i]->adjust_restriction(*R, spCoarseDD, spFineDD);
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

		// check if must be created
		const RevisionCounter& revCnt = spApproxSpace->revision();
		if(m_mProlongation.find(TransferKey(fineGL, coarseGL, revCnt)) != m_mProlongation.end()){
				prolongation(fineGL, coarseGL, spApproxSpace)->apply(uFine, uCoarse);
		} else {
			restriction(coarseGL, fineGL, spApproxSpace)->apply_transposed(uFine, uCoarse);
		}

		// call prolongations due to added constraints (= adjust_restrict, member of class constraint)
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			m_vConstraint[i]->adjust_prolongation(uFine, fineGL, uCoarse, coarseGL);

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

		restriction(coarseGL, fineGL, spApproxSpace)->apply_ignore_zero_rows(uCoarse, m_dampRes, uFine);

		// call restrictions due to added constraints (= adjust_restrict, member of class constraint)
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			m_vConstraint[i]->adjust_restriction(uCoarse, coarseGL, uFine, fineGL);

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
	op->set_debug(m_spDebugWriter);
	op->enable_p1_lagrange_optimization(p1_lagrange_optimization_enabled());
	return op;
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
add_constraint(SmartPtr<IConstraint<TAlgebra> > pp)
{
//	add only once
	if(std::find(m_vConstraint.begin(), m_vConstraint.end(), pp) !=
			m_vConstraint.end()) return;
	m_vConstraint.push_back(pp);
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp)
{
	m_vConstraint.erase(m_vConstraint.begin(),
	                     std::remove(m_vConstraint.begin(), m_vConstraint.end(), pp));
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

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_TRANSFER_IMPL__ */
