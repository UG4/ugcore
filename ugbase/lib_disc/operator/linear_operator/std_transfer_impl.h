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
	std::vector<DoFIndex> vCoarseDoF, vFineDoF;
	for(int si = 0; si < fineDD.num_subsets(); ++si)
	{
		iterBegin = fineDD.template begin<VertexBase>(si);
		iterEnd = fineDD.template end<VertexBase>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			VertexBase* fineVrt = *iter;

		//  get father
			GeometricObject* parent = mg.get_parent(fineVrt);

			if(!parent)
				continue;

		//	type of father
			const ReferenceObjectID roid = parent->reference_object_id();

		//	loop all components
			for(size_t fct = 0; fct < fineDD.num_fct(); fct++)
			{
			//	check that fct defined on subset
				if(!fineDD.is_def_in_subset(fct, si)) continue;

			//  get global indices
				fineDD.inner_dof_indices(fineVrt, fct, vFineDoF);

			//	detect type of father
				switch(roid)
				{
					case ROID_VERTEX:
					{
						VertexBase* vrt = dynamic_cast<VertexBase*>(parent);
						coarseDD.inner_dof_indices(vrt, fct, vCoarseDoF);
						DoFRef(mat, vCoarseDoF[0], vFineDoF[0]) = 1.0;
					}
					break;
					case ROID_EDGE:
					for(int i = 0; i < 2; ++i)
					{
						EdgeBase* edge = dynamic_cast<EdgeBase*>(parent);
						coarseDD.inner_dof_indices(edge->vertex(i), fct, vCoarseDoF);
						DoFRef(mat, vCoarseDoF[0], vFineDoF[0]) = 0.5;
					}
					break;
					case ROID_QUADRILATERAL:
					for(int i = 0; i < 4; ++i)
					{
						Face* face = dynamic_cast<Face*>(parent);
						coarseDD.inner_dof_indices(face->vertex(i), fct, vCoarseDoF);
						DoFRef(mat, vCoarseDoF[0], vFineDoF[0]) = 0.25;
					}
					break;
					case ROID_HEXAHEDRON:
					for(int i = 0; i < 8; ++i)
					{
						Volume* hexaeder = dynamic_cast<Volume*>(parent);
						coarseDD.inner_dof_indices(hexaeder->vertex(i), fct, vCoarseDoF);
						DoFRef(mat, vCoarseDoF[0], vFineDoF[0]) = 0.125;
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
void StdTransfer<TDomain, TAlgebra>::
assemble_restriction_elemwise(matrix_type& mat,
                               const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
                               ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");

//  resize matrix
	mat.resize_and_clear(coarseDD.num_indices(), fineDD.num_indices());

//	check if grid distribution has dofs, otherwise skip creation since father
//	elements may not exist in parallel.
	if(mat.num_rows() == 0 || mat.num_cols() == 0) return;

	std::vector<DoFIndex> vCoarseDoF, vFineDoF;

//	vector of local finite element ids
	const int dim = TDomain::dim;
	std::vector<LFEID> vLFEID(fineDD.num_fct());
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct){
		vLFEID[fct] = fineDD.lfeid(fct);
		if(vLFEID[fct].dim() != dim)
			UG_THROW("StdTransfer: elem-wise transfer currently only implemented"
					" for non-manifold spaces.")
	}

//  iterators
	MultiGrid& mg = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());
	typedef typename DoFDistribution::dim_traits<dim>::const_iterator const_iterator;
	typedef typename DoFDistribution::dim_traits<dim>::geometric_base_object Element;
	typedef typename Element::side Side;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	for(int si = 0; si < coarseDD.num_subsets(); ++si)
	{
		iterBegin = coarseDD.template begin<Element>(si);
		iterEnd = coarseDD.template end<Element>(si);

	//  loop elems on coarse level for subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Element* coarseElem = *iter;

		//  get children
			const size_t numChild = mg.num_children<Element, Element>(coarseElem);
			if(numChild == 0) continue;

			std::vector<Element*> vChild(numChild);
			for(size_t c = 0; c < vChild.size(); ++c)
				vChild[c] = mg.get_child<Element, Element>(coarseElem,  c);

		//	type of coarse element
			const ReferenceObjectID roid = coarseElem->reference_object_id();

		//	get corner coordinates
			std::vector<MathVector<dim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *coarseElem, *spDomain);

		//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map =
					ReferenceMappingProvider::get<dim, dim>(roid, vCornerCoarse);

		//	loop all components
			for(size_t fct = 0; fct < coarseDD.num_fct(); fct++)
			{
			//	check that fct defined on subset
				if(!coarseDD.is_def_in_subset(fct, si)) continue;

			//  get global indices
				coarseDD.dof_indices(coarseElem, fct, vCoarseDoF);

			// 	piecewise constant elements
				if (vLFEID[fct].type() == LFEID::PIECEWISE_CONSTANT)
				{
					//	loop children
					for(size_t c = 0; c < vChild.size(); ++c)
					{
						Element* child = vChild[c];
						fineDD.dof_indices(child, fct, vFineDoF);
						DoFRef(mat, vCoarseDoF[0], vFineDoF[0]) =  1.0;
					}

					continue;
				} // end of piece-wise constant

			//	get local finite element trial spaces
				const LocalShapeFunctionSet<dim>& lsfs =
						LocalFiniteElementProvider::get<dim>(roid, vLFEID[fct]);

			// 	Crouzeix-Raviart elements
				if (vLFEID[fct].type() == LFEID::CROUZEIX_RAVIART)
				{
				//	all fine sides
					std::vector<Side*> vChildSide;
					std::vector<number> vWeight;

				//	loop side children of coarse sides
				//  if there are coarse neighbour elements on a side the interpolation
				//  is weighted with factor 0.5 else with factor 1
					typename Grid::template traits<Side>::secure_container vCoarseSides;
					mg.associated_elements(vCoarseSides,coarseElem);

					for (size_t i = 0; i < vCoarseSides.size() ;i++){
						Side* side = vCoarseSides[i];
						typename Grid::template traits<Element>::secure_container vCoarseNeighbourElem;
						mg.associated_elements(vCoarseNeighbourElem, side);
						const number weight = (1.0)/vCoarseNeighbourElem.size();

						for(size_t c = 0; c < mg.num_children<Side,Side>(side); ++c){
							vChildSide.push_back(mg.get_child<Side,Side>(side,  c));
							vWeight.push_back(weight);
						}
					}

				//	loop inner child sides
					for(size_t c = 0; c < mg.num_children<Side,Element>(coarseElem); ++c){
						vChildSide.push_back(mg.get_child<Side,Element>(coarseElem,c));
						vWeight.push_back(1.0);
					}

				// 	loop fine sides
					for(size_t c = 0; c < vChildSide.size(); ++c)
					{
						Side* childSide = vChildSide[c];

					//	fine dof indices
						fineDD.inner_dof_indices(childSide, fct, vFineDoF);

					//	global positions of fine dofs
						std::vector<MathVector<dim> > vDoFPos;
						DoFPosition(vDoFPos, childSide, *spDomain, vLFEID[fct]);

					//	get local position of DoF
						std::vector<MathVector<dim> > vLocPos(vDoFPos.size(), 0.0);
						map.update(vCornerCoarse);
						map.global_to_local(vLocPos, vDoFPos);

					//	evaluate coarse shape fct at fine local point
						std::vector<std::vector<number> > vvShape;
						lsfs.shapes(vvShape, vLocPos);

					//	add restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								DoFRef(mat, vCoarseDoF[sh], vFineDoF[ip]) += vWeight[c] * vvShape[ip][sh];
					}

					continue;
				} // end of Crouzeix-Raviart transfer

			// 	Lagrange elements
				if (vLFEID[fct].type() == LFEID::LAGRANGE)
				{
				//	loop children
					for(size_t c = 0; c < vChild.size(); ++c)
					{
						Element* child = vChild[c];

					//	fine dof indices
						fineDD.dof_indices(child, fct, vFineDoF);

					//	global positions of fine dofs
						std::vector<MathVector<dim> > vDoFPos;
						DoFPosition(vDoFPos, child, *spDomain, vLFEID[fct]);

					//	get local position of DoF
					//  update map coordinates because they have been changed in dof_indices
						std::vector<MathVector<dim> > vLocPos(vDoFPos.size(), 0.0);
						map.update(vCornerCoarse);
						map.global_to_local(vLocPos, vDoFPos);

					//	evaluate coarse shape fct at fine local point
						std::vector<std::vector<number> > vvShape;
						lsfs.shapes(vvShape, vLocPos);

					//	set restriction
						for(size_t ip = 0; ip < vvShape.size(); ++ip)
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								DoFRef(mat, vCoarseDoF[sh], vFineDoF[ip]) = vvShape[ip][sh];
					}

					continue;
				} // end of Lagrange transfer

				// other finite element types
				UG_THROW("StdTransfer: Local-Finite-Element: "<<vLFEID[fct]<<" is "
				         "not supported by this Transfer.")
			}
		}
	}
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
			set_identity_on_pure_surface(*R, *spCoarseDD, *spFineDD);
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
