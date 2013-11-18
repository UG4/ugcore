/*
 * prolongation_operator_impl.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR_IMPL__

#include "prolongation_operator.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/function_spaces/grid_function_util.h"

namespace ug{


template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
assemble_restriction_p1(typename TAlgebra::matrix_type& mat,
                         const DoFDistribution& coarseDD, const DoFDistribution& fineDD)
{
	PROFILE_FUNC_GROUP("gmg");
// 	allow only lagrange P1 functions
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		if(fineDD.lfeid(fct).type() != LFEID::LAGRANGE ||
			fineDD.lfeid(fct).order() != 1)
			UG_THROW("AssembleStdProlongationForP1Lagrange:"
				"Interpolation only implemented for Lagrange P1 functions.");

//  get subsethandler and grid
	const MultiGrid& grid = *coarseDD.multi_grid();

//  get number of dofs on different levels
	const size_t numFineDoFs = fineDD.num_indices();
	const size_t numCoarseDoFs = coarseDD.num_indices();

//  resize matrix
	mat.resize_and_clear(numCoarseDoFs, numFineDoFs);

//	check if grid distribution has dofs, otherwise skip creation since father
//	elements may not exist in parallel.
	if(numFineDoFs == 0 || numCoarseDoFs == 0) return;

	std::vector<DoFIndex> coarseMultInd, fineMultInd;

//  iterators
	typedef DoFDistribution::traits<VertexBase>::const_iterator const_iterator;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
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
			GeometricObject* parent = grid.get_parent(fineVrt);

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
				fineDD.inner_dof_indices(fineVrt, fct, fineMultInd);

			//	detect type of father
				switch(roid)
				{
					case ROID_VERTEX:
						{
							VertexBase* vrt = dynamic_cast<VertexBase*>(parent);
							coarseDD.inner_dof_indices(vrt, fct, coarseMultInd);
							DoFRef(mat, coarseMultInd[0], fineMultInd[0]) = 1.0;
						}
						break;
					case ROID_EDGE:
						for(int i = 0; i < 2; ++i)
						{
							EdgeBase* edge = dynamic_cast<EdgeBase*>(parent);
							coarseDD.inner_dof_indices(edge->vertex(i), fct, coarseMultInd);
							DoFRef(mat, coarseMultInd[0], fineMultInd[0]) = 0.5;
						}
						break;
					case ROID_QUADRILATERAL:
						for(int i = 0; i < 4; ++i)
						{
							Face* face = dynamic_cast<Face*>(parent);
							coarseDD.inner_dof_indices(face->vertex(i), fct, coarseMultInd);
							DoFRef(mat, coarseMultInd[0], fineMultInd[0]) = 0.25;
						}
						break;
					case ROID_HEXAHEDRON:
						for(int i = 0; i < 8; ++i)
						{
							Volume* hexaeder = dynamic_cast<Volume*>(parent);
							coarseDD.inner_dof_indices(hexaeder->vertex(i), fct, coarseMultInd);
							DoFRef(mat, coarseMultInd[0], fineMultInd[0]) = 0.125;
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
assemble_restriction_elemwise(typename TAlgebra::matrix_type& mat,
                               const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
                               ConstSmartPtr<TDomain> spDomain)
{
	PROFILE_FUNC_GROUP("gmg");
//	dimension
	const int dim = TDomain::dim;

//  get subsethandler and grid
	MultiGrid& grid = *const_cast<MultiGrid*>(coarseDD.multi_grid().get());

//  get number of dofs on different levels
	const size_t numFineDoFs = fineDD.num_indices();
	const size_t numCoarseDoFs = coarseDD.num_indices();

//  resize matrix
	mat.resize_and_clear(numCoarseDoFs, numFineDoFs);

//	check if grid distribution has dofs, otherwise skip creation since father
//	elements may not exist in parallel.
	if(numFineDoFs == 0 || numCoarseDoFs == 0) return;

	std::vector<DoFIndex> vCoarseMultInd, vFineMultInd;

//	vector of local finite element ids
	std::vector<LFEID> vLFEID(fineDD.num_fct());
	for(size_t fct = 0; fct < fineDD.num_fct(); ++fct)
		vLFEID[fct] = fineDD.lfeid(fct);

//  iterators
	typedef typename DoFDistribution::dim_traits<dim>::const_iterator const_iterator;
	typedef typename DoFDistribution::dim_traits<dim>::geometric_base_object Element;
	typedef typename Element::side Side;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on fine level
	for(int si = 0; si < coarseDD.num_subsets(); ++si)
	{
		iterBegin = coarseDD.template begin<Element>(si);
		iterEnd = coarseDD.template end<Element>(si);

	//  loop vertices for fine level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Element* coarseElem = *iter;

		//  get children
			const size_t numChild = grid.num_children<Element, Element>(coarseElem);
			if(numChild == 0) continue;

			std::vector<Element*> vChild(numChild);
			for(size_t c = 0; c < vChild.size(); ++c)
				vChild[c] = grid.get_child<Element, Element>(coarseElem,  c);

		//	type of coarse element
			const ReferenceObjectID roid = coarseElem->reference_object_id();

		//	get corner coordinates
			std::vector<MathVector<dim> > vCornerCoarse;
			CollectCornerCoordinates(vCornerCoarse, *coarseElem, *spDomain);

		//	get Reference Mapping
			DimReferenceMapping<dim, dim>& map = ReferenceMappingProvider::get<dim, dim>(roid, vCornerCoarse);

		//	loop all components
			for(size_t fct = 0; fct < coarseDD.num_fct(); fct++)
			{
				//	check that fct defined on subset
				if(!coarseDD.is_def_in_subset(fct, si)) continue;
				//  get global indices
				coarseDD.dof_indices(coarseElem, fct, vCoarseMultInd);

				// piecewise constant elements
				if (vLFEID[fct].type() == LFEID::PIECEWISE_CONSTANT){
					//	loop children
					for(size_t c = 0; c < vChild.size(); ++c)
					{
						Element* child = vChild[c];
						//	fine dof indices
						fineDD.dof_indices(child, fct, vFineMultInd);
						DoFRef(mat, vCoarseMultInd[0], vFineMultInd[0]) =  1.0;
					}
					continue;
				}

				//	get local finite element trial spaces
				const LocalShapeFunctionSet<dim>& lsfs = LocalFiniteElementProvider::get<dim>(roid, vLFEID[fct]);

				// Crouzeix-Raviart elements
				if (vLFEID[fct].type() == LFEID::CROUZEIX_RAVIART)
				{
					//	loop side children of coarse sides
					//  if there are coarse neighbour elements on a side the interpolation
					//  is weighted with factor 0.5 else with factor 1
					typename Grid::template traits<Side>::secure_container vCoarseSides;
					grid.associated_elements(vCoarseSides,coarseElem);

					for (size_t i=0;i<vCoarseSides.size();i++)
					{
						Side* side = vCoarseSides[i];
						typename Grid::template traits<Element>::secure_container vCoarseNeighbourElem;
						grid.associated_elements(vCoarseNeighbourElem,side);
						number weight = (1.0)/vCoarseNeighbourElem.size();

						const size_t numChild = grid.num_children<Side,Side>(side);
						if(numChild == 0) continue;

						for(size_t c = 0; c < numChild; ++c)
						{
							Side* childSide = grid.get_child<Side,Side>(side,  c);
							//	fine dof indices
							fineDD.inner_dof_indices(childSide, fct, vFineMultInd);
							//	global positions of fine dofs
							std::vector<MathVector<dim> > vDoFPos, vLocPos;
							DoFPosition(vDoFPos, childSide, *spDomain, vLFEID[fct]);

							//	get local position of DoF
							vLocPos.resize(vDoFPos.size());
							for(size_t ip = 0; ip < vLocPos.size(); ++ip) VecSet(vLocPos[ip], 0.0);

							map.update(vCornerCoarse);
							map.global_to_local(vLocPos, vDoFPos);

							//	get all shape functions
							std::vector<std::vector<number> > vvShape;

							//	evaluate coarse shape fct at fine local point
							lsfs.shapes(vvShape, vLocPos);

							for(size_t ip = 0; ip < vvShape.size(); ++ip)
							{
								for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
								{
									DoFRef(mat, vCoarseMultInd[sh], vFineMultInd[ip]) += weight*vvShape[ip][sh];
								}
							}
						}
					}

					// loop fine sides inside coarse element
					const size_t numChild = grid.num_children<Side,Element>(coarseElem);
					for(size_t c = 0; c < numChild; ++c){
						Side* childSide = grid.get_child<Side,Element>(coarseElem,c);
						//	fine dof indices
						fineDD.inner_dof_indices(childSide, fct, vFineMultInd);
						//	global positions of fine dofs
						std::vector<MathVector<dim> > vDoFPos, vLocPos;
						DoFPosition(vDoFPos, childSide, *spDomain, vLFEID[fct]);

						//	get local position of DoF
						vLocPos.resize(vDoFPos.size());
						for(size_t ip = 0; ip < vLocPos.size(); ++ip) VecSet(vLocPos[ip], 0.0);

						map.update(vCornerCoarse);
						map.global_to_local(vLocPos, vDoFPos);

						//	get all shape functions
						std::vector<std::vector<number> > vvShape;

						//	evaluate coarse shape fct at fine local point
						lsfs.shapes(vvShape, vLocPos);

						for(size_t ip = 0; ip < vvShape.size(); ++ip)
						{
							for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
							{
								DoFRef(mat, vCoarseMultInd[sh], vFineMultInd[ip]) = vvShape[ip][sh];
							}
						}
					}
					continue;
				} // end of Crouzeix-Raviart transfer

			// other finite element types

			//	loop children
				for(size_t c = 0; c < vChild.size(); ++c)
				{
					Element* child = vChild[c];

				//	fine dof indices
					fineDD.dof_indices(child, fct, vFineMultInd);

				//	global positions of fine dofs
					std::vector<MathVector<dim> > vDoFPos, vLocPos;
					DoFPosition(vDoFPos, child, *spDomain, vLFEID[fct]);

					UG_ASSERT(vDoFPos.size() == vFineMultInd.size(), "numDoFPos ("
					          <<vDoFPos.size()<<") != numDoFs ("<<vFineMultInd.size()<<").");

				//	get local position of DoF
					vLocPos.resize(vDoFPos.size());
					for(size_t ip = 0; ip < vLocPos.size(); ++ip) VecSet(vLocPos[ip], 0.0);

				//  update map coordinates because they have been changed in dof_indices
					map.update(vCornerCoarse);
					map.global_to_local(vLocPos, vDoFPos);

				//	get all shape functions
					std::vector<std::vector<number> > vvShape;

				//	evaluate coarse shape fct at fine local point
					lsfs.shapes(vvShape, vLocPos);

					for(size_t ip = 0; ip < vvShape.size(); ++ip)
					{
						for(size_t sh = 0; sh < vvShape[ip].size(); ++sh)
						{
							DoFRef(mat, vCoarseMultInd[sh], vFineMultInd[ip]) = vvShape[ip][sh];
						}
					}
				}
			}
		}
	}
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace)
{
	m_spApproxSpace = approxSpace;
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::set_levels(GridLevel coarseLevel, GridLevel fineLevel)
{
	m_fineLevel = fineLevel;
	m_coarseLevel = coarseLevel;

	if(m_fineLevel.level() - m_coarseLevel.level() != 1)
		UG_THROW("StdTransfer<TDomain, TAlgebra>::set_levels:"
				" Can only project between successive level.");

	if(m_fineLevel.type() != GridLevel::LEVEL ||
	   m_coarseLevel.type() != GridLevel::LEVEL)
		UG_THROW("StdTransfer<TDomain, TAlgebra>::set_levels:"
				" Can only project between level dof distributions, but fine="
				<<m_fineLevel<<", coarse="<<m_coarseLevel);
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::init()
{
	PROFILE_FUNC_GROUP("gmg");
	if(!m_spApproxSpace.valid())
		UG_THROW("StdTransfer<TDomain, TAlgebra>::init: "
				"Approximation Space not set. Cannot init Projection.");

	m_Restriction.resize_and_clear(0,0);

// 	check only lagrange P1 functions
	bool P1LagrangeOnly = true;
	for(size_t fct = 0; fct < m_spApproxSpace->num_fct(); ++fct)
		if(m_spApproxSpace->lfeid(fct).type() != LFEID::LAGRANGE ||
			m_spApproxSpace->lfeid(fct).order() != 1)
			P1LagrangeOnly = false;

	try{
		if(P1LagrangeOnly)
		{
			assemble_restriction_p1
			(m_Restriction,
			 *m_spApproxSpace->dof_distribution(m_coarseLevel),
			 *m_spApproxSpace->dof_distribution(m_fineLevel));
		}
		else
		{
			assemble_restriction_elemwise
			(m_Restriction,
			 *m_spApproxSpace->dof_distribution(m_coarseLevel),
			 *m_spApproxSpace->dof_distribution(m_fineLevel),
			 m_spApproxSpace->domain());
		}
	} UG_CATCH_THROW("StdTransfer<TDomain, TAlgebra>::init:"
				"Cannot assemble interpolation matrix.");

	#ifdef UG_PARALLEL
		m_Restriction.set_storage_type(PST_CONSISTENT);
	#endif

	std::stringstream ss; ss<<"Prolongation_"<<m_coarseLevel.level()<<"_"<<m_fineLevel.level();
	write_debug(m_Restriction, ss.str().c_str());

	m_bInit = true;
}

template <typename TDomain, typename TAlgebra>
void StdTransfer<TDomain, TAlgebra>::
prolongate(vector_type& uFine, const vector_type& uCoarse)
{
	PROFILE_FUNC_GROUP("gmg");
//	Check, that operator is initiallized
	if(!m_bInit)
		UG_THROW("StdTransfer<TDomain, TAlgebra>::apply:"
				" Operator not initialized.");

//	Some Assertions
	if(uFine.size() != m_Restriction.num_cols())
		UG_THROW("StdTransfer: Vector ["<<uFine.size()<<"] must be == Cols size "
		         <<m_Restriction.num_cols());
	if(uCoarse.size() != m_Restriction.num_rows())
		UG_THROW("StdTransfer: Vector ["<<uCoarse.size()<<"] must be == Rows size "
		         <<m_Restriction.num_rows());

//	Apply Matrix
	m_Restriction.apply_transposed(uFine, uCoarse);

//	Set dirichlet nodes to zero again
//	todo: We could handle this by eliminating dirichlet rows as well
	try{
	for(size_t i = 0; i < m_vConstraint.size(); ++i){
		if (m_vConstraint[i]->type() & CT_DIRICHLET){
			m_vConstraint[i]->adjust_defect(uFine, uFine, m_spApproxSpace->dof_distribution(m_fineLevel));
		}
	}
	}UG_CATCH_THROW("StdTransfer<TDomain, TAlgebra>::apply: "
					"Error while setting dirichlet defect to zero.");

// call prolongations due to added constraints (= adjust_restrict, member of class constraint)
	try{
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		m_vConstraint[i]->adjust_prolongation(uFine, m_fineLevel, uCoarse, m_coarseLevel);
	} UG_CATCH_THROW("ProjectionOperator::apply_transposed: "
					"Error while setting dirichlet defect to zero.");


// 	check CR functions
#ifdef UG_PARALLEL
	bool bCROnly = true;
	for(size_t fct = 0; fct < m_spApproxSpace->num_fct(); ++fct)
		if(m_spApproxSpace->lfeid(fct).type() != LFEID::CROUZEIX_RAVIART &&
			m_spApproxSpace->lfeid(fct).type() != LFEID::PIECEWISE_CONSTANT)
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
do_restrict(vector_type& uCoarse, const vector_type& uFine)
{
	PROFILE_FUNC_GROUP("gmg");
//	Check, that operator is initialized
	if(!m_bInit)
		UG_THROW("StdTransfer<TDomain, TAlgebra>::apply_transposed:"
				"Operator not initialized.");

//	Some Checks
	if(uFine.size() != m_Restriction.num_cols())
		UG_THROW("StdTransfer: Vector ["<<uFine.size()<<"] must be == Cols size "
		         <<m_Restriction.num_cols());
	if(uCoarse.size() != m_Restriction.num_rows())
		UG_THROW("StdTransfer: Vector ["<<uCoarse.size()<<"] must be == Rows size "
		         <<m_Restriction.num_rows());

//	Apply transposed matrix
	m_Restriction.apply_ignore_zero_rows(uCoarse, m_dampRes, uFine);

//	Set dirichlet nodes to zero again
//	todo: We could handle this by eliminating dirichlet columns as well
	try{
	for(size_t i = 0; i < m_vConstraint.size(); ++i){
		if (m_vConstraint[i]->type() & CT_DIRICHLET){
			m_vConstraint[i]->adjust_defect(uCoarse, uCoarse, m_spApproxSpace->dof_distribution(m_coarseLevel));
		}
	}
	} UG_CATCH_THROW("ProjectionOperator::apply_transposed: "
					"Error while setting dirichlet defect to zero.");

// call restrictions due to added constraints (= adjust_restrict, member of class constraint)
	try{
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		m_vConstraint[i]->adjust_restriction(uCoarse, m_coarseLevel, uFine, m_fineLevel);
	} UG_CATCH_THROW("ProjectionOperator::apply_transposed: "
					"Error while setting dirichlet defect to zero.");

}

template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferOperator<TAlgebra> >
StdTransfer<TDomain, TAlgebra>::clone()
{
	SmartPtr<StdTransfer> op(new StdTransfer);
	op->set_approximation_space(m_spApproxSpace);
	for(size_t i = 0; i < m_vConstraint.size(); ++i)
		op->add_constraint(m_vConstraint[i]);
	op->set_restriction_damping(m_dampRes);
	op->set_debug(m_spDebugWriter);
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
write_debug(const matrix_type& mat, const char* filename)
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
	std::string name(filename); name.append(".mat");

//	write
	GridLevel gridLev = dbgWriter->grid_level();
	dbgWriter->set_grid_levels(m_coarseLevel, m_fineLevel);
	dbgWriter->write_matrix(mat, name.c_str());
	dbgWriter->set_grid_level(gridLev);
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR_IMPL__ */
