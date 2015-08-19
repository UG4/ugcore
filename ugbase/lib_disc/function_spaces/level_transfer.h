/*
 * level_transfer.h
 *
 *  Created on: 06.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__LEVEL_TRANSFER__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__LEVEL_TRANSFER__

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/function_spaces/dof_position_util.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Prolongate
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void ProlongateP1(GridFunction<TDomain, TAlgebra>& uFine,
                  const GridFunction<TDomain, TAlgebra>& uCoarse)
{
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;
	typedef typename TGridFunction::template traits<Vertex>::const_iterator const_iterator;

//  get subsethandler and grid
	SmartPtr<MultiGrid> mg = uFine.domain()->grid();

//	get top level of gridfunctions
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();

//	check
	if(fineTopLevel == GridLevel::TOP || coarseTopLevel == GridLevel::TOP)
		UG_THROW("ProlongateP1: Top Level not supported.")
	if(fineTopLevel < coarseTopLevel)
		UG_THROW("ProlongateP1: fine level must be >= coarse level.");

//	storage
	std::vector<size_t> vFineIndex, vCoarseIndex;

//	loop elements
	const_iterator iterEnd = uFine.template end<Vertex>();
	const_iterator iter = uFine.template begin<Vertex>();
	for(; iter != iterEnd; ++iter)
	{
	//	get vertex
		Vertex* vrt = *iter;
		const int vertexLevel = mg->get_level(vrt);

	//	a) 	if not on the same level as the top level of the fine grid function
	//		and the coarse grid function is already defined on the level
	//		we can simply copy the values, since the coarse grid function and
	//		the fine grid function are covering the identical part here
		if(vertexLevel != fineTopLevel && vertexLevel <= coarseTopLevel)
		{
			uFine.inner_algebra_indices(vrt, vFineIndex);
			uCoarse.inner_algebra_indices(vrt, vCoarseIndex);

			for(size_t i = 0; i < vFineIndex.size(); ++i)
				uFine[ vFineIndex[i] ] = uCoarse[ vCoarseIndex[i] ];

			continue;
		}

	//  get parent and level where coarse grid function is defined
		GridObject* parent = mg->get_parent(vrt);
		const ReferenceObjectID parentBaseObjectID = parent->reference_object_id();
		int parentLevel = mg->get_level(parent);
		while(parentLevel > coarseTopLevel){
			parent = mg->get_parent(parent);
			parentLevel = mg->get_level(parent);
		}

	//	b) 	if the parent, where the coarse grid function is defined and the
	//		fine element are only one level separated, we can use an optimized
	//		interpolation. This case will always apply if the two grid functions
	//		are only one surface level separated.
		if(parentLevel == vertexLevel - 1)
		{
		//	distinguish type of parent
			switch(parentBaseObjectID)
			{
				case ROID_VERTEX:
				{
					Vertex* pParent = static_cast<Vertex*>(parent);
					uFine.inner_algebra_indices(vrt, vFineIndex);
					uCoarse.inner_algebra_indices(pParent, vCoarseIndex);

					for(size_t i = 0; i < vFineIndex.size(); ++i)
						uFine[ vFineIndex[i] ] = uCoarse[ vCoarseIndex[i] ];
				}
				break;
				case ROID_EDGE:
				{
					uFine.inner_algebra_indices(vrt, vFineIndex);
					for(size_t i = 0; i < vFineIndex.size(); ++i)
						uFine[ vFineIndex[i] ] = 0.0;

					Edge* pParent = static_cast<Edge*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						Vertex* edgeVrt = pParent->vertex(i);
						uCoarse.inner_algebra_indices(edgeVrt, vCoarseIndex);

						for(size_t i = 0; i < vFineIndex.size(); ++i)
							VecScaleAdd(uFine[ vFineIndex[i] ],
										1.0, uFine[ vFineIndex[i] ],
										0.5, uCoarse[ vCoarseIndex[i] ]);
					}
				}
				break;
				case ROID_QUADRILATERAL:
				{
					uFine.inner_algebra_indices(vrt, vFineIndex);
					for(size_t i = 0; i < vFineIndex.size(); ++i)
						uFine[ vFineIndex[i] ] = 0.0;

					Face* pParent = static_cast<Face*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						Vertex* faceVrt = pParent->vertex(i);
						uCoarse.inner_algebra_indices(faceVrt, vCoarseIndex);

						for(size_t i = 0; i < vFineIndex.size(); ++i)
							VecScaleAdd(uFine[ vFineIndex[i] ],
										1.0, uFine[ vFineIndex[i] ],
										0.25, uCoarse[ vCoarseIndex[i] ]);
					}
				}
				break;
				case ROID_HEXAHEDRON:
				{
					uFine.inner_algebra_indices(vrt, vFineIndex);
					for(size_t i = 0; i < vFineIndex.size(); ++i)
						uFine[ vFineIndex[i] ] = 0.0;

					Volume* pParent = static_cast<Volume*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						Vertex* hexVrt = pParent->vertex(i);
						uCoarse.inner_algebra_indices(hexVrt, vCoarseIndex);

						for(size_t i = 0; i < vFineIndex.size(); ++i)
							VecScaleAdd(uFine[ vFineIndex[i] ],
										1.0, uFine[ vFineIndex[i] ],
										0.125, uCoarse[ vCoarseIndex[i] ]);
					}
				}
				break;
				case ROID_TRIANGLE:
				case ROID_TETRAHEDRON:
				case ROID_PRISM:
				case ROID_PYRAMID:
				case ROID_OCTAHEDRON: /*nothing to do in those cases */ break;
				default: UG_THROW("Unexpected case appeared.");
			}

			continue;
		}

	//	c) 	we must interpolate the values based on the trial space
		UG_THROW("This case not implemented.");
	}
}



template <typename TDomain, typename TAlgebra>
void ProlongateElemwise(GridFunction<TDomain, TAlgebra>& uFine,
                        const GridFunction<TDomain, TAlgebra>& uCoarse)
{
//	dimension
	const int dim = TDomain::dim;

//  get subsethandler and grid
	SmartPtr<MultiGrid> mg = uFine.domain()->grid();

//	get top level of gridfunctions
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();

//	check
	if(fineTopLevel == GridLevel::TOP || coarseTopLevel == GridLevel::TOP)
		UG_THROW("ProlongateElemwise: Top Level not supported.")
	if(fineTopLevel < coarseTopLevel)
		UG_THROW("ProlongateElemwise: fine level must be >= coarse level.");

//	storage
	std::vector<DoFIndex> vCoarseMI, vFineMI;

//	vector of local finite element ids
	SmartPtr<DoFDistribution> fineDD = uFine.dof_distribution();
	std::vector<LFEID> vFineLFEID(fineDD->num_fct());
	for(size_t fct = 0; fct < fineDD->num_fct(); ++fct)
		vFineLFEID[fct] = fineDD->local_finite_element_id(fct);
	ConstSmartPtr<DoFDistribution> coarseDD = uCoarse.dof_distribution();
	std::vector<LFEID> vCoarseLFEID(coarseDD->num_fct());
	for(size_t fct = 0; fct < coarseDD->num_fct(); ++fct)
		vCoarseLFEID[fct] = coarseDD->local_finite_element_id(fct);

//	check fct
	if(vFineLFEID.size() != vCoarseLFEID.size())
		UG_THROW("ProlongateElemwise: Spaces must contain same number of functions.")

//	get flag if all trial spaces are equal
//	bool bSameLFEID = true;
	for(size_t fct = 0; fct < vFineLFEID.size(); ++fct){
//		if(vFineLFEID[fct] != vCoarseLFEID[fct])
//			bSameLFEID = false;

		if (vCoarseLFEID[fct].type() == LFEID::PIECEWISE_CONSTANT ||
			vFineLFEID[fct].type() == LFEID::PIECEWISE_CONSTANT)
			UG_THROW("Not implemented.")
	}

//  iterators
	typedef typename DoFDistribution::dim_traits<dim>::const_iterator const_iterator;
	typedef typename DoFDistribution::dim_traits<dim>::grid_base_object Element;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	for(int si = 0; si < coarseDD->num_subsets(); ++si)
	{
		iterBegin = coarseDD->template begin<Element>(si);
		iterEnd = coarseDD->template end<Element>(si);

	//  loop elem for coarse level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Element* coarseElem = *iter;

		//  get children where fine grid function is defined
			std::vector<Element*> vChild;
			std::queue<Element*> qElem;
			qElem.push(coarseElem);
			while(!qElem.empty()){
				Element* elem = qElem.front(); qElem.pop();
				if(mg->get_level(elem) == fineTopLevel || !mg->has_children(elem)){
					vChild.push_back(elem);
				} else {
					for(size_t c = 0; c < mg->num_children<Element,Element>(elem); ++c){
						qElem.push(mg->get_child<Element,Element>(elem, c));
					}
				}
			}

		//	type of father
			const ReferenceObjectID coarseROID = coarseElem->reference_object_id();

		//	loop all components
			for(size_t fct = 0; fct < coarseDD->num_fct(); fct++)
			{
			//	check that fct defined on subset
				if(!coarseDD->is_def_in_subset(fct, si)) continue;

			//  get global indices
				coarseDD->dof_indices(coarseElem, fct, vCoarseMI);

			//	get local finite element trial spaces
				const LocalShapeFunctionSet<dim>& lsfs
					= LocalFiniteElementProvider::get<dim>(coarseROID, vCoarseLFEID[fct]);

			//	get corner coordinates
				std::vector<MathVector<dim> > vCornerCoarse;
				CollectCornerCoordinates(vCornerCoarse, *coarseElem, *uFine.domain());

			//	loop children
				for(size_t c = 0; c < vChild.size(); ++c)
				{
					Element* child = vChild[c];

				//	fine dof indices
					fineDD->dof_indices(child, fct, vFineMI);

				//	global positions of fine dofs
					std::vector<MathVector<dim> > vDoFPos, vLocPos;
					DoFPosition(vDoFPos, child, *uFine.domain(), vFineLFEID[fct]);

					UG_ASSERT(vDoFPos.size() == vFineMI.size(), "numDoFPos ("
							  <<vDoFPos.size()<<") != numDoFs ("<<vFineMI.size()<<").");

				//	get Reference Mapping
					DimReferenceMapping<dim, dim>& map
						= ReferenceMappingProvider::get<dim, dim>(coarseROID, vCornerCoarse);


				//	get local position of DoF
					vLocPos.resize(vDoFPos.size());
					for(size_t ip = 0; ip < vLocPos.size(); ++ip) VecSet(vLocPos[ip], 0.0);
					map.global_to_local(vLocPos, vDoFPos);

				//	get all shape functions
					std::vector<std::vector<number> > vvShape;

				//	evaluate coarse shape fct at fine local point
					lsfs.shapes(vvShape, vLocPos);

					for(size_t ip = 0; ip < vvShape.size(); ++ip){
						DoFRef(uFine, vFineMI[ip]) = 0.0;
						for(size_t sh = 0; sh < vvShape[ip].size(); ++sh){
							DoFRef(uFine, vFineMI[ip]) +=
									vvShape[ip][sh] * DoFRef(uCoarse, vCoarseMI[sh]);
						}
					}
				}
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
void Prolongate(GridFunction<TDomain, TAlgebra>& uFine,
                const GridFunction<TDomain, TAlgebra>& uCoarse)
{
//	grid functions must be from same Domain
	if(uFine.domain() != uCoarse.domain())
		UG_THROW("Prolongate: GridFunctions must have same Domain.");

//	grid functions must have same function pattern
	if(uFine.function_pattern().get() != uCoarse.function_pattern().get())
		UG_THROW("Prolongate: GridFunctions must have same Function Pattern.");

//	get grid levels
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();
	if(fineTopLevel == GridLevel::TOP || coarseTopLevel == GridLevel::TOP)
		UG_THROW("Prolongate: Top Level not supported.")
	if(fineTopLevel < coarseTopLevel)
		UG_THROW("Prolongate: fine level must be >= coarse level.");

//	loop functions
	bool bOnlyP1Fct = true;
	for(size_t fct = 0; fct < uFine.num_fct(); ++fct)
		if(uFine.local_finite_element_id(fct).type() != LFEID::LAGRANGE ||
			uFine.local_finite_element_id(fct).order() != 1)
		{
			bOnlyP1Fct = false; break;
		}

	if(bOnlyP1Fct &&
		(fineTopLevel == coarseTopLevel+1 || fineTopLevel == coarseTopLevel)){
		ProlongateP1(uFine, uCoarse);
	}
	else{
		ProlongateElemwise(uFine, uCoarse);
	}

#ifdef UG_PARALLEL
	uFine.set_storage_type(uCoarse.get_storage_mask());
#endif
}

////////////////////////////////////////////////////////////////////////////////
//	Restrict
////////////////////////////////////////////////////////////////////////////////


template <typename TDomain, typename TAlgebra>
void RestrictP1(GridFunction<TDomain, TAlgebra>& uCoarse,
                const GridFunction<TDomain,  TAlgebra>& uFine)
{
	typedef GridFunction<TDomain, TAlgebra> TGridFunction;
	typedef typename TGridFunction::template traits<Vertex>::const_iterator const_iterator;

//  get subsethandler and grid
	SmartPtr<MultiGrid> mg = uCoarse.domain()->grid();

//	get top level of gridfunctions
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();

//	check
	if(fineTopLevel == GridLevel::TOP || coarseTopLevel == GridLevel::TOP)
		UG_THROW("RestrictP1: Top Level not supported.")
	if(fineTopLevel < coarseTopLevel)
		UG_THROW("RestrictP1: fine level must be >= coarse level.");

//	storage
	std::vector<size_t> vFineIndex, vCoarseIndex;

//	loop elements
	const_iterator iterEnd = uCoarse.template end<Vertex>();
	const_iterator iter = uCoarse.template begin<Vertex>();
	for(; iter != iterEnd; ++iter)
	{
	//	get vertex
		Vertex* coarseVrt = *iter;

	//  get children where fine grid function is defined
		Vertex* fineVrt = coarseVrt;
		while(mg->get_level(fineVrt) != fineTopLevel &&
				mg->has_children(fineVrt)){
			fineVrt = mg->get_child<Vertex,Vertex>(fineVrt, 0);
		}

	//	copy values
		uFine.inner_algebra_indices(fineVrt, vFineIndex);
		uCoarse.inner_algebra_indices(coarseVrt, vCoarseIndex);

		for(size_t i = 0; i < vFineIndex.size(); ++i)
			uCoarse[ vCoarseIndex[i] ] = uFine[ vFineIndex[i] ];
	}
}



template <typename TDomain, typename TAlgebra>
void RestrictElemwise(GridFunction<TDomain, TAlgebra>& uCoarse,
                      const GridFunction<TDomain, TAlgebra>& uFine)
{
//	dimension
	const int dim = TDomain::dim;
	const int locDim = TDomain::dim;

//  get subsethandler and grid
	SmartPtr<MultiGrid> mg = uCoarse.domain()->grid();

//	get top level of gridfunctions
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();

//	check
	if(fineTopLevel == GridLevel::TOP || coarseTopLevel == GridLevel::TOP)
		UG_THROW("RestrictElemwise: Top Level not supported.")
	if(fineTopLevel < coarseTopLevel)
		UG_THROW("RestrictElemwise: fine level must be >= coarse level.");

//	storage
	std::vector<DoFIndex> vCoarseMI, vFineMI;

//	vector of local finite element ids
	ConstSmartPtr<DoFDistribution> fineDD = uFine.dof_distribution();
	std::vector<LFEID> vFineLFEID(fineDD->num_fct());
	for(size_t fct = 0; fct < fineDD->num_fct(); ++fct)
		vFineLFEID[fct] = fineDD->local_finite_element_id(fct);
	SmartPtr<DoFDistribution> coarseDD = uCoarse.dof_distribution();
	std::vector<LFEID> vCoarseLFEID(coarseDD->num_fct());
	for(size_t fct = 0; fct < coarseDD->num_fct(); ++fct)
		vCoarseLFEID[fct] = coarseDD->local_finite_element_id(fct);

//	check fct
	if(vFineLFEID.size() != vCoarseLFEID.size())
		UG_THROW("RestrictElemwise: Spaces must contain same number of functions.")

//	get flag if all trial spaces are equal
	bool bSameLFEID = true;
	for(size_t fct = 0; fct < vFineLFEID.size(); ++fct){
		if(vFineLFEID[fct] != vCoarseLFEID[fct])
			bSameLFEID = false;

		if (vCoarseLFEID[fct].type() == LFEID::PIECEWISE_CONSTANT ||
			vFineLFEID[fct].type() == LFEID::PIECEWISE_CONSTANT)
			UG_THROW("Not implemented.")
	}

//  iterators
	typedef typename DoFDistribution::dim_traits<locDim>::const_iterator const_iterator;
	typedef typename DoFDistribution::dim_traits<locDim>::grid_base_object Element;
	const_iterator iter, iterBegin, iterEnd;

//  loop subsets on coarse level
	for(int si = 0; si < coarseDD->num_subsets(); ++si)
	{
		iterBegin = coarseDD->template begin<Element>(si);
		iterEnd = coarseDD->template end<Element>(si);

	//  loop elem for coarse level subset
		for(iter = iterBegin; iter != iterEnd; ++iter)
		{
		//	get element
			Element* coarseElem = *iter;

		//  get children where fine grid function is defined
			std::vector<Element*> vFineElem;
			std::queue<Element*> qElem;
			qElem.push(coarseElem);
			while(!qElem.empty()){
				Element* elem = qElem.front(); qElem.pop();
				if(mg->get_level(elem) == fineTopLevel || !mg->has_children(elem)){
					vFineElem.push_back(elem);
				} else {
					for(size_t c = 0; c < mg->num_children<Element,Element>(elem); ++c){
						qElem.push(mg->get_child<Element,Element>(elem, c));
					}
				}
			}

		//	loop all components
			for(size_t fct = 0; fct < coarseDD->num_fct(); fct++)
			{
			//	check that fct defined on subset
				if(!coarseDD->is_def_in_subset(fct, si)) continue;

			//  get global indices
				coarseDD->inner_dof_indices(coarseElem, fct, vCoarseMI);

			//	global positions of fine dofs
				std::vector<MathVector<dim> > vDoFPos;
				InnerDoFPosition(vDoFPos, coarseElem, *uCoarse.domain(), vCoarseLFEID[fct]);

			//	loop dof points
				for(size_t ip = 0; ip < vDoFPos.size(); ++ip)
				{
				//	loop children
					for(size_t c = 0; c < vFineElem.size(); ++c)
					{
						Element* fineElem = vFineElem[c];

						UG_THROW("This part does not work.");
/*						if(!ContainsPoint(fineElem, vDoFPos[ip], uFine.domain()->position_accessor())){
							if(c == vFineElem.size()-1)
								UG_THROW("Restrict: Cannot find child containing dof.");
							continue;
						}
*/
					//	get corner coordinates
						std::vector<MathVector<dim> > vCornerFine;
						CollectCornerCoordinates(vCornerFine, *fineElem, *uFine.domain());

					//	type of child
						const ReferenceObjectID fineROID = fineElem->reference_object_id();

					//	get local finite element trial spaces
						const LocalShapeFunctionSet<locDim>& lsfs
							= LocalFiniteElementProvider::get<locDim>(fineROID, vFineLFEID[fct]);

					//	get Reference Mapping
						DimReferenceMapping<locDim, dim>& map
							= ReferenceMappingProvider::get<locDim, dim>(fineROID, vCornerFine);

					//	get local position of DoF
						MathVector<locDim> vLocPos;
						VecSet(vLocPos, 0.0);
						map.global_to_local(vLocPos, vDoFPos[ip]);

					//	fine dof indices
						fineDD->dof_indices(fineElem, fct, vFineMI);

					//	get all shape functions
						std::vector<number> vShape;

					//	evaluate coarse shape fct at fine local point
						lsfs.shapes(vShape, vLocPos);

					//	interpolate
						DoFRef(uCoarse, vCoarseMI[ip]) = 0.0;
						for(size_t sh = 0; sh < vShape.size(); ++sh)
						{
							DoFRef(uCoarse, vCoarseMI[ip]) +=
									vShape[sh] * DoFRef(uFine, vFineMI[sh]);
						}
					}
				}
			}
		}
	}
}


template <typename TDomain, typename TAlgebra>
void Restrict(GridFunction<TDomain, TAlgebra>& uCoarse,
              const GridFunction<TDomain, TAlgebra>& uFine)
{
//	grid functions must be from same Domain
	if(uCoarse.domain() != uFine.domain())
		UG_THROW("Restrict: GridFunctions must have same Domain.");

//	grid functions must have same function pattern
	if(uCoarse.function_pattern().get() != uFine.function_pattern().get())
		UG_THROW("Restrict: GridFunctions must have same Function Pattern.");

//	get grid levels
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	if(coarseTopLevel == GridLevel::TOP || fineTopLevel == GridLevel::TOP)
		UG_THROW("Restrict: Top Level not supported.")
	if(coarseTopLevel > fineTopLevel)
		UG_THROW("Restrict: fine level must be >= coarse level.");

//	loop functions
	bool bOnlyP1Fct = true;
	for(size_t fct = 0; fct < uCoarse.num_fct(); ++fct)
		if(uCoarse.local_finite_element_id(fct).type() != LFEID::LAGRANGE ||
				uCoarse.local_finite_element_id(fct).order() != 1)
		{
			bOnlyP1Fct = false; break;
		}

	if(bOnlyP1Fct &&
		(coarseTopLevel+1 == fineTopLevel || coarseTopLevel == fineTopLevel)){
		RestrictP1(uCoarse, uFine);
	}
	else{
		UG_THROW("Restrict: Only P1 implemented.")
		RestrictElemwise(uCoarse, uFine);
	}

#ifdef UG_PARALLEL
	uCoarse.set_storage_type(uFine.get_storage_mask());
#endif
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__LEVEL_TRANSFER__ */
