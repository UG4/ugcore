/*
 * level_transfer.h
 *
 *  Created on: 06.03.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__LEVEL_TRANSFER__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__LEVEL_TRANSFER__

#include "lib_disc/function_spaces/grid_function.h"

namespace ug{

template <typename TDomain, typename TDD, typename TAlgebra>
void ProlongateP1(GridFunction<TDomain, TDD, TAlgebra>& uFine,
                  GridFunction<TDomain, TDD, TAlgebra>& uCoarse)
{
	typedef GridFunction<TDomain, TDD, TAlgebra> TGridFunction;
	typedef typename TGridFunction::template traits<VertexBase>::const_iterator const_iterator;

//  get subsethandler and grid
	SmartPtr<MultiGrid> mg = uFine.domain()->grid();

//	get top level of gridfunctions
	const int fineTopLevel = uFine.dof_distribution()->grid_level().level();
	const int coarseTopLevel = uCoarse.dof_distribution()->grid_level().level();

//	check
	if(fineTopLevel == GridLevel::TOPLEVEL || coarseTopLevel == GridLevel::TOPLEVEL)
		UG_THROW_FATAL("ProlongateP1: Top Level not supported.")
	if(fineTopLevel != coarseTopLevel + 1)
		UG_THROW_FATAL("ProlongateP1: GridFunctions must have one level difference.");

//	storage
	std::vector<size_t> vFineMI, vCoarseMI;

//	loop elements
	const_iterator iterEnd = uFine.template end<VertexBase>();
	const_iterator iter = uFine.template begin<VertexBase>();
	for(; iter != iterEnd; ++iter)
	{
	//	get vertex
		VertexBase* vrt = *iter;

	//	a) 	if not on the same level as the top level of the fine grid function
	//		we can simply copy the values, since the coarse grid function and
	//		the fine grid function are covering the identical part here
		if(mg->get_level(vrt) != fineTopLevel)
		{
			uFine.inner_algebra_indices(vrt, vFineMI);
			uCoarse.inner_algebra_indices(vrt, vCoarseMI);

			for(size_t i = 0; i < vFineMI.size(); ++i)
				uFine[ vFineMI[i] ] = uCoarse[ vCoarseMI[i] ];
		}
		else
		{
		//  get parent
			GeometricObject* parent = mg->get_parent(vrt);

		//	distinguish type of parent
			switch(parent->base_object_type_id())
			{
				case ROID_VERTEX:
				{
					VertexBase* pParent = static_cast<VertexBase*>(parent);
					uFine.inner_algebra_indices(vrt, vFineMI);
					uCoarse.inner_algebra_indices(pParent, vCoarseMI);

					for(size_t i = 0; i < vFineMI.size(); ++i)
						uFine[ vFineMI[i] ] = uCoarse[ vCoarseMI[i] ];
				}
				break;
				case ROID_EDGE:
				{
					uFine.inner_algebra_indices(vrt, vFineMI);
					for(size_t i = 0; i < vFineMI.size(); ++i)
						uFine[ vFineMI[i] ] = 0.0;

					EdgeBase* pParent = static_cast<EdgeBase*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						VertexBase* edgeVrt = pParent->vertex(i);
						uCoarse.inner_algebra_indices(edgeVrt, vCoarseMI);

						for(size_t i = 0; i < vFineMI.size(); ++i)
							VecScaleAdd(uFine[ vFineMI[i] ],
							            1.0, uFine[ vFineMI[i] ],
							            0.5, uCoarse[ vCoarseMI[i] ]);
					}
				}
				break;
				case ROID_QUADRILATERAL:
				{
					uFine.inner_algebra_indices(vrt, vFineMI);
					for(size_t i = 0; i < vFineMI.size(); ++i)
						uFine[ vFineMI[i] ] = 0.0;

					Face* pParent = static_cast<Face*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						VertexBase* faceVrt = pParent->vertex(i);
						uCoarse.inner_algebra_indices(faceVrt, vCoarseMI);

						for(size_t i = 0; i < vFineMI.size(); ++i)
							VecScaleAdd(uFine[ vFineMI[i] ],
							            1.0, uFine[ vFineMI[i] ],
							            0.25, uCoarse[ vCoarseMI[i] ]);
					}
				}
				break;
				case ROID_HEXAHEDRON:
				{
					uFine.inner_algebra_indices(vrt, vFineMI);
					for(size_t i = 0; i < vFineMI.size(); ++i)
						uFine[ vFineMI[i] ] = 0.0;

					Volume* pParent = static_cast<Volume*>(parent);
					for(size_t i = 0; i < pParent->num_vertices(); ++i)
					{
						VertexBase* hexVrt = pParent->vertex(i);
						uCoarse.inner_algebra_indices(hexVrt, vCoarseMI);

						for(size_t i = 0; i < vFineMI.size(); ++i)
							VecScaleAdd(uFine[ vFineMI[i] ],
							            1.0, uFine[ vFineMI[i] ],
							            0.125, uCoarse[ vCoarseMI[i] ]);
					}
				}
				break;
				case ROID_TRIANGLE:
				case ROID_TETRAHEDRON:
				case ROID_PRISM:
				case ROID_PYRAMID: /*nothing to do in those cases */ break;
				default: UG_THROW_FATAL("Unexpected case appeared.");
			}
		}
	}
}


template <typename TDomain, typename TDD, typename TAlgebra>
void Prolongate(GridFunction<TDomain, TDD, TAlgebra>& uFine,
                GridFunction<TDomain, TDD, TAlgebra>& uCoarse)
{
//	grid functions must be from same Domain
	if(uFine.domain().get_impl() != uCoarse.domain().get_impl())
		UG_THROW_FATAL("Prolongate: GridFunctions must have same Domain.");

//	grid functions must have same function pattern
	if(&uFine.function_pattern() != &uCoarse.function_pattern())
		UG_THROW_FATAL("Prolongate: GridFunctions must have same Function Pattern.");

//	loop functions
	bool bOnlyP1Fct = true;
	for(size_t fct = 0; fct < uFine.num_fct(); ++fct)
		if(uFine.local_finite_element_id(fct) != LFEID(LFEID::LAGRANGE, 1))
		{
			bOnlyP1Fct = false; break;
		}

	if(bOnlyP1Fct)
		ProlongateP1(uFine, uCoarse);
	else
	{
		UG_THROW_FATAL("Prolongate: Only implemented for P1.");
	}
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__LEVEL_TRANSFER__ */
