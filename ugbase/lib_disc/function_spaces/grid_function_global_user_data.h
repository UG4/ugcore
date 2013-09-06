/*
 * grid_function_global_user_data.h
 *
 *  Created on: 03.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_GLOBAL_USER_DATA__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_GLOBAL_USER_DATA__

#include "common/common.h"

#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"


namespace ug{



template <typename TGridFunction>
class GlobalGridFunctionNumberData
	: public StdGlobPosData<GlobalGridFunctionNumberData<TGridFunction>, number, TGridFunction::dim>
{
	public:
	///	world dimension of grid function
		static const int dim = TGridFunction::dim;

		private:
	/// grid function
		SmartPtr<TGridFunction> m_spGridFct;

	///	component of function
		size_t m_fct;

	///	local finite element id
		LFEID m_lfeID;

	public:
	/// constructor
		GlobalGridFunctionNumberData(SmartPtr<TGridFunction> spGridFct, const char* cmp)
		: m_spGridFct(spGridFct)
		{
			this->set_functions(cmp);

			//	get function id of name
			m_fct = spGridFct->fct_id_by_name(cmp);

			//	check that function exists
			if(m_fct >= spGridFct->num_fct())
				UG_THROW("GridFunctionNumberData: Function space does not contain"
						" a function with name " << cmp << ".");

			//	local finite element id
			m_lfeID = spGridFct->local_finite_element_id(m_fct);
		};

		virtual bool continuous() const
		{
			return LocalFiniteElementProvider::continuous(m_lfeID);
		}

	///	evaluates the data at a given point and time
		inline void evaluate(number& value, const MathVector<dim>& x, number time, int si) const
		{
			try{
			// \todo: find corresponding element

			GeometricObject* elem;

		//	get corners of element
			std::vector<MathVector<dim> > vCornerCoords;
			CollectCornerCoordinates(vCornerCoords, *elem, *m_spGridFct->domain());

		//	reference object id
			const ReferenceObjectID roid = elem->reference_object_id();

		//	get local position of DoF
			DimReferenceMapping<dim, dim>& map
				= ReferenceMappingProvider::get<dim, dim>(roid, vCornerCoords);
			MathVector<dim> locPos;
			map.global_to_local(locPos, x);

		//	evaluate at shapes at ip
			const LocalShapeFunctionSet<refDim>& rTrialSpace =
					LocalFiniteElementProvider::get<refDim>(roid, m_lfeID);
			std::vector<number> vShape;
			rTrialSpace.shapes(vShape, locPos);

		//	get multiindices of element
			std::vector<MultiIndex<2> > ind;
			m_spGridFct->multi_indices(elem, m_fct, ind);

		// 	compute solution at integration point
			value = 0.0;
			for(size_t sh = 0; sh < vShape.size(); ++sh)
			{
				const number valSH = DoFRef(*m_spGridFct, ind[sh]);
				value += valSH * vShape[sh];
			}

			}
			UG_CATCH_THROW("GlobalGridFunctionNumberData: Evaluation failed.");
		}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_GLOBAL_USER_DATA__ */
