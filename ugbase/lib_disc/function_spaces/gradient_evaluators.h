/*
 * Copyright (c) 2014:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_DISC__GRADIENT_EVALUATORS__
#define __H__UG_DISC__GRADIENT_EVALUATORS__

#include "grid_function.h"
#include "common/util/provider.h"
#include "lib_disc/common/geometry_util.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"

namespace ug{

/**	Provides a function to evaluate the gradient of a given grid function efficiently
 * in each element.*/
template <class TFunction>
class GradientEvaluator_LagrangeP1{
	public:
		static constexpr int dim = TFunction::dim;
		using vector_t = MathVector<dim>;
		using elem_t = typename TFunction::element_type;

		GradientEvaluator_LagrangeP1(TFunction* u, size_t fct)	:
			m_pu(u),
			m_aaPos(u->domain()->position_accessor()),
			m_fct(fct)
		{}

		vector_t evaluate(elem_t* elem)
		{
			TFunction& u = *m_pu;

		//	reference object type
			ReferenceObjectID roid = elem->reference_object_id();

		//	get trial space
			const LocalShapeFunctionSet<dim>& lsfs =
					LocalFiniteElementProvider::get<dim>(roid, LFEID(LFEID::LAGRANGE, dim, 1));

		//	create a reference mapping
			DimReferenceMapping<dim, dim>& map
				= ReferenceMappingProvider::get<dim, dim>(roid);

		//	get local Mid Point
			vector_t localIP = ReferenceElementCenter<dim>(roid);

		//	number of shape functions
			const size_t numSH = lsfs.num_sh();
			vLocalGrad.resize(numSH);
			vGlobalGrad.resize(numSH);

		//	evaluate reference gradient at local midpoint
			lsfs.grads(&vLocalGrad[0], localIP);

		//	get corners of element
			CollectCornerCoordinates(vCorner, *elem, m_aaPos);

		//	update mapping
			map.update(&vCorner[0]);

		//	compute jacobian
			map.jacobian_transposed_inverse(JTInv, localIP);

		//	compute gradient at mid point by summing contributions of all shape fct
			vector_t elemGrad;
			VecSet(elemGrad, 0.0);
			for(size_t sh = 0 ; sh < numSH; ++sh)
			{
			//	get global Gradient
				MatVecMult(vGlobalGrad[sh], JTInv, vLocalGrad[sh]);

			//	get vertex
				Vertex* vert = elem->vertex(sh);

			//	get of of vertex
				std::vector<DoFIndex> ind;
				u.inner_dof_indices(vert, m_fct, ind);

			//	scale global gradient
				vGlobalGrad[sh] *= DoFRef(u, ind[0]);

			//	sum up
				elemGrad += vGlobalGrad[sh];
			}
			return elemGrad;
		}
	private:
		TFunction* m_pu;
		typename TFunction::domain_type::position_accessor_type m_aaPos;
		size_t m_fct;
	//	the following members are declared here so that they can be efficiently reused
		MathMatrix<dim, dim> JTInv;
		std::vector<vector_t > vLocalGrad;
		std::vector<vector_t > vGlobalGrad;
		std::vector<vector_t > vCorner;
};

}//	end of namespace
#endif
