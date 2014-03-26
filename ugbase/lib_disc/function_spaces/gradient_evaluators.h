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
		static const int dim = TFunction::dim;
		typedef MathVector<dim>						vector_t;
		typedef typename TFunction::element_type 	elem_t;

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
