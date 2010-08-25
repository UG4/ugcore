/*
 * constraints_post_process_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__CONSTRAINTS_POST_PROCESS_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__CONSTRAINTS_POST_PROCESS_INTERFACE__

#include "lib_discretization/assemble.h"

namespace ug {

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class IConstraintsPostProcess{
	public:
		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual IAssembleReturn post_process_jacobian(matrix_type& J, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual IAssembleReturn post_process_defect(vector_type& d, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual IAssembleReturn post_process_linear(matrix_type& mat, vector_type& rhs, const discrete_function_type& u)
		{return IAssemble_NOT_IMPLEMENTED;}

		virtual ~IConstraintsPostProcess() {};
};

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class SymP1ConstraintsPostProcess : public IConstraintsPostProcess<TDiscreteFunction, TAlgebra> {
	public:
		// discrete function type
		typedef TDiscreteFunction discrete_function_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		virtual IAssembleReturn post_process_jacobian(matrix_type& J, const discrete_function_type& u)
		{
			typename geometry_traits<ConstrainingEdge>::iterator iter, iterBegin, iterEnd;

			iterBegin = u.template begin<ConstrainingEdge>();
			iterEnd = u.template end<ConstrainingEdge>();

			for(iter = iterBegin; iter != iterEnd; ++iter)
			{
				ConstrainingEdge* cEdge = *iter;

				//ContrainedEdge* Edge =

			}

			return IAssemble_OK;
		}



};


}; // namespace ug



#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__POST_PROCESS__CONSTRAINTS__CONSTRAINTS_POST_PROCESS_INTERFACE__ */
