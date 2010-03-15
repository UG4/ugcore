/*
 * system_discretization_interface.h
 *
 *  Created on: 01.03.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SYSTEM_DISCRETIZATION__SYSTEM_DISCRETIZATION_INTERFACE__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SYSTEM_DISCRETIZATION__SYSTEM_DISCRETIZATION_INTERFACE__

namespace ug{

// interface
template <typename TDomain, typename TAlgebra>
class ISystemDiscretization {
	public:
		// forward types and constants

		// domain type
		typedef TDomain domain_type;

		// algebra type
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename TAlgebra::matrix_type matrix_type;

		// type of algebra vector
		typedef typename TAlgebra::vector_type vector_type;

		// type of discrete function
		typedef DiscreteGridFunction<TDomain, TAlgebra> discrete_function_type;


	public:
		virtual bool prepare_element(GeometricObject* elem) = 0;

		virtual bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution<d>& u, number time, number s_m, number s_a) = 0;
		virtual bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution<d>& u, number time, number s_m, number s_a) = 0;
		virtual bool assemble_defect(GeometricObject* elem, Vector& vec, NumericalSolution<d>& u) = 0;
		virtual bool assemble_jacobian(GeometricObject* elem, Matrix& mat, NumericalSolution<d>& u) = 0;
		virtual bool assemble_linear(GeometricObject* elem, Matrix& mat, Vector& rhs, NumericalSolution<d>& u) = 0;

		virtual bool get_dirichlet_values(SubsetHandler& sh, uint subsetIndex, NumericalSolution<d>& u, DirichletValues<d>& dirVal) = 0;

		virtual ~SystemDiscretization(){}
};


} // namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__SUBSET_DISCRETIZATION__SYSTEM_DISCRETIZATION__SYSTEM_DISCRETIZATION_INTERFACE__ */
