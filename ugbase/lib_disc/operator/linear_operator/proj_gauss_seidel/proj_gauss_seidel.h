/*
 * proj_gauss_seidel.h
 *
 *  Created on: 10.10.2013
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJ_GAUSS_SEIDEL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJ_GAUSS_SEIDEL__

namespace ug{

template <typename TDomain, typename TAlgebra>
class ProjGaussSeidel:
	public ILinearIterator<typename TAlgebra::vector_type>
{

	public:
	///	Domain
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

	public:
		/// constructor
		ProjGaussSeidel(){};

	///////////////////////////////////////////////////////////////////////////
	//	Linear Solver interface methods
	///////////////////////////////////////////////////////////////////////////

	///	name
		virtual const char* name() const {return "Projected GaussSeidel";}

	/// Prepare for Operator J(u) and linearization point u (current solution)
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

	///	Prepare for Linear Operartor L
		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L);

	///	Compute new correction c = B*d
		virtual bool apply(vector_type& c, const vector_type& d);

	///	Compute new correction c = B*d and return new defect d := d - A*c
		virtual bool apply_update_defect(vector_type& c, vector_type& d);

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone();

		///	Destructor
		~ProjGaussSeidel()
		{};

};

} // end namespace ug

// include implementation
#include "proj_gauss_seidel_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJ_GAUSS_SEIDEL__ */
