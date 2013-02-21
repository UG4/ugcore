/*
 * assembled_non_linear_operator.h
 *
 *  Created on: ..
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__

#include "lib_algebra/operator/interface/operator.h"

namespace ug{

template <typename TAlgebra>
class AssembledOperator : public IOperator<typename TAlgebra::vector_type>
{
public:
	/// Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	default constructor
		AssembledOperator()
			: m_pAss(NULL), m_gridLevel() {};

	///	constructor
		AssembledOperator(IAssemble<TAlgebra>& ass)
			: m_pAss(&ass), m_gridLevel(){};

	///	constructor
		AssembledOperator(IAssemble<TAlgebra>& ass, const GridLevel& gl)
			: m_pAss(&ass), m_gridLevel(gl) {};

	///	sets discretization for assembling
		void set_discretization(IAssemble<TAlgebra>& ass) {m_pAss = &ass;}

	///	sets the level used for assembling
		void set_level(const GridLevel& gl) {m_gridLevel = gl;}

	///	returns the level used for assembling
		const GridLevel& level() const {return m_gridLevel;}

	///	Init
		virtual void init() {}

	///	Prepare for apply
		virtual void prepare(vector_type& u);

	/// Compute d = L(u)
		virtual void apply(vector_type& d, const vector_type& u);

	/// return assembling
		IAssemble<TAlgebra>* get_assemble() {return m_pAss;}

	protected:
	///	assembling procedure
		IAssemble<TAlgebra>* m_pAss;

	///	used grid level
		GridLevel m_gridLevel;
};

} // end namepace ug

// include implementation
#include "assembled_non_linear_operator_impl.h"

#endif /*__H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__*/
