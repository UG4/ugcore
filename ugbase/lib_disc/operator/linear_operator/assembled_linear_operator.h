
#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

///	matrix operator based on the assembling of a problem
/**
 * This operator implements the MatrixOperator interface, thus is basically a
 * matrix that can be applied to vectors. In addition the class allows to set
 * an IAssemble object and an the GridLevel. Invoking the init method the
 * matrix is created using the IAssemble routine on the given GridLevel.
 *
 * \tparam	TAlgebra			algebra type
 */
template <typename TAlgebra>
class AssembledLinearOperator :
	public virtual MatrixOperator<	typename TAlgebra::matrix_type,
									typename TAlgebra::vector_type>
{
	public:
	///	Type of Algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Matrix
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Type of base class
		typedef MatrixOperator<matrix_type,vector_type> base_type;

	public:
	///	Default Constructor
		AssembledLinearOperator() :	m_spAss(NULL) {};

	///	Constructor
		AssembledLinearOperator(SmartPtr<IAssemble<TAlgebra> > ass) : m_spAss(ass) {};

	///	Constructor
		AssembledLinearOperator(SmartPtr<IAssemble<TAlgebra> > ass, const GridLevel& gl)
			: m_spAss(ass), m_gridLevel(gl) {};

	///	sets the discretization to be used
		void set_discretization(SmartPtr<IAssemble<TAlgebra> > ass) {m_spAss = ass;}

	///	returns the discretization to be used
		SmartPtr<IAssemble<TAlgebra> > discretization() {return m_spAss;}

	///	sets the level used for assembling
		void set_level(const GridLevel& gl) {m_gridLevel = gl;}

	///	returns the level
		const GridLevel& level() const {return m_gridLevel;}

	///	initializes the operator that may depend on the current solution
		virtual void init(const vector_type& u);

	///	initialize the operator
		virtual void init();

	///	initializes the operator and assembles the passed rhs vector
		void init_op_and_rhs(vector_type& b);

	///	compute d = J(u)*c (here, J(u) is a Matrix)
		virtual void apply(vector_type& d, const vector_type& c);

	///	Compute d := d - J(u)*c
		virtual void apply_sub(vector_type& d, const vector_type& c);

	///	Set Dirichlet values
		void set_dirichlet_values(vector_type& u);

	///	Destructor
		virtual ~AssembledLinearOperator() {};

	protected:
	// 	assembling procedure
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

	// 	DoF Distribution used
		GridLevel m_gridLevel;
};

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/// help function to assemble a linear operator
/**
 * This function initializes the operator, sets the vector b to the computed rhs
 * and sets the dirichlet post processes for the vector u.
 *
 * \param[out]	op		Operator
 * \param[out]	u		Solution
 * \param[out]	b		Rigth-Hand side vector
 *
 * \tparam	TAlgebra			algebra type
 */
template <typename TAlgebra>
void AssembleLinearOperatorRhsAndSolution
		(AssembledLinearOperator<TAlgebra>& op,
		 typename TAlgebra::vector_type& u,
		 typename TAlgebra::vector_type& b);

} // namespace ug

// include implementation
#include "assembled_linear_operator_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */
