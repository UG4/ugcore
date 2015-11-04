
#ifndef __H__UG__LIB_DISC__OPERATOR__ITERATOR_OPERATOR_INVERSE__
#define __H__UG__LIB_DISC__OPERATOR__ITERATOR_OPERATOR_INVERSE__

#include <string>

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioner.h"

#include "common/log.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

/**
 * a LinearIterator which can uses ILinearOperatorInverse to perform B^{-1}
 * this is for the case that some class needs a preconditioner, but we'd like to use a linear solver
 * example: 4x AMG as preconditioner
 * \code
 * linSolver = LinearSolver()
 * linSolver:set_preconditioner(amg)
 * linSolver:set_convergence_check(ConvCheck(4, 0, 0, false) )
 * oii = OperatorInverseIterator(linSolver)
 * someObject:set_preconditioner(oii)
 * \endcode
 */
template <typename TAlgebra>
class OperatorInverseIterator : public ILinearIterator<typename TAlgebra::vector_type>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Matrix Operator type
		typedef typename IPreconditioner<TAlgebra>::matrix_operator_type matrix_operator_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		SmartPtr<ILinearOperatorInverse<vector_type>  >  m_opInv;

	public:
		virtual SmartPtr<ILinearIterator<vector_type, vector_type> > clone()
		{
			UG_ASSERT(0, "not implemented since ILinearOperatorInverse::clone not implemented");
			return SPNULL;
		}
	///	default constructor
		OperatorInverseIterator(SmartPtr<ILinearOperatorInverse<vector_type>  > opInv) : m_opInv(opInv)
		{
			m_name = std::string("OperatorInverseIterator(") + std::string(m_opInv->name()) + std::string(")");
		}
		~OperatorInverseIterator()
		{

		}

		std::string m_name;

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			return m_opInv->supports_parallel();
		}

		virtual const char* name() const
		{
			return m_name.c_str();
		}

		virtual bool init(SmartPtr<ILinearOperator<vector_type> > L)
		{
			if(!m_opInv->init(L))
			{
				UG_LOG("ERROR in '" << name() << "::init'.\n");
				return false;
			}
			return true;
		}

		virtual bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
		{
			if(!m_opInv->init(J, u))
			{
				UG_LOG("ERROR in '" << name() << "::init'.\n");
				return false;
			}
			return true;
		}

		virtual bool apply(vector_type& c, const vector_type& d)
		{
			if(m_opInv->apply(c, d))
			{
				//UG_LOG("ERROR in '" << name() << "::apply'\n");
				return false;
			}
			return true;
		}

		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{
			if(m_opInv->apply_return_defect(c, d))
			{
				//UG_LOG("ERROR in '" << name() << "::apply_update_defect'\n");
				return false;
			}
			return true;
		}

};


} // end namespace ug

#endif // __H__UG__LIB_DISC__OPERATOR__ITERATOR_OPERATOR_INVERSE__
