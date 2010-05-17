#ifndef __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"

namespace ug{

template <typename TDiscreteFunction>
class AssembledDiscreteOperator : public IDiscreteOperator<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// export types:

		// domain function type
		typedef TDiscreteFunction domain_function_type;

		// codomain function type
		typedef TDiscreteFunction codomain_function_type;

		// type of algebra
		typedef typename TDiscreteFunction::algebra_type algebra_type;

	public:
		AssembledDiscreteOperator(IAssemble<algebra_type, domain_function_type>& ass) :
			m_ass(ass)
		{};

		virtual bool init()
		{
			return true;
		}

		virtual bool prepare(domain_function_type& u, codomain_function_type& d)
		{
			return true;
		}

		// compute f = L*u (here, L is a Matrix)
		virtual bool apply(domain_function_type& u, codomain_function_type& d)
		{
			typename codomain_function_type::vector_type& d_vec = d.get_vector();

			// reset vector
			if(d_vec.set(0.0) != true)
			{
				UG_LOG("AssembledDiscreteOperator::apply: Could not reset defect to zero before assembling. Aborting.\n");
				return false;
			}

			// assemble
			if(m_ass.assemble_defect(d_vec, u) != IAssemble_OK)
			{
				UG_LOG("AssembledDiscreteOperator::apply: Could not assemble defect. Aborting.\n");
				return false;
			}

			return true;
		}

		IAssemble<algebra_type, domain_function_type>* get_assemble()
		{
			return &m_ass;
		}

	protected:
		// assembling procedure
		IAssemble<algebra_type, domain_function_type>& m_ass;
};

} // end namepace ug

#endif /*__H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__*/
