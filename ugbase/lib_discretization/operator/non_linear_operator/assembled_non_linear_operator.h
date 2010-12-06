#ifndef __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__

#include "lib_algebra/lib_algebra.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
class AssembledOperator : public IOperator<	typename TAlgebra::vector_type,
											typename TAlgebra::vector_type>
{
public:
	// 	Type of algebra
		typedef TAlgebra algebra_type;

	//	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	//	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	//	Type of DoFDistribution
		typedef TDoFDistribution dof_distribution_type;

	public:
		AssembledOperator() :
			m_bInit(false), m_pAss(NULL), m_pDoFDistribution(NULL)
		{};

		AssembledOperator(IAssemble<dof_distribution_type, algebra_type>& ass) :
			m_bInit(false), m_pAss(&ass), m_pDoFDistribution(NULL)
		{};

		void set_discretization(IAssemble<TDoFDistribution, algebra_type>& ass) {m_pAss = &ass;}

		bool set_dof_distribution(const IDoFDistribution<TDoFDistribution>& dofDistr)
		{
			m_pDoFDistribution = &dofDistr;
			return true;
		}

		const IDoFDistribution<TDoFDistribution>* get_dof_distribution()
		{
			return m_pDoFDistribution;
		}

	//	Init
		virtual bool init()
		{
			if(m_pDoFDistribution == NULL)
			{
				UG_LOG("ERROR in AssembledOperator::init: DoF Distribution not set.\n");
				return false;
			}
			if(m_pAss == NULL)
			{
				UG_LOG("ERROR in AssembledOperator::prepare: Discretization not set.\n");
				return false;
			}

		//	remember that operator has been init
			m_bInit = true;

			return true;
		}

	//	Prepare functions
		virtual bool prepare(vector_type& dOut, vector_type& uIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledOperator::prepare: Operator not initialized.\n");
				return false;
			}

		// 	Set Dirichlet - Nodes to exact values
			if(m_pAss->assemble_solution(uIn, *m_pDoFDistribution) != IAssemble_OK)
				{UG_LOG("AssembledOperator::apply: Cannot set dirichlet values in solution.\n"); return false;}

			return true;
		}

	// 	Compute d = L(u)
		virtual bool apply(vector_type& dOut, const vector_type& uIn)
		{
			if(!m_bInit)
			{
				UG_LOG("ERROR in AssembledOperator::apply: Operator not initialized.\n");
				return false;
			}

			// reset vector
			if(!dOut.set(0.0))
				{UG_LOG("AssembledOperator::apply: Could not reset defect "
						"to zero before assembling. Aborting.\n"); return false;}

			// assemble
			if(m_pAss->assemble_defect(dOut, uIn, *m_pDoFDistribution) != IAssemble_OK)
				{UG_LOG("AssembledOperator::apply: Could not "
						"assemble defect. Aborting.\n"); return false;}

#ifdef UG_PARALLEL
			dOut.set_storage_type(PST_ADDITIVE);
#endif
			return true;
		}

		IAssemble<TDoFDistribution, algebra_type>* get_assemble()
		{
			return m_pAss;
		}

	protected:
		// init flag
		bool m_bInit;

		// assembling procedure
		IAssemble<dof_distribution_type, algebra_type>* m_pAss;

		// DoF Distribution used
		const IDoFDistribution<TDoFDistribution>* m_pDoFDistribution;
};

} // end namepace ug

#endif /*__H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__*/
