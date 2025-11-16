/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Matrix Operator type
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	///	Base type
		using base_type = IPreconditioner<TAlgebra>;

	protected:
		SmartPtr<ILinearOperatorInverse<vector_type>  >  m_opInv;

	public:
		virtual SmartPtr<ILinearIterator<vector_type, vector_type> > clone()
		{
			UG_ASSERT(0, "not implemented since ILinearOperatorInverse::clone not implemented");
			return nullptr;
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
