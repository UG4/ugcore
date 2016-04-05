/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__CONSTRAINED_LINEAR_ITERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__CONSTRAINED_LINEAR_ITERATOR__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"


namespace ug {


/**
 * This class is a template for constraint-respecting versions of
 * any linear iterator derived from ILinearIterator. This includes
 * linear solvers as well as preconditioners.
 * The derived ILinearIterator class must have a constructor
 * without arguments.
 * The handled constraints must implement the two methods:
 *     adjust_correction(),
 *     adjust_defect().
 */
template <typename TDomain, typename TAlgebra, typename TLinIt>
class ConstrainedLinearIterator : public TLinIt
{
	public:
		typedef GridFunction<TDomain, TAlgebra> gf_type;
		typedef typename TAlgebra::vector_type vector_type;
		typedef TLinIt base_type;

		using TLinIt::init;

	public:
		///	constructor
		ConstrainedLinearIterator(SmartPtr<IDomainDiscretization<TAlgebra> > domDisc)
		: base_type(), m_spDomDisc(domDisc), m_time(0.0)
		{}

		/// clone constructor
		ConstrainedLinearIterator(const ConstrainedLinearIterator<TDomain, TAlgebra, TLinIt> &parent)
		: base_type(), m_spDomDisc(parent.m_spDomDisc), m_time(parent.m_time)
		{}

		///	clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ConstrainedLinearIterator<TDomain, TAlgebra, TLinIt>(*this));
		}

		///	destructor
		virtual ~ConstrainedLinearIterator(){}

	protected:
		/// @copydoc ILinearIterator::name()
		virtual const char* name() const
		{
			return base_type::name();
		}

		virtual bool init(SmartPtr<ILinearOperator<vector_type,vector_type> > J, const vector_type& u)
		{
			return base_type::init(J, u);
		}

		virtual bool init(SmartPtr<ILinearOperator<vector_type,vector_type> > L)
		{
			return base_type::init(L);
		}

		virtual bool apply(vector_type& c, const vector_type& d)
		{
			// perform normal step
			if (!base_type::apply(c, d)) return false;

			// apply constraints on correction
			if (!m_spDomDisc.valid())
				UG_THROW("Domain discretization passed to ConstrainedLinearIterator is not valid.\n"
						 "Cannot apply any constraints.")

			gf_type* gf = dynamic_cast<gf_type*>(&c);
			if (gf && gf->dof_distribution().valid())
			{
				size_t nConstr = m_spDomDisc->num_constraints();

				// first Dirichlet (hanging nodes might be constrained by Dirichlet nodes)
				for (size_t i = 0; i < nConstr; i++)
					if (m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), m_time);

				// then hanging nodes and other interpolation
				for (size_t i = 0; i < nConstr; i++)
					if (m_spDomDisc->constraint(i)->type() & CT_CONSTRAINTS)
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), m_time);

				// and Dirichlet again (Dirichlet nodes might also be hanging)
				for (size_t i = 0; i < nConstr; i++)
					if (m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), m_time);
			}
			else UG_THROW("Vector is not a grid function. This is not supported.")

			return true;
		}

		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{
			// perform normal step
			if (!base_type::apply_update_defect(c, d)) return false;

			// apply constraints on correction and defect
			if (!m_spDomDisc.valid())
				UG_THROW("Domain discretization passed to ConstrainedLinearIterator is not valid.\n"
						 "Cannot apply any constraints.")

			gf_type* gf = dynamic_cast<gf_type*>(&c);
			if (gf && gf->dof_distribution().valid())
			{
				size_t nConstr = m_spDomDisc->num_constraints();

				// first Dirichlet (hanging nodes might be constrained by Dirichlet nodes)
				for (size_t i = 0; i < nConstr; i++)
				{
					if (m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
					{
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), m_time);
						m_spDomDisc->constraint(i)->adjust_defect(d, c, gf->dof_distribution(), m_time);
					}
				}

				// then hanging nodes and other interpolation
				for (size_t i = 0; i < nConstr; i++)
				{
					if (m_spDomDisc->constraint(i)->type() & CT_CONSTRAINTS)
					{
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), m_time);
						m_spDomDisc->constraint(i)->adjust_defect(d, c, gf->dof_distribution(), m_time);
					}
				}

				// and Dirichlet again (Dirichlet nodes might also be hanging)
				for (size_t i = 0; i < nConstr; i++)
				{
					if (m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
					{
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), m_time);
						m_spDomDisc->constraint(i)->adjust_defect(d, c, gf->dof_distribution(), m_time);
					}
				}
			}
			else UG_THROW("Vector is not a grid function. This is not supported.")

			return true;
		}

	public:
		/// @copydoc ILinearIterator::supports_parallel()
		virtual bool supports_parallel() const
		{
			return base_type::supports_parallel();
		}

		/// setter for time
		void set_time(number time)
		{
			m_time = time;
		};

	protected:
		///	storage for factorization
		SmartPtr<ILinearIterator<vector_type> > m_spLinIt;
		SmartPtr<IDomainDiscretization<TAlgebra> > m_spDomDisc;
		number m_time;
};


} // end namespace ug


#endif // __H__LIB_ALGEBRA__OPERATOR__INTERFACE__CONSTRAINED_LINEAR_ITERATOR__
