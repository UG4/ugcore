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
// used by plugin: SuperLU
#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__CONSTRAINED_LINEAR_ITERATOR__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__CONSTRAINED_LINEAR_ITERATOR__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/domain_disc_interface.h"
#include "preconditioner.h"

#include <boost/core/enable_if.hpp>
#include <boost/type_traits/is_base_of.hpp>

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
 *
 * @note This class may have lost its use:
 *       Constraints are now usually implemented in such a way that
 *       that both defect and correction is always zero in constrained
 *       DoFs.
 */
template <typename TDomain, typename TAlgebra, typename TLinIt, typename = void>
class ConstrainedLinearIterator : public TLinIt
{
	public:
		using gf_type = GridFunction<TDomain, TAlgebra>;
		using vector_type = typename TAlgebra::vector_type;
		using base_type = TLinIt;

		using TLinIt::init;

	public:
		///	constructor
		ConstrainedLinearIterator(SmartPtr<IDomainDiscretization<TAlgebra> > domDisc)
		: base_type(), m_spDomDisc(domDisc), m_time(0.0)
		{}

		/// clone constructor
		ConstrainedLinearIterator(const ConstrainedLinearIterator<TDomain, TAlgebra, TLinIt> &parent)
		: base_type(parent), m_spDomDisc(parent.m_spDomDisc), m_time(parent.m_time)
		{}

		///	clone
		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ConstrainedLinearIterator<TDomain, TAlgebra, TLinIt>(*this));
		}

		///	destructor
		virtual ~ConstrainedLinearIterator() = default;

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
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_DIRICHLET, m_time);

				// then other interpolation
				for (size_t i = 0; i < nConstr; i++)
					if (m_spDomDisc->constraint(i)->type() & CT_CONSTRAINTS)
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_CONSTRAINTS, m_time);

				// then hanging nodes
				for (size_t i = 0; i < nConstr; i++)
					if (m_spDomDisc->constraint(i)->type() & CT_HANGING)
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_HANGING, m_time);

				// and Dirichlet again (Dirichlet nodes might also be hanging)
				for (size_t i = 0; i < nConstr; i++)
					if (m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
						m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_DIRICHLET, m_time);
			}
			else UG_THROW("Vector is not a grid function. This is not supported.")

			return true;
		}

		virtual bool apply_update_defect(vector_type& c, vector_type& d)
		{
			apply_update_defect_impl<TLinIt> impl(*this);
			return impl(c, d);
		}

	private:
		template <typename S, typename = void>//, typename otherDummy = void>
		struct apply_update_defect_impl
		{
			ConstrainedLinearIterator& cli;

			explicit apply_update_defect_impl(ConstrainedLinearIterator& _cli)
			: cli(_cli) {}

			bool operator () (vector_type& c, vector_type& d)
			{
				// FIXME: This implementation is not correct.
				// As the defect is calculated from possibly erroneous corrections, there is no telling
				// how to adjust the defect. Except for setting it to zero everywhere, which, by the way,
				// is NOT what adjust_defect does, in general. We would need something like
				// adjust_defect_linearSolver() here to distinguish between the non-linear defect and
				// the residual in the linear solver.
				// We CANNOT use our own apply() here and then update the defect ourselves as the
				// ILinearIterator does not have an update method nor any way of applying the underlying
				// linear operator.

				// perform normal step
				if (!cli.base_type::apply_update_defect(c, d)) return false;

				// apply constraints on correction and defect
				if (!cli.m_spDomDisc.valid())
					UG_THROW("Domain discretization passed to ConstrainedLinearIterator is not valid.\n"
							 "Cannot apply any constraints.")

				gf_type* gf = dynamic_cast<gf_type*>(&c);
				if (gf && gf->dof_distribution().valid())
				{
					size_t nConstr = cli.m_spDomDisc->num_constraints();

				// corrections
					// first Dirichlet (hanging nodes might be constrained by Dirichlet nodes)
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
							cli.m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_DIRICHLET, cli.m_time);

					// then other interpolation
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_CONSTRAINTS)
							cli.m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_CONSTRAINTS, cli.m_time);

					// then hanging nodes
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_HANGING)
							cli.m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_HANGING, cli.m_time);

					// and Dirichlet again (Dirichlet nodes might also be hanging)
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
							cli.m_spDomDisc->constraint(i)->adjust_correction(c, gf->dof_distribution(), CT_DIRICHLET, cli.m_time);

				// residuals
					// hanging nodes first
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_HANGING)
							cli.m_spDomDisc->constraint(i)->adjust_linear_residual(d, c, gf->dof_distribution(), CT_HANGING, cli.m_time);

					// then other constraints
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_CONSTRAINTS)
							cli.m_spDomDisc->constraint(i)->adjust_linear_residual(d, c, gf->dof_distribution(), CT_CONSTRAINTS, cli.m_time);

					// and Dirichlet last
					for (size_t i = 0; i < nConstr; i++)
						if (cli.m_spDomDisc->constraint(i)->type() & CT_DIRICHLET)
							cli.m_spDomDisc->constraint(i)->adjust_linear_residual(d, c, gf->dof_distribution(), CT_DIRICHLET, cli.m_time);

				}
				else UG_THROW("Vector is not a grid function. This is not supported.")

				return true;
			}
		};


		// special implementation for IPreconditioner
		template <typename S>
		struct apply_update_defect_impl<S, typename boost::enable_if<boost::is_base_of<IPreconditioner<TAlgebra>, S> >::type>
		{
			ConstrainedLinearIterator& cli;

			explicit apply_update_defect_impl(ConstrainedLinearIterator& _cli)
			: cli(_cli) {}

			bool operator () (vector_type& c, vector_type& d)
			{
				// apply precond and constraints
				if (!cli.apply(c, d)) return false;

				// calculate new defect from preconditioner defect operator
				cli.m_spDefectOperator->apply_sub(d, c);

				return true;
			}
		};

		friend struct apply_update_defect_impl<TLinIt>;


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
		SmartPtr<ILinearIterator<vector_type> > m_spLinIt;
		SmartPtr<IDomainDiscretization<TAlgebra> > m_spDomDisc;
		number m_time;
};




} // end namespace ug


#endif