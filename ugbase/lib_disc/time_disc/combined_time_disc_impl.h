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

#include "common/error.h"  // UG_COND_THROW

namespace ug {


template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::
prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol, number dt)
{
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		m_vTimeDisc[i]->prepare_step(prevSol, dt);
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::prepare_step_elem
(
	SmartPtr<VectorTimeSeries<vector_type> > prevSol,
	number dt,
	const GridLevel& gl
)
{
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		m_vTimeDisc[i]->prepare_step_elem(prevSol, dt, gl);
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::
finish_step(SmartPtr<VectorTimeSeries<vector_type> > currSol)
{
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		m_vTimeDisc[i]->finish_step(currSol);
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::finish_step_elem
(
	SmartPtr<VectorTimeSeries<vector_type> > currSol,
	const GridLevel& gl
)
{
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		m_vTimeDisc[i]->finish_step_elem(currSol, gl);
}

template <typename TAlgebra>
number CombinedTimeDiscretization<TAlgebra>::future_time() const
{
	if (m_vTimeDisc.size())
		return m_vTimeDisc[0]->future_time();
	UG_THROW("At least one time disc must be added to CombinedTimeDiscretization.")
}

template <typename TAlgebra>
size_t CombinedTimeDiscretization<TAlgebra>::num_prev_steps() const
{
	size_t nSteps = 0;
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		nSteps = std::max(nSteps, m_vTimeDisc[i]->num_prev_steps());

	return nSteps;
}

template <typename TAlgebra>
size_t CombinedTimeDiscretization<TAlgebra>::num_stages() const
{
	size_t nStages = 0;
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		nStages = std::max(nStages, m_vTimeDisc[i]->num_stages());

	return nStages;
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::set_stage(size_t stage)
{
	for (size_t i = 0; i < m_vTimeDisc.size(); ++i)
		m_vTimeDisc[i]->set_stage(stage);
}


template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::assemble_jacobian
(
	matrix_type& J,
	const vector_type& u,
	const GridLevel& gl
)
{
	UG_COND_THROW(!m_vTimeDisc.size(),
		"At least one time disc must be added to CombinedTimeDiscretization.")

	m_vTimeDisc[0]->assemble_jacobian(J, u, gl);
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
	{
		// avoid clearing of the matrix before assembling
		m_vTimeDisc[i]->domain_disc()->ass_tuner()->disable_clear_on_resize();

		m_vTimeDisc[i]->assemble_jacobian(J, u, gl);
	}
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::assemble_defect
(
	vector_type& d,
	const vector_type& u,
	const GridLevel& gl
)
{
	UG_COND_THROW(!m_vTimeDisc.size(),
		"At least one time disc must be added to CombinedTimeDiscretization.")

	m_vTimeDisc[0]->assemble_defect(d, u, gl);
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
	{
		// avoid clearing of the matrix before assembling
		m_vTimeDisc[i]->domain_disc()->ass_tuner()->disable_clear_on_resize();

		m_vTimeDisc[i]->assemble_defect(d, u, gl);
	}
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::assemble_linear
(
	matrix_type& A,
	vector_type& b,
	const GridLevel& gl
)
{
	UG_COND_THROW(!m_vTimeDisc.size(),
		"At least one time disc must be added to CombinedTimeDiscretization.")

	m_vTimeDisc[0]->assemble_linear(A, b, gl);
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
	{
		// avoid clearing of the matrix before assembling
		m_vTimeDisc[i]->domain_disc()->ass_tuner()->disable_clear_on_resize();

		m_vTimeDisc[i]->assemble_linear(A, b, gl);
	}
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::assemble_rhs
(
	vector_type& b,
	const vector_type& u,
	const GridLevel& gl
)
{
	UG_COND_THROW(!m_vTimeDisc.size(),
		"At least one time disc must be added to CombinedTimeDiscretization.")

	m_vTimeDisc[0]->assemble_rhs(b, u, gl);
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
	{
		// avoid clearing of the matrix before assembling
		m_vTimeDisc[i]->domain_disc()->ass_tuner()->disable_clear_on_resize();

		m_vTimeDisc[i]->assemble_rhs(b, u, gl);
	}
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::assemble_rhs(vector_type& b, const GridLevel& gl)
{
	UG_COND_THROW(!m_vTimeDisc.size(),
		"At least one time disc must be added to CombinedTimeDiscretization.")

	m_vTimeDisc[0]->assemble_rhs(b, gl);
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
	{
		// avoid clearing of the matrix before assembling
		m_vTimeDisc[i]->domain_disc()->ass_tuner()->disable_clear_on_resize();

		m_vTimeDisc[i]->assemble_rhs(b, gl);
	}
}

template <typename TAlgebra>
void CombinedTimeDiscretization<TAlgebra>::adjust_solution(vector_type& u, const GridLevel& gl)
{
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
		m_vTimeDisc[i]->adjust_solution(u, gl);
}

template <typename TAlgebra>
SmartPtr<AssemblingTuner<TAlgebra> > CombinedTimeDiscretization<TAlgebra>::ass_tuner()
{
	SmartPtr<CombinedAssTuner> sp = make_sp(new CombinedAssTuner());
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
		sp->add_ass_tuner(((IAssemble<TAlgebra>*)m_vTimeDisc[i].get())->ass_tuner());

	return sp;
}

template <typename TAlgebra>
ConstSmartPtr<AssemblingTuner<TAlgebra> > CombinedTimeDiscretization<TAlgebra>::ass_tuner() const
{
	UG_THROW("Unique const AssemblingTuner cannot be provided by CombinedTimeDiscretization.")
}

template <typename TAlgebra>
size_t CombinedTimeDiscretization<TAlgebra>::num_constraints() const
{
	size_t n = 0;
	for (size_t i = 1; i < m_vTimeDisc.size(); ++i)
		n += m_vTimeDisc[i]->num_constraints();

	return n;
}

template <typename TAlgebra>
SmartPtr<IConstraint<TAlgebra> > CombinedTimeDiscretization<TAlgebra>::constraint(size_t i)
{
	UG_COND_THROW(i >= num_constraints(), "Requested constraint " << i << ", but only "
		<< num_constraints() << " constraints available.");

	size_t n = 0;
	size_t k = 0;

	while (n += m_vTimeDisc[k]->num_constraints() <= i)
		++k;

	return m_vTimeDisc[k]->constraint(i - n);
}

} // end namespace ug


