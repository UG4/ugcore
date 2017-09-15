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

#ifndef __H__UG__LIB_DISC__TIME_DISC__COMBINED_TIME_DISC__
#define __H__UG__LIB_DISC__TIME_DISC__COMBINED_TIME_DISC__

#include <cstddef>                                   // for size_t
#include <vector>                                    // for vector

#include "lib_disc/time_disc/time_disc_interface.h"  // for ITimeDiscretization


namespace ug {

/// \ingroup lib_disc_time_assemble
/// @{

/// combine several time discretizations into one
template <typename TAlgebra>
class CombinedTimeDiscretization
	: public IAssemble<TAlgebra>
{
	public:
		typedef TAlgebra algebra_type;
		typedef typename algebra_type::vector_type vector_type;
		typedef typename algebra_type::matrix_type matrix_type;

	public:
		void add_time_disc(SmartPtr<ITimeDiscretization<TAlgebra> > tDisc)
		{m_vTimeDisc.push_back(tDisc);}

	public:
		/**
		 *  @note This is a hack to ensure the time disc can be used with a GMG.
		 *        As the GMG calls IAssemble::ass_tuner(), and then uses
		 *        set_force_regular_grid() on it, this call must be handed down
		 *        to the individual time discs here.
		 */
		class CombinedAssTuner
			: public AssemblingTuner<TAlgebra>
		{
			public:
				void add_ass_tuner(SmartPtr<AssemblingTuner<TAlgebra> > assTuner)
				{m_vAssTuner.push_back(assTuner);}

				void set_force_regular_grid(bool bForce)
				{
					for (size_t i = 0; i < m_vAssTuner.size(); ++i)
						m_vAssTuner[i]->set_force_regular_grid(bForce);
				}

			protected:
				std::vector<SmartPtr<AssemblingTuner<TAlgebra>, FreeDelete> > m_vAssTuner;
		};

	// inherited from ITimeDiscretization
	public:
		/// @copydoc ITimeDiscretization<TAlgebra>::prepare_step()
		virtual void prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol, number dt);

		/// @copydoc ITimeDiscretization<TAlgebra>::prepare_step_elem()
		virtual void prepare_step_elem
		(
			SmartPtr<VectorTimeSeries<vector_type> > prevSol,
			number dt,
			const GridLevel& gl
		);

		/// @copydoc ITimeDiscretization<TAlgebra>::finish_step()
		virtual void finish_step(SmartPtr<VectorTimeSeries<vector_type> > currSol);

		/// @copydoc ITimeDiscretization<TAlgebra>::finish_step_elem()
		virtual void finish_step_elem
		(
			SmartPtr<VectorTimeSeries<vector_type> > currSol,
			const GridLevel& gl
		);

		/// @copydoc ITimeDiscretization<TAlgebra>::future_time()
		virtual number future_time() const;

		/// @copydoc ITimeDiscretization<TAlgebra>::num_prev_steps()
		virtual size_t num_prev_steps() const;

		/// @copydoc ITimeDiscretization<TAlgebra>::num_stages()
		virtual size_t num_stages() const;

		/// @copydoc ITimeDiscretization<TAlgebra>::set_stage()
		virtual void set_stage(size_t stage);

	// inherited from IAssemble
	public:
		/// @copydoc IAssemble::assemble_jacobian
		void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl);

		/// @copydoc IAssemble::assemble_defect
		void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl);

		/// @copydoc IAssemble::assemble_linear
		void assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl);

		/// @copydoc IAssemble::assemble_rhs(vector_type&, vector_type&, GridLevel&)
		void assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl);

		/// @copydoc IAssemble::assemble_rhs(vector_type&, GridLevel&)
		void assemble_rhs(vector_type& b, const GridLevel& gl);

		/// @copydoc IAssemble::adjust_solution
		void adjust_solution(vector_type& u, const GridLevel& gl);

		///\{
		/// @copydoc IAssemble::ass_tuner
		virtual SmartPtr<AssemblingTuner<TAlgebra> > ass_tuner();
		virtual ConstSmartPtr<AssemblingTuner<TAlgebra> > ass_tuner() const;
		///\}

		///	@copydoc IAssemble::num_constraints
		virtual size_t num_constraints() const;

		///	@copydoc IAssemble::constraint
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i);

	protected:
		std::vector<SmartPtr<ITimeDiscretization<TAlgebra> > > m_vTimeDisc;
};

} // end namespace ug


/// }

// include implementation
#include "combined_time_disc_impl.h"

#endif // __H__UG__LIB_DISC__TIME_DISC__COMBINED_TIME_DISC__
