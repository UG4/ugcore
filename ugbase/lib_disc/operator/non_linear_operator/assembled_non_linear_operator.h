/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_disc/assemble_interface.h"

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
			: m_spAss(NULL), m_gridLevel() {};

	///	constructor
		AssembledOperator(SmartPtr<IAssemble<TAlgebra> > ass)
			: m_spAss(ass), m_gridLevel(){};

	///	constructor
		AssembledOperator(SmartPtr<IAssemble<TAlgebra> > ass, const GridLevel& gl)
			: m_spAss(ass), m_gridLevel(gl) {};

	///	sets discretization for assembling
		void set_discretization(SmartPtr<IAssemble<TAlgebra> > ass) {m_spAss = ass;}

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
		SmartPtr<IAssemble<TAlgebra> > discretization() {return m_spAss;}

	protected:
	///	assembling procedure
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

	///	used grid level
		GridLevel m_gridLevel;
};

} // end namepace ug

// include implementation
#include "assembled_non_linear_operator_impl.h"

#endif /*__H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__*/
