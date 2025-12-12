/*
 * Copyright (c) 2013-2017:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PPREPROCESS__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PPREPROCESS__

#include <vector>

#include "common/util/smart_pointer.h"

namespace ug {

///////////////////////////////////////////////////////////////////////////////
// Pre/Post Process for solvers
///////////////////////////////////////////////////////////////////////////////

/// interface for pre- and postprocess functions
/**
 * This is a base (interface) class for pre- and postprocess operations
 * for solvers. A typical application of these operations is the projections
 * eliminating the kernel parts etc.
 *
 * This class defines a generic virtual function for the vectors. If the
 * function requires the geometry, the RTTI should be used to cast the vector
 * to the grid gunction type.
 *
 * \tparam	TVector	the vector type in the algebra
 */
template <typename TVector>
class IPProcessVector
{
public:
	using vector_type = TVector; ///< the vector type
	
///	user-defined pre- or post-process function
	virtual void apply (vector_type& v) = 0;
	
///	virtual destructor
	virtual ~IPProcessVector () = default;
};

///	a chain of pre- or postprocess operations
/**
 * A class for storing a chain of pre- or postprocess operations
 *
 * \tparam	TVector	the vector type in the algebra
 */
template <typename TVector>
class PProcessChain
{
public:
	using p_process_type = IPProcessVector<TVector>;
	using vector_type = TVector; ///< the vector type
	
public:

///	performs all the operations
	void apply (vector_type& v)
	{
		for (size_t i = 0; i < m_pp_chain.size (); i++)
		{
			SmartPtr<p_process_type> op = m_pp_chain[i];
			if (op.valid ()) op->apply (v);
		}
	}

///	adds an operation at the end of the chain
	void add (SmartPtr<p_process_type> p)
	{
		m_pp_chain.push_back (p);
	}
	
///	returns the number of the operations
	[[nodiscard]] size_t size () const {return m_pp_chain.size ();}
	
///	returns the operation #i
	SmartPtr<p_process_type> operator [] (size_t i) {return m_pp_chain [i];}

///	removes an operation from the chain
	void remove (SmartPtr<p_process_type> p)
	{
		size_t i = 0;
		while (i < m_pp_chain.size ())
		{
			if (m_pp_chain[i].get () == p.get ())
				m_pp_chain.erase (m_pp_chain.begin () + i);
			else
				i++;
		}
	}

private:

	std::vector<SmartPtr<p_process_type> > m_pp_chain;
};

} // end namespace ug
#endif