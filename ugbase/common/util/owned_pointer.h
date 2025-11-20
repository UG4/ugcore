/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__owned_pointer__
#define __H__UG__owned_pointer__

namespace ug
{

/// \addtogroup ugbase_common_util
/// \{

///	Holds and automatically deletes a pointer, similar to std::auto_ptr. USE WITH CARE!
/**	WARNING: USE WITH CARE!
 * Class shows uncommon behavior when copied.
 *
 * The behavior of this class is very similar to std::auto_ptr and somewhat
 * similar to boost::scoped_ptr. It takes a pointer in its constructor and
 * automatically deletes the associated object in its destructor or in a call to reset.
 * The most notable property of OwnedPtr is that an OwnedPtr 'owns' the associated
 * object. That means that two OwnedPtrs should never point to the same instance.
 * The copy-constructor thus shows the somewhat irritating behavior of invalidating
 * the original copy and transferring the ownership of the associated object to
 * the new copy of OwnedPtr. This behavior again is the same in std::auto_ptr,
 * however, in contrary to std::auto_ptr, instances of OwnedPtr can be used
 * in standard containers like std::vector. This is explicitly not advised, since
 * many algorithms working with iterators or standard containers will not work
 * with OwnedPtrs.
 * Storing OwnedPtrs in lib_grids attachments may sometimes be useful. Therefore
 * this slightly dangerous class was introduced.
 */
template <typename T>
class OwnedPtr
{
	public:
		using TPtr = T*;
		using TRef = T&;

		OwnedPtr(TPtr p = 0) : m_p(p)	{}

	///	Transfers ownership of the associated object from ap to this.
	/**	Note that ap looses ownership of the associated object.
	 * ap thus is not const at all, despite beeing declared const.*/
		OwnedPtr(const OwnedPtr& op)
		{
			m_p = op.get();
			op.invalidate();
		}

	///	Transfers ownership of the associated object from ap to this.
	/**	Note that ap looses ownership of the associated object.
	 * ap thus is not const at all, despite beeing declared const.*/
		OwnedPtr& operator = (const OwnedPtr& op)
		{
			reset(op.get());
			op.invalidate();
			return *this;
		}

		void reset(TPtr p = 0) {if(m_p) delete m_p; m_p = p;}

		TPtr operator -> () const {return m_p;}
		TRef operator * () const {return *m_p;}

		TPtr& get() const {return m_p;}

		operator bool ()	const {return m_p != 0;}

	private:
	///	This method is a hack, to allow that only one instance owns the pointer
		void invalidate() const				{m_p = 0;}

	private:
		mutable TPtr m_p;
};

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif
