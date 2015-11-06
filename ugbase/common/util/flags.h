/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__flags__
#define __H__UG__flags__

namespace ug{

///	Helps maintaining, activating and deactivating a set of flags from an enum.
/**	Given an enum
 * \code
 * enum SomeEnum{
 * 	E1 = 0,
 * 	E2 = 1,
 * 	E3 = 1 << 1,
 * 	E4 = 1 << 2
 * };
 *
 * \endcode
 * One can use the Flag class as follows:
 *
 * \code
 * typedef Flag<SomeEnum> SomeFlag;
 * ...
 * SomeFlag f1(E2), f2(E4);
 * SomeFlag f3(f1 | f2);
 * if(f3.contains(E2))
 * 	f3.remove(E2);
 * ...
 * \endcode
 *
 */
template <class TEnum, class TStorageType = unsigned int, TStorageType defaultValue = 0>
class Flag{
	public:
		Flag()					: m_value(defaultValue)	{}
		Flag(TStorageType flag)	: m_value(flag)			{}
		Flag(const Flag& flag)	: m_value(flag.m_value)	{}

		bool contains(TStorageType flag) const	{return (m_value & flag) == flag;}
		bool contains(const Flag& flag) const	{return (m_value & flag.m_value) == flag.m_value;}

		bool partially_contains(TStorageType flag) const	{return (m_value & flag) != 0;}
		bool partially_contains(const Flag& flag) const		{return (m_value & flag.m_value) != 0;}

		Flag& set(TStorageType flag)			{m_value = flag; return *this;}
		Flag& add(TStorageType flag)			{m_value |= flag; return *this;}
		Flag& remove(TStorageType flag)			{m_value &= (~flag); return *this;}

		Flag operator& (const Flag& flag) const	{return Flag(m_value & flag.m_value);}
		Flag operator&= (const Flag& flag)		{m_value &= flag.m_value; return *this;}
		Flag operator| (const Flag& flag) const	{return Flag(m_value | flag.m_value);}
		Flag operator|= (const Flag& flag)		{m_value |= flag.m_value; return *this;}
		Flag operator= (const Flag& flag)		{m_value = flag.m_value; return *this;}
		Flag operator= (TStorageType val)		{m_value = val; return *this;}

		TStorageType operator()() const			{return m_value;}
		TStorageType get() const				{return m_value;}

		bool operator== (const Flag& flag) const	{return m_value == flag.m_value;}
		bool operator== (TStorageType val) const	{return m_value == val;}

		bool operator!= (const Flag& flag) const	{return m_value != flag.m_value;}
		bool operator!= (TStorageType val) const	{return m_value != val;}

	private:
		TStorageType	m_value;
};

}// end of namespace

#endif
