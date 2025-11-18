/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__generic_grid_object_iterator__
#define __H__UG__generic_grid_object_iterator__

namespace ug
{

////////////////////////////////////////////////////////////////////////////////////////////////
//	GenericGridObjectIterator
///	Use this class as a tool to create iterators to your own geometric objects.
template <class TValue, class TBaseIterator>
class GenericGridObjectIterator : public TBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		using value_type = TValue;

	public:
		GenericGridObjectIterator()	= default;

		GenericGridObjectIterator(const GenericGridObjectIterator& iter) :
			TBaseIterator(iter)	{}

	///	note that the * operator is read only.
		inline TValue operator* () const	{return static_cast<TValue>(TBaseIterator::operator*());}

	protected:
		explicit GenericGridObjectIterator(const TBaseIterator& iter) :
			TBaseIterator(iter)	{}
};

////////////////////////////////////////////////////////////////////////////////////////////////
//	ConstGenericGridObjectIterator
///	Use this class as a tool to create const_iterators to your own geometric objects.
template <class TValue, class TBaseIterator, class TConstBaseIterator>
class ConstGenericGridObjectIterator : public TConstBaseIterator
{
	friend class Grid;
	template <class TIterDest, class TIterSrc> friend TIterDest iterator_cast(const TIterSrc& iter);

	public:
		using value_type = TValue;

	public:
		ConstGenericGridObjectIterator() = default;

		ConstGenericGridObjectIterator(const ConstGenericGridObjectIterator& iter) :
			TConstBaseIterator(iter)	{}

		ConstGenericGridObjectIterator(const GenericGridObjectIterator<TValue, TBaseIterator>& iter) :
			TConstBaseIterator(iter) {}

	///	note that the * operator is read only.
		inline TValue operator* () const {return static_cast<TValue>(TConstBaseIterator::operator*());}

	protected:
		explicit ConstGenericGridObjectIterator(const TBaseIterator& iter) :
			TConstBaseIterator(iter) {}

		explicit ConstGenericGridObjectIterator(const TConstBaseIterator& iter) :
			TConstBaseIterator(iter) {}
};

////////////////////////////////////////////////////////////////////////
//	iterator_cast
///	You should avoid casting whenever possible!
template <class TIterDest, class TIterSrc>
inline TIterDest
iterator_cast(const TIterSrc& iter)
{
	return TIterDest(iter);
}

}//	end of namespace

#endif
