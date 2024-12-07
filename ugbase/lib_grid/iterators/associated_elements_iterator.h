/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_associated_elements_iterator
#define __H__UG_associated_elements_iterator

#include <iterator>
#include <limits>
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug{

///	Iterator that allows to traverse associated elements of a given element
/** The type of the element whose associated elements shall be traversed
 * has to be specified through the template argument TElem. The type of the
 * associated elements has to be specified through the TAssocElem template
 * argument. Through VRefOrder one can choose whether the traversed associated
 * elements shall be sorted in the order in which they appear in the reference
 * element of given element (false by default).*/
template <class TElem, class TAssocElem, bool VSorted = false>
class AssocElemIter /* todo remove deprecated inheritance
	todo evaluate if simplification of type from TAssocElem* to TAssocElem is possible?
	: public std::iterator<std::input_iterator_tag, TAssocElem*>*/
{
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = TAssocElem*;
			using difference_type = std::ptrdiff_t;
			using pointer = TAssocElem**;
			using reference = TAssocElem*&;

			AssocElemIter(typename Grid::traits<TAssocElem>::callback cbConsiderElem =
							ConsiderAll()) :
				m_i(0),
				m_cbConsiderElem(cbConsiderElem)
			{}

			AssocElemIter(Grid& grid, TElem* elem,
						  typename Grid::traits<TAssocElem>::callback cbConsiderElem =
							ConsiderAll()) :
				m_i(0),
				m_cbConsiderElem(cbConsiderElem)
			{
				init(grid, elem);
			}

			void set_callback(typename Grid::traits<TAssocElem>::callback cbConsiderElem)
			{
				m_cbConsiderElem(cbConsiderElem);
			}

			void reinit(Grid& grid, TElem* elem)
			{
				m_i = 0;
				init(grid, elem);
			}

			void reinit(Grid& grid, TElem* elem,
						typename Grid::traits<TAssocElem>::callback cb)
			{
				m_i = 0;
				m_cbConsiderElem(cb);
				init(grid, elem);
			}

			AssocElemIter end() const
			{
				AssocElemIter e;
				e.m_i = std::numeric_limits<size_t>::max();
				return e;
			}

			bool valid() const	{return m_i < m_assElems.size();}
			bool invalid() const	{return m_i >= m_assElems.size();}


			AssocElemIter& operator ++()			{increment(); return *this;}
			AssocElemIter operator ++(int unused)	{AssocElemIter i = *this; increment(); return i;}

		///	returns true if both iterators are invalid or if both point to the same elemnt.
			bool operator ==(const AssocElemIter& iter) const {return equal(iter);}
		///	returns true if exactly one iterator is invalid or if the iterators point to different elements.
			bool operator !=(const AssocElemIter& iter) const {return !equal(iter);}

			TAssocElem* operator *()	{return dereference();}

		private:
			inline void init(Grid& grid, TElem* elem){
				if(VSorted)
					grid.associated_elements_sorted(m_assElems, elem);
				else
					grid.associated_elements(m_assElems, elem);
				
				if(valid() && (!m_cbConsiderElem(dereference())))
					increment();
			}

		///	returns true if both iterators are invalid or if both point to the same elemnt.
			inline bool equal(AssocElemIter const& other) const{
				if(valid()){
					if(other.valid()){
						if(dereference() == *other)
							return true;
					}
				}
				else{
					if(!other.valid())
						return true;
				}
				return false;
			}

		///	returns next iterator
			void increment(){
				while(1){
					++m_i;
					if(valid()){
						if(m_cbConsiderElem(dereference()))
							break;
					}
					else
						break;
				}
			}

		///	dereference
			inline TAssocElem* dereference() const{
				return m_assElems[m_i];
			}

		private:
			size_t m_i;
			typename Grid::traits<TAssocElem>::secure_container	m_assElems;
			typename Grid::traits<TAssocElem>::callback			m_cbConsiderElem;
	};

}//	end of namespace

#endif	//__H__UG_associated_elements_iterator
