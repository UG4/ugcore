// created by Sebastian Reiter
// s.b.reiter@gmail.com

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
class AssocElemIter : public std::iterator<std::input_iterator_tag, TAssocElem*>
{
		public:
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
