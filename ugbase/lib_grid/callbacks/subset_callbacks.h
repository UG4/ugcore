// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_subset_callbacks
#define __H__UG_subset_callbacks

#include "lib_grid/tools/subset_handler_interface.h"

namespace ug{

/** \ingroup lib_grid_element_callbacks
 * \{ */

///	Element callback that returns true, if an element is contained in a subset
class IsInSubset
{
	public:
		IsInSubset(const ISubsetHandler& sh, int subsetIndex) :
			m_sh(sh),
			m_si(subsetIndex)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}
		bool operator() (Volume* v)		{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return m_sh.get_subset_index(e) == m_si;}

	private:
		const ISubsetHandler& m_sh;
		int m_si;
};

///	Element callback that returns true, if an element is not contained in a subset
class IsNotInSubset
{
	public:
		IsNotInSubset(const ISubsetHandler& sh, int subsetIndex) :
			m_sh(sh),
			m_si(subsetIndex)	{}

		bool operator() (Vertex* v)	{return callback(v);}
		bool operator() (Edge* e)	{return callback(e);}
		bool operator() (Face* f)		{return callback(f);}
		bool operator() (Volume* v)		{return callback(v);}

	private:
		template <class TElem>
		bool callback(TElem* e)			{return m_sh.get_subset_index(e) != m_si;}

	private:
		const ISubsetHandler& m_sh;
		int m_si;
};

/** \} */ //lib_grid_element_callbacks

}//	end of namespace

#endif	//__H__UG_subset_callbacks
