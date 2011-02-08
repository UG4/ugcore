// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SUBSET_HANDLER_CALLBACKS__
#define __H__LIB_GRID__SUBSET_HANDLER_CALLBACKS__

#include "lib_grid/tools/subset_handler_interface.h"
#include "callback_definitions.h"

namespace ug
{

/**	A wrapper that returns whether an object is in a given subset. Instances can
 * be used as callbacks CB_ConsiderVertex, ..., CB_ConsiderVolume.
 */
class IsInSubset
{
	public:
		IsInSubset(const ISubsetHandler& sh, int subsetIndex) :
			m_sh(sh),
			m_si(subsetIndex)	{}

		bool operator() (VertexBase* v)	{return m_sh.get_subset_index(v) == m_si;}
		bool operator() (EdgeBase* e)	{return m_sh.get_subset_index(e) == m_si;}
		bool operator() (Face* f)		{return m_sh.get_subset_index(f) == m_si;}
		bool operator() (Volume* v)		{return m_sh.get_subset_index(v) == m_si;}

	private:
		const ISubsetHandler& m_sh;
		int m_si;
};

/**	A wrapper that returns whether an object is not in a given subset. Instances
 * can be used as callbacks CB_ConsiderVertex, ..., CB_ConsiderVolume.
 */
class IsNotInSubset
{
	public:
		IsNotInSubset(const ISubsetHandler& sh, int subsetIndex) :
			m_sh(sh),
			m_si(subsetIndex)	{}

		bool operator() (VertexBase* v)	{return m_sh.get_subset_index(v) != m_si;}
		bool operator() (EdgeBase* e)	{return m_sh.get_subset_index(e) != m_si;}
		bool operator() (Face* f)		{return m_sh.get_subset_index(f) != m_si;}
		bool operator() (Volume* v)		{return m_sh.get_subset_index(v) != m_si;}

	private:
		const ISubsetHandler& m_sh;
		int m_si;
};

}//	end of namespace

#endif
