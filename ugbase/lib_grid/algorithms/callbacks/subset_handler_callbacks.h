// created by Sebastian Reiter
// y10 m12 d13
// s.b.reiter@googlemail.com

#ifndef __H__LIB_GRID__SUBSET_HANDLER_CALLBACKS__
#define __H__LIB_GRID__SUBSET_HANDLER_CALLBACKS__

#include "lib_grid/tools/subset_handler_interface.h"
#include "callback_definitions.h"
#include "common/serialization.h"

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


///	Serialization callback for subset handlers
/**	Writes subset indices to a binary stream.
 */
/*
class SubsetHandlerSerializer
{
	public:
		SubsetHandlerSerializer(ISubsetHandler& sh) :
			m_sh(sh)	{}

		void operator() (std::ostream& out, VertexBase* v)
		{Serialize(out, m_sh.get_subset_index(v));}

		void operator() (std::ostream& out, EdgeBase* e)
		{Serialize(out, m_sh.get_subset_index(e));}

		void operator() (std::ostream& out, Face* f)
		{Serialize(out, m_sh.get_subset_index(f));}

		void operator() (std::ostream& out, Volume* v)
		{Serialize(out, m_sh.get_subset_index(v));}

	private:
		ISubsetHandler& m_sh;
};
*/
///	Deserialization callback for subset handlers
/**	Reads subset indices from a binary stream and assigns
 * objects accordingly.
 * Make sure that you pass objects in the same order as
 * during serialization.
 */
/*
class SubsetHandlerDeserializer
{
	public:
		SubsetHandlerDeserializer(ISubsetHandler& sh) :
			m_sh(sh)	{}

		void operator() (std::istream& in, VertexBase* v)
		{	int si;
			Deserialize(in, si);
			m_sh.assign_subset(v, si);
		}

		void operator() (std::istream& in, EdgeBase* e)
		{	int si;
			Deserialize(in, si);
			m_sh.assign_subset(e, si);
		}

		void operator() (std::istream& in, Face* f)
		{	int si;
			Deserialize(in, si);
			m_sh.assign_subset(f, si);
		}

		void operator() (std::istream& in, Volume* v)
		{	int si;
			Deserialize(in, si);
			m_sh.assign_subset(v, si);
		}

	private:
		ISubsetHandler& m_sh;
};
*/
}//	end of namespace

#endif
