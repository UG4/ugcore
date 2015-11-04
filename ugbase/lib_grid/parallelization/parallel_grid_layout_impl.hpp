#ifndef __H__LIB_GRID__PARALLEL_GRID_LAYOUT_IMPL__
#define __H__LIB_GRID__PARALLEL_GRID_LAYOUT_IMPL__

#include "parallel_grid_layout.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	GridLayoutMap - implementation
template <class TType>
bool GridLayoutMap::
has_layout(const Key& key) const
{
	const typename Types<TType>::Map& m = get_layout_map<TType>();
	return m.find(key) != m.end();
}

template <class TType>
typename GridLayoutMap::Types<TType>::Layout& GridLayoutMap::
get_layout(const Key& key)
{
	typename Types<TType>::Map& m = get_layout_map<TType>();
	return m[key];
}

template <class TType>
const typename GridLayoutMap::Types<TType>::Layout& GridLayoutMap::
get_layout(const Key& key) const
{
	const typename Types<TType>::Map& m = get_layout_map<TType>();
	typename Types<TType>::Map::const_iterator iter = m.find(key);
	if(iter == m.end()){
		UG_THROW("The specified layout can not be found in the GridLayoutMap.");
	}
	return iter->second;
}

template <class TType>
typename GridLayoutMap::Types<TType>::Map::iterator GridLayoutMap::
layouts_begin()
{
	return get_layout_map<TType>().begin();
}

template <class TType>
typename GridLayoutMap::Types<TType>::Map::const_iterator GridLayoutMap::
layouts_begin() const
{
	return get_layout_map<TType>().begin();
}

template <class TType>
typename GridLayoutMap::Types<TType>::Map::iterator GridLayoutMap::
layouts_end()
{
	return get_layout_map<TType>().end();
}

template <class TType>
typename GridLayoutMap::Types<TType>::Map::const_iterator GridLayoutMap::
layouts_end() const
{
	return get_layout_map<TType>().end();
}

template <class TType>
typename GridLayoutMap::Types<TType>::Map::iterator GridLayoutMap::
erase_layout(typename GridLayoutMap::Types<TType>::Map::iterator iter)
{
	typename Types<TType>::Map& m = get_layout_map<TType>();
	typename Types<TType>::Map::iterator tIter = iter++;
	m.erase(tIter);
	return iter;
}

template <class TType>
void GridLayoutMap::
erase_layout(const Key& key)
{
	typename Types<TType>::Map& m = get_layout_map<TType>();
	typename Types<TType>::Map::iterator iter = m.find(key);
	if(iter != m.end())
		m.erase(iter);
}

inline void GridLayoutMap::clear()
{
	m_vertexLayoutMap = Types<Vertex>::Map();
	m_edgeLayoutMap = Types<Edge>::Map();
	m_faceLayoutMap = Types<Face>::Map();
	m_volumeLayoutMap = Types<Volume>::Map();
}

template <class TType>
inline typename GridLayoutMap::Types<TType>::Map& GridLayoutMap::
get_layout_map()
{
	TType* dummy = NULL;//	set to NULL to avoid warnings
	return get_layout_map(dummy);
}

template <class TType>
inline const typename GridLayoutMap::Types<TType>::Map& GridLayoutMap::
get_layout_map() const
{
	TType* dummy = NULL;//	set to NULL to avoid warnings
	return get_layout_map(dummy);
}

inline GridLayoutMap::Types<Vertex>::Map& GridLayoutMap::
get_layout_map(Vertex*)	
{
	return m_vertexLayoutMap;
}

inline const GridLayoutMap::Types<Vertex>::Map& GridLayoutMap::
get_layout_map(Vertex*) const
{
	return m_vertexLayoutMap;
}

inline GridLayoutMap::Types<Edge>::Map& GridLayoutMap::
get_layout_map(Edge*)
{
	return m_edgeLayoutMap;
}

inline const GridLayoutMap::Types<Edge>::Map& GridLayoutMap::
get_layout_map(Edge*) const
{
	return m_edgeLayoutMap;
}

inline GridLayoutMap::Types<Face>::Map& GridLayoutMap::
get_layout_map(Face*)
{
	return m_faceLayoutMap;
}

inline const GridLayoutMap::Types<Face>::Map& GridLayoutMap::
get_layout_map(Face*) const
{
	return m_faceLayoutMap;
}

inline GridLayoutMap::Types<Volume>::Map& GridLayoutMap::
get_layout_map(Volume*)
{
	return m_volumeLayoutMap;
}

inline const GridLayoutMap::Types<Volume>::Map& GridLayoutMap::
get_layout_map(Volume*) const
{
	return m_volumeLayoutMap;
}

}//	end of namespace

#endif
