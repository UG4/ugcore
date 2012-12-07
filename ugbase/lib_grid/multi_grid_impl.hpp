//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d10

#include "multi_grid.h"

#ifndef __H__LIB_GRID__MULTI_GRID_IMPL__
#define __H__LIB_GRID__MULTI_GRID_IMPL__

namespace ug
{
/*
template <class TChild, class TElem>
int MultiGrid::num_children(TElem* elem)
{

	int sharedPipeSec = elem->container_section();
	assert(sharedPipeSec != -1 && "bad shared pipe section!");

	return get_elem_info(elem)->m_children.num_elements(sharedPipeSec);
}

//	child access
template <class TChild, class TElem>
int MultiGrid::get_children(std::vector<TChild*>& vChildrenOut, TElem* elem)
{

	int sharedPipeSec = geometry_traits<TChild>::CONTAINER_SECTION;
	assert(sharedPipeSec != -1 && "bad shared pipe section!");

	MGElementInfo* elemInfo = get_elem_info(elem);

	typename SectionContainer::iterator iter =
				elemInfo->m_children.section_begin(sharedPipeSec);

	typename SectionContainer::iterator iterEnd =
			elemInfo->m_children.section_end(sharedPipeSec);

	int numChildren = elemInfo->m_children.num_elements(sharedPipeSec);

	if(vChildrenOut.capacity() < numChildren)
		vChildrenOut.reserve(numChildren);

	vChildrenOut.clear();

	for(; iter != iterEnd; ++iter)
	{
		vChildrenOut.push_back((TChild*) *iter);
	}

	return vChildrenOut.size();

}
*/

inline size_t MultiGrid::
top_level() const
{
	if(m_hierarchy.num_subsets() <= 0) return 0;
	return size_t(m_hierarchy.num_subsets() - 1);
}

template <class TElem>
size_t MultiGrid::
num_children_total(TElem* elem)	const
{
	size_t numChildren = num_children<TElem>(elem);
	size_t numChildrenTotal = numChildren;

	for(size_t i = 0; i < numChildren; ++i)
		numChildrenTotal += num_children_total(get_child<TElem>(elem, i));

	return numChildrenTotal;
}


template<class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
MultiGrid::create(size_t level)
{
	typename geometry_traits<TGeomObj>::iterator iter =
										Grid::create<TGeomObj>();
//	put the element into the hierarchy
//	(by default it already was assigned to level 0)
	if(level > 0){
		level_required(level);
		m_hierarchy.assign_subset(*iter, level);
	}
	return iter;
}

template <class TGeomObj>
typename geometry_traits<TGeomObj>::iterator
MultiGrid::create(const typename geometry_traits<TGeomObj>::Descriptor& descriptor,
				size_t level)
{
	typename geometry_traits<TGeomObj>::iterator iter =
										Grid::create<TGeomObj>(descriptor);
//	put the element into the hierarchy
//	(by default it already was assigned to level 0)
	if(level > 0){
		level_required(level);
		m_hierarchy.assign_subset(*iter, level);
	}
	return iter;
}

inline void MultiGrid::level_required(int lvl)
{
	if(m_hierarchy.num_subsets() <= lvl){
		create_levels(lvl - m_hierarchy.num_subsets() + 1);
	}
}


template <class TChild>
size_t MultiGrid::num_children(GeometricObject* elem) const
{
	switch(elem->base_object_id()){
	case VERTEX:	return num_children<TChild>(static_cast<VertexBase*>(elem));
	case EDGE:		return num_children<TChild>(static_cast<EdgeBase*>(elem));
	case FACE:		return num_children<TChild>(static_cast<Face*>(elem));
	case VOLUME:	return num_children<TChild>(static_cast<Volume*>(elem));
	}
	return 0;
}

template <class TChild>
TChild* MultiGrid::get_child(GeometricObject* elem, size_t ind) const
{
	switch(elem->base_object_id()){
	case VERTEX:	return get_child<TChild>(static_cast<VertexBase*>(elem), ind);
	case EDGE:		return get_child<TChild>(static_cast<EdgeBase*>(elem), ind);
	case FACE:		return get_child<TChild>(static_cast<Face*>(elem), ind);
	case VOLUME:	return get_child<TChild>(static_cast<Volume*>(elem), ind);
	}
	return NULL;
}

//	info-access
inline MultiGrid::VertexInfo& MultiGrid::get_info(VertexBase* v)
{
	return m_aaVrtInf[v];
}

inline MultiGrid::EdgeInfo& MultiGrid::get_info(EdgeBase* e)
{
	return m_aaEdgeInf[e];
}

inline MultiGrid::FaceInfo& MultiGrid::get_info(Face* f)
{
	if(FaceInfo* info = m_aaFaceInf[f])
		return *info;
	UG_THROW("MultiGrid::get_info(...): No face info available!");
}

inline MultiGrid::VolumeInfo& MultiGrid::get_info(Volume* v)
{
	if(VolumeInfo* info = m_aaVolInf[v])
		return *info;
	UG_THROW("MultiGrid::get_info(...): No vertex info available!");
}

//	const info-access
inline const MultiGrid::VertexInfo& MultiGrid::get_info(VertexBase* v) const
{
	return m_aaVrtInf[v];
}

inline const MultiGrid::EdgeInfo& MultiGrid::get_info(EdgeBase* e) const
{
	return m_aaEdgeInf[e];
}

inline const MultiGrid::FaceInfo& MultiGrid::get_info(Face* f) const
{
	static FaceInfo	emptyInfo;
	if(FaceInfo* info = m_aaFaceInf[f])
		return *info;
	return emptyInfo;
}

inline const MultiGrid::VolumeInfo& MultiGrid::get_info(Volume* v) const
{
	static VolumeInfo	emptyInfo;
	if(VolumeInfo* info = m_aaVolInf[v])
		return *info;
	return emptyInfo;
}


template <class TElem, class TParent>
void MultiGrid::element_created(TElem* elem, TParent* pParent)
{
//	if hierarchical_insertion is enabled, the element will be put
//	into the next higher level of pParents level.

	int level = 0;
	if(pParent)
	{
	//	the element is inserted into a new layer.
		level = get_level(pParent) + 1;
	}

//	register parent and child
	//typename mginfo_traits<TElem>::info_type& info = get_info(elem);
	//info.m_pParent = pParent;
	set_parent(elem, pParent);
	if(pParent)
	{
	//	make sure that the parent has an info object
		create_child_info(pParent);

	//	add the element to the parents children list
		typename mginfo_traits<TParent>::info_type& parentInfo = get_info(pParent);
		parentInfo.add_child(elem);
	}

//	put the element into the hierarchy
	level_required(level);
	m_hierarchy.assign_subset(elem, level);
}

template <class TElem, class TParent>
void MultiGrid::element_created(TElem* elem, TParent* pParent,
								TElem* pReplaceMe)
{
//	if hierarchical_insertion is enabled, the element will be put
//	into the next higher level of pParents level.

	int level = 0;
	if(pParent)
	{
	//	the element is inserted into a new layer.
		level = get_level(pParent) + 1;
	}

//	register parent and child
	set_parent(elem, pParent);

	if(pParent)
	{
	//	add the element to the parents children list
	//	pParent should have an info object at this time!
		typename mginfo_traits<TParent>::info_type& parentInfo = get_info(pParent);
		parentInfo.replace_child(elem, pReplaceMe);
	}

//	put the element into the hierarchy
	level_required(level);
	m_hierarchy.assign_subset(elem, level);
}

template <class TElem>
void MultiGrid::element_to_be_erased(TElem* elem)
{
//	we have to remove the elements children as well.
	if(has_children(elem)){
		get_info(elem).unregister_from_children(*this);
	}
//	we have to remove the associated info object
	release_child_info(elem);
}

template <class TElem, class TParent>
void MultiGrid::element_to_be_erased(TElem* elem, TParent* pParent)
{
//	unregister the element from its parent.
//	parents always have an info object
	typename mginfo_traits<TParent>::info_type& parentInfo = get_info(pParent);
	parentInfo.remove_child(elem);
	element_to_be_erased(elem);
	if(!parentInfo.has_children()){
		release_child_info(pParent);
	}
}

/*
template <class TElem>
void MultiGrid::element_to_be_replaced(TElem* elemOld, TElem* elemNew)
{
}
*/



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	specialization of wrapper classes

////////////////////////////////////////////////////////////////////////
//	specialization for Grid
template <>
class MGWrapper<Grid>
{
	public:
		MGWrapper(Grid& grid) : m_grid(grid)	{}
		
		inline uint num_levels() const
		{return 1;}

		template <class TElem> inline
		uint num(int level) const
		{return m_grid.num<TElem>();}

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		begin(int level)
		{return m_grid.begin<TElem>();}

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		end(int level)
		{return m_grid.end<TElem>();}

	protected:
		Grid&	m_grid;
};

////////////////////////////////////////////////////////////////////////
//	specialization for MultiGrid
template <>
class MGWrapper<MultiGrid>
{
	public:
		MGWrapper(MultiGrid& grid) : m_grid(grid)	{}
		
		inline uint num_levels() const
		{return m_grid.num_levels();}

		template <class TElem> inline
		uint num(int level) const
		{return m_grid.num<TElem>(level);}

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		begin(int level)
		{return m_grid.begin<TElem>(level);}

		template <class TElem> inline
		typename geometry_traits<TElem>::iterator
		end(int level)
		{return m_grid.end<TElem>(level);}

	protected:
		MultiGrid&	m_grid;
};

}//	end of namespace

#endif
