//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24

#ifndef __H__LIBGRID__SUBSET_HANDLER_INTERFACE__
#define __H__LIBGRID__SUBSET_HANDLER_INTERFACE__

#include <list>
#include <string>
#include <vector>
#include "grid/grid.h"
#include "common_attachments.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SubsetHandlerElements
///	Use these constants to specify which elements shall be supported by a SubsetHandler.
/**
 * You may combine the constants using or-operations.
 */
enum SubsetHandlerElements
{
	SHE_NONE = 0,
	SHE_VERTEX = 1,
	SHE_EDGE = 1<<1,
	SHE_FACE = 1<<2,
	SHE_VOLUME = 1 << 3,
	SHE_ALL = SHE_VERTEX | SHE_EDGE | SHE_FACE | SHE_VOLUME
};

////////////////////////////////////////////////////////////////////////
//	SubsetState
///	The SubsetState is not yet really used inside of libGrid.
/**
 * The main reason why a SubsetState is introduced, is that
 * applications that use libGrid need a mechanism to store
 * information in a subset.
 * It would be a good idea to think about an attachment-like system
 * for subsets.
 */
enum SubsetState
{
	SS_NONE = 0,
	SS_USER_STATE = 1 << 16
};

////////////////////////////////////////////////////////////////////////
//	SubsetInfo
///	a struct that holds information associated with subsets.
/**
 * In the moment a SubsetInfo is a collection of various types.
 * None of them are really required for libGrid (indeed only name and
 * materialIndex are used in the moment).
 * The other variables are introduced mainly for applications that use
 * libGrid. This is not the best way to do this!
 * It would be a good idea to think about an attachment-like system
 * for subsets.
 */
struct SubsetInfo
{
	SubsetInfo();
	std::string	name;
	int			materialIndex;
	vector4		color;
	uint		subsetState;///< an or-combination of SubsetState flags.
};

////////////////////////////////////////////////////////////////////////
//	ISubsetHandler
/**
 * A derived class has to implement the following public methods:
 * <code>
 * virtual void assign_subset(VertexBase* elem, int subsetIndex)
 * virtual void assign_subset(EdgeBase* elem, int subsetIndex)
 * virtual void assign_subset(Face* elem, int subsetIndex)
 * virtual void assign_subset(Volume* elem, int subsetIndex)
 * </code>
 * In those methods 
 * 
 * Derived classes have to store the objects that are selected for
 * a subset in a ISubsetHandler::SectionContainer. Note that multiple
 * SectionContainers per subset may be used.
 */
class ISubsetHandler : public GridObserver
{
	public:
	///	pass an or-combination of SubsetHandlerElements to supportedElements.
	/**	supportedElements define the elements on which the SubsetHandler works.
	 *	Default is SHE_ALL (all element-types).*/
		ISubsetHandler(uint supportedElements = SHE_ALL);

	///	copy-constructor is not implemented correctly in the moment.
		ISubsetHandler(const ISubsetHandler& sh);

	/**	The destructor automatically unregisters the subset-handler from the grid.
	 *	on deregistration erase_subset_lists of the derived class will be called.*/
		virtual ~ISubsetHandler();

	///	returns a pointer to the grid on which the subset-handler works.
	/**	returns NULL if no grid is assigned.*/
		Grid* get_assigned_grid();
		
	///	returns true if the given element-types are supported.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		bool elements_are_supported(uint shElements);

	///	set the type of elements that shall be handled by the SubsetHandler.
	/**	Pass an or-combination of constants enumerated in SubsetHandlerElements.
	 *	\sa SubsetHandler::enable_element_support*/
		void set_supported_elements(uint shElements);

	///	enable support for element-types. Does not invalidate previous settings.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void enable_element_support(uint shElements);

	///	disable support for element-types.
	/**	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void disable_element_support(uint shElements);

	/**	new elements will be automatically assigned to this subset.
	 * 	set this to a negative value to avoid automatic assignment (-1 by default).
	 *	only used if subset_inheritance is disabled or if no parent is specified.*/
		void set_default_subset_index(int subsetIndex);
		inline int get_default_subset_index()	{return m_defaultSubsetIndex;}

	/**	if enabled, newly created elements derive their subset-index from their parents.
	 *	Enabled by default.
	 *	If enabled, the default subset index will be ignored if a parent is specified
	 *	on element creation.*/
		void enable_subset_inheritance(bool bEnable);
		inline bool subset_inheritance_enabled()	{return m_bSubsetInheritanceEnabled;}
		
	///	if the subset with the given index does not yet exist, it will be created.
	/**	All subsets in between num_subsets and index will be created, too.*/
		inline void subset_info_required(int index);

	///	returns the number of subset-infos
		inline uint num_subset_infos() const		{return (uint)m_subsetInfos.size();}
		
	/** if the subset at subsetIndex does not yet exist, it will be created.*/
		void set_subset_info(int subsetIndex, const SubsetInfo& subsetInfo);

	/** if the subset at subsetIndex does not yet exist, it will be created.*/
		SubsetInfo& subset_info(int subsetIndex);

	/** if the subset at subsetIndex does not yet exist, it will be created.*/
		const SubsetInfo& subset_info(int subsetIndex) const;
		
		void clear();
		void clear_subset(int subsetIndex);
		void clear_subsets();


	///	inserts a subset at the given index. Moves all other subsets 1 index higher.
		void insert_subset(int subsetIndex);///< changes subset-indices of other subsets.
	///	erases the subset at the given index. Assigns -1 to all entries. Moves all other subsets 1 index up.
		void erase_subset(int subsetIndex);///< changes subset-indices of other subsets.
	///	Swaps the given subsets,
		void swap_subsets(int subsetIndex1, int subsetIndex2);
	///	Moves the subset from index From to index To. Moves all subsets between indexFrom+1 and indexTo in the opposite direction.
		void move_subset(int indexFrom, int indexTo);///< changes subset indices of other subsets.

		template <class TIterator>
		void assign_subset(TIterator iterBegin, TIterator iterEnd, int subsetIndex);

		int get_subset_index(GeometricObject* elem);
		inline int get_subset_index(VertexBase* elem)	{return m_aaSubsetIndexVRT[elem];}
		inline int get_subset_index(EdgeBase* elem)		{return m_aaSubsetIndexEDGE[elem];}
		inline int get_subset_index(Face* elem)			{return m_aaSubsetIndexFACE[elem];}
		inline int get_subset_index(Volume* elem)		{return m_aaSubsetIndexVOL[elem];}

	//	grid callbacks
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
		virtual void elements_to_be_cleared(Grid* grid);

	//	vertex callbacks
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL);
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt);

	//	edge callbacks
		virtual void edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent = NULL);
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* edge);

	//	face callbacks
		virtual void face_created(Grid* grid, Face* face, GeometricObject* pParent = NULL);
		virtual void face_to_be_erased(Grid* grid, Face* face);

	//	volume callbacks
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL);
		virtual void volume_to_be_erased(Grid* grid, Volume* vol);
		
	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(VertexBase* elem, int subsetIndex) = 0;
	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(EdgeBase* elem, int subsetIndex) = 0;
	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(Face* elem, int subsetIndex) = 0;
	/**	The implementation in a derived class should store the element in a list
	 *	and call subset_assigned with the iterators position and the subset-index.
	 *	The iterator can later be retrieved with get_list_iterator(...).
	 *	The index can be retrieved with get_subset_index(...).*/
		virtual void assign_subset(Volume* elem, int subsetIndex) = 0;
	
	protected:
		typedef SectionContainer<GeometricObject*, std::list<GeometricObject*> >	SectionContainer;
		typedef SectionContainer::iterator iterator;
		
	protected:
	///	set the grid on which the subset-handler shall work.
	/**	The subset-handler can only work on one grid at a time.
	 *	It is cruicial that assign_grid methods of derived classes call
	 *	this method.*/
		void assign_grid(Grid& grid);
		
	///	sets the subset-indices of all elements of m_pGrid to -1.
	/**	Use with care! Only indices are affected. The elements are not
	 *	removed from any lists.
	 *	pass an or-combination of constants enumerated in SubsetHandlerElements.*/
		void reset_subset_indices(uint shElements = SHE_ALL);
		
		inline void subset_assigned(VertexBase* v, iterator iter, int subsetIndex);
		inline void subset_assigned(EdgeBase* e, iterator iter, int subsetIndex);
		inline void subset_assigned(Face* f, iterator iter, int subsetIndex);
		inline void subset_assigned(Volume* v, iterator iter, int subsetIndex);
		
		inline iterator get_list_iterator(VertexBase* v)	{return m_aaIteratorVRT[v];}
		inline iterator get_list_iterator(EdgeBase* e)		{return m_aaIteratorEDGE[e];}
		inline iterator get_list_iterator(Face* f)			{return m_aaIteratorFACE[f];}
		inline iterator get_list_iterator(Volume* v)		{return m_aaIteratorVOL[v];}
		
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/	
		inline void alter_subset_index(VertexBase* v, int subsetIndex)	{m_aaSubsetIndexVRT[v] = subsetIndex;}
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/
		inline void alter_subset_index(EdgeBase* e, int subsetIndex)	{m_aaSubsetIndexEDGE[e] = subsetIndex;}
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/
		inline void alter_subset_index(Face* f, int subsetIndex)		{m_aaSubsetIndexFACE[f] = subsetIndex;}
	/**	alters the subset index only. Suited as a helper for methods like 
	 *	change_subset_indices or reset_subset_indices.
	 *	WARNING: This method only alters the index but does not actually
	 *	move the element to another subset. Use assign_subset instead for this task.*/
		inline void alter_subset_index(Volume* v, int subsetIndex)		{m_aaSubsetIndexVOL[v] = subsetIndex;}
		
		virtual void erase_subset_lists() = 0;
		
		virtual void clear_subset_lists(int index) = 0;
		
		virtual void change_subset_indices(int indOld, int indNew) = 0;
		
		
	///	add a subset if requiered - so that the subset with maxIndex exists.
		virtual void add_required_subset_lists(int maxIndex) = 0;
		
	///	erase the subset-lists but do not touch the subset-indices.
		virtual void erase_subset_lists(int index) = 0;

	///	swap the subset-lists but do not touch the subset-indices.
		virtual void swap_subset_lists(int ind1, int ind2) = 0;

	///	move the subset-lists but do not touch the subset-indices.
		virtual void move_subset_lists(int indexFrom, int indexTo) = 0;

		
	protected:
		typedef AInt	ASubsetIndex;
		typedef Attachment<iterator>	AIterator;
		typedef std::vector<SubsetInfo>	SubsetInfoVec;

	protected:
		Grid*			m_pGrid;
		SubsetInfoVec	m_subsetInfos;
		uint			m_supportedElements;

		ASubsetIndex	m_aSubsetIndex;
		AIterator		m_aIterator;
		
		int				m_defaultSubsetIndex;
		bool			m_bSubsetInheritanceEnabled;
				
		Grid::VertexAttachmentAccessor<ASubsetIndex>	m_aaSubsetIndexVRT;
		Grid::EdgeAttachmentAccessor<ASubsetIndex>		m_aaSubsetIndexEDGE;
		Grid::FaceAttachmentAccessor<ASubsetIndex>		m_aaSubsetIndexFACE;
		Grid::VolumeAttachmentAccessor<ASubsetIndex>	m_aaSubsetIndexVOL;

		Grid::VertexAttachmentAccessor<AIterator>		m_aaIteratorVRT;
		Grid::EdgeAttachmentAccessor<AIterator>			m_aaIteratorEDGE;
		Grid::FaceAttachmentAccessor<AIterator>			m_aaIteratorFACE;
		Grid::VolumeAttachmentAccessor<AIterator>		m_aaIteratorVOL;
};


inline void ISubsetHandler::
subset_assigned(VertexBase* v, iterator iter, int subsetIndex)
{
	m_aaIteratorVRT[v] = iter;
	m_aaSubsetIndexVRT[v] = subsetIndex;
}

inline void ISubsetHandler::
subset_assigned(EdgeBase* e, iterator iter, int subsetIndex)
{
	m_aaIteratorEDGE[e] = iter;
	m_aaSubsetIndexEDGE[e] = subsetIndex;
}

inline void 
ISubsetHandler::
subset_assigned(Face* f, iterator iter, int subsetIndex)
{
	m_aaIteratorFACE[f] = iter;
	m_aaSubsetIndexFACE[f] = subsetIndex;
}

inline void 
ISubsetHandler::
subset_assigned(Volume* v, iterator iter, int subsetIndex)
{
	m_aaIteratorVOL[v] = iter;
	m_aaSubsetIndexVOL[v] = subsetIndex;
}

template <class TIterator>
void ISubsetHandler::
assign_subset(TIterator iterBegin, TIterator iterEnd, int subsetIndex)
{
	typename TIterator::value_type elem;
	while(iterBegin != iterEnd)
	{
		elem = *iterBegin;
		++iterBegin;
		assign_subset(elem, subsetIndex);
	}
}

inline void
ISubsetHandler::
subset_info_required(int index)
{
	if(index >= (int)m_subsetInfos.size())
	{
		m_subsetInfos.resize(index+1);
		add_required_subset_lists(index);
	}
}

}//	end of namespace

////////////////////////////////////////////////
//	include implementation
//#include "subset_handler_base_impl.hpp"

#endif
