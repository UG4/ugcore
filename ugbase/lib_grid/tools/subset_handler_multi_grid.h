//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d24   (reworked y09 m12 d15)

#ifndef __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID__
#define __H__LIBGRID__SUBSET_HANDLER_MULTI_GRID__

#include <vector>
#include <cassert>
#include "lib_grid/multi_grid.h"
#include "common/util/section_container.h"
#include "subset_handler_interface.h"
#include "../lib_grid_messages.h"

namespace ug
{

/** \ingroup lib_grid_tools
 *  \{ */

////////////////////////////////////////////////////////////////////////
//	MultiGridSubsetHandler
/// Handles subsets on a per level basis.
/** The MultiGridSubsetHandler is a specialization of ISubsetHandler for
 * MultiGrids. It allows to access elements given a subset-index and a level index.
 *
 * Note that the number of levels in the MultiGridSubsetHandler always matches
 * the number of levels in the associated multigrid. This is guaranteed through
 * a callback mechanism.
 */
class UG_API MultiGridSubsetHandler : public ISubsetHandler
{
	public:
		using ISubsetHandler::assign_subset;
		
	public:
		MultiGridSubsetHandler(uint supportedElements = SHE_ALL);
		MultiGridSubsetHandler(MultiGrid& mg, uint supportedElements = SHE_ALL);
	/**	WARNING: Don't call the copy-constructor from derived classes,
	  *	Since it calls virtual methods.*/
		MultiGridSubsetHandler(const MultiGridSubsetHandler& sh);
		virtual ~MultiGridSubsetHandler();
				
		void assign_grid(MultiGrid& mg);
		inline MultiGrid* multi_grid()	{return m_pMG;}
		inline const MultiGrid* multi_grid() const {return m_pMG;}

	///	creates the required levels, if they do not yet exist
		inline void level_required(size_t level);

	///	returns the number of levels
		inline uint num_levels() const	{return (uint)m_levels.size();}
		
	///	returns the level in which an element is located
		template <class TGeomObj>
		inline uint get_level(TGeomObj* obj) const	{return m_pMG->get_level(obj);}
		
	////////////////////////////////////////////////
	//	implementation of public virtual methdos of ISubsetHandler.
	///	assigns a vertex to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Vertex* elem, int subsetIndex);

	///	assigns an edge to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(EdgeBase* elem, int subsetIndex);

	///	assigns a face to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Face* elem, int subsetIndex);

	///	assigns a volume to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Volume* elem, int subsetIndex);

	////////////////////////////////////////////////
	//	element-access
	///	returns the begin-iterator for the elements of type TElem in the given subset.
	/**	e.g. begin<Triangle>(0)*/
		template <class TElem>
		typename geometry_traits<TElem>::iterator
		begin(int subsetIndex, int level);

	///	returns the end-iterator for the elements of type TElem in the given subset.
	/**	e.g. end<Triangle>(0)*/
		template <class TElem>
		typename geometry_traits<TElem>::iterator
		end(int subsetIndex, int level);

	///	returns the begin-iterator for the elements of type TElem in the given subset.
	/**	e.g. begin<Triangle>(0)
	 *	Please note that in the const version level < num_levels() has to hold true.*/
		template <class TElem>
		typename geometry_traits<TElem>::const_iterator
		begin(int subsetIndex, int level) const;

	///	returns the end-iterator for the elements of type TElem in the given subset.
	/**	e.g. end<Triangle>(0)
	 *	Please note that in the const version level < num_levels() has to hold true.*/
		template <class TElem>
		typename geometry_traits<TElem>::const_iterator
		end(int subsetIndex, int level) const;
		
	///	returns the total number of elements
		template <class TElem>
		uint num() const;
		
	///	returns the number of elements in the given subset
		template <class TElem>
		uint num(int subsetIndex) const;

	///	returns the number of elements in the given subset on the given level
		template <class TElem>
		uint num(int subsetIndex, int level) const;

	///	removes all elements of type TElem from the specified subset.
		template <class TElem>
		void clear_subset_elements(int subsetIndex);

	///	removes all elements of type TElem from the specified subset on the given level.
		template <class TElem>
		void clear_subset_elements(int subsetIndex, int level);

	///	returns a GridObjectCollection
	/**	the returned GridObjectCollection hold the elements of the
	 *	specified subset on the given level.*/
		GridObjectCollection
		get_grid_objects(int subsetIndex, int level) const;
		
	///	returns a GridObjectCollection with multiple levels
	/**	the returned GridObjectCollection hold the
	 *	elements of the specified subset.*/
		GridObjectCollection
		get_grid_objects_in_subset(int subsetIndex) const;

	///	returns a GridObjectCollection with multiple levels - each representing a subset.
	/**	the returned GridObjectCollection hold the
	 *	elements of the specified level, each level of the collection
	 *	represents a subset.*/
		GridObjectCollection
		get_grid_objects_in_level(int level) const;
		
	///	collects all vertices that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has compexity O(n), where n is the number
	 *	of vertices in the underlying grid.
	 *	\returns number of collected elements.
	 *	\sa begin, end*/
		//virtual size_t collect_subset_elements(std::vector<Vertex*>& vrtsOut, int subsetIndex) const;

	///	collects all edges that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has compexity O(n), where n is the number
	 *	of edges in the underlying grid.
	 *	\returns number of collected elements.
	 *	\sa begin, end*/
		//virtual size_t collect_subset_elements(std::vector<EdgeBase*>& edgesOut, int subsetIndex) const;

	///	collects all faces that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has compexity O(n), where n is the number
	 *	of faces in the underlying grid.
	 *	\returns number of collected elements.
	 *	\sa begin, end*/
		//virtual size_t collect_subset_elements(std::vector<Face*>& facesOut, int subsetIndex) const;

	///	collects all volumes that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has compexity O(n), where n is the number
	 *	of volumes in the underlying grid.
	 *	\returns number of collected elements.
	 *	\sa begin, end*/
		//virtual size_t collect_subset_elements(std::vector<Volume*>& volsOut, int subsetIndex) const;
		
	///	returns true if the subset contains vertices
		virtual bool contains_vertices(int subsetIndex) const	{return num<Vertex>(subsetIndex) > 0;}

	///	returns true if the subset contains edges
		virtual bool contains_edges(int subsetIndex) const		{return num<EdgeBase>(subsetIndex) > 0;}
		
	///	returns true if the subset contains faces
		virtual bool contains_faces(int subsetIndex) const		{return num<Face>(subsetIndex) > 0;}
		
	///	returns true if the subset contains volumes
		virtual bool contains_volumes(int subsetIndex) const	{return num<Volume>(subsetIndex) > 0;}


	///	perform cleanup
		virtual void grid_to_be_destroyed(Grid* grid);

	protected:
	///	returns the number of subsets in the local list
		inline uint num_subsets_in_list() const	{return m_numSubsets;}
		
	///	detaches all attached data.
		void detach_data();

	////////////////////////////////////////////////
	//	implementation of protected virtual methdos of ISubsetHandler.
	///	erases the subsets. Doesn't alter any indices.
		virtual void erase_subset_lists();

	///	non-virtual implementation of erase_subset_lists. Callable from destructor
		void erase_subset_lists_impl();
		
	///	clears the element lists in the given subset. Does not alter any indices.
		void clear_subset_lists(int index);

	///	changes the subset-indices of all elements int the subset.
	/**	WARNING: subsets are not automatically changed accordingly.
	 *	After termination Subset-Indices and Subset-Infos/iterators are asynchronous.
	 *	Make sure to change subset-infos and iterators accordingly.*/		
		void change_subset_indices(int indOld, int indNew);		

		
	///	add a subset
		void add_required_subset_lists(int maxIndex);
		
	///	erases the subset but does not touch the subset-indices.
		void erase_subset_lists(int index);

	///	swaps the subsets but does not touch the subset-indices.
		void swap_subset_lists(int ind1, int ind2);

	///	moves the subset but does not touch the subset-indices.
		void move_subset_lists(int indexFrom, int indexTo);

	///	join the subset-lists but do not touch the subset-indices.
		virtual void join_subset_lists(int target, int src1, int src2);

	///	this method is called by ISubsetHandler when attachment_support has been enabled.
		//void register_subset_elements_at_pipe();

	////////////////////////////////////////////////
	//	protected helper methods
	///	a helper method for the public assign_subset methods.
		template<class TElem>
		void assign_subset_impl(TElem* elem, int subsetIndex);

	///	helper for change_subset_indices
		template<class TElem>
		void change_elem_subset_indices(int indOld, int indNew);
		
	///	Throws an error if the required level does not yet exist
		inline void level_required(size_t level) const;

		void add_level();
		void add_subset_to_all_levels();///< increases m_numSubsets.

	///	helper for collect_subset_elements
		template <class TElem>
		size_t collect_subset_elements_impl(std::vector<TElem*>& elemsOut, int subsetIndex) const;
		
	protected:
		using ISubsetHandler::AttachedVertexList;
		using ISubsetHandler::AttachedEdgeList;
		using ISubsetHandler::AttachedFaceList;
		using ISubsetHandler::AttachedVolumeList;

		using ISubsetHandler::VertexSectionContainer;
		using ISubsetHandler::EdgeSectionContainer;
		using ISubsetHandler::FaceSectionContainer;
		using ISubsetHandler::VolumeSectionContainer;
		
		struct Subset
		{
			VertexSectionContainer	m_vertices;
			EdgeSectionContainer	m_edges;
			FaceSectionContainer	m_faces;
			VolumeSectionContainer	m_volumes;
		};

		typedef std::vector<Subset*>	SubsetVec;
		typedef std::vector<SubsetVec>	LevelVec;

	///	returns the subset with index si on the given level
		inline Subset* subset(int si, int level)	{return m_levels[level][si];}
		inline const Subset* subset(int si, int level)	const {return m_levels[level][si];}

	///	creates a new subset. Caller is responsible for deletion
		Subset* new_subset();

		void cleanup();

	///	returns the iterator at which the given element lies in the section container
	/**	This method may only be called if the element is in a subset != -1.
	 * \{
	 */
		inline VertexSectionContainer::iterator
		get_list_iterator(Vertex* o)
		{
			assert((get_subset_index(o) >= 0) && "invalid subset.");
			return subset(get_subset_index(o), m_pMG->get_level(o))->
					m_vertices.get_container().get_iterator(o);
		}

		inline EdgeSectionContainer::iterator
		get_list_iterator(EdgeBase* o)
		{
			assert((get_subset_index(o) >= 0) && "invalid subset.");
			return subset(get_subset_index(o), m_pMG->get_level(o))->
					m_edges.get_container().get_iterator(o);
		}

		inline FaceSectionContainer::iterator
		get_list_iterator(Face* o)
		{
			assert((get_subset_index(o) >= 0) && "invalid subset.");
			return subset(get_subset_index(o), m_pMG->get_level(o))->
					m_faces.get_container().get_iterator(o);
		}

		inline VolumeSectionContainer::iterator
		get_list_iterator(Volume* o)
		{
			assert((get_subset_index(o) >= 0) && "invalid subset.");
			return subset(get_subset_index(o), m_pMG->get_level(o))->
					m_volumes.get_container().get_iterator(o);
		}
	/**	\}	*/

	///	returns the section container for the given type, subset and level
		template <class TElem> inline
		typename Grid::traits<TElem>::SectionContainer&
		section_container(int si, int lvl);

	///	returns the const section container for the given type, subset and level
		template <class TElem> inline
		const typename Grid::traits<TElem>::SectionContainer&
		section_container(int si, int lvl) const;

	///	callback for multigrid messages
		void multigrid_changed(const GridMessage_MultiGridChanged& gm);

	protected:
		MultiGrid*		m_pMG;
		LevelVec		m_levels;
		int				m_numSubsets;

	//	callback-id (automatically unregisters callback, when the selector is deleted).
		MessageHub::SPCallbackId	m_callbackId;

		AttachedVertexList::AEntry	m_aSharedEntryVRT;
		AttachedEdgeList::AEntry	m_aSharedEntryEDGE;
		AttachedFaceList::AEntry	m_aSharedEntryFACE;
		AttachedVolumeList::AEntry	m_aSharedEntryVOL;
};

typedef MultiGridSubsetHandler MGSubsetHandler;

/** \} */
}//	end of namespace

//	include implementation
#include "subset_handler_multi_grid_impl.hpp"

#endif
