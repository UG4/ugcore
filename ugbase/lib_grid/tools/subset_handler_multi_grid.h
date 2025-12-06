/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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
		explicit MultiGridSubsetHandler(uint supportedElements = SubsetHandlerElements::SHE_ALL);
		explicit MultiGridSubsetHandler(MultiGrid& mg, uint supportedElements = SubsetHandlerElements::SHE_ALL);
	/**	WARNING: Don't call the copy-constructor from derived classes,
	  *	Since it calls virtual methods.*/
		MultiGridSubsetHandler(const MultiGridSubsetHandler& sh);

		~MultiGridSubsetHandler() override;
				
		void assign_grid(MultiGrid& mg);
		inline MultiGrid* multi_grid()	{return m_pMG;}
		[[nodiscard]] inline const MultiGrid* multi_grid() const {return m_pMG;}

	///	creates the required levels, if they do not yet exist
		inline void level_required(int level);

	///	returns the number of levels
		[[nodiscard]] inline uint num_levels() const {return static_cast<uint>(m_levels.size());}
		
	///	returns the level in which an element is located
		template <typename TGeomObj>
		inline uint get_level(TGeomObj* obj) const	{return m_pMG->get_level(obj);}
		
	////////////////////////////////////////////////
	//	implementation of public virtual methods of ISubsetHandler.
	///	assigns a vertex to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Vertex* elem, int subsetIndex) override;

	///	assigns an edge to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Edge* elem, int subsetIndex) override;

	///	assigns a face to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Face* elem, int subsetIndex) override;

	///	assigns a volume to a subset.
	/**	If the subset doesn't already exist, it will be created.*/
		void assign_subset(Volume* elem, int subsetIndex) override;

	////////////////////////////////////////////////
	//	element-access
	///	returns the begin-iterator for the elements of type TElem in the given subset.
	/**	e.g. begin<Triangle>(0)*/
		template <typename TElem>
		typename geometry_traits<TElem>::iterator
		begin(int subsetIndex, int level);

	///	returns the end-iterator for the elements of type TElem in the given subset.
	/**	e.g. end<Triangle>(0)*/
		template <typename TElem>
		typename geometry_traits<TElem>::iterator
		end(int subsetIndex, int level);

	///	returns the begin-iterator for the elements of type TElem in the given subset.
	/**	e.g. begin<Triangle>(0)
	 *	Please note that in the const version level < num_levels() has to hold true.*/
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator
		begin(int subsetIndex, int level) const;

	///	returns the end-iterator for the elements of type TElem in the given subset.
	/**	e.g. end<Triangle>(0)
	 *	Please note that in the const version level < num_levels() has to hold true.*/
		template <typename TElem>
		typename geometry_traits<TElem>::const_iterator
		end(int subsetIndex, int level) const;
		
	///	returns the total number of elements
		template <typename TElem>
		[[nodiscard]] uint num() const;
		
	///	returns the number of elements in the given subset
		template <typename TElem>
		[[nodiscard]] uint num(int subsetIndex) const;

	///	returns the number of elements in the given subset on the given level
		template <typename TElem>
		[[nodiscard]] uint num(int subsetIndex, int level) const;

	///	removes all elements of type TElem from the specified subset.
		template <typename TElem>
		void clear_subset_elements(int subsetIndex);

	///	removes all elements of type TElem from the specified subset on the given level.
		template <typename TElem>
		void clear_subset_elements(int subsetIndex, int level);

	///	returns a GridObjectCollection
	/**	the returned GridObjectCollection hold the elements of the
	 *	specified subset on the given level.*/
		[[nodiscard]] GridObjectCollection
		get_grid_objects(int subsetIndex, int level) const;
		
	///	returns a GridObjectCollection with multiple levels
	/**	the returned GridObjectCollection hold the
	 *	elements of the specified subset.*/
		[[nodiscard]] GridObjectCollection
		get_grid_objects_in_subset(int subsetIndex) const override;

	///	returns a GridObjectCollection with multiple levels - each representing a subset.
	/**	the returned GridObjectCollection hold the
	 *	elements of the specified level, each level of the collection
	 *	represents a subset.*/
		[[nodiscard]] GridObjectCollection
		get_grid_objects_in_level(int level) const;
		
	///	collects all vertices that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has complexity O(n), where n is the number
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
		//virtual size_t collect_subset_elements(std::vector<Edge*>& edgesOut, int subsetIndex) const;

	///	collects all faces that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has complexity O(n), where n is the number
	 *	of faces in the underlying grid.
	 *	\returns number of collected elements.
	 *	\sa begin, end*/
		//virtual size_t collect_subset_elements(std::vector<Face*>& facesOut, int subsetIndex) const;

	///	collects all volumes that are in the given subset.
	/**	Please consider using begin and end methods instead.
	 *	If subset -1 is specified, the method has complexity O(n), where n is the number
	 *	of volumes in the underlying grid.
	 *	\returns number of collected elements.
	 *	\sa begin, end*/
		//virtual size_t collect_subset_elements(std::vector<Volume*>& volsOut, int subsetIndex) const;
		
	///	returns true if the subset contains vertices
		[[nodiscard]] bool contains_vertices(int subsetIndex) const override {return num<Vertex>(subsetIndex) > 0;}

	///	returns true if the subset contains edges
		[[nodiscard]] bool contains_edges(int subsetIndex) const override {return num<Edge>(subsetIndex) > 0;}
		
	///	returns true if the subset contains faces
		[[nodiscard]] bool contains_faces(int subsetIndex) const override {return num<Face>(subsetIndex) > 0;}
		
	///	returns true if the subset contains volumes
		[[nodiscard]] bool contains_volumes(int subsetIndex) const override {return num<Volume>(subsetIndex) > 0;}


	///	perform cleanup
		void grid_to_be_destroyed(Grid* grid) override;

	protected:
	///	returns the number of subsets in the local list
		[[nodiscard]] inline uint num_subsets_in_list() const {return m_numSubsets;}
		
	///	detaches all attached data.
		void detach_data();

	////////////////////////////////////////////////
	//	implementation of protected virtual methods of ISubsetHandler.
	///	erases the subsets. Doesn't alter any indices.
		void erase_subset_lists() override;

	///	non-virtual implementation of erase_subset_lists. Callable from destructor
		void erase_subset_lists_impl();
		
	///	clears the element lists in the given subset. Does not alter any indices.
		void clear_subset_lists(int index) override;

	///	changes the subset-indices of all elements int the subset.
	/**	WARNING: subsets are not automatically changed accordingly.
	 *	After termination Subset-Indices and Subset-Infos/iterators are asynchronous.
	 *	Make sure to change subset-infos and iterators accordingly.*/		
		void change_subset_indices(int indOld, int indNew) override;

		
	///	add a subset
		void add_required_subset_lists(int maxIndex) override;
		
	///	erases the subset but does not touch the subset-indices.
		void erase_subset_lists(int index) override;

	///	swaps the subsets but does not touch the subset-indices.
		void swap_subset_lists(int ind1, int ind2) override;

	///	moves the subset but does not touch the subset-indices.
		void move_subset_lists(int indexFrom, int indexTo) override;

	///	join the subset-lists but do not touch the subset-indices.
		void join_subset_lists(int target, int src1, int src2) override;

	///	this method is called by ISubsetHandler when attachment_support has been enabled.
		//void register_subset_elements_at_pipe();

	////////////////////////////////////////////////
	//	protected helper methods
	///	a helper method for the public assign_subset methods.
		template <typename TElem>
		void assign_subset_impl(TElem* elem, int subsetIndex);

	///	helper for change_subset_indices
		template <typename TElem>
		void change_elem_subset_indices(int indOld, int indNew);
		
	///	Throws an error if the required level does not yet exist
		inline void level_required(int level) const;

		void add_level();
		void add_subset_to_all_levels();///< increases m_numSubsets.

	///	helper for collect_subset_elements
		// template <typename TElem>
		// size_t collect_subset_elements_impl(std::vector<TElem*>& elemsOut, int subsetIndex) const;
		
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

		using SubsetVec = std::vector<Subset*>;
		using LevelVec = std::vector<SubsetVec>;

	///	returns the subset with index si on the given level
		inline Subset* subset(int si, int level)	{return m_levels[level][si];}
		[[nodiscard]] inline const Subset* subset(int si, int level)	const {return m_levels[level][si];}

	///	creates a new subset. Caller is responsible for deletion
		[[nodiscard]] Subset* new_subset() const;

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
		get_list_iterator(Edge* o)
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
		template <typename TElem> inline
		typename Grid::traits<TElem>::SectionContainer&
		section_container(int si, int lvl);

	///	returns the const section container for the given type, subset and level
		template <typename TElem> inline
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

using MGSubsetHandler = MultiGridSubsetHandler;

/** \} */
}//	end of namespace

//	include implementation
#include "subset_handler_multi_grid_impl.hpp"

#endif
