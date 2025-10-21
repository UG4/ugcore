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

#ifndef __H__LIB_GRID__DISTRIBUTED_GRID__
#define __H__LIB_GRID__DISTRIBUTED_GRID__

#include <map>
#include <vector>
#include "parallel_grid_layout.h"
#include "lib_grid/multi_grid.h"
#include "common/util/owned_pointer.h"
#include "distro_adjuster.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization
/// @{

///	the states with which elements are marked in ug::DistributedGridManager
/**	Note that the constants are directly related to the constants enumerated
 * in InterfaceNodeTypes.
 * Please also note that the constants are currently stored in bytes. No
 * value using more than 8 bits is thus allowed.
 */
enum ElementStatusTypes
{
	ES_NONE = INT_NONE,
	ES_H_MASTER = INT_H_MASTER,
	ES_H_SLAVE = INT_H_SLAVE,
	ES_V_MASTER = INT_V_MASTER,
	ES_V_SLAVE = INT_V_SLAVE,

	//ES_GHOST = 1 << 5,//currently unused
	ES_SCHEDULED_FOR_INTERFACE = 1 << 6,
	ES_IN_INTERFACE = 1 << 7
};


///	manages the layouts and interfaces which are associated with a distributed grid.
/** The DistributedGridManager is a grid observer, which manages the GridLayoutMap
 * of a distributed grid. New elements are automatically added to interfaces as
 * required and erased elements are automatically removed from interfaces.
 *
 * Note that while a DistributedGridManager observes a grid, one may only
 * create new elements in the associated grid between calls to
 * begin_ordered_element_insertion() and end_ordered_element_insertion().
 *
 * Similarly elements of the associated grid may only be erased during calls to
 * begin_element_deletion() and end_element_deletion().
 *
 * Those barriers are important, so that distributed grid managers can insert
 * associated elements on different processes in the same order without communicating.
 *
 * Note that only one layer of elements may be created / erased between calls to
 * begin_ordered_element_insertion / end_ordered_... etc...
 */
class DistributedGridManager : public GridObserver
{	
	public:
		DistributedGridManager();
		DistributedGridManager(MultiGrid& grid);
		virtual ~DistributedGridManager();
		
	//	assignment
		void assign(MultiGrid& grid);
			
		inline MultiGrid* get_assigned_grid()	{return m_pGrid;}
		inline const MultiGrid* get_assigned_grid()	const {return m_pGrid;}

	//	layout access
	/**	if you change the layout externally, be sure to call
	 *	DistributedGrid::layout_changed() afterwards.*/
		inline GridLayoutMap& grid_layout_map()				{return m_gridLayoutMap;}
		inline const GridLayoutMap& grid_layout_map() const	{return m_gridLayoutMap;}	
		
	///	call this method if you altered the layout externally.
	/**	This should be done as seldom as possible.
	 *	If you only added elements you may set addedElemsOnly to true.
	 *	The complexity in this case is proportional to the number of elements
	 *	in the layout.
	 *	If you removed elements or if you are unsure what operations have been
	 *	performed on the layout, you have to set addedElemsOnly to false
	 *	(the default value). Complexity in this case is proportional to the
	 *	number of elements in the underlying grid (or numer of elements in
	 *	the layout - whichever is higher).*/
	 	void grid_layouts_changed(bool addedElemsOnly = false);
		
	///	returns the status of the given object.
	/**	The status gives information about the type of interfaces in which a node
	 * lies. The returned value is an or combination of the constants enumerated
	 * in InterfaceNodeTypes and ElementStatusTypes.
	 * \sa contains_status
	 * \{ */
		unsigned char get_status(GridObject* go) const;
		inline unsigned char get_status(Vertex* vrt)	const	{return elem_info(vrt).get_status();}
		inline unsigned char get_status(Edge* edge) const	{return elem_info(edge).get_status();}
		inline unsigned char get_status(Face* face) const		{return elem_info(face).get_status();}
		inline unsigned char get_status(Volume* vol) const		{return elem_info(vol).get_status();}
	/**	\} */

	///	returns true if the status of the given object contains the given status.
	/**	status can be an or-combination of constants enumerated in InterfaceNodeTypes
	 * and ElementStatusTypes.*/
		template <class TGeomObj>
		bool contains_status(TGeomObj* o, unsigned char status) const	{return (get_status(o) & status) == status;}

	///	returns true if the element is a ghost
	/**	ghost elements are vertical masters that are in no other interfaces.
	 *	Those elements shouldn't be refined.*/
	 	template<class TElem>
		inline bool is_ghost(TElem* elem) const;

	///	returns true if the element is contained in a horizontal interface
	 	template<class TElem>
		inline bool is_in_horizontal_interface(TElem* elem) const;

	///	returns true if the element is contained in a vertical interface
	 	template<class TElem>
		inline bool is_in_vertical_interface(TElem* elem) const;

	//	element creation
	///	call this method before you start creating new elements in the associated grid.
	/** You shouldn't add new interfaces to the associated communication-set
	 *  between begin_ and end_element_creation.*/
		void begin_ordered_element_insertion();
		
	///	call this method when you're done with element creation.
	/**	Elements will not be added to the associated \sa GridCommunicationSet
	 *  until this method is called.*/
		void end_ordered_element_insertion();
		
	//	element deletion
	///	call this method before you start deleting elements in the associated grid
		void begin_element_deletion();

	///	call this method after you're done deleting elements from the associated grid
		void end_element_deletion();

	///	returns true if an element is in one or more interfaces
		template <class TElem>
		bool is_interface_element(TElem* elem);
		
	/**	returns a list of pairs (procID, index) that tells for each element
	 *	where in which interfaces it lies.
	 *	\param	statusType	may be any or-combination of values enumerated in ElementStatusTypes.*/
	 	template <class TElem>
	 	void collect_interface_entries(
						std::vector<std::pair<int, size_t> >& vEntriesOut,
						TElem* elem, unsigned char statusType, bool clearContainer = true);


	///	Enables or disables interface managment. Use with care!
	/**	Interface managment is enabled by default. If you intend to completly
	 * restructure the grid and its interfaces, it may be beneficial to
	 * disable interface management before doing so. You should then use the method
	 * grid_layouts_changed to inform the DistributedGridManager that you
	 * modified the interfaces externally.*/
		void enable_interface_management(bool bEnable)	{m_interfaceManagementEnabled = bEnable;}


	/// set a global distribution adjuster
	/** A distribution adjuster may be used to manually adjust how elements are distributed
	 *  by the domain distribution. For example, in the circumstances that led to the implementation
	 *  of this feature, it was necessary to distribute certain vertices to all processes which held
	 *  elements of a specific subset (regardless of whether those processes had any elements
	 *  connected to these vertices or not).
	 *  During the distribution process, the adjuster's adjust() method is called at the end of
	 *  SelectElementsForTargetPartition() in distribution.cpp.
	 *
	 * @warning This is a highly experimental feature and as such not guaranteed to work properly.
	 */
		void set_distro_adjuster(SmartPtr<DistroAdjuster> adj) {m_spDistroAdjuster = adj;}
		
	/// get the distribution adjuster
		SmartPtr<DistroAdjuster> distro_adjuster() {return m_spDistroAdjuster;}

	////////////////////////////////
	//	grid callbacks
		virtual void grid_to_be_destroyed(Grid* grid);
		
	//	vertex callbacks
		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void edge_created(Grid* grid, Edge* e,
									GridObject* pParent = NULL,
									bool replacesParent = false);
		
		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void volume_created(Grid* grid, Volume* v,
									GridObject* pParent = NULL,
									bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, Vertex* vrt,
										 Vertex* replacedBy = NULL);

		virtual void edge_to_be_erased(Grid* grid, Edge* e,
										 Edge* replacedBy = NULL);

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL);

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL);

	protected:
	///	performs registration and deregistration at a grid.
	/**	call set_grid(NULL) to unregister the observer from a grid.*/
		void set_grid(Grid* grid);
		
	///	free's all grid related data
		void free_grid_data();

		template <class TGeomObj>
		void reset_elem_infos();
		
	/*	Currently unused. See implementation for more details.
		template <class TElem>
		void set_preliminary_ghost_states();

		void update_ghost_states();
	 */

		template <class TGeomObj, class TLayoutMap>
		void update_elem_info(TLayoutMap& layoutMap, int nodeType,
							  unsigned char newStatus, bool addStatus = false);

		template <class TGeomObj>
		void update_all_elem_infos();
							  
	///	vertex_created, edge_created, ... callbacks call this method.
		template <class TElem>
		void handle_created_element(TElem* pElem, GridObject* pParent,
									bool replacesParent);
		
		template <class TElem, class TScheduledElemMap, class TParent>
		void schedule_element_for_insertion(TScheduledElemMap& elemMap,
												TElem* elem,
												TParent* pParent);

		void clear_scheduled_elements();
		
		template <class TScheduledElemMap>
		void perform_ordered_element_insertion(TScheduledElemMap& elemMap);
		
		template <class TElem>
		void add_element_to_interface(TElem* pElem, int procID);

		template <class TElem>
		void element_to_be_erased(TElem* elem);

	/**	Note that the content of the given vector may be extended during this method.*/
		template <class TElem>
		void create_missing_constrained_h_interfaces(std::vector<TElem*>& newConstrainedElems);

	protected:
	///	Be careful when creating copies of ElementInfo.
	/**	Ownership of the internal data-object is transfered to the new copy.
	 * The old instance will thus point to a NULL pointer instead of the data
	 * object after the copy operation.
	 */
		template <class TGeomObj>
		class ElementInfo
		{
			public:
			//	types
				typedef typename GridLayoutMap::template Types<TGeomObj>
						::Interface		Interface;
				typedef typename Interface::iterator InterfaceElemIter;
				//typedef std::pair<Interface*, InterfaceElemIter> Entry;

				struct Entry{
					Entry()	{}
					Entry(Interface* intfc, InterfaceElemIter intfcElemIter, int intfcType) :
						m_interface(intfc), m_interfaceElemIter(intfcElemIter),
						m_interfaceType(intfcType)	{}

					Interface* 			m_interface;
					InterfaceElemIter	m_interfaceElemIter;
					int					m_interfaceType;
				};

				typedef std::list<Entry>				EntryList;
				typedef typename EntryList::iterator	EntryIterator;
				typedef typename EntryList::const_iterator	ConstEntryIterator;
				
			//	methods
				ElementInfo()	{}

				~ElementInfo()	{if(has_data()) m_data.reset();}

				void reset()
					{
						if(has_data()){
						//todo: reuse m_data
							m_data.reset();
						}
					}
				
				void add_entry(Interface* interface,
								InterfaceElemIter iter,
								int intfcType)				{data().m_entries.push_back(Entry(interface, iter, intfcType));}
				
				void remove_entry(Interface* interface)		{data().m_entries.erase(find_entry(interface));}
				
			///	Note: This method may only be called if is_interface_entry() returns true.
			/**	\{ */
				inline EntryIterator entries_begin()		{assert(has_data()); return m_data->m_entries.begin();}
				inline EntryIterator entries_end()			{assert(has_data()); return m_data->m_entries.end();}
				
				inline ConstEntryIterator entries_begin() const	{assert(has_data()); return m_data->m_entries.begin();}
				inline ConstEntryIterator entries_end() const	{assert(has_data()); return m_data->m_entries.end();}
			/**	\} */

				size_t get_local_id(EntryIterator iter) const	{return iter->m_interface->get_local_id(iter->m_interfaceElemIter);}
				size_t get_local_id(ConstEntryIterator iter) const	{return iter->m_interface->get_local_id(iter->m_interfaceElemIter);}
				int get_target_proc(EntryIterator iter) const	{return iter->m_interface->get_target_proc();}
				int get_target_proc(ConstEntryIterator iter) const	{return iter->m_interface->get_target_proc();}
				Interface* get_interface(EntryIterator iter)	{return iter->m_interface;}
				int get_interface_type(EntryIterator iter) const	{return iter->m_interfaceType;}
				int get_interface_type(ConstEntryIterator iter) const	{return iter->m_interfaceType;}

			///	Note: This method may only be called if is_interface_entry() returns true.
				EntryIterator find_entry(Interface* interface)
					{	assert(has_data());
						for(EntryIterator iter = entries_begin(); iter != entries_end(); ++iter){
							if(iter->m_interface == interface)
								return iter;
						}
						return entries_end();
					}
				
				void set_status(unsigned char status)
				{
					if(!has_data() && (status == ES_NONE))
						return;
					data().m_status = status;
				}
				unsigned char get_status() const
				{
					if(!has_data()) return ES_NONE;
					return m_data->m_status;
				}
				
				bool is_interface_element()
				{
					if(!has_data()) return false;
					return !m_data->m_entries.empty();
				}

			protected:
				struct Data{
					Data() : m_status(ES_NONE) {}
					EntryList	m_entries;
					unsigned char		m_status;
				};

			///	returns the data object. Creates it if necessary.
				inline Data& data()
				{
					if(!has_data())
						m_data.get() = new Data;
					return *m_data;
				}

				inline bool has_data() const	{return m_data.get() != NULL;}

			///	OwnedPtr is required to transfer ownership of the data-ptr during copy-operations.
			/**	Since ElementInfo objects are stored in attachments, they will be copied
			 * from time to time. Ownership is thereby transfered to the new copy.
			 */
				OwnedPtr<Data>	m_data;
		};
		
		typedef ElementInfo<Vertex>	ElemInfoVrt;
		typedef ElementInfo<Edge>	ElemInfoEdge;
		typedef ElementInfo<Face>		ElemInfoFace;
		typedef ElementInfo<Volume>		ElemInfoVol;
		
		typedef Attachment<ElemInfoVrt>		AElemInfoVrt;
		typedef Attachment<ElemInfoEdge>	AElemInfoEdge;
		typedef Attachment<ElemInfoFace>	AElemInfoFace;
		typedef Attachment<ElemInfoVol>		AElemInfoVol;
		
	///	Used to schedule an element for insertion during ordered-insertion-mode.
	/**	For each interface in which the parent lies, an instance of
	 *	this class is added to scheduledElementMap of the parents type.
	 *	The local ID of the parent will be used as key.
	 *	The type of layout into which the element shall be inserted can
	 *	be retrieved from the status of the associated geomObj.*/
		struct ScheduledElement
		{
			ScheduledElement(GridObject* obj, int procID) :
				geomObj(obj), connectedProcID(procID)				{}

			GridObject*	geomObj;
			int					connectedProcID;
		};
		
		typedef std::multimap<size_t, ScheduledElement>	ScheduledElemMap;

	protected:
		inline ElemInfoVrt& elem_info(Vertex* ele)	{return m_aaElemInfoVRT[ele];}
		inline ElemInfoEdge& elem_info(Edge* ele)	{return m_aaElemInfoEDGE[ele];}
		inline ElemInfoFace& elem_info(Face* ele)		{return m_aaElemInfoFACE[ele];}
		inline ElemInfoVol& elem_info(Volume* ele)		{return m_aaElemInfoVOL[ele];}

		inline const ElemInfoVrt& elem_info(Vertex* ele) const	{return m_aaElemInfoVRT[ele];}
		inline const ElemInfoEdge& elem_info(Edge* ele) const	{return m_aaElemInfoEDGE[ele];}
		inline const ElemInfoFace& elem_info(Face* ele) const		{return m_aaElemInfoFACE[ele];}
		inline const ElemInfoVol& elem_info(Volume* ele) const		{return m_aaElemInfoVOL[ele];}

		inline void got_new_constrained_vertical(Vertex* v)	{m_newConstrainedVerticalVrts.push_back(v);}
		inline void got_new_constrained_vertical(Edge* e)	{m_newConstrainedVerticalEdges.push_back(e);}
		inline void got_new_constrained_vertical(Face* f)		{m_newConstrainedVerticalFaces.push_back(f);}
		inline void got_new_constrained_vertical(Volume*)		{UG_THROW("There are no constrained volumes!");}


	protected:
		MultiGrid*		m_pGrid;
		GridLayoutMap	m_gridLayoutMap;
		
		bool m_interfaceManagementEnabled;///<only for debug purposes
		
		bool m_bOrderedInsertionMode;
		bool m_bElementDeletionMode;
		
		AElemInfoVrt	m_aElemInfoVrt;
		AElemInfoEdge	m_aElemInfoEdge;
		AElemInfoFace	m_aElemInfoFace;
		AElemInfoVol	m_aElemInfoVol;
		
		Grid::VertexAttachmentAccessor<AElemInfoVrt>	m_aaElemInfoVRT;
		Grid::EdgeAttachmentAccessor<AElemInfoEdge>		m_aaElemInfoEDGE;
		Grid::FaceAttachmentAccessor<AElemInfoFace>		m_aaElemInfoFACE;
		Grid::VolumeAttachmentAccessor<AElemInfoVol>	m_aaElemInfoVOL;
		
		ScheduledElemMap	m_vrtMap;	///< holds all elements that were scheduled by vertices
		ScheduledElemMap	m_edgeMap;	///< holds all elements that were scheduled by edges
		ScheduledElemMap	m_faceMap;	///< holds all elements that were scheduled by faces
		ScheduledElemMap	m_volMap;	///< holds all elements that were scheduled by volumes

		std::vector<Vertex*>	m_newConstrainedVerticalVrts;
		std::vector<Edge*>		m_newConstrainedVerticalEdges;
		std::vector<Face*>			m_newConstrainedVerticalFaces;

		SmartPtr<DistroAdjuster> m_spDistroAdjuster;
};

/// @}

}// end of namespace

////////////////////////////////
//	include implementation
#include "distributed_grid_impl.hpp"

#endif
