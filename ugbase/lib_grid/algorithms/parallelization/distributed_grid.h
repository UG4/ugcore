// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y09 m08 d17

#ifndef __H__LIB_GRID__DISTRIBUTED_GRID__
#define __H__LIB_GRID__DISTRIBUTED_GRID__

#include <map>
#include <vector>
#include "lib_grid/lg_base.h"
#include "parallel_grid_layout.h"

namespace libGrid
{
enum ElementStatus
{
	ES_NONE = 0,
	ES_SCHEDULED_FOR_INTERFACE = 1 << 1,
	ES_IN_INTERFACE = 1 << 2,
	ES_MASTER = 1 << 3,
	ES_SLAVE = 1 << 4
};
////////////////////////////////////////////////////////////////////////
///	Helps to create new interface-elements in the correct order.
/**
 * Between calls to begin_ordered_element_insertion()
 * and end_ordered_element_insertion(),
 * instances of this class will collect all created elements that have
 * been created from parent-elements that lie on an interface.
 * On end_element_creation() those elements will be added to interfaces
 * of the associated \sa GridCommunicationSet. The order in which those
 * elements are inserted is the same as the order that their parents
 * have in their interfaces.
 */
class DistributedGrid : public GridObserver
{		
	public:
		DistributedGridObserver();
		DistributedGridObserver(Grid& grid);
		
		virtual ~DistributedGridObserver();
		
	//	assignment
		void assign(Grid& grid);
			
	//	layout access
	/**	if you change the layout externally, be sure to call
	 *	DistributedGrid::layout_changed() afterwards.*/
		inline ParallelGridLayout& grid_layout()				{return m_gridLayout;}
		inline const ParallelGridLayout& grid_layout() const	{return m_gridLayout;}	
		
	///	call this method if you altered the layout externally.
	/**	This should be done as seldom as possible.
	 *	If you only added elements you may set addedElemsOnly to true.
	 *	The complexity in this case is proportional to the number of elements
	 *	in the layout.
	 *	If you removed elements or if you are unsure what operations have been
	 *	performed on the layout, you have to set addedElemsOnly to false
	 *	(the default value). Complexity in this case is proportial to the
	 *	number of elements in the underlying grid (or numer of elements in
	 *	the layout - whichever is higher).*/
	 	void grid_layout_changed(bool addedElemsOnly = false);
		
	//	element creation
	///	call this method before you start creating new elements in the associated grid.
	/** You shouldn't add new interfaces to the associated communication-set
	 *  between begin_ and end_element_creation.*/
		void begin_ordered_element_insertion();
		
	///	call this method when you're done with element creation.
	/**	Elements will not be added to the associated \sa GridCommunicationSet
	 *  until this method is called.*/
		void end_ordered_element_insertion();
		
	//	element-status
		inline bool check_status(VertexBase* vrt, byte status)
			{return ((m_aaStatusVRT[vrt] & status) == status);}
		
		inline bool check_status(EdgeBase* edge, byte status)
			{return ((m_aaStatusEDGE[edge] & status) == status);}

		inline bool check_status(Face* face, byte status)
			{return ((m_aaStatusFACE[face] & status) == status);}

		inline bool check_status(Volume* vol, byte status)
			{return ((m_aaStatusVOL[vol] & status) == status);}
			
		inline byte get_status(VertexBase* vrt)	{return m_aaStatusVRT[vrt];}
		inline byte get_status(EdgeBase* edge)	{return m_aaStatusEDGE[edge];}
		inline byte get_status(Face* face)		{return m_aaStatusFACE[face];}
		inline byte get_status(Volume* vol)		{return m_aaStatusVOL[vol];}
		
	//	grid callbacks
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);
		virtual void elements_to_be_cleared(Grid* grid);
		
	//	vertex callbacks
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL);
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt);
		virtual void vertex_to_be_replaced(Grid* grid, VertexBase* vrtOld, VertexBase* vrtNew);

	//	edge callbacks
		virtual void edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent = NULL);
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* edge);
		virtual void edge_to_be_replaced(Grid* grid, EdgeBase* edgeOld, EdgeBase* edgeNew);

	//	face callbacks
		virtual void face_created(Grid* grid, Face* face, GeometricObject* pParent = NULL)				{}
		virtual void face_to_be_erased(Grid* grid, Face* face)			{}
		virtual void face_to_be_replaced(Grid* grid, Face* faceOld, Face* faceNew)	{}

	//	volume callbacks
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL)			{}
		virtual void volume_to_be_erased(Grid* grid, Volume* vol)		{}
		virtual void volume_to_be_replaced(Grid* grid, Volume* volOld, Volume* volNew)	{}
		
	protected:
		
		class ScheduledElement
		{
			public:
				ScheduledElement()	{}
				ScheduledElement(GeometricObject* obj, InterfaceNodeTypes nt,
								int procID) : pObj(obj),
											nodeType(nt),
											connectedProcID(procID)	{}
				GeometricObject*	pObj;
				InterfaceNodeTypes	nodeType;
				int					connectedProcID;
		};
		
		typedef std::multimap<int, ScheduledElement> ScheduledElemMap;
/*
		typedef VertexCommunicationGroup::HNODE HVertex;
		typedef EdgeCommunicationGroup::HNODE HEdge;
		typedef FaceCommunicationGroup::HNODE HFace;
		typedef VolumeCommunicationGroup::HNODE HVolume;
*/		
		typedef util::Attachment<byte>	AStatus;

	protected:	
		void clear_scheduled_elements();
		void resize_scheduled_element_map_vecs();
		
		byte get_geom_obj_status(GeometricObject* go);
		
		inline void set_status(VertexBase* vrt, byte status)
			{m_aaStatusVRT[vrt] = status;}
		
		inline void set_status(EdgeBase* edge, byte status)
			{m_aaStatusEDGE[edge] = status;}

		inline void set_status(Face* face, byte status)
			{m_aaStatusFACE[face] = status;}

		inline void set_status(Volume* vol, byte status)
			{m_aaStatusVOL[vol] = status;}
			
		const std::vector<int>& get_interface_indices(GeometricObject* go);
		const std::vector<int>& get_interface_entry_indices(GeometricObject* go);
		
		template <class TParallelElemLayout, class TAttachmentAccessor>
		void set_elem_statuses(TParallelElemLayout& pel,
								TAttachmentAccessor& aaStatus, byte newStatus)

	///	vertex_created, edge_created, ... callbacks call this method.
		template <class TElem>
		void handle_created_element(TElem* pElem,
									GeometricObject* pParent);
		
	///	vertex_to_be_erased, edge..., ... callbacks call this method.
		template <class TElem, class TCommGrp>
		void handle_erased_element(TElem* e, TCommGrp& commGrp);
	
	///	vertex_to_be_replaced, edge..., ... callbacks call this method.	
		template <class TElem, class TCommGrp>
		void handle_replaced_element(TElem* eOld, TElem* eNew,
									TCommGrp& commGrp);
									
		template <class TScheduledElemMap>
		void perform_ordered_element_insertion(TScheduledElemMap& elemMap);
		
		template <class TElem, class TCommGrp>
		void add_element_to_interface(TElem* pElem, TCommGrp& commGrp,
										int procID, InterfaceNodeTypes nodeType);
						
	protected:
		Grid*					m_pGrid;
		ParallelGridLayout		m_gridLayout;
		
		ScheduledElemMap		m_vrtMap; ///< holds all elements that were scheduled by vertices
		ScheduledElemMap		m_edgeMap; ///< holds all elements that were scheduled by edges
		ScheduledElemMap		m_faceMap; ///< holds all elements that were scheduled by faces
		ScheduledElemMap		m_volMap; ///< holds all elements that were scheduled by volumes
		
		bool	m_bOrderedInsertionMode;
		AStatus	m_aStatus;

		Grid::VertexAttachmentAccessor<AStatus>	m_aaStatusVRT;
		Grid::EdgeAttachmentAccessor<AStatus>	m_aaStatusEDGE;
		Grid::FaceAttachmentAccessor<AStatus>	m_aaStatusFACE;
		Grid::VolumeAttachmentAccessor<AStatus>	m_aaStatusVOL;
};

}// end of namespace

#endif
