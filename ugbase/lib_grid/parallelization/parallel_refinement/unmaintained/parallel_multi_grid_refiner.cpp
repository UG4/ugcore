#include <vector>
#include "parallel_multi_grid_refiner.h"
#include "pcl/pcl.h"

using namespace std;

namespace ug
{

/**
 * Uses Grid::mark during collect (begins in begin_layout_collection and
 * ends in end_layout_collection).
 *
 * The RefinementMarkDistributor thus only works for layout-collections.
 *
 * TODO: optimize the distribution process!
 */
template <class TLayout>
class RefinementMarkDistributor : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;
		
		RefinementMarkDistributor(Grid& grid,
								ParallelMultiGridRefiner& refiner,
								std::vector<Element>& vMarkedInterfaceElems)
			 :	m_grid(grid),
			 	m_refiner(refiner),
				m_vMarkedInterfaceElems(vMarkedInterfaceElems)
		{
		}
		
		virtual bool
		begin_layout_collection(const Layout* pLayout)
		{
		//	begin marking and mark all elements in m_vMarkedInterfaceElems
			m_grid.begin_marking();
			for(size_t i = 0; i < m_vMarkedInterfaceElems.size(); ++i)
				m_grid.mark(m_vMarkedInterfaceElems[i]);
				
			return true;
		}

	///	signals the end of a layout collection
	/**	the default implementation returns true and does nothing else.*/
		virtual bool
		end_layout_collection(const Layout*)
		{
			m_grid.end_marking();
			return true;
		}
		
	///	writes entries for marked interface elements
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			if(!m_vMarkedInterfaceElems.empty())
			{
				int counter = 0;
				for(InterfaceIter iter = interface.begin();
					iter != interface.end(); ++iter)
				{
					Element elem = interface.get_element(iter);
					if(m_grid.is_marked(elem)){
						buff.write((char*)&counter, sizeof(int));
						int rule = m_refiner.get_rule(elem);
						buff.write((char*)&rule, sizeof(int));
					}
					++counter;
				}
			}
		//	write an end-sign.
			int endMark = -1;
			buff.write((char*)&endMark, sizeof(int));
			return true;
		}
		
		virtual bool
		begin_layout_extraction(const Layout* pLayout)
		{
			m_vNewMarks.clear();
			m_vNewRules.clear();
			return true;
		}
		
	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	iterate through interface elements.
		//	if indices match then mark it.
			InterfaceIter iter = interface.begin();
			int counter = 0;
			int markIndex = -1;
			buff.read((char*)&markIndex, sizeof(int));
			while(markIndex != -1){
				if(markIndex == counter){
					Element e = interface.get_element(iter);
				//	the entry is marked on the other side of the interface
					m_vNewMarks.push_back(e);
				//	get the associated rule
					int rule;
					buff.read((char*)&rule, sizeof(int));
					m_vNewRules.push_back(rule);
				//	get next entry					
					buff.read((char*)&markIndex, sizeof(int));
				}
				++iter;
				++counter;
				
				if(counter > (int)interface.size()){
					UG_LOG("ERROR in RefinementMarkDistributor: invalid element index in extract.\n");
					return false;
				}
			}
			
			return true;
		}
		
		const vector<Element>& get_new_marks() const
		{
			return m_vNewMarks;
		}

		const vector<int>& get_new_rules() const
		{
			return m_vNewRules;
		}
		
	protected:
		Grid&						m_grid;
		ParallelMultiGridRefiner&	m_refiner;
		std::vector<Element>& 		m_vMarkedInterfaceElems;
		vector<Element>				m_vNewMarks;
		vector<int>				m_vNewRules;
};
/*
ParallelMultiGridRefiner::
ParallelMultiGridRefiner() : m_pLayoutMap(NULL)
{
}
*/
ParallelMultiGridRefiner::
ParallelMultiGridRefiner(DistributedGridManager& distGridMgr) :
	MultiGridRefiner(*distGridMgr.get_assigned_grid()),
	m_distGridMgr(distGridMgr),
	m_bCommunicationEnabled(true)
{
}

ParallelMultiGridRefiner::~ParallelMultiGridRefiner()
{
}

template <class TDistributor, class TCommunicator>
void ParallelMultiGridRefiner::
exchange_data(TDistributor& distributor,
				TCommunicator& communicator,
				std::vector<typename TDistributor::Element>* pvReceivedElemsOut)
{
	if(communication_is_enabled())
	{
		GridLayoutMap& layoutMap = m_distGridMgr.grid_layout_map();

	//	send data SLAVE -> MASTER
		communicator.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
										distributor);
		communicator.communicate();
		
	//	mark received elements
		mark_received_elements(distributor);
		
	//	store received elements in the given vector
		if(pvReceivedElemsOut)
		{
			const vector<typename TDistributor::Element>& vElems = distributor.get_new_marks();
			for(size_t i = 0; i < vElems.size(); ++i)
				pvReceivedElemsOut->push_back(vElems[i]);
		}
				
	//	send data MASTER -> SLAVE
		communicator.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
										distributor);
		communicator.communicate();
		
		mark_received_elements(distributor);
		
	//	store received elements in the given vector
		if(pvReceivedElemsOut)
		{
			const vector<typename TDistributor::Element>& vElems = distributor.get_new_marks();
			for(size_t i = 0; i < vElems.size(); ++i)
				pvReceivedElemsOut->push_back(vElems[i]);
		}
	}
}
/*
template <class TGeomObj, class TIterator>
void ParallelMultiGridRefiner::
mark_fixed_elements(TIterator iterBegin, TIterator iterEnd)
{
	if(m_distGridMgr.grid_layout_map().has_layout<TGeomObj>(INT_V_MASTER))
	{
		for(TIterator iter = iterBegin; iter != iterEnd; ++iter)
		{
			byte status = distGridMgr.get_status(*iter);
			if(( status & ES_V_MASTER)
				&! (status & ES_H_MASTER)
				&! (status & ES_H_SLAVE))
			{
				set_rule(*iter, RM_FIXED);
			}
		}
	}
}
*/
void ParallelMultiGridRefiner::
collect_objects_for_refine()
{
/*
Algorithm Layout:
mark elements that may not be refined as fixed

initial mark distribution
	- adjust_initial_selection
	- distribute marks for edges (and faces).
	- mark associated vertices.

closure selection
	- select_closure
	- distribute closure-marks for edges (and faces)
		* differentiate between refinement and copy-marks.
	- mark associated vertices
	
copy element selection
	- iterate from 0 to copy_range
		- select_copy_elements
		- distribute marks for copy-elements
		- mark associated vertices
		
'behind the scenes'
There are vectors for vertices, edges, faces and volumes, that store all elements
whose rule changed during the collect_objects_for_refie method (m_vNewlyMarked...).
This is realised via protected virtual set_rule methods in MultiGridRefiner.
Since those methods are only called during collect_objects_for_refine,
it is clear that those vectors always hold valid and existing elements.
*/

//	mark fixed elements
	//mark_fixed_elements<Vertex>(m_pMG->begin<Vertex>(), m_pMG->end<Vertex>());
	
//	The refinement distributor classes
	RefinementMarkDistributor<VertexLayout> vertexMarkDistributor(
												*m_distGridMgr.get_assigned_grid(),
												*this,
												m_vNewlyMarkedInterfaceVrts);
	
	RefinementMarkDistributor<EdgeLayout> edgeMarkDistributor(
												*m_distGridMgr.get_assigned_grid(),
												*this,
												m_vNewlyMarkedInterfaceEdges);
												
	pcl::InterfaceCommunicator<VertexLayout>	vertexCommunicator;	
	pcl::InterfaceCommunicator<EdgeLayout>	edgeCommunicator;
	
	//GridLayoutMap& layoutMap = m_distGridMgr.grid_layout_map();
	
//	element buffers
	vector<Vertex*> vVrts;
	
	//bool bCommunicate = true;
	//while(bCommunicate)
	{	
	////////////////////////////////
	//	initial selection
		adjust_initial_selection();
		
	//	exchange data
		exchange_data(vertexMarkDistributor, vertexCommunicator);
		exchange_data(edgeMarkDistributor, edgeCommunicator);
	}
	
////////////////////////////////
//	closure selection
	clear_newly_marked_element_buffers();
	select_closure(vVrts);
	
	exchange_data(vertexMarkDistributor, vertexCommunicator, &vVrts);
	exchange_data(edgeMarkDistributor, edgeCommunicator);
	
////////////////////////////////
//	collect copy-elements
	int iFirst = 0;
	for(int i = 0; i < get_copy_range(); ++i)
	{
		clear_newly_marked_element_buffers();
		int oldSize = vVrts.size();
		clear_newly_marked_element_buffers();
		select_copy_elements(vVrts, iFirst, 1);
		iFirst = oldSize;
		
		exchange_data(vertexMarkDistributor, vertexCommunicator, &vVrts);
		exchange_data(edgeMarkDistributor, edgeCommunicator);
	}

	clear_newly_marked_element_buffers();
}

void ParallelMultiGridRefiner::
refinement_step_begins()
{
	m_distGridMgr.begin_ordered_element_insertion();
}

void ParallelMultiGridRefiner::
refinement_step_ends()
{
	m_distGridMgr.end_ordered_element_insertion();
}

void ParallelMultiGridRefiner::
set_rule(Vertex* e, RefinementMark mark)
{
	MultiGridRefiner::set_rule(e, mark);
	if(m_distGridMgr.is_interface_element(e))
		m_vNewlyMarkedInterfaceVrts.push_back(e);
}

void ParallelMultiGridRefiner::
set_rule(Edge* e, RefinementMark mark)
{
	MultiGridRefiner::set_rule(e, mark);
	if(m_distGridMgr.is_interface_element(e))
		m_vNewlyMarkedInterfaceEdges.push_back(e);
}

void ParallelMultiGridRefiner::
set_rule(Face* e, RefinementMark mark)
{
	MultiGridRefiner::set_rule(e, mark);
	if(m_distGridMgr.is_interface_element(e))
		m_vNewlyMarkedInterfaceFaces.push_back(e);
}

void ParallelMultiGridRefiner::
set_rule(Volume* e, RefinementMark mark)
{
	MultiGridRefiner::set_rule(e, mark);
	if(m_distGridMgr.is_interface_element(e))
		m_vNewlyMarkedInterfaceVols.push_back(e);
}

template <class TMarkDistributor>
void ParallelMultiGridRefiner::
mark_received_elements(TMarkDistributor& distributor)
{
	const vector<typename TMarkDistributor::Element>& vElems =
									distributor.get_new_marks();
	const vector<int>& vRules = distributor.get_new_rules();
	
	for(size_t i = 0; i < vElems.size(); ++i){
		m_selMarks.select(vElems[i]);
		set_rule(vElems[i], (RefinementMark)vRules[i]);
	}
}

void ParallelMultiGridRefiner::
clear_newly_marked_element_buffers()
{
	m_vNewlyMarkedInterfaceVrts.clear();
	m_vNewlyMarkedInterfaceEdges.clear();
	m_vNewlyMarkedInterfaceFaces.clear();
	m_vNewlyMarkedInterfaceVols.clear();
}

void ParallelMultiGridRefiner::
adjust_parallel_selection(const std::vector<Vertex*>* pvVrts,
							const std::vector<Edge*>* pvEdges,
							const std::vector<Face*>* pvFaces,
							const std::vector<Volume*>* pvVols)
{
//	the grid on which we're operating
	Grid& grid = *m_distGridMgr.get_assigned_grid();
	
//	we'll use those vectors to collect associated elements
	vector<Edge*> vAssEdges;
	vector<Face*> vAssFaces;
	
//	in the moment only edges are supported
//	TODO: support all elements
	if(pvVrts)
		UG_ASSERT(pvVrts->empty(), "not yet supported");
	if(pvFaces)
		UG_ASSERT(pvFaces->empty(), "not yet supported");
	if(pvVols)
		UG_ASSERT(pvVols->empty(), "not yet supported");

//	select associated elements of edges
	if(pvEdges)
	{
		const vector<Edge*>& vEdges = *pvEdges;
		
		for(size_t i = 0; i < vEdges.size(); ++i){
			Edge* e = vEdges[i];
			
		//	make sure that associated vertices of the edge are selected
			for(size_t j = 0; j < 2; ++j){
				m_selMarks.select(e->vertex(j));
			}
			
		//	if the edge is a refine-edge, we have to mark associated
		//	faces and volumes. we have to mark their associated
		//	faces, edges and vertices as copy-elements - if they are
		//	not already selected.
//TODO:	add support for volumes
			
			int rule = get_rule(e);
			switch(rule){
				case RM_REFINE:{
				//	mark associated faces
					CollectFaces(vAssFaces, grid, e);
					
					for(size_t j = 0; j < vAssFaces.size(); ++j)
					{
						Face* f = vAssFaces[j];
						if(!m_selMarks.is_selected(f)){
						//	mark the face
							m_selMarks.select(f);
							
						//	mark associated edges as copy-elements
							CollectEdges(vAssEdges, grid, f);
							
							bool refineRegular = true;
							for(size_t k = 0; k < vAssEdges.size(); ++k){
								Edge* assEdge = vAssEdges[k];
								if(m_selMarks.is_selected(assEdge)){
									if(get_rule(assEdge) != RM_REFINE)
										refineRegular = false;
								}
								else{
									m_selMarks.select(assEdge);
									set_rule(assEdge, RM_COPY);
									refineRegular = false;
								}
							}
							
						//	set the refinement rule of the face
							if(refineRegular)
								set_rule(f, RM_REFINE);
							else
								set_rule(f, RM_IRREGULAR);
								
						//	make sure that all vertices are selected
							for(size_t k = 0; k < f->num_vertices(); ++k){
								if(!m_selMarks.is_selected(f->vertex(k)))
									m_selMarks.select(f->vertex(k));
							}
						}
					}
				}break;
				
				default: break;
			}
		}
	}
	
}

}//	end of namespace
