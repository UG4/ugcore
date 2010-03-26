//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m03 d23

#include <vector>
#include "parallel_multi_grid_refiner.h"
#include "pcl/pcl.h"

using namespace std;

namespace ug
{

template <class TLayout>
class RefinementMarkDistributor : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout							Layout;
		typedef typename Layout::Type			GeomObj;
		typedef typename Layout::Element		Element;
		typedef typename Layout::Interface		Interface;
		typedef typename Interface::iterator	InterfaceIter;
		
		RefinementMarkDistributor(Selector& selMarks) : m_pSelMarks(&selMarks)
		{
		}
		
	///	writes entries for marked interface elements
		virtual bool
		collect(std::ostream& buff, Interface& interface)
		{
		//	write the entry indices of marked elements.
			int counter = 0;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				if(m_pSelMarks->is_selected(interface.get_element(iter))){
					buff.write((char*)&counter, sizeof(int));
					++counter;
				}
			}
			
		//	write an end-sign.
			int endMark = -1;
			buff.write((char*)&endMark, sizeof(int));
		}
		
	///	reads marks from the given stream
		virtual bool
		extract(std::istream& buff, Interface& interface)
		{
		//	iterate through interface elements.
		//	if an indices match then mark it.
			InterfaceIter iter = interface.begin();
			int counter = 0;
			int markIndex = -1;
			buff.read((char*)&markIndex, sizeof(int));
			while(markIndex != -1){
				if(markIndex == counter){
					Element e = interface.get_element(iter);
				//	the entry is marked on the other side of the interface
					if(!m_pSelMarks->is_selected(e)){
						m_pSelMarks->select(e);
						m_vNewMarks.push_back(e);
					}
				//	get next entry
					buff.read((char*)&markIndex, sizeof(int));
				}
				++iter;
				++counter;
			}
			
		}
		
		const vector<Element>& get_new_marks()
		{
			return m_vNewMarks;
		}
		
		void clear_new_marks()
		{
			m_vNewMarks.clear();
		}
		
	protected:
		Selector* m_pSelMarks;
		vector<Element>	m_vNewMarks;
};
/*
ParallelMultiGridRefiner::
ParallelMultiGridRefiner() : m_pLayoutMap(NULL)
{
}
*/
ParallelMultiGridRefiner::
ParallelMultiGridRefiner(DistributedGridManager& distGridMgr) :
	m_distGridMgr(distGridMgr),
	MultiGridRefiner(*distGridMgr.get_assigned_grid())
{
}

ParallelMultiGridRefiner::~ParallelMultiGridRefiner()
{
}

void ParallelMultiGridRefiner::
collect_objects_for_refine()
{
	RefinementMarkDistributor<EdgeLayout>	edgeMarkDistributor(m_selMarks);
	pcl::ParallelCommunicator<EdgeLayout>	edgeCommunicator;
	
	GridLayoutMap& layoutMap = m_distGridMgr.grid_layout_map();
	
//	vertices stored in here are used during copy-element-selection.
	vector<VertexBase*> vVrts;
	
	//bool bCommunicate = true;
	//while(bCommunicate)
	{	
	////////////////////////////////
	//	initial selection
		adjust_initial_selection();
		
	//	clear the marks-dump in the distributor
		edgeMarkDistributor.clear_new_marks();
		
	//	exchange data
		edgeCommunicator.exchange_data(layoutMap, INT_SLAVE, INT_MASTER,
										edgeMarkDistributor);
		edgeCommunicator.communicate();
		
		edgeCommunicator.exchange_data(layoutMap, INT_MASTER, INT_SLAVE,
										edgeMarkDistributor);
		edgeCommunicator.communicate();
		
	//	readjust initial selection
		adjust_initial_selection();
	}
	
////////////////////////////////
//	closure selection
	//while(!m_vNewlyMarkedEdges.empty())
	{
		m_vNewlyMarkedEdges.clear();
		m_vNewlyMarkedVertices.clear();
		
		select_closure(vVrts);
/*
	//	clear the marks-dump in the distributor
		edgeMarkDistributor.clear_new_marks();
		
	//	exchange data
		edgeCommunicator.exchange_data(layoutMap, INT_SLAVE, INT_MASTER,
										edgeMarkDistributor);
		edgeCommunicator.communicate();
		
		edgeCommunicator.exchange_data(layoutMap, INT_MASTER, INT_SLAVE,
										edgeMarkDistributor);
		edgeCommunicator.communicate();

		//adjust_initial_selection();
		
		for(size_t i = 0; i < m_vNewlyMarkedEdges.size(); ++i)
			for(size_t j = 0; j < 2; ++j)
				mark_for_refinement(m_vNewlyMarkedEdges[i]->vertex(j));

		for(size_t i = 0; i < m_vNewlyMarkedVertices.size(); ++i)
			vVrts.push_back(m_vNewlyMarkedVertices[i]);
			
		m_vNewlyMarkedVertices.clear();
	*/
	}

////////////////////////////////
//	collect copy-elements
	select_copy_elements(vVrts);

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
element_marked(VertexBase* elem)
{
	m_vNewlyMarkedVertices.push_back(elem);
}

void ParallelMultiGridRefiner::
element_marked(EdgeBase* elem)
{
	m_vNewlyMarkedEdges.push_back(elem);
}

void ParallelMultiGridRefiner::
element_marked(Face* elem)
{
}

void ParallelMultiGridRefiner::
element_marked(Volume* elem)
{
}

}//	end of namespace
