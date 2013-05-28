// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 16.11.2012 (d,m,y)

#include "parallel_hnode_adjuster.h"
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/algorithms/debug_util.h"
#include "common/error.h"

namespace ug{

template <class TLayout>
class ComPol_BroadcastRefineMarks : public pcl::ICommunicationPolicy<TLayout>
{
	public:
		typedef TLayout								Layout;
		typedef typename Layout::Type				GeomObj;
		typedef typename Layout::Element			Element;
		typedef typename Layout::Interface			Interface;
		typedef typename Interface::const_iterator	InterfaceIter;

		ComPol_BroadcastRefineMarks(IRefiner& ref, byte consideredMarks)
			 :	m_ref(ref), m_consideredMarks(consideredMarks)
		{}

		virtual ~ComPol_BroadcastRefineMarks()	{}
		virtual int
		get_required_buffer_size(const Interface& interface)
		{return interface.size() * sizeof(byte);}

	///	writes writes the selection states of the interface entries
		virtual bool
		collect(ug::BinaryBuffer& buff, const Interface& interface)
		{
		//	write the entry indices of marked elements.
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				byte refMark = m_ref.get_mark(elem);
				buff.write((char*)&refMark, sizeof(byte));
			}

			return true;
		}

	///	reads marks from the given stream
		virtual bool
		extract(ug::BinaryBuffer& buff, const Interface& interface)
		{
			byte val;
			for(InterfaceIter iter = interface.begin();
				iter != interface.end(); ++iter)
			{
				Element elem = interface.get_element(iter);
				buff.read((char*)&val, sizeof(byte));

				val &= m_consideredMarks;

			//	check the current status and adjust the mark accordingly
				byte curVal = m_ref.get_mark(elem);

				if(val > curVal){
					if(val & RM_COARSEN)
						m_ref.mark(elem, RM_COARSEN);
					if(val & RM_ANISOTROPIC)
						m_ref.mark(elem, RM_ANISOTROPIC);
					if(val & RM_REFINE)
						m_ref.mark(elem, RM_REFINE);
				}
			}
			return true;
		}

		IRefiner& m_ref;
		byte m_consideredMarks;
};



template <class TStdVector>
static bool ContainsInterfaceElem(const TStdVector& elems,
								  DistributedGridManager& distGridMgr)
{
	for(size_t i = 0; i < elems.size(); ++i){
		if(distGridMgr.is_interface_element(elems[i]))
			return true;
	}
	return false;
}


void ParallelHNodeAdjuster::
ref_marks_changed(IRefiner& ref,
			   	  const std::vector<VertexBase*>& vrts,
			   	  const std::vector<EdgeBase*>& edges,
			   	  const std::vector<Face*>& faces,
			   	  const std::vector<Volume*>& vols)
{
	UG_DLOG(LIB_GRID, 1, "refMarkAdjuster-start: ParallelHNodeAdjuster::ref_marks_changed\n");
	UG_ASSERT(ref.grid(), "A refiner has to operate on a grid, before marks can be adjusted!");
	if(!ref.grid()){
		UG_DLOG(LIB_GRID, 1, "refMarkAdjuster-stop: ParallelHNodeAdjuster::ref_marks_changed\n");
		return;
	}
	
	Grid& grid = *ref.grid();
	if(!grid.is_parallel()){
		UG_DLOG(LIB_GRID, 1, "refMarkAdjuster-stop: ParallelHNodeAdjuster::ref_marks_changed\n");
		return;
	}

	DistributedGridManager& distGridMgr = *grid.distributed_grid_manager();
	GridLayoutMap& layoutMap = distGridMgr.grid_layout_map();

//	check whether new interface elements have been selected
	bool newInterfaceVrtsMarked = ContainsInterfaceElem(vrts, distGridMgr);
	bool newInterfaceEdgeMarked = ContainsInterfaceElem(edges, distGridMgr);
	bool newInterfaceFacesMarked = ContainsInterfaceElem(faces, distGridMgr);
	bool newInterfaceVolsMarked = ContainsInterfaceElem(vols, distGridMgr);

	bool newlyMarkedElems = newInterfaceVrtsMarked ||
							newInterfaceEdgeMarked ||
							newInterfaceFacesMarked ||
							newInterfaceVolsMarked;

	bool exchangeFlag = pcl::OneProcTrue(newlyMarkedElems);

	if(exchangeFlag){
		const byte consideredMarks = RM_REFINE | RM_ANISOTROPIC;
		ComPol_BroadcastRefineMarks<VertexLayout> compolRefVRT(ref, consideredMarks);
		ComPol_BroadcastRefineMarks<EdgeLayout> compolRefEDGE(ref, consideredMarks);
		ComPol_BroadcastRefineMarks<FaceLayout> compolRefFACE(ref, consideredMarks);

	//	send data SLAVE -> MASTER
		m_intfComVRT.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
									compolRefVRT);

		m_intfComEDGE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
									compolRefEDGE);

		m_intfComFACE.exchange_data(layoutMap, INT_H_SLAVE, INT_H_MASTER,
									compolRefFACE);

		m_intfComVRT.communicate();
		m_intfComEDGE.communicate();
		m_intfComFACE.communicate();

	//	and now MASTER -> SLAVE (the selection has been adjusted on the fly)
		m_intfComVRT.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
									compolRefVRT);

		m_intfComEDGE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
									compolRefEDGE);

		m_intfComFACE.exchange_data(layoutMap, INT_H_MASTER, INT_H_SLAVE,
									compolRefFACE);

		m_intfComVRT.communicate();
		m_intfComEDGE.communicate();
		m_intfComFACE.communicate();

		UG_DLOG(LIB_GRID, 1, "refMarkAdjuster-stop (force continue): ParallelHNodeAdjuster::ref_marks_changed\n");
	}

	UG_DLOG(LIB_GRID, 1, "refMarkAdjuster-stop: ParallelHNodeAdjuster::ref_marks_changed\n");
}
}// end of namespace
