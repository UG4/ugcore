/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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
					if(val & RM_CLOSURE)
						m_ref.mark(elem, RM_CLOSURE);
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
			   	  const std::vector<Vertex*>& vrts,
			   	  const std::vector<Edge*>& edges,
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
