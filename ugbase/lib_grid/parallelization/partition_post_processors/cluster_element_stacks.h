/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_cluster_element_stacks
#define __H__UG_cluster_element_stacks

#include <stack>

#include "../partitioner.h"
// #include "lib_grid/iterators/lg_for_each.h"
#include "lib_grid/algorithms/normal_calculation.h"

namespace ug {

template <typename elem_t, typename vector_t>
class ClusterElementStacks : public IPartitionPostProcessor
{
	public:
		ClusterElementStacks()
		{
			VecSet(m_stackingDir, 0);
		}

		ClusterElementStacks(
				const Attachment<vector_t>& aPos,
				const vector_t& stackingDir) :
			m_aPos(aPos),
			m_stackingDir(stackingDir)
		{}

		~ClusterElementStacks() override = default;

		void set_position_attachment(const Attachment<vector_t>& aPos)
		{
			m_aPos = aPos;
		}

		void set_stacking_direction(const vector_t& stackingDir)
		{
			m_stackingDir = stackingDir;
		}

		void init_post_processing(MultiGrid* mg, SubsetHandler* partitions) override {
			m_mg = mg;
			m_partitions = partitions;
			mg->attach_to<elem_t>(m_aProcessed);
		}

		void post_process(int partitionLvl) override {
			using namespace std;

			UG_COND_THROW(VecLengthSq(m_stackingDir) == 0,
						  "Please specify a valid stacking direction. "
						  "Current stacking direction has length 0.");
			UG_COND_THROW(!m_mg, "no valid grid was specified during 'init_post_processing'");
			UG_COND_THROW(!m_partitions, "no valid partition map was specified during "
						  "'init_post_processing'");
			UG_COND_THROW(!m_mg->has_vertex_attachment(m_aPos), "Please specify "
						  "a valid position attachment (attached to the underlying grid) "
						  "before performing partitioning through "
						  "'ClusterElementStacks::set_position_attachment'.");

			MultiGrid& mg = *m_mg;
			SubsetHandler& sh = *m_partitions;

			Grid::VertexAttachmentAccessor<Attachment<vector_t> > aaPos(mg, m_aPos);
			Grid::AttachmentAccessor<elem_t, ABool> aaProcessed(mg, m_aProcessed);

		//	set all processed flags to false in this level
			SetAttachmentValues(aaProcessed, mg.begin<elem_t>(partitionLvl),
								mg.end<elem_t>(partitionLvl), false);

			vector_t stackingDir;
			VecNormalize(stackingDir, m_stackingDir);

		//	iterate over all elements. If one was already processed, ignore it.
		//	For each unprocessed element: Init the current stack with this element.
		//	Then traverse all elements which are neighbored to a stack element
		//	and which are connected by a side which is not orthogonal to the
		//	stacking direction. Add those elements to the stack.
		
			stack<elem_t*>	stack;
			typename Grid::traits<elem_t>::secure_container	elems;
			typename Grid::traits<side_t>::secure_container sides;

			for(auto _feI = mg.begin<elem_t>(partitionLvl); _feI != mg.end<elem_t>(partitionLvl); ++_feI){
				elem_t* rootElem = *_feI;
				if(aaProcessed[rootElem])
					continue;

				int newSi = sh.get_subset_index(rootElem);
				
				stack.push(rootElem);
				aaProcessed[rootElem] = true;

				while(!stack.empty()){
					elem_t* elem = stack.top();
					stack.pop();

					mg.associated_elements(sides, elem);
					for(size_t _vfeI = 0; _vfeI < sides.size(); ++_vfeI){
						side_t* s = sides[_vfeI];
						vector_t n = CalculateNormal(s, aaPos);
						if(fabs(VecDot(n, stackingDir)) > SMALL){
						//	cluster neighbors of that side
							mg.associated_elements(elems, s);
							for(size_t _vfeI = 0; _vfeI < elems.size(); ++_vfeI){
								elem_t* nbr = elems[_vfeI];
								if(!aaProcessed[nbr]){
									stack.push(nbr);
									aaProcessed[nbr] = true;
									sh.assign_subset(nbr, newSi);
								}	
							}
						}
					}
				}
			}
		}

		void partitioning_done() override {
			m_mg->detach_from<elem_t>(m_aProcessed);
		}

	private:
		using side_t = typename elem_t::side;
		using a_position_t = Attachment<vector_t>;

		MultiGrid*		m_mg;
		SubsetHandler*	m_partitions;
		ABool			m_aProcessed;
		a_position_t	m_aPos;
		vector_t		m_stackingDir;
};


}//	end of namespace

#endif