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

#ifndef __H__UG_smooth_partition_bounds
#define __H__UG_smooth_partition_bounds

#include "../partitioner.h"
#include "lib_grid/iterators/lg_for_each.h"

namespace ug{

///	early draft. Currently only useful for prism-geometries in the d3f-wipp setting
template <class elem_t>/*, int dim>*/
class SmoothPartitionBounds : public IPartitionPostProcessor
{
	public:
		SmoothPartitionBounds() :
			m_mg(nullptr),
			m_partitions(nullptr)
			// m_consideringVerticalSidesOnly(false)
		{}

		~SmoothPartitionBounds() override {
			if(m_mg && m_mg->has_attachment<side_t>(m_aSubsetNbrs))
				m_mg->detach_from<side_t>(m_aSubsetNbrs);
		}

		void init_post_processing(MultiGrid* mg, SubsetHandler* partitions) override {
			m_mg = mg;
			m_partitions = partitions;
			mg->attach_to<side_t>(m_aSubsetNbrs);
		}

		// void set_position_attachment(position_attachment_t aPos);

		// void consider_vertical_sides_only(bool enable)
		// {
		// 	m_consideringVerticalSidesOnly = enable;
		// 	if(enable){
		// 		UG_COND_THROW(!mg->has_vertex_attachment(m_aPos),
		// 					  "consider_vertical_sides_only may only be enabled if "
		// 					  "a valid position attachment for the underlying grid "
		// 					  "has been set previously.");
		// 		m_aaPos.access(*mg, m_aPos);
		// 	}
		// }
		
		// bool considering_vertical_sides_only() const	{return m_consideringVerticalSidesOnly;}

		void post_process(int partitionLvl) override {
			using namespace std;
		//	we'll regularize the partition to reduce h-interface sizes
		//	(also important to improve gmg-smoother efficiency).
		//todo:	PARALLEL IMPLEMENTATION, EQUALLY DISTRIBUTE SWAP-ELEMENTS
			MultiGrid& mg = *m_mg;
			SubsetHandler& sh = *m_partitions;

			Grid::AttachmentAccessor<side_t, a_subset_pair_t>	aaSubsetNbrs(mg, m_aSubsetNbrs);

			typename MultiGrid::traits<side_t>::secure_container	sides;
			
			vector<int> nbrSubs;
			bool assignToSmallerSubsetsOnly = true;

			for(int main_iteration = 0; main_iteration < 2; ++main_iteration){
				SetAttachmentValues(aaSubsetNbrs, mg.begin<side_t>(partitionLvl),
									mg.end<side_t>(partitionLvl), pair<int, int>(-1, -1));

			//	regularization step 1: assign element-subset-indices to sides
				for(typename Grid::traits<elem_t>::iterator _feI = mg.begin<elem_t>(partitionLvl); _feI != mg.end<elem_t>(partitionLvl); ++_feI){ elem_t* e = *_feI;{
					int si = sh.get_subset_index(e);
					if(si >= 0){
						mg.associated_elements(sides, e);
						for(size_t _vfeI = 0; _vfeI < sides.size(); ++_vfeI){ side_t* s = sides[_vfeI];{
							if(aaSubsetNbrs[s].first == -1)
								aaSubsetNbrs[s].first = si;
							else
								aaSubsetNbrs[s].second = si;
						}};
					}
				}};

			//todo:	regularization step 1.5: communicate subsetNbrs
			//	...
			
			//	regularization step 2: check for elements whose neighbors belong to other partitions
				for(typename Grid::traits<elem_t>::iterator _feI = mg.begin<elem_t>(partitionLvl); _feI != mg.end<elem_t>(partitionLvl); ++_feI){ elem_t* e = *_feI;{
					int si = sh.get_subset_index(e);
					if(si >= 0){
						nbrSubs.clear();
						mg.associated_elements(sides, e);
						int numOwnSubsetNbrs = 0;
						for(size_t _vfeI = 0; _vfeI < sides.size(); ++_vfeI){ side_t* s = sides[_vfeI];{
							if(s->reference_object_id() == ROID_QUADRILATERAL){
								const subset_pair_t& nbrs = aaSubsetNbrs[s];
								if(nbrs.first == si && nbrs.second == si)
									++numOwnSubsetNbrs;
								else if(nbrs.first != -1 && nbrs.first != si)
									nbrSubs.push_back(nbrs.first);
								else if(nbrs.second != -1 && nbrs.second != si)
									nbrSubs.push_back(nbrs.second);
							}
						}};

					//	we'll now search for the subset which shares the most sides.
						int newSubsetNbrs = numOwnSubsetNbrs;
						int newSub = si;

						for(size_t i_nbrSub = 0; i_nbrSub < nbrSubs.size(); ++i_nbrSub){
							int nbrSub = nbrSubs[i_nbrSub];
							if(nbrSub != -1){// otherwise the subset has already been processed
								int count = 1;
								for(size_t i = i_nbrSub + 1; i < nbrSubs.size(); ++i){
									if(nbrSubs[i] == nbrSub){
										nbrSubs[i] = -1;
										++count;
									}
								}

								if(count > newSubsetNbrs){
									newSubsetNbrs = count;
									newSub = nbrSub;
								}
							}
						}

					//	assign the new subset.
					//todo:	use a more elaborate decider to which subset an element shall be assigned.
						if(newSub != -1
						   && ((assignToSmallerSubsetsOnly && newSub < si)
						   	   || (!assignToSmallerSubsetsOnly && newSub != si)))
						{
							sh.assign_subset(e, newSub);
						}
					}
				}};

			//	in the second iteration we'll assign to all neighbors
				assignToSmallerSubsetsOnly = false;
			}
		}

		void partitioning_done() override {
			m_mg->detach_from<side_t>(m_aSubsetNbrs);
		}

	private:
		using side_t = typename elem_t::side;
		using subset_pair_t = std::pair<int, int>;
		using a_subset_pair_t = Attachment<subset_pair_t>;
		
		MultiGrid*		m_mg;
		SubsetHandler*	m_partitions;
		a_subset_pair_t	m_aSubsetNbrs;
		// bool			m_consideringVerticalSidesOnly;
};

}//	end of namespace

#endif