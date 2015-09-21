// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_cluster_element_stacks
#define __H__UG_cluster_element_stacks

#include <stack>
#include "../load_balancing.h"
#include "lib_grid/iterators/lg_for_each.h"
#include "lib_grid/algorithms/normal_calculation.h"

namespace ug{

template <class elem_t, class vector_t>
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

		virtual ~ClusterElementStacks()
		{}

		void set_position_attachment(const Attachment<vector_t>& aPos)
		{
			m_aPos = aPos;
		}

		void set_stacking_direction(const vector_t& stackingDir)
		{
			m_stackingDir = stackingDir;
		}

		void init_post_processing(MultiGrid* mg, SubsetHandler* partitions)
		{
			m_mg = mg;
			m_partitions = partitions;
			mg->attach_to<elem_t>(m_aProcessed);
		}

		void post_process(int partitionLvl)
		{
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

			lg_for_each_in_lvl_template(elem_t, rootElem, mg, partitionLvl){
				if(aaProcessed[rootElem])
					continue;

				int newSi = sh.get_subset_index(rootElem);
				
				stack.push(rootElem);
				aaProcessed[rootElem] = true;

				while(!stack.empty()){
					elem_t* elem = stack.top();
					stack.pop();

					mg.associated_elements(sides, elem);
					for_each_in_vec(side_t* s, sides){
						vector_t n = CalculateNormal(s, aaPos);
						if(fabs(VecDot(n, stackingDir)) > SMALL){
						//	cluster neighbors of that side
							mg.associated_elements(elems, s);
							for_each_in_vec(elem_t* nbr, elems){
								if(!aaProcessed[nbr]){
									stack.push(nbr);
									aaProcessed[nbr] = true;
									sh.assign_subset(nbr, newSi);
								}	
							}end_for;
						}
					}end_for;
				}
			}lg_end_for;
		}

		void partitioning_done()
		{
			m_mg->detach_from<elem_t>(m_aProcessed);
		}

	private:
		typedef typename elem_t::side		side_t;
		typedef Attachment<vector_t>		a_position_t;

		MultiGrid*		m_mg;
		SubsetHandler*	m_partitions;
		ABool			m_aProcessed;
		a_position_t	m_aPos;
		vector_t		m_stackingDir;
};


}//	end of namespace

#endif	//__H__UG_cluster_element_stacks
