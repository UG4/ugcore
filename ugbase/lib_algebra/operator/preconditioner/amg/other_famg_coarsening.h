

template<typename neighborstruct>
bool other_UpdateRating(size_t node, std::vector<neighborstruct> &PN, std::vector<famg_nodeinfo> &nodes, std::vector<bool> &prolongation_calculated)
{
	if(prolognation_calculated[node])
		return UpdateRating(node, PN, nodes);
	else
	{
		int coarseNeighbors=0;
		for(cgraph::cRowIterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
		{
			if(rating[(*conn)].is_coarse())
				coarseNeighbors++;
		}
		int r = min(0, 2-coarseNeighbors);
		if(r != rating[node])
		{
			rating[node] = r;
			return true;
		}
	}
}


template<typename neighborstruct>
void other_UpdateNeighbors(const cgraph &SymmNeighGraph, size_t node, std::vector<std::vector<neighborstruct> > &possible_neighbors,
		std::vector<famg_nodeinfo> &nodes, maxheap<famg_nodeinfo> &heap, std::vector<bool> &prolongation_calculated)
{
	for(cgraph::cRowIterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
	{
		size_t neigh = (*conn);
		UG_ASSERT(node != neigh, "");
		if(!nodes[neigh].is_valid_rating())
			continue;
		if(other_UpdateRating(neigh, possible_neighbors[neigh], nodes, prolongation_calculated))
		{
			if(nodes[neigh].is_uninterpolateable())
			{
				FAMG_LOG(2, " remove from heap! ");
				if(heap.is_in(neigh)) heap.remove(neigh);
			}
			else
			{
				FAMG_LOG(2, " update in heap! ");
				if(heap.is_in(neigh)) heap.update(neigh);
				else heap.insert(neigh);
			}
		}
	}
}

void other_coarsening()
{
	famg_log_level = 4;
	size_t N = A.num_rows();

	size_t i;
	for(i=0; i<N; i++)
	{
		if(A.get_num_connections(i) > 2)
			break;
	}

	std::vector<bool> prolongation_calculated;
	prolongation_calculated.resize(N, false);

	FAMGParentNodesCalculator calculator(A);

	calculator.GetPossibleParentPairs(i, possible_neighbors[i], nodes[i]);
	prolongation_calculated[i] = true;

	while(1)
	{
		neighborstruct2 &n = possible_neighbors[i][0];

		FAMG_LOG(2, "\n\n\nSelect next node...\n");
		FAMG_LOG(2, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)
			FAMG_LOG(2, n.parents[j].from << " ");

		FAMG_LOG(2, "\nUpdate Neighbors of " << i << "\n");
		// node i gets fine. update neighbors.
		rating[i].set_fine();
		unassigned--;
		other_UpdateNeighbors(SymmNeighGraph, i, possible_neighbors, rating, heap, prolongation_calculated);

		FAMG_LOG(2, "Set coarse parents:\n");
		// get parent pair, set as coarse (if not already done), update neighbors.

		for(size_t j=0; j < n.parents.size(); j++)
		{
			size_t node = n.parents[j].from;

			if(rating[node].is_coarse()) { FAMG_LOG(2, "\nnode " << node << " is already coarse\n"); }
			else { FAMG_LOG(2, "\nnode " << node << " has rating " << rating[node] << ". now gets coarse.\nUpdate Neighbors of " << node << "\n"); }

			if(!rating[node].is_coarse())
			{
				if(heap.is_in(node)) heap.remove(node);
				rating[node].set_coarse();
				unassigned--;
				newIndex[node] = iNrOfCoarse++;
				P(node, newIndex[node]) = 1.0;
				other_UpdateNeighbors(SymmNeighGraph, node, possible_neighbors, rating, heap, prolongation_calculated);
			}
			P(i, newIndex[node]) = n.parents[j].value;
		}

		if(heap.height() == 0) break;
		while(1)
		{
			i = heap.get_max();
			if(prolongation_calculated[i] == false)
			{
				calculator.GetPossibleParentPairs(i, possible_neighbors[i], nodes[i]);
				prolongation_calculated[i] = true;
				heap.update(i);
			}
			else
			{
				i = heap.remove_max();
				break;
			}
		}
	}
}
