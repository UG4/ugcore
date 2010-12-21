
template<typename neighborstruct, typename matrix_type>
bool other_UpdateRating(size_t node, std::vector<neighborstruct> &PN, famg_nodes &nodes,
		std::vector<bool> &prolongation_calculated, cgraph &SymmNeighGraph,	FAMGInterpolationCalculator<matrix_type> &calculator)
{
	if(prolongation_calculated[node])
		return nodes.update_rating(node, PN);
	else
	{
#if 0
		int coarseNeighbors=0;
		for(cgraph::cRowIterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
		{
			if(rating[(*conn)].is_coarse())
				coarseNeighbors++;
		}
		int r = min(0, 2-coarseNeighbors);
		if(r == rating[node])
			return false;

		rating[node] = r;
		if(r == 2)
			return true;
#else
		int hasCoarseNeighbors=false;
		for(cgraph::cRowIterator conn = SymmNeighGraph.begin_row(node); conn != SymmNeighGraph.end_row(node); ++conn)
		{
			if(nodes[(*conn)].is_coarse())
			{
				hasCoarseNeighbors = true;
				break;
			}
		}
		if(!hasCoarseNeighbors)
			return false;

#endif
		FAMG_LOG(2, "node " << node << " has coarse neighbors, but rating is not calculated yet. calculating rating now\n");
		calculator.get_possible_parent_pairs(node, PN, nodes[node]);
		prolongation_calculated[node] = true;
		nodes.update_rating(node, PN);
		return true;


	}
}

template<typename neighborstruct, typename matrix_type>
void other_Update(size_t node, std::vector<std::vector<neighborstruct> > &possible_parents, famg_nodes &nodes,
		maxheap<famg_nodeinfo> &heap,
		std::vector<bool> &prolongation_calculated,	cgraph &SymmNeighGraph, FAMGInterpolationCalculator<matrix_type> &calculator)
{
	if(!nodes[node].is_valid_rating())
		return;
	if(other_UpdateRating(node, possible_parents[node], nodes, prolongation_calculated, SymmNeighGraph, calculator))
	{
		if(nodes[node].is_uninterpolateable())
		{
			if(heap.is_in(node))
			{
				FAMG_LOG(2, " remove node " << node << " from heap! ");
				heap.remove(node);
			}
			else
				FAMG_LOG(2, " node " << node << " uninterpolateable, but not in heap anyway! ");
		}
		else
		{

			if(heap.is_in(node))
			{
				heap.update(node);
				FAMG_LOG(2, " updated node " << node << " in heap! ");
			}
			else
			{
				heap.insert_item(node);
				FAMG_LOG(2, " inserted node " << node << " into heap!");
			}
		}
	}
}

void AddUnmarkedNeighbors(cgraph &SymmNeighGraph, size_t i, std::vector<bool> &mark, std::vector<size_t> &list)
{
	for(cgraph::cRowIterator conn = SymmNeighGraph.begin_row(i); conn != SymmNeighGraph.end_row(i); ++conn)
		if(mark[(*conn)] == false) list.push_back((*conn));
}



template<typename matrix_type, typename neighborstruct>
void other_coarsening(const matrix_type &A, cgraph &SymmNeighGraph, matrix_type &P,
		std::vector<std::vector<neighborstruct> > &possible_parents,
		famg_nodes &rating, maxheap<famg_nodeinfo> &heap, std::vector<int> &newIndex,
		size_t &iNrOfCoarse, size_t &unassigned)
{
	size_t N = A.num_rows();

	possible_parents.clear();
	possible_parents.resize(A.num_rows());

	for(size_t j=0; j<N; j++)
	{
		if(A.is_isolated(j))
		{
			rating[j].set_fine();
			unassigned--;
		}
	}

	size_t i;
	for(i=0; i<N; i++)
	{
		if(IsCloseToBoundary(A, i, 2) == false)
			break;
	}

	FAMG_LOG(2, "\nStarting with node " << i << "\n");

	std::vector<bool> prolongation_calculated;
	prolongation_calculated.resize(N, false);

	std::vector<bool> bvisited;
	bvisited.resize(N, false);
	std::vector<size_t> neighborsToUpdate;

	FAMGInterpolationCalculator<matrix_type> calculator(A);

	calculator.get_possible_parent_pairs(i, possible_parents[i], rating[i]);
	prolongation_calculated[i] = true;
	rating.update_rating(i, possible_parents[i]);

	while(1)
	{
		neighborstruct2 &n = possible_parents[i][0];

		FAMG_LOG(2, "\n\n\nSelect next node...\n");
		FAMG_LOG(2, "======================================\n");
		FAMG_LOG(2, "node " << i << " has rating " << rating[i] << ". now gets fine. parents: ");
		for(size_t j=0; j < n.parents.size(); j++)
			FAMG_LOG(2, n.parents[j].from << " ");

		// node i gets fine. update neighbors.
		rating[i].set_fine();
		bvisited[i] = true;
		unassigned--;

		for(size_t j=0; j < n.parents.size(); j++)
			bvisited[i] = true;

		AddUnmarkedNeighbors(SymmNeighGraph, i, bvisited, neighborsToUpdate);

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
				AddUnmarkedNeighbors(SymmNeighGraph, node, bvisited, neighborsToUpdate);
			}
			P(i, newIndex[node]) = n.parents[j].value;
		}

		IF_FAMG_LOG(4) print_vector(neighborsToUpdate, "neighborsToUpdate");

		// update neighbors
		for(size_t j=0; j<neighborsToUpdate.size(); j++)
		{
			other_Update(neighborsToUpdate[j], possible_parents, rating, heap, prolongation_calculated, SymmNeighGraph, calculator);
			if(rating[j].is_valid_rating())
				bvisited[j] = false;
		}
		neighborsToUpdate.clear();


		FAMG_LOG(2, "\n\nSearching for next node...\n");
		if(heap.height() == 0)
		{

			FAMG_LOG(2, "heap is empty!\n");
			break;
			for(i=0; i<N; i++)
			{
				if(A.num_connections(i) > 2 && rating[i].is_valid_rating())
					break;
			}
			if(i == N) break;
			FAMG_LOG(2, "\n\nRESTARTING WITH NODE " << i << "!!!\n\n");

			calculator.get_possible_parent_pairs(i, possible_parents[i], rating[i]);
			prolongation_calculated[i] = true;
			rating.update_rating(i, possible_parents[i]);
		}
		else
			while(1)
			{
				i = heap.get_max();
				FAMG_LOG(2, "node " << i << " has best rating");
				if(prolongation_calculated[i] == false)
				{
					FAMG_LOG(2, ", but prolongation not calculated. update.\n");
					calculator.get_possible_parent_pairs(i, possible_parents[i], rating[i]);
					prolongation_calculated[i] = true;
					rating.update_rating(i, possible_parents[i]);
					heap.update(i);
				}
				else
				{
					FAMG_LOG(2, ", and prolongation calculated, take this node.\n");
					i = heap.remove_max();
					break;
				}
			}
	}

	for(size_t i=0; i<N; i++)
	{
		if(rating[i].rating == 0.0)
			rating[i].set_uninterpolateable();
	}
}
