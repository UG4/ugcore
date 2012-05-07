/**
 * \file agglomeration.cpp
 *
 * \author Martin Rupp
 *
 * \date 03.08.2011
 *
 * Goethe-Center for Scientific Computing 2011.
 */

#include "common/common.h"
#include <vector>
#include <map>
#include <algorithm>

namespace ug
{

struct supernode
{
	std::vector<int> components; // [0] is master
	int size;
	std::vector<int> connections;
	int father;

	supernode()
	{
		size = 0;
		father = -1;
	}

	void print()
	{
		UG_DLOG(LIB_ALG_AMG, 4, " size: " << size << "\n components: ")
		for(size_t i=0; i<components.size(); i++) UG_DLOG(LIB_ALG_AMG, 4, components[i] << " ");
		UG_DLOG(LIB_ALG_AMG, 4, "\n");
	}
	void print(std::vector<supernode> &supernodes)
	{
		print();
		UG_DLOG(LIB_ALG_AMG, 4, " connections to: ");
		for(size_t i=0; i<connections.size(); i++)
		{
			size_t j = connections[i];
			UG_ASSERT(supernodes[j].size != -1, j);
			UG_DLOG(LIB_ALG_AMG, 4, "{");
			for(size_t k=0; k<supernodes[j].components.size(); k++)
				UG_DLOG(LIB_ALG_AMG, 4, supernodes[j].components[k] << " ");
			UG_DLOG(LIB_ALG_AMG, 4, "} ");
		}
		UG_DLOG(LIB_ALG_AMG, 4, "\n");
	}
};


void EasyAgglomeration(const std::vector<size_t> sizes,
		const std::vector<std::map<int, size_t> > connections,
		std::vector<std::vector<int> > &mergeWith,
		size_t minimalSize, size_t preferredSize)
{
	// while(there are processors which have not enough elements)
	// get the smallest of them and merge it with a neighbor so that the least amount of new interface is produced

	preferredSize = std::max(minimalSize, preferredSize);

	std::vector<supernode> supernodes;
	supernodes.resize(sizes.size());
	for(size_t i=0; i<sizes.size(); i++)
	{
		supernodes[i].size = (int) sizes[i];
		supernodes[i].components.clear();
		supernodes[i].components.push_back(i);
		supernodes[i].connections.clear();
		for(std::map<int, size_t>::const_iterator it = connections[i].begin(); it != connections[i].end(); ++it)
			supernodes[i].connections.push_back(it->first);
		supernodes[i].father = -1;
	}

	size_t cnt;
	size_t nrAgglomerated=0;
	size_t smallest;
	do
	{
		UG_DLOG(LIB_ALG_AMG, 4, "========================\n");
		// get smallest supernode
		size_t ismallest = 0;
		smallest = 0;
		cnt = 0;
		for(size_t i=0; i<sizes.size(); i++)
		{
			if(supernodes[i].size == -1) continue; // node deleted
			UG_DLOG(LIB_ALG_AMG, 4, i << ":\n");
			supernodes[i].print(supernodes);

			if(smallest == 0 || supernodes[i].size < (int)smallest)
			{
				smallest = supernodes[i].size;
				ismallest = i;
			}
			cnt++;
		}


		if(cnt == 1) break;
		if(nrAgglomerated == 0 && smallest > minimalSize) break;
		else if(nrAgglomerated > 0 && smallest > preferredSize) break;

		// agglomerate with the smallest neighbor connection
		supernode &snSmallest = supernodes[ismallest];

		UG_DLOG(LIB_ALG_AMG, 4, "\nsmallest supernode is " << ismallest << "\n");
		supernodes[ismallest].print(supernodes);

		int iAgglo = snSmallest.connections[0];
		for(size_t i=1; i<snSmallest.connections.size(); i++)
		{
			UG_ASSERT(supernodes[snSmallest.connections[i]].size != -1, snSmallest.connections[i]);
			if(supernodes[snSmallest.connections[i]].size < supernodes[iAgglo].size)
				iAgglo = snSmallest.connections[i];
		}
		supernode &snAgglo = supernodes[iAgglo];

		UG_DLOG(LIB_ALG_AMG, 4, "\nagglomerating with node " << iAgglo << "\n");
		supernodes[iAgglo].print(supernodes);

		nrAgglomerated++;

		// put all components of snSmallest in snAgglo
		for(size_t i=0; i<snSmallest.components.size(); i++)
		{
			snAgglo.components.push_back(snSmallest.components[i]);
			supernodes[snSmallest.components[i]].father = iAgglo;
		}
		// increase size of snAgglo
		snAgglo.size += snSmallest.size;
		snSmallest.size = -1;
		snSmallest.father = iAgglo;

		// put all connections of snSmallest in snAgglo
		for(size_t i=0; i<snSmallest.connections.size(); i++)
		{
			int k = snSmallest.connections[i];
			if(k == iAgglo) continue;
			if(std::find(snAgglo.connections.begin(), snAgglo.connections.end(), k) == snAgglo.connections.end())
				snAgglo.connections.push_back(k);
		}


		// all nodes which are connected to snSmallest are now connected to snAgglo
		for(size_t i=0; i<snSmallest.connections.size(); i++)
		{
			int j = snSmallest.connections[i];
			supernode &s = supernodes[j];

			std::vector<int>::iterator it = std::find(s.connections.begin(), s.connections.end(), (int)ismallest);
			if(it != s.connections.end())
				s.connections.erase(it);

			if(j != iAgglo)
				if(std::find(s.connections.begin(), s.connections.end(), iAgglo) == s.connections.end())
					s.connections.push_back(iAgglo);
		}

		// delete snSmallest
		snSmallest.components.clear();
		snSmallest.connections.clear();

		// remove node ismallest

	} while (1);


	//UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 4);
	if(cnt==1)
		{ UG_DLOG(LIB_ALG_AMG, 4, "\nOnly one node left. Result\n"); }
	else
		{ UG_DLOG(LIB_ALG_AMG, 4, "\nSmallest supernode has size " << smallest << " < "
				<< preferredSize << " nodes. Result:\n"); }
	for(size_t i=0; i<sizes.size(); i++)
	{
		if(supernodes[i].size == -1) continue; // node deleted
		supernodes[i].print();
	}
	//UG_SET_DEBUG_LEVEL(LIB_ALG_AMG, 0);

	mergeWith.resize(sizes.size());
	for(size_t i=0; i<mergeWith.size(); i++)
	{
		mergeWith[i].clear();
		if(supernodes[i].father != -1)
			mergeWith[i].push_back(supernodes[i].father);
		else
		{
			UG_ASSERT(supernodes[i].components[0] == (int)i, i << " != " << supernodes[i].components[0]);
			mergeWith[i] = supernodes[i].components;
		}
	}
}

}
