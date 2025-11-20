/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef __H__UG__LIB_GRID__PARALLELIZATION__UTIL__PARTITION_WEIGHTING_CALLBACKS__
#define __H__UG__LIB_GRID__PARALLELIZATION__UTIL__PARTITION_WEIGHTING_CALLBACKS__


#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/tools.h"
#include "common/error.h"
#include <vector>

namespace ug
{

/**
*	Base class for partition weighting callbacks.
**/

class PartitionWeighting
{
	public:
		PartitionWeighting() : m_sh(nullptr), m_hWeight(1), m_vWeight(1) {};
		virtual ~PartitionWeighting() = default;

		virtual int operator () (Edge* e1, Edge* e2) {return weigh(e1,e2);};
		virtual int operator () (Face* f1, Face* f2) {return weigh(f1,f2);};
		virtual int operator () (Volume* v1, Volume* v2) {return weigh(v1,v2);};

		void set_subset_handler(MGSubsetHandler* sh)
		{
			m_sh = sh;
		}

		void set_default_weights(int hWeight, int vWeight)
		{
			m_hWeight = hWeight;
			m_vWeight = vWeight;
		}

	private:
		template <typename TElem>
		int weigh(TElem* e1, TElem* e2)
		{
			if (!this->m_sh)
				UG_THROW("Subset handler must be assigned to InterSubsetPartitionWeighting before it is used!");

			if (this->m_sh->get_level(e1) == this->m_sh->get_level(e2))
				return m_hWeight;

			return m_vWeight;
		}


	protected:
		MGSubsetHandler* m_sh;

		int m_hWeight;	// horizontal weight
		int m_vWeight;	// vertical weight
};


/**
 * First simple implementation: Preventing division between two distinct subsets.
**/

class InterSubsetPartitionWeighting : public PartitionWeighting
{
	public:
		InterSubsetPartitionWeighting() : PartitionWeighting(), m_vsi1(0), m_vsi2(0), m_vweights(0) {};
		virtual ~InterSubsetPartitionWeighting() = default;

	public:
		void set_inter_subset_weight(int si1, int si2, int weight)
		{
			m_vsi1.push_back(si1);
			m_vsi2.push_back(si2);
			m_vweights.push_back(weight);
		}

		int operator () (Edge* e1, Edge* e2) override {return weigh(e1,e2);};
		int operator () (Face* f1, Face* f2) override {return weigh(f1,f2);};
		int operator () (Volume* v1, Volume* v2) override {return weigh(v1,v2);};

	private:
		template <typename TElem>
		int weigh(TElem* e1, TElem* e2)
		{
			if (!this->m_sh)
				UG_THROW("Subset handler must be assigned to InterSubsetPartitionWeighting before it is used!");

			if (this->m_sh->get_level(e1) == this->m_sh->get_level(e2))
			{
				// check whether elems fulfill one of the indivisibility conditions
				for (size_t i = 0; i < m_vsi1.size(); i++)
				{
					if (this->m_sh->get_subset_index(e1) == m_vsi1[i]
						&& this->m_sh->get_subset_index(e2) == m_vsi2[i])
					{
						return m_vweights[i];
					}
				}
				return m_hWeight;
			}
			return m_vWeight;
		}

	private:
		std::vector<int> m_vsi1;
		std::vector<int> m_vsi2;
		std::vector<int> m_vweights;
};


/**
 * This PartitionWeighting sets out to protect specific subsets from having
 * vertices in any process boundaries.
 * To that purpose, in the weighing function, it checks both elements for
 * a common vertex in one of the specified subsets to protect. If such a vertex
 * exists, the weight for the division along the two elements is set to the
 * value defined via set_weight().
**/

class ProtectSubsetPartitionWeighting : public PartitionWeighting
{
	public:
		using vertex_list = MultiGrid::traits<Vertex>::secure_container;


	public:
		ProtectSubsetPartitionWeighting() : m_vSi(0), m_vWeights(0) {};
		~ProtectSubsetPartitionWeighting() override = default;

	public:
		void set_weight(int si, int weight)
		{
			m_vSi.push_back(si);
			m_vWeights.push_back(weight);
		}

		int operator () (Edge* e1, Edge* e2) override {return weigh(e1,e2);};
		int operator () (Face* f1, Face* f2) override {return weigh(f1,f2);};
		int operator () (Volume* v1, Volume* v2) override {return weigh(v1,v2);};

	private:
		template <typename TElem>
		int weigh(TElem* e1, TElem* e2)
		{
			if (!this->m_sh)
				UG_THROW("Subset handler must be assigned to InterSubsetPartitionWeighting before it is used!")

			if (this->m_sh->get_level(e1) == this->m_sh->get_level(e2))
			{
				// check whether elems fulfill one of the indivisibility conditions
				vertex_list vl1, vl2;
				this->m_sh->grid()->associated_elements(vl1, e1);
				this->m_sh->grid()->associated_elements(vl2, e2);

				for (size_t i = 0; i < vl1.size(); i++)
				{
					for (size_t j = 0; j < m_vSi.size(); j++)
					{
						// check if vertex is in one of the restricted subsets
						if (this->m_sh->get_subset_index(vl1[i]) != m_vSi[j])
							continue;

						// if so, check that it is a shared vertex of both elems
						for (size_t k = 0; k < vl2.size(); k++)
						{
							if (vl1[i] == vl2[k])
								return m_vWeights[j];
						}
					}
				}
				return m_hWeight;
			}
			return m_vWeight;
		}

	private:
		std::vector<int> m_vSi;
		std::vector<int> m_vWeights;
};


} //	end of namespace


#endif