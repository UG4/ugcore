/*
 * edge_weighting_callbacks.h
 *
 *  Created on: 29.01.2013
 *      Author: mbreit
 */

#ifndef __H__UG__LIB_GRID__PARALLELIZATION__UTIL__EDGE_WEIGHTING_CALLBACKS__
#define __H__UG__LIB_GRID__PARALLELIZATION__UTIL__EDGE_WEIGHTING_CALLBACKS__


#include "lib_grid/multi_grid.h"
#include "lib_grid/tools/tools.h"
#include "common/error.h"
#include <vector>

namespace ug
{

/**
*	Base class for edge weighting callbacks.
**/

class EdgeWeighting
{
	public:
		EdgeWeighting() : m_sh(NULL) {};
		virtual ~EdgeWeighting() {};

		// virtual operator() should never be called
		// if it is, my best guess is that a boost::function was constructed using the copy of an object
		// derived from EdgeWeighting instead of a reference
		// DerivedFromEdgeWeighting d;
		// EdgeWeighting& e = d;
		// boost::function<int (int, int)> f = e; 				<-- BAD
		// boost::function<int (int, int)> f = boost::ref(e);	<-- GOOD
		virtual int operator() (EdgeBase* e1, EdgeBase* e2)
		{UG_WARNING("This method must never be called. Only derived implementations may be."); return 0;};
		virtual int operator() (Volume* v1, Volume* v2)
		{UG_WARNING("This method must never be called. Only derived implementations may be."); return 0;};
		virtual int operator() (Face* f1, Face* f2)
		{UG_WARNING("This method must never be called. Only derived implementations may be."); return 0;};

		void set_subset_handler(SmartPtr<MGSubsetHandler> sh)
		{
			m_sh = sh.get_nonconst();
		}

	protected:
		MGSubsetHandler* m_sh;
};


/**
 * First simple implementation: Preventing division between two distinct subsets.
**/

class InterSubsetEdgeWeighting : public EdgeWeighting
{
	public:
		InterSubsetEdgeWeighting() : EdgeWeighting(), m_vsi1(0), m_vsi2(0), m_vweights(0), m_hWeight(1), m_vWeight(1) {};
		virtual ~InterSubsetEdgeWeighting() {};

	public:
		void setIndivisibleBndBetweenSubsets(int si1, int si2, int weight)
		{
			m_vsi1.push_back(si1);
			m_vsi2.push_back(si2);
			m_vweights.push_back(weight);
		}

		void setDefaultWeights(int hWeight, int vWeight)
		{
			m_hWeight = hWeight;
			m_vWeight = vWeight;
		}

		virtual int operator() (EdgeBase* e1, EdgeBase* e2) {return weigh(e1,e2);};
		virtual int operator() (Face* f1, Face* f2) {return weigh(f1,f2);};
		virtual int operator() (Volume* v1, Volume* v2) {return weigh(v1,v2);};

	private:
		template <class TElem>
		int weigh(TElem* e1, TElem* e2)
		{
			if (this->m_sh)
			{
				if (this->m_sh->get_level(e1) == this->m_sh->get_level(e2))
				{
					// check whether elems fullfill one of the indivisibility conditions
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
			else UG_THROW("Subset handler must be assigned to InterSubsetEdgeWeighting before it is used!")
		}

	private:
		std::vector<int> m_vsi1;
		std::vector<int> m_vsi2;
		std::vector<int> m_vweights;
		int m_hWeight;
		int m_vWeight;

};


} //	end of namespace


#endif /* __H__UG__LIB_GRID__PARALLELIZATION__UTIL__EDGE_WEIGHTING_CALLBACKS__ */
