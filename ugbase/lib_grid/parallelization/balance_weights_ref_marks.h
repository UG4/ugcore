/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__balance_weights_ref_marks__
#define __H__UG__balance_weights_ref_marks__

#include "lib_grid/refinement/refiner_interface.h"
#include "partitioner.h"

namespace ug{

class BalanceWeightsRefMarks : public IBalanceWeights
{
	public:
		BalanceWeightsRefMarks(IRefiner* refiner) : m_refiner(refiner)	{};

		number get_refined_weight(Vertex* e) override {
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE || m == RM_ANISOTROPIC)
				return get_weight(e);
			return 0;
		}

		number get_refined_weight(Edge* e) override {
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE || m == RM_ANISOTROPIC)
				return 2. * get_weight(e);
			return 0;
		}

		number get_refined_weight(Face* e) override {
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE)
				return 4. * get_weight(e);
			if(m == RM_ANISOTROPIC)
				return 2. * get_weight(e);// ok - that isn't always true...
			return 0;
		}

		number get_refined_weight(Volume* e) override {
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE)
				return 8. * get_weight(e);	// ok - that isn't always true...
			if(m == RM_ANISOTROPIC)
				return 2. * get_weight(e);	// yep - that isn't always true...
			return 0;
		}

		bool has_level_offsets() override {return true;}

	///	Indicator in which level the specified elements should be partitioned.
	/** This will only have effect if the given element does not have any children
	 * \{ */
		bool consider_in_level_above(Vertex* e) override {return consider_in_level_above_impl(e);}
		bool consider_in_level_above(Edge* e) override {return consider_in_level_above_impl(e);}
		bool consider_in_level_above(Face* e) override {return consider_in_level_above_impl(e);}
		bool consider_in_level_above(Volume* e) override {return consider_in_level_above_impl(e);}
	/** \} */

	private:
		template <typename TElem>
		bool consider_in_level_above_impl(TElem* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			return (m == RM_REFINE) || (m == RM_ANISOTROPIC);
		}

		IRefiner*	m_refiner;
};

}// end of namespace

#endif
