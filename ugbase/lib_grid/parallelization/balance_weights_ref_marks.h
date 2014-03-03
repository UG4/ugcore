// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Mar 3, 2014

#ifndef __H__UG__balance_weights_ref_marks__
#define __H__UG__balance_weights_ref_marks__

#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "load_balancer.h"

namespace ug{

class BalanceWeightsRefMarks : public IBalanceWeights
{
	public:
		BalanceWeightsRefMarks(IRefiner* refiner) : m_refiner(refiner)	{};

		virtual number get_weight(Vertex* e){
			return 1;
		}

		virtual number get_weight(Edge* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE || m == RM_ANISOTROPIC)
				return 2;
			return 1;
		}

		virtual number get_weight(Face* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE)
				return 4;
			if(m == RM_ANISOTROPIC)
				return 2;// ok - that isn't always true...
			return 1;
		}

		virtual number get_weight(Volume* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE)
				return 8;	// ok - that isn't always true...
			if(m == RM_ANISOTROPIC)
				return 2;	// yep - that isn't always true...
			return 1;
		}

	private:
		IRefiner*	m_refiner;
};

}// end of namespace

#endif
