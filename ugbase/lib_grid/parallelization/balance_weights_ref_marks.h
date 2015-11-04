#ifndef __H__UG__balance_weights_ref_marks__
#define __H__UG__balance_weights_ref_marks__

#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "load_balancer.h"

namespace ug{

class BalanceWeightsRefMarks : public IBalanceWeights
{
	public:
		BalanceWeightsRefMarks(IRefiner* refiner) : m_refiner(refiner)	{};

		virtual number get_refined_weight(Vertex* e){
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE || m == RM_ANISOTROPIC)
				return get_weight(e);
			return 0;
		}

		virtual number get_refined_weight(Edge* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE || m == RM_ANISOTROPIC)
				return 2. * get_weight(e);
			return 0;
		}

		virtual number get_refined_weight(Face* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE)
				return 4. * get_weight(e);
			if(m == RM_ANISOTROPIC)
				return 2. * get_weight(e);// ok - that isn't always true...
			return 0;
		}

		virtual number get_refined_weight(Volume* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			if(m == RM_REFINE)
				return 8. * get_weight(e);	// ok - that isn't always true...
			if(m == RM_ANISOTROPIC)
				return 2. * get_weight(e);	// yep - that isn't always true...
			return 0;
		}

		virtual bool has_level_offsets()		{return true;}

	///	Indicator in which level the specifed elements should be partitioned.
	/** This will only have effect if the given element does not have any children
	 * \{ */
		virtual bool consider_in_level_above(Vertex* e)	{return consider_in_level_above_impl(e);}
		virtual bool consider_in_level_above(Edge* e) 	{return consider_in_level_above_impl(e);}
		virtual bool consider_in_level_above(Face* e) 	{return consider_in_level_above_impl(e);}
		virtual bool consider_in_level_above(Volume* e)	{return consider_in_level_above_impl(e);}	
	/** \} */

	private:
		template <class TElem>
		bool consider_in_level_above_impl(TElem* e)
		{
			RefinementMark m = m_refiner->get_mark(e);
			return (m == RM_REFINE) || (m == RM_ANISOTROPIC);
		}

		IRefiner*	m_refiner;
};

}// end of namespace

#endif
