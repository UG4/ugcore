// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 25, 2013 (d,m,y)

#ifndef __H__UG__load_balancer_util__
#define __H__UG__load_balancer_util__

#include "load_balancer.h"

namespace ug{

///	Creates a process-hierarchy that fullfills the given conditions.
SPProcessHierarchy
CreateProcessHierarchy(size_t* numElemsOnLvl, size_t numLvls,
					   size_t minNumElemsPerProcPerLvl, size_t maxNumRedistProcs,
					   size_t maxNumProcs, int maxLvlsWithoutRedist);


template <int dim>
class StdConnectionWeights : public ConnectionWeights<dim>{
	public:
		typedef ConnectionWeights<dim>			base_class;
		typedef typename base_class::elem_type	elem_type;

		StdConnectionWeights() : m_wgt(1.0)					{}
		StdConnectionWeights(number wgt) : m_wgt(wgt)		{}
		virtual ~StdConnectionWeights()						{}

		virtual void set_weight(number wgt)					{m_wgt = wgt;}

		virtual void set_grid(MultiGrid*, Attachment<MathVector<dim> >)	{}
		virtual void refresh_weights(int)					{}
		virtual number get_weight(elem_type*, elem_type*)	{return m_wgt;}

	private:
		number m_wgt;
};


/**	If a level-factor > 0 is specified, then the get_weight method returns
 * for an element e with on level l:
 * \code
 * (l+1) * lvlFac * wgt.
 * \endcode
 * If level-factor <= 0, then m_wgt is simply returned.
 * The default level-factor is 0, the default weight is 1.
 */
template <int dim>
class StdBalanceWeights : public BalanceWeights<dim>{
	public:
		typedef BalanceWeights<dim>				base_class;
		typedef typename base_class::elem_type	elem_type;

		StdBalanceWeights() :
			m_mg(NULL), m_wgt(1.0), m_lvlFac(0)				{}
		StdBalanceWeights(number wgt, number lvlFac) :
			m_mg(NULL), m_wgt(wgt), m_lvlFac(lvlFac)		{}

		virtual ~StdBalanceWeights()						{}

		virtual void set_weight(number wgt)					{m_wgt = wgt;}
		virtual void set_level_fac(number fac)				{m_lvlFac = fac;}

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> >)	{m_mg = mg;}
		virtual void refresh_weights(int)					{}
		virtual number get_weight(elem_type* e)
		{
			if(m_lvlFac > 0){
				UG_ASSERT(m_mg, "A mult-grid has to be set before get_weight is called!");
				return number(m_mg->get_level(e) + 1) * m_lvlFac * m_wgt;
			}
			else
				return m_wgt;
		}

	private:
		MultiGrid*	m_mg;
		number		m_wgt;
		number		m_lvlFac;
};


}// end of namespace

#endif
