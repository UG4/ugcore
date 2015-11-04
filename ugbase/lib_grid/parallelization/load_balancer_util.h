#ifndef __H__UG__load_balancer_util__
#define __H__UG__load_balancer_util__

#include "load_balancer.h"
#include "lib_grid/algorithms/volume_calculation.h"
#include "lib_grid/lib_grid.h"

namespace ug{

///	Creates a process-hierarchy that fullfills the given conditions.
SPProcessHierarchy
CreateProcessHierarchy(size_t* numElemsOnLvl, size_t numLvls,
					   size_t minNumElemsPerProcPerLvl, size_t maxNumRedistProcs,
					   size_t maxNumProcs, int minDistLvl,
					   int maxLvlsWithoutRedist);

//
//template <int dim>
//class StdConnectionWeights : public ConnectionWeights<dim>{
//	public:
//		typedef ConnectionWeights<dim>			base_class;
//		typedef typename base_class::elem_type	elem_type;
//
//		StdConnectionWeights() : m_wgt(1.0)					{}
//		StdConnectionWeights(number wgt) : m_wgt(wgt)		{}
//		virtual ~StdConnectionWeights()						{}
//
//		virtual void set_weight(number wgt)					{m_wgt = wgt;}
//
//		virtual void set_grid(MultiGrid*, Attachment<MathVector<dim> >)	{}
//		virtual void refresh_weights(int)					{}
//		virtual number get_weight(elem_type*, elem_type*)	{return m_wgt;}
//
//	private:
//		number m_wgt;
//};


/**	If a level-factor > 0 is specified, then the get_weight method returns
 * for an element e with on level l:
 * \code
 * (l+1) * lvlFac * wgt.
 * \endcode
 * If level-factor <= 0, then m_wgt is simply returned.
 * The default level-factor is 0, the default weight is 1.
 */
class StdBalanceWeights : public IBalanceWeights{
	public:
		StdBalanceWeights() :
			m_wgt(1.0)				{}

		virtual ~StdBalanceWeights()						{}

		virtual void set_weight(number wgt)					{m_wgt = wgt;}
		virtual void refresh_weights(int)					{}

		virtual number get_weight(Vertex* e)	{return m_wgt;}
		virtual number get_weight(Edge* e)		{return m_wgt;}
		virtual number get_weight(Face* e)		{return m_wgt;}
		virtual number get_weight(Volume* e)	{return m_wgt;}

	private:
		number		m_wgt;
};


///	The higher the volume, the higher the weight when anisotropic refinement is used...
template <int dim>
class AnisotropicBalanceWeights : public IBalanceWeights{
	public:
		typedef Attachment<MathVector<dim> >	position_attachment_t;
		typedef typename GeomObjBaseTypeByDim<dim>::base_obj_type elem_t;
		AnisotropicBalanceWeights() : m_weightFactor(1)	{}
		virtual ~AnisotropicBalanceWeights()	{}

		virtual void set_weight_factor(number weightFactor)
		{
			m_weightFactor = weightFactor;
		}

		virtual number weight_factor() const	{return m_weightFactor;}

		virtual void set_grid(MultiGrid* mg, Attachment<MathVector<dim> > aPos)
		{
			m_aaPos.access(*mg, aPos);
		}

		virtual void refresh_weights(int baseLevel)	{}

		virtual number get_weight(Vertex* e)	{return 1;}
		virtual number get_weight(Edge* e)		{return get_weight_impl(e);}
		virtual number get_weight(Face* e)		{return get_weight_impl(e);}
		virtual number get_weight(Volume* e)	{return get_weight_impl(e);}

	private:
		number get_weight_impl(elem_t* e){
			return CalculateVolume(e, m_aaPos) * m_weightFactor;
		}

		number get_weight_impl(GridObject* e){
			return 1;
		}

		number m_weightFactor;
		Grid::VertexAttachmentAccessor<position_attachment_t>	m_aaPos;

};


}// end of namespace

#endif
