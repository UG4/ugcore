/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__load_balancer_util__
#define __H__UG__load_balancer_util__

#include "partitioner.h"
#include "lib_grid/algorithms/volume_calculation.h"

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
//		using base_class = ConnectionWeights<dim>;
//		using elem_type = typename base_class::elem_type;
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
		using position_attachment_t = Attachment<MathVector<dim> >;
		using elem_t = typename GeomObjBaseTypeByDim<dim>::base_obj_type;
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
