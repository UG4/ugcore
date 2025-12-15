/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__SUBDIVISION_RULES_PIECEWISE_LOOP__
#define __H__UG__SUBDIVISION_RULES_PIECEWISE_LOOP__

#include <vector>
#include <cassert>

#include "lib_grid/lg_base.h"

namespace ug {

///	\addtogroup lib_grid_algorithms_refinement_subdivision
///	@{

///	A singleton that stores all rules for a piecewise-loop subdivision surface.
class SubdivRules_PLoop
{
	public:
		struct NeighborInfo{
			NeighborInfo() = default;
			NeighborInfo(Vertex* n, size_t cval) :
				nbr(n), creaseValence(cval)	{}
				
			Vertex* nbr;
		/**	0 means that the neighbor is not a crease vertex.
		 *	> 0: The valence of the crease regarding only the
		 *	part on the side of the center-vertex.*/
			size_t creaseValence;
		};
		
	public:
	///	returns the only instance to this singleton.
		static SubdivRules_PLoop& inst()
		{
			static SubdivRules_PLoop subdivRules;
			return subdivRules;
		}

	////////////////////////////////
	//	WEIGHTS
		number ref_even_inner_center_weight(size_t valence) const
			{return 1. - static_cast<number>(valence) * get_beta(valence);}
		
		number ref_even_inner_nbr_weight(size_t valence) const
			{return get_beta(valence);}
		
	///	returns weights for center, nbr1 and nbr2.
		vector3 ref_even_crease_weights() const
			{return vector3(0.75, 0.125, 0.125);}
			
		vector4 ref_odd_inner_weights() const
			{return vector4(0.375, 0.375, 0.125, 0.125);}
		
	///	weights of an odd vertex on an inner edge that is connected to a crease.
	/**	The weight for the vertex on the crease is in v.x(), the inner edge vertex
	 *	in v.y() and the two indirectly connected vertex weights are in v.z() and v.w.
	 *	creaseValence specifies the number of associated edges of the crease vertex.
	 *
	 *	Rules are taken from:
	 *	"Piecewise Smooth Subdivision Surfaces with Normal Control",
	 *	H. Biermann, Adi Levin, Denis Zorin.*/
		vector4 ref_odd_inner_weights(size_t creaseValence) const
		{
			assert(creaseValence > 2 && "Bad crease valence. Underlying grid is not a surface grid.");
			if(creaseValence == 4)
				return ref_odd_inner_weights();
			number gamma = 0.5 - 0.25 * cos(PI / static_cast<number>(creaseValence - 1));
			return vector4(0.75 - gamma, gamma, 0.125, 0.125);
			
			//number c = 0.25 * cos((2.*PI) / (number)(creaseValence - 1));
			//return vector4(0.5 - c, 0.25 + c, 0.125, 0.125);
		}
		
		vector2 ref_odd_crease_weights() const
			{return vector2(0.5, 0.5);}

		number proj_inner_center_weight(size_t valence) const
			{return 1.0 - static_cast<number>(valence) / (0.375 / get_beta(valence) + valence);}
		
		number proj_inner_nbr_weight(size_t valence) const
			{return 1.0 / (0.375 / get_beta(valence) + valence);}
			
		vector3 proj_crease_weights() const
			{return vector3(2./3., 1./6., 1./6.);}

	/**	nbrInfos have to be specified in clockwise or counter-clockwise order.*/
		void proj_inner_crease_nbr_weights(number& centerWgtOut, number* nbrWgtsOut,
										   NeighborInfo* nbrInfos, size_t numNbrs) const
		{
			number wcntrProj = proj_inner_center_weight(numNbrs);
			number wnbrProj = proj_inner_nbr_weight(numNbrs);
			
		//	initialize all weights with 0
			centerWgtOut = 0;
			for(size_t i = 0; i < numNbrs; ++i)
				nbrWgtsOut[i] = 0;
				
			for(size_t i = 0; i < numNbrs; ++i){
				NeighborInfo& nbrInfo = nbrInfos[i];
				vector4 oddWeights(0.5, 0.5, 0, 0);
				
				if(nbrInfo.creaseValence == 0){
				//	Calculate the weights for odd-refinement on the given edge
					oddWeights = ref_odd_inner_weights();
				}
				else{
					oddWeights = ref_odd_inner_weights(nbrInfo.creaseValence);
				}
				
				nbrWgtsOut[i] += oddWeights.x() * wnbrProj;
				centerWgtOut += oddWeights.y() * wnbrProj;
				nbrWgtsOut[next_ind(i, numNbrs)] += oddWeights.z() * wnbrProj;
				nbrWgtsOut[prev_ind(i, numNbrs)] += oddWeights.w() * wnbrProj;
			}
			
			
		//	add scaled weights for evem refinement mask
			centerWgtOut += wcntrProj * ref_even_inner_center_weight(numNbrs);
			for(size_t i = 0; i < numNbrs; ++i)
				nbrWgtsOut[i] += wcntrProj * ref_even_inner_nbr_weight(numNbrs);
		}
/*			
		number proj_inner_crease_nbr_center_weight(size_t valence);
		number proj_inner_crease_nbr_nbr_weight(size_t valence);
		number proj_inner_crease_nbr_nbr_weight(size_t valence, size_t cValence);
*/
	///	returns beta as it is used in the subdivision masks.
	/**	performs a lookup if the valency is small enough.
	 *	calculates a fresh beta else.*/
		number get_beta(size_t valency) const;
/*
	////////////////////////////////
	//	EVEN MASKS
		template <typename TAAPos>
		typename TAAPos::ValueType
		apply_even_mask(Grid& grid, Vertex* center,
						TAAPos& aaPos);

		template <typename TAAPos>
		typename TAAPos::ValueType
		apply_even_crease_mask(Vertex* center, Vertex* nbr1,
							   Vertex* nbr2, TAAPos& aaPos);

	////////////////////////////////
	//	ODD MASKS
		template <typename TAAPos>
		typename TAAPos::ValueType
		apply_odd_mask(Vertex* vrt, Edge* parent,
					   TAAPos& aaPos);

		template <typename TAAPos>
		typename TAAPos::ValueType
		apply_odd_crease_mask(Vertex* vrt, Edge* parent,
							  TAAPos& aaPos);

		template <typename TAAPos>
		typename TAAPos::ValueType
		apply_odd_crease_nbr_mask(Grid& grid, Vertex* vrt,
								  Edge* parent, TAAPos& aaPos);

	////////////////////////////////
	//	PROJECTION
		template <typename TAAPos>
		typename TAAPos::ValueType
		project_inner_vertex(Grid& grid, Vertex* vrt,
							 TAAPos& aaPos);

		template <typename TAAPos>
		typename TAAPos::ValueType
		project_inner_vertex(Vertex* vrt, Vertex* nbrs,
							 int* nbrCreaseValencies, int numNbrs,
							 TAAPos& aaPos);

		template <typename TAAPos>
		typename TAAPos::ValueType
		project_crease_vertex(Vertex* vrt, Vertex* nbr1,
							  Vertex* nbr2, TAAPos& aaPos);
*/
	private:
	///	calculates beta as it is used in the subdivision masks.
		number calculate_beta(size_t valency) const;
		
	private:
	///	private constructor prohibits multiple instantiation.
		SubdivRules_PLoop();
		
	///	private copy constructor prohibits copy-construction.
		SubdivRules_PLoop(const SubdivRules_PLoop& src);
		
	///	private assignment operator prohibits assignment.
		SubdivRules_PLoop& operator = (const SubdivRules_PLoop& src);
		
	///	returns the next index in a cyclic index set
		inline size_t next_ind(size_t ind, size_t numInds) const	{return (ind + 1) % numInds;}
	///	returns the previous index in a cyclic index set
		inline size_t prev_ind(size_t ind, size_t numInds) const	{return (ind + numInds - 1) % numInds;}
		
		
	private:
		std::vector<number>	m_betas;//< precalculated betas.
};

/// @}	//	end of add_to_group

}//	end of namespace

#endif
