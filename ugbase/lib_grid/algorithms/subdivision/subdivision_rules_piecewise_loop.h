//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m10 d9

#ifndef __H__UG__SUBDIVISION_RULES_PIECEWISE_LOOP__
#define __H__UG__SUBDIVISION_RULES_PIECEWISE_LOOP__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

///	A singleton that stores all rules for a piecewise-loop subdivision surface.
class SubdivRules_PLoop
{
	public:
	///	returns the only instance to this singleton.
		static SubdivRules_PLoop& inst()
		{
			static SubdivRules_PLoop subdivRules;
			return subdivRules;
		}

	////////////////////////////////
	//	WEIGHTS
		number ref_even_inner_center_weight(size_t valence)
			{return 1. - (number)valence * get_beta(valence);}
		
		number ref_even_inner_nbr_weight(size_t valence)
			{return get_beta(valence);}
		
	///	returns weights for center, nbr1 and nbr2.
		vector3 ref_even_crease_weights()
			{return vector3(0.75, 0.125, 0.125);}
			
		vector4 ref_odd_inner_weight()
			{return vector4(0.375, 0.375, 0.125, 0.125);}
			
		//vector4 ref_odd_inner_weight(size_t creaseValence);
		vector2 ref_odd_crease_weight()
			{return vector2(0.5, 0.5);}
/*
		number proj_inner_center_weight(size_t valence);
		number proj_inner_nbr_weight(size_t valence);
		number proj_crease_center_weight(size_t valence);
		number proj_crease_nbr_weight(size_t valence);
		number proj_inner_crease_nbr_center_weight(size_t valence);
		number proj_inner_crease_nbr_nbr_weight(size_t valence);
		number proj_inner_crease_nbr_nbr_weight(size_t valence, size_t cValence);
*/		
	///	returns beta as it is used in the subdivision masks.
	/**	performs a lookup if the valency is small enough.
	 *	calculates a fresh beta else.*/
		number get_beta(size_t valency);
/*
	////////////////////////////////
	//	EVEN MASKS
		template <class TAAPos>
		typename TAAPos::ValueType
		apply_even_mask(Grid& grid, VertexBase* center,
						TAAPos& aaPos);

		template <class TAAPos>
		typename TAAPos::ValueType
		apply_even_crease_mask(VertexBase* center, VertexBase* nbr1,
							   VertexBase* nbr2, TAAPos& aaPos);

	////////////////////////////////
	//	ODD MASKS
		template <class TAAPos>
		typename TAAPos::ValueType
		apply_odd_mask(VertexBase* vrt, EdgeBase* parent,
					   TAAPos& aaPos);

		template <class TAAPos>
		typename TAAPos::ValueType
		apply_odd_crease_mask(VertexBase* vrt, EdgeBase* parent,
							  TAAPos& aaPos);

		template <class TAAPos>
		typename TAAPos::ValueType
		apply_odd_crease_nbr_mask(Grid& grid, VertexBase* vrt,
								  EdgeBase* parent, TAAPos& aaPos);

	////////////////////////////////
	//	PROJECTION
		template <class TAAPos>
		typename TAAPos::ValueType
		project_inner_vertex(Grid& grid, VertexBase* vrt,
							 TAAPos& aaPos);

		template <class TAAPos>
		typename TAAPos::ValueType
		project_inner_vertex(VertexBase* vrt, VertexBase* nbrs,
							 int* nbrCreaseValencies, int numNbrs,
							 TAAPos& aaPos);

		template <class TAAPos>
		typename TAAPos::ValueType
		project_crease_vertex(VertexBase* vrt, VertexBase* nbr1,
							  VertexBase* nbr2, TAAPos& aaPos);
*/
	private:
	///	calculates beta as it is used in the subdivision masks.
		number calculate_beta(size_t valency);
		
	private:
	///	private constructor prohibits multiple instantiation.
		SubdivRules_PLoop();
		
	///	private copy constructor prohibits copy-construction.
		SubdivRules_PLoop(const SubdivRules_PLoop& src);
		
	///	private assignment operator prohibits assignment.
		SubdivRules_PLoop& operator=(const SubdivRules_PLoop& src);
		
	private:
		std::vector<number>	m_betas;//< precalculated betas.
};

}//	end of namespace

#endif