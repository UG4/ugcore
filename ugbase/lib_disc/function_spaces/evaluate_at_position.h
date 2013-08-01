/*
 * evaluate_at_closest_vertex.h
 *
 *  Created on: 01.08.2013
 *      Author: imuha
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACES__EVALUATE_AT_CLOSEST_VERTEX__
#define __H__UG__LIB_DISC__FUNCTION_SPACES__EVALUATE_AT_CLOSEST_VERTEX__

#include "common/common.h"
#include "common/util/smart_pointer.h"

#include "lib_disc/domain_util.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/reference_element/reference_mapping.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug{

template <typename TGridFunction>
number Evaluate_at_vertex2D(number pos_x, number pos_y,
                           SmartPtr<TGridFunction> spGridFct,
                           size_t fct,
                           const SubsetGroup& ssGrp)
{

//	domain type and position_type
	typedef typename TGridFunction::domain_type domain_type;
	typedef typename domain_type::position_type position_type;

// get position accessor
	const typename domain_type::position_accessor_type& aaPos
										= spGridFct->domain()->position_accessor();

	std::vector<MultiIndex<2> > ind;
	typename TGridFunction::template dim_traits<0>::const_iterator iterEnd, iter,chosen;
	double min_distance;

	for(size_t i = 0; i < ssGrp.size(); ++i)
	{
	//	get subset index
		const int si = ssGrp[i];

	//	skip if function is not defined in subset
		if(!spGridFct->is_def_in_subset(fct, si)) continue;

	// 	iterate over all elements
		iterEnd = spGridFct->template end<VertexBase>(si);
		iter = spGridFct->template begin<VertexBase>(si);
		for(; iter != iterEnd; ++iter)
		{
		//	get element
			VertexBase* vrt = *iter;

		//	global position
			position_type glob_pos = aaPos[vrt];
			if(i==0 && iter == spGridFct->template begin<VertexBase>(si))
			{
				min_distance = 0;
				min_distance += (glob_pos[0]-pos_x)*(glob_pos[0]-pos_x);
				min_distance += (glob_pos[1]-pos_y)*(glob_pos[1]-pos_y);
				min_distance = sqrt(min_distance);
				chosen = iter;
			}
			else
			{
				double buffer = 0;
				min_distance += (glob_pos[0]-pos_x)*(glob_pos[0]-pos_x);
				min_distance += (glob_pos[1]-pos_y)*(glob_pos[1]-pos_y);
				buffer = sqrt(min_distance);
				if(buffer<min_distance)
				{
					min_distance = buffer;
					chosen = iter;
				}
			}
		}


	}
	VertexBase* vrt = *chosen;
	spGridFct->multi_indices(vrt, fct, ind);
		if(ind.size()>0)
		{
			UG_THROW("Only one component supported.");
		}
	return 	BlockRef((*spGridFct)[ind[0][0]], ind[0][1]);

}


template <typename TGridFunction>
number Evaluate_at_closest_vertex(number pos_x,number pos_y,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{

//	get function id of name
	const size_t fct = spGridFct->fct_id_by_name(cmp);

//	check that function found
	if(fct > spGridFct->num_fct())
		UG_THROW("Evaluate: Name of component '"<<cmp<<"' not found.");


//	create subset group
	SubsetGroup ssGrp(spGridFct->domain()->subset_handler());
	if(subsets != NULL)
	{
		ssGrp.add(TokenizeString(subsets));
	}
	else
	{
	//	add all subsets and remove lower dim subsets afterwards
		ssGrp.add_all();
	}

	return Evaluate_at_vertex2D<TGridFunction>(pos_x,pos_y, spGridFct, fct, ssGrp);

}

/*
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 number time)
{Interpolate(spInterpolFunction, spGridFct, cmp, NULL, time);}
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp,
                 const char* subsets)
{Interpolate(spInterpolFunction, spGridFct, cmp, subsets, 0.0);}
template <typename TGridFunction>
void Interpolate(SmartPtr<UserData<number, TGridFunction::dim> > spInterpolFunction,
                 SmartPtr<TGridFunction> spGridFct, const char* cmp)
{Interpolate(spInterpolFunction, spGridFct, cmp, NULL, 0.0);}
*/


} // namespace ug

#endif /*EVALUATE_AT_CLOSEST_VERTEX*/
