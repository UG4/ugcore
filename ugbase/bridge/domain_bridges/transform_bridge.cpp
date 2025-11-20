/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

using namespace std;

namespace ug{

/**
 * \defgroup transform_bridge Transformation Bridge
 * \ingroup domain_bridge
 * \{
 */

////////////////////////////////////////////////////////////////////////////////
///	Translates (moves) selected elements by the given offset
template <typename TDomain>
void TranslateDomain(TDomain& dom, ISelector& sel, const vector3& offset)

{
	using pos_t = typename TDomain::position_type;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t o;
	VecCopy(o, offset, 0);

//	perform the translation
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

	for(size_t i = 0; i < vrts.size(); ++i){
		pos_t& v = aaPos[vrts[i]];
		VecAdd(v, v, o);
	}
}


////////////////////////////////////////////////////////////////////////////////
///	Scales selected elements around the given center.
template <typename TDomain>
void ScaleDomain(TDomain& dom, const vector3& center, const vector3& scale)

{
	using pos_t = typename TDomain::position_type;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	using grid_t = typename TDomain::grid_type;
	using vrt_iterator_t = typename grid_t::template traits<Vertex>::iterator;

	grid_t& g = *dom.grid();
	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//	perform the scaling
	for(vrt_iterator_t iter = g.template begin<Vertex>();
		iter != g.template end<Vertex>(); ++iter)
	{
		pos_t& v = aaPos[*iter];
		for(size_t j = 0; j < pos_t::Size; ++j)
			v[j] = c[j] + (v[j] - c[j]) * s[j];
	}
}


////////////////////////////////////////////////////////////////////////////////
///	Scales selected elements around the given center.
template <typename TDomain>
void ScaleDomain(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	using pos_t = typename TDomain::position_type;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//	perform the scaling
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

    if(vrts.empty()){
        return;
	}

	for(size_t i = 0; i < vrts.size(); ++i){
		pos_t& v = aaPos[vrts[i]];

		for(size_t j = 0; j < pos_t::Size; ++j)
			v[j] = c[j] + (v[j] - c[j]) * s[j];
	}
}

template <typename TDomain>
void ScaleDomainSquaredWeighting(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	using pos_t = typename TDomain::position_type;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//  prepare data
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);
    if(vrts.empty())
        return;

    number dMaxSq = -1;
    for(size_t i = 0; i < vrts.size(); ++i){
        number distSq = VecDistanceSq(c,  aaPos[vrts[i]]);
        if( distSq > dMaxSq )
                dMaxSq = distSq;
    }
    number dMax = sqrt(dMaxSq);

//  perform the scaling
    for(size_t i = 0; i < vrts.size(); ++i){
        pos_t& v = aaPos[vrts[i]];
        number dist = VecDistance(c, v);
        number sqrdWeight = (dist/dMax)*(dist/dMax);

        for(size_t j = 0; j < pos_t::Size; ++j)
            v[j] = c[j] + (v[j] - c[j]) * ((s[j]-1.0)*sqrdWeight + 1.0);
    }

}

template <typename TDomain>
void ScaleDomainWeighting(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	using pos_t = typename TDomain::position_type;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//  prepare data
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);
    if(vrts.empty())
        return;

    number dMaxSq = -1;
    for(size_t i = 0; i < vrts.size(); ++i){
        number distSq = VecDistanceSq(c,  aaPos[vrts[i]]);
        if( distSq > dMaxSq )
                dMaxSq = distSq;
    }
    number dMax = sqrt(dMaxSq);

//  perform the scaling
    for(size_t i = 0; i < vrts.size(); ++i){
        pos_t& v = aaPos[vrts[i]];
        number dist = VecDistance(c, v);
		number weight = dist/dMax;

        for(size_t j = 0; j < pos_t::Size; ++j)
            v[j] = c[j] + (v[j] - c[j]) * ((s[j]-1.0)*weight + 1.0);
    }
}


template <typename TDomain>
void ScaleDomainSqrtWeighting(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	using pos_t = typename TDomain::position_type;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//  prepare data
	vector<Vertex*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);
    if(vrts.empty())
        return;

    number dMaxSq = -1;
    for(size_t i = 0; i < vrts.size(); ++i){
        number distSq = VecDistanceSq(c,  aaPos[vrts[i]]);
        if( distSq > dMaxSq )
                dMaxSq = distSq;
    }
    number dMax = sqrt(dMaxSq);

//  perform the scaling
    for(size_t i = 0; i < vrts.size(); ++i){
        pos_t& v = aaPos[vrts[i]];
        number dist = VecDistance(c, v);
		number weight = sqrt(dist/dMax);

        for(size_t j = 0; j < pos_t::Size; ++j)
            v[j] = c[j] + (v[j] - c[j]) * ((s[j]-1.0)*weight + 1.0);
    }

}

// end group transform_bridge
/// \}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace bridge{
namespace Transform{

/// \addtogroup transform_bridge
/// \{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain dependent parts.
 * All Functions and Classes depending on the Domain
 * are to be placed here when registering. The method is called for all
 * available Domain types, based on the current build options.
 *
 * @param reg				registry
 * @param grp				group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	using domain_type = TDomain;

	reg.add_function("TranslateDomain", &TranslateDomain<domain_type>, grp, "", "dom#sel#offset");
	reg.add_function("ScaleDomain",
			static_cast<void (*)(domain_type&, const vector3&, const vector3&)>
				(&ScaleDomain<domain_type>), grp, "", "dom#center#scale");
	reg.add_function("ScaleDomain",
			static_cast<void (*)(domain_type&, ISelector&, const vector3&, const vector3&)>
				(&ScaleDomain<domain_type>), grp, "", "dom#sel#center#scale");
	reg.add_function("ScaleDomainSqrtWeighting", &ScaleDomainSqrtWeighting<domain_type>, grp, "", "dom#sel#center#scale");
	reg.add_function("ScaleDomainWeighting", &ScaleDomainWeighting<domain_type>, grp, "", "dom#sel#center#scale");
	reg.add_function("ScaleDomainSquaredWeighting", &ScaleDomainSquaredWeighting<domain_type>, grp, "", "dom#sel#center#scale");

}

}; // end Functionality

// end group transform_bridge
/// \}

}// end Refinement

/// \addtogroup transform_bridge
void RegisterBridge_Transform(Registry& reg, string grp)
{
	grp.append("/Transform");
	using Functionality = Transform::Functionality;

	try{
		RegisterDomainDependent<Functionality>(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace bridge
}// end of namespace ug
