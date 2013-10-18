// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.07.2011 (m,d,y)
 
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
template <class TDomain>
void TranslateDomain(TDomain& dom, ISelector& sel, const vector3& offset)

{
	typedef typename TDomain::position_type pos_t;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t o;
	VecCopy(o, offset, 0);

//	perform the translation
	vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

	for(size_t i = 0; i < vrts.size(); ++i){
		pos_t& v = aaPos[vrts[i]];
		VecAdd(v, v, o);
	}
}


////////////////////////////////////////////////////////////////////////////////
///	Scales selected elements around the given center.
template <class TDomain>
void ScaleDomain(TDomain& dom, const vector3& center, const vector3& scale)

{
	typedef typename TDomain::position_type pos_t;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();
	typedef typename TDomain::grid_type	grid_t;
	typedef typename grid_t::template traits<VertexBase>::iterator	vrt_iterator_t;

	grid_t& g = *dom.grid();
	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//	perform the scaling
	for(vrt_iterator_t iter = g.template begin<VertexBase>();
		iter != g.template end<VertexBase>(); ++iter)
	{
		pos_t& v = aaPos[*iter];
		for(size_t j = 0; j < pos_t::Size; ++j)
			v[j] = c[j] + (v[j] - c[j]) * s[j];
	}
}


////////////////////////////////////////////////////////////////////////////////
///	Scales selected elements around the given center.
template <class TDomain>
void ScaleDomain(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	typedef typename TDomain::position_type pos_t;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//	perform the scaling
	vector<VertexBase*> vrts;
	CollectVerticesTouchingSelection(vrts, sel);

    if(vrts.empty())
        return;

	for(size_t i = 0; i < vrts.size(); ++i){
		pos_t& v = aaPos[vrts[i]];

		for(size_t j = 0; j < pos_t::Size; ++j)
			v[j] = c[j] + (v[j] - c[j]) * s[j];
	}
}

template <class TDomain>
void ScaleDomainSquaredWeighting(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	typedef typename TDomain::position_type pos_t;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//  prepare data
	vector<VertexBase*> vrts;
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

template <class TDomain>
void ScaleDomainWeighting(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	typedef typename TDomain::position_type pos_t;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//  prepare data
	vector<VertexBase*> vrts;
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


template <class TDomain>
void ScaleDomainSqrtWeighting(TDomain& dom, ISelector& sel, const vector3& center,
				 const vector3& scale)

{
	typedef typename TDomain::position_type pos_t;
	typename TDomain::position_accessor_type& aaPos = dom.position_accessor();

	pos_t c, s;
	VecCopy(c, center, 0);
	VecCopy(s, scale, 1);

//  prepare data
	vector<VertexBase*> vrts;
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
 * @param parentGroup		group for sorting of functionality
 */
template <typename TDomain>
static void Domain(Registry& reg, string grp)
{
	typedef TDomain 							domain_type;
	typedef typename TDomain::position_type		pos_type;

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
	typedef Transform::Functionality Functionality;

	try{
		RegisterDomainDependent<Functionality>(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace bridge
}// end of namespace ug
