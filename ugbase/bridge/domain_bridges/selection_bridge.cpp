// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.07.2011 (m,d,y)
 
#include <vector>
#include <string>
#include <sstream>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_domain_dependent.h"

#include "lib_disc/domain.h"
#include "lib_grid/lib_grid.h"

using namespace std;

namespace ug{

/**
 * \defgroup selection_bridge Selection Bridge
 * \ingroup domain_bridge
 * \{
 */

////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects all elements
void SelectDomainElements(ISelector& sel, bool bSelect, bool selectVrts,
				  bool selectEdges, bool selectFaces, bool selectVolumes)
{
	Grid* g = sel.grid();

	if(!g)
		return;

	ISelector::status_t status = bSelect ? ISelector::SELECTED
										 : ISelector::DESELECTED;

	if(selectVrts)
		sel.select(g->begin<Vertex>(), g->end<Vertex>(), status);
	if(selectEdges)
		sel.select(g->begin<EdgeBase>(), g->end<EdgeBase>(), status);
	if(selectFaces)
		sel.select(g->begin<Face>(), g->end<Face>(), status);
	if(selectVolumes)
		sel.select(g->begin<Volume>(), g->end<Volume>(), status);
}

////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects associated (low dimensional) elements.
void SelectAssociatedElements(ISelector& sel, bool bSelect, bool selectVrts,
				  	  	  	  bool selectEdges, bool selectFaces)
{
	ISelector::status_t status = bSelect ? ISelector::SELECTED
										 : ISelector::DESELECTED;


	GridObjectCollection goc = sel.get_grid_objects();

	for(size_t lvl = 0; lvl < goc.num_levels(); ++lvl){
		if(selectVrts){
			SelectAssociated<Vertex>(sel, goc.begin<EdgeBase>(lvl),
										 goc.end<EdgeBase>(lvl), status);
			SelectAssociated<Vertex>(sel, goc.begin<Face>(lvl),
										 goc.end<Face>(lvl), status);
			SelectAssociated<Vertex>(sel, goc.begin<Volume>(lvl),
										 goc.end<Volume>(lvl), status);
		}
		if(selectEdges){
			SelectAssociated<EdgeBase>(sel, goc.begin<Face>(lvl),
										 goc.end<Face>(lvl), status);
			SelectAssociated<EdgeBase>(sel, goc.begin<Volume>(lvl),
									   goc.end<Volume>(lvl), status);
		}

		if(selectFaces){
			SelectAssociated<Face>(sel, goc.begin<Volume>(lvl),
								   goc.end<Volume>(lvl), status);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects elements in the given subset of the given domain
template <class TDomain>
void SelectDomainSubset(ISelector& sel, TDomain& dom, int subsetIndex,
						bool bSelect, bool selectVrts, bool selectEdges,
						bool selectFaces, bool selectVolumes)
{
	ISubsetHandler& sh = *dom.subset_handler();

	ISelector::status_t status = bSelect ? ISelector::SELECTED
										 : ISelector::DESELECTED;

	if(selectVrts)
		SelectSubsetElements<Vertex>(sel, sh, subsetIndex, status);
	if(selectEdges)
		SelectSubsetElements<EdgeBase>(sel, sh, subsetIndex, status);
	if(selectFaces)
		SelectSubsetElements<Face>(sel, sh, subsetIndex, status);
	if(selectVolumes)
		SelectSubsetElements<Volume>(sel, sh, subsetIndex, status);
}

// end group selection_bridge
/// \}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace bridge{
namespace Selection{

/// \addtogroup selection_bridge
/// \{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
//	register domain independent mark methods
	reg.add_function("SelectDomainElements", &SelectDomainElements, grp, "", "sel#bSelect#selectVrts#selectEdges#selectFaces#selectVolumes");
	reg.add_function("SelectAssociatedElements", &SelectAssociatedElements, grp, "", "sel#bSelect#selectVrts#selectEdges#selectFaces");
}

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

	reg.add_function("SelectDomainSubset",
					 &SelectDomainSubset<domain_type>, grp,
					 "", "sel#dom#subsetIndex#bSelect#selectVrts#selectEdges#selectFaces#selectVolumes");
}

}; // end Functionality

// end group selection_bridge
/// \}

}// end Refinement

/// \addtogroup selection_bridge
void RegisterBridge_Selection(Registry& reg, string grp)
{
	grp.append("/Selection");
	typedef Selection::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg, grp);
		RegisterDomainDependent<Functionality>(reg, grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

}// end of namespace bridge
}// end of namespace ug
