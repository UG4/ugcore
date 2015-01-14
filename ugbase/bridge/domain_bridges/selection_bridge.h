// author: stephan
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
				  bool selectEdges, bool selectFaces, bool selectVolumes);


////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects associated (low dimensional) elements.
void SelectAssociatedElements(ISelector& sel, bool bSelect, bool selectVrts,
				  	  	  	  bool selectEdges, bool selectFaces);


////////////////////////////////////////////////////////////////////////////////
///	Selects / Deselects elements in the given subset of the given domain
template <class TDomain>
void SelectDomainSubset(ISelector& sel, TDomain& dom, int subsetIndex,
						bool bSelect, bool selectVrts, bool selectEdges,
						bool selectFaces, bool selectVolumes);
}
