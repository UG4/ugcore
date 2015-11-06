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
		sel.select(g->begin<Edge>(), g->end<Edge>(), status);
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
			SelectAssociated<Vertex>(sel, goc.begin<Edge>(lvl),
										 goc.end<Edge>(lvl), status);
			SelectAssociated<Vertex>(sel, goc.begin<Face>(lvl),
										 goc.end<Face>(lvl), status);
			SelectAssociated<Vertex>(sel, goc.begin<Volume>(lvl),
										 goc.end<Volume>(lvl), status);
		}
		if(selectEdges){
			SelectAssociated<Edge>(sel, goc.begin<Face>(lvl),
										 goc.end<Face>(lvl), status);
			SelectAssociated<Edge>(sel, goc.begin<Volume>(lvl),
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
		SelectSubsetElements<Edge>(sel, sh, subsetIndex, status);
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
