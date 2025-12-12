/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#include "local_transfer_interface.h"

#include "local_transfer.h"

namespace ug {

template <typename TDomain>
SmartPtr<IElemProlongation<TDomain> >
GetStandardElementProlongation(const LFEID& lfeid)
{
	switch(lfeid.type()){
		case LFEID::LAGRANGE:
			if(lfeid.order() == 1)
				return make_sp(new P1LagrangeElemTransfer<TDomain>(lfeid));
			else
				return make_sp(new StdLagrangeElemTransfer<TDomain>(lfeid));
		case LFEID::PIECEWISE_CONSTANT:
			return make_sp(new PiecewiseConstantElemTransfer<TDomain>(lfeid));

		case LFEID::CROUZEIX_RAVIART:
			return make_sp(new CrouzeixRaviartElemTransfer<TDomain>(lfeid));

		default: UG_THROW("No Standard Element Prolongation found for "<<lfeid);
	}
}


template <typename TDomain>
SmartPtr<IElemRestriction<TDomain> >
GetStandardElementRestriction(const LFEID& lfeid)
{
	switch(lfeid.type()){
		case LFEID::LAGRANGE:
			if(lfeid.order() == 1)
				return make_sp(new P1LagrangeElemTransfer<TDomain>(lfeid));
			else
				return make_sp(new StdLagrangeElemTransfer<TDomain>(lfeid));

		case LFEID::PIECEWISE_CONSTANT:
			return make_sp(new PiecewiseConstantElemTransfer<TDomain>(lfeid));

		case LFEID::CROUZEIX_RAVIART:
			return make_sp(new CrouzeixRaviartElemTransfer<TDomain>(lfeid));

		default: UG_THROW("No Standard Element Restriction found for "<<lfeid);
	}
}


#ifdef UG_DIM_1
template SmartPtr<IElemProlongation<Domain1d> >
GetStandardElementProlongation<Domain1d>(const LFEID& lfeid);

template SmartPtr<IElemRestriction<Domain1d> >
GetStandardElementRestriction<Domain1d>(const LFEID& lfeid);
#endif
#ifdef UG_DIM_2
template SmartPtr<IElemProlongation<Domain2d> >
GetStandardElementProlongation<Domain2d>(const LFEID& lfeid);

template SmartPtr<IElemRestriction<Domain2d> >
GetStandardElementRestriction<Domain2d>(const LFEID& lfeid);
#endif
#ifdef UG_DIM_3
template SmartPtr<IElemProlongation<Domain3d> >
GetStandardElementProlongation<Domain3d>(const LFEID& lfeid);

template SmartPtr<IElemRestriction<Domain3d> >
GetStandardElementRestriction<Domain3d>(const LFEID& lfeid);
#endif

} // end namespace ug
