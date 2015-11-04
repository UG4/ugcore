
#include "local_transfer_interface.h"
#include "local_transfer.h"

namespace ug{

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
