
#ifndef LAGRANGE_MULTIPLIER_DISC_INTERFACE_H_
#define LAGRANGE_MULTIPLIER_DISC_INTERFACE_H_

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/common/multi_index.h"

namespace ug{

template <typename TDomain, typename TGridFunction>
class ILagrangeMultiplierDisc
{
	private:
	///	own type
		typedef ILagrangeMultiplierDisc<TDomain, TGridFunction> this_type;

	public:
	///	Domain type
		typedef TDomain domain_type;

	///	World dimension
		static const int dim = TDomain::dim;

	public:
		ILagrangeMultiplierDisc(){};

	/// Virtual destructor
		virtual ~ILagrangeMultiplierDisc() {}

		virtual void lagrange_multiplier(TGridFunction& lagMult, const TGridFunction& u,
				std::vector<DoFIndex> vActiveSet, std::vector<int> vActiveSubsets) = 0;
};

} //end namespace ug

#endif /* LAGRANGE_MULTIPLIER_DISC_INTERFACE_H_ */
