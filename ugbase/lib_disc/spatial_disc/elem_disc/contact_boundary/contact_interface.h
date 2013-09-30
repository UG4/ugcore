/*
 * contact_interface.h
 *
 *  Created on: 06.05.2013
 *      Author: raphaelprohl
 */

#ifndef CONTACT_INTERFACE_H_
#define CONTACT_INTERFACE_H_

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/common/multi_index.h"

namespace ug{

template <typename TDomain, typename TGridFunction>
class IContactDisc
{
	private:
	///	own type
		typedef IContactDisc<TDomain, TGridFunction> this_type;

	public:
	///	Domain type
		typedef TDomain domain_type;

	///	World dimension
		static const int dim = TDomain::dim;

	public:
		IContactDisc(){};

	/// Virtual destructor
		virtual ~IContactDisc() {}

		virtual void contactForces(TGridFunction& force, const TGridFunction& u,
				std::vector<DoFIndex> vActiveSet, std::vector<int> vActiveSubsets) = 0;
};

} //end namespace ug

#endif /* CONTACT_INTERFACE_H_ */
