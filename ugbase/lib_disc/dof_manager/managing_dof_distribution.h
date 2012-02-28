/*
 * managing_dof_distribution.h
 *
 *  Created on: 21.02.2012
 *      Author: andreasvogel
 */

#ifndef MANAGING_DOF_DISTRIBUTION_H_
#define MANAGING_DOF_DISTRIBUTION_H_

#include "lib_disc/function_spaces/grid_function.h"

namespace ug{

class ManagingDoFDistribution
{
	public:
	///	registers a grid function for adaptation management
		void manage_grid_function(IGridFunction& gridFct);

	///	unregisters a grid function for adaptation management
		void unmanage_grid_function(IGridFunction& gridFct);

	protected:
	///	permutes values in managed functions, if indices permuted
		void permute_values(const std::vector<size_t>& vIndNew);

	///	swaps values in managed functions, if indices swapped
		void copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,
		                                         bool bDisjunct);

	///	changes values in managed functions, number of indices changed
		void resize_values(size_t newSize);

	protected:
	///	managed grid functions
		std::vector<IGridFunction*> m_vpGridFunction;
};

} // end namespace ug

#endif /* MANAGING_DOF_DISTRIBUTION_H_ */
