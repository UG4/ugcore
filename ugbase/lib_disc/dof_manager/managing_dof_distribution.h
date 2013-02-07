/*
 * managing_dof_distribution.h
 *
 *  Created on: 21.02.2012
 *      Author: andreasvogel
 */

#ifndef MANAGING_DOF_DISTRIBUTION_H_
#define MANAGING_DOF_DISTRIBUTION_H_

namespace ug{

///	predeclaration
class IGridFunction;
class LevelMGDoFDistribution;

class ManagingDoFDistribution
{
	friend class LevelMGDoFDistribution;

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
	/**
	 * This vector holds a pointer to all grid functions, that should be
	 * managed (i.e. adapted), when the dof distribution is changed.
	 * NOTE:	No SmartPtr is used here on purpose. The GridFunction stores a
	 * 			SmartPtr to the DoFDistribution. If we would use a SmartPtr
	 * 			here those objects would reference each other and would never
	 * 			be deleted.
	 */
		std::vector<IGridFunction*> m_vpGridFunction;
};

} // end namespace ug

#endif /* MANAGING_DOF_DISTRIBUTION_H_ */
