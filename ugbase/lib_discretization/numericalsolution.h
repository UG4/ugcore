/*
 * numericalsolution.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__NUMERICALSOLUTION__
#define __H__LIBDISCRETIZATION__NUMERICALSOLUTION__

// extern libraries
#include <string>
#include <algorithm>

// other ug4 modules
#include "../lib_grid/lib_grid.h"
#include "../lib_algebra/lib_algebra.h"

// library intern headers
#include "trialspace.h"
#include "referenceelement.h"
#include "domain.h"
#include "dofpattern.h"

namespace ug {


template <int d>
class NumericalSolution {
	public:
		/// constructor
		NumericalSolution(std::string shortname, std::string longname, std::string description, DoFPattern& pattern, Domain<d>& domain);

		/// sets the value of the function
		bool set_values(bool (*fct)(MathVector<d>, number&), uint nr_fct, ISubsetHandler& sh, uint subsetIndex, uint level = 0);

		/// returns the associated GridVector
		Vector* GridVector(uint level = 0);

		/// returns the associated pattern
		DoFPattern* get_pattern();

		/// returns the trial space type of numerical solution 'nr_fct'
		TrialSpaceType get_TrialSpaceType(uint nr_func);

		/// returns the trial space of numerical solution 'nr_fct'
		//template <class TElem>
		//const TrialSpace<TElem>& get_TrialSpace(uint nr_fct);

		/// returns the associated domain
		Domain<d>* get_domain();

		/// returns the number of dofs
		uint num_dofs(uint level = 0);

		/// returns the dof indices of a given element for the solution 'nr_fct'
		template<typename TElem>
		bool get_indices(TElem* elem, uint nr_fct, int* ind);

		/// returns the dof values of a given element for the solution 'nr_fct'
		template <class TElem>
		bool get_local_DoFValues(TElem* elem, uint nr_func, number* DoFValues);

		/// print information
		bool print();

		/// assignment
		bool assign(NumericalSolution<d>& v);

		/// destructor
		~NumericalSolution();

		/// name of numerical solution 'nr_fct'
		std::string get_name(uint nr_fct);

	protected:
		std::string _shortname;
		std::string _longname;
		std::string _description;

		std::vector<Vector*> _GridVector;
		DoFPattern* _pattern;
		Domain<d>* _domain;
};


} /* end namespace ug */

#include "numericalsolution_impl.h"

#endif /* __H__LIBDISCRETIZATION__NUMERICALSOLUTION__ */
