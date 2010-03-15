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
#include <boost/any.hpp>

// other ug4 modules
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"
#include "lib_discretization/reference_element/reference_elements.h"
#include "lib_discretization/domain.h"
#include "lib_discretization/dof_manager/dof_manager.h"

namespace ug {
/*
template<int dimDomain>
class IDiscreteFunction{
	public:
		virtual void evaluate(MathVector<dimDomain> x) = 0;
};

template<int dimDomain>
class ILinearOperator{
	// init this function
	virtual void init() = 0;

	// L(u) = v
	virtual void apply(const IDiscreteFunction<dimDomain>& u, IDiscreteFunction<dimDomain>& v) = 0;
};

template <int dimDomain>
class DiscreteFunction{
	public:
		typedef DoFPattern::subvec_type subvec_type;
	public:
		DiscreteFunction(const Domain<dimDomain>& domain, const DoFPattern& pattern) : m_domain(domain), m_pattern(pattern)
		{
		};

		bool add_sub_vector(subvec_type subvec)
		{
			for(std::size_t i = 0; i < subvec.size(); ++i)
			{
				const subvec_type::index_type& index = subvec.index(i);
				subvec_type::index_type::single_index_type subsetIndex = index[0];
				subvec_type::index_type::single_index_type DoFGroup = index[1];

			}
		}

	private:
		template <typename BlockType>
		Vector<BlockType>* get_sparse_matrix(int subsetIndex, uint DoFGroup)
		{
			return *(boost::any_cast<Vector<BlockType> >(&m_vector[subsetIndex][DoFGroup]));
		}

		template <typename BlockType>
		bool set_sparse_matrix(int subsetIndex, uint DoFGroup, Vector<BlockType>* vec)
		{
			m_vector[subsetIndex][DoFGroup] = *vec;
			return true;
		}


	private:
		std::vector<std::vector<boost::any> > m_vector;
		std::vector<std::vector<int> > m_blockSize;
		Domain<dimDomain>& m_domain;
		DoFPattern& m_pattern;
};
*/

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
