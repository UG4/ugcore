/*
 * domain_discretization.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__
#define __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__

// other ug4 modules
#include "common/common.h"
#include "common/string_util.h"
#include "lib_algebra/lib_algebra.h"

// library intern headers
#include "./domain_discretization_interface.h"
#include "./elem_disc/elem_disc_assemble_util.h"
#include "./coupled_elem_disc/coupled_elem_disc_assemble_util.h"
#include "./post_process/post_process_interface.h"
#include "./subset_assemble_util.h"
#include "lib_discretization/common/function_group.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// domain discretization implementing the interface
/**
 * This class is an implementation of the IDomainDiscretization interface. It
 * is designed to simply group several discreizations on different subsets and
 * perform element based assemblings and post processes in the same order.
 */
template <	typename TDoFDistribution,
			typename TAlgebra>
class DomainDiscretization :
	public IDomainDiscretization<TDoFDistribution, TAlgebra>
{
	protected:
	///	Type of DoF Distribution
		typedef TDoFDistribution dof_distribution_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	Empty Constructor
		DomainDiscretization()
		{};

	///////////////////////////
	// Time independent part
	///////////////////////////

	/// \copydoc IAssemble::assemble_jacobian()
	IAssembleReturn assemble_jacobian(matrix_type& J, const vector_type& u,
	                                  const dof_distribution_type& dofDistr);

	/// \copydoc IAssemble::assemble_defect()
	IAssembleReturn assemble_defect(vector_type& d, const vector_type& u,
	                                const dof_distribution_type& dofDistr);

	/// \copydoc IAssemble::assemble_linear()
	IAssembleReturn assemble_linear(matrix_type& A, vector_type& b,
	                                const vector_type& u,
	                                const dof_distribution_type& dofDistr);

	/// \copydoc IAssemble::assemble_solution()
	IAssembleReturn assemble_solution(vector_type& u,
	                                  const dof_distribution_type& dofDistr);

	///////////////////////
	// Time dependent part
	///////////////////////

	/// \copydoc IDomainDiscretization::assemble_jacobian()
	virtual
	IAssembleReturn assemble_jacobian(matrix_type& J, const vector_type& u,
	                                  const dof_distribution_type& dofDistr,
	                                  number time, number s_m, number s_a);

	/// \copydoc IDomainDiscretization::assemble_defect()
	virtual
	IAssembleReturn assemble_defect(vector_type& d, const vector_type& u,
	                                const dof_distribution_type& dofDistr,
									number time, number s_m, number s_a);

	/// \copydoc IDomainDiscretization::assemble_linear()
	virtual
	IAssembleReturn assemble_linear(matrix_type& A, vector_type& b,
	                                const vector_type& u,
	                                const dof_distribution_type& dofDistr,
									number time, number s_m, number s_a);

	/// \copydoc IDomainDiscretization::assemble_solution()
	virtual
	IAssembleReturn assemble_solution(vector_type& u,
	                                  const dof_distribution_type& dofDistr,
	                                  number time);

	public:
	/// adds an element discretization to the assembling process
	/**
	 * This function adds an Element Discretization to the assembling. During
	 * the assembling process the elem disc will assemble a given group of
	 * subsets and add the output to the passed vectors and matrices.
	 * For each Elem Disc one loop over all elements of the subset will be
	 * performed.
	 *
	 * \param[in] 	elem		Element Discretization to be added
	 */
		bool add_elem_disc(IElemDisc<TAlgebra>& elem)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(m_vElemDisc[i] == &elem)
					return true;

		//	add it
			m_vElemDisc.push_back(&elem);
			return true;
		}

	protected:
	///	vector holding all registered elem discs
		std::vector<IElemDisc<TAlgebra>*> m_vElemDisc;

	public:
	/// adds a Coupled Element Discretization that will be performed on given subsets
	/**
	 * This function adds a Coupled Element Discretization for a given group
	 * of functions. The Functions must match the requiries of the
	 * Discretization and must be defined for the Subsets.
	 *
	 * \param[in] 	coupledSystem	coupled element discretization
	 * \param[in]	functionGroup	Function used in discretization
	 * \param[in]	subsetGroup		Subsets, where discretization should be performed
	 */
		bool add(CoupledSystem<TAlgebra>& coupledSystem,
		         const FunctionGroup& functionGroup,
		         const SubsetGroup& subsetGroup);

	/// adds a Coupled Element Discretization that will be performed on given subsets
	/**
	 * This function adds a Coupled Element Discretization for a given group of
	 * functions. The Functions must match the requiries of the Discretization
	 * and must be defined for the Subsets.
	 *
	 * \param[in] 	coupledSystem	coupled element discretization
	 * \param[in]	pattern			underlying function pattern to be used
	 * \param[in]	functions		Function Names
	 * \param[in]	subsets			Subset Names
	 */
		bool add(CoupledSystem<TAlgebra>& coupledSystem,
		         const FunctionPattern& pattern,
		         const char* functions,
		         const char* subsets);

	protected:
		struct CoupledDisc
		{
			CoupledDisc(CoupledSystem<TAlgebra>& disc_,
			            const FunctionGroup& fcts_,
			            const SubsetGroup& subsetGroup_) :
				disc(&disc_), functionGroup(fcts_), subsetGroup(subsetGroup_)
			{};

			CoupledSystem<TAlgebra>* disc;
			FunctionGroup functionGroup;
			SubsetGroup subsetGroup;
		};
		std::vector<CoupledDisc> m_vCoupledDisc;

	public:
	/// adds a post process to the assembling process
	/**
	 * This function adds a Post Process to the assembling. The post process is
	 * called when all element-wise assembling have been performed.
	 *
	 * \param[in] 	pp		Post Process to be added
	 */
		bool add_post_process(IPostProcess<TDoFDistribution, TAlgebra>& pp)
		{
		// 	get type of post process
			const int type = pp.type();

		//	check that not already registered
			for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
				if(m_vvPostProcess[type][i] == &pp)
					return true;

		//	add post process
			m_vvPostProcess[type].push_back(&pp);
			return true;
		}

	protected:
	//	vector holding all registered post processes
		std::vector<IPostProcess<TDoFDistribution, TAlgebra>*>
			m_vvPostProcess[PPT_NUM_POST_PROCESS_TYPES];

	protected:
	///	returns number of registered post processes
		virtual size_t num_post_process() const
		{
			return m_vvPostProcess[PPT_DIRICHLET].size();
		}

	///	returns the i'th post process
		virtual IPostProcess<TDoFDistribution, TAlgebra>* get_post_process(size_t i)
		{
			return m_vvPostProcess[PPT_DIRICHLET].at(i);
		}

	// todo: What is this function used for???? Do we have to include it
	// todo: Remove if possible.
		virtual size_t num_fct() const
		{
			size_t sum = 0;
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				sum += m_vElemDisc[i]->get_function_group().num_fct();

			return sum;
		}
};

/// @}

} // end namespace ug

// inlcude documentation
#include "domain_discretization_impl.h"

#endif /* __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__ */
