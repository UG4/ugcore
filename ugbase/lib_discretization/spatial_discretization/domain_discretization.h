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
#include "common/util/string_util.h"

// library intern headers
#include "./domain_discretization_interface.h"
#include "./elem_disc/elem_disc_assemble_util.h"
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
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	Empty Constructor
		DomainDiscretization() : m_bForceRegGrid(false)
		{
			m_vvPostProcess.resize(PPT_NUM_POST_PROCESS_TYPES);
		};

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

	///////////////////////////
	// Mass and Stiffness Matrix
	///////////////////////////

	/// assembles the mass matrix
	IAssembleReturn assemble_mass_matrix(matrix_type& M, const vector_type& u,
	                                     const dof_distribution_type& dofDistr);

	/// assembles the stiffness matrix
	IAssembleReturn assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
	                                          const dof_distribution_type& dofDistr);

	/// forces the assembling to consider the grid as regular
	virtual void force_regular_grid(bool bForce) {m_bForceRegGrid = bForce;}

	protected:
	/// forces the assembling to regard the grid as regular
		bool m_bForceRegGrid;

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
		std::vector<std::vector<IPostProcess<TDoFDistribution, TAlgebra>*> >
			m_vvPostProcess;

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
};

/// @}

} // end namespace ug

// inlcude documentation
#include "domain_discretization_impl.h"

#endif /* __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__ */
