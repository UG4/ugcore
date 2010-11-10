/*
 * domain_discretization.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__

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

template <	typename TDoFDistribution,
			typename TAlgebra>
class DomainDiscretization :
	public IDomainDiscretization<TDoFDistribution, TAlgebra>
{
	protected:
	// 	Type of DoF Distribution
		typedef TDoFDistribution dof_distribution_type;

	// 	Type of algebra
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		DomainDiscretization()
		{};

		///////////////////////////
		// Time independent part
		///////////////////////////
		IAssembleReturn assemble_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr);

		IAssembleReturn assemble_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr);

		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr);

		IAssembleReturn assemble_solution(vector_type& u, const dof_distribution_type& dofDistr);

		///////////////////////
		// Time dependent part
		///////////////////////
		IAssembleReturn assemble_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr,
											number time, number s_m, number s_a);

		IAssembleReturn assemble_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr,
										number time, number s_m, number s_a);

		IAssembleReturn assemble_linear(matrix_type& mat, vector_type& rhs, const vector_type& u, const dof_distribution_type& dofDistr,
										number time, number s_m, number s_a);

		IAssembleReturn assemble_solution(vector_type& u, const dof_distribution_type& dofDistr, number time);

	public:
		/** adds an element discretization to the assembling process
		 * This function adds an Element Discretization to the assembling. During the
		 * assembling process the elem disc will assemble a given group of subsets and
		 * add the output to the passed vectors and matrices. Per Elem Disc one loop
		 * over all elements of the subset will be performed.
		 *
		 * \param[in] 	elem		Element Discretization to be added
		 */
		bool add_elem_disc(IElemDisc<TAlgebra>& elem)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				if(m_vElemDisc[i] == &elem)
					return true;

			m_vElemDisc.push_back(&elem);
			return true;
		}

	protected:
		std::vector<IElemDisc<TAlgebra>*> m_vElemDisc;


	public:
		/** adds a Coupled Element Discretization that will be performed on given subsets
		 * This function adds a Coupled Element Discretization for a given group of functions. The Functions
		 * must match the requiries of the Discretization and must be defined for the Subsets.
		 *
		 * \param[in] 	coupledSystem	coupled element discretization
		 * \param[in]	functionGroup	Function used in discretization
		 * \param[in]	subsetGroup		Subsets, where discretization should be performed
		 */
		bool add(CoupledSystem<TAlgebra>& coupledSystem, const FunctionGroup& functionGroup, const SubsetGroup& subsetGroup);

		/** adds a Coupled Element Discretization that will be performed on given subsets
		 * This function adds a Coupled Element Discretization for a given group of functions. The Functions
		 * must match the requiries of the Discretization and must be defined for the Subsets.
		 *
		 * \param[in] 	coupledSystem	coupled element discretization
		 * \param[in]	pattern			underlying function pattern to be used
		 * \param[in]	functions		string of Function Names used in discretization, separated by ','
		 * \param[in]	subsets			string of Subset Names, where discretization should be performed, separated by ','
		 */
		bool add(CoupledSystem<TAlgebra>& coupledSystem, const FunctionPattern& pattern, const char* functions, const char* subsets);

	protected:
		struct CoupledDisc
		{
			CoupledDisc(CoupledSystem<TAlgebra>& disc_, const FunctionGroup& fcts_, const SubsetGroup& subsetGroup_) :
				disc(&disc_), functionGroup(fcts_), subsetGroup(subsetGroup_) {};

			CoupledSystem<TAlgebra>* disc;
			FunctionGroup functionGroup;
			SubsetGroup subsetGroup;
		};
		std::vector<CoupledDisc> m_vCoupledDisc;

	public:
		/** adds a post process to the assembling process
		 * This function adds a Post Process to the assembling. The post process is called when all
		 * element-wise assembling have been performed.
		 *
		 * \param[in] 	pp		Post Process to be added
		 */
		bool add_post_process(IPostProcess<TDoFDistribution, TAlgebra>& pp)
		{
			const int type = pp.type();

		//	check that not already registered
			for(size_t i = 0; i < m_vvPostProcess[type].size(); ++i)
				if(m_vvPostProcess[type][i] == &pp)
					return true;

			m_vvPostProcess[type].push_back(&pp);
			return true;
		}

	protected:
		std::vector<IPostProcess<TDoFDistribution, TAlgebra>*> m_vvPostProcess[PPT_NUM_POST_PROCESS_TYPES];

	protected:
		virtual size_t num_post_process() const {return m_vvPostProcess[PPT_DIRICHLET].size();}
		virtual IPostProcess<TDoFDistribution, TAlgebra>* get_post_process(size_t i)  {return m_vvPostProcess[PPT_DIRICHLET].at(i);};

		// todo: What is this function used for???? Do we have to include it
		virtual size_t num_fct() const
		{
			size_t sum = 0;
			for(size_t i = 0; i < m_vElemDisc.size(); ++i)
				sum += m_vElemDisc[i]->get_function_group().num_fct();

			return sum;
		}
};

} // end namespace ug

#include "domain_discretization_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__ */
