/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
 *  (main parts are based on the structure of
 *  	newton.h and some ideas of Sebastian Reiter & Andreas Vogel)
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_GAUSS_SEIDEL__NL_GAUSS_SEIDEL_H_
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_GAUSS_SEIDEL__NL_GAUSS_SEIDEL_H_

#include "lib_algebra/operator/interface/operator_inverse.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/spatial_disc/local_to_global/local_to_global_mapper.h"

namespace ug {

template <typename TAlgebra>
class LocalToGlobalMapperNLGS : public ILocalToGlobalMapper<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	///	default constructor
		LocalToGlobalMapperNLGS() {}

	///	adds a local vector to the global rhs
		void add_local_vec_to_global(vector_type& vec, const LocalVector& lvec,
				ConstSmartPtr<DoFDistribution> dd);

	///	adds a local matrix to the global matrix
		void add_local_mat_to_global(matrix_type& mat, const LocalMatrix& lmat,
				ConstSmartPtr<DoFDistribution> dd);

	///	modifies local solution vector for adapted defect computation
		 void modify_LocalSol(LocalVector& vecMod, const LocalVector& lvec, ConstSmartPtr<DoFDistribution> dd){}

	/// sets assembling index
		void set_assembling_index(const size_t assIndex){ m_assemblingIndex = assIndex;}

	///	destructor
		~LocalToGlobalMapperNLGS() {};

	private:
		size_t m_assemblingIndex;
};

/// Nonlinear GaussSeidel-method
/**
 * 	Let L(u) denote a nonlinear functional of n components (l_1,...,l_n).
 * 	Then the basic step of the nonlinear GaussSeidel method is to solve the
 * 	i-th equation
 *
 * 	l_i(u_1^{k+1},...,u_{i-1}^{k+1},u_i,u_{i+1}^{k},...,u_{n}^{k}) = 0
 *
 * 	for u_i and to set u_i^{k+1} = u_i. Here k denotes the iteration-index.
 * 	Note, that the already computed, updated values (.)^{k+1} are used in this
 * 	method.
 * 	Thus, in order to obtain u^{k+1} from u^k, we solve successively the n
 * 	dimensional nonlinear equations for i = 1,...,n. Here this is done
 * 	by a scalar newton step for every i. But every other scalar nonlinear method
 * 	could be applied as well.
 *
 * 	Using a damped version of the nonlinear GaussSeidel method (= nonlinear
 * 	SOR-method) results in the following update of the variables
 *
 * 	u_i^{k+1} = u_i^k + damp * (u_i -u_i^k).
 *
 * References:
 * <ul>
 * <li> J. M. Ortega and W. C. Rheinbolt. Iterative Solution of nonlinear equations in several variables.(1970)
 * </ul>
 *
 *  \tparam 	TDomain			Domain type
 *  \tparam 	TAlgebra		Algebra type
 */
template <typename TDomain, typename TAlgebra>
class NLGaussSeidelSolver
	: public IOperatorInverse<typename TAlgebra::vector_type>,
	  public DebugWritingObject<TAlgebra>
{
	private:
	///	own type
		typedef NLGaussSeidelSolver<TDomain, TAlgebra> this_type;

	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Value type
		typedef typename vector_type::value_type value_type;

	///	Domain type
		typedef TDomain domain_type;

	///	Type of approximation space
		typedef ApproximationSpace<domain_type>	approx_space_type;

	///	Type of geometric base object
		typedef typename domain_traits<TDomain::dim>::grid_base_object grid_base_object;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;

	public:
	///	constructor
		NLGaussSeidelSolver();

	///	constructor
		NLGaussSeidelSolver(SmartPtr<approx_space_type> spApproxSpace,
					SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_approximation_space(SmartPtr<approx_space_type> spApproxSpace)
		{m_spApproxSpace = spApproxSpace;}

		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

		void set_damp(const number damp) {m_damp = damp;}

		void set_constraint(const vector_type& cons) {m_ConsVec = cons; m_bProjectedGS = true;}

	///////////////////////////////////////////////////////////////////////////
	//	OperatorInverse interface methods
	///////////////////////////////////////////////////////////////////////////

	/// This operator inverts the Operator op: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > op);

	/// prepare Operator
		virtual bool prepare(vector_type& u);

	/// apply Operator, i.e. op^{-1}(0) = u
		virtual bool apply(vector_type& u);

	///	returns information about configuration parameters
		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << "NonlinearGaussSeidelSolver( damp = " << m_damp << ")\n";
			ss << " ConvergenceCheck: ";
			if(m_spConvCheck.valid())	ss << ConfigShift(m_spConvCheck->config_string()) << "\n";
			else							ss << " NOT SET!\n";

			return ss.str();
		}

	private:
	///	help functions for debug output
	///	\{
		void write_debug(const vector_type& vec, const char* filename);
		void write_debug(const matrix_type& mat, const char* filename);
	/// \}

	private:
		///	Approximation Space
		SmartPtr<approx_space_type> m_spApproxSpace;

		///	DoF distribution pointer
		ConstSmartPtr<DoFDistribution> m_spLevDD;
		ConstSmartPtr<DoFDistribution> m_spSurfDD;

		/// DoF Distribution used
		GridLevel m_gridLevel;

		SmartPtr<IConvergenceCheck<vector_type> > m_spConvCheck;

		SmartPtr<AssembledOperator<algebra_type> > m_spAssOp;
		SmartPtr<AssembledLinearOperator<algebra_type> > m_spJBlock;
		SmartPtr<IAssemble<TAlgebra> > m_spAss;

		///	damping factor
		number m_damp;

		///	vector describing a constraint
		vector_type m_ConsVec;
		bool m_bProjectedGS;

		///	call counter
		int m_dgbCall;

		/* TODO: hier alle typen von grid_base_object zulassen?
		typedef Attachment<vector<grid_base_object*> > AElemList; 	//attachment type: attachment of ElemDatas
		AElemList m_aElemList;						//the instance of the attachment type
		typedef Grid::VertexAttachmentAccessor<AElemList>	ElemListAccessor;
		ElemListAccessor m_aaElemList;
		///	use this method to make sure that all required attachments are attached
		void attach_attachments()
		{
			typename TDomain::grid_type& grid = *this->domain().grid();
			grid.attach_to_vertices(m_aElemList);
			m_aaElemList.access(grid, m_aElemList);
		}*/

		// TODO: hier alle typen von grid_base_object zulassen?
		typedef std::vector<grid_base_object*> elemList;
		std::vector<elemList> m_vElemList;

		///	selector of elements with contributions to a specific DoF
		Selector m_sel;

		///	LocalToGlobal-Mapper which accounts for index-wise assembling
		LocalToGlobalMapperNLGS<TAlgebra> m_map;
};

}

#include "nl_gauss_seidel_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NL_GAUSS_SEIDEL__NL_GAUSS_SEIDEL_H_ */
