/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/algebra_common/core_smoothers.h"
#include "lib_algebra/algebra_common/sparsematrix_util.h"
#ifdef UG_PARALLEL
	//#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_algebra/parallelization/matrix_overlap.h"
	#include "lib_algebra/parallelization/parallel_matrix_overlap.h"
#endif

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/algebra_common/permutation_util.h"

namespace ug {

template<typename TAlgebra>
class GaussSeidelBase : public IPreconditioner<TAlgebra>
{
	public:
		using algebra_type = TAlgebra;
		using vector_type = typename TAlgebra::vector_type;
		using matrix_type = typename TAlgebra::matrix_type;
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;
		using base_type = IPreconditioner<TAlgebra>;

	///	Ordering type
		using ordering_container_type = std::vector<size_t>;
		using ordering_algo_type = IOrderingAlgorithm<TAlgebra, ordering_container_type>;

	protected:
		using base_type::set_debug;
		using base_type::debug_writer;
		using base_type::write_debug;

	public:
	//	Constructor
		GaussSeidelBase() :
			m_relax(1.0),
			m_bConsistentInterfaces(false),
			m_useOverlap(false) {};

	/// clone constructor
		GaussSeidelBase( const GaussSeidelBase &parent )
			: base_type(parent),
			  m_bConsistentInterfaces(parent.m_bConsistentInterfaces),
			  m_useOverlap(parent.m_useOverlap),
			  m_spOrderingAlgo(parent.m_spOrderingAlgo)
		{
			set_sor_relax(parent.m_relax);
		}

	///	set relaxation parameter to define a SOR-method
		void set_sor_relax(number relaxFactor){ m_relax = relaxFactor;}

	///	activates the new parallelization approach (disabled by default)
		void enable_consistent_interfaces(bool enable) {m_bConsistentInterfaces = enable;}

		void enable_overlap (bool enable) {m_useOverlap = enable;}

	/// 	sets an ordering algorithm
		void set_ordering_algorithm(SmartPtr<ordering_algo_type> ordering_algo){
			m_spOrderingAlgo = ordering_algo;
		}

		[[nodiscard]] const char* name() const override = 0;
	protected:
		[[nodiscard]] bool supports_parallel() const override {return true;}

	//	Preprocess routine
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override {
			PROFILE_BEGIN_GROUP(GaussSeidel_preprocess, "algebra gaussseidel");
			matrix_type *pA;
#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				if(m_useOverlap){
					m_A = *pOp;
					CreateOverlap(m_A);
					m_oD.set_layouts(m_A.layouts());
					m_oC.set_layouts(m_A.layouts());
					m_oD.resize(m_A.num_rows(), false);
					m_oC.resize(m_A.num_rows(), false);
				}
				else if (m_bConsistentInterfaces)
				{
					m_A = *pOp;
					MatMakeConsistentOverlap0(m_A);
				}
				else
				{
					//	copy original matrix
					MakeConsistent(*pOp, m_A);
					//	set zero on slaves
					std::vector<IndexLayout::Element> vIndex;
					CollectUniqueElements(vIndex,  m_A.layouts()->slave());
					SetDirichletRow(m_A, vIndex);
				}
				pA = &m_A;
			}
			else
				pA = &(*pOp);
#else
			pA = &(*pOp);
#endif
			THROW_IF_NOT_EQUAL(pA->num_rows(), pA->num_cols());
//			UG_ASSERT(CheckDiagonalInvertible(A), "GS: A has noninvertible diagonal");
			UG_COND_THROW(CheckDiagonalInvertible(*pA) == false, name() << ": A has noninvertible diagonal");

			return true;
		}

	//	Postprocess routine
		bool postprocess() override {return true;}

		virtual void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax) = 0;

	//	Stepping routine
		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override {
			PROFILE_BEGIN_GROUP(GaussSeidel_step, "algebra gaussseidel");

#ifdef UG_PARALLEL
			if(pcl::NumProcs() > 1)
			{
				if(m_useOverlap){
					for(size_t i = 0; i < d.size(); ++i)
						m_oD[i] = d[i];
					for(size_t i = d.size(); i < m_oD.size(); ++i)
						m_oD[i] = 0;
					m_oD.set_storage_type(PST_ADDITIVE);
					m_oD.change_storage_type(PST_CONSISTENT);

					step(m_A, m_oC, m_oD, m_relax);

					for(size_t i = 0; i < c.size(); ++i)
						c[i] = m_oC[i];
					SetLayoutValues(&c, c.layouts()->slave(), typename TAlgebra::vector_type::value_type(0));
					c.set_storage_type(PST_UNIQUE);
				}
				else if (m_bConsistentInterfaces)
				{
					// make defect consistent
					SmartPtr<vector_type> spDtmp = d.clone();
					spDtmp->change_storage_type(PST_CONSISTENT);

					THROW_IF_NOT_EQUAL_3(c.size(), spDtmp->size(), m_A.num_rows());
					step(m_A, c, *spDtmp, m_relax);

					// declare c unique to enforce that only master correction is used
					// when it is made consistent below
					c.set_storage_type(PST_UNIQUE);

					//UG_COND_THROW(!d.has_storage_type(PST_ADDITIVE), "Additive or unique defect expected.");
					//step(m_A, c, d, m_relax);
					//c.set_storage_type(PST_ADDITIVE);
				}
				else
				{
					// todo: do not clone every time
				//	make defect unique
					SmartPtr<vector_type> spDtmp = d.clone();
					spDtmp->change_storage_type(PST_UNIQUE);

					THROW_IF_NOT_EQUAL_3(c.size(), spDtmp->size(), m_A.num_rows());
					step(m_A, c, *spDtmp, m_relax);
					c.set_storage_type(PST_UNIQUE);
				}

				// make correction consistent
				c.change_storage_type(PST_CONSISTENT);

				return true;
			}
			else
#endif
			{
				matrix_type &A = *pOp;
				THROW_IF_NOT_EQUAL_4(c.size(), d.size(), A.num_rows(), A.num_cols());

				step(A, c, d, m_relax);
#ifdef UG_PARALLEL
				c.set_storage_type(PST_CONSISTENT);
#endif
				return true;
			}
		}

	protected:
#ifdef UG_PARALLEL
		matrix_type m_A;
#endif
	protected:
	///	relaxation parameter
		number m_relax;

	///	for overlaps only
		vector_type m_oD;
		vector_type m_oC;

		bool m_bConsistentInterfaces;
		bool m_useOverlap;


	/// for ordering algorithms
		SmartPtr<ordering_algo_type> m_spOrderingAlgo;
#ifdef NOT_YET
		ordering_container_type m_ordering, m_old_ordering;
		std::vector<size_t> m_newIndex, m_oldIndex;
		bool m_bSortIsIdentity;
#endif
};

/// Gauss-Seidel preconditioner for the 'forward' ordering of the dofs
/**
 * This class implements the Gauss-Seidel preconditioner (and smoother) for the
 * 'forward' ordering of the dofs. When a relaxation parameter is set by the method
 * 'set_sor_relax', the resulting preconditioner is better known as (forward) 'SOR'-method.
 * References:
 * <ul>
 * <li> W. Hackbusch. Iterative solution of large sparse systems of equations. New York: Springer, 1994
 * </ul>
 *
 * \tparam	TAlgebra	Algebra type
 */
template <typename TAlgebra>
class GaussSeidel : public GaussSeidelBase<TAlgebra>
{
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;
	using base_type = GaussSeidelBase<TAlgebra>;

public:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "Gauss-Seidel";}

	/// constructor
		GaussSeidel() = default;

	/// clone constructor
		GaussSeidel( const GaussSeidel &parent )
			: base_type(parent)
		{	}

		~GaussSeidel() override = default;

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new GaussSeidel(*this));
		}

	//	Stepping routine
		void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax) override {
			gs_step_LL(A, c, d, relax);
		}
};

/// Gauss-Seidel preconditioner for the 'backward' ordering of the dofs
/**
 * This class implements the Gauss-Seidel preconditioner (and smoother) for the
 * 'backward' ordering of the dofs. When a relaxation parameter is set by the method
 * 'set_sor_relax', the resulting preconditioner is better known as backward 'SOR'-method.
 * References:
 * <ul>
 * <li> W. Hackbusch. Iterative solution of large sparse systems of equations. New York: Springer, 1994
 * </ul>
 *
 * \tparam	TAlgebra	Algebra type
 */
template <typename TAlgebra>
class BackwardGaussSeidel : public GaussSeidelBase<TAlgebra>
{
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;
	using base_type = GaussSeidelBase<TAlgebra>;

public:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "Backward Gauss-Seidel";}

	/// constructor
		BackwardGaussSeidel() = default;

	/// clone constructor
		BackwardGaussSeidel( const BackwardGaussSeidel &parent )
			: base_type(parent)
		{	}

		~BackwardGaussSeidel() override = default;

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new BackwardGaussSeidel(*this));
		}

	//	Stepping routine
		void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax) override {
			gs_step_UR(A, c, d, relax);
		}
};


template <typename TAlgebra>
class SymmetricGaussSeidel : public GaussSeidelBase<TAlgebra>
{
	using algebra_type = TAlgebra;
	using vector_type = typename TAlgebra::vector_type;
	using matrix_type = typename TAlgebra::matrix_type;
	using base_type = GaussSeidelBase<TAlgebra>;

public:
	//	Name of preconditioner
		[[nodiscard]] const char* name() const override {return "Symmetric Gauss-Seidel";}

	/// constructor
		SymmetricGaussSeidel() = default;

	/// clone constructor
		SymmetricGaussSeidel( const SymmetricGaussSeidel &parent )
			: base_type(parent)
		{	}

		~SymmetricGaussSeidel() override = default;
	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new SymmetricGaussSeidel(*this));
		}

	//	Stepping routine
		void step(const matrix_type &A, vector_type &c, const vector_type &d, const number relax) override {
			sgs_step(A, c, d, relax);
		}
};

} // end namespace ug

#endif
