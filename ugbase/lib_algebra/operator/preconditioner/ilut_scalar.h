/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT_SCALAR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__ILUT_SCALAR__

#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/preconditioner.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "common/progress.h"
#include "common/util/ostream_util.h"

#include "lib_algebra/algebra_common/vector_util.h"
#include "lib_algebra/operator/linear_solver/linear_solver.h"
#include "lib_algebra/cpu_algebra_types.h"
#include "ilut.h"

namespace ug{

template <typename TAlgebra>
class ILUTScalarPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
		using algebra_type = TAlgebra;
		using vector_type = typename TAlgebra::vector_type;
		using matrix_type = typename TAlgebra::matrix_type;
		using matrix_operator_type = typename IPreconditioner<TAlgebra>::matrix_operator_type;

	private:
		using block_type = typename matrix_type::value_type;
		using base_type = IPreconditioner<TAlgebra>;

	public:
	//	Constructor
		ILUTScalarPreconditioner(double eps=1e-6) :
			m_eps(eps), m_info(false), m_show_progress(true), m_bSort(true)
		{};

	/// clone constructor
		ILUTScalarPreconditioner( const ILUTScalarPreconditioner &parent )
			: base_type(parent)
		{
			m_eps = parent.m_eps;
			set_info(parent.m_info);
			set_show_progress(parent.m_show_progress);
			set_sort(parent.m_bSort);
		}

	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone() override {
			return make_sp(new ILUTScalarPreconditioner(*this));
		}

	// 	Destructor
		~ILUTScalarPreconditioner() override = default;

	///	returns if parallel solving is supported
		bool supports_parallel() const override {return true;}

	///	sets threshold for incomplete LU factorisation (added 01122010ih)
		void set_threshold(number thresh)
		{
			m_eps = thresh;
		}
		
	///	sets storage information output to true or false
		void set_info(bool info)
		{
			m_info = info;
		}
		
	///	sets the indication of the progress to true or false
		void set_show_progress(bool b)
		{
			m_show_progress = b;
		}
		
		void set_sort(bool b)
		{
			m_bSort = b;
		}

	public:
		void preprocess(const matrix_type &M)
		{
#ifdef 	UG_PARALLEL
			matrix_type A;
			A = M;

			MatAddSlaveRowsToMasterRowOverlap0(A);

		//	set zero on slaves
			std::vector<IndexLayout::Element> vIndex;
			CollectUniqueElements(vIndex, M.layouts()->slave());
			SetDirichletRow(A, vIndex);
#else
			const matrix_type &A = M;
#endif

			STATIC_ASSERT(matrix_type::rows_sorted, Matrix_has_to_have_sorted_rows);

			ilut = make_sp(new ILUTPreconditioner<CPUAlgebra>(m_eps));
			ilut->set_info(m_info);
			ilut->set_show_progress(m_show_progress);
			ilut->set_sort(m_bSort);

			mo = make_sp(new MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type>);
			CPUAlgebra::matrix_type &mat = mo->get_matrix();

			#ifdef UG_PARALLEL
				mat.set_storage_type(PST_ADDITIVE);
				mat.set_layouts(CreateLocalAlgebraLayouts());
			#endif

			m_size = GetDoubleSparseFromBlockSparse(mat, A);

			SmartPtr convCheck(new StdConvCheck<CPUAlgebra::vector_type>);
			convCheck->set_maximum_steps(100);
			convCheck->set_minimum_defect(1e-50);
			convCheck->set_reduction(1e-20);
			convCheck->set_verbose(false);

			linearSolver.set_preconditioner(ilut);
			linearSolver.set_convergence_check(convCheck);
			linearSolver.init(mo);
		}
	protected:
	//	Name of preconditioner
		const char* name() const override {return "ILUTScalar";}


	//	Preprocess routine
		bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp) override {
			matrix_type &A = *pOp;
			preprocess(A);

			return true;
		}

		void get_vector(CPUAlgebra::vector_type &v_scalar, const vector_type& v)
		{
			for(size_t i=0, k=0; i<v.size(); i++)
			{
				for(size_t j=0; j<GetSize(v[i]); j++)
					v_scalar[k++] = BlockRef(v[i],j);
			}
		}

		void set_vector(CPUAlgebra::vector_type &v_scalar, vector_type& v)
		{
			for(size_t i=0, k=0; i<v.size(); i++)
			{
				for(size_t j=0; j<GetSize(v[i]); j++)
					BlockRef(v[i],j) = v_scalar[k++];
			}
		}

	//	Stepping routine
		bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d) override {
#ifdef UG_PARALLEL
			SmartPtr<vector_type> spDtmp = d.clone();
			spDtmp->change_storage_type(PST_UNIQUE);
			bool b = apply_double(c, *spDtmp);

			c.set_storage_type(PST_ADDITIVE);
			c.change_storage_type(PST_CONSISTENT);
			return b;
#else
			return apply_double(c, d);
#endif
			return true;
		}

		bool apply_double(vector_type &c, const vector_type &d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);
			m_c.set_layouts(CreateLocalAlgebraLayouts());
			m_d.set_layouts(CreateLocalAlgebraLayouts());
#endif
			get_vector(m_d, d);
			m_c.set(0.0);

			//ApplyLinearSolver(mo, m_u, m_b, linearSolver);
			ilut->apply(m_c, m_d);
			//ilut->apply(m_c, m_d);

			set_vector(m_c, c);
			return true;
		}



	public:
		bool solve(vector_type &c, const vector_type &d)
		{
			m_c.resize(m_size);
			m_d.resize(m_size);

#ifdef UG_PARALLEL
			m_d.set_storage_type(PST_ADDITIVE);
			m_c.set_storage_type(PST_CONSISTENT);
			m_c.set_layouts(CreateLocalAlgebraLayouts());
			m_d.set_layouts(CreateLocalAlgebraLayouts());
#endif
			get_vector(m_d, d);
			m_c.set(0.0);

			//ApplyLinearSolver(mo, m_u, m_b, linearSolver);
			linearSolver.apply_return_defect(m_c, m_d);
			//ilut->apply(m_c, m_d);

			set_vector(m_c, c);
			return true;
		}

	public:
		std::string config_string() const override {
			std::stringstream ss ; ss << "ILUTScalar(threshold = " << m_eps << ", sort = " << (m_bSort?"true":"false") << ")";
			if(m_eps == 0.0) ss << " = Sparse LU";
			return ss.str();
		}

	protected:
	//	Postprocess routine
		bool postprocess() override {return true;}

protected:
	SmartPtr<ILUTPreconditioner<CPUAlgebra> > ilut;
	SmartPtr<MatrixOperator<CPUAlgebra::matrix_type, CPUAlgebra::vector_type> > mo;
	LinearSolver<CPUAlgebra::vector_type> linearSolver;
	CPUAlgebra::vector_type m_c, m_d;
	double m_eps;
	bool m_info, m_show_progress, m_bSort;
	size_t m_size;
};

} // end namespace ug

#endif
