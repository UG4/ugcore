/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__MYCOLOR_GS_PRECOND__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__MYCOLOR_GS_PRECOND__




#include <iostream>
#include <sstream>
#include <string>
#include <set>

#include "lib_algebra/operator/interface/preconditioner.h"
#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/interface/preconditioned_linear_operator_inverse.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/interface/matrix_operator_inverse.h"
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallelization.h"
#include "pcl/pcl.h"
#endif

#include "lib_algebra/operator/debug_writer.h"

#include "common/log.h"

#include "lib_algebra/adapter/slicing.h"


namespace ug{

extern DebugID SchurDebug;



namespace ColorGS {


typedef std::vector<int> ColoringVector;
static const size_t MAX_COLORS = 256;

typedef SlicingData<ColoringVector, MAX_COLORS> SlicingData;


/// operator implementation of the DD Schur complement solver
/**
 * This operator implements a Schur complement solver */
//template <typename TAlgebra>
//class SchurSolver : public IMatrixOperatorInverse<	typename TAlgebra::matrix_type,
//													typename TAlgebra::vector_type>,
//	public DebugWritingObject<TAlgebra>

template <typename TAlgebra>
class ColorPrecond: public IPreconditioner<TAlgebra>
{
	public:
	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Base type
		typedef IPreconditioner<TAlgebra> base_type;

	protected:
		//using base_type::set_debug;
		using base_type::write_debug;
		using base_type::debug_writer;

	public:
	///	constructor
		ColorPrecond();


		ColorPrecond(const ColorPrecond<TAlgebra> &parent) : base_type(parent)
		{
			// m_spDirichletSolver = parent.m_spDirichletSolver;
		}


		void set_slicing(SlicingData &s)
		{ m_slicing = s; }

		const SlicingData &slicing() const
		{return m_slicing;}


		virtual SmartPtr<ILinearIterator<vector_type> > clone()
		{
			return make_sp(new ColorPrecond<algebra_type>(*this));
		}

	protected:
	///	name of solver
		virtual const char* name() const {return "Schur complement";}

		//	Preprocess routine
		virtual bool preprocess(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp);

		//	Stepping routine
		virtual bool step(SmartPtr<MatrixOperator<matrix_type, vector_type> > pOp, vector_type& c, const vector_type& d);

		//	Postprocess routine
		virtual bool postprocess();

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{
			return true;
		}

	public:
		///	set debug output
		void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter)
		{
			base_type::set_debug(spDebugWriter);
		}


		virtual std::string config_string() const
		{
			std::stringstream ss; ss << name() << "\n";
			return ss.str();
		}



	private:

		bool check_requirements();

		void create_aux_vectors(const vector_type& d);
		void clear_aux_vectors();

		void create_diag_vectors(const matrix_type& A);
		void clear_diag_vectors();


	protected:

		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_pA;

		// Slicing data.
		SlicingData m_slicing;

		// temporary vectors for correction/defect
		SmartPtr<vector_type> m_aux_rhs[MAX_COLORS];
		SmartPtr<vector_type> m_aux_sol[MAX_COLORS];
		std::vector<typename matrix_type::value_type> m_diag[MAX_COLORS];

};

}  // namespace ColorGS


} // namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__COLOR_GS_PRECOND__ */
