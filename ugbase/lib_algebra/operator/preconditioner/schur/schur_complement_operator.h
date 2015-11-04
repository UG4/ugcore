
#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_COMPLEMENT_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_COMPLEMENT_OPERATOR__


#ifdef UG_PARALLEL

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
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_algebra/operator/debug_writer.h"
#include "schur_complement_inverse_interface.h"
#include "lib_algebra/operator/algebra_debug_writer.h"
#include "pcl/pcl.h"

#include "common/log.h"

#include "slicing.h"


namespace ug{




template <typename TAlgebra>
class SchurComplementOperator
	: public ILinearOperator<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>
{
	public:

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	// 	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;


	public:
	///	constructor
	SchurComplementOperator(SmartPtr<MatrixOperator<matrix_type, vector_type> > Alocal,
							SlicingData::slice_desc_type_vector &sdv)
	: m_spOperator(Alocal),
	  m_slicing(sdv)
	{
		m_op[0][0] = make_sp(new MatrixOperator<matrix_type, vector_type>());
		m_op[0][1] = make_sp(new MatrixOperator<matrix_type, vector_type>());
		m_op[1][0] = make_sp(new MatrixOperator<matrix_type, vector_type>());
		m_op[1][1] = make_sp(new MatrixOperator<matrix_type, vector_type>());
	}

	// destructor
	virtual ~SchurComplementOperator() {};

	///	name of solver
	virtual const char* name() const {return "My local Schur complement Solver";}


	/// implementation of the operator for the solution dependent initialization.
	void init(const vector_type& u) {init();}

	///	initializes the solver for operator A
	virtual void init();

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := S times u'
	virtual void apply(vector_type& f, const vector_type& u);

	///	applies the Schur complement built from matrix operator set via 'set_matrix()'
	/// to 'u' and returns the result 'f := f - S times u'
	virtual void apply_sub(vector_type& f, const vector_type& u);

	//	save current operator
	void set_matrix(SmartPtr<MatrixOperator<matrix_type, vector_type> > A)
	{ m_spOperator = A; }

	///	sets a Dirichlet solver
	void set_dirichlet_solver(SmartPtr<ILinearOperatorInverse<vector_type> > dirichletSolver)
	{ m_spDirichletSolver = dirichletSolver; }


	matrix_type &sub_matrix(int r, int c)
	{return sub_operator(r,c)->get_matrix();}

	SmartPtr<MatrixOperator<matrix_type, vector_type> > sub_operator(int r, int c)
	{return m_op[r][c];}

	size_t sub_size(slice_desc_type type)
	{return m_slicing.get_num_elems(type);}

	const SlicingData &slicing() const
	{return m_slicing;}


	// for debugging: computes schur operator
	void debug_compute_matrix();

	void compute_matrix(matrix_type &schur_matrix, double threshold=0.0);
	virtual void set_debug(SmartPtr<IDebugWriter<algebra_type> > spDebugWriter);

	template<typename T>
	void set_inner_debug(SmartPtr<T> op)
	{
		if(!m_spDebugWriterInner.valid()) return;
		SmartPtr<VectorDebugWritingObject<vector_type> > dvwo
				= op.template cast_dynamic<VectorDebugWritingObject<vector_type> >();
		if(dvwo.valid())
			dvwo->set_debug(m_spDebugWriterInner);
	}

	template<typename T>
	void set_skeleton_debug(SmartPtr<T> op)
	{
		if(!m_spDebugWriterSkeleton.valid()) return;
		SmartPtr<VectorDebugWritingObject<vector_type> > dvwo
				= op.template cast_dynamic<VectorDebugWritingObject<vector_type> >();
		if(dvwo.valid())
			dvwo->set_debug(m_spDebugWriterSkeleton);
	}

protected:
	// 	Operator that is inverted by this Inverse Operator
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_spOperator;

	// slices from matrix
	const SlicingData m_slicing;

	// 	Linear Solver to invert the local Dirichlet problems
	SmartPtr<ILinearOperatorInverse<vector_type> > m_spDirichletSolver;



	// sub matrices/operator (will be set by init)
	SmartPtr<MatrixOperator<matrix_type,vector_type> > m_op[2][2];


	template<int dim> void set_debug_dim();

	SmartPtr<AlgebraDebugWriter<algebra_type> > m_spDebugWriterInner;
	SmartPtr<AlgebraDebugWriter<algebra_type> > m_spDebugWriterSkeleton;
	SmartPtr<IDebugWriter<algebra_type> > m_spDebugWriter;

	int m_applyCnt;

};



} // end namespace ug

#endif /* UG_PARALLEL */

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR_COMPLEMENT_OPERATOR__ */
