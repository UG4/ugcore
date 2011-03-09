/*
 * gauss_seidel.h
 *
 *  Created on: 14.07.2010
 *      Author: Martin Rupp
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__

#include "lib_algebra/operator/operator_interface.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TAlgebra>
class GSPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		GSPreconditioner() {};

	// 	Clone
		virtual ILinearIterator<vector_type,vector_type>* clone()
		{
			GSPreconditioner<algebra_type>* clone = new GSPreconditioner<algebra_type>();
			return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat) {return true;}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
			// 	set the computed correction to additive
			//	todo: handle parallel case
#ifdef UG_PARALLEL
				c.set_storage_type(PST_ADDITIVE);
#endif
			return gs_step_LL(mat, c, d);
		}
};

template <typename TAlgebra>
class BGSPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		BGSPreconditioner() {};

	// 	Clone
		virtual ILinearIterator<vector_type,vector_type>* clone()
		{
			BGSPreconditioner<algebra_type>* clone = new BGSPreconditioner<algebra_type>();
			return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Backward Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat) {return true;}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
			// 	set the computed correction to additive
			//	todo: handle parallel case
#ifdef UG_PARALLEL
				c.set_storage_type(PST_ADDITIVE);
#endif
			return gs_step_UR(mat, c, d);
		}
};

template <typename TAlgebra>
class SGSPreconditioner : public IPreconditioner<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	//	Constructor
		SGSPreconditioner() {};

	// 	Clone
		virtual ILinearIterator<vector_type,vector_type>* clone()
		{
			SGSPreconditioner<algebra_type>* clone = new SGSPreconditioner<algebra_type>();
			return dynamic_cast<ILinearIterator<vector_type,vector_type>* >(clone);
		}

	protected:
	//	Name of preconditioner
		virtual const char* name() const {return "Symmetric Gauss-Seidel";}

	//	Preprocess routine
		virtual bool preprocess(matrix_type& mat) {return true;}

	//	Postprocess routine
		virtual bool postprocess() {return true;}

	//	Stepping routine
		virtual bool step(matrix_type& mat, vector_type& c, const vector_type& d)
		{
			// 	set the computed correction to additive
			//	todo: handle parallel case
#ifdef UG_PARALLEL
				c.set_storage_type(PST_ADDITIVE);
#endif
			return sgs_step(mat, c, d);
		}
};

} // end namespace ug

#endif // __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__GAUSS_SEIDEL__
