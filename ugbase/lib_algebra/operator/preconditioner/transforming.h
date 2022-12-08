/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel, Arne Nägel
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

#ifndef __H__UG__LIB_ALGEBRA__OPERATOR__LINEAR_OPERATOR__TRANSFORMING_ITER__
#define __H__UG__LIB_ALGEBRA__OPERATOR__LINEAR_OPERATOR__TRANSFORMING_ITER__

#include "lib_algebra/operator/interface/linear_iterator.h"
#include "lib_algebra/operator/interface/matrix_operator.h"
#include "lib_algebra/operator/debug_writer.h"

#include "lib_disc/assemble_interface.h"
#include "lib_disc/function_spaces/grid_function.h"

namespace ug{

/**
 * Abstract base class for transforming iterations.
 * Supporting both left- and right transformations
 *
 * Given
 * \f{eqnarray*}{
 *  A = T_L^{-1} \hat{A} T_R
 * \f}
 * this implements a subspace correction based on a defect correction:
 * \f{eqnarray*}{
 *  x := x +  T_R^{-1} {\hat{A}}^{-1} T_L (b-Ax)
 * \f}
 * If inversion is to expensive, we replace may replace this by a (single step) iterative solver.
 *
 * In order
 *
 *
 */
template<typename TAlgebra, typename TDerived>
class ITransformingIteration :
	public ILinearIterator<typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
private:

	/// CRTP operator
	TDerived& derived()
	{ return *static_cast<TDerived*>(this);}

public:

	///	Algebra type
	typedef TAlgebra algebra_type;

	///	Vector type
	typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
	typedef typename algebra_type::matrix_type matrix_type;


	ITransformingIteration() : DebugWritingObject<TAlgebra>() {}
	ITransformingIteration(const ITransformingIteration &ti)
	{ UG_THROW("Do you really want to call a copy constructor of this abstract base class?"); }


	// Begin CRTP

	/// implementation of init for non-linear (CRTP)
	bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
	{ return derived().init(J,u);};

	/// implementation of init for linear (CRTP)
	bool init(SmartPtr<ILinearOperator<vector_type> > L)
	{ return derived().init(L);};

	/// original operator (CRTP)
	SmartPtr<ILinearOperator<vector_type> > original_operator()
	{ return derived().original_operator();}

	/// map: d -> dtilde (CRTP)
	bool transform_defect(vector_type &c, const vector_type& d)
	{ return derived().transform_defect(c,d); }

	/// map: dtilde -> ctilde (CRTP)
	bool apply_transformed(vector_type &c, const vector_type &d)
	{ return derived().apply_transformed(c,d); }

	/// map: ctilde -> c (CRTP)
	bool untransform_correction(vector_type& c, const vector_type& d)
	{ return derived().untransform_correction(c,d); }


	// End CRTP

	///  implementation of apply (final, non-CRTP)
	virtual bool apply(vector_type& c, const vector_type& d)
	{
		return (derived().transform_defect(c,d) &&
				derived().apply_transformed(c,d) &&
				derived().untransform_correction(c,d));
	}


	///  implementation of apply (final, non-CRTP)
	virtual bool apply_update_defect(vector_type &c, vector_type& d)
	{
		// compute correction
		if (!apply(c,d)) return false;

		//	compute updated defect
		original_operator()->apply_sub(d, c);

		return true;
	}



};


/***
 * Implementation of Transforming-/DGS smoothers, following Wittum, 1989.
 *
 *
 * A) For DGS type, assemble
 *
 * \hat A = \begin{pmatrix} Q& W\\
 * 			-div_h & -\triangle_h^N  \end{pmatrix}
 *
 * T_R**{-1} =  \begin{pmatrix} I & grad_h\\
 * 								0 & -Q_h**l  \end{pmatrix}
 *
 * and specify a preconditioner for $\hat A$.
 *
 *
 * B) For Wittum-type transforming smoothers, assemble
 *
 * \hat A = \begin{pmatrix} Q& W\\
 * 			-div_h & -\triangle_h^N  \end{pmatrix}
 *
 * T_R**{-1} =  \begin{pmatrix} A & B \\
 * 								0 & I \end{pmatrix}
 *
 *
 * and specify a preconditioner for both $\hat A$ and T_R.
 *
 *  */
template <typename TDomain, typename TAlgebra>
class AssembledTransformingSmoother :
	//public ILinearIterator<typename TAlgebra::vector_type>,
	//public DebugWritingObject<TAlgebra>
	public ITransformingIteration<TAlgebra, AssembledTransformingSmoother<TDomain,TAlgebra> >
{
	public:
	///	Domain
		typedef TDomain domain_type;

	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename algebra_type::vector_type vector_type;

	///	Matrix type
		typedef typename algebra_type::matrix_type matrix_type;

		using DebugWritingObject<TAlgebra>::write_debug;
		using ILinearIterator<typename TAlgebra::vector_type>::damping;

	public:
	/// constructor setting approximation space
		AssembledTransformingSmoother(SmartPtr<IAssemble<TAlgebra> > spAuxSystemAss,
		                              SmartPtr<ILinearIterator<vector_type> > spAuxSmoother,
		                              SmartPtr<IAssemble<TAlgebra> > spRightTrafoAss,
		                              SmartPtr<ILinearIterator<vector_type> > spRightTrafoSmoother=SPNULL)
			: m_spAuxSystemAss(spAuxSystemAss),
			  m_spAuxMatrixOperator(new MatrixOperator<matrix_type, vector_type>),
			  m_spAuxSmoother(spAuxSmoother),
			  m_spRightTrafoAss(spRightTrafoAss),
			  m_spRightTrafoMat(new MatrixOperator<matrix_type, vector_type>),
			  m_spRightTrafoSmoother(spRightTrafoSmoother)
		{}

	///	name
		virtual const char* name() const {return "Assembled Transform Smoother";}

	///	returns if parallel solving is supported
		virtual bool supports_parallel() const
		{return (m_spAuxSmoother.valid()) ? m_spAuxSmoother->supports_parallel() : false;}

		// Begin CRTP functions
		bool init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u);

		/// Since we need grid information, linear operators are not supported...
		bool init(SmartPtr<ILinearOperator<vector_type> > L)
		{ UG_THROW("Not Implemented!!"); }

		bool transform_defect(vector_type& c, const vector_type& d) {return true;}
		bool apply_transformed(vector_type &c, const vector_type &d);
		bool untransform_correction(vector_type& c, const vector_type& d);


	///	Clone
		SmartPtr<ILinearIterator<vector_type> > clone();




		SmartPtr<ILinearOperator<vector_type> > original_operator()
		{return m_spOriginalSystemOp;}

	protected:

	///	matrix for original system
		SmartPtr<ILinearOperator<vector_type> > m_spOriginalSystemOp;

	/// auxiliary correction vector
		SmartPtr<vector_type> m_spAuxCorrection;

	///	assembling the transformed system
		SmartPtr<IAssemble<TAlgebra> > m_spAuxSystemAss;

	///	matrix (operator) of transformed system
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spAuxMatrixOperator;

	///	smoother used on transformed system
		SmartPtr<ILinearIterator<vector_type> > m_spAuxSmoother;


	///	assembling the right-transformation
		SmartPtr<IAssemble<TAlgebra> > m_spRightTrafoAss;

	///	matrix (operator) of right-transformation
		SmartPtr<MatrixOperator<matrix_type, vector_type> > m_spRightTrafoMat;

	///	smoother used on transformed system
		SmartPtr<ILinearIterator<vector_type> > m_spRightTrafoSmoother;




};



template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
init(SmartPtr<ILinearOperator<vector_type> > J, const vector_type& u)
{
	//m_spOriginalSystemOp = J;

	GridLevel gridLevel;

	const GridFunction<TDomain, TAlgebra>* pSol =
		dynamic_cast<const GridFunction<TDomain, TAlgebra>*>(&u);
	if(pSol){
		gridLevel = pSol->dof_distribution()->grid_level();
	}

	// init aux operator
	try{
		m_spAuxSystemAss->assemble_jacobian(*m_spAuxMatrixOperator, u, gridLevel);
	}
	UG_CATCH_THROW("AssembledTransformingSmoother: "
					" Cannot assemble transformed system matrix.");

	write_debug(*m_spAuxMatrixOperator, "TrafoSystem");

	if(!m_spAuxSmoother->init(m_spAuxMatrixOperator, u))
		UG_THROW("AssembledTransformingSmoother: "
				" Cannot init smoother for transformed system matrix.");


	// init right transform
	try{
			m_spRightTrafoAss->assemble_jacobian(*m_spRightTrafoMat, u, gridLevel);
		}
		UG_CATCH_THROW("AssembledTransformingSmoother: "
				" Cannot assemble right-transformation matrix.");

	write_debug(*m_spRightTrafoMat, "RightTrafo");


	return true;
};




/**
 *
 * */
template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
apply_transformed(vector_type &c, const vector_type& d)
{
	/**
//	temporary vector for defect
	vector_type dTmp; dTmp.resize(d.size());

//	copy defect
	dTmp = d;

//	work on copy
	return apply_update_defect(c, dTmp);

	*/

	// obtain aux correction vector
	m_spAuxCorrection = c.clone_without_values();

	//	apply smoother of transformed system
	if(!m_spAuxSmoother->apply(*m_spAuxCorrection, d))
	{
		UG_LOG("AssembledTransformingSmoother: Smoother applied incorrectly.\n");
		return false;
	}

	return true;
}

template <typename TDomain, typename TAlgebra>
bool AssembledTransformingSmoother<TDomain, TAlgebra>::
untransform_correction(vector_type& c, const vector_type& d)
{

	//	apply right-transformation
	if (m_spRightTrafoSmoother.valid())
	{
		// e.g., for Wittum-type smoothers
		if(!m_spRightTrafoSmoother->apply(*m_spAuxCorrection, d))
		{
			UG_LOG("AssembledTransformingSmoother: Right Trafo smoother failed!\n");
			return false;
		}
	}
	else
	{
		// e.g., for Brandt-Dinar-DGS-type smoothers
		m_spRightTrafoMat->apply(c, *m_spAuxCorrection);
	}


	// release auy correction vecto
	m_spAuxCorrection = SPNULL;

	#ifdef UG_PARALLEL
		c.change_storage_type(PST_CONSISTENT);
	#endif

	const number damp = this->damping()->damping();
	c *= damp;

	return true;
}


template <typename TDomain, typename TAlgebra>
SmartPtr<ILinearIterator<typename TAlgebra::vector_type> >
AssembledTransformingSmoother<TDomain, TAlgebra>::
clone()
{
	SmartPtr<AssembledTransformingSmoother<TDomain, TAlgebra> > clone(
		new AssembledTransformingSmoother<TDomain, TAlgebra>(
										m_spAuxSystemAss, m_spAuxSmoother,
										m_spRightTrafoAss, m_spRightTrafoSmoother));

	return clone;
}


} // end namespace ug

#endif /* __H__UG__PLUGINS__NAVIER_STOKES__INCOMPRESSIBLE__TRANSFORMING_SMOOTHER__ */
