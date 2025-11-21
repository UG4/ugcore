/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__

// extern libraries
#include <cassert>
#include <map>

// other ug4 modules
#include "common/math/ugmath.h"

// library intern headers
#include "local_finite_element_id.h"
#include "lib_disc/local_finite_element/local_dof_set.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
//	Interface for local shape function sets
////////////////////////////////////////////////////////////////////////////////

/// \ingroup lib_disc_local_finite_elements
/// @{

/// virtual base class for local shape function sets
/**
 * This class is a base class for the supply of local shape functions on finite
 * elements. The class provides evaluation of the shape functions and the
 * gradients at arbitrary points in the interior of a reference element.
 *
 * \tparam 	TDim	Reference Element Dimension
 * \tparam	TShape	type of Range of Shape Functions
 * \tparam	TGrad	type of gradient of shape functions
 */
template <	int TDim,
			typename TShape = number,
			typename TGrad = MathVector<TDim> >
class LocalShapeFunctionSet : public DimLocalDoFSet<TDim>
{
	public:
	///	Dimension, where shape functions are defined
		static constexpr int dim = TDim;

	///	Domain position type
	using position_type = MathVector<dim>;

	///	Shape type
	using shape_type = TShape;

	///	Gradient type
	using grad_type = TGrad;

	public:
	///	returns if space constructs continuous functions
		virtual bool continuous() const = 0;

	/// evaluates the shape function
	/**
	 * This function returns the value of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return shape function value at point
	 */
		virtual shape_type shape(size_t i, const MathVector<dim>& x) const = 0;

	/// evaluates the shape function
	/**
	 * This function returns the value of Shape Function i at
	 * an element-local evaluation point.
	 * \param[out]	shape	the shape function
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
		virtual void shape(shape_type& shape, size_t i, const MathVector<dim>& x) const = 0;

	/// returns all shape functions evaluated at a point
	/**
	 * This function returns the values of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	vShape	Vector of Shapes
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
	///	\{
		virtual void shapes(shape_type* vShape, const MathVector<dim>& x) const = 0;
		virtual void shapes(std::vector<shape_type>& vShape, const MathVector<dim>& x) const = 0;
	///	\}

	/// returns all shape functions evaluated at several point
	/**
	 * This function returns the values of all Shape Functions at
	 * several element-local evaluation point in an array.
	 * \param[out]	vvShape	Vector of Shapes
	 * \param[in]	vLocPos	Vector of Position on reference element
	 */
		virtual void shapes(std::vector<std::vector<shape_type> >& vvShape,
		                    const std::vector<MathVector<dim> >& vLocPos) const = 0;

	/// evaluates the gradient of the shape function
	/** This function returns the gradient of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return gradient at point
	 */
		virtual void grad(grad_type& g, size_t i, const MathVector<dim>& x) const = 0;

	/// returns all gradients evaluated at a point
	/**
	 * This function returns the gradients of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	vGrad	Vector of gradients
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
	///	\{
		virtual void grads(grad_type* vGrad, const MathVector<dim>& x) const = 0;
		virtual void grads(std::vector<grad_type>& vGrad, const MathVector<dim>& x) const = 0;
	///	\}

	/// returns all gradients evaluated at a several points
	/**
	 * This function returns the gradients of all Shape Functions at
	 * several element-local evaluation point in an array.
	 * \param[out]	vvGrad	Vector of gradients
	 * \param[in]	vLocPos	Vector of Position on reference element
	 */
		virtual void grads(std::vector<std::vector<grad_type> >& vvGrad,
		                   const std::vector<MathVector<dim> >& vLocPos) const = 0;

	///	virtual destructor
	~LocalShapeFunctionSet() override = default;
};

/// @}

////////////////////////////////////////////////////////////////////////////////
//	Common base class for local shape function sets to ease implementation
////////////////////////////////////////////////////////////////////////////////

/// static interface for trial spaces
template <typename TImpl, int TDim,
		  typename TShape = number,
		  typename TGrad = MathVector<TDim> >
class BaseLSFS
{
	public:
	///	type of implementation
		using ImplType = TImpl;

	///	dimension of reference element
		static constexpr int dim = TDim;

	///	Shape type
		using shape_type = TShape;

	///	Gradient type
		using grad_type = TGrad;

	//////////////////////////////////////////
	//	methods implemented by derived class
	//////////////////////////////////////////

	public:
	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline shape_type shape(size_t i, const MathVector<dim>& x) const
		{
			return getImpl().shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, size_t i, const MathVector<dim>& x) const
		{
			 getImpl().grad(g, i, x); return;
		}

	//////////////////////////////////////////
	//	methods generated generically
	//////////////////////////////////////////

	public:
	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void shape(shape_type& sh, size_t i, const MathVector<dim>& x) const
		{
			sh = shape(i, x);
		}


	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		inline void shapes(shape_type* vShape, const MathVector<dim>& x) const
		{
			for(size_t sh = 0; sh < getImpl().num_sh(); ++sh)
				vShape[sh] = shape(sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		inline void shapes(std::vector<shape_type>& vShape, const MathVector<dim>& x) const
		{
			vShape.resize(getImpl().num_sh()); shapes(&vShape[0], x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		inline void shapes(std::vector<std::vector<shape_type> >& vvShape,
		                    const std::vector<MathVector<dim> >& vLocPos) const
		{
			vvShape.resize(vLocPos.size());
			for(size_t ip = 0; ip < vLocPos.size(); ++ip)
				shapes(vvShape[ip], vLocPos[ip]);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		inline void grads(grad_type* vGrad, const MathVector<dim>& x) const
		{
			for(size_t sh = 0; sh < getImpl().num_sh(); ++sh)
				grad(vGrad[sh], sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		inline void grads(std::vector<grad_type>& vGrad, const MathVector<dim>& x) const
		{
			vGrad.resize(getImpl().num_sh()); grads(&vGrad[0], x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		inline void grads(std::vector<std::vector<grad_type> >& vvGrad,
		                   const std::vector<MathVector<dim> >& vLocPos) const
		{
			vvGrad.resize(vLocPos.size());
			for(size_t ip = 0; ip < vLocPos.size(); ++ip)
				grads(vvGrad[ip], vLocPos[ip]);
		}

	protected:
	///	access to implementation
		ImplType& getImpl() {return static_cast<ImplType&>(*this);}

	///	const access to implementation
		const ImplType& getImpl() const {return static_cast<const ImplType&>(*this);}

};

} // namespace ug

#endif