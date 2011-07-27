/*
 * local_shape_function_set.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__
#define __H__LIBDISCRETIZATION__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__

// extern libraries
#include <cassert>
#include <map>

// other ug4 modules
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_discretization/reference_element/reference_element.h"
#include "local_finite_element_id.h"

namespace ug {

// Exception
struct UG_ERROR_InvalidShapeFunctionIndex
{
	UG_ERROR_InvalidShapeFunctionIndex(size_t i)
		: index(i)
	{}
	size_t index;
};

/// \ingroup lib_disc_local_finite_elements
/// @{

/// base class for local shape function sets
/**
 * This class is a base class for the supply of local shape functions on finite
 * elements. The class provides evaluation of the shape functions and the
 * gradients at arbitrary points in the interior of a reference element.
 *
 * \tparam 	tDim	Reference Element Dimension
 */
template <int TDim>
class DimLocalShapeFunctionSet
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = TDim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	public:
	///	Number of DoFs (shapes) on finite element
		virtual size_t num_sh() const = 0;

	///	local position of DoF i
	/**
	 * This function returns the local position of a DoF if possible.
	 * \param[in] 	i		number of DoF
	 * \param[out]	pos		Position of DoF
	 * \retval		true 	if position exists
	 * \retval		false 	if no meaningful position available
	 */
		virtual bool position(size_t i, position_type& pos) const = 0;

	/// evaluates the shape function
	/**
	 * This function returns the value of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return shape function value at point
	 */
		virtual shape_type shape(size_t i, const position_type& x) const = 0;

	/// returns all shape functions evaluated at a point
	/**
	 * This function returns the values of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	sOut	Vector of Shapes
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
		virtual void shapes(shape_type* sOut, const position_type& x) const = 0;

	/// evaluates the gradient of the shape function
	/** This function returns the gradient of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return gradient at point
	 */
		virtual grad_type grad(size_t i, const position_type& x) const = 0;

	/// returns all gradients evaluated at a point
	/**
	 * This function returns the gradients of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	gOut	Vector of gradients
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
		virtual void grads(grad_type* gOut, const position_type& x) const = 0;

	///	virtual destructor
		virtual ~DimLocalShapeFunctionSet() {};
};

template <typename TRefElem>
class LocalShapeFunctionSet
	: public DimLocalShapeFunctionSet<TRefElem::dim>
{
	public:
	///	Reference Element type
		typedef TRefElem reference_element_type;
};

/// @}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// wrapper class implementing the LocalShapeFunctionSet interface
/**
 * This class wrappes a class passed by the template argument into the
 * virtual ILocalShapeFunctionSet interface and makes it thus usable in that
 * context on the price of virtual functions.
 *
 * \tparam 	TImpl		Implementation of a Local Shape Function Set
 */
template <typename TImpl>
class LocalShapeFunctionSetWrapper
	: public LocalShapeFunctionSet<typename TImpl::reference_element_type>,
	  public TImpl
{
	/// Implementation
		typedef TImpl ImplType;

	public:
	///	Reference Element type
		typedef typename ImplType::reference_element_type reference_element_type;

	///	Order of Shape functions
		static const size_t order = ImplType::order;

	///	Dimension, where shape functions are defined
		static const int dim = ImplType::dim;

	///	Domain position type
		typedef typename ImplType::position_type position_type;

	///	Shape type
		typedef typename ImplType::shape_type shape_type;

	///	Gradient type
		typedef typename ImplType::grad_type grad_type;

	/// Number of shape functions
		static const size_t nsh = ImplType::nsh;

	public:
	///	constructor
		LocalShapeFunctionSetWrapper(){}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		virtual size_t num_sh() const {return nsh;}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		virtual bool position(size_t i, position_type& pos) const
		{
			return ImplType::position(i, pos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		virtual shape_type shape(size_t i, const position_type& x) const
		{
			return ImplType::shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* sOut, const position_type& x) const
		{
			ImplType::shapes(sOut, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual grad_type grad(size_t i, const position_type& x) const
		{
			return ImplType::grad(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* gOut, const position_type& x) const
		{
			ImplType::grads(gOut, x);
		}
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


/// Exception thrown when local shape function set not found
struct UG_ERROR_LocalShapeFunctionSetNotRegistered
	: public UGFatalError
{
	UG_ERROR_LocalShapeFunctionSetNotRegistered(int dim_, ReferenceObjectID roid_, LFEID lfeid_)
	: UGFatalError(""), dim(dim_), roid(roid_), lfeid(lfeid_)
	{
		std::stringstream ss; ss << "Local Shape Function Set not found for "
							<<roid<<" (dim="<<dim<<") and type = "<<lfeid;
		UGFatalError::set_msg(ss.str());
	}
	int dim;
	ReferenceObjectID roid;
	LFEID lfeid;
};

// LocalShapeFunctionSetProvider
/** class to provide local shape function sets
 *
 *	This class provides references to Local Shape functions sets.
 *	It is implemented as a Singleton.
 */
class LocalShapeFunctionSetProvider {
	private:
	// 	disallow private constructor
		LocalShapeFunctionSetProvider();

	// disallow copy and assignment (intentionally left unimplemented)
		LocalShapeFunctionSetProvider(const LocalShapeFunctionSetProvider&);
		LocalShapeFunctionSetProvider& operator=(const LocalShapeFunctionSetProvider&);

	// 	private destructor
		~LocalShapeFunctionSetProvider(){};

	// 	Singleton provider
		static LocalShapeFunctionSetProvider& inst()
		{
			static LocalShapeFunctionSetProvider myInst;
			return myInst;
		};

	private:
	// 	initialize the standard trialspaces (called during construction)
		template <typename TRefElem>
		bool init_standard_sets();

	// 	clears all maps
		template <typename TRefElem>
		static void clear_maps();

	// 	return a map of element_trial_spaces
		template <typename TRefElem>
		static std::map<LFEID, const LocalShapeFunctionSet<TRefElem>* >&
		get_map();

	// 	return a map of element_trial_spaces
		template <int dim>
		static std::vector<std::map<LFEID, const DimLocalShapeFunctionSet<dim>* > >&
		get_dim_map();

	public:
	/// register a local shape function set for a given reference element type
	/**
	 * This function is used to register a Local Shape Function set for an element
	 * type and the corresponding local shape function set id.
	 *
	 * \param[in]		id 		Identifier for local shape function set
	 * \param[in]		set		Local Shape Function Set to register
	 * \return			bool	true iff registration successful
	 */
		template <typename TRefElem>
		static bool
		register_set(LFEID id, const LocalShapeFunctionSet<TRefElem>& set);

	/// unregister a local shape function set for a given reference element type
	/**
	 *  This function is used to unregister a Local Shape Function set for an element
	 * type and the corresponding local shape function set id from this Provider.
	 *
	 * \param[in]		id 		Identifier for local shape function set
	 * \return			bool	true iff removal successful
	 */
		template <typename TRefElem>
		static bool unregister_set(LFEID id);

	///	returns the Local Shape Function Set
	/**
	 * This function returns the Local Shape Function Set for a reference element
	 * type and an Identifier if a set has been registered for the identifier.
	 * Else an exception is thrown.
	 *
	 * \param[in]	id		Identifier for local shape function set
	 * \return 		set		A const reference to the shape function set
	 */
		// get the local shape function set for a given reference element and id
		template <typename TRefElem>
		static const LocalShapeFunctionSet<TRefElem>& get(LFEID id);

	///	returns the Local Shape Function Set
	/**
	 *  This function returns the Local Shape Function Set for a reference element
	 * type and an Identifier if a set has been registered for the identifier.
	 * Else an exception is thrown.
	 *
	 * \param[in]	id		Identifier for local shape function set
	 * \return 		set		A const reference to the shape function set
	 */
		template <int dim>
		static const DimLocalShapeFunctionSet<dim>& get(ReferenceObjectID roid,
		                                             LFEID id);
};

} // namespace ug

#include "local_shape_function_set_impl.h"

#endif /* __H__LIBDISCRETIZATION__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__ */
