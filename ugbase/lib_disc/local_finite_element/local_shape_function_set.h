/*
 * local_shape_function_set.h
 *
 *  Created on: 12.05.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__

// extern libraries
#include <cassert>
#include <map>

// other ug4 modules
#include "common/math/ugmath.h"
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/reference_element/reference_element.h"
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

/// static interface for trial spaces
template <typename TImpl, int TDim>
class BaseLocalShapeFunctionSet
{
	public:
	///	type of implementation
		typedef TImpl ImplType;

	///	dimension of reference element
		static const int dim = TDim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef number shape_type;

	///	Gradient type
		typedef MathVector<dim> grad_type;

	//////////////////////////////////////////
	//	methods implemented by derived class
	//////////////////////////////////////////

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline LFEID type() const {return getImpl().type();}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return getImpl().num_sh();}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, position_type& pos) const
		{
			return getImpl().position(i, pos);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline shape_type shape(size_t i, const position_type& x) const
		{
			return getImpl().shape(i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		inline void grad(grad_type& g, size_t i, const position_type& x) const
		{
			return getImpl().grad(g, i, x);
		}

	//////////////////////////////////////////
	//	methods generated generically
	//////////////////////////////////////////

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		inline void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				sOut[sh] = shape(sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		inline void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(gOut[sh], sh, x);
		}

	protected:
	///	access to implementation
		ImplType& getImpl() {return static_cast<ImplType&>(*this);}

	///	const access to implementation
		const ImplType& getImpl() const {return static_cast<const ImplType&>(*this);}

};


/// virtual base class for local shape function sets
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
	///	type of shape functions
		virtual LFEID type() const = 0;

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
	///	\{
		virtual void shapes(shape_type* sOut, const position_type& x) const = 0;

		inline void shapes(std::vector<shape_type>& vShapeOut, const position_type& x) const
		{
			vShapeOut.resize(num_sh()); shapes(&vShapeOut[0], x);
		}
	///	\}

	/// evaluates the gradient of the shape function
	/** This function returns the gradient of Shape Function i at
	 * an element-local evaluation point.
	 * \param[in]	i		number of DoF
	 * \param[in]	x		Position on reference element (evaluation point)
	 * \return gradient at point
	 */
		virtual void grad(grad_type& g, size_t i, const position_type& x) const = 0;

	/// returns all gradients evaluated at a point
	/**
	 * This function returns the gradients of all Shape Functions at
	 * an element-local evaluation point in an array.
	 * \param[out]	gOut	Vector of gradients
	 * \param[in]	x		Position on reference element (evaluation point)
	 */
	///	\{
		virtual void grads(grad_type* gOut, const position_type& x) const = 0;

		inline void grads(std::vector<grad_type>& vGradOut, const position_type& x) const
		{
			vGradOut.resize(num_sh()); grads(&vGradOut[0], x);
		}
	///	\}

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
	///	Domain position type
		typedef typename ImplType::position_type position_type;

	///	Shape type
		typedef typename ImplType::shape_type shape_type;

	///	Gradient type
		typedef typename ImplType::grad_type grad_type;

	public:
	///	constructor
		LocalShapeFunctionSetWrapper(){}

	///	\copydoc ug::LocalShapeFunctionSet::type()
		virtual LFEID type() const {return ImplType::type();}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		virtual size_t num_sh() const {return ImplType::num_sh();}

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

	///	\copydoc ug::LocalShapeFunctionSet::grad()
		virtual void grad(grad_type& g, size_t i, const position_type& x) const
		{
			ImplType::grad(g, i, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		virtual void shapes(shape_type* sOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < ImplType::num_sh(); ++sh)
				sOut[sh] = ImplType::shape(sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		virtual void grads(grad_type* gOut, const position_type& x) const
		{
			for(size_t sh = 0; sh < ImplType::num_sh(); ++sh)
				ImplType::grad(gOut[sh], sh, x);
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
		UGFatalError::push_msg(ss.str());
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
		~LocalShapeFunctionSetProvider();

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

	// 	initialize the standard trialspaces (called during construction)
		template <typename TRefElem>
		static bool init_flex_lagrange(size_t order);

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

	//	vector of dynamically created spaces
		template <int dim>
		static std::vector<DimLocalShapeFunctionSet<dim>*>&
		get_dynamic_allocated_vector();

	//	creates new set at runtime if available
		static void dynamically_create_set(ReferenceObjectID roid, LFEID id);

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
		static const LocalShapeFunctionSet<TRefElem>& get(LFEID id, bool bCreate = true);

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
		                                                LFEID id, bool bCreate = true);
};

} // namespace ug

#include "local_shape_function_set_impl.h"

#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__ */
