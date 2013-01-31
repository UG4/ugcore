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

// library intern headers
#include "lib_disc/reference_element/reference_element.h"
#include "local_finite_element_id.h"

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
class LocalShapeFunctionSet
{
	public:
	///	Dimension, where shape functions are defined
		static const int dim = TDim;

	///	Domain position type
		typedef MathVector<dim> position_type;

	///	Shape type
		typedef TShape shape_type;

	///	Gradient type
		typedef TGrad grad_type;

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
		virtual bool position(size_t i, MathVector<dim>& pos) const = 0;

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
		virtual ~LocalShapeFunctionSet() {};
};

/// @}

////////////////////////////////////////////////////////////////////////////////
//	Common base class for local shape function sets to ease implementation
////////////////////////////////////////////////////////////////////////////////

/// static interface for trial spaces
template <typename TImpl, int TDim,
		  typename TShape = number,
		  typename TGrad = MathVector<TDim> >
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
		typedef TShape shape_type;

	///	Gradient type
		typedef TGrad grad_type;

	//////////////////////////////////////////
	//	methods implemented by derived class
	//////////////////////////////////////////

	///	\copydoc ug::LocalShapeFunctionSet::type()
		inline LFEID type() const {return getImpl().type();}

	///	\copydoc ug::LocalShapeFunctionSet::num_sh()
		inline size_t num_sh() const {return getImpl().num_sh();}

	///	\copydoc ug::LocalShapeFunctionSet::position()
		inline bool position(size_t i, MathVector<dim>& pos) const
		{
			return getImpl().position(i, pos);
		}

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

	///	\copydoc ug::LocalShapeFunctionSet::shape()
		inline void shape(shape_type& sh, size_t i, const MathVector<dim>& x) const
		{
			sh = shape(i, x);
		}


	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		inline void shapes(shape_type* vShape, const MathVector<dim>& x) const
		{
			for(size_t sh = 0; sh < num_sh(); ++sh)
				vShape[sh] = shape(sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::shapes()
		inline void shapes(std::vector<shape_type>& vShape, const MathVector<dim>& x) const
		{
			vShape.resize(num_sh()); shapes(&vShape[0], x);
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
			for(size_t sh = 0; sh < num_sh(); ++sh)
				grad(vGrad[sh], sh, x);
		}

	///	\copydoc ug::LocalShapeFunctionSet::grads()
		inline void grads(std::vector<grad_type>& vGrad, const MathVector<dim>& x) const
		{
			vGrad.resize(num_sh()); grads(&vGrad[0], x);
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


////////////////////////////////////////////////////////////////////////////////
//	Provider for local shape function sets
////////////////////////////////////////////////////////////////////////////////

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
		void init_standard_sets();

	// 	initialize the standard trialspaces (called during construction)
		template <typename TRefElem>
		static void init_flex_lagrange(size_t order);

	// 	return a map of element_trial_spaces
		template <int dim>
		static std::map<LFEID, const LocalShapeFunctionSet<dim>* >*
		get_map();

	//	vector of dynamically created spaces
		template <int dim>
		static std::vector<LocalShapeFunctionSet<dim>*>&
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
	 * \param[in]		roid	Reference Object id
	 * \param[in]		set		Local Shape Function Set to register
	 */
		template <int dim>
		static void	register_set(LFEID id,
		           	             const ReferenceObjectID roid,
		           	             const LocalShapeFunctionSet<dim>& set);

	/// unregister a local shape function set for a given reference element type
	/**
	 *  This function is used to unregister a Local Shape Function set for an element
	 * type and the corresponding local shape function set id from this Provider.
	 *
	 * \param[in]		id 		Identifier for local shape function set
	 * \return			bool	true iff removal successful
	 */
		template <int dim>
		static bool unregister_set(LFEID id);

	///	returns the Local Shape Function Set
	/**
	 *  This function returns the Local Shape Function Set for a reference element
	 * type and an Identifier if a set has been registered for the identifier.
	 * Else an exception is thrown.
	 *
	 * \param[in]	roid	Reference object id
	 * \param[in]	id		Identifier for local shape function set
	 * \return 		set		A const reference to the shape function set
	 */
		template <int dim>
		static const LocalShapeFunctionSet<dim>& get(ReferenceObjectID roid,
		                                             LFEID id, bool bCreate = true);
};

} // namespace ug

#include "local_shape_function_set_impl.h"

#endif /* __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_SHAPE_FUCNTION_SET__ */
