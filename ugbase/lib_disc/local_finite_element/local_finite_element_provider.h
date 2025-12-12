/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_FINITE_ELEMENT_PROVIDER__
#define __H__UG__LIB_DISC__LOCAL_FINITE_ELEMENT__LOCAL_FINITE_ELEMENT_PROVIDER__

// extern libraries
//#include <cassert>
#include <map>

// other ug4 modules
#include "common/math/ugmath.h"

// library intern headers
#include "local_finite_element_id.h"
#include "lib_disc/local_finite_element/local_shape_function_set.h"
#include "lib_disc/local_finite_element/local_dof_set.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
//	Provider for local finite element infos
////////////////////////////////////////////////////////////////////////////////

// LocalFiniteElementProvider
/** class to provide infos on local finite elements
 *
 *	This class provides references to Local Shape functions and Local DoF Sets.
 *	It is implemented as a Singleton.
 */
class LocalFiniteElementProvider {
	private:
	// 	disallow private constructor
		LocalFiniteElementProvider();

	// disallow copy and assignment (intentionally left unimplemented)
		LocalFiniteElementProvider(const LocalFiniteElementProvider&);
		LocalFiniteElementProvider& operator = (const LocalFiniteElementProvider&);

	// 	private destructor
		~LocalFiniteElementProvider();

	// 	Singleton provider
		inline static LocalFiniteElementProvider& inst()
		{
			static LocalFiniteElementProvider myInst;
			return myInst;
		};

	private:
	/// create the standard lagrange space
	///	\{
		template <typename TRefElem>
		static void create_lagrange_set(const LFEID& id);
		static void create_lagrange_set(ReferenceObjectID roid, const LFEID& id);
	///	\}

	/// create the mini bubble space
	///	\{
		template <typename TRefElem>
		static void create_mini_bubble_set(const LFEID& id);
		static void create_mini_bubble_set(ReferenceObjectID roid, const LFEID& id);
	///	\}

	/// creates nedelec space
	///	\{
		template <typename TRefElem>
		static void create_nedelec_set(const LFEID& id);
		static void create_nedelec_set(ReferenceObjectID roid, const LFEID& id);
	///	\}

	/// creates piecewise constant space
	///	\{
		template <typename TRefElem>
		static void create_piecewise_constant_set(const LFEID& id);
		static void create_piecewise_constant_set(ReferenceObjectID roid, const LFEID& id);
	///	\}

	/// creates crouxeiz-raviart space
	///	\{
		template <typename TRefElem>
		static void create_crouxeiz_raviart_set(const LFEID& id);
		static void create_crouxeiz_raviart_set(ReferenceObjectID roid, const LFEID& id);
	///	\}

	///	creates new set at runtime if available
		static void create_set(ReferenceObjectID roid, const LFEID& id);

	///	creates new set at runtime if available
		static void create_set(const LFEID& id);

	///	creates dof set on sub elements
		template <int rdim, int dim>
		static void create_sub_dof_set(ReferenceObjectID roid, const LFEID& id);
		static void create_dof_set(ReferenceObjectID roid, const LFEID& id);

	private:
		template <int dim, typename TShape, typename TGrad>
		struct LocalShapeFunctionSets{
			ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> >& operator [] (size_t i)  {return ptr[i];}
			const ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> >& operator [] (size_t i) const {return ptr[i];}
			ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > ptr[NUM_REFERENCE_OBJECTS];
		};

	// 	return a map of element_trial_spaces
		template <int dim, typename TShape, typename TGrad>
		static std::map<LFEID, LocalShapeFunctionSets<dim, TShape, TGrad> >&
		lsfs_map();

	//	returns the continuous information
		static std::map<LFEID, bool> m_mContSpace;

	private:
		struct LocalDoFSets{
			ConstSmartPtr<LocalDoFSet>& operator [] (size_t i)  {return ptr[i];}
			const ConstSmartPtr<LocalDoFSet>& operator [] (size_t i) const {return ptr[i];}
			ConstSmartPtr<LocalDoFSet> ptr[NUM_REFERENCE_OBJECTS];
		};

		template <int dim>
		struct DimLocalDoFSets{
			ConstSmartPtr<DimLocalDoFSet<dim> >& operator [] (size_t i)  {return ptr[i];}
			const ConstSmartPtr<DimLocalDoFSet<dim> >& operator [] (size_t i) const {return ptr[i];}
			ConstSmartPtr<DimLocalDoFSet<dim> > ptr[NUM_REFERENCE_OBJECTS];
		};


	//	map holding common dof set for roid of same dimension
		static std::map<LFEID, CommonLocalDoFSet> m_mCommonDoFSet;

	//	map holding dof sets for a reference object id
		static std::map<LFEID, LocalDoFSets> m_mLocalDoFSets;

	//	returns map for dof set with local pos
		template <int dim>
		static std::map<LFEID, DimLocalDoFSets<dim> >& lds_map();

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
		template <int dim, typename TShape, typename TGrad>
		static void	register_set(const LFEID& id,
		           	             ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> > set);

	/** register a local DoF set for a given reference element type
	 * This function is used to register a Local Shape Function set for an element
	 * type and the corresponding local DoF set id.
	 *
	 * \param[in]		id 		Identifier for local DoF set
	 * \param[in]		set		Local Shape Function Set to register
	 */
		template <int dim>
		static void register_set(const LFEID& id, ConstSmartPtr<DimLocalDoFSet<dim> > set);

	/** register a local DoF set for a given reference element type
	 * This function is used to register a Local Shape Function set for an element
	 * type and the corresponding local DoF set id.
	 *
	 * \param[in]		id 		Identifier for local DoF set
	 * \param[in]		set		Local Shape Function Set to register
	 */
		static void register_set(const LFEID& id, ConstSmartPtr<LocalDoFSet> set);

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
	///\{
		template <int dim, typename TShape, typename TGrad>
		static const LocalShapeFunctionSet<dim, TShape, TGrad>&
		get(ReferenceObjectID roid, const LFEID& id, bool bCreate = true);

		template <int dim>
		static const LocalShapeFunctionSet<dim>&
		get(ReferenceObjectID roid, const LFEID& id, bool bCreate = true)
			{return get<dim,number,MathVector<dim> >(roid, id, bCreate);}

		template <int dim, typename TShape, typename TGrad>
		static ConstSmartPtr<LocalShapeFunctionSet<dim, TShape, TGrad> >
		getptr(ReferenceObjectID roid, const LFEID& id, bool bCreate = true);

		template <int dim>
		static ConstSmartPtr<LocalShapeFunctionSet<dim> >
		getptr(ReferenceObjectID roid, const LFEID& id, bool bCreate = true)
			{return getptr<dim,number,MathVector<dim> >(roid, id, bCreate);}
	///\}


	///	returns the common dof set for all reference objects of a dimension
		static const CommonLocalDoFSet& get_dofs(const LFEID& id, bool bCreate = true);

	///	returns the local DoF set base for an id
		static const LocalDoFSet& get_dofs(ReferenceObjectID roid, const LFEID& id, bool bCreate = true);

	///	returns the local DoF set base for an id
	/// \{
		template <int dim>
		static const DimLocalDoFSet<dim>&
		get_dofs(ReferenceObjectID roid, const LFEID& id, bool bCreate = true);

		template <int dim>
		static ConstSmartPtr<DimLocalDoFSet<dim> >
		get_dof_ptr(ReferenceObjectID roid, const LFEID& id, bool bCreate = true);
	/// \}

	///returns if a Local Shape Function Set is continuous
		static bool continuous(const LFEID& id, bool bCreate = true);
};

} // namespace ug

#include "local_finite_element_provider_impl.h"

#endif