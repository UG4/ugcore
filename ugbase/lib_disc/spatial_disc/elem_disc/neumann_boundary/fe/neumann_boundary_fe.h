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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FE__

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "../neumann_boundary_base.h"

namespace ug{

template<typename TDomain>
class NeumannBoundaryFE
	: public NeumannBoundaryBase<TDomain>
{
	private:
	///	Base class type
		typedef NeumannBoundaryBase<TDomain> base_type;

	///	Base class type
		typedef NeumannBoundaryFE<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	default constructor
		NeumannBoundaryFE(const char* function);

	///	add a boundary value
	///	\{
		void add(SmartPtr<CplUserData<number, dim> > data, 			const char* BndSubsets, const char* InnerSubsets);
		void add(SmartPtr<CplUserData<number, dim, bool> > user, 		const char* BndSubsets, const char* InnerSubsets);
		void add(SmartPtr<CplUserData<MathVector<dim>, dim> > user, 	const char* BndSubsets, const char* InnerSubsets);
	/// \}

	protected:
		using typename base_type::Data;

	///	Unconditional scalar user data
		struct NumberData : public base_type::Data
		{
			NumberData(SmartPtr<CplUserData<number, dim> > data,
					   std::string BndSubsets, std::string InnerSubsets,
					   NeumannBoundaryFE* this_)
				: base_type::Data(BndSubsets, InnerSubsets), This(this_)
			{
				import.set_data(data);
			}

			template<typename TElem, typename TFEGeom>
			void extract_bip(const TFEGeom& geo);

			template <typename TElem, typename TFEGeom>
			void lin_def(const LocalVector& u,
						 std::vector<std::vector<number> > vvvLinDef[],
						 const size_t nip);

			DataImport<number, dim> import;
			std::vector<MathVector<dim> > vLocIP;
			std::vector<MathVector<dim> > vGloIP;
			NeumannBoundaryFE* This;
		};
		friend struct NumberData;

	///	Conditional scalar user data
		struct BNDNumberData : public base_type::Data
		{
			BNDNumberData(SmartPtr<CplUserData<number, dim, bool> > functor_,
						  std::string BndSubsets, std::string InnerSubsets)
				: base_type::Data(BndSubsets, InnerSubsets), functor(functor_) {}

			SmartPtr<CplUserData<number, dim, bool> > functor;
		};

	///	Unconditional vector user data
		struct VectorData : public base_type::Data
		{
			VectorData(SmartPtr<CplUserData<MathVector<dim>, dim> > functor_,
					   std::string BndSubsets, std::string InnerSubsets)
			: base_type::Data(BndSubsets, InnerSubsets), functor(functor_) {}

			SmartPtr<CplUserData<MathVector<dim>, dim> > functor;
		};

		std::vector<NumberData> m_vNumberData;
		std::vector<BNDNumberData> m_vBNDNumberData;
		std::vector<VectorData> m_vVectorData;

		void update_subset_groups();

	public:
	///	type of trial space for each function used
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid);

	protected:
	///	current order of disc scheme
		int m_order;
		int m_quadOrder;

	///	current shape function set
		LFEID m_lfeID;

	///	current inner subset
		int m_si;

	protected:
	///	assembling functions for fv1
	///	\{
		template<typename TElem, typename TFEGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);
		template<typename TElem, typename TFEGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);
		template<typename TElem, typename TFEGeom>
		void finish_elem_loop();
		template<typename TElem, typename TFEGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	/// \}

		static const int _C_ = 0;

	protected:
		void register_all_funcs(int order);
		template<typename TElem, typename TFEGeom> void register_func();
};

} // end namespac ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FE1__*/
