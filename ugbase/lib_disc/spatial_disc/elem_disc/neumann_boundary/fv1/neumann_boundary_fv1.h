/*
 * neumann_boundary_fv1.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "../neumann_boundary_base.h"

namespace ug{

template<typename TDomain>
class NeumannBoundaryFV1
	: public NeumannBoundaryBase<TDomain>
{
	private:
	///	Base class type
		typedef NeumannBoundaryBase<TDomain> base_type;

	///	Base class type
		typedef NeumannBoundaryFV1<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	default constructor
		NeumannBoundaryFV1(const char* function);

	///	add a boundary value
	///	\{
		void add(SmartPtr<UserData<number, dim> > data, 			const char* BndSubsets, const char* InnerSubsets);
		void add(SmartPtr<UserData<number, dim, bool> > user, 		const char* BndSubsets, const char* InnerSubsets);
		void add(SmartPtr<UserData<MathVector<dim>, dim> > user, 	const char* BndSubsets, const char* InnerSubsets);
	/// \}

	protected:
		using base_type::Data;

	///	Unconditional scalar user data
		struct NumberData : public base_type::Data
		{
			NumberData(SmartPtr<UserData<number, dim> > data,
			           std::string BndSubsets, std::string InnerSubsets)
				: base_type::Data(BndSubsets, InnerSubsets)
			{
				import.set_data(data);
			}

			template<typename TElem, typename TFVGeom>
			void extract_bip(const TFVGeom& geo);

			template <typename TElem, typename TFVGeom>
			void lin_def(const LocalVector& u,
			             std::vector<std::vector<number> > vvvLinDef[],
			             const size_t nip);

			DataImport<number, dim> import;
			std::vector<MathVector<dim> > vLocIP;
			std::vector<MathVector<dim> > vGloIP;
		};

	///	Conditional scalar user data
		struct BNDNumberData : public base_type::Data
		{
			BNDNumberData(SmartPtr<UserData<number, dim, bool> > functor_,
						  std::string BndSubsets, std::string InnerSubsets)
				: base_type::Data(BndSubsets, InnerSubsets), functor(functor_) {}

			SmartPtr<UserData<number, dim, bool> > functor;
		};

	///	Unconditional vector user data
		struct VectorData : public base_type::Data
		{
			VectorData(SmartPtr<UserData<MathVector<dim>, dim> > functor_,
			           std::string BndSubsets, std::string InnerSubsets)
			: base_type::Data(BndSubsets, InnerSubsets), functor(functor_) {}

			SmartPtr<UserData<MathVector<dim>, dim> > functor;
		};

		std::vector<NumberData> m_vNumberData;
		std::vector<BNDNumberData> m_vBNDNumberData;
		std::vector<VectorData> m_vVectorData;

		void update_subset_groups();
		using base_type::update_subset_groups;

	///	current inner subset
		int m_si;

	public:
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID);

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular);

	protected:
	///	assembling functions for fv1
	///	\{
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);
		template<typename TElem, typename TFVGeom>
		void prep_elem(TElem* elem, const LocalVector& u);
		template<typename TElem, typename TFVGeom>
		void finish_elem_loop();
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d);
	/// \}

		static const int _C_ = 0;

	private:
		void register_all_funcs(bool bHang);
		template <typename TElem, typename TFVGeom>	void register_func();
};

} // end namespac ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__*/
