/*
 * neumann_boundary_fv.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV__

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "../neumann_boundary_base.h"

namespace ug{

template<typename TDomain>
class NeumannBoundaryFV
	: public NeumannBoundaryBase<TDomain>
{
	private:
	///	Base class type
		typedef NeumannBoundaryBase<TDomain> base_type;

	///	Base class type
		typedef NeumannBoundaryFV<TDomain> this_type;

	public:
	///	World dimension
		static const int dim = base_type::dim;

	public:
	///	default constructor
		NeumannBoundaryFV(const char* function);

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
					   NeumannBoundaryFV* this_)
				: base_type::Data(BndSubsets, InnerSubsets), This(this_)
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
			NeumannBoundaryFV* This;
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

	///	current shape function set
		LFEID m_lfeID;

	///	current inner subset
		int m_si;

	protected:
	///	assembling functions for fv1
	///	\{
		template<typename TElem, typename TFVGeom>
		void prep_elem_loop(const ReferenceObjectID roid, const int si);
		template<typename TElem, typename TFVGeom>
		void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);
		template<typename TElem, typename TFVGeom>
		void finish_elem_loop();
		template<typename TElem, typename TFVGeom>
		void add_rhs_elem(LocalVector& d, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	/// \}

		static const int _C_ = 0;

	protected:
		void register_all_funcs(int order);
		template<typename TElem, typename TFVGeom> void register_func();
};

} // end namespac ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY___NEUMANN_BOUNDARY_FV1__*/
