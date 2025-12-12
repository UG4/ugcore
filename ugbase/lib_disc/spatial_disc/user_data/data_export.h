/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__

#include "data_export.h"

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"
// #include "common/util/string_util.h"

// #include "lib_disc/reference_element/reference_mapping_provider.h"
// #include "lib_disc/local_finite_element/local_finite_element_provider.h"

#include "std_user_data.h"

namespace ug {


////////////////////////////////////////////////////////////////////////////////
// ValueDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class ValueDataExport
	: public StdDependentUserData<ValueDataExport<dim>,number,dim>
{
	public:
		ValueDataExport(const char* functions){this->set_functions(functions);}

		template <int refDim>
		void eval_and_deriv(number vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<number> > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = nullptr) const;

		virtual void check_setup() const;

		virtual bool continuous() const;
};

////////////////////////////////////////////////////////////////////////////////
// GradientDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class GradientDataExport
	: public StdDependentUserData<GradientDataExport<dim>, MathVector<dim>,dim>
{
	public:
		GradientDataExport(const char* functions){this->set_functions(functions);}

		template <int refDim>
		void eval_and_deriv(MathVector<dim> vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = nullptr) const;

		virtual void check_setup() const;

		virtual bool continuous() const{return false;}
};

////////////////////////////////////////////////////////////////////////////////
// VectorDataExport
////////////////////////////////////////////////////////////////////////////////

template <int dim>
class VectorDataExport
	: public StdDependentUserData<VectorDataExport<dim>, MathVector<dim>,dim>
{
	public:
		VectorDataExport(const char* functions){this->set_functions(functions);}

		template <int refDim>
		void eval_and_deriv(MathVector<dim> vValue[],
							const MathVector<dim> vGlobIP[],
							number time, int si,
							GridObject* elem,
							const MathVector<dim> vCornerCoords[],
							const MathVector<refDim> vLocIP[],
							const size_t nip,
							LocalVector* u,
							bool bDeriv,
							int s,
							std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
							const MathMatrix<refDim, dim>* vJT = nullptr) const;

		virtual void check_setup() const;

		virtual bool continuous() const;
};

////////////////////////////////////////////////////////////////////////////////
// Data Export
////////////////////////////////////////////////////////////////////////////////

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim>
class DataExport :
	public StdDependentUserData<DataExport<TData,dim>, TData, dim>
{
	public:
	///	default constructor
		DataExport(const char* functions);

		template <int refDim>
		void eval_and_deriv(TData vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<TData> > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = nullptr) const;

	///	register evaluation of export function
		template <typename TClass, int refDim>
		void set_fct(ReferenceObjectID id, TClass* obj,
		             void (TClass::*func)(	TData vValue[],
		      								const MathVector<dim> vGlobIP[],
		    								number time, int si,
		    								const LocalVector& u,
		      								GridObject* elem,
		      								const MathVector<dim> vCornerCoords[],
		      								const MathVector<refDim> vLocIP[],
		      								const size_t nip,
		      								bool bDeriv,
		      								std::vector<std::vector<TData> > vvvDeriv[]));

	///	register evaluation of export function
		template <int refDim>
		void set_fct(ReferenceObjectID id,
		             void (*func)(	TData vValue[],
									const MathVector<dim> vGlobIP[],
									number time, int si,
									const LocalVector& u,
									GridObject* elem,
									const MathVector<dim> vCornerCoords[],
									const MathVector<refDim> vLocIP[],
									const size_t nip,
									bool bDeriv,
									std::vector<std::vector<TData> > vvvDeriv[]));

	///	clears all export functions
		void clear_fct();

	///	clear dependent data
		void clear() {m_vDependData.clear();}

	///	add data dependency
		void add_needed_data(SmartPtr<ICplUserData<dim> > data);

	///	remove needed data
		void remove_needed_data(SmartPtr<ICplUserData<dim> > data);

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual SmartPtr<ICplUserData<dim> > needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return true;}

	protected:
		/* The following classes are used to implement the functors to support
		 * free functions and member functions. We do not use boost::bind or
		 * loki here, since they usually do not support that many arguments. In
		 * addition, arguments are known and can simply be hardcoded.
		 */

		// base class
		template <int refDim>
		class FunctorBase{
			public:
				virtual void operator () (TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GridObject* elem,
				                const MathVector<dim> vCornerCoords[],
				                const MathVector<refDim> vLocIP[],
				                const size_t nip,
				                bool bDeriv,
				                std::vector<std::vector<TData> > vvvDeriv[]) const = 0;
				virtual ~FunctorBase() = default;
		};

		// free function functor
		template <int refDim>
		class FreeFunctionFunctor : public FunctorBase<refDim>{
			using FreeFunc = void(*)(TData vValue[],
			                         const MathVector<dim> vGlobIP[],
			                         number time, int si,
			                         const LocalVector& u,
			                         GridObject* elem,
			                         const MathVector<dim> vCornerCoords[],
			                         const MathVector<refDim> vLocIP[],
			                         const size_t nip,
			                         bool bDeriv,
			                         std::vector<std::vector<TData> > vvvDeriv[]);

			public:
				FreeFunctionFunctor(FreeFunc f) : m_f(f) {}

				void operator ()(TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GridObject* elem,
				                const MathVector<dim> vCornerCoords[],
				                const MathVector<refDim> vLocIP[],
				                const size_t nip,
				                bool bDeriv,
				                std::vector<std::vector<TData> > vvvDeriv[]) const
				{
					m_f(vValue, vGlobIP, time, si, u, elem, vCornerCoords, vLocIP, nip, bDeriv, vvvDeriv);
				}

			protected:
				FreeFunc m_f;
		};

		// member function functor
		template <typename TClass, int refDim>
		class MemberFunctionFunctor : public FunctorBase<refDim>{
			using MemFunc = void(TClass::*)(TData vValue[],
			                                const MathVector<dim> vGlobIP[],
			                                number time, int si,
			                                const LocalVector& u,
			                                GridObject* elem,
			                                const MathVector<dim> vCornerCoords[],
			                                const MathVector<refDim> vLocIP[],
			                                const size_t nip,
			                                bool bDeriv,
			                                std::vector<std::vector<TData> > vvvDeriv[]);

			public:
				MemberFunctionFunctor(TClass* obj, MemFunc f) : m_pObj(obj), m_mf(f) {}

				void operator () (TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GridObject* elem,
				                const MathVector<dim> vCornerCoords[],
				                const MathVector<refDim> vLocIP[],
				                const size_t nip,
				                bool bDeriv,
				                std::vector<std::vector<TData> > vvvDeriv[]) const override {
					((m_pObj)->*m_mf)(vValue, vGlobIP, time, si, u, elem, vCornerCoords, vLocIP, nip, bDeriv, vvvDeriv);
				}

			protected:
				TClass* m_pObj;
				MemFunc m_mf;
		};

		// generic functor for both types
		template <int refDim>
		class Functor{
			public:
				Functor() : m_spImpl(nullptr) {}

				template <typename FreeFunc>
				Functor(FreeFunc f) : m_spImpl(new FreeFunctionFunctor<refDim>(f)) {}

				template <typename TClass, typename MemFunc>
				Functor(TClass* obj, MemFunc f) : m_spImpl(new MemberFunctionFunctor<TClass, refDim>(obj, f)) {}

				void operator () (TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GridObject* elem,
				                const MathVector<dim> vCornerCoords[],
				                const MathVector<refDim> vLocIP[],
				                const size_t nip,
				                bool bDeriv,
				                std::vector<std::vector<TData> > vvvDeriv[]) const
				{
					(*m_spImpl)(vValue, vGlobIP, time, si, u, elem, vCornerCoords, vLocIP, nip, bDeriv, vvvDeriv);
				}

				bool valid() const {return m_spImpl.valid();}
				bool invalid() const {return m_spImpl.invalid();}
				void invalidate() {m_spImpl = SmartPtr<FunctorBase<refDim> >();}

			protected:
				ConstSmartPtr<FunctorBase<refDim> > m_spImpl;
		};

		template <int refDim>
		Functor<refDim>& eval_fct(ReferenceObjectID id) {return eval_fct(id, Int2Type<refDim>());}
		template <int refDim>
		const Functor<refDim>& eval_fct(ReferenceObjectID id) const {return eval_fct(id, Int2Type<refDim>());}

		Functor<1>& eval_fct(ReferenceObjectID id, Int2Type<1>) {return m_vCompFct1d[id];}
		Functor<2>& eval_fct(ReferenceObjectID id, Int2Type<2>) {return m_vCompFct2d[id];}
		Functor<3>& eval_fct(ReferenceObjectID id, Int2Type<3>) {return m_vCompFct3d[id];}
		const Functor<1>& eval_fct(ReferenceObjectID id, Int2Type<1>) const {return m_vCompFct1d[id];}
		const Functor<2>& eval_fct(ReferenceObjectID id, Int2Type<2>) const {return m_vCompFct2d[id];}
		const Functor<3>& eval_fct(ReferenceObjectID id, Int2Type<3>) const {return m_vCompFct3d[id];}
		Functor<1> m_vCompFct1d[NUM_REFERENCE_OBJECTS];
		Functor<2> m_vCompFct2d[NUM_REFERENCE_OBJECTS];
		Functor<3> m_vCompFct3d[NUM_REFERENCE_OBJECTS];

		bool eval_fct_set(ReferenceObjectID id) const;

	protected:
	///	data the export depends on
		std::vector<SmartPtr<ICplUserData<dim> > > m_vDependData;
};

} // end namespace ug

#include "data_export_impl.h"

#endif