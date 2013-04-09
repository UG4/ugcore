/*
 * data_export.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__

#include "lib_disc/common/function_group.h"
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/common/groups_util.h"
#include "common/util/string_util.h"

#include "user_data.h"
#include "std/std_user_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Data Export
////////////////////////////////////////////////////////////////////////////////

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim>
class DataExport :
	public StdUserData<	DataExport<TData,dim>, DependentUserData<TData, dim>, TData,dim>
{
	public:
	///	default constructor
		DataExport();

	///	implement compute() method of IUserData
		virtual void compute(LocalVector* u, GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false);

		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			UG_THROW("DataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		inline void evaluate (TData vValue[],
		                      const MathVector<dim> vGlobIP[],
		                      number time, int si, const size_t nip) const
		{
			UG_THROW("DataExport: Solution, element and local ips required "
					"for evaluation, but not passed. Cannot evaluate.");
		}

		template <int refDim>
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
			evaluate<refDim>(&value,&globIP,time,si,u,elem,
			                 vCornerCoords,&locIP,1,NULL);
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			const Functor<refDim>& func = eval_fct<refDim>(m_id);
			(func)(vValue,vGlobIP,time,si,u,elem,
					vCornerCoords,vLocIP,nip, false, NULL);
		}

	///	sets the geometric object type
		virtual void set_roid(ReferenceObjectID id);

	///	register evaluation of export function
		template <typename TClass, int refDim>
		void set_fct(ReferenceObjectID id, TClass* obj,
		             void (TClass::*func)(	TData vValue[],
		      								const MathVector<dim> vGlobIP[],
		    								number time, int si,
		    								const LocalVector& u,
		      								GeometricObject* elem,
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
									GeometricObject* elem,
									const MathVector<dim> vCornerCoords[],
									const MathVector<refDim> vLocIP[],
									const size_t nip,
									bool bDeriv,
									std::vector<std::vector<TData> > vvvDeriv[]));

	///	clears all export functions
		void clear_fct();

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return false;}

	///	clear dependent data
		void clear() {m_vDependData.clear();}

	///	add data dependency
		void add_needed_data(SmartPtr<IUserData<dim> > data);

	///	remove needed data
		void remove_needed_data(SmartPtr<IUserData<dim> > data);

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual SmartPtr<IUserData<dim> > needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if the dependent data is ready for evaluation
		virtual void check_setup() const;

	protected:
		/* The following classes are used to implement the functors to support
		 * free functions and member functions. We do not use boost::bind or
		 * loki here, since they usually do not support that many arguments. In
		 * addition arguments are known and can simply be hardcoded.
		 */

		// base class
		template <int refDim>
		class FunctorBase{
			public:
				virtual void operator()(TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GeometricObject* elem,
				                const MathVector<dim> vCornerCoords[],
				                const MathVector<refDim> vLocIP[],
				                const size_t nip,
				                bool bDeriv,
				                std::vector<std::vector<TData> > vvvDeriv[]) const = 0;
				virtual ~FunctorBase() {}
		};

		// free function functor
		template <int refDim>
		class FreeFunctionFunctor : public FunctorBase<refDim>{
			typedef void (*FreeFunc)(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								const LocalVector& u,
								GeometricObject* elem,
								const MathVector<dim> vCornerCoords[],
								const MathVector<refDim> vLocIP[],
								const size_t nip,
								bool bDeriv,
								std::vector<std::vector<TData> > vvvDeriv[]);

			public:
				FreeFunctionFunctor(FreeFunc f) : m_f(f) {}

				void operator()(TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GeometricObject* elem,
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
			typedef void (TClass::*MemFunc)(TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si,
								const LocalVector& u,
								GeometricObject* elem,
								const MathVector<dim> vCornerCoords[],
								const MathVector<refDim> vLocIP[],
								const size_t nip,
								bool bDeriv,
								std::vector<std::vector<TData> > vvvDeriv[]);

			public:
				MemberFunctionFunctor(TClass* obj, MemFunc f) : m_pObj(obj), m_mf(f) {}

				void operator()(TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GeometricObject* elem,
				                const MathVector<dim> vCornerCoords[],
				                const MathVector<refDim> vLocIP[],
				                const size_t nip,
				                bool bDeriv,
				                std::vector<std::vector<TData> > vvvDeriv[]) const
				{
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
				Functor() : m_spImpl(NULL) {}

				template <typename FreeFunc>
				Functor(FreeFunc f) : m_spImpl(new FreeFunctionFunctor<refDim>(f)) {}

				template <typename TClass, typename MemFunc>
				Functor(TClass* obj, MemFunc f) : m_spImpl(new MemberFunctionFunctor<TClass, refDim>(obj, f)) {}

				void operator()(TData vValue[],
				                const MathVector<dim> vGlobIP[],
								number time, int si,
				                const LocalVector& u,
				                GeometricObject* elem,
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

		template <int refDim>
		void comp(const LocalVector& u, GeometricObject* elem,
		          const MathVector<dim> vCornerCoords[], bool bDeriv);

	protected:
	/// current Geom Object
		ReferenceObjectID m_id;

	///	data the export depends on
		std::vector<SmartPtr<IUserData<dim> > > m_vDependData;
};

} // end namespace ug

#include "data_export_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__ */
