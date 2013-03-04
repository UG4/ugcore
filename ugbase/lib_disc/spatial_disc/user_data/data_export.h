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

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Data Export
////////////////////////////////////////////////////////////////////////////////

/// Data export
/**
 * A DataExport is user data produced by an element discretization.
 */
template <typename TData, int dim>
class DataExport : 	public DependentUserData<TData, dim>
{
	public:
	///	default constructor
		DataExport();

	///	implement compute() method of IUserData
		virtual void compute(LocalVector* u, GeometricObject* elem, bool bDeriv = false);

	///	sets the geometric object type
		virtual void set_roid(ReferenceObjectID id);

	///	register evaluation of export function
		template <typename TClass, int refDim>
		void set_fct(ReferenceObjectID id, TClass* obj,
		             void (TClass::*func)(
		            		 const LocalVector& u,
		            		 const MathVector<dim> vGlobIP[],
		            		 const MathVector<refDim> vLocIP[],
		            		 const size_t nip,
		            		 TData vValue[],
		            		 bool bDeriv,
		            		 std::vector<std::vector<TData> > vvvDeriv[]));

	///	register evaluation of export function
		template <int refDim>
		void set_fct(ReferenceObjectID id,
		             void (*func)(
		            		 const LocalVector& u,
		            		 const MathVector<dim> vGlobIP[],
		            		 const MathVector<refDim> vLocIP[],
		            		 const size_t nip,
		            		 TData vValue[],
		            		 bool bDeriv,
		            		 std::vector<std::vector<TData> > vvvDeriv[]));

	///	clears all export functions
		void clear_fct();

	///	returns if data depends on solution
		virtual bool zero_derivative() const {return false;}

	///	clear dependent data
		void clear() {m_vDependData.clear();}

	///	add data dependency
		void add_needed_data(SmartPtr<IUserData> data);

	///	remove needed data
		void remove_needed_data(SmartPtr<IUserData> data);

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vDependData.size();}

	///	return needed data
		virtual SmartPtr<IUserData> needed_data(size_t i) {return m_vDependData.at(i);}

	///	returns if the dependent data is ready for evaluation
		virtual void check_setup() const;

	protected:
		template <int refDim>
		struct traits{
			typedef boost::function<void (const LocalVector& u,
			                              const MathVector<dim> vGlobIP[],
			                              const MathVector<refDim> vLocIP[],
			                              const size_t nip,
			                              TData vValue[],
			                              bool bDeriv,
			                              std::vector<std::vector<TData> > vvvDeriv[])>
			EvalFunc;
		};

		template <int refDim>
		typename traits<refDim>::EvalFunc& eval_fct(ReferenceObjectID id) {return eval_fct(id, Int2Type<refDim>());}
		template <int refDim>
		const typename traits<refDim>::EvalFunc& eval_fct(ReferenceObjectID id) const {return eval_fct(id, Int2Type<refDim>());}

		typename traits<1>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<1>) {return m_vCompFct1d[id];}
		typename traits<2>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<2>) {return m_vCompFct2d[id];}
		typename traits<3>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<3>) {return m_vCompFct3d[id];}
		const typename traits<1>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<1>) const {return m_vCompFct1d[id];}
		const typename traits<2>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<2>) const {return m_vCompFct2d[id];}
		const typename traits<3>::EvalFunc& eval_fct(ReferenceObjectID id, Int2Type<3>) const {return m_vCompFct3d[id];}
		typename traits<1>::EvalFunc m_vCompFct1d[NUM_REFERENCE_OBJECTS];
		typename traits<2>::EvalFunc m_vCompFct2d[NUM_REFERENCE_OBJECTS];
		typename traits<3>::EvalFunc m_vCompFct3d[NUM_REFERENCE_OBJECTS];

		bool eval_fct_set(ReferenceObjectID id) const;

		template <int refDim>
		void comp(const LocalVector& u, bool bDeriv);

	protected:
	/// current Geom Object
		ReferenceObjectID m_id;

	///	data the export depends on
		std::vector<SmartPtr<IUserData> > m_vDependData;
};

} // end namespace ug

#include "data_export_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DATA_EXPORT__ */
