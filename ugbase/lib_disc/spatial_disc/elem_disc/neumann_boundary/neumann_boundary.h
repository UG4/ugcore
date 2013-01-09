/*
 * neumann_boundary.h
 *
 *  Created on: 26.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__

#include <boost/function.hpp>

// other ug4 modules
#include "common/common.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"
#include "common/util/provider.h"

namespace ug{

template<typename TDomain>
class NeumannBoundary
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	Base class type
		typedef NeumannBoundary<TDomain> this_type;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	public:
	///	default constructor
		NeumannBoundary(const char* subsets);

	///	adds a lua callback (cond and non-cond)
#ifdef UG_FOR_LUA
		void add(const char* name, const char* function, const char* subsets);
#endif

	///	add a boundary value
	///	\{
		void add(number val, const char* function, const char* subsets);
		void add(SmartPtr<UserData<number, dim> > data, const char* function, const char* subsets);
		void add(SmartPtr<UserData<number, dim, bool> > user, const char* function, const char* subsets);
		void add(SmartPtr<UserData<MathVector<dim>, dim> > user, const char* function, const char* subsets);
		void add(const std::vector<number>& val, const char* function, const char* subsets);
	/// \}

	private:
	///	base class for user data
		struct Data
		{
			Data(std::string fctName_, std::string ssName_)
							: fctName(fctName_), ssNames(ssName_) {}
			size_t locFct;
			std::string fctName;
			SubsetGroup ssGrp;
			std::string ssNames;
		};

	///	Unconditional scalar user data
		struct NumberData : public Data
		{
			NumberData(SmartPtr<UserData<number, dim> > data,
			           std::string fctName_, std::string ssName_)
				: Data(fctName_, ssName_)
			{
				import.set_data(data);
			}

			template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
			void extract_bip(const TFVGeom<TElem,dim>& geo);

			template <typename TElem, template <class Elem, int  Dim> class TFVGeom>
			void lin_def_fv1(const LocalVector& u,
			                 std::vector<std::vector<number> > vvvLinDef[],
			                 const size_t nip);

			DataImport<number, dim> import;
			std::vector<MathVector<dim> > vLocIP;
			std::vector<MathVector<dim> > vGloIP;
		};

	///	Conditional scalar user data
		struct BNDNumberData : public Data
		{
			BNDNumberData(SmartPtr<UserData<number, dim, bool> > functor_,
						  std::string fctName_, std::string ssName_)
				: Data(fctName_, ssName_), functor(functor_) {}

			SmartPtr<UserData<number, dim, bool> > functor;
		};

	///	Unconditional vector user data
		struct VectorData : public Data
		{
			VectorData(SmartPtr<UserData<MathVector<dim>, dim> > functor_,
			           std::string fctName_, std::string ssName_)
			: Data(fctName_, ssName_), functor(functor_) {}

			SmartPtr<UserData<MathVector<dim>, dim> > functor;
		};

		std::vector<NumberData> m_vNumberData;
		std::vector<BNDNumberData> m_vBNDNumberData;
		std::vector<VectorData> m_vVectorData;

	///	method used to extract subsets and function id
		void extract_data();
		void extract_data(Data& userData, FunctionGroup& commonFctGrp, std::string& fctNames);

	///	callback invoked, when approximation space is changed
		virtual void approximation_space_changed() {extract_data();}

	public:
	///	 returns the type of elem disc
		virtual int type() const {return EDT_BND;}

	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID);

	///	switches between non-regular and regular grids
		virtual bool request_non_regular_grid(bool bNonRegular);

	private:
		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void prep_elem_loop_fv1();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void prep_elem_fv1(TElem* elem, const LocalVector& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void finish_elem_loop_fv1();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void add_JA_elem_fv1(LocalMatrix& J, const LocalVector& u) {}

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void add_JM_elem_fv1(LocalMatrix& J, const LocalVector& u) {}

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void add_dA_elem_fv1(LocalVector& d, const LocalVector& u) {}

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void add_dM_elem_fv1(LocalVector& d, const LocalVector& u) {}

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		void add_rhs_elem_fv1(LocalVector& d);

	private:
		void register_all_fv1_funcs(bool bHang);

		template <template <class Elem, int WorldDim> class TFVGeom>
		struct RegisterFV1 {
				RegisterFV1(this_type* pThis) : m_pThis(pThis){}
				this_type* m_pThis;
				template< typename TElem > void operator()(TElem&)
				{m_pThis->register_fv1_func<TElem, TFVGeom>();}
		};

		template <typename TElem, template <class Elem, int WorldDim> class TFVGeom>
		void register_fv1_func();

};

} // end namespac ug

#endif /*__H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__NEUMANN_BOUNDARY__FV1__NEUMANN_BOUNDARY__*/
