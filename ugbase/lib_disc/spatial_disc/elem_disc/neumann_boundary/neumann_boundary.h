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
#include "lib_grid/lg_base.h"

// library intern headers
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug{

template<typename TDomain>
class FV1NeumannBoundaryElemDisc
	: public IDomainElemDisc<TDomain>
{
	private:
	///	Base class type
		typedef IDomainElemDisc<TDomain> base_type;

	///	Base class type
		typedef FV1NeumannBoundaryElemDisc<TDomain> this_type;

	///	explicitly forward function
		using base_type::time;

	public:
	///	Domain type
		typedef typename base_type::domain_type domain_type;

	///	World dimension
		static const int dim = base_type::dim;

	///	Position type
		typedef typename base_type::position_type position_type;

	protected:
	///	type of bnd number
		typedef boost::function<bool (number& value, const MathVector<dim>& x, number time)> BNDNumberFunctor;
		typedef boost::function<void (MathVector<dim>& value, const MathVector<dim>& x, number time)> VectorFunctor;

	public:
	///	default constructor
		FV1NeumannBoundaryElemDisc(const char* subsets);

	///	add a boundary value
		void add(BNDNumberFunctor& user, const char* function, const char* subsets);
		void add(VectorFunctor& user, const char* function, const char* subsets);

	private:
	///	Functor, function grouping
		struct BNDNumberData
		{
			BNDNumberData(BNDNumberFunctor functor_,
			              std::string fctName_, std::string ssName_)
				: functor(functor_), fctName(fctName_), ssNames(ssName_) {}

			BNDNumberFunctor functor;
			size_t loc_fct;
			std::string fctName;
			SubsetGroup ssGrp;
			std::string ssNames;
		};

	///	Functor, function grouping
		struct VectorData
		{
			VectorData(VectorFunctor functor_,
			           std::string fctName_, std::string ssName_)
			: functor(functor_), fctName(fctName_), ssNames(ssName_) {}

			VectorFunctor functor;
			size_t loc_fct;
			std::string fctName;
			SubsetGroup ssGrp;
			std::string ssNames;
		};

		std::vector<BNDNumberData> m_vBNDNumberData;
		std::vector<VectorData> m_vVectorData;

		std::map<int, std::vector<BNDNumberData*> > m_mBNDNumberBndSegment;
		std::map<int, std::vector<VectorData*> > m_mVectorBndSegment;

		template <typename TUserData>
		void extract_data(std::map<int, std::vector<TUserData*> >& mvUserDataBndSegment,
		                  std::vector<TUserData>& vUserData,
		                  FunctionGroup& commonFctGrp, std::string& fctNames);

		void extract_data();

		virtual void approximation_space_changed() {extract_data();}

	public:
	///	type of trial space for each function used
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID)
		{
		//	check that Lagrange 1st order
			for(size_t i = 0; i < vLfeID.size(); ++i)
				if(vLfeID[i] != LFEID(LFEID::LAGRANGE, 1)) return false;
			return true;
		}

	///	switches between non-regular and regular grids
		virtual bool treat_non_regular_grid(bool bNonRegular)
		{
		//	switch, which assemble functions to use.
			if(bNonRegular)
			{
				UG_LOG("ERROR in 'FVNeumannBoundaryElemDisc::treat_non_regular_grid':"
						" Non-regular grid not implemented.\n");
				return false;
			}

		//	this disc supports regular grids
			return true;
		}


	private:
		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool prepare_element_loop();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool prepare_element(TElem* elem, const LocalVector& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool finish_element_loop();

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_JA(LocalMatrix& J, const LocalVector& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_JM(LocalMatrix& J, const LocalVector& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_A(LocalVector& d, const LocalVector& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_M(LocalVector& d, const LocalVector& u);

		template<typename TElem, template <class Elem, int  Dim> class TFVGeom>
		bool assemble_f(LocalVector& d);

	private:
	// 	position access
		const position_type* m_vCornerCoords;

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
