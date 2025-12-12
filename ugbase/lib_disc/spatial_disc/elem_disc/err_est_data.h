/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Dmitry Logashenko, Markus Breit
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

/*
 *	Data shared by element discretizations for a-posteriori error estimation
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__

// extern headers
#include <vector>
#include <string>
// #include <limits>

// intern headers
#include "lib_grid/tools/surface_view.h"
// #include "lib_grid/algorithms/multi_grid_util.h"
#include "lib_disc/function_spaces/integrate.h"

#ifdef UG_PARALLEL
 	//#include "lib_grid/parallelization/util/compol_attachment_reduce.h"
 	//#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#endif


namespace ug {

/// Base class for error estimator data
/**
 * This virtual class should be the base of any particular error estimator
 * implemented in the elem_disc's. Every elem_disc class (not object!) should
 * declare its own derived class for keeping intermediate information that
 * should be accumulated in the computation of the local error estimators.
 * Several objects of the elem_disc class may share the same object of the
 * derived class for a consistent computation of the error estimator.
 */
template <typename TDomain>
class IErrEstData
{
	public:
	/// world dimension
		static constexpr int dim = TDomain::dim;
	
	///	class constructor
		IErrEstData () : m_consider(true), m_scale(1.0) {};
		
	///	virtual class destructor
		virtual ~IErrEstData () {};

	///	virtual function to allocate data structures for the error estimator
		virtual void alloc_err_est_data (ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl) = 0;
		
	///	virtual function called after the computation of the error estimator data in all the elements
		virtual void summarize_err_est_data (SmartPtr<TDomain> spDomain) = 0;

	/// calculate L2 integrals
		virtual number get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[]) = 0;
		
	///	virtual function to release data structures for the error estimator
		virtual void release_err_est_data () = 0;

	/// virtual function granting get access to the m_consider member
		bool consider_me() const {return m_consider;};

	/// whether or not this instance is to be considered by domainDisc
	/**
	 * The domainDisc calls alloc_err_est_data(), summarize_err_est_data(),
	 * get_elem_error_indicator() etc. only for ErrEstData objects that have
	 * consider_me() == true.
	 * This is useful when using a MultipleErrEstData object combining ErrEstData objects
	 * which are already set to some ElemDisc. In this case, set_consider_me(false).
	 * The default value is true.
	 */
		void set_consider_me(bool b) {m_consider = b;};

	/// set scaling factor for final error calculation
		void set_scaling_factor(number scale) {m_scale = scale;};

	/// get scaling factor
	/**
	 * 	This factor is used in the element-wise calculation of error indicators.
	 * 	After calculation of indicators for each IErrEstData object in the domain
	 * 	discretization (typically one per equation), the final overall indicator
	 * 	is calculated as weighted sum (with the scaling factors as weights).
	 * 	The default value (if not set) is 1.0.
	 *
	 * @return scale scaling factor
	 */
		number scaling_factor() {return m_scale;}

	private:
	bool m_consider;
	number m_scale;
};

/// Error estimator data class storing one scalar number per side
/**
 * This class allocates an attachment keeping one number per full-dimensional
 * element side. Furthermore, the data are collected at the boundaries of the
 * patches (in the case of the adaptive refinement).
 *
 * \tparam TDomain	domain type
 */
template <typename TDomain>
class SideFluxErrEstData : public IErrEstData<TDomain>
{
public:
	///	domain type
		using domain_type = TDomain;
		
	/// world dimension
		static constexpr int dim = TDomain::dim;
		
	///	type of the sides (face, edge) and the elems (volume, face)
		using side_type = typename domain_traits<dim>::side_type;
	
public:
	/// constructor
		SideFluxErrEstData() : IErrEstData<TDomain>() {};

	///	virtual class destructor
		virtual ~SideFluxErrEstData() = default;

//	Functions to access data

	///	get the data reference for a given side
		number& operator ()
		(
			side_type* pSide ///< pointer to the side
		)
		{
			return m_aaFluxJump[pSide];
		};
		
	///	get the surface view
		ConstSmartPtr<SurfaceView>& surface_view() {return m_spSV;};

//	Interface virtual functions inherited from IErrEstData

	///	virtual function to allocate data structures for the error estimator
		virtual void alloc_err_est_data (ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl);
		
	///	virtual function called after the computation of the error estimator data in all the elements
		virtual void summarize_err_est_data (SmartPtr<TDomain> spDomain);

	/// calculate L2 integrals
		virtual number get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[]) {return 0;};
		
	///	virtual function to release data structures of the error estimator
		virtual void release_err_est_data ();
	
private:
	///	Flux jumps for the error estimator
		ANumber m_aFluxJumpOverSide;
		
	///	Attachment accessor
		MultiGrid::AttachmentAccessor<side_type, ANumber> m_aaFluxJump;
		
	///	Grid for the attachment
		ConstSmartPtr<SurfaceView> m_spSV;
		
	///	Finest grid level
		GridLevel m_errEstGL;
};



/// Error estimator data class storing a number vector per side and per element.
/**
 * This class represents an H1 error estimator.
 * It can integrate expressions on elements and their sides with arbitrary order.
 * A vector of values at defined integration points is attached to any element and
 * side to that end.
 *
 * RECOMMENDED (INTENDED) USAGE
 * The data will typically consist of the values of certain functions at integration
 * points (IP) on the sides and the element.
 * A pointer to an object of this class can be handed to any element discretization
 * involved in a discretization. They will access the attachments in their method
 * compute_err_est_elem and add their respective parts of the function to be
 * integrated for the error estimator. Exactly one of them (or maybe some other object)
 * then has to do the actual integration using the given values at the IPs and add up
 * side and element terms according to the error estimator formula used.
 *
 * \tparam TDomain	domain type
 */
template <typename TDomain>
class SideAndElemErrEstData : public IErrEstData<TDomain>
{
public:
	///	domain type
	using domain_type = TDomain;

	/// world dimension
		static constexpr int dim = TDomain::dim;

	///	type of the sides (face, edge) and the elems (volume, face)
	using side_type = typename domain_traits<dim>::side_type;
	using elem_type = typename domain_traits<dim>::element_type;

	/// attachment type
	using attachment_type = Attachment<std::vector<number> >;

	/// this class
	using this_type = SideAndElemErrEstData<TDomain>;

	/// maximal number of sides of any element
		static constexpr int MAX_NUM_SIDES = 8;

public:
	/// constructors
		SideAndElemErrEstData(size_t _sideOrder, size_t _elemOrder, const char* subsets);
		SideAndElemErrEstData(size_t _sideOrder, size_t _elemOrder,
							  std::vector<std::string> subsets = std::vector<std::string>(0));

	///	virtual class destructor
		~SideAndElemErrEstData() override = default;

	//	Functions to access data

	/// getting the side integration order
		size_t side_order() const {return sideOrder;}

	/// getting the elem integration order
		size_t elem_order() const {return elemOrder;}

	///	get the data reference for a given side and ip
		number& operator ()
		(
			side_type* pSide, 	///< pointer to the side
			size_t ip		///< integration point id on the side
		);

	///	get the data reference for a given elem and ip
		number& operator ()
		(
			elem_type* pElem, 	///< pointer to the elem
			size_t ip		///< integration point id on the elem
		);

	/// get the local side integration points for a specific roid
		template <int refDim>
		const MathVector<refDim>* side_local_ips(const ReferenceObjectID roid);

	/// get the local elem integration points for a specific roid
		template <int refDim>
		const MathVector<refDim>* elem_local_ips(const ReferenceObjectID roid);

	/// get all global side integration points
		MathVector<TDomain::dim>* all_side_global_ips(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// get the global side integration points for a specific side roid
		MathVector<TDomain::dim>* side_global_ips(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// get the global elem integration points for a specific roid
		MathVector<TDomain::dim>* elem_global_ips(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// get number of side IPs of a specific side
		size_t num_side_ips(const side_type* pSide);

	/// get number of side IPs of a specific side type
		size_t num_side_ips(const ReferenceObjectID roid);

	/// get number of first IP belonging to a specific side
		size_t first_side_ips(const ReferenceObjectID roid, const size_t side);

	/// get number of side IPs
		size_t num_all_side_ips(const ReferenceObjectID roid);

	/// get number of elem IPs
		size_t num_elem_ips(const ReferenceObjectID roid);

	/// get index of specific side IP in sideIP array returned by side_local_ips
		size_t side_ip_index(const ReferenceObjectID roid, const size_t side, const size_t ip);

	///	get the surface view
		ConstSmartPtr<SurfaceView>& surface_view () {return m_spSV;};

	//	virtual functions inherited from IErrEstData
	///	virtual function to allocate data structures for the error estimator
		void alloc_err_est_data (ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl) override;

	///	virtual function called after the computation of the error estimator data in all the elements
		void summarize_err_est_data (SmartPtr<TDomain> spDomain) override;

	/// calculate L2 integrals
		number get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[]) override;

	///	virtual function to release data structures of the error estimator
		void release_err_est_data () override;

	/// select L2/H1 Estimator
		void set_type(int type) {m_type = (type<=H1_ERROR_TYPE) ? type : H1_ERROR_TYPE;}

protected:
	/// initialization of quadrature (to be called during construction)
		void init_quadrature();

	/// helper struct for getting quadrature rules by use of mpl::lists
		template<int refDim>
		struct GetQuadRules
		{
				GetQuadRules(QuadratureRule<refDim>** ppQuadRule, size_t quadOrder) :
					m_ppQuadRule(ppQuadRule), m_quadOrder(quadOrder) {}
				QuadratureRule<refDim>** m_ppQuadRule;
				size_t m_quadOrder;
				template< typename TElem > void operator () (TElem&)
				{
					const ReferenceObjectID roid = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
					m_ppQuadRule[roid] =
						const_cast<QuadratureRule<refDim>*>(&QuadratureRuleProvider<refDim>::get(roid, m_quadOrder));
				}
		};

private:
	/// order of side and elem function approximations for integrating
		size_t sideOrder;
		size_t elemOrder;

	/// the subsets this error estimator will produce values for
		std::vector<std::string> m_vSs;
		SubsetGroup m_ssg;

	/// storage for integration rules
		QuadratureRule<dim-1>* quadRuleSide[NUM_REFERENCE_OBJECTS];
		QuadratureRule<dim>* quadRuleElem[NUM_REFERENCE_OBJECTS];

	/// extra storage for local side IPs (elem IPs are contained in elem quad rules)
		std::vector<MathVector<TDomain::dim> > m_SideIPcoords[NUM_REFERENCE_OBJECTS];

	/// storage for global elem and side IPs
		std::vector<MathVector<TDomain::dim> > m_sideGlobalIPcoords;
		std::vector<MathVector<TDomain::dim> > m_singleSideGlobalIPcoords;	// not the most elegant solution...
		std::vector<MathVector<TDomain::dim> > m_elemGlobalIPcoords;

	/// the first index for IPs of a specific side in the sideIP series for a roid
		size_t m_sideIPsStartIndex[NUM_REFERENCE_OBJECTS][MAX_NUM_SIDES];

	///	vector of attachments for sides
		attachment_type m_aSide;

	///	vector of attachments for elems
		attachment_type m_aElem;

	///	vector of side attachment accessors
		MultiGrid::AttachmentAccessor<side_type, attachment_type > m_aaSide;

	///	vector of elem attachment accessors
		MultiGrid::AttachmentAccessor<elem_type, attachment_type > m_aaElem;

	///	Grid for the attachment
		ConstSmartPtr<SurfaceView> m_spSV;

	///	Finest grid level
		GridLevel m_errEstGL;

	/// set type
		enum type {L2_ERROR_TYPE=0, H1_ERROR_TYPE};
		int m_type;
};


/// Error estimator data class for discretizations with more than one unknown.
/**
 * This class is a kind of wrapper for a bundle of error estimator objects.
 * It can be useful if a discretization depends on more than one unknown and needs to
 * compute error estimators for both of them:
 * One can only pass one error estimator to any ElemDisc, but at the same time
 * it is necessary to calculate the errors for different unknowns separately!
 *
 * This class will not actually do anything, but pass any request on to the underlying
 * objects.
 * As they most probably figure in some other ElemDisc of their corresponding unknown,
 * the considerMe property is set to false by default.
 *
 * Make sure that the error estimation routines of any elem disc that is given this
 * MultipleErrEstData object write to the correct sub-objects (i.e. ErrEstData objects)!
 * There has to be some kind of mapping between the unknowns of the elem disc and
 * the order in which the single-function ErrEstData objects are added to this object.
 * This will typically be the exact same order. So try not to add the same object of
 * MultipleErrEstData to elem discs with unknowns defined in a different order!
 *
 * The template parameter TErrEstData must implement the IErrEstData interface.
 *
 * \todo Maybe find a better way to deal with the orders of ErrEstData objects here
 *       and unknowns in the elem discs.
 * \tparam TDomain	domain type
 */
template <typename TDomain, typename TErrEstData>
class MultipleErrEstData : public IErrEstData<TDomain>
{
	public:
	/// world dimension
		static constexpr int dim = TDomain::dim;

	///	class constructor
		MultipleErrEstData(ConstSmartPtr<ApproximationSpace<TDomain> > approx)
		: IErrEstData<TDomain>(), m_spApprox(approx),
		  m_fctGrp(m_spApprox->dof_distribution_info())
		{
			this->set_consider_me(false);
		}

	///	virtual class destructor
		~MultipleErrEstData() override = default;

	/// adding error estimator data objects
		virtual void add(SmartPtr<TErrEstData> spEed, const char* fct)
		{
			// check that fct is not already contained in fctGrp
			size_t uid = m_spApprox->fct_id_by_name(fct);
			if (m_fctGrp.contains(uid))
			{
				UG_THROW("Error estimator for function '" << fct << "' can not be added\n"
						 "as another error estimator object for the same function is already\n"
						 "contained here.")
			}

			// add to function group
			try
			{
				m_fctGrp.add(fct);
			}
			UG_CATCH_THROW("Error estimator data object for function '"
							<< fct << "' could not be added.");

			m_vEed.push_back(spEed.get());
		}

	/// getting the number of underlying error estimator data objects
		size_t num() const {return m_vEed.size();};

	/// accessing the underlying error estimator data objects via function id
		TErrEstData* get(size_t uid)
		{
			if (!m_fctGrp.contains(uid))
			{
				std::string name = m_spApprox->name(uid);
				UG_THROW("Trying to access error estimator data object "
						 "for unique function index " << uid << "(aka '" << name << "')\n"
						 "which is not present in this collection.")
			}

			return m_vEed[m_fctGrp.local_index(uid)];
		}

	//	inherited from IErrEstData
	///	virtual function to allocate data structures for the error estimator
		void alloc_err_est_data(ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl) override;

	///	virtual function called after the computation of the error estimator data in all the elements
		void summarize_err_est_data(SmartPtr<TDomain> spDomain) override;

	/// calculate L2 integrals
		number get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[]) override;

	///	virtual function to release data structures for the error estimator
		void release_err_est_data() override;

	protected:
		std::vector<TErrEstData*> m_vEed;

		/// approx space
		ConstSmartPtr<ApproximationSpace<TDomain> > m_spApprox;

		/// function group (in order to map fcts to error estimator objects)
		FunctionGroup m_fctGrp;
};


template <typename TDomain>
class MultipleSideAndElemErrEstData
	: public MultipleErrEstData<TDomain, SideAndElemErrEstData<TDomain> >
{
	public:
	/// world dimension
		static constexpr int dim = TDomain::dim;

	///	type of the sides (face, edge) and the elems (volume, face)
	using side_type = typename SideAndElemErrEstData<TDomain>::side_type;
	using elem_type = typename SideAndElemErrEstData<TDomain>::elem_type;

	/// constructor
		MultipleSideAndElemErrEstData(ConstSmartPtr<ApproximationSpace<TDomain> > approx)
			: MultipleErrEstData<TDomain, SideAndElemErrEstData<TDomain> >(approx),
			  m_bEqSideOrder(false), m_bEqElemOrder(false) {};

	/// destructor
		~MultipleSideAndElemErrEstData() override = default;

	/// adding error estimator data objects
	/// overrides parent add method; performs check for equal order after adding
		void add(SmartPtr<SideAndElemErrEstData<TDomain> > spEed, const char* fct) override;

	/// returns whether all underlying err ests have the same elem and side integration orders
		bool equal_side_order() const {return m_bEqSideOrder;}

	/// returns whether all underlying err ests have the same elem and side integration orders
		bool equal_elem_order() const {return m_bEqElemOrder;}

	protected:
		/// find out whether all underlying err_ests have the same integration orders
		/// (makes assembling easier)
			void check_equal_order();

		/// find out whether all underlying err_ests have the same side integration orders
		/// (makes assembling easier)
			void check_equal_side_order();

		/// find out whether all underlying err_ests have the same elem integration orders
		/// (makes assembling easier)
			void check_equal_elem_order();

	private:
		bool m_bEqSideOrder;
		bool m_bEqElemOrder;

};

} // end of namespace ug

#include "err_est_data_impl.h"

#endif
