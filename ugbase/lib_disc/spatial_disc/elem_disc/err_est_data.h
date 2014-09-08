/*
 * err_est_data.h
 *
 *	Data shared by element discretizations for a-posteriori error estimation
 *
 *  Created on: 25.02.2014
 *     Authors: Dmitriy Logashenko, Markus Breit
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__

// extern headers
#include <vector>
#include <string>
#include <limits>

// intern headers
#include "lib_grid/tools/surface_view.h"
#include "lib_grid/algorithms/multi_grid_util.h"
#include "lib_disc/function_spaces/integrate.h"
//#include "lib_disc/common/subset_util.h"	// for DimensionOfSubset

#ifdef UG_PARALLEL
 	#include "lib_grid/parallelization/util/compol_attachment_reduce.h"
 	#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#endif

#include <boost/mpl/for_each.hpp>


namespace ug{

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
		static const int dim = TDomain::dim;
	
	///	class constructor
		IErrEstData () : m_consider(true) {};
		
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

	/// virtual function granting set access to the m_consider member
		void set_consider_me(bool b) {m_consider = b;};

	private:
	///	whether or not the instance of this class is to be considered when the domainDisc
	///	calls alloc, summarize, get_elem etc. (true by default)
	bool m_consider;
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
		typedef TDomain domain_type;
		
	/// world dimension
		static const int dim = TDomain::dim;
		
	///	type of the sides (face, edge) and the elems (volume, face)
		typedef typename domain_traits<dim>::side_type side_type;
	
public:
	/// constructor
		SideFluxErrEstData() : IErrEstData<TDomain>() {};

	///	virtual class destructor
		virtual ~SideFluxErrEstData() {};

//	Functions to access data

	///	get the data reference for a given side
		number& operator()
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
 * This class represents an L2 error estimator.
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
		typedef TDomain domain_type;

	/// world dimension
		static const int dim = TDomain::dim;

	///	type of the sides (face, edge) and the elems (volume, face)
		typedef typename domain_traits<dim>::side_type side_type;
		typedef typename domain_traits<dim>::element_type elem_type;

	/// attachment type
		typedef Attachment<std::vector<number> > attachment_type;

	/// this class
		typedef SideAndElemErrEstData<TDomain> this_type;

	/// maximal number of sides of any element
		static const int MAX_NUM_SIDES = 6;

public:
	/// constructors
		SideAndElemErrEstData(std::size_t _sideOrder, std::size_t _elemOrder, const char* subsets);
		SideAndElemErrEstData(std::size_t _sideOrder, std::size_t _elemOrder,
							  std::vector<std::string> subsets = std::vector<std::string>(0));

	///	virtual class destructor
		virtual ~SideAndElemErrEstData() {};

	//	Functions to access data

	/// getting the side integration order
		std::size_t side_order() const {return sideOrder;}

	/// getting the elem integration order
		std::size_t elem_order() const {return elemOrder;}

	///	get the data reference for a given side and ip
		number& operator()
		(
			side_type* pSide, 	///< pointer to the side
			std::size_t ip		///< integration point id on the side
		);

	///	get the data reference for a given elem and ip
		number& operator()
		(
			elem_type* pElem, 	///< pointer to the elem
			std::size_t ip		///< integration point id on the elem
		);

	/// get the local side integration points for a specific roid
		template <int refDim>
		const MathVector<refDim>* side_local_ips(const ReferenceObjectID roid);

	/// get the local elem integration points for a specific roid
		template <int refDim>
		const MathVector<refDim>* elem_local_ips(const ReferenceObjectID roid);

	/// get all global side integration points
	//	globIPs MUST be of the size num_side_ips()!
		MathVector<dim>* all_side_global_ips(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// get the global side integration points for a specific side roid
	//	globIPs MUST be of the size num_side_ips()!
		MathVector<dim>* side_global_ips(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// get the global elem integration points for a specific roid
	//	globIPs MUST be of the size num_elem_ips()!
		MathVector<dim>* elem_global_ips(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// get number of side IPs of a specific side
		std::size_t num_side_ips(const side_type* pSide);

	/// get number of side IPs of a specific side type
		std::size_t num_side_ips(const ReferenceObjectID roid);

	/// get number of first IP belonging to a specific side
		std::size_t first_side_ips(const ReferenceObjectID roid, const std::size_t side);

	/// get number of side IPs
		std::size_t num_all_side_ips(const ReferenceObjectID roid);

	/// get number of elem IPs
		std::size_t num_elem_ips(const ReferenceObjectID roid);

	/// get index of specific side IP in sideIP array returned by side_local_ips
		std::size_t side_ip_index(const ReferenceObjectID roid, const std::size_t side, const std::size_t ip);

	///	get the surface view
		ConstSmartPtr<SurfaceView>& surface_view () {return m_spSV;};

	//	virtual functions inherited from IErrEstData
	///	virtual function to allocate data structures for the error estimator
		virtual void alloc_err_est_data (ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl);

	///	virtual function called after the computation of the error estimator data in all the elements
		virtual void summarize_err_est_data (SmartPtr<TDomain> spDomain);

	/// calculate L2 integrals
		virtual number get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	virtual function to release data structures of the error estimator
		virtual void release_err_est_data ();

protected:
	/// initialization of quadrature (to be called during construction)
		void init_quadrature();

	/// helper struct for getting quadrature rules by use of mpl::lists
		template<int refDim>
		struct GetQuadRules
		{
				GetQuadRules(QuadratureRule<refDim>** ppQuadRule, std::size_t quadOrder) :
					m_ppQuadRule(ppQuadRule), m_quadOrder(quadOrder) {}
				QuadratureRule<refDim>** m_ppQuadRule;
				std::size_t m_quadOrder;
				template< typename TElem > void operator()(TElem&)
				{
					const ReferenceObjectID roid = geometry_traits<TElem>::REFERENCE_OBJECT_ID;
					m_ppQuadRule[roid] =
						const_cast<QuadratureRule<refDim>*>(&QuadratureRuleProvider<refDim>::get(roid, m_quadOrder));
				}
		};

private:
	/// order of side and elem function approximations for integrating
		std::size_t sideOrder;
		std::size_t elemOrder;

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
		std::size_t m_sideIPsStartIndex[NUM_REFERENCE_OBJECTS][MAX_NUM_SIDES];

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
 * The template parameter TErrEstData must implement the IErrEstData interface.
 *
 * \tparam TDomain	domain type
 */
template <typename TDomain, typename TErrEstData>
class MultipleErrEstData : public IErrEstData<TDomain>
{
	public:
	/// world dimension
		static const int dim = TDomain::dim;

	///	class constructor
		MultipleErrEstData() : IErrEstData<TDomain>() {this->set_consider_me(false);};

	///	virtual class destructor
		virtual ~MultipleErrEstData() {};

	/// adding error estimator data objects
		virtual void add(SmartPtr<TErrEstData> spEed) {m_vEed.push_back(spEed.get());};

	/// getting the number of underlying error estimator data objects
		std::size_t num() const {return m_vEed.size();};

	/// accessing the underlying error estimator data objects
		TErrEstData* get(std::size_t eed)
		{
			if (eed >= num()) UG_THROW("Trying to access an index that does not exist.")
			return m_vEed[eed];
		}

	//	inherited from IErrEstData
	///	virtual function to allocate data structures for the error estimator
		virtual void alloc_err_est_data(ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl);

	///	virtual function called after the computation of the error estimator data in all the elements
		virtual void summarize_err_est_data(SmartPtr<TDomain> spDomain);

	/// calculate L2 integrals
		virtual number get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	virtual function to release data structures for the error estimator
		virtual void release_err_est_data();

	protected:
		std::vector<TErrEstData*> m_vEed;
};


template <typename TDomain>
class MultipleSideAndElemErrEstData
	: public MultipleErrEstData<TDomain, SideAndElemErrEstData<TDomain> >
{
	public:
	/// world dimension
		static const int dim = TDomain::dim;

	/// constructor
		MultipleSideAndElemErrEstData()
			: MultipleErrEstData<TDomain, SideAndElemErrEstData<TDomain> >(),
			  m_bEqSideOrder(false), m_bEqElemOrder(false) {};

	/// destructor
		virtual ~MultipleSideAndElemErrEstData() {};

	/// adding error estimator data objects
	/// overrides parent add method; performs check for equal order after adding
		virtual void add(SmartPtr<SideAndElemErrEstData<TDomain> > spEed);

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

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ERR_EST_DATA__ */

/* End of File */
