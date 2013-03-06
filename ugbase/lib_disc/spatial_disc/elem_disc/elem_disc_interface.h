/*
 * elem_disc_interface.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE__
#define __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE__

// extern headers
#include <vector>
#include <string>

// intern headers
#include "lib_disc/common/local_algebra.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "lib_disc/function_spaces/approximation_space.h"
#include "lib_disc/local_finite_element/local_finite_element_id.h"
#include "lib_disc/reference_element/reference_element_traits.h"
#include "lib_disc/spatial_disc/user_data/data_import.h"
#include "common/util/provider.h"
#include "lib_disc/domain_util.h"

namespace ug{

/// Types of elem disc
enum ElemDiscType
{
	EDT_NONE = 0,
	EDT_ELEM = 1 << 0,
	EDT_SIDE = 1 << 1,
	EDT_BND = 1 << 2,
	EDT_ALL = EDT_NONE | EDT_SIDE | EDT_ELEM | EDT_BND
};

/**
 * Element Discretizations
 *
 * \defgroup lib_disc_elem_disc Elem Disc
 * \ingroup lib_discretization
 */

/// \ingroup lib_disc_elem_disc
/// @{

///	base class for all element-wise discretizations
/**
 * This class is the base class for element-wise discretizations. An
 * implementation of this class must provide local stiffness/mass-matrix
 * contribution of one element to the global jacobian and local contributions
 * of one element to the local defect.
 */
class IElemDisc
{
	public:
	///	Constructor
		IElemDisc(const char* functions = NULL, const char* subsets = NULL);

	///	Constructor
		IElemDisc(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);

	////////////////////////////
	// Functions and Subsets

	///	sets functions by name list, divided by ','
		void set_functions(std::string functions);

	/// sets functions by vector of names
		void set_functions(const std::vector<std::string>& functions) {m_vFct = functions;};

	///	sets subset(s) by name list, divided by ','
		void set_subsets(std::string subsets);

	///	sets subset(s) by name list, divided by ','
		void set_subsets(const std::vector<std::string>& subsets) {m_vSubset = subsets;}

	/// number of functions this discretization handles
		size_t num_fct() const {return m_vFct.size();}

	///	returns the symbolic functions
		const std::vector<std::string>& symb_fcts() const {return m_vFct;}

	/// number of subsets this discretization handles
		size_t num_subsets() const {return m_vSubset.size();}

	///	returns the symbolic subsets
		const std::vector<std::string>& symb_subsets() const {return m_vSubset;}

	protected:
	///	vector holding name of all symbolic functions
		std::vector<std::string> m_vFct;

	///	vector holding name of all symbolic subsets
		std::vector<std::string> m_vSubset;

	////////////////////////////
	// UserData and Coupling
	////////////////////////////
	public:
	///	registers a data import
		void register_import(IDataImport& Imp);

	///	registers a data export
		void register_export(SmartPtr<IUserData> Exp);

	///	returns number of imports
		size_t num_imports() const {return m_vIImport.size();}

	/// returns an import
		IDataImport& get_import(size_t i);

	///	removes all imports
		void clear_imports() {m_vIImport.clear();}

	///	returns number of exports
		size_t num_exports() const {return m_vIExport.size();}

	/// returns an export
		SmartPtr<IUserData> get_export(size_t i);

	///	removes all exports
		void clear_exports() {m_vIExport.clear();}

	protected:
	/// data imports
		std::vector<IDataImport*> m_vIImport;

	///	data exports
		std::vector<SmartPtr<IUserData> > m_vIExport;

	public:
	////////////////////////////
	// Assembling functions
	////////////////////////////
	///	 returns the type of elem disc
		virtual int type() const {return EDT_ELEM | EDT_SIDE;}

	/// requests assembling for a finite element id
	/**
	 * This function is called before the assembling starts. In the vector
	 * exactly this->num_fct() Local Finite Element IDs must be passed. The
	 * IElemDisc-Implementation checks if it can assemble the set of LFEID and
	 * registers the corresponding assembling functions. If this is not the
	 * case instead false is returned.
	 *
	 * \param[in]		vLfeID		vector of Local Finite Element IDs
	 * \returns			true		if assemble routines are present and selected
	 * 					false		if no assembling for those Spaces available
	 */
		virtual bool request_finite_element_id(const std::vector<LFEID>& vLfeID) = 0;

	///	informs the assembling, that hanging nodes must be taken into account
	/**
	 * This method is called before the assembling of elements, that may have
	 * hanging nodes, constrained edges, etc. Thus, if the assembling must take
	 * special care, it can prepare for such needs. Typically other assembling
	 * functions are used then and registered.
	 *
	 * \param[in]		bNonRegular 	true iff non-regular grid used
	 * \returns			bool			true  if successful
	 * 									false if cannot be handled by disc
	 */
		virtual bool request_non_regular_grid(bool bNonRegular) = 0;

	///	returns if discretization acts on hanging nodes if present
	/**
	 * This function returns if a discretization really needs the hanging nodes
	 * in case of non-regular grid. This may not be the case for e.g. finite
	 * element assemblings but is needed for finite volumes
	 */
		virtual bool use_hanging() const {return false;}

	///	sets if assembling should be time-dependent and the local time series
	/**
	 * This function specifies if the assembling is time-dependent. If NULL is
	 * passed, the assembling is assumed to be time-independent. If a local
	 * time series is passed, this series is used as previous solution.
	 *
	 * \param[in]	locTimeSeries	Time series of previous solutions
	 */
		void set_time_dependent(const LocalVectorTimeSeries& locTimeSeries,
		        				const std::vector<number>& vScaleMass,
		        				const std::vector<number>& vScaleStiff);

	///	sets that the assembling is time independent
		void set_time_independent();

	///	returns if assembling is time-dependent
		bool is_time_dependent() const {return (m_pLocalVectorTimeSeries != NULL) && !m_bStationaryForced;}

	///	sets that the assembling is always stationary (even in instationary case)
		void set_stationary(bool bStationaryForced = true) {m_bStationaryForced = bStationaryForced;}
		void set_stationary() {set_stationary(true);}

	///	returns if assembling is forced to be stationary
		bool is_stationary() const {return m_bStationaryForced;}

	///	returns if local time series needed by assembling
	/**
	 * This callback must be implemented by a derived Elem Disc in order to handle
	 * time-dependent data. As return the derived Elem Disc can specify, if
	 * it really needs data from previous time steps for the (spatial) disc. The
	 * default is false.
	 *
	 * \returns 	if elem disc needs time series local solutions
	 */
		virtual bool requests_local_time_series() {return false;}

	///	sets the current time point
		void set_time_point(const size_t timePoint) {m_timePoint = timePoint;}

	///	returns the currently considered time point of the time-disc scheme
		size_t time_point() const {return m_timePoint;}

	///	returns currently set timepoint
		number time() const {if(m_pLocalVectorTimeSeries)
								return m_pLocalVectorTimeSeries->time(m_timePoint);
							else return 0.0;}

	///	returns the local time solutions
	/**
	 * This function returns the local time solutions. This a type of vector,
	 * that holds the local unknowns for each time point of a time series.
	 * Note, that the first sol is always the current (iterated, unknown)
	 * time point, while all remaining sols are the already computed time steps
	 * used e.g. in a time stepping scheme.
	 *
	 * \returns vLocalTimeSol		vector of local time Solutions
	 */
		const LocalVectorTimeSeries* local_time_solutions() const
			{return m_pLocalVectorTimeSeries;}

	///	returns the weight factors of the time-disc scheme
	///	\{
		const std::vector<number>& mass_scales() const {return m_vScaleMass;}
		const std::vector<number>& stiff_scales() const {return m_vScaleStiff;}

		number mass_scale(const size_t timePoint) const {return m_vScaleMass[timePoint];}
		number stiff_scale(const size_t timePoint) const {return m_vScaleStiff[timePoint];}

		number mass_scale() const {return m_vScaleMass[m_timePoint];}
		number stiff_scale() const {return m_vScaleStiff[m_timePoint];}
	///	\}

	/// prepare the timestep
	/**
	 * This function prepares a timestep (iff timedependent). This function is
	 * called once for every element before the spatial assembling procedure
	 * begins.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		template <typename TElem>
		void fast_prep_timestep_elem(TElem* elem, const LocalVector& u);

	/// prepare the timestep
		virtual void prep_timestep_elem(GeometricObject* elem, const LocalVector& u) {}

	///	prepares the loop over all elements of one type
	/**
	 * This function should prepare the element loop for elements of one fixed
	 * type. This function is called before e.g. the loop over all geometric
	 * objects of a chosen type is performed.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_prep_elem_loop(const ReferenceObjectID roid, const int si)
		{UG_ASSERT(m_vPrepareElemLoopFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vPrepareElemLoopFct[m_id])(roid, si);}

	///	virtual prepares the loop over all elements of one type
		virtual void prep_elem_loop(const ReferenceObjectID roid, const int si) {}

	///	prepare one elements for assembling
	/**
	 * This function prepares one Geometric object, that will be assembled in
	 * the next step.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 *
	 * \param[in]		elem		The geometric object
	 * \param[in]		u			The current local solution
	 */
		template <typename TElem>
		void fast_prep_elem(TElem* elem, const LocalVector& u);

	///	virtual prepare one elements for assembling
		virtual void prep_elem(GeometricObject* elem, const LocalVector& u) {}

	///	postprocesses the loop over all elements of one type
	/**
	 * This function should post process the element loop for elements of one fixed
	 * type. This function is called after e.g. the loop over all geometric
	 * objects of a chosen type has been performed.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_fsh_elem_loop()
		{UG_ASSERT(m_vFinishElemLoopFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vFinishElemLoopFct[m_id])();}

	///	virtual postprocesses the loop over all elements of one type
		virtual void fsh_elem_loop() {}

	/// finish the timestep
	/**
	 * This function finishes the timestep (iff timedependent). This function is
	 * called in the PostProcess of a timestep.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		template <typename TElem>
		void fast_fsh_timestep_elem(TElem* elem, const number time, const LocalVector& u);

	/// virtual finish the timestep
		virtual void fsh_timestep_elem(GeometricObject* elem, const number time, const LocalVector& u) {}

	/// Assembling of Jacobian (Stiffness part)
	/**
	 * This function assembles the local (stiffness) jacobian for the current
	 * solution u.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_add_jac_A_elem(LocalMatrix& J, const LocalVector& u)
		{UG_ASSERT(m_vElemJAFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vElemJAFct[m_id])(J, u);}

	/// Assembling of Jacobian (Stiffness part)
		virtual void add_jac_A_elem(GeometricObject* elem, LocalMatrix& J, const LocalVector& u) {}

	/// Assembling of Jacobian (Mass part)
	/**
	 * This function assembles the local (mass) jacobian for the current
	 * solution u.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_add_jac_M_elem(LocalMatrix& J, const LocalVector& u)
		{UG_ASSERT(m_vElemJMFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vElemJMFct[m_id])(J, u);}

	/// Assembling of Jacobian (Mass part)
		virtual void add_jac_M_elem(GeometricObject* elem, LocalMatrix& J, const LocalVector& u) {}

	/// Assembling of Defect (Stiffness part)
	/**
	 * This function assembles the local (stiffness) defect for the current
	 * solution u.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_add_def_A_elem(LocalVector& d, const LocalVector& u)
		{UG_ASSERT(m_vElemdAFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vElemdAFct[m_id])(d, u);}

		// explicit reaction, reaction_rate and source
   	    void fast_add_def_A_elem_explicit(LocalVector& d, const LocalVector& u)
   	    {
   	    	if(this->m_vElemdAFct_explicit[m_id] != NULL)
   	    		(this->*m_vElemdAFct_explicit[m_id])(d, u);
   	    }


	/// virtual Assembling of Defect (Stiffness part)
		virtual void add_def_A_elem(GeometricObject* elem, LocalVector& d, const LocalVector& u) {}

        // explicit defect for reaction, reaction_rate and source
		virtual void add_def_A_elem_explicit(GeometricObject* elem, LocalVector& d, const LocalVector& u) {}

	/// Assembling of Defect (Mass part)
	/**
	 * This function assembles the local (mass) defect for the current
	 * solution u.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_add_def_M_elem(LocalVector& d, const LocalVector& u)
		{UG_ASSERT(m_vElemdMFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vElemdMFct[m_id])(d, u);}

	/// virtual Assembling of Defect (Mass part)
		virtual void add_def_M_elem(GeometricObject* elem, LocalVector& d, const LocalVector& u) {}

	/// Assembling of Right-Hand Side
	/**
	 * This function assembles the local rhs.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_roid' must have been called to set the elem type.
	 */
		void fast_add_rhs_elem(LocalVector& rhs)
		{UG_ASSERT(m_vElemRHSFct[m_id]!=NULL, "Fast-Assemble Method missing.");
			(this->*m_vElemRHSFct[m_id])(rhs);}

	/// virtual Assembling of Right-Hand Side
		virtual void add_rhs_elem(GeometricObject* elem, LocalVector& rhs) {}

	/// Virtual destructor
		virtual ~IElemDisc() {}

	protected:
	///	number of functions
		size_t m_numFct;

	///	time point
		size_t m_timePoint;

	///	list of local vectors for all solutions of the time series
		const LocalVectorTimeSeries* m_pLocalVectorTimeSeries;

	///	weight factors for time dependent assembling
	/// \{
		std::vector<number> m_vScaleMass;
		std::vector<number> m_vScaleStiff;
	/// \}

	///	flag if stationary assembling is to be used even in instationary assembling
		bool m_bStationaryForced;

	private:
	//	abbreviation for own type
		typedef IElemDisc T;

	// 	types of timestep function pointers
		typedef void (T::*PrepareTimestepElemFct)(const LocalVector& u);
		typedef void (T::*FinishTimestepElemFct)(const LocalVector& u);

	// 	types of loop function pointers
		typedef void (T::*PrepareElemLoopFct)(ReferenceObjectID roid, int si);
		typedef void (T::*PrepareElemFct)(GeometricObject* obj, const LocalVector& u);
		typedef void (T::*FinishElemLoopFct)();

	// 	types of Jacobian assemble functions
		typedef void (T::*ElemJAFct)(LocalMatrix& J, const LocalVector& u);
		typedef void (T::*ElemJMFct)(LocalMatrix& J, const LocalVector& u);

	// 	types of Defect assemble functions
		typedef void (T::*ElemdAFct)(LocalVector& d, const LocalVector& u);
		typedef void (T::*ElemdMFct)(LocalVector& d, const LocalVector& u);

	// 	types of right hand side assemble functions
		typedef void (T::*ElemRHSFct)(LocalVector& d);

	protected:
	// 	register the functions
		template <typename TAssFunc> void set_prep_timestep_elem_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_fsh_timestep_elem_fct(ReferenceObjectID id, TAssFunc func);

		template <typename TAssFunc> void set_prep_elem_loop_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_prep_elem_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_fsh_elem_loop_fct(ReferenceObjectID id, TAssFunc func);

		template <typename TAssFunc> void set_add_jac_A_elem_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_add_jac_M_elem_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_add_def_A_elem_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_add_def_A_elem_fct_explicit(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_add_def_M_elem_fct(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_add_rhs_elem_fct(ReferenceObjectID id, TAssFunc func);

	///	sets usage of fast assemble functions
		void enable_fast_add_elem(bool bEnable) {m_bFastAssembleEnabled = bEnable;}

	///	sets all assemble functions to NULL
		void clear_add_fct();

	public:
	///	returns if fast assembling for elememts is used
		bool fast_add_elem_enabled() const {return m_bFastAssembleEnabled;}

	/// sets the geometric object type
	/**
	 * This functions set the geometric object type of the object, that is
	 * assembled next. The user has to call this function before most of the
	 * assembling routines can be called. Keep in mind, that the elements are
	 * looped type by type, thus this function has to be called very few times.
	 */
		void set_roid(ReferenceObjectID id, int discType);

	private:
	///	flag if fast assemble is used
		bool m_bFastAssembleEnabled;

	// 	timestep function pointers
		PrepareTimestepElemFct 		m_vPrepareTimestepElemFct[NUM_REFERENCE_OBJECTS];
		FinishTimestepElemFct 		m_vFinishTimestepElemFct[NUM_REFERENCE_OBJECTS];

	// 	loop function pointers
		PrepareElemLoopFct 	m_vPrepareElemLoopFct[NUM_REFERENCE_OBJECTS];
		PrepareElemFct 		m_vPrepareElemFct[NUM_REFERENCE_OBJECTS];
		FinishElemLoopFct 	m_vFinishElemLoopFct[NUM_REFERENCE_OBJECTS];

	// 	Jacobian function pointers
		ElemJAFct 	m_vElemJAFct[NUM_REFERENCE_OBJECTS];
		ElemJMFct 	m_vElemJMFct[NUM_REFERENCE_OBJECTS];

	// 	Defect function pointers
		ElemdAFct 	m_vElemdAFct[NUM_REFERENCE_OBJECTS];
		ElemdAFct 	m_vElemdAFct_explicit[NUM_REFERENCE_OBJECTS];
		ElemdMFct 	m_vElemdMFct[NUM_REFERENCE_OBJECTS];

	// 	Rhs function pointers
		ElemRHSFct 	m_vElemRHSFct[NUM_REFERENCE_OBJECTS];

	protected:
	/// current Geometric Object
		ReferenceObjectID m_id;
};

template <typename TDomain>
class IDomainElemDisc : public IElemDisc
{
	private:
	///	base class type
		typedef IElemDisc base_type;

	public:
	///	Domain type
		typedef TDomain domain_type;

	///	World dimension
		static const int dim = TDomain::dim;

	///	Position type
		typedef typename TDomain::position_type position_type;

	public:
	///	Constructor
		IDomainElemDisc(const char* functions = NULL, const char* subsets = NULL)
			: IElemDisc(functions, subsets), m_spApproxSpace(NULL) {};
		
	///	Constructor
		IDomainElemDisc(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset)
			: IElemDisc(vFct, vSubset), m_spApproxSpace(NULL) {};

	///	sets the approximation space
	/**	Calls protected virtual 'approximation_space_changed', when a new approximation space
	 * has been set. Note that 'approximation_space_changed' is only called once if the
	 * same approximation space is set multiple times.*/
		void set_approximation_space(SmartPtr<ApproximationSpace<domain_type> > approxSpace)
		{
		//	check whether the approximation space has already been set
			bool newApproxSpace = (m_spApproxSpace != approxSpace);

		//	remember approx space
			m_spApproxSpace = approxSpace;

		//	invoke callback
			if(newApproxSpace)
				approximation_space_changed();
		}

	///	returns approximation space
		SmartPtr<ApproximationSpace<domain_type> > approx_space() {return m_spApproxSpace;}

	///	returns approximation space
		ConstSmartPtr<ApproximationSpace<domain_type> > approx_space() const {return m_spApproxSpace;}

	///	returns the domain
		domain_type& domain()
		{
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set.");
			return *m_spApproxSpace->domain();
		}

	///	returns the domain
		const domain_type& domain() const
		{
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set.");
			return *m_spApproxSpace->domain();
		}

	///	returns the function pattern
		const FunctionPattern& function_pattern() const {return *m_spApproxSpace;}

	///	returns if function pattern set
		bool fct_pattern_set() const {return m_spApproxSpace.valid();}

	///	returns the subset handler
		typename domain_type::subset_handler_type& subset_handler()
		{
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set.");
			return *m_spApproxSpace->domain()->subset_handler();
		}

	///	returns the subset handler
		const typename domain_type::subset_handler_type& subset_handler() const
		{
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set.");
			return *m_spApproxSpace->domain()->subset_handler();
		}

	///	returns the corner coordinates of an Element in a C-array
		template<typename TElem>
		const position_type* element_corners(TElem* elem)
		{
		//	check domain
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set");
		
		//	get and update the provider
			ElemGlobCornerCoords<TDomain, TElem>& co_coord = Provider<ElemGlobCornerCoords<TDomain, TElem> >::get();
			co_coord.update((m_spApproxSpace->domain()).get(), elem);
			return co_coord.vGlobalCorner();
		}

	protected:
	///	callback invoked, when approximation space is changed
		virtual void approximation_space_changed() {}

	protected:
	///	Approximation Space
		SmartPtr<ApproximationSpace<domain_type> > m_spApproxSpace;

};
/// @}

} // end namespace ug

#include "elem_disc_interface_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__ELEM_DISC__ELEM_DISC_INTERFACE__ */
