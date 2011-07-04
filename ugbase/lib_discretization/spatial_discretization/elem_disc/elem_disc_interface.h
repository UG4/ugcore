/*
 * elem_disc_interface.h
 *
 *  Created on: 07.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__

// extern headers
#include <vector>

// intern headers
#include "lib_discretization/common/local_algebra.h"
#include "lib_discretization/time_discretization/solution_time_series.h"
#include "lib_discretization/function_spaces/approximation_space.h"

namespace ug{

// predeclaration
class IDataExport;
class IDataImport;

/**
 * Element Discretizations
 *
 * \defgroup lib_disc_elem_disc Elem Disc
 * \ingroup lib_discretization
 */

/// \ingroup lib_disc_elem_disc
/// @{


class IElemDisc{
	public:
	///	Local matrix type
		typedef LocalMatrix local_matrix_type;

	/// Local vector type
		typedef LocalVector local_vector_type;

	/// Local index type
		typedef LocalIndices local_index_type;

	public:
	///	Constructor
		IElemDisc()
			: 	m_bTimeDependent(false), m_time(0.0),
			  	m_pLocalVectorTimeSeries(NULL), m_id(-1)
		{}

	////////////////////////////
	// Functions and Subsets
	////////////////////////////

	///	sets functions by name list, divided by ','
		virtual bool set_functions(const char* functions);

	///	sets subset(s) by name list, divided by ','
		virtual bool set_subsets(const char* subsets);

	///	returns the symbolic functions
		const std::vector<std::string>& symb_fcts() const {return m_vFct;}

	///	returns the symbolic subsets
		const std::vector<std::string>& symb_subsets() const {return m_vSubset;}

	protected:
	///	vector holding name of all symbolic functions
		std::vector<std::string> m_vFct;

	///	vector holding name of all symbolic subsets
		std::vector<std::string> m_vSubset;

	////////////////////////////
	// IP Data and Coupling
	////////////////////////////
	public:
	///	registers a data import
		void register_import(IDataImport& Imp);

	///	registers a data export
		void register_export(IDataExport& Exp);

	///	returns number of imports
		size_t num_imports() const {return m_vIImport.size();}

	/// returns an import
		IDataImport& get_import(size_t i);

	///	returns number of exports
		size_t num_exports() const {return m_vIExport.size();}

	/// returns an export
		IDataExport& get_export(size_t i);

	protected:
	/// data imports
		std::vector<IDataImport*> m_vIImport;

	///	data exports
		std::vector<IDataExport*> m_vIExport;

	public:
	////////////////////////////
	// Assembling functions
	////////////////////////////

	/// number of functions this discretization handles
		virtual size_t num_fct() = 0;

	/// shape function set of the functions handled by this discretization
		virtual LocalShapeFunctionSetID local_shape_function_set_id(size_t i) = 0;

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
		virtual bool treat_non_regular_grid(bool bNonRegular) = 0;

	///	returns if discretization acts on hanging nodes if present
	/**
	 * This function returns if a discretization really needs the hanging nodes
	 * in case of non-regular grid. This may not be the case for e.g. finite
	 * element assemblings but is needed for finite volumes
	 */
		virtual bool use_hanging() const {return false;}

	/// sets the geometric object type
	/**
	 * This functions set the geometric object type of the object, that is
	 * assembled next. The user has to call this function before most of the
	 * assembling routines can be called. Keep in mind, that the elements are
	 * looped type by type, thus this function has to be called very few times.
	 */
		bool set_geometric_object_type(int id);

	///	sets if assembling should be time-dependent (and the time point iff)
	/**
	 * This function specifies if the assembling is time-dependent. In that case
	 * the time point is set and can be accessed by the method time(). Iff the
	 * problem is time dependent, the time_point_changed() callback is invoked
	 * with the new time point.
	 *
	 * \param[in]	bTimeDependent	flag if time-dependent
	 * \param[in]	time			time point
	 * \param[in]	time	time point
	 *
	 * \returns 	if elem disc needs time series local solutions
	 */
		bool set_time_dependent(bool bTimeDependent,
		                        number time = 0.0,
		                        const LocalVectorTimeSeries* locTimeSeries = NULL)
		{
			m_bTimeDependent = bTimeDependent;
			m_time = time;
			m_pLocalVectorTimeSeries = locTimeSeries;
			if(is_time_dependent()) return time_point_changed(m_time);
			else return false;
		}

	///	returns if assembling is time-dependent
		bool is_time_dependent() const {return m_bTimeDependent;}

	///	callback, invoked when new time point is set
	/**
	 * This callback must be implemented by a derived Elem Disc in order to handle
	 * time-dependent data. The callback is invoked, whenever the time point
	 * of the disc changes. As return the derived Elem Disc can specify, if
	 * it really needs data from previous time steps for the (spatial) disc. The
	 * default is false.
	 *
	 * \param[in]	time		new time point
	 * \returns 	if elem disc needs time series local solutions
	 */
		virtual bool time_point_changed(number time) {return false;}

	///	returns the current time point
	/**
	 * This function returns the current time point.
	 *
	 * \returns 	time 	time point
	 */
		inline number time() const {return m_time;}

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

	///	prepares the loop over all elements of one type
	/**
	 * This function should prepare the element loop for elements of one fixed
	 * type. This function is called before e.g. the loop over all geometric
	 * objects of a chosen type is performed.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool prepare_element_loop()
			{return (this->*(m_vPrepareElemLoopFunc[m_id]))();}

	///	prepare one elements for assembling
	/**
	 * This function prepares one Geometric object, that will be assembled in
	 * the next step.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 *
	 * \param[in]		obj			The geometric object
	 * \param[in]		u			The current local solution
	 * \param[in]		glob_ind	The global indices of the local solution
	 */
		bool prepare_element(GeometricObject* obj,
		                     const local_vector_type& u,
		                     const local_index_type& glob_ind)
			{return (this->*(m_vPrepareElemFunc[m_id]))(obj, u, glob_ind);}

	///	postprocesses the loop over all elements of one type
	/**
	 * This function should post process the element loop for elements of one fixed
	 * type. This function is called after e.g. the loop over all geometric
	 * objects of a chosen type has been performed.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool finish_element_loop()
			{return (this->*(m_vFinishElemLoopFunc[m_id]))();}

	/// Assembling of Jacobian (Stiffness part)
	/**
	 * This function assembles the local (stiffness) jacobian for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u)
			{return (this->*(m_vAssJAFunc[m_id]))(J, u);}

	/// Assembling of Jacobian (Mass part)
	/**
	 * This function assembles the local (mass) jacobian for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u)
			{return (this->*(m_vAssJMFunc[m_id]))(J, u);}

	/// Assembling of Defect (Stiffness part)
	/**
	 * This function assembles the local (stiffness) defect for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_A(local_vector_type& d, const local_vector_type& u)
			{return (this->*(m_vAssAFunc[m_id]))(d, u);}


	/// Assembling of Defect (Mass part)
	/**
	 * This function assembles the local (mass) defect for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_M(local_vector_type& d, const local_vector_type& u)
			{return (this->*(m_vAssMFunc[m_id]))(d, u);}

	/// Assembling of Right-Hand Side
	/**
	 * This function assembles the local rhs.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_f(local_vector_type& rhs)
			{return (this->*(m_vAssFFunc[m_id]))(rhs);}

	/// Virtual destructor
		virtual ~IElemDisc() {}

	protected:
	///	flag if time dependent
		bool m_bTimeDependent;

	///	time point
		number m_time;

	///	list of local vectors for all solutions of the time series
		const LocalVectorTimeSeries* m_pLocalVectorTimeSeries;

	protected:
	// 	register the functions
		template <typename TAssFunc> void reg_prepare_vol_loop_fct(int id, TAssFunc func);
		template <typename TAssFunc> void reg_prepare_vol_fct(int id, TAssFunc func);
		template <typename TAssFunc> void reg_finish_vol_loop_fct(int id, TAssFunc func);

		template <typename TAssFunc> void reg_ass_JA_vol_fct(int id, TAssFunc func);
		template <typename TAssFunc> void reg_ass_JM_vol_fct(int id, TAssFunc func);
		template <typename TAssFunc> void reg_ass_dA_vol_fct(int id, TAssFunc func);
		template <typename TAssFunc> void reg_ass_dM_vol_fct(int id, TAssFunc func);
		template <typename TAssFunc> void reg_ass_rhs_vol_fct(int id, TAssFunc func);

	protected:
	// 	checks if the needed functions are registered for the id type
		bool function_registered(int id);

	// 	checks if the functions are present
		bool prepare_element_loop_fct_registered(int id);
		bool prepare_element_fct_registered(int id);
		bool finish_element_loop_fct_registered(int id);

	// 	checks if the functions are present
		bool ass_JA_vol_fct_registered(int id);
		bool ass_JM_vol_fct_registered(int id);
		bool ass_dA_vol_fct_registered(int id);
		bool ass_dM_vol_fct_registered(int id);
		bool ass_rhs_vol_fct_registered(int id);

	private:
	//	abbreviation for own type
		typedef IElemDisc T;

	// 	types of loop function pointers
		typedef bool (T::*PrepareElementLoopFunc)();
		typedef bool (T::*PrepareElementFunc)(	GeometricObject* obj,
												const local_vector_type& u,
												const local_index_type& glob_ind);
		typedef bool (T::*FinishElementLoopFunc)();

	// 	types of Jacobian assemble functions
		typedef bool (T::*AssembleJAFunc)(	local_matrix_type& J,
											const local_vector_type& u);
		typedef bool (T::*AssembleJMFunc)(	local_matrix_type& J,
											const local_vector_type& u);

	// 	types of Defect assemble functions
		typedef bool (T::*AssembleAFunc)( 	local_vector_type& d,
											const local_vector_type& u);
		typedef bool (T::*AssembleMFunc)(	local_vector_type& d,
											const local_vector_type& u);

	// 	types of right hand side assemble functions
		typedef bool (T::*AssembleFFunc)(local_vector_type& d);

	private:
	// 	loop function pointers
		std::vector<PrepareElementLoopFunc> m_vPrepareElemLoopFunc;
		std::vector<PrepareElementFunc> 	m_vPrepareElemFunc;
		std::vector<FinishElementLoopFunc> 	m_vFinishElemLoopFunc;

	// 	Jacobian function pointers
		std::vector<AssembleJAFunc> 	m_vAssJAFunc;
		std::vector<AssembleJMFunc> 	m_vAssJMFunc;

	// 	Defect function pointers
		std::vector<AssembleAFunc> 	m_vAssAFunc;
		std::vector<AssembleMFunc> 	m_vAssMFunc;

	// 	Rhs function pointers
		std::vector<AssembleFFunc> 	m_vAssFFunc;

	protected:
		// current Geometric Object
		int m_id;
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

	///	Local matrix type
		typedef typename base_type::local_matrix_type local_matrix_type;

	/// Local vector type
		typedef typename base_type::local_vector_type local_vector_type;

	/// Local index type
		typedef typename base_type::local_index_type local_index_type;

	public:
		IDomainElemDisc() : m_pDomain(NULL) {};

	///	sets the approximation space
		void set_approximation_space(IApproximationSpace<domain_type>& approxSpace)
		{
		//	remember approx space
			m_pApproxSpace = &approxSpace;

		//	remember domain
			set_domain(approxSpace.get_domain());
		}

	///	sets the domain
		void set_domain(domain_type& domain)
		{
		//	remember domain
			m_pDomain = &domain;

		//	remember position accessor
			m_aaPos = m_pDomain->get_position_accessor();
		}

	///	returns the domain
		domain_type& get_domain()
		{
			UG_ASSERT(m_pDomain != NULL, "Domain not set.");
			return *m_pDomain;
		}

	///	returns the domain
		const domain_type& get_domain() const
		{
			UG_ASSERT(m_pDomain != NULL, "Domain not set.");
			return *m_pDomain;
		}

	///	returns the function pattern
		const FunctionPattern& get_fct_pattern() const {return *m_pApproxSpace;}

	///	returns the subset handler
		typename domain_type::subset_handler_type& get_subset_handler()
		{
			UG_ASSERT(m_pDomain != NULL, "Domain not set.");
			return m_pDomain->get_subset_handler();
		}

	///	returns the subset handler
		const typename domain_type::subset_handler_type& get_subset_handler() const
		{
			UG_ASSERT(m_pDomain != NULL, "Domain not set.");
			return m_pDomain->get_subset_handler();
		}

	///	returns the corner coordinates of an Element in a C++-array
		template<typename TElem>
		const position_type* get_element_corners(TElem* elem)
		{
			typedef typename reference_element_traits<TElem>::reference_element_type
							ref_elem_type;

		//	resize corners
			m_vCornerCoords.resize(ref_elem_type::num_corners);

		//	check domain
			UG_ASSERT(m_pDomain != NULL, "Domain not set");

		//	extract corner coordinates
			for(size_t i = 0; i < m_vCornerCoords.size(); ++i)
			{
				VertexBase* vert = elem->vertex(i);
				m_vCornerCoords[i] = m_aaPos[vert];
			}

		//	return corner coords
			return &m_vCornerCoords[0];
		}

	protected:
	///	Domain
		TDomain* m_pDomain;

	///	Position access
		typename TDomain::position_accessor_type m_aaPos;

	///	Corner Coordinates
		std::vector<position_type> m_vCornerCoords;

	///	Approximation Space
		IApproximationSpace<domain_type>* m_pApproxSpace;
};
/// @}

} // end namespace ug

#include "elem_disc_interface_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__ */
