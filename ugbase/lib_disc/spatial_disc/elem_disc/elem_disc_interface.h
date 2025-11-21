/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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
#include "lib_disc/domain_traits.h"
#include "elem_modifier.h"
#include "lib_disc/spatial_disc/elem_disc/err_est_data.h"
#include "bridge/util_algebra_dependent.h"
#include "lib_disc/common/multi_index.h"

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


/// Proxy struct for generic passing of any vector type
struct VectorProxyBase
{
	virtual ~VectorProxyBase() {};
	virtual number evaluate(const DoFIndex& di) const = 0;
};

template <typename TVector>
struct VectorProxy : public VectorProxyBase
{
	VectorProxy(const TVector& v) : m_v(v) {}

	virtual number evaluate(const DoFIndex& di) const {return DoFRef(m_v, di);}

	const TVector& m_v;

};

/**
 * Element Discretizations
 *
 * \defgroup lib_disc_elem_disc Elem Disc
 * \ingroup lib_discretization
 */

/// \ingroup lib_disc_elem_disc
/// @{
/*
template <typename TDomain>
class IElemDiscBaseData
{
public:
	///	Domain type
	using domain_type = TDomain;

	///	World dimension
	static constexpr int dim = TDomain::dim;
};
*/
/// This class encapsulates all functions related to error estimation
template <typename TLeaf, typename TDomain>
class IElemAssembleFuncs
{
public:
	/// constructor
	IElemAssembleFuncs() { set_default_add_fct(); }

	/// Virtual destructor
	virtual ~IElemAssembleFuncs()= default;

	/// Barton Nackman trick (TODO: needed?)
	using leaf_type = TLeaf;

	TLeaf& asLeaf()
	{ return static_cast<TLeaf&>(*this); }

	///	Domain type
	using domain_type = TDomain;

	///	World dimension
	static constexpr int dim = TDomain::dim;

	////////////////////////////
	// assembling functions
	////////////////////////////
public:
	///	virtual prepares the loop over all elements of one type
	virtual void prep_assemble_loop() {}

	///	virtual prepares the loop over all elements of one type
	virtual void post_assemble_loop() {}

	/// prepare the time step
	virtual void prep_timestep(number future_time, number time, VectorProxyBase* u);

	/// prepare the time step element-wise
	virtual void prep_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	virtual prepares the loop over all elements of one type
	virtual void prep_elem_loop(const ReferenceObjectID roid, const int si);

	///	virtual prepare one elements for assembling
	virtual void prep_elem(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	///	virtual postprocesses the loop over all elements of one type
	virtual void fsh_elem_loop();

	/// finish the time step
	virtual void fsh_timestep(number time, VectorProxyBase* u);

	/// virtual finish the time step element-wise
	virtual void fsh_timestep_elem(const number time, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// Assembling of Jacobian (Stiffness part)
	virtual void add_jac_A_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// Assembling of Jacobian (Mass part)
	virtual void add_jac_M_elem(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// virtual Assembling of Defect (Stiffness part)
	virtual void add_def_A_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// defect for explicit terms
	virtual void add_def_A_expl_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// virtual Assembling of Defect (Mass part)
	virtual void add_def_M_elem(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	/// virtual Assembling of Right-Hand Side
	virtual void add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[]);


	///	function dispatching call to implementation
	/// \{
	void do_prep_timestep(number future_time, const number time, VectorProxyBase* u, size_t algebra_id);
	void do_prep_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_prep_elem_loop(const ReferenceObjectID roid, const int si);
	void do_prep_elem(LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);
	void do_fsh_elem_loop();
	void do_fsh_timestep(const number time, VectorProxyBase* u, size_t algebra_id);
	void do_fsh_timestep_elem(const number time, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_add_jac_A_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_add_jac_M_elem(LocalMatrix& J, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_add_def_A_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_add_def_A_expl_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_add_def_M_elem(LocalVector& d, LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_add_rhs_elem(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[]);




protected:
	// 	register the functions
	template <typename TAssFunc> void set_prep_timestep_fct(size_t algebra_id, TAssFunc func);
	template <typename TAssFunc> void set_prep_timestep_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_fsh_timestep_fct(size_t algebra_id, TAssFunc func);
	template <typename TAssFunc> void set_fsh_timestep_elem_fct(ReferenceObjectID id, TAssFunc func);

	template <typename TAssFunc> void set_prep_elem_loop_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_prep_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_fsh_elem_loop_fct(ReferenceObjectID id, TAssFunc func);

	template <typename TAssFunc> void set_add_jac_A_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_add_jac_M_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_add_def_A_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_add_def_A_expl_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_add_def_M_elem_fct(ReferenceObjectID id, TAssFunc func);
	template <typename TAssFunc> void set_add_rhs_elem_fct(ReferenceObjectID id, TAssFunc func);



	//	unregister functions
	void remove_prep_timestep_fct(size_t algebra_id);
	void remove_prep_timestep_elem_fct(ReferenceObjectID id);
	void remove_fsh_timestep_fct(size_t algebra_id);
	void remove_fsh_timestep_elem_fct(ReferenceObjectID id);

	void remove_prep_elem_loop_fct(ReferenceObjectID id);
	void remove_prep_elem_fct(ReferenceObjectID id);
	void remove_fsh_elem_loop_fct(ReferenceObjectID id);

	void remove_add_jac_A_elem_fct(ReferenceObjectID id);
	void remove_add_jac_M_elem_fct(ReferenceObjectID id);
	void remove_add_def_A_elem_fct(ReferenceObjectID id);
	void remove_add_def_A_expl_elem_fct(ReferenceObjectID id);
	void remove_add_def_M_elem_fct(ReferenceObjectID id);
	void remove_add_rhs_elem_fct(ReferenceObjectID id);

protected:
	///	sets all assemble functions to the corresponding virtual ones
	void set_default_add_fct();

	///	sets all assemble functions to nullptr for a given ReferenceObjectID
	void clear_add_fct(ReferenceObjectID id);

	///	sets all assemble functions to nullptr (for all ReferenceObjectID's)
	void clear_add_fct();

private:
//	abbreviation for own type
	using T = IElemAssembleFuncs<TLeaf, TDomain>;

// 	types of timestep function pointers
	using PrepareTimestepFct = void(T::*)(number, number, VectorProxyBase*);

	using PrepareTimestepElemFct = void(T::*)(number, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	using FinishTimestepFct = void(T::*)(number, VectorProxyBase*);

	using FinishTimestepElemFct = void(T::*)(number, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

// 	types of loop function pointers
	using PrepareElemLoopFct = void(T::*)(ReferenceObjectID roid, int si);

	using PrepareElemFct = void(T::*)(const LocalVector& u, GridObject* elem, const ReferenceObjectID roid, const MathVector<dim> vCornerCoords[]);

	using FinishElemLoopFct = void(T::*)();

// 	types of Jacobian assemble functions
	using ElemJAFct = void(T::*)(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	using ElemJMFct = void(T::*)(LocalMatrix& J, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

// 	types of Defect assemble functions
	using ElemdAFct = void(T::*)(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	using ElemdMFct = void(T::*)(LocalVector& d, const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

// 	types of right hand side assemble functions
	using ElemRHSFct = void(T::*)(LocalVector& rhs, GridObject* elem, const MathVector<dim> vCornerCoords[]);


private:
// 	timestep function pointers
	PrepareTimestepFct			m_vPrepareTimestepFct[bridge::NUM_ALGEBRA_TYPES];
	PrepareTimestepElemFct 		m_vPrepareTimestepElemFct[NUM_REFERENCE_OBJECTS];
	FinishTimestepFct			m_vFinishTimestepFct[bridge::NUM_ALGEBRA_TYPES];
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
	ElemdAFct 	m_vElemdAExplFct[NUM_REFERENCE_OBJECTS];
	ElemdMFct 	m_vElemdMFct[NUM_REFERENCE_OBJECTS];

// 	Rhs function pointers
	ElemRHSFct 	m_vElemRHSFct[NUM_REFERENCE_OBJECTS];

public:
/// sets the geometric object type
/**
 * This functions set the geometric object type of the object, that is
 * assembled next. The user has to call this function before most of the
 * assembling routines can be called. Keep in mind, that the elements are
 * looped type by type, thus this function has to be called very few times.
 */
	void set_roid(ReferenceObjectID id, int discType);

	/// check, if all inputs have been set
	void check_roid(ReferenceObjectID roid, int discType);

protected:
/// current Geometric Object
	ReferenceObjectID m_roid;
};




/// This class encapsulates all functions related to error estimation
template <typename TLeaf, typename TDomain>
class IElemEstimatorFuncs
{
public:
	/// constructor
	IElemEstimatorFuncs() : m_bDoErrEst(false), m_spErrEstData(nullptr)
	{ set_default_add_fct(); }

	/// Virtual destructor
	virtual ~IElemEstimatorFuncs()= default;

	/// Barton Nackman trick (TODO: needed?)
	using leaf_type = TLeaf;

	TLeaf& asLeaf()
	{ return static_cast<TLeaf&>(*this); }

	///	Domain type
	using domain_type = TDomain;

	///	World dimension
	static constexpr int dim = TDomain::dim;

	void do_prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);
	void do_prep_err_est_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);
	void do_compute_err_est_A_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);
	void do_compute_err_est_M_elem(LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);
	void do_compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);
	void do_fsh_err_est_elem_loop();

public:
	///	virtual prepares the loop over all elements of one type for the computation of the error estimator
		virtual void prep_err_est_elem_loop(const ReferenceObjectID roid, const int si);

	///	virtual prepares the loop over all elements of one type for the computation of the error estimator
		virtual void prep_err_est_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	///	virtual compute the error estimator (stiffness part) contribution for one element
		virtual void compute_err_est_A_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	virtual compute the error estimator (mass part) contribution for one element
		virtual void compute_err_est_M_elem(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	virtual compute the error estimator (rhs part) contribution for one element
		virtual void compute_err_est_rhs_elem(GridObject* elem, const MathVector<dim> vCornerCoords[], const number& scale);

	///	virtual postprocesses the loop over all elements of one type in the computation of the error estimator
		virtual void fsh_err_est_elem_loop();

protected:
		template <typename TAssFunc> void set_prep_err_est_elem_loop(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_prep_err_est_elem(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_compute_err_est_A_elem(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_compute_err_est_M_elem(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_compute_err_est_rhs_elem(ReferenceObjectID id, TAssFunc func);
		template <typename TAssFunc> void set_fsh_err_est_elem_loop(ReferenceObjectID id, TAssFunc func);

		void remove_prep_err_est_elem_loop(ReferenceObjectID id);
		void remove_prep_err_est_elem(ReferenceObjectID id);
		void remove_compute_err_est_A_elem(ReferenceObjectID id);
		void remove_compute_err_est_M_elem(ReferenceObjectID id);
		void remove_compute_err_est_rhs_elem(ReferenceObjectID id);
		void remove_fsh_err_est_elem_loop(ReferenceObjectID id);

		///	sets all assemble functions to nullptr for a given ReferenceObjectID
		void clear_add_fct(ReferenceObjectID id);

		///	sets all assemble functions to nullptr (for all ReferenceObjectID's)
		void clear_add_fct();

		///	sets all assemble functions to the corresponding virtual ones
		void set_default_add_fct();



private:
	//	abbreviation for own type
	using T = IElemEstimatorFuncs<TLeaf, TDomain>;

	//	types of the error estimator assembler
	using PrepareErrEstElemLoopFct = void(T::*)(ReferenceObjectID roid, int si);

	using PrepareErrEstElemFct = void(T::*)(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[]);

	using ElemComputeErrEstAFct = void(T::*)(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number&);

	using ElemComputeErrEstMFct = void(T::*)(const LocalVector& u, GridObject* elem, const MathVector<dim> vCornerCoords[], const number&);

	using ElemComputeErrEstRhsFct = void(T::*)(GridObject* elem, const MathVector<dim> vCornerCoords[], const number&);

	using FinishErrEstElemLoopFct = void(T::*)();

	//	Error estimator functions
	PrepareErrEstElemLoopFct	m_vPrepareErrEstElemLoopFct[NUM_REFERENCE_OBJECTS];
	PrepareErrEstElemFct		m_vPrepareErrEstElemFct[NUM_REFERENCE_OBJECTS];
	ElemComputeErrEstAFct		m_vElemComputeErrEstAFct[NUM_REFERENCE_OBJECTS];
	ElemComputeErrEstMFct		m_vElemComputeErrEstMFct[NUM_REFERENCE_OBJECTS];
	ElemComputeErrEstRhsFct		m_vElemComputeErrEstRhsFct[NUM_REFERENCE_OBJECTS];
	FinishErrEstElemLoopFct		m_vFinishErrEstElemLoopFct[NUM_REFERENCE_OBJECTS];


	// //////////////////////////
	// Error estimator
	// //////////////////////////
public:
	///	sets the pointer to an error estimator data object (or nullptr)
	/**
	 * This function sets the pointer to an error estimator data object
	 * that should be used for this discretization. Note that the ElemDisc
	 * object must use RTTI to try to convert this pointer to the type
	 * of the objects accepted by it for this purpose. If the conversion
	 * fails than an exception must be thrown since this situation is not
	 * allowed.
	 */
	void set_error_estimator(SmartPtr<IErrEstData<TDomain> > ee) {m_spErrEstData = ee; m_bDoErrEst = true;}

	/// find out whether or not a posteriori error estimation is to be performed for this disc
	bool err_est_enabled() const {return m_bDoErrEst;}

	///	returns the pointer to the error estimator data object (or nullptr)
	virtual SmartPtr<IErrEstData<TDomain> > err_est_data() {return m_spErrEstData;}

private:
	/// flag indicating whether or not a posteriori error estimation is to be performed for this disc
	bool m_bDoErrEst;

protected:
	/// error estimation object associated to the element discretization
	SmartPtr<IErrEstData<TDomain> > m_spErrEstData;

public:
/// sets the geometric object type
/**
 * This functions set the geometric object type of the object, that is
 * assembled next. The user has to call this function before most of the
 * assembling routines can be called. Keep in mind, that the elements are
 * looped type by type, thus this function has to be called very few times.
 */
	void set_roid(ReferenceObjectID id, int discType);

	/// check, if all inputs have been set
	void check_roid(ReferenceObjectID roid, int discType);

protected:
/// current Geometric Object
	ReferenceObjectID m_roid;
};

///	base class for all element-wise discretizations
/**
 * This class is the base class for element-wise discretizations. An
 * implementation of this class must provide local stiffness/mass-matrix
 * contribution of one element to the global jacobian and local contributions
 * of one element to the local defect.
 */

template <typename TDomain>
class IElemDiscBase
{
	public:
	///	Domain type
		using domain_type = TDomain;

	///	Position type
		using position_type = typename TDomain::position_type;

	///	World dimension
		static constexpr int dim = TDomain::dim;
		
	public:
	///	Constructor
		IElemDiscBase(const char* functions = "", const char* subsets = "");

	///	Constructor
		IElemDiscBase(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset);

	/// Virtual destructor
		virtual ~IElemDiscBase()= default;

	public:
	///	sets the approximation space
	/**	Calls protected virtual 'approximation_space_changed', when a new approximation space
	 * has been set. Note that 'approximation_space_changed' is only called once if the
	 * same approximation space is set multiple times.*/
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	///	returns approximation space
		SmartPtr<ApproximationSpace<TDomain> > approx_space() {return m_spApproxSpace;}

	///	returns approximation space
		ConstSmartPtr<ApproximationSpace<TDomain> > approx_space() const {return m_spApproxSpace;}

	///	returns the domain
		SmartPtr<TDomain> domain(){return m_spApproxSpace->domain();}

	///	returns the domain
		ConstSmartPtr<TDomain> domain() const{return m_spApproxSpace->domain();}

	///	returns the subset handler
		typename TDomain::subset_handler_type& subset_handler()
		{
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set.");
			return *m_spApproxSpace->domain()->subset_handler();
		}

	///	returns the subset handler
		const typename TDomain::subset_handler_type& subset_handler() const
		{
			UG_ASSERT(m_spApproxSpace.valid(), "ApproxSpace not set.");
			return *m_spApproxSpace->domain()->subset_handler();
		}
/*
		void add_elem_modifier(SmartPtr<IElemDiscModifier<TDomain> > elemModifier )
		{
			m_spElemModifier.push_back(elemModifier);
			elemModifier->set_elem_disc(this);
		}
		std::vector<SmartPtr<IElemDiscModifier<TDomain> > >& get_elem_modifier(){ return m_spElemModifier;}
		*/

	protected:
	///	callback invoked, when approximation space is changed
		virtual void approximation_space_changed() {}

	///	Approximation Space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	Approximation Space
	//	std::vector<SmartPtr<IElemDiscModifier<TDomain> > > m_spElemModifier;

	////////////////////////////
	// Functions and Subsets
	////////////////////////////
	public:
	///	sets functions by name list, divided by ','
		void set_functions(const std::string& functions);

	/// sets functions by vector of names
		void set_functions(const std::vector<std::string>& functions);

	///	sets subset(s) by name list, divided by ','
		void set_subsets(const std::string& subsets);

	///	sets subset(s) by name list, divided by ','
		void set_subsets(const std::vector<std::string>& subsets);

	/// number of functions this discretization handles
		size_t num_fct() const {return m_vFct.size();}

	///	returns the symbolic functions
		const std::vector<std::string>& symb_fcts() const {return m_vFct;}

	/// number of subsets this discretization handles
		size_t num_subsets() const {return m_vSubset.size();}

	///	returns the symbolic subsets
		const std::vector<std::string>& symb_subsets() const {return m_vSubset;}

	///	returns the current function pattern
		ConstSmartPtr<FunctionPattern> function_pattern() const {return m_spFctPattern;}

	///	returns the current function group
		const FunctionGroup& function_group() const {return m_fctGrp;}

	///	returns the current function index mapping
		const FunctionIndexMapping& map() const {return m_fctIndexMap;}

	///	checks the setup of the elem disc
		void check_setup(bool bNonRegularGrid);

	protected:
	///	vector holding name of all symbolic functions
		std::vector<std::string> m_vFct;

	///	vector holding name of all symbolic subsets
		std::vector<std::string> m_vSubset;

	///	current function pattern
		ConstSmartPtr<FunctionPattern> m_spFctPattern;

	///	current function group
		FunctionGroup m_fctGrp;

	///	current function index mapping
		FunctionIndexMapping m_fctIndexMap;

	///	sets current function pattern
		void set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt);

	///	updates the function index mapping
		void update_function_index_mapping();

	////////////////////////////
	// UserData and Coupling
	////////////////////////////
	public:
	///	registers a data import
		void register_import(IDataImport<dim>& Imp);

	///	returns number of imports
		size_t num_imports() const {return m_vIImport.size();}

	/// returns an import
		IDataImport<dim>& get_import(size_t i)
		{
			UG_ASSERT(i < num_imports(), "Invalid index");
			return *m_vIImport[i];
		}

	///	removes all imports
		void clear_imports() {m_vIImport.clear();}

	protected:
	/// data imports
		std::vector<IDataImport<dim>*> m_vIImport;

	////////////////////////////
	// time handling
	////////////////////////////
	public:
	///	sets if assembling should be time-dependent and the local time series
	/**
	 * This function specifies if the assembling is time-dependent. If nullptr is
	 * passed, the assembling is assumed to be time-independent. If a local
	 * time series is passed, this series is used as previous solution.
	 *
	 * \param[in]	locTimeSeries	Time series of previous solutions
	 */
		void set_time_dependent(LocalVectorTimeSeries& locTimeSeries,
		        				const std::vector<number>& vScaleMass,
		        				const std::vector<number>& vScaleStiff);

	///	sets that the assembling is time independent
		void set_time_independent();

	///	returns if assembling is time-dependent
		bool is_time_dependent() const {return (m_pLocalVectorTimeSeries != nullptr) && !m_bStationaryForced;}

	///	sets that the assembling is always stationary (even in instationary case)
		void set_stationary(bool bStationaryForced = true) {m_bStationaryForced = bStationaryForced;}
		void set_stationary() {set_stationary(true);}

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

		bool local_time_series_needed() {return is_time_dependent() && requests_local_time_series();}

	///	sets the current time point
		void set_time_point(const size_t timePoint) {m_timePoint = timePoint;}

	///	returns the currently considered time point of the time-disc scheme
		size_t time_point() const {return m_timePoint;}

	///	returns currently set timepoint
		number time() const
		{
			if(m_pLocalVectorTimeSeries) return m_pLocalVectorTimeSeries->time(m_timePoint);
			else return 0.0;
		}

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

	protected:
	///	time point
		size_t m_timePoint;

	///	list of local vectors for all solutions of the time series
		LocalVectorTimeSeries* m_pLocalVectorTimeSeries;

	///	weight factors for time dependent assembling
	/// \{
		std::vector<number> m_vScaleMass;
		std::vector<number> m_vScaleStiff;
	/// \}

	///	flag if stationary assembling is to be used even in instationary assembling
		bool m_bStationaryForced;

	

	////////////////////////////
	// general info
	////////////////////////////
	public:
	///	 returns the type of elem disc
		virtual int type() const {return EDT_ELEM | EDT_SIDE;}

	/// requests assembling for trial spaces and grid type
	/**
	 * This function is called before the assembling starts. The
	 * IElemDisc-Implementation is supposed to checks if it can assemble the set
	 * of LFEID and the grid type. It may register corresponding assembling
	 * functions or perform other initialization.
	 * If the ElemDisc does not support the setting it should throw an exception.
	 *
	 * \param[in] vLfeID			vector of Local Finite Element IDs
	 * \param[in] bNonRegularGrid	regular grid type
	 */
		virtual void prepare_setting(const std::vector<LFEID>& vLfeID, bool bNonRegularGrid) = 0;

	///	returns if discretization acts on hanging nodes if present
	/**
	 * This function returns if a discretization really needs the hanging nodes
	 * in case of non-regular grid. This may not be the case for e.g. finite
	 * element assemblings but is needed for finite volumes
	 */
		virtual bool use_hanging() const {return false;}
};


template <typename TDomain>
class IElemError :
	public IElemDiscBase<TDomain>,
	public IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>
{
public:
	using domain_type = TDomain;
	static constexpr int dim = TDomain::dim;

	friend class IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>;
	using estimator_base_type = IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>;

	IElemError(const char* functions, const char* subsets)
	: IElemDiscBase<TDomain>(functions, subsets), estimator_base_type()  {}

	IElemError(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset)
	: IElemDiscBase<TDomain>(vFct, vSubset), estimator_base_type()  {}


protected:

	///	sets all assemble functions to nullptr for a given ReferenceObjectID
	void clear_add_fct(ReferenceObjectID id)
	{ estimator_base_type::clear_add_fct(id); }

	///	sets all assemble functions to nullptr (for all ReferenceObjectID's)
	void clear_add_fct()
	{ estimator_base_type::clear_add_fct(); }

	///	sets all assemble functions to the corresponding virtual ones
	//void set_default_add_fct();
	using estimator_base_type::set_default_add_fct;

};


/**
 * Element discretization (including error indicator)
 * TODO: Should be separated!!!
 * */
template <typename TDomain>
class IElemDisc :
		public IElemAssembleFuncs<IElemDisc<TDomain>, TDomain>,
		public IElemError<TDomain>
{
public:
	using domain_type = TDomain;
	static constexpr int dim = TDomain::dim;

	/// real base class
	using base_type = IElemError<TDomain>;
	using estimator_base_type = IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>;
	using assemble_base_type = IElemAssembleFuncs<IElemDisc<TDomain>, TDomain>;

	friend class IElemEstimatorFuncs<IElemDisc<TDomain>, TDomain>;
	friend class IElemAssembleFuncs<IElemDisc<TDomain>, TDomain>;


	IElemDisc(const char* functions, const char* subsets)
	: assemble_base_type(), IElemError<TDomain>(functions, subsets) {}

	IElemDisc(const std::vector<std::string>& vFct, const std::vector<std::string>& vSubset)
	: assemble_base_type(), IElemError<TDomain>(vFct, vSubset) {}

protected:

	///	sets all assemble functions to nullptr for a given ReferenceObjectID
		void clear_add_fct(ReferenceObjectID id)
		{
			base_type::clear_add_fct(id);
			assemble_base_type::clear_add_fct(id);
		}

	///	sets all assemble functions to nullptr (for all ReferenceObjectID's)
		void clear_add_fct()
		{
			base_type::clear_add_fct();
			assemble_base_type::clear_add_fct();
		}

	///	sets all assemble functions to the corresponding virtual ones
		void set_default_add_fct()
		{
			base_type::set_default_add_fct();
			assemble_base_type::set_default_add_fct();
		}


public:
	void add_elem_modifier(SmartPtr<IElemDiscModifier<TDomain> > elemModifier )
	{
			m_spElemModifier.push_back(elemModifier);
			elemModifier->set_elem_disc(this);
	}

	std::vector<SmartPtr<IElemDiscModifier<TDomain> > >& get_elem_modifier()
	{ return m_spElemModifier;}

protected:
	///	Approximation Space
	std::vector<SmartPtr<IElemDiscModifier<TDomain> > > m_spElemModifier;


};



/// @}

} // end namespace ug

#include "elem_disc_interface_impl.h"

#endif