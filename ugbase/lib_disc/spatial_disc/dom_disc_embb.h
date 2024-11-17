/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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
 * Global assembling of the problems with the embedded boundary.
 */
#ifndef __H__UG__PLUGINS__D3F__EMBASS__
#define __H__UG__PLUGINS__D3F__EMBASS__

#include <vector>

// ug4 headers
#include "common/common.h"
#include "common/util/smart_pointer.h"
#include "lib_disc/domain_traits.h"
#include "lib_disc/spatial_disc/domain_disc.h"
#include "lib_disc/operator/linear_operator/std_injection.h"
#include "bridge/util_algebra_dependent.h"
#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug {

///	Base class for the extrapolation over an embedded boundary
/**
 * This class provides an interface for the access to the extrapolation over
 * an embedded boundary.
 *
 * \tparam TDomain	domain type
 * \tparam TAlgebra	algebra type for the functions to extrapolate
 */
template <typename TDomain, typename TAlgebra>
class IInterfaceExtrapolation
{
public:

///	domain type
	typedef TDomain domain_type;
	
///	algebra type for the functions to extrapolate
	typedef TAlgebra algebra_type;
	
///	vector type (for the functions to extrapolate)
	typedef typename algebra_type::vector_type vector_type;
	
///	matrix type
	typedef typename algebra_type::matrix_type matrix_type;
	
///	dimensionality (the World dimension)
	static const int dim = domain_type::dim;
	
///	Constructor
	IInterfaceExtrapolation () {}
	
///	Destructor
	virtual ~IInterfaceExtrapolation () {}

///	checks whether the element is intersected by the interface, or what, and prepares the data
	virtual int check_elem_lsf
	(
		size_t n_co, ///< number of the corners of the element
		GridObject * pElem, ///< the element to process
		int si, ///< subset of the element
		int g_level, ///< grid level of the element
		bool use_hanging, ///< if there can be hanging nodes
		const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
		number time ///< the phisical time
	)
	{UG_THROW ("IInterfaceExtrapolation: Virtual functions are not implemented in the base class.");}
	
///	(slower version) checks whether the element is intersected by the interface, or what, and prepares the data
	virtual int check_elem_lsf
	(
		size_t n_co, ///< number of the corners of the element
		GridObject * pElem, ///< the element to process
		int si, ///< subset of the element
		bool use_hanging, ///< if there can be hanging nodes
		const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
		number time ///< the phisical time
	)
	{UG_THROW ("IInterfaceExtrapolation: Virtual functions are not implemented in the base class.");}
	
///	extrapolates a component of the solution to the vertices behind the interface (w.r.t. a base corner)
	virtual void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		size_t base_co, ///< the base corner
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const
	{UG_THROW ("IInterfaceExtrapolation: Virtual functions are not implemented in the base class.");}
	
///	extrapolates a component of the solution to the vertices behind the interface (by averaging)
	virtual void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const
	{UG_THROW ("IInterfaceExtrapolation: Virtual functions are not implemented in the base class.");}
	
///	returns true if the corner is "inside" (use after check_elem_lsf)
	virtual bool corner_inside
	(
		size_t co ///< the corner
	) const
	{UG_THROW ("IInterfaceExtrapolation: Virtual functions are not implemented in the base class.");}
	
///	returns the effective value of the LSF at a corner (use after check_elem_lsf)
	virtual number lsf_at
	(
		size_t co ///< the corner
	) const
	{UG_THROW ("IInterfaceExtrapolation: Virtual functions are not implemented in the base class.");}

}; // class IInterfaceExtrapolation

/// Global assembler based on the ghost-fluid method with a level-set function
/**
 * Template class of the global assembler based on the ghost-fluid method
 * with a piecewise linear level-set function.
 *
 * \tparam TDomain			domain type
 * \tparam TAlgebra			algebra type
 * \tparam TExtrapolation	extrapolation class
 */
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
class LSGFGlobAssembler
{
public:

///	Domain type
	typedef TDomain domain_type;
	
///	Algebra type
	typedef TAlgebra algebra_type;
	
///	type of approximation space
	typedef ApproximationSpace<domain_type> approx_space_type;
	
///	Vector type in the algebra
	typedef typename algebra_type::vector_type vector_type;
	
///	Matrix type in the algebra
	typedef typename algebra_type::matrix_type matrix_type;
	
///	Extrapolation type
	typedef TExtrapolation extrapolation_type;
	
///	Grid function type for the LSF
	typedef typename extrapolation_type::ls_grid_func_type ls_grid_func_type;
	
///	world dimension
	static const int dim = TDomain::dim;
	
////////////////////////////////////////////////////////////////////////////////
// Constructor/destructor
////////////////////////////////////////////////////////////////////////////////

public:

///	class constructor (may not have any arguments!)
	LSGFGlobAssembler () : m_bAssembleOnlyCut(false) {};
	
///	virtual destructor
	virtual ~LSGFGlobAssembler () {};

////////////////////////////////////////////////////////////////////////////////
// Assembling tools
////////////////////////////////////////////////////////////////////////////////

public:

	template <typename TElem, typename TIterator>
	void
	AssembleStiffnessMatrix(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
								ConstSmartPtr<domain_type> spDomain,
								ConstSmartPtr<DoFDistribution> dd,
								TIterator iterBegin,
								TIterator iterEnd,
								int si, bool bNonRegularGrid,
								matrix_type& A,
								const vector_type& u,
								ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleMassMatrix( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<domain_type> spDomain,
						ConstSmartPtr<DoFDistribution> dd,
						TIterator iterBegin,
						TIterator iterEnd,
						int si, bool bNonRegularGrid,
						matrix_type& M,
						const vector_type& u,
						ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleJacobian(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<domain_type> spDomain,
						ConstSmartPtr<DoFDistribution> dd,
						TIterator iterBegin,
						TIterator iterEnd,
						int si, bool bNonRegularGrid,
						matrix_type& J,
						const vector_type& u,
						ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleJacobian(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
						ConstSmartPtr<domain_type> spDomain,
						ConstSmartPtr<DoFDistribution> dd,
						TIterator iterBegin,
						TIterator iterEnd,
						int si, bool bNonRegularGrid,
						matrix_type& J,
						ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
						number s_a0,
						ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleDefect( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& d,
					const vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleDefect( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& d,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					const std::vector<number>& vScaleMass,
					const std::vector<number>& vScaleStiff,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleLinear( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					matrix_type& A,
					vector_type& rhs,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	AssembleLinear( const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					matrix_type& A,
					vector_type& rhs,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					const std::vector<number>& vScaleMass,
					const std::vector<number>& vScaleStiff,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

////////////////////////////////////////////////////////////////////////////////
// Assemble Rhs: it cannot be done for the ghost-fluid method independently of the matrix
////////////////////////////////////////////////////////////////////////////////

public:

	template <typename TElem, typename TIterator>
	static void
	AssembleRhs(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& rhs,
					const vector_type& u,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
		UG_THROW ("LSGFGlobAssembler::AssembleRhs: Cannot assemble the RHS in GF independently of the matrix");
	}

	template <typename TElem, typename TIterator>
	static void
	AssembleRhs(	const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					vector_type& rhs,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					const std::vector<number>& vScaleMass,
					const std::vector<number>& vScaleStiff,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner)
	{
		UG_THROW ("LSGFGlobAssembler::AssembleRhs: Cannot assemble the RHS in GF independently of the matrix");
	}

////////////////////////////////////////////////////////////////////////////////
// Prepare and Finish Timestep: these version merely skip the outer elements
////////////////////////////////////////////////////////////////////////////////

public:
	/**
	 * This function prepares the global discretization for a time-stepping scheme
	 * by calling the "prepare_timestep" methods of all passed element
	 * discretizations.
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		vSol			current and previous solutions
	 * \param[in]		spAssTuner		assemble adapter
	 */
	void PrepareTimestep
	(
		const std::vector<IElemDisc<domain_type>*>& vElemDisc,
		ConstSmartPtr<DoFDistribution> dd,
		bool bNonRegularGrid,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		number future_time,
		ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner
	);

	template <typename TElem, typename TIterator>
	void
	PrepareTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
					ConstSmartPtr<domain_type> spDomain,
					ConstSmartPtr<DoFDistribution> dd,
					TIterator iterBegin,
					TIterator iterEnd,
					int si, bool bNonRegularGrid,
					ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
					ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);
	/**
	 * This function finishes the global discretization for a time-stepping scheme
	 * by calling the "finish_timestep" methods of all passed element
	 * discretizations.
	 *
	 * \param[in]		vElemDisc		element discretizations
	 * \param[in]		dd				DoF Distribution
	 * \param[in]		bNonRegularGrid flag to indicate if non regular grid is used
	 * \param[in]		vSol			current and previous solutions
	 * \param[in]		spAssTuner		assemble adapter
	 */
	void FinishTimestep
	(
		const std::vector<IElemDisc<domain_type>*>& vElemDisc,
		ConstSmartPtr<DoFDistribution> dd,
		bool bNonRegularGrid,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner
	);

	template <typename TElem, typename TIterator>
	void
	FinishTimestepElem(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				   ConstSmartPtr<domain_type> spDomain,
				   ConstSmartPtr<DoFDistribution> dd,
				   TIterator iterBegin,
				   TIterator iterEnd,
				   int si, bool bNonRegularGrid,
				   ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
				   ConstSmartPtr<AssemblingTuner<TAlgebra> > spAssTuner);

	template <typename TElem, typename TIterator>
	void
	InitAllExports(const std::vector<IElemDisc<domain_type>*>& vElemDisc,
				   ConstSmartPtr<DoFDistribution> dd,
				   TIterator iterBegin,
				   TIterator iterEnd,
				   int si, bool bNonRegularGrid);

////////////////////////////////////////////////////////////////////////////////
// Error estimators: Not implemented for the ghost-fluid method
////////////////////////////////////////////////////////////////////////////////

public:

	template <typename TElem, typename TIterator>
	static void
	AssembleErrorEstimator
	(
		const std::vector<IElemError<domain_type>*>& vElemDisc,
		ConstSmartPtr<domain_type> spDomain,
		ConstSmartPtr<DoFDistribution> dd,
		TIterator iterBegin,
		TIterator iterEnd,
		int si,
		bool bNonRegularGrid,
		const vector_type& u
	)
	{
		UG_THROW ("AssembleErrorEstimator: No error estimator implemented for the Ghost-Fluid method.");
	}

	template <typename TElem, typename TIterator>
	static void
	AssembleErrorEstimator
	(
		const std::vector<IElemError<domain_type>*>& vElemDisc,
		ConstSmartPtr<domain_type> spDomain,
		ConstSmartPtr<DoFDistribution> dd,
		TIterator iterBegin,
		TIterator iterEnd,
		int si,
		bool bNonRegularGrid,
		std::vector<number> vScaleMass,
		std::vector<number> vScaleStiff,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol
	)
	{
		UG_THROW ("AssembleErrorEstimator: No error estimator implemented for the Ghost-Fluid method.");
	}
	
////////////////////////////////////////////////////////////////////////////////
// Data
////////////////////////////////////////////////////////////////////////////////

private:

	extrapolation_type m_extrapol; ///< the extrapolation at the interface
	
	/**
	 * The following flag is used only for research purposes (measuring the volume
	 * sources/sinks at the embedded interface etc.). It switches off assembing
	 * in all inner elements (i.e. makes the procedures to assemble only in the cut
	 * elements:
	 */
	 bool m_bAssembleOnlyCut;

public:

///	set the level-set function and check it
	void set_LSF
	(
		SmartPtr<ls_grid_func_type> spLSF ///< the function to set
	)
	{
		m_extrapol.set_LSF (spLSF);
	}
	
///	prepares the boundary conditions at the interface: sets all them to Dirichlet-0
	void prepare_interface_bc
	(
		SmartPtr<approx_space_type> spApproxSpace ///< the approximation space of the domain discretization
	)
	{
		m_extrapol.prepare_interface_bc (spApproxSpace);
	}
	
///	adds a Dirichlet BC with a given value on the interface
	void set_Dirichlet_on_if_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		number value ///< the Dirichlet value
	)
	{
		m_extrapol.set_Dirichlet_for (spApproxSpace, fct_name, value);
	}
	
///	adds a Dirichlet BC with a given value on the interface
	void set_Dirichlet_on_if_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		SmartPtr<CplUserData<number, dim> > func ///< the Dirichlet function
	)
	{
		m_extrapol.set_Dirichlet_for (spApproxSpace, fct_name, func);
	}
	
///	adds a "plain" Dirichlet BC with a given value on the interface
	void set_plain_Dirichlet_on_if_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		SmartPtr<CplUserData<number, dim> > func ///< the Dirichlet function
	)
	{
		m_extrapol.set_plain_Dirichlet_for (spApproxSpace, fct_name, func);
	}
	
///	adds a "plain" Dirichlet BC with a given value on the interface
	void set_plain_Dirichlet_on_if_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name, ///< function to impose the condition for
		number value ///< the Dirichlet value
	)
	{
		m_extrapol.set_plain_Dirichlet_for (spApproxSpace, fct_name, value);
	}
	
///	adds a Neumann-0 with on the interface
	void set_Neumann0_on_if_for
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* fct_name ///< function to impose the condition for
	)
	{
		m_extrapol.set_Neumann0_for (spApproxSpace, fct_name);
	}
	
///	excludes a (boundary) subsets from the extrapolation
	void exclude_subsets
	(
		SmartPtr<approx_space_type> spApproxSpace, ///< the approximation space of the domain discretization
		const char* subset_names ///< names of the subsets to exclude
	)
	{
		m_extrapol.exclude_subsets (spApproxSpace, subset_names);
	}
	
///	set the "assemble only in cut elements" flag
	void set_assemble_only_cut (bool b)
	{
		m_bAssembleOnlyCut = b;
	}
	
///	project the level-set function to the coarse levels
	void project_LSF ()
	{
		m_extrapol.project_LSF ();
	}
	
///	returns the extrapolation
	extrapolation_type & extrapolation () {return m_extrapol;}
	
///	checks whether the element is intersected by the interface, or what, and prepares the data
	int check_elem_lsf
	(
		size_t n_co, ///< number of the corners of the element
		GridObject * pElem, ///< the element to process
		int si, ///< subset of the element
		int g_level, ///< grid level of the element
		bool use_hanging, ///< if there can be hanging nodes
		const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
		number time ///< the phisical time
	)
	{
		return m_extrapol.check_elem_lsf (n_co, pElem, si, g_level, use_hanging, vCornerCoords, time);
	}
	
///	extrapolates a component the solution to the vertices behind the interface (w.r.t. a base corner)
	void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		size_t base_co, ///< the base corner
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const
	{
		m_extrapol.extrapolate_by_lsf (num_co, base_co, u, fct);
	}
	
///	extrapolates a component the solution to the vertices behind the interface (by averaging)
	void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const
	{
		m_extrapol.extrapolate_by_lsf (num_co, u, fct);
	}
	
///	returns true if the corner is "inside" (use after check_elem_lsf)
	bool corner_inside
	(
		size_t co ///< the corner
	) const
	{return m_extrapol.corner_inside (co);}
	
///	returns the effective value of the LSF at a corner (use after check_elem_lsf)
	number lsf_at
	(
		size_t co ///< the corner
	) const
	{return m_extrapol.lsf_at (co);}
	
/// sets the values at the outer vertices to 0
	virtual void clear_outer_values
	(
		vector_type & d, ///< the vector where to set
		const DoFDistribution * dd ///< dof distribution of the grid function to reset
	) const
	{m_extrapol.clear_outer_values (d, dd);}
	
/// sets the values at the outer vertices to given values
	virtual void set_outer_values
	(
		vector_type & u, ///< the vector where to set
		const DoFDistribution * dd, ///< dof distribution of the grid function to reset
		number time ///< the physical time
	)
	{m_extrapol.set_outer_values (u, dd, time);}
	
/// sets the matrices at outer vertices to identity
	virtual void set_outer_matrices
	(
		matrix_type & A, ///< the vector where to set
		const DoFDistribution * dd ///< dof distribution of the grid function to reset
	) const
	{m_extrapol.set_outer_matrices (A, dd);}
	
}; // class LSGFGlobAssembler

///	a special constraint that sets functions and matrices in the outer subdomain to given values
/**
 * This class implements a constraint for the LSF based ghost-fluid method
 * to set grid functions and matrices in the outer subdomain to given values.
 *
 * \tparam TDomain			domain type
 * \tparam TAlgebra			algebra type
 * \tparam TExtrapolation	extrapolation class
 */
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
class LSGFConstraint
:	public IDomainConstraint<TDomain, TAlgebra>
{
public:

///	Domain type
	typedef TDomain domain_type;
	
///	Algebra type
	typedef TAlgebra algebra_type;
	
///	Vector type in the algebra
	typedef typename algebra_type::vector_type vector_type;
	
///	Matrix type in the algebra
	typedef typename algebra_type::matrix_type matrix_type;
	
///	Extrapolation type
	typedef TExtrapolation extrapolation_type;
	
///	class constructor
	LSGFConstraint
	(
		extrapolation_type & rExtrapolation ///< the GF extrapolation
	)
	:	m_rExtrapolation (rExtrapolation)
	{}

/// virtual destructor
	virtual ~LSGFConstraint () {};

/// sets a unity row for all conductor indices
	void adjust_jacobian
	(
		matrix_type & J,
		const vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		int type,
		number time = 0.0,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
		const number s_a0 = 1.0
	)
	{
		m_rExtrapolation.set_outer_matrices (J, dd.get ());
	}

/// sets a zero value in the defect for all conductor indices
	void adjust_defect
	(
		vector_type & d,
		const vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		int type,
		number time = 0.0,
		ConstSmartPtr<VectorTimeSeries<vector_type> > vSol = SPNULL,
		const std::vector<number> * vScaleMass = NULL,
		const std::vector<number> * vScaleStiff = NULL
	)
	{
		m_rExtrapolation.clear_outer_values (d, dd.get ());
	}

/// sets the value in the solution for all conductor indices
	void adjust_solution
	(
		vector_type	& u,
		ConstSmartPtr<DoFDistribution> dd,
		int type,
		number time = 0.0
	)
	{
		m_rExtrapolation.set_outer_values (u, dd.get (), time);
	}

///	sets unity rows in A and dirichlet values in right-hand side b
	void adjust_linear
	(
		matrix_type & A,
		vector_type & b,
		ConstSmartPtr<DoFDistribution> dd,
		int type,
		number time = 0.0
	)
	{ // Note that this function is not really used, so it needs not to be optimal.
		m_rExtrapolation.set_outer_matrices (A, dd.get ());
		m_rExtrapolation.clear_outer_values (b, dd.get ());
	}

///	sets the dirichlet value in the right-hand side
	void adjust_rhs
	(
		vector_type & b,
		const vector_type & u,
		ConstSmartPtr<DoFDistribution> dd,
		int type,
		number time = 0.0
	)
	{
		m_rExtrapolation.clear_outer_values (b, dd.get ());
	}

///	returns the type of the constraints
	int type () const {return CT_DIRICHLET;}
	
private:
	
/// Extrapolation in the GF method
	extrapolation_type & m_rExtrapolation;
};

/// domain discretization for the Level-Set Ghost-Fluid method
/**
 * This class template is an implementation of the IDomainDiscretization
 * interface for the Level-Set Ghost-Fluid method.
 *
 * \tparam TDomain          domain type
 * \tparam TAlgebra         algebra type
 */
template <typename TDomain, typename TAlgebra, typename TExtrapolation>
class LSGFDomainDiscretization
:	public DomainDiscretizationBase<TDomain, TAlgebra, LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation> >,
	public IInterfaceExtrapolation<TDomain, TAlgebra>
{
/// Type of the global assembler
	typedef LSGFGlobAssembler<TDomain, TAlgebra, TExtrapolation> gass_type;
	
///	Type of the base class
	typedef DomainDiscretizationBase<TDomain, TAlgebra, gass_type> base_type;

/// Type of the constraint
	typedef LSGFConstraint<TDomain, TAlgebra, TExtrapolation> ls_constraint_type;
	
public:
///	Type of Domain
	typedef TDomain domain_type;

///	Type of algebra
	typedef TAlgebra algebra_type;

///	Type of the grid
	typedef typename TDomain::grid_type grid_type;

///	Type of algebra matrix
	typedef typename algebra_type::matrix_type matrix_type;

///	Type of algebra vector
	typedef typename algebra_type::vector_type vector_type;

///	Type of approximation space
	typedef ApproximationSpace<TDomain>	approx_space_type;
	
///	Type of the LSF grid functions
	typedef typename gass_type::ls_grid_func_type ls_grid_func_type;
	
///	world dimension
	static const int dim = TDomain::dim;
	
public:
///	default Constructor
	LSGFDomainDiscretization (SmartPtr<approx_space_type> pApproxSpace)
	: 	DomainDiscretizationBase<domain_type, algebra_type, gass_type> (pApproxSpace),
		m_spLSFGFConstraint (new ls_constraint_type (this->extrapolation ()))
	{
	//	register the constraint
		this->add (SmartPtr<IDomainConstraint<domain_type, algebra_type> > (m_spLSFGFConstraint));
	//	set the default the boundary condtions at the interface
		gass_type::prepare_interface_bc (pApproxSpace);
	}

/// virtual destructor
	virtual ~LSGFDomainDiscretization () {};
	
///	set the level-set function and check it
	void set_LSF
	(
		SmartPtr<ls_grid_func_type> spLSF ///< the function to set
	)
	{
		gass_type::set_LSF (spLSF);
	}
	
///	sets the Dirichlet boundary condition at the interface for a component of the solution
	void set_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		number value ///< the Dirichlet value
	)
	{
		gass_type::set_Dirichlet_on_if_for (base_type::m_spApproxSpace, fct_name, value);
	}
	
///	sets the Dirichlet boundary condition at the interface for a component of the solution
	void set_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		SmartPtr<CplUserData<number, dim> > func ///< the Dirichlet value (as a function)
	)
	{
		gass_type::set_Dirichlet_on_if_for (base_type::m_spApproxSpace, fct_name, func);
	}
	
///	sets the 'plain' Dirichlet boundary condition at the interface for a component of the solution
	void set_plain_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		number value ///< the Dirichlet value
	)
	{
		gass_type::set_plain_Dirichlet_on_if_for (base_type::m_spApproxSpace, fct_name, value);
	}
	
///	sets the 'plain' Dirichlet boundary condition at the interface for a component of the solution
	void set_plain_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		SmartPtr<CplUserData<number, dim> > func ///< the Dirichlet value (as a function)
	)
	{
		gass_type::set_plain_Dirichlet_on_if_for (base_type::m_spApproxSpace, fct_name, func);
	}
	
#ifdef UG_FOR_LUA
	
///	sets the Dirichlet boundary condition at the interface for a component of the solution
	void set_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		const char* func_name ///< the Dirichlet value (as a name of the LUA function)
	)
	{
		set_Dirichlet_on_if_for (fct_name, LuaUserDataFactory<number,dim>::create(func_name));
	}
	
///	sets the Dirichlet boundary condition at the interface for a component of the solution
	void set_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		LuaFunctionHandle func ///< the Dirichlet value (as a LUA function)
	)
	{
		set_Dirichlet_on_if_for (fct_name, make_sp(new LuaUserData<number,dim>(func)));
	}
	
///	sets the Dirichlet boundary condition at the interface for a component of the solution
	void set_plain_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		const char* func_name ///< the Dirichlet value (as a name of the LUA function)
	)
	{
		set_plain_Dirichlet_on_if_for (fct_name, LuaUserDataFactory<number,dim>::create(func_name));
	}
	
///	sets the Dirichlet boundary condition at the interface for a component of the solution
	void set_plain_Dirichlet_on_if_for
	(
		const char* fct_name, ///< function to impose the condition for
		LuaFunctionHandle func ///< the Dirichlet value (as a LUA function)
	)
	{
		set_plain_Dirichlet_on_if_for (fct_name, make_sp(new LuaUserData<number,dim>(func)));
	}
	
#endif
	
///	sets the Neumann-0 boundary condition at the interface for a component of the solution
	void set_Neumann0_on_if_for
	(
		const char* fct_name ///< function to impose the condition for
	)
	{
		gass_type::set_Neumann0_on_if_for (base_type::m_spApproxSpace, fct_name);
	}
	
///	excludes a (boundary) subsets from the extrapolation
	void exclude_subsets
	(
		const char* subset_names ///< names of the subsets to exclude
	)
	{
		gass_type::exclude_subsets (base_type::m_spApproxSpace, subset_names);
	}
	
///	set the "assemble only in cut elements" flag
	void set_assemble_only_cut (bool b)
	{
		gass_type::set_assemble_only_cut (b);
	}
	
///	project the level-set function to the coarse levels
	void project_LSF ()
	{
		gass_type::project_LSF ();
	}
	
///	checks whether the element is intersected by the interface, or what, and prepares the data
	virtual int check_elem_lsf
	(
		size_t n_co, ///< number of the corners of the element
		GridObject * pElem, ///< the element to process
		int si, ///< subset of the element
		int g_level, ///< grid level of the element
		bool use_hanging, ///< if there can be hanging nodes
		const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
		number time ///< the phisical time
	)
	{
		return gass_type::check_elem_lsf (n_co, pElem, si, g_level, use_hanging, vCornerCoords, time);
	}
	
///	(slow version) checks whether the element is intersected by the interface, or what, and prepares the data
	virtual int check_elem_lsf
	(
		size_t n_co, ///< number of the corners of the element
		GridObject * pElem, ///< the element to process
		int si, ///< subset of the element
		bool use_hanging, ///< if there can be hanging nodes
		const MathVector<dim> vCornerCoords [], ///< coordinates of the corners of the element
		number time ///< the phisical time
	)
	{
		ConstSmartPtr<grid_type> mg = base_type::m_spApproxSpace->domain()->grid ();
		int g_level = mg->get_level (pElem);
		UG_ASSERT (g_level >= 0, "LSGFDomainDiscretization: Grid element without grid level.");
		return gass_type::check_elem_lsf (n_co, pElem, si, g_level, use_hanging, vCornerCoords, time);
	}
	
///	extrapolates a component of the solution to the vertices behind the interface (w.r.t. a base corner)
	virtual void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		size_t base_co, ///< the base corner
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const
	{
		gass_type::extrapolate_by_lsf (num_co, base_co, u, fct);
	}
	
///	extrapolates a component of the solution to the vertices behind the interface (by averaging)
	virtual void extrapolate_by_lsf
	(
		size_t num_co, ///< number of the corners
		number * u, ///< nodal values to extrapolate
		size_t fct ///< index of the function (to identify to type of the extrapolation)
	) const
	{
		gass_type::extrapolate_by_lsf (num_co, u, fct);
	}
	
///	returns true if the corner is "inside" (use after check_elem_lsf)
	virtual bool corner_inside
	(
		size_t co ///< the corner
	) const
	{return gass_type::corner_inside (co);}
	
///	returns the effective value of the LSF at a corner (use after check_elem_lsf)
	virtual number lsf_at
	(
		size_t co ///< the corner
	) const
	{return gass_type::lsf_at (co);}
	
protected:
///	the Level-Set Function constraint
	SmartPtr<ls_constraint_type> m_spLSFGFConstraint;
};

} // end namespace ug

#include "dom_disc_embb_impl.h"

#endif // __H__UG__PLUGINS__D3F__EMBASS__

/* End of File */
