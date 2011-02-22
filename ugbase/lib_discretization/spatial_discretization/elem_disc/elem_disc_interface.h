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
#include "lib_discretization/assemble.h"
#include "lib_discretization/common/local_algebra.h"

namespace ug{

enum IElemDiscNeed {
	IEDN_NONE = 0,
	IEDN_DEFECT = 1 << 0,
	IEDN_JACOBIAN = 1 << 1,
	IEDN_LINEAR = 1 << 2,
	IEDN_STIFFNESS = 1 << 3,
	IEDN_MASS = 1 << 4
};

// predeclaration
template <typename TAlgebra> class IDataExport;
template <typename TAlgebra> class IDataImport;

/**
 * Element Discretizations
 *
 * \defgroup lib_disc_elem_disc Elem Disc
 * \ingroup lib_discretization
 */

/// \ingroup lib_disc_elem_disc
/// @{


template <typename TAlgebra>
class IElemDisc{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Local matrix type
		typedef LocalMatrix<typename TAlgebra::matrix_type::value_type>
				local_matrix_type;

	/// Local vector type
		typedef LocalVector<typename TAlgebra::vector_type::value_type>
				local_vector_type;

	/// Local index type
		typedef LocalIndices local_index_type;

	public:
	///	Constructor
		IElemDisc() : m_pPattern(NULL), m_id(-1) {}

	////////////////////////////
	// Functions and Subsets
	////////////////////////////

	///	set underlying function pattern
		void set_pattern(const FunctionPattern& pattern)
		{
			m_pPattern = &pattern;
			m_FunctionGroup.set_function_pattern(pattern);
			m_SubsetGroup.set_subset_handler(*pattern.get_subset_handler());
		}

	///	sets functions by name list, divided by ','
		virtual bool set_functions(const char* functions);

	///	sets functions using function group
		virtual bool set_functions(const FunctionGroup& funcGroup);

	///	sets subset(s) by name list, divided by ','
		virtual bool set_subsets(const char* subsets);

	///	sets subset(s) by subset group
		virtual bool set_subsets(const SubsetGroup& subsetGroup);

	///	returns the funtion group
		const FunctionGroup& get_function_group() const {return m_FunctionGroup;}

	///	returns the subset group
		const SubsetGroup& get_subset_group() const {return m_SubsetGroup;}

	protected:
	///	underlying function pattern
		const FunctionPattern* m_pPattern;

	/// function group
		FunctionGroup m_FunctionGroup;

	/// subset group
		SubsetGroup m_SubsetGroup;

	public:
	////////////////////////////
	// IP Data and Coupling
	////////////////////////////

	///	registers a data import
		void register_import(IDataImport<TAlgebra>& Imp)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vIImport.size(); ++i)
				if(m_vIImport[i] == &Imp)
					throw(UGFatalError("Trying to register import twice."));

		//	set function group
			Imp.set_function_group(m_FunctionGroup);

		//	add it
			m_vIImport.push_back(&Imp);
		}

	///	registers a data export
		void register_export(IDataExport<TAlgebra>& Exp)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vIExport.size(); ++i)
				if(m_vIExport[i] == &Exp)
					throw(UGFatalError("Trying to register export twice."));

		//	set function group
			Exp.set_function_group(m_FunctionGroup);

		//	add it
			m_vIExport.push_back(&Exp);
		}

	///	returns number of imports
		size_t num_imports() const {return m_vIImport.size();}

	/// returns an import
		IDataImport<TAlgebra>* import(size_t i)
		{
			UG_ASSERT(i < num_imports(), "Invalid index");
			return m_vIImport[i];
		}

	///	returns number of imports
		size_t num_exports() const {return m_vIExport.size();}

	/// returns an import
		IDataExport<TAlgebra>* get_export(size_t i)
		{
			UG_ASSERT(i < num_exports(), "Invalid index");
			return m_vIExport[i];
		}

protected:
	///	sets the function group in all registered import/exports
		void set_function_group_for_ipdata(const FunctionGroup& fctGrp)
		{
			for(size_t i = 0; i < m_vIImport.size(); ++i)
				m_vIImport[i]->set_function_group(fctGrp);

			for(size_t i = 0; i < m_vIExport.size(); ++i)
				m_vIExport[i]->set_function_group(fctGrp);
		}

	/// data imports
		std::vector<IDataImport<TAlgebra>*> m_vIImport;

	///	data exports
		std::vector<IDataExport<TAlgebra>*> m_vIExport;

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
		bool set_geometric_object_type(int id, IElemDiscNeed need);

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
		bool assemble_JA(local_matrix_type& J, const local_vector_type& u, number time=0.0)
			{return (this->*(m_vAssJAFunc[m_id]))(J, u, time);}

	/// Assembling of Jacobian (Mass part)
	/**
	 * This function assembles the local (mass) jacobian for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_JM(local_matrix_type& J, const local_vector_type& u, number time=0.0)
			{return (this->*(m_vAssJMFunc[m_id]))(J, u, time);}

	/// Assembling of Defect (Stiffness part)
	/**
	 * This function assembles the local (stiffness) defect for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_A(local_vector_type& d, const local_vector_type& u, number time=0.0)
			{return (this->*(m_vAssAFunc[m_id]))(d, u, time);}


	/// Assembling of Defect (Mass part)
	/**
	 * This function assembles the local (mass) defect for the current
	 * solution u and the time (iff timedependent).
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_M(local_vector_type& d, const local_vector_type& u, number time=0.0)
			{return (this->*(m_vAssMFunc[m_id]))(d, u, time);}

	/// Assembling of Right-Hand Side
	/**
	 * This function assembles the local rhs.
	 * <b>NOTE:</b>Before this method can be used, the method
	 * 'set_geometric_object_type must have been called to set the elem type.
	 */
		bool assemble_f(local_vector_type& rhs, number time=0.0)
			{return (this->*(m_vAssFFunc[m_id]))(rhs, time);}

	/// Virtual destructor
		virtual ~IElemDisc() {}

	protected:
		// register the functions
		template <typename TAssFunc> void register_prepare_element_loop_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_prepare_element_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_finish_element_loop_function(int id, TAssFunc func);

		template <typename TAssFunc> void register_assemble_JA_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_JM_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_A_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_M_function(int id, TAssFunc func);
		template <typename TAssFunc> void register_assemble_f_function(int id, TAssFunc func);

	protected:
		// checks if the needed functions are registered for the id type
		bool function_registered(int id, IElemDiscNeed need);

		// checks if the functions are present
		bool prepare_element_loop_function_registered(int id);
		bool prepare_element_function_registered(int id);
		bool finish_element_loop_function_registered(int id);

		// checks if the functions are present
		bool assemble_JA_function_registered(int id);
		bool assemble_JM_function_registered(int id);
		bool assemble_A_function_registered(int id);
		bool assemble_M_function_registered(int id);
		bool assemble_f_function_registered(int id);

	private:
		// abbreviation for own type
		typedef IElemDisc<TAlgebra> T;

		// types of loop function pointers
		typedef bool (T::*PrepareElementLoopFunc)();
		typedef bool (T::*PrepareElementFunc)(	GeometricObject* obj,
												const local_vector_type& u,
												const local_index_type& glob_ind);
		typedef bool (T::*FinishElementLoopFunc)();

		// types of Jacobian assemble functions
		typedef bool (T::*AssembleJAFunc)(	local_matrix_type& J,
											const local_vector_type& u, number time);
		typedef bool (T::*AssembleJMFunc)(	local_matrix_type& J,
											const local_vector_type& u, number time);

		// types of Defect assemble functions
		typedef bool (T::*AssembleAFunc)( 	local_vector_type& d,
											const local_vector_type& u, number time);
		typedef bool (T::*AssembleMFunc)(	local_vector_type& d,
											const local_vector_type& u, number time);

		// types of right hand side assemble functions
		typedef bool (T::*AssembleFFunc)(local_vector_type& d, number time);

	private:
		// loop function pointers
		std::vector<PrepareElementLoopFunc> m_vPrepareElemLoopFunc;
		std::vector<PrepareElementFunc> 	m_vPrepareElemFunc;
		std::vector<FinishElementLoopFunc> 	m_vFinishElemLoopFunc;

		// Jacobian function pointers
		std::vector<AssembleJAFunc> 	m_vAssJAFunc;
		std::vector<AssembleJMFunc> 	m_vAssJMFunc;

		// Defect function pointers
		std::vector<AssembleAFunc> 	m_vAssAFunc;
		std::vector<AssembleMFunc> 	m_vAssMFunc;

		// Rhs function pointers
		std::vector<AssembleFFunc> 	m_vAssFFunc;

	protected:
		// current Geometric Object
		int m_id;
};

/// @}

} // end namespace ug

#include "elem_disc_interface_impl.h"

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__ELEM_DISC__ELEM_DISC_INTERFACE__ */
