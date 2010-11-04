/*
 * grid_function_space.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__

#include "lib_algebra/operator/operator_base_interface.h"

namespace ug{

// predeclaration
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
class ApproximationSpace;


// A grid function brings approximation space and algebra together. For a given DoFManager and level, a grid function
// represents the solutions on the level 'level'
template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
class GridFunction :	public TAlgebra::vector_type,
						public virtual IFunctionBase
{
	public:
	// 	This type
		typedef GridFunction<TDomain, TDoFDistribution, TAlgebra> this_type;

	// 	Type of Approximation space
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra> approximation_space_type;

	// 	Domain
		typedef TDomain domain_type;

	// 	Algebra type
		typedef TAlgebra algebra_type;

	// 	Vector type used to store dof values
		typedef typename algebra_type::vector_type vector_type;

	// 	Local vector type
		typedef LocalVector<typename vector_type::value_type> local_vector_type;

	// 	Local index type
		typedef LocalIndices local_index_type;

	//	Multi index
		typedef typename TDoFDistribution::multi_index_vector_type multi_index_vector_type;

	// 	Algebra index
		typedef typename TDoFDistribution::algebra_index_vector_type algebra_index_vector_type;

		// DOF DISTRIBUTION
		// dof manager used for this approximation space
		typedef TDoFDistribution dof_distribution_type;

	public:
		// Default constructor
		GridFunction() :
			m_name(""), m_pApproxSpace(NULL), m_pDoFDistribution(NULL), m_pNonConstDoFDistribution(NULL)
			{}

		// Constructor
		GridFunction(	std::string name, approximation_space_type& approxSpace,
						const dof_distribution_type& DoFDistr, bool allocate = true) :
			m_name(name), m_pApproxSpace(&approxSpace)
		{
			assign_dof_distribution(DoFDistr);

			if(allocate)
				if(!create_storage())
					UG_ASSERT(0, "Cannot create vector memory.\n");
		};

		// copy constructor
		GridFunction(const this_type& v) :
			m_name(v.m_name), m_pApproxSpace(v.m_pApproxSpace),
			m_pDoFDistribution(v.m_pDoFDistribution)
		{
			assign(v);
		};

		// sets grid function
		this_type& operator=(const this_type& v)
			{assign(v); return *this;}

		//	set a vector
		bool assign(const vector_type& v)
		{
			if(v.size() != get_vector().size()) return false;
			get_vector() = v;
			return true;
		}

		// destructor
		virtual ~GridFunction()
		{
			release_storage();
		}

		// clone
		this_type& clone()
		{
			return *(new this_type(*this));
		}

		// copies the GridFunction v, except that the values are copied.
		virtual bool clone_pattern(const this_type& v)
		{
		// 	delete memory if vector storage exists already
			release_storage();

		// 	copy informations
			m_name = v.m_name;
			m_pApproxSpace = v.m_pApproxSpace;
			m_pDoFDistribution = v.m_pDoFDistribution;

		// 	create new vector
			if(!create_storage())
				{UG_LOG("Cannot create pattern.\n"); return false;}

			return true;
		};

	/////////////////////////////////
	// help functions
	/////////////////////////////////
	protected:
		// create storage vector. DoFDistrution must be set.
		bool create_storage()
		{
		//	DoFDistribution is mandatory
			if(m_pDoFDistribution == NULL)
				{UG_LOG("Cannot create vector without DoFDistribution.\n"); return false;}

		//	read number of algebra dofs
			size_t num_dofs = m_pDoFDistribution->num_dofs();

		//	check if vector has correct size. Otherwise resize
			if(get_vector().size() != num_dofs)
			{
				if(!get_vector().destroy()) return false;
				if(!get_vector().create(num_dofs)) return false;
			}
			return true;
		}

		// deletes the memory
		bool release_storage()
		{
		//	delete vector
			get_vector().destroy();
			return true;
		}

		// sets the values of GridFunction 'v' to this GridFunction
		// DofManager and level must be the same
		bool assign(const this_type& v)
		{
			// check that approximation space is equal
			if(m_pApproxSpace != v.m_pApproxSpace)
				return false;

			// check that Grid functions are of same type
			if(	m_pDoFDistribution != v.m_pDoFDistribution)
				return false;

			if(!create_storage())
			{UG_LOG("Cannot create storage for assignment.\n"); return false;}

			if (v.get_vector().size() != get_vector().size())
				{UG_LOG("Size of discrete function does not match."); return false;}

			// copy values
			get_vector() = v.get_vector();
			return true;
		}

		/////////////////////////////////
		// General informations
		/////////////////////////////////
	public:
		// name of grid function
		std::string name()
			{return m_name;}

		// dimension
		int get_dim() const {return domain_type::dim;}

		void assign_dof_distribution(const dof_distribution_type& DoFDistr, bool allocate = true)
		{
			m_pDoFDistribution = &DoFDistr;
			m_pNonConstDoFDistribution = const_cast<dof_distribution_type*>(&DoFDistr);
			if(allocate)
				create_storage();
			else
				release_storage();
		}

		void assign_approximation_space(approximation_space_type& ApproxSpace)
		{
			m_pApproxSpace = &ApproxSpace;
		}

		void assign_surface_approximation_space(approximation_space_type& ApproxSpace)
		{
		//	if already assigned, do nothing
			if(m_pApproxSpace != NULL) return;

		//	assign
			assign_approximation_space(ApproxSpace);
			assign_dof_distribution(ApproxSpace.get_surface_dof_distribution(), true);
		}

		// get dof distribution
		const dof_distribution_type& get_dof_distribution() const
			{return *m_pDoFDistribution;}

		// get approximation space
		const approximation_space_type& get_approximation_space() const {return *m_pApproxSpace;}
		approximation_space_type& get_approximation_space() {return *m_pApproxSpace;}

		/////////////////////////////////
		// DoF Distribution requirements
		/////////////////////////////////

		/// number of discrete functions
		size_t num_fct() const {return m_pDoFDistribution->num_fct();}

		/// number of discrete functions on subset si
		size_t num_fct(int si) const {return m_pDoFDistribution->num_fct(si);}

		/// returns the trial space of the discrete function fct
		LocalShapeFunctionSetID local_shape_function_set_id(size_t fct) const  {return m_pDoFDistribution->local_shape_function_set_id(fct);}

		/// returns the name of the discrete function nr_fct
		std::string name(size_t fct) const {return m_pDoFDistribution->name(fct);}

		/// returns the dimension in which solution lives
		int dim(size_t fct) const {return m_pDoFDistribution->dim(fct);}

		/// returns true if the discrete function nr_fct is defined on subset s
		bool is_def_in_subset(size_t fct, int si) const {return m_pDoFDistribution->is_def_in_subset(fct, si);}


		/// number of subsets
		inline int num_subsets() const {return m_pDoFDistribution->num_subsets();}

		/// return the number of dofs distributed
		inline size_t num_dofs() const {return m_pDoFDistribution->num_dofs();}

		/// return the number of dofs distributed on subset si
		inline size_t num_dofs(int si) const {return m_pDoFDistribution->num_dofs(si);}

		// number of elements of this type for a subset
		template <typename TElem>
		inline size_t num() const
			{return m_pDoFDistribution->template num<TElem>();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::const_iterator begin() const
			{return m_pDoFDistribution->template begin<TElem>();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::const_iterator end() const
			{return m_pDoFDistribution->template end<TElem>();}

		// number of elements of this type for a subset
		template <typename TElem>
		inline size_t num(int si) const
			{return m_pDoFDistribution->template num<TElem>(si);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::const_iterator begin(int si) const
			{return m_pDoFDistribution->template begin<TElem>(si);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::const_iterator end(int si) const
			{return m_pDoFDistribution->template end<TElem>(si);}

		////////// Local Algebra ////////////

		/// number of algebra indices on an element
		size_t num_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
			{return m_pDoFDistribution->num_indices(refID, si, funcGroup);}

		/// number of algebra indices on an element
		size_t num_inner_indices(ReferenceObjectID refID, int si, const FunctionGroup& funcGroup) const
			{return m_pDoFDistribution->num_inner_indices(refID, si, funcGroup);}

		/// fill local informations in LocalIndex
		bool prepare_indices(ReferenceObjectID refID, int si, LocalIndices& ind, bool useHanging = false) const
			{return m_pDoFDistribution->prepare_indices(refID, si, ind, useHanging);}

		/// fill local informations in LocalIndex
		bool prepare_inner_indices(ReferenceObjectID refID, int si, LocalIndices& ind) const
			{return m_pDoFDistribution->prepare_inner_indices(refID, si, ind);}

		/// fill the global algebra indices in LocalIndex
		template<typename TElem>
		void update_indices(TElem* elem, LocalIndices& ind, bool useHanging = false) const
			{return m_pDoFDistribution->update_indices(elem, ind, useHanging);}

		/// fill the global algebra indices in LocalIndex
		template<typename TElem>
		void update_inner_indices(TElem* elem, LocalIndices& ind) const
			{return m_pDoFDistribution->update_inner_indices(elem, ind);}

		////////// Multi indices ////////////

		// number of multi indices on an finite element in canonical order
		template <typename TElem>
		inline size_t num_multi_indices(TElem* elem, size_t fct) const
			{return m_pDoFDistribution->num_multi_indices(elem, fct);}

		// number of multi indices on an geometric object in canonical order
		template <typename TGeomObj>
		inline size_t num_inner_multi_indices(TGeomObj* elem, size_t fct) const
			{return m_pDoFDistribution->num_inner_multi_indices(elem, fct);}

		// get multi indices on an finite element in canonical order
		template <typename TElem>
		inline size_t get_multi_indices(TElem* elem, size_t fct, multi_index_vector_type& ind) const
			{return m_pDoFDistribution->get_multi_indices(elem, fct, ind);}

		// get multi indices on an geometric object in canonical order
		template <typename TGeomObj>
		inline size_t get_inner_multi_indices(TGeomObj* elem, size_t fct,	multi_index_vector_type& ind) const
			{return m_pDoFDistribution->get_inner_multi_indices(elem, fct, ind);}

		////////// Algebra indices ////////////

		// number of algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		inline size_t num_algebra_indices(TGeomObj* elem, size_t fct) const
			{return m_pDoFDistribution->num_algebra_indices(elem, fct);}

		// number of algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		inline size_t num_inner_algebra_indices(TGeomObj* elem, size_t fct) const
			{return m_pDoFDistribution->num_inner_algebra_indices(elem, fct);}

		// get algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		inline void get_algebra_indices(TGeomObj* elem, algebra_index_vector_type& ind) const
			{m_pDoFDistribution->get_algebra_indices(elem, ind);}

		// get algebra indices on an geometric object in canonical order
		template <typename TGeomObj>
		inline void get_inner_algebra_indices(TGeomObj* elem, algebra_index_vector_type& ind) const
			{m_pDoFDistribution->get_inner_algebra_indices(elem, ind);}

		////////// DoF Values ////////////

		// get dof values
		inline number get_dof_value(size_t i, size_t comp) const
			{return BlockRef(((get_vector())[i]), comp);}

		////////////////////////////
		// Algebra requirements
		////////////////////////////
		// TODO: Since a GridFunction is a Vector now, we can replace get_vector everywhere
		// 			and use GridFunction directly
		// export the dof storage of this vector
		vector_type& get_vector()
			{return *(dynamic_cast<vector_type*>(this));}

		// export the dof storage of this vector
		const vector_type& get_vector() const
			{return *(dynamic_cast<const vector_type*>(this));}

	protected:
		// name
		std::string m_name;

		// Approximation Space
		approximation_space_type* m_pApproxSpace;

		// dof manager of this discrete function
		const dof_distribution_type* m_pDoFDistribution;

		// Todo: remove this
		dof_distribution_type* m_pNonConstDoFDistribution;
};

template <typename TDomain, typename TDoFDistribution, typename TAlgebra>
inline std::ostream& operator<< (std::ostream& outStream, const GridFunction<TDomain, TDoFDistribution, TAlgebra>& v)
{
	outStream << v.get_vector();
	return outStream;
}

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__ */
