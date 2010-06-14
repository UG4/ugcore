/*
 * grid_function_space.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__
#define __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_discretization/parallelization/parallelization.h"
#endif

namespace ug{

/**
 * The storage type of a grid function is used in parallel applications.
 * We assume that the dofs are distributed to the processes in the way that
 * each dof is master on exactly one process and can be a slave (i.e. a local
 * copy) on several other processes. Given the real values of the dofs the different
 * storage type are defined as follows:
 *  - GFST_UNDEFINED: no information given
 *  - GFST_CONSISTENT: The real value is saved in the master and every slave
 *  - GFST_ADDITIVE: The sum over the values in the master and all slaves gives the exact value
 *  - GFST_UNIQUE: Same as GFST_ADDITIV, but value is zero in all slaves (i.e. master has exact value)
 *
 *  Note, that a GridFunction can have more than one type. E.g. every unique GridFunction is additive. Moreover,
 *  the GridFunction being zero everywhere is consistent, additive and unique at the same time. Therefore,
 *  the information is given bitwise.
 *  To check, whether a function is in a type one may ask:
 *  	 -- v.has_storage_type(GFST_CONSISTENT), etc.
 *  To change into another type
 *       -- v.change_storage_type(GFST_ADDITIVE)
 */
enum GridFunctionStorageType {
	GFST_UNDEFINED = 0,
	GFST_CONSISTENT = 1 << 0,
	GFST_ADDITIVE = 1 << 1,
	GFST_UNIQUE = 1 << 2
};

// A grid function brings approximation space and algebra together. For a given DoFManager and level, a grid function
// represents the solutions on the level 'level'
template <typename TApproximationSpace, typename TDoFManager, typename TAlgebra>
class GridFunction{
	public:
		typedef GridFunction<TApproximationSpace, TDoFManager, TAlgebra> this_type;

		// DOMAIN
		// domain type
		typedef typename TApproximationSpace::domain_type domain_type;

		// global coordinate type
		typedef typename domain_type::position_type position_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::grid_type grid_type;

		// subset handler, where DoF Manager is defined
		typedef typename domain_type::subset_handler_type subset_handler_type;


		// ALGEBRA
		// algebra type
		typedef TAlgebra algebra_type;

		// vector type used to store dof values
		typedef typename algebra_type::vector_type vector_type;

		// local vector type
		typedef typename vector_type::local_vector_type local_vector_type;

		// local index type
		typedef typename vector_type::local_index_type local_index_type;


		// APPROXIMATIONSPACE
		// Approximation Space
		typedef TApproximationSpace approximation_space_type;

		// DOFMANAGER
		// dof manager used for this approximation space
		typedef TDoFManager dof_manager_type;

	public:
		// Default constructor
		GridFunction() :
			m_pApproxSpace(NULL), m_pDoFManager(NULL), m_level(0),
			m_name(""), m_pVector(NULL), m_type(GFST_UNDEFINED)
			{}

		// Constructor
		GridFunction(	std::string name, approximation_space_type& ApproxSpace,
						dof_manager_type& DoFManager, size_t level, bool allocate = true) :
			m_pApproxSpace(&ApproxSpace), m_pDoFManager(&DoFManager), m_level(level),
			m_name(name), m_pVector(NULL), m_type(GFST_UNDEFINED)
		{
			UG_ASSERT(level < m_pDoFManager->num_levels(), "Accessing level that does not exist");
			if(allocate){
				size_t num_dofs = m_pDoFManager->num_dofs(level);
				if(create_storage(num_dofs) != true)
					UG_ASSERT(0, "Cannot create vector memory for " << num_dofs << " dofs.\n");
			}
		};

		// copy constructor
		GridFunction(const this_type& v) :
			m_pApproxSpace(v.m_pApproxSpace), m_pDoFManager(v.m_pDoFManager), m_level(v.m_level),
			m_name(v.m_name), m_pVector(NULL), m_type(GFST_UNDEFINED)
		{
			assign(v);
		};

		// destructor
		~GridFunction()
		{
			release_storage();
		}

		// copies the GridFunction v, except that the values are copied.
		bool clone_pattern(const this_type& v)
		{
			// delete memory if vector storage exists already
			if(m_pVector != NULL) release_storage();

			// copy informations
			m_pApproxSpace = v.m_pApproxSpace;
			m_pDoFManager = v.m_pDoFManager;
			m_level = v.m_level;
			m_name = v.m_name;

			// create new vector
			if(create_storage(v.num_dofs()) != true)
				{UG_LOG("Cannot create pattern.\n"); return false;}

			return true;
		};

		// projects a surface function to this grid function
		// currently this is only implemented for a full refinement
		// (surface level == full refinement level).
		// Then, only the pointer to the dof storage is copied to avoid unnecessary copy-work.
		template <typename TSurfaceDoFManager>
		bool project_surface(GridFunction<TApproximationSpace, TSurfaceDoFManager, TAlgebra>& v)
		{
			if(m_pApproxSpace != v.m_pApproxSpace) {return false;}

			if(m_pDoFManager == v.m_pDoFManager && m_level == v.m_level)
			{
				// set pointer to dof storage
				m_pVector = v.m_pVector;
				m_type = v.m_type;
				return true;
			}
			else{UG_ASSERT(0, "Not implemented.");}
			return false;
		}

		// inverse operation to project surface
		template <typename TSurfaceDoFManager>
		bool release_surface(GridFunction<TApproximationSpace, TSurfaceDoFManager, TAlgebra>& v)
		{
			if(m_pApproxSpace != v.m_pApproxSpace) {return false;}

			if(m_pDoFManager == v.m_pDoFManager && m_level == v.m_level)
			{
				if(m_pVector != v.m_pVector) return false;

				// give information about storage type to surface
				v.m_type = m_type;

				// forget about memory without deleting it.
				m_pVector = NULL;
				set_storage_type(GFST_UNDEFINED);
				return true;
			}
			else
			{UG_ASSERT(0, "Not implemented.");}
			return false;
		}

		// number of dofs
		inline size_t num_dofs() const
			{return m_pDoFManager->num_dofs(m_level);}

		// number of functions
		inline size_t num_fct() const
			{return m_pDoFManager->num_fct();}

		// name of function
		inline std::string get_name(size_t fct) const
			{return m_pDoFManager->get_name(fct);}

		// dim of function
		inline size_t dim_fct(size_t fct) const
			{return m_pDoFManager->dim_fct(fct);}

		// returns if function 'fct' is defined in subset 'si'
		inline bool fct_is_def_in_subset(size_t fct, int si) const
			{return m_pDoFManager->fct_is_def_in_subset(fct, si);}

		// number of subsets
		inline int num_subsets() const
			{return m_pDoFManager->num_subsets();}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator begin(int si) const
			{return m_pDoFManager->template begin<TElem>(si, m_level);}

		// iterator for elements where this grid function is defined
		template <typename TElem>
		inline typename geometry_traits<TElem>::iterator end(int si) const
			{return m_pDoFManager->template end<TElem>(si, m_level);}

		// number of elements of this type for a subset
		template <typename TElem>
		inline size_t num(int si) const
			{return m_pDoFManager->template num<TElem>(si, m_level);}

/* NEEDED ???
		// evaluate Grid function on an element
		template <typename TElem>
		number evaluate(TElem* elem, size_t fct, size_t comp,
						const MathVector<reference_element_traits<TElem>::dim>& loc_pos) const;

		// evaluate Grid function for a global position
		number evaluate(size_t fct, size_t comp, const position_type& glob_pos) const;
*/
		// get dof values
		inline bool get_dof_values(local_vector_type& val, local_index_type& ind) const
			{m_pVector->get(val, ind); return true;}

		// get multiindices on an finite element in canonical order
		template <typename TElem>
		inline size_t get_multi_indices(TElem* elem, size_t fct,
										local_index_type& ind, size_t offset = 0) const
			{return m_pDoFManager->get_multi_indices(elem, fct, ind, offset);}

		// get multiindices on an geometric object in canonical order
		template <typename TGeomObj>
		inline size_t get_multi_indices_of_geom_obj(TGeomObj* elem, size_t fct,
													local_index_type& ind, size_t offset = 0) const
			{return m_pDoFManager->get_multi_indices_of_geom_obj(elem, fct, ind, offset);}

		// get local dof values
		template <typename TElem>
		inline int get_dof_values(TElem* elem, size_t fct, local_vector_type& val) const
		{
			const int num_ind = m_pDoFManager->num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			m_pDoFManager->get_multi_indices(elem, fct, ind);

			m_pVector->get(val, ind);

			return num_ind;
		}

		template <typename TGeomObj>
		inline int get_dof_values_of_geom_obj(TGeomObj* elem, size_t fct, local_vector_type& val) const
		{
			int num_ind = m_pDoFManager->num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: This is a quick hack. Needed: function num_multi_indices_of_geom_obj
			num_ind = m_pDoFManager->get_multi_indices_of_geom_obj(elem, fct, ind);

			ind.resize(num_ind);

			m_pVector->get(val, ind);

			return num_ind;
		}

		template <typename TElem>
		inline bool set_dof_values(TElem* elem, size_t fct, local_vector_type& val)
		{
			const int num_ind = m_pDoFManager->num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			m_pDoFManager->get_multi_indices(elem, fct, ind);

			m_pVector->set(val, ind);

			UG_DEBUG_BEGIN(LIB_DISC_DISCRETE_FUNCTION, 2);
				UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "\n " << num_ind << " Values set: ");
				for(size_t i = 0; i < ind.size(); ++i)
					UG_DLOG(LIB_DISC_DISCRETE_FUNCTION, 2, "(" << ind[i] << ", " << val[i] << ") ");
			UG_DEBUG_END(LIB_DISC_DISCRETE_FUNCTION, 2);

			return true;
		}

		template <typename TGeomObj>
		inline bool set_dof_values_of_geom_obj(TGeomObj* elem, size_t fct, local_vector_type& val)
		{
			int num_ind = m_pDoFManager->num_multi_indices(elem, fct);
			local_index_type ind(num_ind);

			// TODO: This is a quick hack. Needed: function num_multi_indices_of_geom_obj
			num_ind = m_pDoFManager->get_multi_indices_of_geom_obj(elem, fct, ind);

			ind.resize(num_ind);

			m_pVector->set(val, ind);

			return true;
		}

		// approximation space of this gridfunction
		approximation_space_type& get_approximation_space()
			{return *m_pApproxSpace;}

		// returns the dofmanager used for this discrete function
		dof_manager_type& get_dof_manager()
			{return *m_pDoFManager;}

		// type of fundamental function
		LocalShapeFunctionSetID get_local_shape_function_set_id(size_t fct) const
			{return m_pDoFManager->get_local_shape_function_set_id(fct);}

		// returns the level of the dofmanager, that is used for this grid function
		size_t get_level() const
			{return m_level;}

		// export the dof storage of this vector
		vector_type& get_vector()
			{return *m_pVector;}

		// export the dof storage of this vector
		const vector_type& get_vector() const
			{return *m_pVector;}

		// return the domain
		const domain_type& get_domain() const
			{return m_pApproxSpace->get_domain();}

		// return the domain
		domain_type& get_domain()
			{return m_pApproxSpace->get_domain();}

		// this function finalizes the dof pattern.
		// Afterwards the pattern can only be changed by
		// destroying the vector and creating a new one.
		bool finalize()
			{return m_pVector->finalize();}

		// sets grid function
		this_type& operator=(const this_type& v)
			{assign(v); return *this;}

		// add a grid function to this grid function
		this_type& operator+=(const this_type& v)
			{*m_pVector += *(v.m_pVector); return *this;}

		// subtract a grid function from this grid function
		this_type& operator-=(const this_type& v)
			{*m_pVector -= *(v.m_pVector); return *this;}

		// multiply grid function by scalar
		this_type& operator*=(number w)
			{*m_pVector *= w; return *this;}

		// set all dofs on level 'level' to value 'w'
		bool set(number w, GridFunctionStorageType type = GFST_CONSISTENT)
		{
			if(type & GFST_UNDEFINED) return false;

			if(m_pVector->set(w) != true) return false;
			set_storage_type(GFST_CONSISTENT);

			// if w == 0.0, we have all types
			if(w == 0.0){
				add_storage_type(GFST_ADDITIVE);
				add_storage_type(GFST_UNIQUE);
				return true;
			}
			// consistent required
			if(type & GFST_CONSISTENT) return true;

			// additive or additive unique
			return change_storage_type(GFST_UNIQUE);
		}

		// name of grid function
		std::string name()
			{return m_name;}

		// changes to the requested storage type if possible
		bool change_storage_type(GridFunctionStorageType type)
		{
			// can only change if current state is defined
			if(has_storage_type(GFST_UNDEFINED)) return false;

			// if already in that type
			if(has_storage_type(type)) return true;

			// else switch to that type
			switch(type)
			{
			case GFST_CONSISTENT:
					 if(has_storage_type(GFST_UNIQUE)){
						 UniqueToConsistent(	m_pVector,
												m_pDoFManager->get_master_layout(m_level),
												m_pDoFManager->get_slave_layout(m_level));
						 set_storage_type(GFST_CONSISTENT);
						 break;
					 }
					 else if(has_storage_type(GFST_ADDITIVE)){
						AdditiveToConsistent(	m_pVector,
												m_pDoFManager->get_master_layout(m_level),
												m_pDoFManager->get_slave_layout(m_level));
						set_storage_type(GFST_CONSISTENT);
						break;
					}
					else return false;
			case GFST_ADDITIVE:
					if(has_storage_type(GFST_UNIQUE)){
						add_storage_type(GFST_ADDITIVE);
						break;
					}
					else if(has_storage_type(GFST_CONSISTENT)){
						ConsistentToUnique( m_pVector,
											m_pDoFManager->get_slave_layout(m_level));
						set_storage_type(GFST_ADDITIVE);
						add_storage_type(GFST_UNIQUE);
						break;
					}
					else return false;
			case GFST_UNIQUE:
					if(has_storage_type(GFST_ADDITIVE)){
						AdditiveToUnique(	m_pVector,
											m_pDoFManager->get_master_layout(m_level),
											m_pDoFManager->get_slave_layout(m_level));
						add_storage_type(GFST_UNIQUE);
						break;
					}
					else if(has_storage_type(GFST_CONSISTENT)){
						ConsistentToUnique( m_pVector,
											m_pDoFManager->get_slave_layout(m_level));
						set_storage_type(GFST_ADDITIVE);
						add_storage_type(GFST_UNIQUE);
						break;
					}
					else return false;
			default: return false;
			}
			return true;
		}

		// returns if the current storage type has a given representation
		bool has_storage_type(GridFunctionStorageType type)
			{return (bool)(m_type & type);}

		// sets the storage type
		void set_storage_type(GridFunctionStorageType type)
			{m_type = type;}

		// adds the storage type
		void add_storage_type(GridFunctionStorageType type)
			{m_type |= type;}

		// removes the storage type
		void remove_storage_type(GridFunctionStorageType type)
			{m_type &= ~type;}

		// copies the storage type from another vector
		void copy_storage_type(const this_type& v)
			{m_type = v.m_type;}

		// two norm
		inline number two_norm()
		{
			number norm;

			// step 1: make vector d additive unique
#ifdef UG_PARALLEL
			if(!change_storage_type(GFST_UNIQUE)) {
				UG_ASSERT(0, "Cannot change to unique representation.");
				return -1;
			}
#endif

			// step 2: compute new defect norm
			double tNormLocal = (double)local_two_norm();
#ifdef UG_PARALLEL
			tNormLocal *= tNormLocal;
			double tNormGlobal;
			pcl::AllReduce(&tNormLocal, &tNormGlobal, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
			norm = sqrt((number)tNormGlobal);
#else
			norm = tNormLocal;
#endif
			return norm;
		}

		// square norm for of the part of grid function on this process
		inline number local_two_norm() const
			{return m_pVector->two_norm();}

	protected:
		// sets the values of GridFunction 'v' to this GridFunction
		// DofManager and level must be the same
		bool assign(const this_type& v)
		{
			// check that Grid functions are of same type
			if(	m_pApproxSpace != v.m_pApproxSpace ||
				m_pDoFManager != v.m_pDoFManager ||
				m_level != v.m_level)
					return false;

			// allocate memory if needed
			if(m_pVector == NULL) create_storage(v.m_pVector->size());
			else if (v.m_pVector->size() != m_pVector->size())
				{UG_LOG("Size of discrete function does not match."); return false;}

			// copy values
			*m_pVector = *v.m_pVector;
			m_type = v.m_type;
			return true;
		}

		// creates storage for dofs
		bool create_storage(size_t num_dofs)
		{
			m_pVector = new vector_type;
			if(m_pVector->create(num_dofs) != true) return false;
			else return true;
		}

		// deletes the memory
		bool release_storage()
		{
			if(m_pVector != NULL){m_pVector->destroy(); delete m_pVector;}
			return true;
		}

	protected:
		// Approximation Space
		approximation_space_type* m_pApproxSpace;

		// dof manager of this discrete function
		dof_manager_type* m_pDoFManager;

		// level of Approximation Space, where this Grid Function lives
		size_t m_level;

		// name
		std::string m_name;

		// vector storage, to store values of local degrees of freedom
		vector_type* m_pVector;

		// type of storage in parallel (i.e. consistent, additiv, additiv unique)
		int m_type;
};

template <typename TApproximationSpace, typename TDoFManager, typename TAlgebra>
inline std::ostream& operator<< (std::ostream& outStream, const GridFunction<TApproximationSpace, TDoFManager, TAlgebra>& v)
{
	outStream << v.get_vector();
	return outStream;
}


// Surface Grid function lives on the surface grid
template <typename TApproximationSpace, typename TAlgebra>
class SurfaceGridFunction : public GridFunction<TApproximationSpace, typename TApproximationSpace::surface_dof_manager_type, TAlgebra>{
	public:
		// Approximation Space
		typedef TApproximationSpace approximation_space_type;

		// DOFMANAGER
		// dof manager used for this approximation space
		typedef typename TApproximationSpace::surface_dof_manager_type dof_manager_type;

	public:
		// Constructor
		SurfaceGridFunction(approximation_space_type& approximationSpace, dof_manager_type& dof_manager) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(approximationSpace, dof_manager, dof_manager.num_levels() - 1)
		{};

		SurfaceGridFunction(const SurfaceGridFunction& v) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(v)
		{}
};

// Level Grid function lives on a level grid
template <typename TApproximationSpace, typename TAlgebra>
class LevelGridFunction : public GridFunction<TApproximationSpace, typename TApproximationSpace::level_dof_manager_type, TAlgebra>{
public:
	// Approximation Space
	typedef TApproximationSpace approximation_space_type;

	// DOFMANAGER
	// dof manager used for this approximation space
	typedef typename TApproximationSpace::level_dof_manager_type dof_manager_type;

	public:
		// Constructor
		LevelGridFunction(approximation_space_type& approximationSpace, dof_manager_type& dof_manager, size_t level) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(approximationSpace, dof_manager, level)
		{};

		LevelGridFunction(const LevelGridFunction& v) :
			GridFunction<TApproximationSpace, dof_manager_type, TAlgebra>(v)
		{}


};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__FUNCTION_SPACE__GRID_FUNCTION__ */
