/*
 * parallel_grid_function.h
 *
 *  Created on: 19.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_GRID_FUNCTION__
#define __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_GRID_FUNCTION__

#include "pcl/pcl.h"
#include "lib_algebra/parallelization/parallelization.h"
#include "lib_discretization/parallelization/parallelization.h"

namespace ug
{

template <typename TGridFunction>
class ParallelGridFunction : public TGridFunction
{
	public:
		typedef ParallelGridFunction<TGridFunction> this_type;

		// approximation space type
		typedef typename TGridFunction::approximation_space_type approximation_space_type;

		// Domain
		typedef typename TGridFunction::domain_type domain_type;

		// algebra type
		typedef typename TGridFunction::algebra_type algebra_type;

		// vector type used to store dof values
		typedef typename TGridFunction::vector_type vector_type;

		// local vector type
		typedef typename TGridFunction::local_vector_type local_vector_type;

		// local index type
		typedef typename TGridFunction::local_index_type local_index_type;

		// dof distribution
		typedef typename TGridFunction::dof_distribution_type dof_distribution_type;

	public:
		// Constructor
		ParallelGridFunction() : TGridFunction() {};

		// Constructor
		ParallelGridFunction(std::string name, approximation_space_type& approxSpace, const dof_distribution_type& DoFDistr, bool allocate = true)
			: TGridFunction(name, approxSpace, DoFDistr, allocate)
		{
			if(allocate)
			{
				set_layouts();
				set_storage_type(PST_UNDEFINED);
			}
		};


		///////////////////////////////
		// overloaded functions
		///////////////////////////////
		void assign_dof_distribution(const dof_distribution_type& DoFDistr, bool allocate = true)
		{
			TGridFunction::assign_dof_distribution(DoFDistr, allocate);
			set_layouts();
			set_storage_type(PST_UNDEFINED);
		}

		virtual bool clone_pattern(const this_type& v)
		{
			if(!TGridFunction::clone_pattern(v)) return false;
			set_layouts();
			set_storage_type(PST_UNDEFINED);
			return true;
		}

		// clone
		this_type& clone()
		{return *(new this_type(*this));}

		//	set a vector
		bool assign(const vector_type& v)
		{
			if(v.size() != TGridFunction::get_vector().size()) return false;
			TGridFunction::get_vector() = v;
			set_layouts();
			TGridFunction::get_vector().copy_storage_type(v);
			return true;
		}

		bool assign(const this_type& v)
		{
			if(!TGridFunction::assign(v))
			{
				UG_ASSERT(0, "Assigning failed.");
				return false;
			}
			set_layouts();
			copy_storage_type(v);
			return true;
		}

		// sets grid function
		this_type& operator=(const this_type& v)
			{assign(v); return *this;}

		// for Template Expressions
		/*template<typename T>
		this_type& operator =(const T &t)
		{
			VectorAssign(*this, t);
			return *this;
		}
		template<typename T>
		this_type& operator -=(const T &t)
		{
			VectorSub(*this, t);
			return *this;
		}
		template<typename T>
		this_type& operator +=(const T &t)
		{
			VectorAdd(*this, t);
			return *this;
		}*/

		// set all dofs on level 'level' to value 'w'
		// removed this because it is the same in ParallelVector and Vector (MR)
		/*bool set(number w, ParallelStorageType type)
			{return TGridFunction::get_vector().set(w, type);}

		bool set(number w)
			{return TGridFunction::get_vector().set(w, PST_CONSISTENT);}*/

		////////////////////////////
		// Storage type
		////////////////////////////

		bool change_storage_type_by_string(std::string type)
		{
			if(type == "consistent")
				return change_storage_type(PST_CONSISTENT);
			else if(type == "additive")
				return change_storage_type(PST_ADDITIVE);
			else if(type == "unique")
				return change_storage_type(PST_UNIQUE);
			else
			{
				UG_LOG("Could not change storage type by name.\n");
			}
			return false;
		}

		void set_storage_type_by_string(std::string type)
		{
			if(type == "consistent")
				{set_storage_type(PST_CONSISTENT); return;}
			else if(type == "additive")
				{set_storage_type(PST_ADDITIVE); return;}
			else if(type == "unique")
				{set_storage_type(PST_UNIQUE); return;}
			else
			{
				UG_LOG("Could not set storage type by name.\n");
			}
			return;
		}

		// changes to the requested storage type if possible
		bool change_storage_type(ParallelStorageType type)
			{return TGridFunction::get_vector().change_storage_type(type);}

		// returns if the current storage type has a given representation
		bool has_storage_type(ParallelStorageType type) const
			{return TGridFunction::get_vector().has_storage_type(type);}

		// returns storage type mask
		ParallelStorageType get_storage_mask() const { return TGridFunction::get_vector().get_storage_mask(); }

		// sets the storage type
		void set_storage_type(ParallelStorageType type)
			{TGridFunction::get_vector().set_storage_type(type);}

		// adds the storage type
		void add_storage_type(ParallelStorageType type)
			{TGridFunction::get_vector().add_storage_type(type);}

		// removes the storage type
		void remove_storage_type(ParallelStorageType type)
			{TGridFunction::get_vector().remove_storage_type(type);}

		// copies the storage type from another vector
		void copy_storage_type(const this_type& v)
			{TGridFunction::get_vector().copy_storage_type(*v.get_vector());}

		///////////////////////////////
		// index layouts
		///////////////////////////////

		inline IndexLayout& get_slave_layout()	{return TGridFunction::m_pNonConstDoFDistribution->get_slave_layout();}
		inline IndexLayout& get_master_layout()	{return TGridFunction::m_pNonConstDoFDistribution->get_master_layout();}
		inline IndexLayout& get_vertical_slave_layout()		{return TGridFunction::m_pNonConstDoFDistribution->get_vertical_slave_layout();}
		inline IndexLayout& get_vertical_master_layout()	{return TGridFunction::m_pNonConstDoFDistribution->get_vertical_master_layout();}

		inline pcl::ParallelCommunicator<IndexLayout>& get_communicator() {return TGridFunction::m_pNonConstDoFDistribution->get_communicator();;}
		inline pcl::ProcessCommunicator& get_process_communicator()	{return TGridFunction::m_pNonConstDoFDistribution->get_process_communicator();}

		///////////////////////////////
		// help functions
		///////////////////////////////

	protected:
		void set_layouts()
		{
			TGridFunction::get_vector().set_slave_layout(get_slave_layout());
			TGridFunction::get_vector().set_master_layout(get_master_layout());
			TGridFunction::get_vector().set_vertical_slave_layout(get_vertical_slave_layout());
			TGridFunction::get_vector().set_vertical_master_layout(get_vertical_master_layout());

			TGridFunction::get_vector().set_communicator(get_communicator());
			TGridFunction::get_vector().set_process_communicator(get_process_communicator());
		}


};

// for template expressions
// d = alpha*v1
template<typename T>
inline void VecScale(ParallelGridFunction<T> &dest,
		double alpha1, const ParallelGridFunction<T> &v1)
{
	dest.copy_storage_type(v1);
	VecScale(dest.get_vector(), alpha1, v1.get_vector());
}


// dest = alpha1*v1 + alpha2*v2
template<typename T>
inline void VecScaleAdd(ParallelGridFunction<T>  &dest,
		double alpha1, const ParallelGridFunction<T> &v1,
		double alpha2, const ParallelGridFunction<T> &v2)
{
	ParallelStorageType mask = v1.get_storage_mask() & v2.get_storage_mask();
	UG_ASSERT(mask != 0, "VecScaleAdd: cannot add vectors v1 and v2");
	dest.set_storage_type(mask);

	VecScaleAdd(dest.get_vector(), alpha1, v1.get_vector(), alpha2, v2.get_vector());
}

// dest = alpha1*v1 + alpha2*v2 + alpha3*v3
template<typename T>
inline void VecScaleAdd(ParallelGridFunction<T> &dest,
		double alpha1, const ParallelGridFunction<T> &v1,
		double alpha2, const ParallelGridFunction<T> &v2,
		double alpha3, const ParallelGridFunction<T> &v3)
{
	ParallelStorageType mask = v1.get_storage_mask() & v2.get_storage_mask() & v3.get_storage_mask();
	UG_ASSERT(mask != 0, "VecScaleAdd: cannot add vectors v1 and v2");
	dest.set_storage_type(mask);

	VecScaleAdd(dest.get_vector(), alpha1, v1.get_vector(), alpha2, v2.get_vector(), alpha3, v3.get_vector());
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_GRID_FUNCTION__ */
