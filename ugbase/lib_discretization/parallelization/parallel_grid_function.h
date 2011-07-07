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
	/// Default Constructor
		ParallelGridFunction() : TGridFunction() {};

	/// Initializing Constructor
		ParallelGridFunction(approximation_space_type& approxSpace, dof_distribution_type& DoFDistr)
			: TGridFunction(approxSpace, DoFDistr)
		{
			copy_layouts_into_vector();
			set_storage_type(PST_UNDEFINED);
		};


		///////////////////////////////
		// overloaded functions
		///////////////////////////////

	///	assign dof distribution
		void assign_dof_distribution(dof_distribution_type& DoFDistr)
		{
		//	assign distribution
			TGridFunction::assign_dof_distribution(DoFDistr);

		//	copy layouts into vector
			copy_layouts_into_vector();

		//	set storage type to undefined
			set_storage_type(PST_UNDEFINED);
		}

		virtual void clone_pattern(const this_type& v)
		{
		//	clone process local grid function
			TGridFunction::clone_pattern(v);

		//	copy layouts
			copy_layouts_into_vector();

		//	set storage type to undefined
			set_storage_type(PST_UNDEFINED);
		}

	/// clone
		this_type& clone() {return *(new this_type(*this));}

	///	assigns a vector
		bool assign(const vector_type& v)
		{
		//	assign values on own process
			TGridFunction::assign(v);

		//	set layouts
			copy_layouts_into_vector();

		//	copy storage type
			get_vector().copy_storage_type(v);

		//	we're done
			return true;
		}

		bool assign(const this_type& v)
		{
		//	assign own grid function
			if(!TGridFunction::assign(v))
			{
				UG_ASSERT(0, "Assigning failed.");
				return false;
			}

		//	copy layouts
			copy_layouts_into_vector();

		//	copy storage type
			copy_storage_type(v);

		//	we're done
			return true;
		}

		// sets grid function
		this_type& operator=(const this_type& v)
			{assign(v); return *this;}

		////////////////////////////
		// Storage type
		////////////////////////////

	///	changes storage type
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

	///	sets the storage type
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

	/// changes to the requested storage type if possible
		bool change_storage_type(ParallelStorageType type)
			{return get_vector().change_storage_type(type);}

	/// returns if the current storage type has a given representation
		bool has_storage_type(uint type) const
			{return get_vector().has_storage_type(type);}

	/// returns storage type mask
		uint get_storage_mask() const
			{return get_vector().get_storage_mask(); }

	/// sets the storage type
		void set_storage_type(uint type)
			{get_vector().set_storage_type(type);}

	/// adds the storage type
		void add_storage_type(uint type)
			{get_vector().add_storage_type(type);}

	/// removes the storage type
		void remove_storage_type(uint type)
			{get_vector().remove_storage_type(type);}

	/// copies the storage type from another vector
		void copy_storage_type(const this_type& v)
			{get_vector().copy_storage_type((v.get_vector()));}

		///////////////////////////////
		// index layouts
		///////////////////////////////

		IndexLayout& get_slave_layout()	{return this->get_dof_distribution().get_slave_layout();}
		IndexLayout& get_master_layout()	{return this->get_dof_distribution().get_master_layout();}

		IndexLayout& get_vertical_slave_layout()		{return this->get_dof_distribution().get_vertical_slave_layout();}
		IndexLayout& get_vertical_master_layout()	{return this->get_dof_distribution().get_vertical_master_layout();}

		pcl::ParallelCommunicator<IndexLayout>& get_communicator() {return this->get_dof_distribution().get_communicator();}
		pcl::ProcessCommunicator& get_process_communicator()	{return this->get_dof_distribution().get_process_communicator();}

	protected:
	///	copies references of the layouts from the underlying dof distribution into the vector
		void copy_layouts_into_vector()
		{
		//	copy all horizontal layouts (for all domain decomps)
			get_vector().set_layouts(get_master_layout(), get_slave_layout());

		//	copy vertical layouts
			get_vector().set_vertical_layouts(get_vertical_master_layout(),
			                                 get_vertical_slave_layout());

		//	copy communicator
			get_vector().set_communicator(get_communicator());
			get_vector().set_process_communicator(get_process_communicator());
		}

	/// get own vector
		vector_type& get_vector() {return *(dynamic_cast<vector_type*>(this));}

	///	get own const vector
		const vector_type& get_vector() const {return *(dynamic_cast<const vector_type*>(this));}
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
