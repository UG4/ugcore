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

#include "lib_algebra/martin_algebra/template_operations/operations.h"

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
		ParallelGridFunction(	std::string name, approximation_space_type& approxSpace,
								dof_distribution_type& DoFDistr, bool allocate = true)
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


		// set all dofs on level 'level' to value 'w'
		bool set(number w, ParallelStorageType type = PST_CONSISTENT)
			{return TGridFunction::get_vector().set(w, type);}

		////////////////////////////
		// Storage type
		////////////////////////////

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
			{TGridFunction::get_vector().copy_storage_type(v.TGridFunction::get_vector());}

	///////////////////////////////
	// index layouts
	///////////////////////////////

		inline IndexLayout& get_slave_layout()	{return TGridFunction::m_pDoFDistribution->get_slave_layout();}
		inline IndexLayout& get_master_layout()	{return TGridFunction::m_pDoFDistribution->get_master_layout();}
		inline IndexLayout& get_vertical_slave_layout()		{return TGridFunction::m_pDoFDistribution->get_vertical_slave_layout();}
		inline IndexLayout& get_vertical_master_layout()	{return TGridFunction::m_pDoFDistribution->get_vertical_master_layout();}

		inline pcl::ParallelCommunicator<IndexLayout>& get_communicator() {return TGridFunction::m_pDoFDistribution->get_communicator();;}
		inline pcl::ProcessCommunicator& get_process_communicator()	{return TGridFunction::m_pDoFDistribution->get_process_communicator();}

	///////////////////////////////
	// help functions
	///////////////////////////////
	protected:
		// forward all layouts from DoFDistribution to algebra Vector
		void set_layouts()
		{
		//	since we inherite from a Parallel<TVector>, the following functions are inherited
			set_slave_layout(get_slave_layout());
			set_master_layout(get_master_layout());
			set_vertical_slave_layout(get_vertical_slave_layout());
			set_vertical_master_layout(get_vertical_master_layout());

			set_communicator(get_communicator());
			set_process_communicator(get_process_communicator());
		}


};




} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__PARALLELIZATION__PARALLEL_GRID_FUNCTION__ */
