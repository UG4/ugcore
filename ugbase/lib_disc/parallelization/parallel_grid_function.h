/*
 * parallel_grid_function.h
 *
 *  Created on: 19.7.2010
 *      Author: A. Vogel
 */

#ifndef __H__UG__LIB_DISC__PARALLELIZATION__PARALLEL_GRID_FUNCTION__
#define __H__UG__LIB_DISC__PARALLELIZATION__PARALLEL_GRID_FUNCTION__

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

		// dof distribution type
		typedef typename TGridFunction::dof_distribution_type dof_distribution_type;

		using TGridFunction::clone_pattern;

	public:
	/// Default Constructor
		ParallelGridFunction(SmartPtr<approximation_space_type> approxSpace)
			: TGridFunction(approxSpace)
		{
			copy_layouts_into_vector();
			this->set_storage_type(PST_UNDEFINED);
		};

	/// Initializing Constructor
		ParallelGridFunction(SmartPtr<approximation_space_type> approxSpace,
		                     SmartPtr<dof_distribution_type> DoFDistr)
			: TGridFunction(approxSpace, DoFDistr)
		{
			copy_layouts_into_vector();
			this->set_storage_type(PST_UNDEFINED);
		};


		///////////////////////////////
		// overloaded functions
		///////////////////////////////
		virtual void clone_pattern(const this_type& v)
		{
		//	clone process local grid function
			TGridFunction::clone_pattern(v);

		//	copy layouts
			copy_layouts_into_vector();

		//	set storage type to undefined
			this->set_storage_type(PST_UNDEFINED);
		}

	/// clone
		this_type& clone() {return *(new this_type(*this));}

	///	assigns a vector
		void assign(const vector_type& v)
		{
		//	assign values on own process
			TGridFunction::assign(v);

		//	set layouts
			copy_layouts_into_vector();

		//	copy storage type
			vector_type::copy_storage_type(v);
		}

		void assign(const this_type& v)
		{
		//	assign own grid function
			TGridFunction::assign(v);

		//	copy layouts
			copy_layouts_into_vector();

		//	copy storage type
			copy_storage_type(v);
		}

		// sets grid function
		this_type& operator=(const this_type& v)
			{assign(v); return *this;}

		////////////////////////////
		// Storage type
		////////////////////////////

	/// copies the storage type from another vector
		void copy_storage_type(const this_type& v)
			{vector_type::copy_storage_type((v.get_vector()));}

	protected:
	///	copies references of the layouts from the underlying dof distribution into the vector
		void copy_layouts_into_vector()
		{
		//	copy all horizontal layouts (for all domain decomps)
			vector_type::set_layouts(this->m_spDD->master_layout(), this->m_spDD->slave_layout());

		//	copy vertical layouts
			vector_type::set_vertical_layouts(this->m_spDD->vertical_master_layout(),
			                                  this->m_spDD->vertical_slave_layout());

		//	copy communicator
			vector_type::set_communicator(this->m_spDD->communicator());
			vector_type::set_process_communicator(this->m_spDD->process_communicator());
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

#endif /* __H__UG__LIB_DISC__PARALLELIZATION__PARALLEL_GRID_FUNCTION__ */
