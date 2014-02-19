/*
 * grid_function_impl.h
 *
 *  Created on: 13.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__

#include "grid_function.h"

#include "lib_algebra/algebra_type.h"
#include "adaption_surface_grid_function.h"
#include "lib_grid/algorithms/serialization.h"
#include "common/profiler/profiler.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
	#include "lib_algebra/parallelization/parallelization.h"
	#include "lib_grid/parallelization/util/compol_copy_attachment.h"
#endif

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// GridFunction : init
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
             SmartPtr<DoFDistribution> spDoFDistr, bool bManage)
{
	init(spApproxSpace, spDoFDistr, bManage);
};

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, bool bManage)
{
	init(spApproxSpace, spApproxSpace->dof_distribution(GridLevel(GridLevel::TOP, GridLevel::SURFACE, false)), bManage);
};

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace, int level, bool bManage)
{
	init(spApproxSpace, spApproxSpace->dof_distribution(GridLevel(level, GridLevel::SURFACE, false)), bManage);
};

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>::
GridFunction(SmartPtr<approximation_space_type> spApproxSpace, const GridLevel& gl, bool bManage)
{
	init(spApproxSpace, spApproxSpace->dof_distribution(gl), bManage);
};

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
init(SmartPtr<ApproximationSpace<TDomain> > spApproxSpace,
     SmartPtr<DoFDistribution> spDoFDistr, bool bManage)
{
	m_spApproxSpace = spApproxSpace;
	m_spDD = spDoFDistr;
	m_bManaged = bManage;
	this->set_dof_distribution_info(m_spApproxSpace->dof_distribution_info());
	m_spAdaptGridFct = SPNULL;

//	check correct passings
	if(m_spDD.invalid()) UG_THROW("GridFunction: DoF Distribution is null.");
	if(m_spApproxSpace.invalid()) UG_THROW("GridFunction: ApproxSpace is null.");

//	check correct choice of compile-time algebra
	check_algebra();

//	resize the vector to correct size
	resize_values(num_indices());

	if(bManage) {
		//	registered as managed by dof distribution
		m_spDD->manage_grid_function(*this);

		//	register to observe grid
		register_at_adaption_msg_hub();
	}


#ifdef UG_PARALLEL
//	set layouts
	this->set_layouts(m_spDD->layouts());

//	set storage type
	this->set_storage_type(PST_UNDEFINED);
#endif
};

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::check_algebra()
{
//	get blocksize of algebra
	const int blockSize = algebra_type::blockSize;

//	a)	If blocksize fixed and > 1, we need grouping.
	if(blockSize > 1 && !this->m_spDD->grouped())
	{
		UG_THROW("Fixed block algebra needs grouped dofs.");
	}
//	b) 	If blocksize flexible, we group
	else if (blockSize == AlgebraType::VariableBlockSize
			&& !this->m_spDD->grouped())
	{
		UG_THROW("Variable block algebra needs grouped dofs.");
	}
//	c)	If blocksize == 1, we do not group. This will allow us to handle
//		this case for any problem.
	else if (blockSize == 1 && this->m_spDD->grouped())
	{
		UG_THROW("block 1x1 algebra needs non-grouped dofs.");
	}
}

////////////////////////////////////////////////////////////////////////////////
// GridFunction : cloning
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
clone_pattern(const this_type& v)
{
//	init normally
	init(v.m_spApproxSpace, v.m_spDD, v.m_bManaged);

#ifdef UG_PARALLEL
//	copy storage type
	this->set_storage_type(v.get_storage_mask());
	this->set_layouts(v.layouts());
#endif
};


template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::assign(const vector_type& v)
{
//	check size
	if(v.size() != vector_type::size())
		UG_THROW("GridFunction: Assigned vector has incorrect size.");

//	assign vector
	*(dynamic_cast<vector_type*>(this)) = v;

#ifdef UG_PARALLEL
//	copy storage type
	this->set_storage_type(v.get_storage_mask());
#endif
}

template <typename TDomain, typename TAlgebra>
void GridFunction<TDomain, TAlgebra>::assign(const this_type& v)
{
//	clone pattern
	clone_pattern(v);

//  copy values
	*(dynamic_cast<vector_type*>(this)) = *dynamic_cast<const vector_type*>(&v);
}

template <typename TDomain, typename TAlgebra>
GridFunction<TDomain, TAlgebra>*
GridFunction<TDomain, TAlgebra>::virtual_clone_without_values() const
{
	GridFunction<TDomain, TAlgebra>* p =
		new GridFunction<TDomain, TAlgebra>(m_spApproxSpace, m_spDD, m_bManaged);
	if(p->size() != this->size())
		p->resize(this->size());
#ifdef UG_PARALLEL
	p->set_layouts(this->layouts());
#endif

	return p;
}

////////////////////////////////////////////////////////////////////////////////
// GridFunction : dof distribution callbacks
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
resize_values(size_t s, number defaultValue)
{
//	remember old values
	const size_t oldSize = vector_type::size();

//	resize vector
	vector_type::resize_sloppy(s);

//	set vector to zero-values
	for(size_t i = oldSize; i < s; ++i)
		this->operator[](i) = defaultValue;
}

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
permute_values(const std::vector<size_t>& vIndNew)
{
//	check sizes
	if(vIndNew.size() != this->size())
		UG_THROW("GridFunction::permute_values: For a permutation the"
				" index set must have same cardinality as vector.");

// \todo: avoid tmp vector, only copy values into new vector and use that one
//	create tmp vector
	vector_type vecTmp; vecTmp.resize(this->size());
#ifdef UG_PARALLEL
//	copy storage type
	vecTmp.set_storage_type(this->get_storage_mask());
#endif

//	loop indices and copy values
	for(size_t i = 0; i < vIndNew.size(); ++i)
		vecTmp[vIndNew[i]] = this->operator[](i);

//	copy tmp vector into this vector
	this->assign(vecTmp);
}

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
copy_values(const std::vector<std::pair<size_t, size_t> >& vIndexMap,bool bDisjunct)
{
//	disjunct case
	if(bDisjunct)
		for(size_t i = 0; i < vIndexMap.size(); ++i)
			this->operator[](vIndexMap[i].second)
				= this->operator[](vIndexMap[i].first);
	else {
		typedef typename vector_type::value_type value_type;
		std::vector<value_type> values;
		values.resize(vIndexMap[vIndexMap.size()-1].first);
		for(size_t i = 0; i < vIndexMap.size(); ++i){ 
			const size_t index = vIndexMap[i].first;
			if (index>=values.size()) values.resize(index+1);
			values[index] = this->operator[](index);
		}
		for(size_t i = 0; i < vIndexMap.size(); ++i)
			this->operator[](vIndexMap[i].second)
				= values[vIndexMap[i].first];
	}
}

////////////////////////////////////////////////////////////////////////////////
// GridFunction : grid adaption
////////////////////////////////////////////////////////////////////////////////

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
register_at_adaption_msg_hub()
{
//	register function for grid adaption
	SPMessageHub msgHub = domain()->grid()->message_hub();
	m_spGridAdaptionCallbackID =
		msgHub->register_class_callback(this,
		&this_type::grid_changed_callback);

	m_spGridDistributionCallbackID =
		msgHub->register_class_callback(this,
		&this_type::grid_distribution_callback);
}

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
grid_changed_callback(const GridMessage_Adaption& msg)
{
	// before adaption begins: copy values into grid attachments
	if(msg.adaption_begins()){
		// prepare
		m_spAdaptGridFct = SmartPtr<AdaptionSurfaceGridFunction<TDomain> >(
							new AdaptionSurfaceGridFunction<TDomain>(this->domain()));
		m_spAdaptGridFct->copy_from_surface(*this);
	}

	// before coarsening: restrict values
	if(msg.coarsening() && msg.step_begins()){
		#ifdef UG_PARALLEL
		//	since ghosts may exist in a parallel environment and since those ghosts
		//	may be removed during coarsening, we have to make sure, that the correct
		//	values are stored in those ghosts before restriction is performed.
			Grid& grid = *domain()->grid();
			typedef typename AdaptionSurfaceGridFunction<TDomain>::AValues AValues;
			if(m_spDDI->max_dofs(VERTEX)){
				ComPol_CopyAttachment<VertexLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<VertexLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
			}
			if(m_spDDI->max_dofs(EDGE)){
				ComPol_CopyAttachment<EdgeLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<EdgeLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
			}
			if(m_spDDI->max_dofs(FACE)){
				ComPol_CopyAttachment<FaceLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<FaceLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
			}
			if(m_spDDI->max_dofs(VOLUME)){
				ComPol_CopyAttachment<VolumeLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<VolumeLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
			}
		#endif
		m_spAdaptGridFct->do_restrict(msg);
	}

	// after refinement: prolongate values
	if(msg.refinement() && msg.step_ends()){
		m_spAdaptGridFct->prolongate(msg);
	}

	// at end of adaption: copy values back into algebra vector
	if(msg.adaption_ends())
	{
		// all grid functions must resize to the current number of dofs
		resize_values(num_indices());

		#ifdef UG_PARALLEL
		//	set layouts
		this->set_layouts(m_spDD->layouts());
		#endif


		m_spAdaptGridFct->copy_to_surface(*this);
		m_spAdaptGridFct = SPNULL;
	}
}

template <typename TDomain, typename TAlgebra>
void
GridFunction<TDomain, TAlgebra>::
grid_distribution_callback(const GridMessage_Distribution& msg)
{
	PROFILE_FUNC();

	#ifdef UG_PARALLEL
	GridDataSerializationHandler& sh = msg.serialization_handler();

	switch(msg.msg()){
		case GMDT_DISTRIBUTION_STARTS:{
			m_preDistStorageType = this->get_storage_mask();
			if(!(this->has_storage_type(PST_CONSISTENT) || this->has_storage_type(PST_UNDEFINED))){
				this->change_storage_type(PST_CONSISTENT);
			}

			m_spAdaptGridFct = SmartPtr<AdaptionSurfaceGridFunction<TDomain> >(
						new AdaptionSurfaceGridFunction<TDomain>(this->domain(), false));
			m_spAdaptGridFct->copy_from_surface(*this);
			Grid& grid = *domain()->grid();

			typedef typename AdaptionSurfaceGridFunction<TDomain>::AValues AValues;

			if(m_spDDI->max_dofs(VERTEX)){
				ComPol_CopyAttachment<VertexLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<VertexLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
				sh.add(GeomObjAttachmentSerializer<Vertex, AValues>::
							create(grid, m_spAdaptGridFct->value_attachment()));
			}
			if(m_spDDI->max_dofs(EDGE)){
				ComPol_CopyAttachment<EdgeLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<EdgeLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
				sh.add(GeomObjAttachmentSerializer<EdgeBase, AValues>::
							create(grid, m_spAdaptGridFct->value_attachment()));
			}
			if(m_spDDI->max_dofs(FACE)){
				ComPol_CopyAttachment<FaceLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<FaceLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				com.communicate();
				sh.add(GeomObjAttachmentSerializer<Face, AValues>::
							create(grid, m_spAdaptGridFct->value_attachment()));
			}
			if(m_spDDI->max_dofs(VOLUME)){
				ComPol_CopyAttachment<VolumeLayout, AValues> compol(grid, m_spAdaptGridFct->value_attachment());
				pcl::InterfaceCommunicator<VolumeLayout> com;
				com.exchange_data(grid.distributed_grid_manager()->grid_layout_map(),
								  INT_V_SLAVE, INT_V_MASTER, compol);
				sh.add(GeomObjAttachmentSerializer<Volume, AValues>::
							create(grid, m_spAdaptGridFct->value_attachment()));
			}
		}break;

		case GMDT_DISTRIBUTION_STOPS:
			{
				PROFILE_BEGIN(grid_func_distribution_stops)
				// all grid functions must resize to the current number of dofs
				resize_values(num_indices());

				#ifdef UG_PARALLEL
				//	set layouts
				this->set_layouts(m_spDD->layouts());
				#endif

				m_spAdaptGridFct->copy_to_surface(*this);
				m_spAdaptGridFct = SPNULL;

				if(m_preDistStorageType != this->get_storage_mask()){
					if((m_preDistStorageType & PST_ADDITIVE) == PST_ADDITIVE)
						this->change_storage_type(PST_ADDITIVE);
					else if((m_preDistStorageType & PST_UNIQUE) == PST_UNIQUE)
						this->change_storage_type(PST_UNIQUE);
					else{
						UG_THROW("Can't reestablish storage type!");
					}
				}

				PROFILE_END();
			}break;

		default:
			break;
	}
	#endif
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_IMPL__ */
