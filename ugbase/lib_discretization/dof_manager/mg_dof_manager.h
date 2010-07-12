/*
 * mg_dof_manager.h
 *
 *  Created on: 12.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__
#define __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__

#include "lib_grid/lib_grid.h"
#include "lib_discretization/geom_object_container/surface_view.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"
#include "lib_algebra/multi_index/multi_indices.h"

namespace ug{

class P1ConformDoFDistributor
{
	public:
		// type of multiindex used
		typedef MultiIndex<1> index_type;

		// value container for element local indices
		typedef std::vector<index_type> local_index_type;

	public:
		P1ConformDoFDistributor(GeometricObjectCollection& goc, ISubsetHandler& sh)
		: m_pGoc(&goc), m_pISubsetHandler(&sh), m_numDoFs(0)
		{
			m_vFunctionInfo.clear();
			m_vSubsetInfo.clear();
		}

		bool add_discrete_function(std::string name, LocalShapeFunctionSetID id, int dim)
		{
			// for a P1 dof manager only Lagrange P1 function space is permitted
			if(id != LSFS_LAGRANGEP1)
				{UG_LOG("P1ConformDoFDistributor: Only LSFS_LAGRANGEP1 functions are supported.\n"); return false;}

			// add to function list
			m_vFunctionInfo.push_back(FunctionInfo(name, dim));

			return true;
		}

		/// number of subsets
		inline size_t num_subsets() const {return m_pISubsetHandler->num_subset_infos();}

		/// number of discrete functions this dof distributor handles
		inline size_t num_fct() const {return m_vFunctionInfo.size();}

		/// returns the name of the discrete function nr_fct
		std::string get_name(size_t fct) const {return m_vFunctionInfo[fct].name;}

		/// return the number of dofs distributed
		inline size_t num_dofs() const {return m_numDoFs;}

		// distribute dofs for given goc
		bool distribute_dofs()
		{
			// Attach indices
			subset_required(num_subsets());

			// iterators
			geometry_traits<VertexBase>::iterator iter, iterBegin, iterEnd;

			// loop subsets
			size_t i = 0;

			for(size_t si = 0; si < num_subsets(); ++si)
			{
				iterBegin = m_pGoc->begin<VertexBase>(si);
				iterEnd =  m_pGoc->end<VertexBase>(si);

				for(iter = iterBegin; iter != iterEnd; ++iter)
				{
					VertexBase* vrt = *iter;
					(m_vSubsetInfo[si].aaDoF)[vrt] = i;
					i += num_fct();
				}
			}
			m_numDoFs = i;

			UG_LOG( std::setw(8) << i << " DoF indices distributed [each carrying " << 1 << " component(s)]" << std::endl);

			return true;
		}

		template<typename TElem>
		size_t get_multi_indices(TElem* elem, size_t fct, local_index_type& ind, size_t offset) const
		{
			UG_ASSERT(fct < num_fct(), "Function not defined in DoF Manager");
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;

			for(int i = 0; i < ref_elem_type::num_corners; ++i)
			{
				VertexBase* vrt = elem->vertex(i);
				int si = m_pISubsetHandler->get_subset_index(vrt);

				ind[offset + i][0] = (m_vSubsetInfo[si].aaDoF)[vrt] + fct;
			}

			return ref_elem_type::num_corners;
		}

		template<typename TGeomObj>
		size_t get_multi_indices_of_geom_obj(TGeomObj* obj, size_t fct, local_index_type& ind, size_t offset) const
		{
			if(geometry_traits<TGeomObj>::BASE_OBJECT_TYPE_ID == VERTEX)
			{
				VertexBase* vrt = (VertexBase*)obj;
				int si = m_pISubsetHandler->get_subset_index(vrt);
				ind[offset + 0][0] = (m_vSubsetInfo[si].aaDoF)[vrt] + fct;
				return 1;
			}
			else return 0;
		}


	protected:
		inline void subset_required(size_t si)
		{
			m_pISubsetHandler->enable_subset_attachments(true);

			// Create level dof distributors
			for(size_t s = m_vSubsetInfo.size(); s <= si; ++s)
			{
				m_vSubsetInfo.push_back(SubsetInfo());
				m_pISubsetHandler->attach_to<VertexBase>(m_vSubsetInfo[si].aDoF, si);
				m_vSubsetInfo[si].aaDoF.access(*m_pISubsetHandler, m_vSubsetInfo[si].aDoF, si);
			}
		}


	protected:
		struct FunctionInfo
		{
			FunctionInfo(std::string name_, int dim_) : name(name_), dim(dim_){};
			std::string name;
			int dim;
		};

		struct SubsetInfo
		{
			typedef ug::Attachment<size_t> ADoF;
			typedef ISubsetHandler::AttachmentAccessor<VertexBase, ADoF> attachment_accessor_type;

			attachment_accessor_type aaDoF;
			ADoF aDoF;
		};

	protected:
		// geometric object collection for this Distributor
		GeometricObjectCollection* m_pGoc;

		// subset handler for this distributor
		ISubsetHandler* m_pISubsetHandler;

		// informations about Functions
		std::vector<FunctionInfo> m_vFunctionInfo;

		// informations about Subsets
		std::vector<SubsetInfo> m_vSubsetInfo;

		// number of distributed dofs
		size_t m_numDoFs;

};



template <typename TDoFDistributor>
class MGDoFManager
{
	public:
		MGDoFManager() : m_pMGSubsetHandler(NULL), m_pMultiGrid(NULL), m_pSurfaceView(NULL)
		{
			m_vLevelDoFDistributor.clear();
		};

		MGDoFManager(MultiGridSubsetHandler& mgsh)
		: m_pMGSubsetHandler(&mgsh), m_pMultiGrid(m_pMGSubsetHandler->get_assigned_multi_grid())
		{
			m_vLevelDoFDistributor.clear();
			m_pSurfaceView = new SurfaceView(*m_pMultiGrid);
		};

		// set multi grid subset handler
		void assign_multi_grid_subset_handler(MultiGridSubsetHandler& mgsh)
		{
			m_pMGSubsetHandler = &mgsh;
			m_pMultiGrid = m_pMGSubsetHandler->get_assigned_multi_grid();
			m_pSurfaceView = new SurfaceView(*m_pMultiGrid);
			// TODO: init surface view
		}

		// number of levels
		inline size_t num_levels() const {return m_pMultiGrid->num_levels();}

		// distribute dofs on all levels + surface level
		bool distribute_dofs()
		{
			level_required(num_levels());

			// distribute on level grids
			for(size_t l = 0; l < num_levels(); ++l)
			{
				if(!m_vLevelDoFDistributor[l]->distribute_dofs()) return false;
			}

			// distribute on surface grid
			if(!m_pSurfaceDoFDistributor->distribute_dofs()) return false;

			return true;
		}

		~MGDoFManager()
		{
			if(m_pSurfaceView != NULL) delete m_pSurfaceView;
			for(size_t l = 0; l < m_vLevelDoFDistributor.size(); ++l)
			{
				delete m_vLevelDoFDistributor[l];
			}
		}

	protected:
		inline void level_required(size_t level)
		{
			// Create surface dof distributor
			if(m_pSurfaceDoFDistributor == NULL)
				m_pSurfaceDoFDistributor = new TDoFDistributor(m_pSurfaceView->get_geometric_object_collection(), m_pSurfaceView);
			else

			// Create level dof distributors
			for(size_t l = m_vLevelDoFDistributor.size(); l <= level; ++l)
			{
				m_vLevelDoFDistributor.push_back(new TDoFDistributor(m_pMGSubsetHandler->get_goc_by_level(l), m_pMGSubsetHandler));
			}
		}


	protected:
		// MultiGridSubsetHandler this DofManager works on
		MultiGridSubsetHandler* m_pMGSubsetHandler;

		// MultiGrid associated to the SubsetHandler
		MultiGrid* m_pMultiGrid;

		// Surface View
		SurfaceView* m_pSurfaceView;

		// Level DoF Distributors
		std::vector<TDoFDistributor*> m_vLevelDoFDistributor;

		// Surface Grid DoF Distributor
		TDoFDistributor* m_pSurfaceDoFDistributor;
};









} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__DOF_MANAGER__MG_DOF_MANAGER__ */
