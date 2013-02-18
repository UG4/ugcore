/*
 * level_dof_distribution.h
 *
 *  Created on: 24.01.2012
 *      Author: andreasvogel
 */

#ifndef LEVEL_DOF_DISTRIBUTION_H_
#define LEVEL_DOF_DISTRIBUTION_H_

#include "mg_dof_distribution.h"
#include "dof_distribution.h"
#include "managing_dof_distribution.h"
#include "lib_disc/domain_traits.h"

namespace ug{

class LevelMGDoFDistribution : public MGDoFDistribution
{
	friend class LevelDoFDistribution;

	public:
	///	constructor
		LevelMGDoFDistribution(SmartPtr<MultiGrid> spMG,
		                       SmartPtr<MGSubsetHandler> spMGSH,
							   const DoFDistributionInfo& rDDInfo,
		                       bool bGrouped);

	///	removes holes in the index set
	/**
	 * This method removes holes in the index set such that the index set is
	 * contiguous. Therefore, free indices are replaced by those at the end
	 * of the index set. The replacement is stored in the vReplaced vector (and
	 * may be used to adjust associated data, e.g. a grid vector).
	 *
	 * \param[in,out]	vReplaced	vector with all pairs of replacements
	 */
		void defragment(std::vector<std::pair<size_t,size_t> >& vReplaced, int lev);

	///	redistributes all dofs. Resizes associated level vectors afterwards.
		virtual void redistribute_dofs();

	///	register a ManagingDoFDistribution for a given level
	/**	Those dof-distributions will be informed, whenever redistribute_dofs is
	 * executed.
	 * This method is e.g. used by LevelDoFDistribution.*/
		void register_managing_dof_distribution(ManagingDoFDistribution* mdd, int lvl);

	protected:
	///	initializes the indices
		void init();

	///	initializes the indices
		template <typename TBaseElem>
		void init();

	///	removes holes in the index set
	/**
	 * This method removes holes in the index set such that the index set is
	 * contiguous. Therefore, free indices are replaced by those at the end
	 * of the index set. The replacement is stored in the vReplaced vector (and
	 * may be used to adjust associated data, e.g. a grid vector).
	 *
	 * \param[in,out]	vReplaced	vector with all pairs of replacements
	 */
		template <typename TBaseElem>
		void defragment(std::vector<std::pair<size_t,size_t> >& vReplaced, int lev);

	///	adds indices to created objects
	/**
	 * When an element is inserted into the grid, this function is called an
	 * adds needed indices to the grid object.
	 */
		template <typename TBaseElem>
		inline void obj_created(TBaseElem* obj, GeometricObject* pParent = NULL,
		                        bool replacesParent = false);

	///	removes indices, when a grid element is removed
	/**
	 * When a grid element is removed from the grid, this function is called
	 * and takes care about the indices. All indices associated with the element
	 * are removed and stored in a free index container, counters are adjusted.
	 * In general the removal of a grid element will lead to holes in the index
	 * set. Those can be removed by calling defragment.
	 *
	 * \param[in]		obj		grid object that will be removed
	 */
		template <typename TBaseElem>
		inline void obj_to_be_erased(TBaseElem* obj, TBaseElem* replacedBy = NULL);

	public:
		/// grid observer callbacks
		/// \{
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL, bool replacesParent = false);
		virtual void edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent = NULL, bool replacesParent = false);
		virtual void face_created(Grid* grid, Face* f, GeometricObject* pParent = NULL, bool replacesParent = false);
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL, bool replacesParent = false);

		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt, VertexBase* replacedBy = NULL);
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* e, EdgeBase* replacedBy = NULL);
		virtual void face_to_be_erased(Grid* grid, Face* f, Face* replacedBy = NULL);
		virtual void volume_to_be_erased(Grid* grid, Volume* vol, Volume* replacedBy = NULL);
		/// \}

	protected:
		///	returns the number of indices on whole level
		size_t num_indices(const GridLevel& gl) const {
			UG_ASSERT(gl.type() == GridLevel::LEVEL, "Not level.");
			return m_vLev[gl.level()].numIndex;
		}

		///	returns the number of indices on a level and a subset
		size_t num_indices(const GridLevel& gl, const int si) const {
			UG_ASSERT(gl.type() == GridLevel::LEVEL, "Not level.");
			return m_vLev[gl.level()].vNumIndexOnSubset[si];
		}

		///	permutes the indices on a grid level
		void permute_indices(const std::vector<size_t>& vNewInd, int lev);

		/// permutes the indices on a grid level for an base element type
		template <typename TBaseElem>
		void permute_indices(const std::vector<size_t>& vNewInd, int lev);

	///	adjusts storage for requested level
		void level_required(int level);

	///	informations for each level
		std::vector<LevInfo<> > m_vLev;

		LevInfo<>& lev_info(int lev) {return m_vLev[lev];}
		const LevInfo<>& lev_info(int lev) const {return m_vLev[lev];}

		std::vector<ManagingDoFDistribution*>	m_managingDoFDists;

#ifdef UG_PARALLEL
		void create_layouts_and_communicator(int l);

		void create_index_layout(IndexLayout& layout, InterfaceNodeTypes keyType, int l);

		template <typename TBaseElem>
		void add_indices_from_layouts(IndexLayout& indexLayout, InterfaceNodeTypes keyType, int l);

	protected:
		DistributedGridManager* m_pDistGridMgr;
#endif
};



class LevelDoFDistribution : public DoFDistributionInfoProvider,
							 public ManagingDoFDistribution,
							 public DoFDistribution
{
	public:
	///	constructor
		LevelDoFDistribution(SmartPtr<LevelMGDoFDistribution> spLevMGDD,
		                     SmartPtr<SurfaceView> spSurfView,
		                     int level);

		virtual ~LevelDoFDistribution();

	///	returns multigrid
		const MultiGrid& multi_grid() const {return *m_spSurfView->subset_handler()->multi_grid();}


		///////////////////////////////////////
		// Index Access
		///////////////////////////////////////

	/// return the number of dofs distributed
		size_t num_indices() const {return m_spMGDD->num_indices(grid_level());}

	/// return the number of dofs distributed on subset si
		size_t num_indices(int si) const {return m_spMGDD->num_indices(grid_level(), si);}

	///	returns adjacency graph if available
		virtual bool get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

	///	renames the indices
		virtual void permute_indices(const std::vector<size_t>& vIndNew);

	///	removes wholes in index set
		void defragment();

	protected:
		template <typename TBaseElem>
		void get_connections(std::vector<std::vector<size_t> >& vvConnection) const;

		LevInfo<>& lev_info() {return m_spMGDD->lev_info(grid_level().level());}
		const LevInfo<>& lev_info() const {return m_spMGDD->lev_info(grid_level().level());}

#ifdef UG_PARALLEL
	public:
	///	returns the algebra layouts
		const AlgebraLayouts& layouts() const {return lev_info().algebraLayouts;}

	// \TODO: Non-const access should be private or be removed
	public:
	///	returns the algebra layouts
		AlgebraLayouts& layouts() {return lev_info().algebraLayouts;}
#endif

	protected:
	///	MultiGrid Level DoF Distribution
		SmartPtr<LevelMGDoFDistribution> m_spMGDD;

	///	MultiGrid Subset Handler
		SmartPtr<SurfaceView> m_spSurfView;
};

} // end namespace ug

#endif /* LEVEL_DOF_DISTRIBUTION_H_ */
