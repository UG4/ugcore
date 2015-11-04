#ifndef __H__UG_regular_refiner_multi_grid
#define __H__UG_regular_refiner_multi_grid

#include "lib_grid/tools/selector_multi_grid.h"
#include "refiner_interface.h"

namespace ug{

class UG_API RegularRefiner_MultiGrid : public IRefiner{
	public:
		RegularRefiner_MultiGrid ();
		RegularRefiner_MultiGrid (MultiGrid* pmg);

		void 		set_grid (MultiGrid* pmg);
		MultiGrid*	multi_grid () const;

		virtual Grid* get_associated_grid ();
		virtual Grid* grid ();

		virtual bool adaptivity_supported () const;
		virtual bool coarsening_supported () const;

		virtual bool mark (Vertex* v, RefinementMark refMark = RM_REFINE);
		virtual bool mark (Edge* e, RefinementMark refMark = RM_REFINE);
		virtual bool mark (Face* f, RefinementMark refMark = RM_REFINE);
		virtual bool mark (Volume* v, RefinementMark refMark = RM_REFINE);

		virtual void mark_neighborhood (size_t numIterations,
										RefinementMark refMark,
										bool sideNbrsOnly);

		virtual RefinementMark get_mark (Vertex* v);
		virtual RefinementMark get_mark (Edge* e);
		virtual RefinementMark get_mark (Face* f);
		virtual RefinementMark get_mark (Volume* v);

		virtual bool save_marks_to_file (const char* filename);

		template <class TElem>
		bool is_closure (TElem* elem);

	protected:
		enum Marks{
			NONE = RM_NONE,
			COPY = RM_COPY,
			REGULAR = RM_REFINE,
			ANISOTROPIC = RM_ANISOTROPIC,
			LIFT = REGULAR | ANISOTROPIC | COPY
		};

		virtual void perform_refinement ();
		virtual bool perform_coarsening ();

	///	returns the number of locally marked edges on all levels of the hierarchy
		virtual void num_marked_edges_local (std::vector<int>& numMarkedEdgesOut);
	///	returns the number of locally marked faces on all levels of the hierarchy
		virtual void num_marked_faces_local (std::vector<int>& numMarkedFacesOut);
	///	returns the number of locally marked volumes on all levels of the hierarchy
		virtual void num_marked_volumes_local (std::vector<int>& numMarkedVolsOut);


		template <class TElem>
		void adjust_side_states (
				size_t lvl,
				uint considerElemMarks,
				uint ignoreSideMarks,
				RefinementMark newSideMark,
				bool closure);

		template <class TElem>
		void copy_state_to_sides (
				size_t lvl,
				uint considerElemMarks,
				bool closure);

		template <class TSide>
		void adjust_side_of_states (
				size_t lvl,
				uint considerSideMarks,
				uint ignoreElemMarks,
				RefinementMark newElemMark,
				bool closure);

		template <class TElem>
		void clear_dummies ();

		template <class TElem>
		void mark_by_level_discrepancy (
				int lvl,
				Grid::VertexAttachmentAccessor<AInt> aaLvl);

		void collect_objects_for_refine ();

		bool refinement_is_allowed(Vertex* v);
		bool refinement_is_allowed(Edge* e);
		bool refinement_is_allowed(Face* f);
		bool refinement_is_allowed(Volume* v);

	private:
		RegularRefiner_MultiGrid (const RegularRefiner_MultiGrid&)	{}

		MGSelector						m_marks;
		MultiGrid*						m_pMG;

		ABool									m_aClosure;
		MultiElementAttachmentAccessor<AByte>	m_aaClosure;
};

}//	end of namespace

#endif	//__H__UG_regular_refiner_multi_grid
