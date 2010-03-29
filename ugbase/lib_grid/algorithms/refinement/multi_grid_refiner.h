//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d09

#ifndef __H__LIB_GRID__MULTI_GRID_REFINER__
#define __H__LIB_GRID__MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"

namespace ug
{

class MultiGridRefiner : public GridObserver
{
	public:
	/**	Those marks are used to mark the status of an element.
	 *	constants must be in the range 0 - 0xFD*/
		enum StatusMark
		{
			SM_NONE = 0,
			SM_REGULAR = 1,
			SM_COPY = 2,
			SM_IRREGULAR = 3,
		};

	public:
		MultiGridRefiner();
		MultiGridRefiner(MultiGrid& mg);
		~MultiGridRefiner();
		
		virtual void registered_at_grid(Grid* grid);
		virtual void unregistered_from_grid(Grid* grid);

		void assign_grid(MultiGrid& mg);

	////////////////////////////////
	//	marks
		void clear_marks();

//TODO:	Add an optional parameter anisotropic = false, that determines whether an
//		element is part of an anisotropic refinement.
//		If all marked lower-dimensional elements are marked anisotropic (at least two), then
//		the resulting elements will be marked anisotropic, too.
		template <class TElem>
		inline void mark_for_refinement(TElem* elem)
		{
			m_selMarks.select(elem);
		}

	///	the value-type of TIterator has to be a pointer to a type derived from either EdgeBase, Face or Volume.
		template <class TIterator>
		inline void mark_for_refinement(const TIterator& iterBegin,
										const TIterator& iterEnd)
					{m_selMarks.select(iterBegin, iterEnd);}

		template <class TElem>
		inline bool is_marked(TElem* elem)	{return m_selMarks.is_selected(elem);}

	////////////////////////////////
	//	refine
	///	performs refinement on the marked elements.
		void refine();

	////////////////////////////////
	//	settings
	///	determines how many unrefined neighbours will be copied to the next level in subsequent refinement-steps.
		inline void set_copy_range(int range)	{m_copyRange = range;}
	///	returns how many unrefined neighbours are copied to the next level in subsequent refinement-steps.
		inline int get_copy_range()				{return m_copyRange;}

	////////////////////////////////
	//	element-status
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.
		inline int get_status(VertexBase* e)	{return m_aaIntVRT[e] & MR_STATUS;}
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.	
		inline int get_status(EdgeBase* e)		{return m_aaIntEDGE[e] & MR_STATUS;}
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.
		inline int get_status(Face* e)			{return m_aaIntFACE[e] & MR_STATUS;}
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.	
		inline int get_status(Volume* e)		{return m_aaIntVOL[e] & MR_STATUS;}

	protected:
	/**	Those marks are used to mark the refinement rule that will be
	 *	applied to an element.
	 *	constants must be in the range 0x00FF - 0xFF00*/
		enum RefinementMark
		{
			RM_NONE = 0x00FF,
			RM_REGULAR = 1 << 8,
			RM_ANISOTROPIC = 1 << 9,	//not yet used
			RM_COPY = 1 << 10,
			RM_IRREGULAR = 1 << 11,
			RM_UNKNOWN = 1 << 12
		};

	/**	Those constants define the range in which associated marks lie.*/
		enum MarkRanges
		{
			MR_STATUS = 0xFF,
			MR_REFINEMENT = 0xFF00
		};

	protected:
		virtual void collect_objects_for_refine();
	///	this method helps derived classes to perform operations directly before actual element refinment is performed.
	/**	Called from the refine() method in each refinement-iteration after
	 *	collect_objects_for_refine().
	 *	Default implementation is empty.*/
		virtual void refinement_step_begins()	{};

	///	this method helps derived classes to perform operations directly after actual element refinment took place.
	/**	Called from the refine() method in each refinement-iteration after
	 *	all scheduled elements had been refined.
	 *	The refine process will either terminate after this method or will
	 *	start a new iteration, if new elements had been marked during refine.
	 *	Default implementation is empty.*/
		virtual void refinement_step_ends()		{};
		
		void adjust_initial_selection();
		void select_closure(std::vector<VertexBase*>& vVrts);
		void select_copy_elements(std::vector<VertexBase*>& vVrts);
		
		inline void set_status(VertexBase* e, StatusMark mark)	{m_aaIntVRT[e] = (m_aaIntVRT[e] & ~MR_STATUS) | mark;}
		inline void set_status(EdgeBase* e, StatusMark mark)	{m_aaIntEDGE[e] = (m_aaIntEDGE[e] & ~MR_STATUS) | mark;}
		inline void set_status(Face* e, StatusMark mark)		{m_aaIntFACE[e] = (m_aaIntFACE[e] & ~MR_STATUS) | mark;}
		inline void set_status(Volume* e, StatusMark mark)		{m_aaIntVOL[e] = (m_aaIntVOL[e] & ~MR_STATUS) | mark;}

		virtual void set_rule(VertexBase* e, RefinementMark mark)	{m_aaIntVRT[e] = (m_aaIntVRT[e] & ~MR_REFINEMENT) | mark;}
		virtual void set_rule(EdgeBase* e, RefinementMark mark)		{m_aaIntEDGE[e] = (m_aaIntEDGE[e] & ~MR_REFINEMENT) | mark;}
		virtual void set_rule(Face* e, RefinementMark mark)			{m_aaIntFACE[e] = (m_aaIntFACE[e] & ~MR_REFINEMENT) | mark;}
		virtual void set_rule(Volume* e, RefinementMark mark)		{m_aaIntVOL[e] = (m_aaIntVOL[e] & ~MR_REFINEMENT) | mark;}

		inline int get_rule(VertexBase* e)	{return m_aaIntVRT[e] & MR_REFINEMENT;}
		inline int get_rule(EdgeBase* e)	{return m_aaIntEDGE[e] & MR_REFINEMENT;}
		inline int get_rule(Face* e)		{return m_aaIntFACE[e] & MR_REFINEMENT;}
		inline int get_rule(Volume* e)		{return m_aaIntVOL[e] & MR_REFINEMENT;}

	protected:
		MultiGrid*	m_pMG;

	//	selection-marks
		Selector	m_selMarks;

	//	copy-range
		int			m_copyRange;

	//	status-marks
		AInt		m_aInt;
		Grid::VertexAttachmentAccessor<AInt>	m_aaIntVRT;
		Grid::EdgeAttachmentAccessor<AInt>		m_aaIntEDGE;
		Grid::FaceAttachmentAccessor<AInt>		m_aaIntFACE;
		Grid::VolumeAttachmentAccessor<AInt>	m_aaIntVOL;
};

}//	end of namespace

#endif
