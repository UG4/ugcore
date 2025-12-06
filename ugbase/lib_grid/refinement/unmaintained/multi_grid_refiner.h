øunused
/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__LIB_GRID__MULTI_GRID_REFINER__
#define __H__LIB_GRID__MULTI_GRID_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{
/**	This class is only a prototype for an adaptive multigrid refiner
 *  with regular closure and copy elements.
 *
 *  It is currently not used in the simulation system.
 *  Instead, on can use e.g. the HangingNodeRefiner_MultiGrid.
 *
 *  Before this class should be used in simulations, one should
 *  derive it from IRefiner and adjust the mark methods as well
 *  as the refinement-callbacks.
 *
 *  Please note that the functionality of this class is limited.
 *  Only elements in the highest grid hierarchy may be refined.
 */
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
		~MultiGridRefiner() override;
		
		virtual void grid_to_be_destroyed(Grid* grid);

		void assign_grid(MultiGrid& mg);

	////////////////////////////////
	//	marks
		void clear_marks();

//TODO:	Add an optional parameter anisotropic = false, that determines whether an
//		element is part of an anisotropic refinement.
//		If all marked lower-dimensional elements are marked anisotropic (at least two), then
//		the resulting elements will be marked anisotropic, too.
		template <typename TElem>
		inline void mark_for_refinement(TElem* elem)
		{
			m_selMarks.select(elem);
		}

	///	the value-type of TIterator has to be a pointer to a type derived from either Edge, Face or Volume.
		template <typename TIterator>
		inline void mark_for_refinement(const TIterator& iterBegin,
										const TIterator& iterEnd)
					{m_selMarks.select(iterBegin, iterEnd);}

		template <typename TElem>
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
		inline int get_status(Vertex* e) {return m_aaIntVRT[e] & MR_STATUS;}
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.	
		inline int get_status(Edge* e) {return m_aaIntEDGE[e] & MR_STATUS;}
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.
		inline int get_status(Face* e) {return m_aaIntFACE[e] & MR_STATUS;}
	///	returns a constant enumerated in MultiGridRefiner::StatusMark.	
		inline int get_status(Volume* e) {return m_aaIntVOL[e] & MR_STATUS;}

	protected:
	/**	Those marks are used to mark the refinement rule that will be
	 *	applied to an element.
	 *	constants must be in the range 0x00FF - 0xFF00*/
		enum RefinementMark
		{
			RM_NONE = 0x00FF,
			RM_REFINE = 1 << 8,
			RM_ANISOTROPIC = 1 << 9,	//not yet used
			RM_COPY = 1 << 10,
			RM_IRREGULAR = 1 << 11,
			RM_UNKNOWN = 1 << 12
//			RM_FIXED = 1 << 13
		};

	/**	Those constants define the range in which associated marks lie.*/
		enum MarkRanges
		{
			MR_STATUS = 0xFF,
			MR_REFINEMENT = 0xFF00
		};

	protected:
	///	performs registration and deregistration at a grid.
	/**	call set_grid(nullptr) to unregister the observer from a grid.*/
		void set_grid(Grid* grid);
		
		virtual void collect_objects_for_refine();
	///	this method helps derived classes to perform operations directly before actual element refinment is performed.
	/**	Called from the refine() method in each refinement-iteration after
	 *	collect_objects_for_refine().
	 *	Default implementation is empty.*/
		virtual void refinement_step_begins() {};

	///	this method helps derived classes to perform operations directly after actual element refinment took place.
	/**	Called from the refine() method in each refinement-iteration after
	 *	all scheduled elements had been refined.
	 *	The refine process will either terminate after this method or will
	 *	start a new iteration, if new elements had been marked during refine.
	 *	Default implementation is empty.*/
		virtual void refinement_step_ends() {};
		
		void adjust_initial_selection();
		void select_closure(std::vector<Vertex*>& vVrts);
		void select_copy_elements(std::vector<Vertex*>& vVrts,
								  int iFirst = 0, int copyRange = -1);
		
		inline void set_status(Vertex* e, StatusMark mark) {m_aaIntVRT[e] = (m_aaIntVRT[e] & ~MR_STATUS) | mark;}
		inline void set_status(Edge* e, StatusMark mark) {m_aaIntEDGE[e] = (m_aaIntEDGE[e] & ~MR_STATUS) | mark;}
		inline void set_status(Face* e, StatusMark mark) {m_aaIntFACE[e] = (m_aaIntFACE[e] & ~MR_STATUS) | mark;}
		inline void set_status(Volume* e, StatusMark mark) {m_aaIntVOL[e] = (m_aaIntVOL[e] & ~MR_STATUS) | mark;}

		virtual void set_rule(Vertex* e, RefinementMark mark) {m_aaIntVRT[e] = (m_aaIntVRT[e] & ~MR_REFINEMENT) | mark;}
		virtual void set_rule(Edge* e, RefinementMark mark) {m_aaIntEDGE[e] = (m_aaIntEDGE[e] & ~MR_REFINEMENT) | mark;}
		virtual void set_rule(Face* e, RefinementMark mark) {m_aaIntFACE[e] = (m_aaIntFACE[e] & ~MR_REFINEMENT) | mark;}
		virtual void set_rule(Volume* e, RefinementMark mark) {m_aaIntVOL[e] = (m_aaIntVOL[e] & ~MR_REFINEMENT) | mark;}

		inline int get_rule(Vertex* e) {return m_aaIntVRT[e] & MR_REFINEMENT;}
		inline int get_rule(Edge* e) {return m_aaIntEDGE[e] & MR_REFINEMENT;}
		inline int get_rule(Face* e) {return m_aaIntFACE[e] & MR_REFINEMENT;}
		inline int get_rule(Volume* e) {return m_aaIntVOL[e] & MR_REFINEMENT;}
/*
		template <typename TIterator>
		void mark_fixed_elements(TIterator iterBegin, TIterator iterEnd);
*/
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

/// @}
}//	end of namespace

#endif
