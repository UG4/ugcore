/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__REFINER_INTERFACE__
#define __H__UG__REFINER_INTERFACE__

#include <string>
#include "projectors/refinement_projector.h"

namespace ug
{
///	\addtogroup lib_grid_algorithms_refinement
///	@{

///	refinement-marks allow to specify how an element shall be processed during refinement.
//	Make sure not to use refinement marks with a value of 128 or higher! Those

///< Fully refines an element and all associated sides and edges
enum RefinementMark{
	RM_NONE = 0,					///< no refinement is performed
	RM_CLOSURE = 1,					///< Refines elements according to associated marked edges
	RM_COPY = RM_CLOSURE,			///< DEPRECATED. Use RM_CLOSURE or RM_LOCAL with localMark = 0 instead.
	RM_ANISOTROPIC = RM_CLOSURE,	///< DEPRECATED. Use RM_CLOSURE instead.
	RM_LOCAL = 1 << 1,				///< Refines elements according to a local refinement mark (use 'mark_local')
	RM_FULL = 1 << 2,				///< Fully refines an element and all associated sides and edges
	RM_REFINE = RM_FULL,			///< DEPRECATED. Use RM_FULL instead.
	RM_COARSEN = 1 << 3,			///< the element is coarsened (only valid for adaptive multi-grid refinement)
	RM_DUMMY = 1 << 4,				///< used internally during mark-adjustment
	RM_MAX							///< the highest constant in RefinementMark. Should always be smaller than 128!
};

///	The refiner interface allows to mark elements for refinement and to call refine.
/**	A refiner always operates on a grid. A grid thus has to be assigned
 * before refinement starts. Please take a look at the specializations
 * of IRefiner, for more information.
 */
class IRefiner
{
	public:
		IRefiner(SPRefinementProjector projector = SPNULL) :
			m_msgIdAdaption(-1), m_projector(projector),
			m_adaptionIsActive(false), m_debuggingEnabled(false)	{}

		virtual ~IRefiner()	{}

		void set_projector(SPRefinementProjector projector)
			{m_projector = projector;}

		SPRefinementProjector projector()
			{return m_projector;}

	///	DEPRECIATED! Use grid(). Has to return the associated grid. Pure virtual
		virtual Grid* get_associated_grid() = 0;
	///	Returns the grid associated with the refiner
	/**	Pure virtual. Specify this method in derived classes!*/
		virtual Grid* grid() = 0;

	///	clears all marks. Default implementation is empty
		virtual void clear_marks()	{}

	///	returns whether the refiner is able to perform adaptive refinement
	/**	pure virtual!*/
		virtual bool adaptivity_supported() const = 0;

	///	returns true, if the refiner supports coarsening.
	/**	pure virtual!*/
		virtual bool coarsening_supported() const = 0;

	///	returns true, if the refiner supports local marks.
		virtual bool local_marks_supported() const 	{return false;}

	///	Marks an element for refinement. Default implementation is empty
	/**	\{ */
		virtual bool mark(Vertex* v, RefinementMark refMark = RM_REFINE)	{return false;}
		virtual bool mark(Edge* e, RefinementMark refMark = RM_REFINE)		{return false;}
		virtual bool mark(Face* f, RefinementMark refMark = RM_REFINE)		{return false;}
		virtual bool mark(Volume* v, RefinementMark refMark = RM_REFINE)	{return false;}
	/**	\} */

	///	marks the specified geometric object
	/**	The default implementation casts the object to a more concrete type
	 * (Vertex, Edge, Face, Volume) and calls the appropriate mark method.*/
		virtual bool mark(GridObject* o, RefinementMark refMark = RM_REFINE);


	///	Marks a face or volume for local refinement.
	/**	The passed mark is an or combination. If the i-th edge of the element
	 * shall be refined, it should hold true: 'mark & 1<<i != 0'.
	 * The passed element will also receive the RM_LOCAL flag.
	 * \note	local-marks differ from mark(e, RM_CLOSURE). The former will
	 *			refine an element according to the marks of associated edges.
	 *			Elements marked with mark_local, however, will only be refined
	 *			according to their local mark.
	 * \{ */
		virtual void mark_local(Face* e, int mark)		{UG_THROW("mark_local not supported by this refiner!");}
		virtual void mark_local(Volume* e, int mark)	{UG_THROW("mark_local not supported by this refiner!");}
	/** \} */

	///	returns the local mark of the specified face or volume.
	/** If the i-th edge of the element shall be refined, it holds true:
	 * 'get_local_mark(e) & 1<<i != 0'
	 * \{ */
		virtual int get_local_mark(Face* e) const	{return 0;}
		virtual int get_local_mark(Volume* e) const	{return 0;}
	/** \} */


	///	marks the neighborhood of currently marked elements.
	/**	In each step direct neighbors of currently marked elements are selected.
	 * The number of iterations thus specifies the width of the neighborhood which
	 * will be marked.
	 * Calls mark_neighborhood(numIterations, RM_NONE, false)*/
	 	void mark_neighborhood(size_t numIterations)
	 	{mark_neighborhood(numIterations, RM_NONE, false);}

	///	marks the neighborhood of currently marked elements.
	/**	In each step direct neighbors of currently marked elements are also marked.
	 * You may specify the refinement mark that will be applied to newly mared elements
	 * - elements which already were marked will be ignored.
	 * By passing RM_NONE as refMark, the refinement-mark will be derived from
	 * neighbored elements.
	 * If sideNbrsOnly is set to true, only elements which are connected to
	 * sides of marked elements are also marked. Otherwise all elements which are
	 * connected to vertices of marked elements are marked.*/
		virtual void mark_neighborhood(
						size_t numIterations,
						RefinementMark refMark,
						bool sideNbrsOnly)			{}

	///	Returns the mark of a given element. Default returns RM_REFINE
	/**	\{ */
		virtual RefinementMark get_mark(Vertex* v)	{return RM_REFINE;}
		virtual RefinementMark get_mark(Edge* e)	{return RM_REFINE;}
		virtual RefinementMark get_mark(Face* f)		{return RM_REFINE;}
		virtual RefinementMark get_mark(Volume* v)		{return RM_REFINE;}
	/**	\} */

	///	returns the mark of the specified geometric object
	/**	The default implementation casts the object to a more concrete type
	 * (Vertex, Edge, Face, Volume) and calls the appropriate get_mark method.*/
		virtual RefinementMark get_mark(GridObject* o);

	///	marks all elements between iterBegin and iterEnd.
	/**	the value-type of TIterator has to be a pointer to a type derived
	 * 	from either Edge, Face or Volume.*/
		template <class TIterator>
		void mark(const TIterator& iterBegin, const TIterator& iterEnd,
				  RefinementMark refMark = RM_REFINE)
			{
				TIterator iter = iterBegin;
				while(iter != iterEnd){
					mark(*iter, refMark);
					++iter;
				}
			}


	///	notifies all listeners of the associated message-hub, that adaption begins / ends.
	/**	While this message is not important to the refiner itself, it may be important
	 * to listeners of the associated grid's message-hub.
	 * \{ */
		void adaption_begins();
		void adaption_ends();
	/**	\} */

	/// Performs refinement on the marked elements.
	/**	internally calls the virtual method 'perform_refinement'*/
		void refine();

	///	Performs coarsening on the elements marked RM_COARSEN.
	/**	Note that coarsening is not supported by all refiners. Normally only
	 * MultiGrid-Refiner do support coarsening.
	 *
	 * coarsen returns false, if no elements have been coarsened, true if at
	 * least one has been coarsened.
	 *
	 * Internally calls the virtual method 'perform_coarsening'
	 */
		bool coarsen();


	///	returns the number of (globally) marked edges on all levels of the hierarchy
		size_t num_marked_edges(std::vector<int>& numMarkedEdgesOut);
	///	returns the number of (globally) marked faces on all levels of the hierarchy
		size_t num_marked_faces(std::vector<int>& numMarkedFacesOut);
	///	returns the number of (globally) marked volumes on all levels of the hierarchy
		size_t num_marked_volumes(std::vector<int>& numMarkedVolsOut);
	///	returns the number of (globally) marked grid-objects of highest dimension
		size_t num_marked_elements(std::vector<int>& numMarkedElemsOut);

	///	returns the number of (globally) marked edges on all levels of the hierarchy
		size_t num_marked_edges()		{std::vector<int> t; return num_marked_edges(t);}
	///	returns the number of (globally) marked faces on all levels of the hierarchy
		size_t num_marked_faces()		{std::vector<int> t; return num_marked_faces(t);}
	///	returns the number of (globally) marked volumes on all levels of the hierarchy
		size_t num_marked_volumes()		{std::vector<int> t; return num_marked_volumes(t);}
	///	returns the number of (globally) marked grid-objects of highest dimension
		size_t num_marked_elements()	{std::vector<int> t; return num_marked_elements(t);}

	///	Writes the associated grid and marks to a file. Pure virtual.
	/**	Elements should be assigned to subsets depending on their current
	 * refinement-mark.*/
		virtual bool save_marks_to_file(const char* filename) = 0;

	///	sets a filename to which adjusted marks are saved during refinement / coarsening
	/**	If no filename is set, then no marks are being saved during refinement / coarsening.
	 * If you want to unset the file, either pass a NULL pointer or an empty string.*/
		void set_adjusted_marks_debug_filename(const char* filename);

		void enable_debugging(bool enable)	{m_debuggingEnabled = enable;}
		bool debugging_enabled() const		{return m_debuggingEnabled;}

	protected:
	///	sets the message hub.
	/**	A message hub is required, since it is used transmit messages regarding
	 * adaption, refinement and coarsening.*/
		void set_message_hub(SPMessageHub msgHub);

	///	called by refine(). Derived classes should implement their refinement algorithm here.
		virtual void perform_refinement() = 0;

	///	Called by coarsen(). Derived classes sould implement their coarsen algorithm here.
	/** Since the default implementation does not perform coarsening, it returns false.*/
		virtual bool perform_coarsening()		{return false;}

	///	returns the number of locally marked edges on all levels of the hierarchy
		virtual void num_marked_edges_local(std::vector<int>& numMarkedEdgesOut) = 0;
	///	returns the number of locally marked faces on all levels of the hierarchy
		virtual void num_marked_faces_local(std::vector<int>& numMarkedFacesOut) = 0;
	///	returns the number of locally marked volumes on all levels of the hierarchy
		virtual void num_marked_volumes_local(std::vector<int>& numMarkedVolsOut) = 0;

	protected:
		SPMessageHub			m_messageHub;
		int						m_msgIdAdaption;
		SPRefinementProjector	m_projector;
		bool					m_adaptionIsActive;
		bool					m_debuggingEnabled;
		std::string				m_adjustedMarksDebugFilename;
};

/// @}	// end of add_to_group command

}//	end of namespace

#endif
