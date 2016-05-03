/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__GLOBAL_FRACTURED_MEDIA_REFINER__
#define __H__LIB_GRID__GLOBAL_FRACTURED_MEDIA_REFINER__

#include <vector>
#include "lib_grid/lg_base.h"
#include "lib_grid/multi_grid.h"
#include "refinement_projectors_old/refinement_callbacks.h"
#include "hanging_node_refiner_multi_grid.h"
#include "lib_grid/tools/bool_marker.h"


namespace ug
{

class ISubsetHandler;

///	\addtogroup lib_grid_algorithms_refinement
///	@{
///	Performs global grid refinement, while performig anisotropic refinement in the given fractures.
/** If a subset shall be refined as a fracture, you may tell the refiner through
 * its method mark_as_fracture. Make sure that the subset specified in this
 * call, consists of two layers of elements. Pleas have a look at the literature
 * for a more detailed description of a fractures topology.
 * If you don't want to treat a subset as a fracture anymore (i.e. since the elements
 * have a good aspect ratio after a couple of refinements), you may tell the
 * refiner by setting the fracture mark to false using the method
 * mark_as_fracture once again.
 */
//template <class TAPosition>
class GlobalFracturedMediaRefiner : public IRefiner, public GridObserver
{
	public:
		GlobalFracturedMediaRefiner(IRefinementCallback* refCallback = NULL);
		GlobalFracturedMediaRefiner(MultiGrid& mg,
							   	     IRefinementCallback* refCallback = NULL);
							   
		virtual ~GlobalFracturedMediaRefiner();

		virtual void grid_to_be_destroyed(Grid* grid);
		
		void assign_grid(MultiGrid& mg);
		void assign_grid(MultiGrid* mg);

		void set_subset_handler(ISubsetHandler* sh)	{m_pSH = sh;}
		void set_subset_handler(ISubsetHandler& sh)	{m_pSH = &sh;}

	///	sets the position attachment
	/**	If you don't explicitly set the position attachment through this method,
	 * the default position attachment is used, as returned by
	 * ug::GetDefaultPositionAttachment.*/
		//void set_position_attachment(TAPosition& aPos)	{m_aPos = aPos;}

	///	if enabled, a subset will be regarded as a fracture.
	/**	If a subset is regarded as a fracture, it will be refined appropriately.
	 * Please make sure, that the specified subset has a valid topology, as
	 * described in ug::GlobalFracturedMediaRefiner.
	 * @param subInd	The index of the subset whose property is set
	 * @param isFracture	true or false, indicating whether the subset shall
	 * 						be regarded as a fracture or not.*/
		void mark_as_fracture(int subInd, bool isFracture);

	///	returns whether the specified subset is regarded as a fracture.
		bool is_fracture(int subInd);


		virtual Grid* get_associated_grid()		{return m_pMG;}
		virtual Grid* grid()					{return m_pMG;}
		virtual MultiGrid* multi_grid()			{return m_pMG;}

		virtual bool adaptivity_supported() const	{return false;}
		virtual bool coarsening_supported() const	{return false;}

		virtual bool save_marks_to_file(const char* filename);

	protected:
	///	returns the number of (globally) marked edges on this level of the hierarchy
		virtual void num_marked_edges_local(std::vector<int>& numMarkedEdgesOut);
	///	returns the number of (globally) marked faces on this level of the hierarchy
		virtual void num_marked_faces_local(std::vector<int>& numMarkedFacesOut);
	///	returns the number of (globally) marked volumes on this level of the hierarchy
		virtual void num_marked_volumes_local(std::vector<int>& numMarkedVolsOut);

		template <class TElem>
		void num_marked_elems(std::vector<int>& numMarkedElemsOut);
		
	////////////////////////////////
	///	performs refinement on the marked elements.
		virtual void perform_refinement();

	///	called by perform_refinement to adjust the marks
	/**	Everything that shall be refined, should be marked in m_marker.
	 * Note that the method calls communicate_marks twice, to allow derived
	 * classes to e.g. communicate marks in a parallel environment. There thus
	 * shouldn't be a need to reimplement adjust_marks in a derived class.*/
		virtual void adjust_marks();

	///	Called by adjust_marks. Default implementation does nothing.
	/**	If you communicate marks (using an or operation) in this method, then the
	 * GlobalFracturedMediaRefiner should run fine in a parallel environment, too.
	 * The default implementation does nothing (that's fine for a serial environment).*/
		virtual void communicate_marks(BoolMarker& marker)		{}

	///	performs the actual marking
	/**	This class is specialized for Face and Volume.*/
		template <class TElem>
		void assign_elem_and_side_marks();

	///	recursively marks sides of all marked top level elements of the given type
		template <class TElem>
		void mark_sides_of_marked_top_level_elements();

	///	returns the number of marked entries
		template <class TElem>
		size_t num_marked(const std::vector<TElem*>& elems) const;

	///	a callback that allows to deny refinement of special vertices
		virtual bool refinement_is_allowed(Vertex* elem)	{return true;}
	///	a callback that allows to deny refinement of special edges
		virtual bool refinement_is_allowed(Edge* elem)		{return true;}
	///	a callback that allows to deny refinement of special faces
		virtual bool refinement_is_allowed(Face* elem)			{return true;}
	///	a callback that allows to deny refinement of special volumes
		virtual bool refinement_is_allowed(Volume* elem)		{return true;}
		
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
		
	///	returns true if the specified element is a fracture element.
	/**	Note that this method does not check whether the underlying subset-handler
	 * is valid. Make sure to check that beforehand.*/
		template <class TElem>
		bool is_fracture_element(TElem* e)	{return is_fracture(m_pSH->get_subset_index(e));}

	protected:
		BoolMarker			m_marker;
		std::vector<bool>	m_subsetIsFracture;
		//TAPosition			m_aPos;
		MultiGrid*			m_pMG;
		ISubsetHandler*		m_pSH;
};

/// @}
}//	end of namespace

#endif
