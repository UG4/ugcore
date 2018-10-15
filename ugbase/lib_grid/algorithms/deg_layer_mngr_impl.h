/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Dmitry Logashenko
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

/*
 * Implementation of the degenerated layer subset manager.
 */
// ug4 headers
#include "common/util/string_util.h"
#include "lib_grid/grid/grid.h"

namespace ug {

/**
 * Class constructor: Parses the subset names, attaches itself to the grid
 * to be updated if the grid is refined, and writes itself to the properties
 * of the subsets.
 */
template <int dim>
DegeneratedLayerManager<dim>::DegeneratedLayerManager
(
	SmartPtr<MultiGridSubsetHandler> spSH ///< [in] subset handler of the grid
)
:	m_spSH (spSH), m_bClosed (false)
{
// Check the subset handler:
	if (m_spSH.invalid ())
		UG_THROW ("DegeneratedLayerManager: Invalid subset handler specified!");
	
//	Initialize the data:
	m_layerSsGrp.set_subset_handler (m_spSH);
	
//	Initialize the attachment:
	MultiGrid * pMG = (MultiGrid *) (m_spSH->multi_grid ());
	pMG->template attach_to_dv<Vertex> (m_aVertexMarks, D_LAYER_UNKNOWN);
	m_aaVertMarks.access (*pMG, m_aVertexMarks);
	
//	Register the callbacks in the message hub:
	m_spGridAdaptionCallbackID =
		pMG->message_hub()->register_class_callback
			(this, & DegeneratedLayerManager<dim>::grid_adaption_callback);
	m_spGridDistributionCallbackID =
		pMG->message_hub()->register_class_callback
			(this, & DegeneratedLayerManager<dim>::grid_distribution_callback);
}

/**
 * Class destructor: Detaches the attachment.
 */
template <int dim>
DegeneratedLayerManager<dim>::~DegeneratedLayerManager ()
{
	MultiGrid * pMG = (MultiGrid *) (m_spSH->multi_grid ());
	m_aaVertMarks.invalidate ();
	pMG->template detach_from<Vertex> (m_aVertexMarks);
	m_bClosed = false; // (to be on the very safe side...)
}

/**
 * Adds subsets to the manager
 */
template <int dim>
void DegeneratedLayerManager<dim>::add
(
	const char * ss_names ///< [in] subset names of the degenerated layers
)
{
	if (m_bClosed)
		UG_LOG ("DegeneratedLayerManager: Adding subsets to a closed manager.");
	
//	If extended - then not closed any more
	m_bClosed = false;
	
//	Parse the subset names:
	std::vector<std::string> vssNames;
	TokenizeString (ss_names, vssNames);
	
//	Add the subsets to the group:
	m_layerSsGrp.add (vssNames);
}

/**
 * Removes subsets from the manager
 */
template <int dim>
void DegeneratedLayerManager<dim>::remove
(
	const char * ss_names ///< [in] subset names of the fractures
)
{
//	If extended - then not closed any more
	m_bClosed = false;

//	Parse the subset names:
	std::vector<std::string> vssNames;
	TokenizeString (ss_names, vssNames);

//	Add the subsets to the group:
	m_layerSsGrp.remove (vssNames);
}

/**
 * Closes the manager: Computes the vertex marks
 */
template <int dim>
void DegeneratedLayerManager<dim>::close ()
{
	size_t init_n_subsets = m_layerSsGrp.size ();
	
//	Check the dimensionality of the subsets:
	size_t i = 0;
	while (i < m_layerSsGrp.size ())
	if (DimensionOfSubset (*m_spSH, m_layerSsGrp [i]) != dim)
	{
		UG_LOG ("DegeneratedLayerManager: Low-dimensional subset '"
			<< m_spSH->get_subset_name (m_layerSsGrp [i])
				<< "' specified as a degenerated layer. Removed from the list.\n");
		m_layerSsGrp.remove (m_layerSsGrp [i]);
		i = 0;
	}
	else i++;
	
//	Check if all the subsets are removed due to the illegal dimensionality
	if (init_n_subsets != 0 && m_layerSsGrp.size () == 0)
		UG_THROW ("DegeneratedLayerManager: All the specified subsets have illegal dimensionality!\n");
	
//	Mark the vertices:
	mark_vertices ();
	
//	Done:
	m_bClosed = true;
}

/**
 * Adds the registered subsets to the refiner
 */
template <int dim>
void DegeneratedLayerManager<dim>::init_refiner
(
	SmartPtr<GlobalFracturedMediaRefiner> refiner,
	bool as_low_dim
)
{
	if (! is_closed ())
		UG_THROW ("DegeneratedLayerManager: The manager must be closed before use.");
	for (size_t i = 0; i <  m_layerSsGrp.size (); i++)
		refiner->mark_as_fracture (m_layerSsGrp[i], as_low_dim);
}

/**
 * Marks the vertices in such a way that the inner layer vertices are
 * marked with D_LAYER_INNER, whereas all other vertices (at the fractures
 * and in the bulk medium) are marked with D_LAYER_OUTER. After the call,
 * there should be no vertices marked with D_LAYER_UNKNOWN.
 */
template <int dim>
void DegeneratedLayerManager<dim>::mark_vertices ()
{
	typedef typename geometry_traits<element_type>::const_iterator t_elem_iter;
	
	MultiGrid * pMG = (MultiGrid *) (m_spSH->multi_grid ());
	
//	Mark all vertices of the fracture elements:
	for (size_t i = 0; i < m_layerSsGrp.size (); i++)
	{
		int si = m_layerSsGrp [i];
		for (size_t lev = 0; lev < pMG->num_levels (); lev++)
		{
			t_elem_iter list_end = m_spSH->template end<element_type> (si, lev);
			for (t_elem_iter iElem = m_spSH->template begin<element_type> (si, lev);
					iElem != list_end; ++iElem)
			{
				element_type * elem = *iElem;
				for (size_t i = 0; i < elem->num_vertices (); i++)
				{
					Vertex * vert = elem->vertex(i);
					m_aaVertMarks [vert] = (! vert->is_constrained ())? D_LAYER_INNER : D_LAYER_OUTER;
					// REMARK: Constrained vertices cannot be inner for degenerated layers!
					// We cannot catch this situation somewhere else, so we do it here.
				}
			}
		}
	}
	
//	Unmark all vertices in the other (not-degenerated-layer) elements:
	SubsetGroup bulkSsGrp (m_spSH);
	for (int si = 0; si < m_spSH->num_subsets (); si++)
	if (DimensionOfSubset (*m_spSH, si) == dim && ! m_layerSsGrp.contains (si))
		for (size_t lev = 0; lev < pMG->num_levels (); lev++)
		{
			t_elem_iter list_end = m_spSH->template end<element_type> (si, lev);
			for (t_elem_iter iElem = m_spSH->template begin<element_type> (si, lev);
				iElem != list_end; ++iElem)
			{
				element_type * elem = *iElem;
				for (size_t i = 0; i < elem->num_vertices (); i++)
					m_aaVertMarks [elem->vertex(i)] = D_LAYER_OUTER;
			}
		}
}

/**
 * The grid adaption callback: Catches the messages about grid refinements, ...
 */
template <int dim>
void DegeneratedLayerManager<dim>::grid_adaption_callback
(
	const GridMessage_Adaption & msg ///< the message from the message hub
)
{
	if (m_bClosed && msg.adaption_ends ())
		mark_vertices ();
}

/**
 * The grid adaption callback: Catches the messages about the distribution of the grid
 */
template <int dim>
void DegeneratedLayerManager<dim>::grid_distribution_callback
(
	const GridMessage_Distribution & msg ///< the message from the message hub
)
{
	if (m_bClosed && msg.msg () == GMDT_DISTRIBUTION_STOPS)
		mark_vertices ();
}

/**
 * For a given element, finds its inner and outer layer sides. For these
 * sides, gets the correspondence of the vertices in the degenerated element
 * (i.e. considering all other sides as degenerated).
 */
template <int dim>
void DegeneratedLayerManager<dim>::get_layer_sides
(
	element_type * elem, ///< [in] the element
	size_t & num_side_co, ///< [out] number of corners of the inner/outer sides
	side_type * & inner_side, ///< [out] its fracture inner side
	size_t & inner_side_idx, ///< [out] index of the inner side in the reference element
	size_t inner_side_corners [], ///< [out] inner side corner idx -> elem. corner idx (maxLayerSideCorners elements)
	side_type * & outer_side, ///< [out] its fracture outer side
	size_t & outer_side_idx, ///< [out] index of the outer side in the reference element
	size_t outer_side_corners [], ///< [out] outer side corner idx -> elem. corner idx (maxLayerSideCorners elements)
	size_t ass_co [] ///< [out] correspondence of the corners of the sides (2 * maxLayerSideCorners elements or NULL)
)
{
	size_t n_inner, n_outer;
	size_t n_co = elem->num_vertices ();
	MultiGrid * pMG = (MultiGrid *) (m_spSH->multi_grid ());
	typename Grid::traits<Edge>::secure_container edge_list;
	
//	First, we check, whether the inner and outer sides are well defined.
//	Loop over the sides: We look for sides with only inner or only outer corners
	inner_side = NULL; n_inner = 0;
	outer_side = NULL; n_outer = 0;
	for (size_t i = 0; i < elem->num_sides (); i++)
	{
		side_type * side = pMG->get_side (elem, i);
		bool has_inner = false, has_outer = false;
		size_t side_n_co = side->num_vertices ();
		for (size_t co = 0; co < side_n_co; co++)
		{
			int mark = vert_mark (side->vertex (co));
			if (mark == D_LAYER_OUTER)
				has_outer = true;
			else if (mark == D_LAYER_INNER)
				has_inner = true;
			else
				UG_THROW ("DegeneratedLayerManager: Some vertices are not marked!");
		}
		
		if (has_inner && has_outer) continue; // this is a 'degenerated' side
		
		if (has_inner)
		{	// this is the inner side
			if (inner_side != NULL)
				UG_THROW ("DegeneratedLayerManager: Two inner sides of a degenerated layer element found!");
			inner_side = side;
			inner_side_idx = i;
			n_inner = side_n_co;
		}
		else
		{	// this is the outer side
			if (outer_side != NULL)
				UG_THROW ("DegeneratedLayerManager: Two outer sides of a degenerated layer element found!");
			outer_side = side;
			outer_side_idx = i;
			n_outer = side_n_co;
		}
	}
	if (n_inner != n_outer || 2 * n_inner != n_co)
		UG_THROW ("DegeneratedLayerManager: Unrecognized degenerated layer element!");
	num_side_co = n_inner;
	
//	Now: Corner of inner side are exactly those marked as inner, and the same for the outer side
	Vertex * inner_vert_ptr [maxLayerSideCorners], * outer_vert_ptr [maxLayerSideCorners];
	size_t inner_corners [maxLayerSideCorners], outer_corners [maxLayerSideCorners];
	size_t inner_i = 0, outer_i = 0;
	for (size_t co = 0; co < n_co; co++)
	{
		Vertex * vrt = elem->vertex (co);
		if (vert_mark (vrt) == D_LAYER_INNER)
		{
			inner_vert_ptr [inner_i] = vrt;
			inner_corners [inner_i] = co;
			inner_i++;
		}
		else
		{
			outer_vert_ptr [outer_i] = vrt;
			outer_corners [outer_i] = co;
			outer_i++;
		}
	}
	if (inner_i != n_inner || outer_i != n_outer)
		UG_THROW ("DegeneratedLayerManager: Failed to get the vertex-side correspondence!");
	
//	Order the corner indices according to the corners of the sides
	for (size_t i = 0; i < n_inner; i++)
	{
		size_t co;
		
		co = 0;
		while (inner_side->vertex (co) != inner_vert_ptr [i])
			if (++co == n_inner)
				UG_THROW ("DegeneratedLayerManager: Failed! Missing inner vertex?");
		inner_side_corners [co] = inner_corners [i];
		
		co = 0;
		while (outer_side->vertex (co) != outer_vert_ptr [i])
			if (++co == n_outer)
				UG_THROW ("DegeneratedLayerManager: Failed! Missing outer vertex?");
		outer_side_corners [co] = outer_corners [i];
	}
	
//	Loop over the edges: Find out the associated corners
	if (ass_co != NULL)
	{
		pMG->associated_elements (edge_list, elem);
		for (size_t i = 0; i < edge_list.size (); i++)
		{
			Edge * e = edge_list [i];
			Vertex * v_1 = e->vertex (0); int mark_1 = vert_mark (v_1);
			Vertex * v_2 = e->vertex (1); int mark_2 = vert_mark (v_2);
			Vertex * t;
			size_t co_1, co_2;
			
			if (mark_1 == mark_2)
				continue; // an edge of the inner or the outer size
			
			if (mark_1 == D_LAYER_OUTER)
			{	// Vertex # 1 must be inner, vertex # 2 must be outer
				t = v_1; v_1 = v_2; v_2 = t;
			}
		
			co_1 = 0;
			while (inner_vert_ptr[co_1] != v_1)
				if (++co_1 == n_inner)
					UG_THROW ("DegeneratedLayerManager: Cannot find an inner node by vertex!");
			
			co_2 = 0;
			while (outer_vert_ptr[co_2] != v_2)
				if (++co_2 == n_outer)
					UG_THROW ("DegeneratedLayerManager: Cannot find an outer node by vertex!");
			
			ass_co [inner_corners [co_1]] = outer_corners [co_2];
			ass_co [outer_corners [co_2]] = inner_corners [co_1];
		}
	}
}

/**
 * Assigns a different subset index to the inner sides of a layer. (Only to
 * the sides, not vertices, or edges in 3d!) If -1 passed instead of the index
 * of the subset to assign, new subset is created. The function returns the
 * index of the assigned subset.
 */
template <int dim>
int DegeneratedLayerManager<dim>::assign_middle_subset
(
	int layer_si, ///< subset index of the layer
	int middle_si ///< the subset index to assign (or -1 to create a new subset)
)
{
	typedef typename geometry_traits<element_type>::const_iterator t_elem_iter;
	
	MultiGrid * pMG = (MultiGrid *) (m_spSH->multi_grid ());
	
//	Assign the default value of the subset (if necessary)
	if (middle_si < 0) middle_si = m_spSH->num_subsets ();
	
//	Loop the elements of the layer in all the grid levels
	pMG->message_hub()->post_message (GridMessage_Creation (GMCT_CREATION_STARTS, 0));
	for (size_t lev = 0; lev < pMG->num_levels (); lev++)
	{
		t_elem_iter list_end = m_spSH->template end<element_type> (layer_si, lev);
		for (t_elem_iter iElem = m_spSH->template begin<element_type> (layer_si, lev);
				iElem != list_end; ++iElem)
		{
			element_type * elem = *iElem;
			size_t num_side_co;
			side_type * inner_side, * outer_side;
			size_t inner_side_idx, outer_side_idx;
			size_t inner_side_corners [maxLayerSideCorners], outer_side_corners [maxLayerSideCorners];
			
		//	get the inner side
			get_layer_sides (elem, num_side_co,
				inner_side, inner_side_idx, inner_side_corners,
				outer_side, outer_side_idx, outer_side_corners);
			
		//	assign the subset
			m_spSH->assign_subset (inner_side, middle_si);
		}
	}
	pMG->message_hub()->post_message (GridMessage_Creation (GMCT_CREATION_STOPS, 0));
	
	return middle_si;
}

/**
 * Assigns a different subset to the inner sides of a layer. (Only to
 * the sides, not vertices, or edges in 3d!) If the given subset does not exist,
 * it is created. The function returns the index of the assigned subset.
 */
template <int dim>
int DegeneratedLayerManager<dim>::assign_middle_subset
(
	int layer_si, ///< subset index of the layer
	const char* middle_ss_name ///< name of the subset to assign
)
{
//	Look for the subset to assign
	int middle_si = -1;
	bool new_ss = true;
	for (int si = 0; si < m_spSH->num_subsets (); si++)
	if (strcmp (middle_ss_name, m_spSH->get_subset_name (si)) == 0)
	{
		middle_si = si;
		new_ss = false;
		break;
	}
	
//	Assign the subset
	middle_si = assign_middle_subset (layer_si, middle_si);
	
//	Set the name
	if (new_ss)
		m_spSH->set_subset_name (middle_ss_name, middle_si);
	
	return middle_si;
}

/**
 * Assigns a different subset to the inner sides of a layer. (Only to
 * the sides, not vertices, or edges in 3d!) If the given subset does not exist,
 * it is created. The function returns the index of the assigned subset.
 */
template <int dim>
int DegeneratedLayerManager<dim>::assign_middle_subset
(
	const char* layer_ss_name, ///< subset name of the layer
	const char* middle_ss_name ///< name of the subset to assign
)
{
//	Look for the subset of the layer
	for (int si = 0; si < m_spSH->num_subsets (); si++)
	if (strcmp (layer_ss_name, m_spSH->get_subset_name (si)) == 0)
		return assign_middle_subset (si, middle_ss_name);
	UG_THROW ("DegeneratedLayerManager: Subset '" << layer_ss_name << "' not found");
}

} // end namespace ug

/* End of File */
