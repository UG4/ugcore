/*
 * A manager for the degenerated layer subsets. It helps to distinguish between the sides
 * of the degenerated elements and to find out the correspondence of the nodes
 * in them.
 *
 * Created on: Mar. 17, 2014
 * Author: D. Logashenko
 */
#ifndef __H__UG__PLUGINS__D3F__DEGENERATED_LAYER_MANAGER__
#define __H__UG__PLUGINS__D3F__DEGENERATED_LAYER_MANAGER__

#include <map>
#include <vector>

// ug4 headers
#include "common/common.h"
#include "common/util/smart_pointer.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_grid/tools/subset_handler_multi_grid.h"
#include "lib_grid/algorithms/subset_dim_util.h"
#include "lib_grid/grid_objects/grid_dim_traits.h"
#include "lib_grid/algorithms/refinement/global_fractured_media_refiner.h"

namespace ug {

/// Gegenerated layer subset manager
/**
 * Class for the manager of the degenerated layer (e.g. fracture) subsets. Note
 * that the object of the class should get ALL the subsets that belong to ALL
 * the degenerated layers in the domain. This object gets updated every time
 * the grid is refined.
 *
 * Usage instructions:
 * <ul>
 *  <li> Create the object of the DegeneratedLayerManager class </li>
 *  <li> Use the method DegeneratedLayerManager::add to register a degenerated
 *       layer subsets in the object </li>
 *  <li> Call the DegeneratedLayerManager::close method to mark the vertices
 *       as 'inner' or 'outer'. </li>
 * </ul>
 *
 * Remarks:
 * <ul>
 *  <li> The grid must be created (loaded) before the object is created. </li>
 *  <li> All the degenerated layer subsets must be registered in the same degenerated
 *       layer manager. There should exist only one such the manager at all. Alternatively,
 *       the different managers must hold groups of the layers that are not intersect or
 *       connected to each other. But every degenerated subset must be registered
 *       in one of the degenerated layer managers. </li>
 * </ul>
 *
 * References:
 * <ul>
 *  <li> S. Reiter, D. Logashenko, A. Grillo, G. Wittum, Preparation of grids
 *       for simulations of groundwater flow in fractured porous media,
 *       Computing and Visualization in Science, Vol. 15, No. 4 (2012),
 *       pp. 209-225, DOI: 10.1007/s00791-013-0210-7
 *  </li>
 * </ul>
 *
 * \tparam dim	(topological) dimensionality of the grid
 */
template <int dim>
class DegeneratedLayerManager
{
public:
	/// type of the attachment for the marks
		typedef Attachment<signed char> mark_attachment_type;
		
	///	base grid element object type
		typedef typename grid_dim_traits<dim>::grid_base_object element_type;
		
	///	grid element's side base object type
		typedef typename grid_dim_traits<dim>::side_type side_type;
		
	///	max. number of corners of the elements
		static const size_t maxElemCorners = grid_dim_traits<dim>::MaxNumVerticesOfElem;
		
	///	max. number of corners of non-degenerated sides
		static const size_t maxLayerSideCorners = maxElemCorners / 2;
	
	/// Marks for the grid vertices
		enum t_grid_object_mark
		{
			D_LAYER_UNKNOWN = -1,
			D_LAYER_OUTER = 0,
			D_LAYER_INNER = 1
		};

public:
	/// Constructor
		DegeneratedLayerManager
		(
			SmartPtr<MultiGridSubsetHandler> spSH ///< [in] subset handler of the grid
		);
		
	/// Destructor
		virtual ~DegeneratedLayerManager ();
		
	///	Adds a fracture subdomain
		void add
		(
			const char * ss_names ///< [in] subset names of the fractures
		);
		
	///	Removes a fracture subdomain (e.g. for dimension-adaptive method)
		void remove
		(
			const char * ss_names ///< [in] subset names of the fractures
		);

	///	Closes the manager, i.e. computes all the data, ...
		void close ();
		
	///	Initializes a refiner with the fracture subsets
		void init_refiner
		(
			SmartPtr<GlobalFracturedMediaRefiner> refiner, ///< the refiner
			bool as_low_dim ///< whether it should consider the fractures as low-dimentional
		);
		
	///	Returns true if the manager is closed (and can be used) or false otherwise
		bool is_closed () {return m_bClosed;};
		
	///	Whether a subset is registered in the manager
		bool contains
		(
			int si ///< [in] subset index
		)
		{return m_layerSsGrp.contains (si);};
		
	///	Returns the subset group of the fracture network
		const SubsetGroup & subset_grp () {return m_layerSsGrp;};
		
	///	Returs the mark of a vertex
		int vert_mark (Vertex * vrt) {return m_aaVertMarks [vrt];};
		
	///	Gets the inner and the outer fracture sides of an element
		void get_layer_sides
		(
			element_type * elem, ///< [in] the element
			size_t & num_fract_co, ///< [out] number of corners of the inner/outer sides
			side_type * & inner_side, ///< [out] its fracture inner side
			size_t & inner_side_idx, ///< [out] index of the inner side in the reference element
			size_t inner_side_corners [], ///< [out] inner side corner idx -> elem. corner idx
			side_type * & outer_side, ///< [out] its fracture outer side
			size_t & outer_side_idx, ///< [out] index of the outer side in the reference element
			size_t outer_side_corners [], ///< [out] outer side corner idx -> elem. corner idx
			size_t ass_co [] ///< [out] correspondence of the corners of the sides
		);
		
protected:
	/// Marks the inner fracture vertices
		void mark_vertices ();
	
	///	Called when a grid adaption has been performed
		void grid_adaption_callback (const GridMessage_Adaption& msg);

	///	Called when a grid has been distributed between different processes
		void grid_distribution_callback (const GridMessage_Distribution& msg);
	
private:
	/// Subset handler to use
		SmartPtr<MultiGridSubsetHandler> m_spSH;
	
	///	Subset group of the fractures
		SubsetGroup m_layerSsGrp;
	
	/// Attachment keeping the grid object marks for the vertices
		mark_attachment_type m_aVertexMarks;
	///	Attachment accessor
		MultiGrid::AttachmentAccessor<Vertex, mark_attachment_type> m_aaVertMarks;
		
	///	'closed'-flag
		bool m_bClosed;
		
	//	Message hub callback id's for the notifications of the changes of the grid:
		MessageHub::SPCallbackId m_spGridAdaptionCallbackID;
		MessageHub::SPCallbackId m_spGridDistributionCallbackID;
};

} // end namespace ug

#include "deg_layer_mngr_impl.h"

#endif // __H__UG__PLUGINS__D3F__DEGENERATED_LAYER_MANAGER__

/* End of File */
