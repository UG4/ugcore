/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__DOMAIN__
#define __H__UG__LIB_DISC__DOMAIN__

#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/refinement/projectors/refinement_projector.h"

#include <map>

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#endif

namespace ug{

/**
 * Domain
 *
 * \defgroup lib_disc_domain Domain
 * \ingroup lib_discretization
 * \{
 */

///	Describes the contents of a domain.
/**	In a parallel environment those are the contents of the global distributed domain.
*/
class DomainInfo
{
	public:
		inline int element_type()	const							{return m_elementType;}
		inline size_t num_levels() const							{return m_numElems.size();}
	///	returns the global number of elements on the given level (excluding ghosts...)
		inline int num_elements_on_level(size_t lvl) const			{return m_numElems[lvl];}
	///	returns the local number of elements on the given level (excluding ghosts...)
		inline int num_local_elements_on_level(size_t lvl) const	{return m_numLocalElems[lvl];}
	///	returns the minimum number of elements a process has on a given leven (excluding ghosts)
		inline int min_num_local_elements_on_level(size_t lvl) const	{return m_minNumLocalElems[lvl];}
	///	returns the maximum number of elements a process has on a given leven (excluding ghosts)
		inline int max_num_local_elements_on_level(size_t lvl) const	{return m_maxNumLocalElems[lvl];}
	///	returns the local number of ghosts on the given level
		inline int num_local_ghosts_on_level(size_t lvl) const		{return m_numLocalGhosts[lvl];}

		inline int num_subsets() const
			{return (int)m_subsetDims.size();}
		inline int subset_dim(int si) const
			{return m_subsetDims[si];}

		inline int num_elements() const
			{
				int total = 0;
				for(size_t i = 0; i < m_numElems.size(); ++i){
					total += m_numElems[i];
				}
				return total;
			}

		inline void set_info(GridBaseObjectId elemType,
								const std::vector<int>& numElems,
								const std::vector<int>& numLocalElems,
								const std::vector<int>& minNumLocalElems,
								const std::vector<int>& maxNumLocalElems,
								const std::vector<int>& numLocalGhosts,
								const std::vector<int>& subsetDims)
			{m_elementType = elemType;
			 m_numElems = numElems;
			 m_numLocalElems = numLocalElems;
			 m_minNumLocalElems = minNumLocalElems;
			 m_maxNumLocalElems = maxNumLocalElems;
			 m_numLocalGhosts = numLocalGhosts;
			 m_subsetDims = subsetDims;}

		std::string to_string() const;

	private:
		GridBaseObjectId	m_elementType;
		std::vector<int>	m_numElems;
		std::vector<int>	m_numLocalElems;///< local number of elements excluding ghosts.
		std::vector<int>	m_minNumLocalElems;
		std::vector<int>	m_maxNumLocalElems;
		std::vector<int>	m_numLocalGhosts;
		std::vector<int>	m_subsetDims;
};


/// describes a physical domain
/**
 * A Domain collects and exports relevant information about the
 * physical domain, that is intended to be discretized. It will be used as
 * a template Parameter in several classes to distinguish at compile-time
 * between needed types and parameters. It mainly has a grid and a subset
 * handler for the grid to define different physical subsets. In addition
 * to a grid, that may only contain topological informations a domain always
 * has a Position Attachment holding the physical coordinates.
 *
 * \tparam	TGrid			Grid type
 * \tparam	TSubsetHandler	Subset Handler type
 */
template <typename TGrid = MultiGrid, typename TSubsetHandler = MGSubsetHandler>
class IDomain
{
	public:
	///	Grid type
		typedef TGrid grid_type;

	///	Subset Handler type
		typedef TSubsetHandler subset_handler_type;

	public:
	///	Default constructor
	/**
	 * creates an empty domain. Grid and Subset Handler are set up. The
	 * Distributed Grid Manager is set in the parallel case.
	 * \param[in]	options		Grid Options (optinal)*/
		IDomain(bool isAdaptive = true);

	///	Destructor
		virtual ~IDomain();

	///	World Dimension
		virtual int get_dim() const = 0;
		
	///	returns Grid
		inline SmartPtr<TGrid> grid() {return m_spGrid;};

	///	const access to Grid
		inline const ConstSmartPtr<TGrid> grid() const {return m_spGrid;};

	///	returns Subset Handler
		inline SmartPtr<TSubsetHandler> subset_handler() {return m_spSH;};

	///	const access to Subset Handler
		inline const ConstSmartPtr<TSubsetHandler> subset_handler() const {return m_spSH;};

	///	returns the message hub of the grid
		SPMessageHub message_hub() {return m_spGrid->message_hub();}

	///	returns whether the domain may be used for adaptive refinement
		bool is_adaptive() const		{return m_isAdaptive;}

	///	returns whether the associated grid is empty
	/**	Note that one vertex is enough to consider the grid as non-empty.*/
		bool empty() const			{return m_spGrid->num_vertices() == 0;}

	///	updates and broadcasts subset names and dimensions from the given rootProc to all other processes.
		void update_subset_infos(int rootProc);

	///	returns information on the current domain
	/**	In a parallel environment, this information relates to the global
	 * (distributed) domain.*/
		const DomainInfo& domain_info() const		{return m_domainInfo;}

	///	updates the internal domain-info object.
	/**	This method is called automatically each time the associated grid has changed.*/
		void update_domain_info();

	///	returns whether the domain can be used for parallel computations
	/**	If ug was build with support for parallelism, this method will always return true
	 * and always false, if ug was build for serial environments.*/
		bool is_parallel()											{return m_spGrid->is_parallel();}

	///	returns Distributed Grid Manager
		inline DistributedGridManager* distributed_grid_manager()	{return m_spGrid->distributed_grid_manager();}

	///	creates an additional subset-handler with the given name
	/**	If this subset-handler is created before the domain is loaded from a file,
	 * the contents of the handler will be filled with the contents of the
	 * file's subset-handler-node with the same name.*/
		bool create_additional_subset_handler(std::string name);

	///	returns a list with the names of additional subset handlers
		std::vector<std::string> additional_subset_handler_names() const;

	///	returns an additional subset handler Subset Handler
		SmartPtr<TSubsetHandler> additional_subset_handler(std::string name);

	///	const access to Subset Handler
		const ConstSmartPtr<TSubsetHandler> additional_subset_handler(std::string name) const;

	///	sets the ug::RefinementProjector which can be used by refiners during refinement
		void set_refinement_projector(SPRefinementProjector proj);

	///	returns the domain's ug::RefinementProjector. The pointer may be invalid.
		SPRefinementProjector refinement_projector() const;

	///	returns the geometry of the domain
		virtual SPIGeometry3d geometry3d() const = 0;

	protected:
		#ifdef UG_PARALLEL
		/// helper method to broadcast ug::RefinementProjectors to different processes
			SPRefinementProjector
			broadcast_refinement_projector (
					int rootProc,
					pcl::ProcessCommunicator& procCom,
					SPIGeometry3d geometry,
					SPRefinementProjector projector = SPNULL);

			void serialize_projector (
					BinaryBuffer& bufOut,
					SPRefinementProjector proj);

			SPRefinementProjector deserialize_projector (BinaryBuffer& buf);
		#endif
			
		SmartPtr<TGrid> m_spGrid;			///< Grid
		SmartPtr<TSubsetHandler> m_spSH;	///< Subset Handler
		std::map<std::string, SmartPtr<TSubsetHandler> >	m_additionalSH; ///< additional subset handlers

		SPRefinementProjector		m_refinementProjector;

		MessageHub::SPCallbackId 	m_spGridAdaptionCallbackID;
		MessageHub::SPCallbackId 	m_spGridCreationCallbackID;
		MessageHub::SPCallbackId 	m_spGridDistributionCallbackID;

		DomainInfo	m_domainInfo;

		bool	m_isAdaptive;
		bool	m_adaptionIsActive;

	/**	this callback is called by the message hub, when a grid adaption has been
	 * performed. It will call all necessary actions in order to keep the grid
	 * correct for computations. */
		inline void grid_adaption_callback(const GridMessage_Adaption& msg);

	///	Called when a domain has been loaded and during domain distribution
		inline void grid_creation_callback(const GridMessage_Creation& msg);

	/**	this callback is called by the message hub, when a grid has been distributed
	 * between different processes.*/
		inline void grid_distribution_callback(const GridMessage_Distribution& msg);

#ifdef UG_PARALLEL
	protected:
	/**	make sure that elements of the given type are contained in at most one
	 * vmaster interface. This is always the case for highest dimensional elements.*/
		template <class TElem>
		void count_ghosts(std::vector<int>& numGhostsOnLvlOut);
#endif
};


template <int d, typename TGrid = MultiGrid, typename TSubsetHandler = MGSubsetHandler>
class Domain : public IDomain<TGrid, TSubsetHandler>
{
	private:
	/// base type
		typedef IDomain<TGrid, TSubsetHandler> base_type;

	public:
	///	World dimension
		static const int dim = d;

	///	Type of position coordinates
		typedef MathVector<dim> position_type;

	///	Type of Position Attachment
		typedef Attachment<position_type> position_attachment_type;

	///	Type of Accessor to the Position Data Attachment
		typedef Grid::VertexAttachmentAccessor<position_attachment_type>
					position_accessor_type;

	public:
	///	Grid type
		typedef typename base_type::grid_type grid_type;

	///	Subset Handler type
		typedef typename base_type::subset_handler_type subset_handler_type;

	public:
	///	Default constructor
	/**
	 * creates an empty domain. Grid and Subset Handler are set up. The
	 * Distributed Grid Manager is set in the parallel case.
	 * \param[in]	options		Grid Options (optional)
	 */
		Domain(bool isAdaptive = true);
		virtual ~Domain() {};

	///	World Dimension
		virtual int get_dim() const {return dim;}

	///	returns Position Attachment
		inline position_attachment_type& position_attachment(){return m_aPos;}

	///	const access to Position Attachment
		inline const position_attachment_type& position_attachment() const {return m_aPos;}

	///	get Position Accessor
		inline position_accessor_type& position_accessor() {return m_aaPos;}

	///	const access to Position Accessor
		inline const position_accessor_type& position_accessor() const{return m_aaPos;}

		virtual SPIGeometry3d geometry3d() const	{return m_geometry3d;}

	protected:
		position_attachment_type m_aPos;	///<Position Attachment
		position_accessor_type	m_aaPos;		///<Accessor
		SPIGeometry3d			m_geometry3d;
};

typedef Domain<1, MultiGrid, MGSubsetHandler> Domain1d;
typedef Domain<2, MultiGrid, MGSubsetHandler> Domain2d;
typedef Domain<3, MultiGrid, MGSubsetHandler> Domain3d;

} // end namespace ug

// end group lib_disc_domain
/// \}

// include implementation
#include "domain_impl.h"

#endif /* __H__UG__LIB_DISC__DOMAIN__ */
