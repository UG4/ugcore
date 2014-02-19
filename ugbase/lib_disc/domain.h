/*
 * domain.h
 *
 *  Created on: 17.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOMAIN__
#define __H__UG__LIB_DISC__DOMAIN__

#include "lib_grid/algorithms/subset_util.h"
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
		inline GridBaseObjectId element_type()	const			{return m_elementType;}
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

	protected:
		SmartPtr<TGrid> m_spGrid;			///< Grid
		SmartPtr<TSubsetHandler> m_spSH;	///< Subset Handler
		std::map<std::string, SmartPtr<TSubsetHandler> >	m_additionalSH; ///< additional subset handlers

		MessageHub::SPCallbackId m_spGridAdaptionCallbackID;
		MessageHub::SPCallbackId m_spGridCreationCallbackID;
		MessageHub::SPCallbackId m_spGridDistributionCallbackID;

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
	 * \param[in]	options		Grid Options (optinal)
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

	protected:
		position_attachment_type m_aPos;	///<Position Attachment
		position_accessor_type m_aaPos;		///<Accessor
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
