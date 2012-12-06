/*
 * domain.h
 *
 *  Created on: 17.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__DOMAIN__
#define __H__UG__LIB_DISC__DOMAIN__

#include "lib_grid/algorithms/subset_util.h"

#ifdef UG_PARALLEL
#include "lib_grid/parallelization/distributed_grid.h"
#include "lib_grid/parallelization/util/parallel_subset_util.h"
#endif

namespace ug{

/**
 * Domain
 *
 * \defgroup lib_disc_domain Domain
 * \ingroup lib_discretization
 */

/// \ingroup lib_disc_domain
/// @{

///	Describes the contents of a domain.
/**	In a parallel environment those are the contents of the global distributed domain.
*/
class DomainInfo
{
	public:
		inline GeometricBaseObject element_type()	const		{return m_elementType;}
		inline size_t num_levels() const						{return m_numElemsOnLvl.size();}
		inline int num_elements_on_level(size_t lvl) const	{return m_numElemsOnLvl[lvl];}

		inline int num_subsets() const
			{return (int)m_subsetDims.size();}
		inline int subset_dim(int si) const
			{return m_subsetDims[si];}

		inline void set_info(GeometricBaseObject elemType,
								const std::vector<int>& numElemsOnLvl,
								const std::vector<int>& subsetDims)
			{m_elementType = elemType; m_numElemsOnLvl = numElemsOnLvl; m_subsetDims = subsetDims;}

	private:
		GeometricBaseObject	m_elementType;
		std::vector<int>	m_numElemsOnLvl;
		std::vector<int>	m_subsetDims;
};


/// describes a physical domain
/**
 * A Domain collects and exports relevant informations about the
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

	///	updates the subsets dimension property "dim" (integer) locally.
	/** This method normally only has to be called after the subset-handler
	 * has been changed. In most cases a call after the domain has been loaded
	 * should be enough. This is done automatically by e.g. LoadDomain.
	 *
	 * \todo	calling this method after coarsening could be useful.
	 * 			Think about registering a callback at the message-hub.*/
		void update_local_subset_dim_property()	{UpdateMaxDimensionOfSubset(*this->m_spSH, "dim");}


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

	protected:
		SmartPtr<TGrid> m_spGrid;			///< Grid
		SmartPtr<TSubsetHandler> m_spSH;	///< Subset Handler
		MessageHub::SPCallbackId m_spGridAdaptionCallbackID; ///< SmartPointer to grid adaption callback id

		DomainInfo	m_domainInfo;

		bool	m_isAdaptive;
		bool	m_adaptionIsActive;

	/**	this callback is called by the message hub, when a grid change has been
	 * performed. It will call all necessary actions in order to keep the grid
	 * correct for computations. */
		inline void grid_changed_callback(int, const GridMessage_Adaption* msg);

#ifdef UG_PARALLEL
	public:
	///	updates the subsets dimension property "dim" (integer) globally.
	/** This method normally only has to be called after the subset-handler
	 * has been changed. In most cases a call after the domain has been loaded
	 * should be enough. This is done automatically by e.g. LoadDomain.
	 *
	 * You may optionally specify a process communicator, which will be used
	 * to perform the communication involved.
	 *
	 * \todo	calling this method after coarsening could be useful.
	 * 			Think about registering a callback at the message-hub.*/
		void update_global_subset_dim_property(pcl::ProcessCommunicator procCom =
													pcl::ProcessCommunicator())
		{UpdateGlobalMaxDimensionOfSubset(*this->m_spSH, "dim", procCom);}

	protected:
	///	updates the number of levels on each block to global maximum
	/**	This method must be invoked when the grid has been changed and there
	 * might be diffent numbers of levels on different processors. This method
	 * will compute the global maximum of levels (involving an allreduce) and
	 * adapt all local grids to the maximum number by inserting empty levels.
	 * After this method has been invoked, all local grids will return the same
	 * number of levels, when invoking num_levels().*/
		void update_local_multi_grid();
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

/// @}

// include implementation
#include "domain_impl.h"

#endif /* __H__UG__LIB_DISC__DOMAIN__ */
