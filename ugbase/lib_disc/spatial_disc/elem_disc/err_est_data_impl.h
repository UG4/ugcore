/*
 * err_est_data_impl.h
 *
 *	Implementation of classes from err_est_data.h
 *
 *  Created on: 26.03.2014
 *     Authors: Dmitriy Logashenko, Markus Breit
 */

#include "err_est_data.h"

namespace ug{

// ******** class SideFluxErrEstData ********

/// Allocates data structures for the error estimator
template <typename TDomain>
void SideFluxErrEstData<TDomain>::alloc_err_est_data
(
	ConstSmartPtr<SurfaceView> spSV,
	const GridLevel& gl
)
{
//	Get and check the grid level:
	if (gl.type () != GridLevel::SURFACE)
		UG_THROW("SideFluxErrEstData::alloc_err_est_data:"
			" The error estimator can work only with grid functions of the SURFACE type.");
	
//	Copy the parameters to the object:
	m_errEstGL = gl;
	m_spSV = spSV;
	
//	Prepare the attachment for the jumps of the fluxes over the sides:
	typedef typename domain_traits<dim>::side_type side_type;
	MultiGrid * pMG = (MultiGrid *) (spSV->subset_handler()->multi_grid());
	pMG->template attach_to_dv<side_type,ANumber>(m_aFluxJumpOverSide, 0);
	m_aaFluxJump.access(*pMG, m_aFluxJumpOverSide);
};

/// Called after the computation of the error estimator data in all the elements
template <typename TDomain>
void SideFluxErrEstData<TDomain>::summarize_err_est_data (SmartPtr<TDomain> spDomain)
{
	typedef typename domain_traits<dim>::side_type side_type;
	
	const MultiGrid * pMG = m_spSV->subset_handler()->multi_grid();
	
//	Loop the rim sides and add the jumps
	typedef typename SurfaceView::traits<side_type>::const_iterator t_iterator;
	t_iterator end_rim_side_iter = m_spSV->template end<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
	for (t_iterator rim_side_iter = m_spSV->template begin<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
		rim_side_iter != end_rim_side_iter; ++rim_side_iter)
	{
	//	Get the sides on both the levels
		side_type * c_rim_side = *rim_side_iter;
		if (pMG->template num_children<side_type>(c_rim_side) != 1) // we consider _NONCOPY only, no hanging nodes
			UG_THROW ("SideFluxErrEstData::summarize_error_estimator:"
					" The error estimator does not accept hanging nodes");
		side_type * f_rim_side = pMG->template get_child<side_type> (c_rim_side, 0);
		
	//	Compute the total jump and save it for both the sides:
		number & c_rim_flux = m_aaFluxJump[c_rim_side];
		number & f_rim_flux = m_aaFluxJump[f_rim_side];
		number flux_jump = f_rim_flux + c_rim_flux;
		c_rim_flux = f_rim_flux = flux_jump;
	}
};

/// Releases data structures of the error estimator
template <typename TDomain>
void SideFluxErrEstData<TDomain>::release_err_est_data ()
{
//	Release the attachment
	typedef typename domain_traits<dim>::side_type side_type;
	MultiGrid * pMG = (MultiGrid *) (m_spSV->subset_handler()->multi_grid());
	m_aaFluxJump.invalidate ();
	pMG->template detach_from<side_type> (m_aFluxJumpOverSide);
	//m_spSV = ConstSmartPtr<SurfaceView> (NULL);	// this raises a rte
};




// ******** class SideAndElemErrEstData ********

inline void check_subset_strings(std::vector<std::string> s)
{
	//	remove white space
	for (size_t i = 0; i < s.size(); i++)
		RemoveWhitespaceFromString(s[i]);

	//	if no subset passed, clear subsets
	if (s.size() == 1 && s[0].empty()) s.clear();

	//	if subsets passed with separator, but not all tokens filled, throw error
	for (size_t i = 0; i < s.size(); i++)
	{
		if (s.empty())
		{
			UG_THROW("Error while setting subsets in SideAndElemErrEstData: Passed subset string lacks a "
					 "subset specification at position " << i << "(of " << s.size()-1 << ")");
		}
	}
}


template <typename TDomain>
SideAndElemErrEstData<TDomain>::SideAndElemErrEstData
(
	std::size_t _sideOrder,
	std::size_t _elemOrder,
	const char* subsets
) :
	IErrEstData<TDomain>(),
	sideOrder(_sideOrder), elemOrder(_elemOrder),
	m_aSide(attachment_type("errEstSide")), m_aElem(attachment_type("errEstElem")),
	m_aaSide(MultiGrid::AttachmentAccessor<side_type, attachment_type >()),
	m_aaElem(MultiGrid::AttachmentAccessor<elem_type, attachment_type >()),
	m_spSV(SPNULL), m_errEstGL(GridLevel())
{
	m_vSs = TokenizeString(subsets);
	check_subset_strings(m_vSs);

	if (m_vSs.size() == 0)
	{
		UG_LOG("Warning: SideAndElemErrEstData ist constructed without definition of subsets. This is likely not to work.\n"
			   "Please specify a subset of the same dimension as your domain that the error estimator is supposed to work on.\n");
	}

	init_quadrature();
}


template <typename TDomain>
SideAndElemErrEstData<TDomain>::SideAndElemErrEstData
(
	std::size_t _sideOrder,
	std::size_t _elemOrder,
	std::vector<std::string> subsets
) :
	IErrEstData<TDomain>(),
	sideOrder(_sideOrder), elemOrder(_elemOrder),
	m_aSide(attachment_type("errEstSide")), m_aElem(attachment_type("errEstElem")),
	m_aaSide(MultiGrid::AttachmentAccessor<side_type, attachment_type >()),
	m_aaElem(MultiGrid::AttachmentAccessor<elem_type, attachment_type >()),
	m_spSV(SPNULL), m_errEstGL(GridLevel())
{
	m_vSs = subsets;
	check_subset_strings(m_vSs);

	if (m_vSs.size() == 0)
	{
		UG_LOG("Warning: SideAndElemErrEstData ist constructed without definition of subsets. This is likely not to work.\n"
					   "Please specify a subset of the same dimension as your domain that the error estimator is supposed to work on.\n");
	}

	init_quadrature();
}

template <typename TDomain>
void SideAndElemErrEstData<TDomain>::init_quadrature()
{
	// get quadrature rules for sides and elems
	boost::mpl::for_each<typename domain_traits<dim>::ManifoldElemList>(GetQuadRules<dim-1>(&quadRuleSide[0], sideOrder));
	boost::mpl::for_each<typename domain_traits<dim>::DimElemList>(GetQuadRules<dim>(&quadRuleElem[0], elemOrder));

	// fill in values for local side IPs (must be transformed from side local coords to elem local coords)
	// and fill IP indexing structure along the way
	for (ReferenceObjectID roid = ROID_VERTEX; roid != NUM_REFERENCE_OBJECTS; roid++)
	{
		// no reference elements for octahedra so far
		if (roid == ROID_OCTAHEDRON) continue;

		// get reference element for roid
		const ReferenceElement& re = ReferenceElementProvider::get(roid);
		int ref_dim = re.dimension();
		if (ref_dim != dim) continue;

		// make a DimReferenceElement from roid (we need access to corner coords)
		const DimReferenceElement<dim>& ref_elem = ReferenceElementProvider::get<dim>(roid);

		// clear IP coords
		m_SideIPcoords[roid].clear();

		// loop sides of ref elem
		for (std::size_t side = 0; side < ref_elem.num(ref_dim-1); side++)
		{
			// get side roid
			ReferenceObjectID side_roid = ref_elem.roid(ref_dim-1, side);

			// get number of IPs for this roid from quad rules
			if (!quadRuleSide[side_roid])
				UG_THROW("Requesting side IPs for roid " << roid << ", but no quadrature rule has been created for it.");
			std::size_t nIPs = quadRuleSide[side_roid]->size();

			// save start index for this side's IPs
			if (side == 0) m_sideIPsStartIndex[roid][side] = 0;
			else m_sideIPsStartIndex[roid][side] = m_sideIPsStartIndex[roid][side-1] + nIPs;

			// get the side-local IPs
			const MathVector<dim-1>* sideloc_IPs = quadRuleSide[side_roid]->points();

			// fill vector of side corners (in element-dimensional coords)
			std::size_t nCo = ref_elem.num(dim-1, side, 0);
			std::vector<MathVector<dim> > side_corners(nCo);
			for (std::size_t co = 0; co < ref_elem.num(dim-1, side, 0); co++)
			{
				std::size_t co_id = ref_elem.id(dim-1, side, 0, co);
				side_corners[co] = ref_elem.corner(co_id);
			}

			// get reference mapping
			DimReferenceMapping<dim-1,dim>& ref_map = ReferenceMappingProvider::get<dim-1,dim>(side_roid, &side_corners[0]);

			// map IPs
			for (std::size_t ip = 0; ip < nIPs; ip++)
			{
				m_SideIPcoords[roid].push_back(MathVector<TDomain::dim>());
				ref_map.local_to_global(m_SideIPcoords[roid].back(), sideloc_IPs[ip]);
			}
		}
	}
}


///	get the data reference for a given side and ip
template <typename TDomain>
number& SideAndElemErrEstData<TDomain>::operator() (side_type* pSide, std::size_t ip)
{
	try
	{
		return m_aaSide[pSide].at(ip);
	}
	catch (const std::out_of_range& oor)
	{
		UG_THROW("Requested attachment for side integration point " << ip <<
				 ", which does not appear to be allocated.");

	}
	UG_CATCH_THROW("Could not access error estimator side attachment for IP " << ip << ".");

	// silence no return warning; this code is never reached
	UG_THROW("Reached unreachable! This should be impossible. Check out why it is not!");
	return m_aaSide[pSide].at(ip);
}

template <typename TDomain>
number& SideAndElemErrEstData<TDomain>::operator() (elem_type* pElem, std::size_t ip)
{
	try
	{
		return m_aaElem[pElem].at(ip);
	}
	catch (const std::out_of_range& oor)
	{
		UG_THROW("Requested attachment for elem integration point " << ip <<
				 ", which does not appear to be allocated.");

	}
	UG_CATCH_THROW("Could not access error estimator elem attachment for IP " << ip << ".");

	// silence no return warning; this code is never reached
	UG_THROW("Reached unreachable! This should be impossible. Check out why it is not!");
	return m_aaElem[pElem].at(ip);
}


template <typename TDomain>
template <int refDim>
const MathVector<refDim>* SideAndElemErrEstData<TDomain>::side_local_ips(const ReferenceObjectID roid)
{
	// the usual case: return all side IPs of an element (belonging to all sides)
	if (TDomain::dim == refDim)
	{
	// check that IP series exists
		if (m_SideIPcoords[roid].size() == 0)
			UG_THROW("No side IP series available for roid " << roid << ".");

		// cast is necessary, since TDomain::dim might be != refDim,
		// but in that case, this return is not reached
		return reinterpret_cast<const MathVector<refDim>*>(&m_SideIPcoords[roid][0]);
	}
	// special case (assembling over manifold, needed for inner_boundary): only IPs of one side
	else if (TDomain::dim == refDim+1)
	{
		// check that quad rule exists
		if (!quadRuleSide[roid])
			UG_THROW("Requesting side IPs for roid " << roid << ", but no quadrature rule has been created for it.");

		return reinterpret_cast<const MathVector<refDim>*>(quadRuleSide[roid]->points());
	}

	UG_THROW("Local IPs requested with the wrong dimension refDim." << std::endl
			 << "Either call with refDim == TDomain::dim and a TDomain::dim-dimensional roid "
			 "for the local side IPs of all of its sides" << std::endl
			 << "or with refDim == TDomain::dim-1 and a (TDomain::dim-1)-dimensional roid"
			 "for the local side IPs for a side of this roid.");

	return NULL;
}

template <typename TDomain>
template <int refDim>
const MathVector<refDim>* SideAndElemErrEstData<TDomain>::elem_local_ips(const ReferenceObjectID roid)
{
//	return NULL if dim is not fitting (not meaningful for the purpose of error estimation)
	if (TDomain::dim != refDim)
		return NULL;

//	check that quad rule exists
	if (!quadRuleElem[roid])
		UG_THROW("Requesting side IPs for roid " << roid << ", but no quadrature rule has been created for it.");
//	check that IP series exists
	if (quadRuleElem[roid]->size() == 0)
		UG_THROW("No elem IP series available for roid " << roid << ".");

	// cast is necessary, since TDomain::dim might be != refDim,
	// but in that case, this return is not reached
	return reinterpret_cast<const MathVector<refDim>*>(quadRuleElem[roid]->points());
}

template <typename TDomain>
void SideAndElemErrEstData<TDomain>::all_side_global_ips
(
	MathVector<dim>* globIPs,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get reference object ID
	ReferenceObjectID roid = elem->reference_object_id();

	// get reference mapping
	DimReferenceMapping<dim,dim>& ref_map =
		ReferenceMappingProvider::get<dim,dim>(roid, &vCornerCoords[0]);

	// map IPs
	try
	{
		ref_map.local_to_global(globIPs, &m_SideIPcoords[roid][0], m_SideIPcoords[roid].size());
	}
	catch (std::exception& e)
	{
		UG_THROW("Encountered exception while trying to fill array of global IPs: "
				 << std::endl << "'" << e.what() << "'");
	}
}

template <typename TDomain>
void SideAndElemErrEstData<TDomain>::side_global_ips
(
	MathVector<dim>* globIPs,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get reference object ID
	ReferenceObjectID roid = elem->reference_object_id();

	// get reference mapping
	DimReferenceMapping<dim-1,dim>& ref_map =
		ReferenceMappingProvider::get<dim-1,dim>(roid, &vCornerCoords[0]);

	// map IPs
	try
	{
		ref_map.local_to_global(globIPs, quadRuleSide[roid]->points(), quadRuleSide[roid]->size());
	}
	catch (std::exception& e)
	{
		UG_THROW("Encountered exception while trying to fill array of global IPs: "
				 << std::endl << "'" << e.what() << "'");
	}
}

template <typename TDomain>
void SideAndElemErrEstData<TDomain>::elem_global_ips
(
	MathVector<dim>* globIPs,
	GridObject* elem,
	const MathVector<dim> vCornerCoords[]
)
{
	// get reference object ID
	ReferenceObjectID roid = elem->reference_object_id();

	// get reference mapping
	DimReferenceMapping<dim,dim>& ref_map =
		ReferenceMappingProvider::get<dim,dim>(roid, &vCornerCoords[0]);

	// map IPs
	try
	{
		for (std::size_t ip = 0; ip < quadRuleElem[roid]->size(); ip++)
			ref_map.local_to_global(globIPs[ip], quadRuleElem[roid]->point(ip));
	}
	catch (std::exception e)
	{
		UG_THROW("Encountered exception while trying to fill array of global IPs: "
				 << std::endl << "'" << e.what() << "'");
	}
}

template <typename TDomain>
std::size_t SideAndElemErrEstData<TDomain>::num_side_ips(const side_type* pSide)
{
	return m_aaSide[pSide].size();
}

template <typename TDomain>
std::size_t SideAndElemErrEstData<TDomain>::num_side_ips(const ReferenceObjectID roid)
{
	// check that quad rule exists
	if (!quadRuleSide[roid])
		UG_THROW("Requesting number of side IPs for roid " << roid << ", but no quadrature rule has been created for it.");

	return quadRuleSide[roid]->size();
}

template <typename TDomain>
std::size_t SideAndElemErrEstData<TDomain>::first_side_ips(const ReferenceObjectID roid, const std::size_t side)
{
	return m_sideIPsStartIndex[roid][side];
}


template <typename TDomain>
std::size_t SideAndElemErrEstData<TDomain>::num_all_side_ips(const ReferenceObjectID roid)
{
	return m_SideIPcoords[roid].size();
}

template <typename TDomain>
std::size_t SideAndElemErrEstData<TDomain>::num_elem_ips(const ReferenceObjectID roid)
{
	// check that quad rule exists
	if (!quadRuleElem[roid])
		UG_THROW("Requesting elem IPs for roid " << roid << ", but no quadrature rule has been created for it.");

	return quadRuleElem[roid]->size();
}

template <typename TDomain>
std::size_t SideAndElemErrEstData<TDomain>::side_ip_index
(	const ReferenceObjectID roid,
	const std::size_t side,
	const std::size_t ip
)
{
	// TODO: check validity of side index
	return m_sideIPsStartIndex[roid][side];
}


/// Allocates data structures for the error estimator
template <typename TDomain>
void SideAndElemErrEstData<TDomain>::alloc_err_est_data
(
	ConstSmartPtr<SurfaceView> spSV,
	const GridLevel& gl
)
{
//	get and check the grid level
	if (gl.type () != GridLevel::SURFACE)
		UG_THROW("SideFluxErrEstData::alloc_err_est_data:"
			" The error estimator can work only with grid functions of the SURFACE type.");

//	get the subset handler
	ConstSmartPtr<MGSubsetHandler> ssh = spSV->subset_handler();

//	copy the parameters to the object
	m_errEstGL = gl;
	m_spSV = spSV;

//	prepare the attachments and their accessors
	MultiGrid* pMG = (MultiGrid*) (ssh->multi_grid());

//	sides
	pMG->template attach_to_dv<side_type, attachment_type >(m_aSide, std::vector<number>(0));
	m_aaSide.access(*pMG, m_aSide);

//	elems
	pMG->template attach_to_dv<elem_type, attachment_type >(m_aElem, std::vector<number>(0));
	m_aaElem.access(*pMG, m_aElem);

//	construct subset group from subset info
	m_ssg.set_subset_handler(ssh);
	if (m_vSs.size() == 0) m_ssg.add_all();
	else m_ssg.add(m_vSs);

//	find out whether we work on an elem subset or a side subset
	int ssDimMax = 0;
	int ssDimMin = 4;
	for (size_t si = 0; si < m_ssg.size(); si++)
	{
		const int dim_si = DimensionOfSubset(*ssh, m_ssg[si]);
		if (dim_si > ssDimMax) ssDimMax = dim_si;
		if (dim_si < ssDimMin) ssDimMin = dim_si;
	}

	if (ssDimMax != ssDimMin || ssDimMax != TDomain::dim)
	{
		UG_THROW("Subsets passed to an instance of SideAndElemErrEstData have varying or inadmissable dimension.\n"
				 "(NOTE: Only sets of subsets of the same dimension as the domain are allowed.)");
	}

	// iterate over elems and associated sides and resize the IP value vectors
	for (size_t si = 0; si < m_ssg.size(); si++)
	{
		typedef typename SurfaceView::traits<elem_type>::const_iterator elem_iterator_type;
		elem_iterator_type elem_iter_end = m_spSV->template end<elem_type>(m_ssg[si], m_errEstGL, SurfaceView::ALL);
		for (elem_iterator_type elem_iter = m_spSV->template begin<elem_type>(m_ssg[si], m_errEstGL, SurfaceView::ALL);
			 elem_iter != elem_iter_end; ++elem_iter)
		{
			elem_type* elem = *elem_iter;

		//	get roid of elem
			ReferenceObjectID roid = elem->reference_object_id();

		//	get number of IPs from quadrature rule for the roid and specified side quadrature order
			std::size_t size = quadRuleElem[roid]->size();

		//	resize attachment accordingly
			m_aaElem[elem].resize(size, 0.0);

		// loop associated sides
			typename MultiGrid::traits<side_type>::secure_container side_list;
			pMG->associated_elements(side_list, elem);

			for (std::size_t side = 0; side < side_list.size(); side++)
			{
				side_type* pSide = side_list[side];

			//	get roid of side
				ReferenceObjectID roid = pSide->reference_object_id();

			//	get number of IPs from quadrature rule for the roid and specified side quadrature order
				std::size_t size = quadRuleSide[roid]->size();

			//	resize attachment accordingly
				if (m_aaSide[pSide].size() != size)
					m_aaSide[pSide].resize(size, 0.0);
			}
		}
	}

	// equalize sizes at horizontal interface
	// done by summing up, which, in a first step, will resize the smaller of the two vectors
	// to the same size as the bigger one (cf. std_number_vector_attachment_reduce_traits)
	// the actual summing does not hurt, since it performs 0 = 0 + 0;
#ifdef UG_PARALLEL
	typedef typename GridLayoutMap::Types<side_type>::Layout layout_type;
	DistributedGridManager& dgm = *pMG->distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();
	pcl::InterfaceCommunicator<layout_type> icom;

	// sum all copies at the h-master attachment
	ComPol_AttachmentReduce<layout_type, attachment_type> compolSumSideErr(*pMG, m_aSide, PCL_RO_SUM);
	icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSumSideErr);
	icom.communicate();
	icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolSumSideErr);
	icom.communicate();

	icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolSumSideErr);
	icom.communicate();
#endif

	// for the case where subset border goes through SHADOW_RIM:
	// resize SHADOW_RIM children if parent is resized and vice-versa
	typedef typename Grid::traits<side_type>::const_iterator side_iter_type;
	side_iter_type end_side = pMG->template end<side_type>();
	for (side_iter_type side_iter = pMG->template begin<side_type>(); side_iter != end_side; ++side_iter)
	{
		// get the sides on both the levels (coarse and fine)
		side_type* c_rim_side = *side_iter;

		// if side is not SHADOW_RIM: do nothing
		if (!m_spSV->surface_state(c_rim_side).partially_contains(SurfaceView::SHADOW_RIM)) continue;

		// loop rim side children
		bool resize_parent = false;
		const size_t num_children = pMG->template num_children<side_type>(c_rim_side);
		for (std::size_t ch = 0; ch < num_children; ch++)
		{
			side_type* child_side = pMG->template get_child<side_type>(c_rim_side, ch);

			// get roid of side
			ReferenceObjectID child_roid = child_side->reference_object_id();

			// get number of IPs from quadrature rule for the roid and specified side quadrature order
			std::size_t child_size = quadRuleSide[child_roid]->size();

			// resize attachment accordingly
			if (m_aaSide[child_side].size() != child_size && num_side_ips(c_rim_side) > 0)
				m_aaSide[child_side].resize(child_size, 0.0);
			else if (m_aaSide[child_side].size() == child_size && num_side_ips(c_rim_side) == 0)
				resize_parent = true;
		}

		if (resize_parent)
		{
			// get roid of side
			ReferenceObjectID roid = c_rim_side->reference_object_id();

			// get number of IPs from quadrature rule for the roid and specified side quadrature order
			std::size_t size = quadRuleSide[roid]->size();

			// resize attachment accordingly
			m_aaSide[c_rim_side].resize(size, 0.0);
		}
	}

#ifdef UG_PARALLEL
	icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolSumSideErr);
	icom.communicate();
#endif
};

/// Called after the computation of the error estimator data in all the elements
/**	Because of the multigrid hierarchy, the sides at the rim of a multigrid level
 * 	only have one of the two flux terms (cis or trans). In order to calculate a
 * 	jump, we need to add up the resp. attachments on parent and child for that side.
 */
template <typename TDomain>
void SideAndElemErrEstData<TDomain>::summarize_err_est_data(SmartPtr<TDomain> spDomain)
{
	MultiGrid* pMG = spDomain->subset_handler()->multi_grid();

	// STEP 1: Ensure consistency over horizontal interfaces.
	// Add the IP values of horizontal interfaces. Care has to be taken in the case where a side
	// is SHADOW_RIM, then it only has values from one side of the interface. This will be checked
	// by calling size() of the IP value vector.
#ifdef UG_PARALLEL
	typedef typename GridLayoutMap::Types<side_type>::Layout layout_type;
	DistributedGridManager& dgm = *pMG->distributed_grid_manager();
	GridLayoutMap& glm = dgm.grid_layout_map();
	pcl::InterfaceCommunicator<layout_type> icom;

	// sum all copies at the h-master attachment
	ComPol_AttachmentReduce<layout_type, attachment_type> compolSumSideErr(*pMG, m_aSide, PCL_RO_SUM);
	icom.exchange_data(glm, INT_H_SLAVE, INT_H_MASTER, compolSumSideErr);
	icom.communicate();

	// copy the sum from the master to all of its slave-copies
	ComPol_CopyAttachment<layout_type, attachment_type> compolCopySideErr(*pMG, m_aSide);
	icom.exchange_data(glm, INT_H_MASTER, INT_H_SLAVE, compolCopySideErr);
	icom.communicate();

	// STEP 2: Ensure correct values in vertical master rim side children.
	// since we're copying from vmasters to vslaves later on, we have to make
	// sure, that also all v-masters contain the correct values.
	// todo: communication to vmasters may not be necessary here...
	//       it is performed to make sure that all surface-rim-children
	//       contain their true value.
	icom.exchange_data(glm, INT_V_SLAVE, INT_V_MASTER, compolCopySideErr);
	icom.communicate();
#endif

	// STEP 3: Calculate values on rim sides.
	// loop the rim sides and add the jumps
	// TODO: not sure whether we cannot simply loop with
	//       m_spSV->template end<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM)

	//typedef typename SurfaceView::traits<side_type>::const_iterator t_iterator;
	//t_iterator end_rim_side_iter = m_spSV->template end<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
	//for (t_iterator rim_side_iter = m_spSV->template begin<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
	//	rim_side_iter != end_rim_side_iter; ++rim_side_iter)
	typedef typename Grid::traits<side_type>::const_iterator side_iter_type;
	side_iter_type end_side = pMG->template end<side_type>();
	for (side_iter_type side_iter = pMG->template begin<side_type>(); side_iter != end_side; ++side_iter)
	{
		// get the sides on both the levels (coarse and fine)
		side_type* c_rim_side = *side_iter;

		// if side is not SHADOW_RIM: do nothing
		if (!m_spSV->surface_state(c_rim_side).partially_contains(SurfaceView::SHADOW_RIM)) continue;

		// if no ips registered on this side: do nothing
		if (num_side_ips(c_rim_side) == 0) continue;

		// distinguish hanging nodes and clozure elements by number of children
		const size_t num_children = pMG->template num_children<side_type>(c_rim_side);

		// closure elements
		if (num_children == 1)
		{
			side_type* f_rim_side = pMG->template get_child<side_type>(c_rim_side, 0);

			// add up total jump and save it for both sides
			for (std::size_t i = 0; i < m_aaSide[c_rim_side].size(); i++)
			{
				number& c_rim_flux = m_aaSide[c_rim_side][i];
				number& f_rim_flux = m_aaSide[f_rim_side][i];
				number flux_jump = f_rim_flux + c_rim_flux;
				c_rim_flux = f_rim_flux = flux_jump;
			}
		}

		// hanging node
		else
		{
			// TODO: This is a somehow dirty hack (but will converge with h->0),
			// a correct interpolation with all given points would be preferable.
			// And even if it is not: Check whether it could be beneficial to use structures like k-d tree,
			// since this is the most naive implementation possible. I guess, the number of IPs
			// is likely to be very low in most of the applications.

			// The IPs are not in the same locations on both sides:
			// Interpolate values for both sides, then add up and distribute.
			// Interpolation is piecewise constant (use value of nearest neighbour).

			// map coarse side local IPs to global
			std::vector<MathVector<dim> > c_coCo;
			CollectCornerCoordinates(c_coCo, c_rim_side, spDomain->position_accessor(), false);
			std::vector<MathVector<dim> > c_gloIPs(num_side_ips(c_rim_side));

			side_global_ips(&c_gloIPs[0], c_rim_side, &c_coCo[0]);

			// interpolate fine IPs on coarse side
			for (std::size_t ch = 0; ch < num_children; ch++)
			{
				// get fine side
				side_type* f_rim_side = pMG->template get_child<side_type>(c_rim_side, ch);

				// map fine side local IPs to global
				std::vector<MathVector<dim> > f_coCo;
				CollectCornerCoordinates(f_coCo, f_rim_side, spDomain->position_accessor(), false);
				std::vector<MathVector<dim> > f_gloIPs(num_side_ips(f_rim_side));
				side_global_ips(&f_gloIPs[0], f_rim_side, &f_coCo[0]);

				// for each fine side IP, find nearest coarse side IP
				for (std::size_t fip = 0; fip < m_aaSide[f_rim_side].size(); fip++)
				{
					if (m_aaSide[c_rim_side].size() < 1)
						{UG_THROW("No IP defined for coarse side.");}

					std::size_t nearest = 0;
					MathVector<dim> diff;
					VecSubtract(diff, f_gloIPs[fip], c_gloIPs[0]);
					number dist = VecLengthSq(diff);

					// loop coarse IPs
					for (std::size_t cip = 1; cip < m_aaSide[c_rim_side].size(); cip++)
					{
						VecSubtract(diff, f_gloIPs[fip], c_gloIPs[cip]);
						if (VecLengthSq(diff) < dist)
						{
							dist = VecLengthSq(diff);
							nearest = cip;
						}
					}
					m_aaSide[f_rim_side][fip] += m_aaSide[c_rim_side][nearest];
				}
			}

			// interpolate coarse IPs on fine side:
			for (std::size_t cip = 0; cip < m_aaSide[c_rim_side].size(); cip++)
			{
				MathVector<dim> diff;
				number dist = std::numeric_limits<number>::infinity();
				number val = 0.0;

				// we have to loop all child sides
				for (std::size_t ch = 0; ch < pMG->template num_children<side_type>(c_rim_side); ch++)
				{
					// get fine side
					side_type* f_rim_side = pMG->template get_child<side_type>(c_rim_side, ch);

					// map fine side local IPs to global
					std::vector<MathVector<dim> > f_coCo;
					CollectCornerCoordinates(f_coCo, f_rim_side, spDomain->position_accessor(), false);
					std::vector<MathVector<dim> > f_gloIPs(num_side_ips(f_rim_side));
					side_global_ips(&f_gloIPs[0], f_rim_side, &f_coCo[0]);

					// loop coarse IPs
					for (std::size_t fip = 1; fip < m_aaSide[f_rim_side].size(); fip++)
					{
						VecSubtract(diff, f_gloIPs[fip], c_gloIPs[cip]);
						if (VecLengthSq(diff) < dist)
						{
							dist = VecLengthSq(diff);
							val = m_aaSide[f_rim_side][fip];
						}
					}
				}

				// the fine grid already contains the correct values
				m_aaSide[c_rim_side][cip] = val;
			}
		}
	}

	// STEP 4: Ensure consistency in slave rim sides.
	// Copy from v-masters to v-slaves, since there may be constrained sides which locally
	// do not have a constraining element. Note that constrained V-Masters always have a local
	// constraining element and thus contain the correct value.
#ifdef UG_PARALLEL
	icom.exchange_data(glm, INT_V_MASTER, INT_V_SLAVE, compolCopySideErr);
	icom.communicate();
#endif
};


template <typename TDomain>
number SideAndElemErrEstData<TDomain>::get_elem_error_indicator(GridObject* pElem, const MathVector<dim> vCornerCoords[])
{

// the indicator
	number etaSq = 0.0;

// elem terms
	// info about reference element type
	ReferenceObjectID roid = pElem->reference_object_id();
	const DimReferenceElement<dim>& refElem = ReferenceElementProvider::get<dim>(roid);

	// only take into account elem contributions of the subsets defined in the constructor
	int ssi = surface_view()->subset_handler()->get_subset_index(pElem);
	if (!m_ssg.contains(ssi)) return 0.0;

	// check number of integration points
	std::size_t nIPs = quadRuleElem[roid]->size();
	std::vector<number>& integrand = m_aaElem[dynamic_cast<elem_type*>(pElem)];
	if (nIPs != integrand.size())
		UG_THROW("Element attachment vector does not have the required size for integration!");

	// get reference element mapping
	DimReferenceMapping<dim,dim>& mapping = ReferenceMappingProvider::get<dim,dim>(roid);
	mapping.update(&vCornerCoords[0]);

	//	compute det of jacobian at each IP
	std::vector<number> det(nIPs);
	mapping.sqrt_gram_det(&det[0], quadRuleElem[roid]->points(), nIPs);

	// integrate
	number sum = 0.0;
	for (std::size_t ip = 0; ip < nIPs; ip++)
		sum += quadRuleElem[roid]->weight(ip) * std::pow(integrand[ip], 2.0) * det[ip];

	// scale by diam^2(elem)
	// c* vol(elem) >= diam^3(elem) >= vol(elem)
	// therefore, up to a constant, error estimator can calculate diam(elem) as (vol(elem))^(1/3)
	number diamSq = std::pow(ElementSize<dim>(roid, &vCornerCoords[0]), 2./dim);

	// add to error indicator
	etaSq += diamSq * sum;

// side terms
	// get the sides of the element
	MultiGrid* pErrEstGrid = (MultiGrid*) (surface_view()->subset_handler()->multi_grid());
	typename MultiGrid::traits<side_type>::secure_container side_list;
	pErrEstGrid->associated_elements(side_list, pElem);

	// loop sides
	number diam;
	for (std::size_t side = 0; side < side_list.size(); side++)
	{
		side_type* pSide = side_list[side];

		// info about reference side type
		ReferenceObjectID side_roid = pSide->reference_object_id();

		// check number of integration points
		std::size_t nsIPs = quadRuleSide[side_roid]->size();
		if (nsIPs != m_aaSide[pSide].size())
			UG_THROW("Side attachment vector does not have the required size for integration!");

		// get side corners
		std::vector<MathVector<dim> > vSideCornerCoords(0);
		std::size_t nsCo = refElem.num(dim-1, side, 0);
		for (std::size_t co = 0; co < nsCo; co++)
			vSideCornerCoords.push_back(vCornerCoords[refElem.id(dim-1, side, 0, co)]);

		// get reference element mapping
		DimReferenceMapping<dim-1,dim>& mapping = ReferenceMappingProvider::get<dim-1,dim>(side_roid);
		mapping.update(&vSideCornerCoords[0]);

		//	compute det of jacobian at each IP
		det.resize(nsIPs);
		mapping.sqrt_gram_det(&det[0], quadRuleSide[side_roid]->points(), nsIPs);

		// integrate
		number sum = 0.0;
		for (std::size_t ip = 0; ip < nsIPs; ip++)
			sum += quadRuleSide[side_roid]->weight(ip) * std::pow(m_aaSide[pSide][ip], 2.0) * det[ip];

		// scale by diam(side)
		// c* vol(side) >= diam^2(side) >= vol(side)
		// therefore, up to a constant, error estimator can calculate diam as sqrt(vol(side))
		if (dim == 1)
			diam = 1.0;
		else
		{
			number exponent = dim-1;
			diam = std::pow(ElementSize<dim>(side_roid, &vSideCornerCoords[0]), 1.0/exponent);
		}
		// add to error indicator
		etaSq += diam * sum;
	}

	return etaSq;
}


/// Releases data structures of the error estimator
template <typename TDomain>
void SideAndElemErrEstData<TDomain>::release_err_est_data ()
{
//	release the attachments
	MultiGrid * pMG = (MultiGrid *) (m_spSV->subset_handler()->multi_grid());

	m_aaSide.invalidate();
	pMG->template detach_from<side_type>(m_aSide);

	m_aaElem.invalidate();
	pMG->template detach_from<elem_type>(m_aElem);
	//m_spSV = ConstSmartPtr<SurfaceView> (NULL);	// this raises a rte
};



// ******** class MultipleErrEstData ********

template <typename TDomain, typename TErrEstData>
void MultipleErrEstData<TDomain,TErrEstData>::
alloc_err_est_data(ConstSmartPtr<SurfaceView> spSV, const GridLevel& gl)
{
	// only called if consider_me()
	for (std::size_t eed = 0; eed < num(); eed++)
		m_vEed[eed]->alloc_err_est_data(spSV, gl);
}

template <typename TDomain, typename TErrEstData>
void MultipleErrEstData<TDomain,TErrEstData>::
summarize_err_est_data(SmartPtr<TDomain> spDomain)
{
	// only called if consider_me()
	for (std::size_t eed = 0; eed < num(); eed++)
		m_vEed[eed]->summarize_err_est_data(spDomain);
}

template <typename TDomain, typename TErrEstData>
number MultipleErrEstData<TDomain,TErrEstData>::
get_elem_error_indicator(GridObject* elem, const MathVector<dim> vCornerCoords[])
{
	// only called if consider_me()
	number sum = 0.0;
	for (std::size_t eed = 0; eed < num(); eed++)
		sum += m_vEed[eed]->get_elem_error_indicator(elem, vCornerCoords);

	return sum;
}

template <typename TDomain, typename TErrEstData>
void MultipleErrEstData<TDomain,TErrEstData>::
release_err_est_data()
{
	// only called if consider_me()
	for (std::size_t eed = 0; eed < num(); eed++)
		m_vEed[eed]->release_err_est_data();
}



// ******** class MultipleSideAndElemErrEstData ********

template <typename TDomain>
void MultipleSideAndElemErrEstData<TDomain>::add(SmartPtr<SideAndElemErrEstData<TDomain> > spEed)
{
	this->m_vEed.push_back(spEed.get());
	check_equal_order();
}

template <typename TDomain>
void MultipleSideAndElemErrEstData<TDomain>::check_equal_order()
{
	check_equal_side_order();
	check_equal_elem_order();
}

template <typename TDomain>
void MultipleSideAndElemErrEstData<TDomain>::check_equal_side_order()
{
	m_bEqSideOrder = false;

	if (this->m_vEed.size() == 0)
	{
		m_bEqSideOrder = true;
		return;
	}

	std::size_t side_order = this->m_vEed[0]->side_order();

	for (std::size_t ee = 1; ee < this->m_vEed.size(); ee++)
		if (this->m_vEed[ee]->side_order() != side_order) return;

	m_bEqSideOrder = true;
}

template <typename TDomain>
void MultipleSideAndElemErrEstData<TDomain>::check_equal_elem_order()
{
	m_bEqElemOrder = false;

	if (this->m_vEed.size() == 0)
	{
		m_bEqElemOrder = true;
		return;
	}

	std::size_t elem_order = this->m_vEed[0]->elem_order();

	for (std::size_t ee = 1; ee < this->m_vEed.size(); ee++)
		if (this->m_vEed[ee]->elem_order() != elem_order) return;

	m_bEqElemOrder = true;
}


} // end of namespace ug

/* End of File */
