/*
 * err_est_data_impl.h
 *
 *	Implementation of classes from err_est_data.h
 *
 *  Created on: 26.03.2014
 *     Authors: Dmitriy Logashenko, Markus Breit
 */

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
void SideFluxErrEstData<TDomain>::summarize_err_est_data ()
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

/// Constructor
template <typename TDomain>
SideAndElemErrEstData<TDomain>::SideAndElemErrEstData(std::size_t _sideOrder, std::size_t _elemOrder) :
	IErrEstData<TDomain>(),
	sideOrder(_sideOrder), elemOrder(_elemOrder),
	m_aSide(attachment_type("errEstSide")), m_aElem(attachment_type("errEstElem")),
	m_aaSide(MultiGrid::AttachmentAccessor<side_type, attachment_type >()),
	m_aaElem(MultiGrid::AttachmentAccessor<elem_type, attachment_type >()),
	m_spSV(SPNULL), m_errEstGL(GridLevel())
{
	// get quadrature rules for sides and elems
	boost::mpl::for_each<typename domain_traits<dim>::ManifoldElemList>(GetQuadRules<dim-1>(&quadRuleSide[0], sideOrder));
	boost::mpl::for_each<typename domain_traits<dim>::DimElemList>(GetQuadRules<dim>(&quadRuleElem[0], elemOrder));

	// fill in values for local side IPs (must be transformed from side local coords to elem local coords)
	// and fill IP indexing structure along the way
	for (ReferenceObjectID roid = ROID_VERTEX; roid != NUM_REFERENCE_OBJECTS; roid++)
	{
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

//	copy the parameters to the object
	m_errEstGL = gl;
	m_spSV = spSV;

//	prepare the attachments and their accessors
	MultiGrid* pMG = (MultiGrid*) (spSV->subset_handler()->multi_grid());

//	sides
	pMG->template attach_to_dv<side_type, attachment_type >(m_aSide, std::vector<number>(0));
	m_aaSide.access(*pMG, m_aSide);

//	elems
	pMG->template attach_to_dv<elem_type, attachment_type >(m_aElem, std::vector<number>(0));
	m_aaElem.access(*pMG, m_aElem);

//	now iterate over the grid and get the number of IPs from stored quadrature rules
//	sides first
	typedef typename SurfaceView::traits<side_type>::const_iterator side_iterator_type;
	side_iterator_type side_iter_end = m_spSV->template end<side_type> (m_errEstGL, SurfaceView::ALL);
	for (side_iterator_type side_iter = m_spSV->template begin<side_type> (m_errEstGL, SurfaceView::ALL);
		 side_iter != side_iter_end; ++side_iter)
	{
		side_type* side = *side_iter;

	//	get roid of side
		ReferenceObjectID roid = side->reference_object_id();

	//	get number of IPs from quadrature rule for the roid and specified side quadrature order
		std::size_t size = quadRuleSide[roid]->size();

	//	resize attachment accordingly
		m_aaSide[side].resize(size, 0.0);
	}
//	then elems
	typedef typename SurfaceView::traits<elem_type>::const_iterator elem_iterator_type;
	elem_iterator_type elem_iter_end = m_spSV->template end<elem_type> (m_errEstGL, SurfaceView::ALL);
	for (elem_iterator_type elem_iter = m_spSV->template begin<elem_type> (m_errEstGL, SurfaceView::ALL);
		 elem_iter != elem_iter_end; ++elem_iter)
	{
		elem_type* elem = *elem_iter;

	//	get roid of elem
		ReferenceObjectID roid = elem->reference_object_id();

	//	get number of IPs from quadrature rule for the roid and specified side quadrature order
		std::size_t size = quadRuleElem[roid]->size();

	//	resize attachment accordingly
		m_aaElem[elem].resize(size, 0.0);
	}
};

/// Called after the computation of the error estimator data in all the elements
/**	Because of the multigrid hierarchy, the sides at the rim of a multigrid level
 * 	only have one of the two flux terms (cis or trans). In order to calculate a
 * 	jump, we need to add up the resp. attachments on parent and child for that side.
 */
template <typename TDomain>
void SideAndElemErrEstData<TDomain>::summarize_err_est_data()
{
	const MultiGrid* pMG = m_spSV->subset_handler()->multi_grid();

//	loop the rim sides and add the jumps
	typedef typename SurfaceView::traits<side_type>::const_iterator t_iterator;
	t_iterator end_rim_side_iter = m_spSV->template end<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
	for (t_iterator rim_side_iter = m_spSV->template begin<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
		rim_side_iter != end_rim_side_iter; ++rim_side_iter)
	{
	//	get the sides on both the levels (coarse and fine)
		side_type* c_rim_side = *rim_side_iter;
		if (pMG->template num_children<side_type>(c_rim_side) != 1) // we consider _NONCOPY only, no hanging nodes // isn't that COPY?
			UG_THROW ("SideFluxErrEstData::summarize_error_estimator:"
					" The error estimator does not accept hanging nodes");
		side_type* f_rim_side = pMG->template get_child<side_type>(c_rim_side, 0);

	//	compute the total jump and save it for both the sides
		for (std::size_t i = 0; i < m_aaSide[c_rim_side].size(); i++)
		{
			number& c_rim_flux = m_aaSide[c_rim_side][i];
			number& f_rim_flux = m_aaSide[f_rim_side][i];
			number flux_jump = f_rim_flux + c_rim_flux;
			c_rim_flux = f_rim_flux = flux_jump;
		}

	// TODO: something similar for the parallel case!?
	}
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

	// check number of integration points
	std::size_t nIPs = quadRuleElem[roid]->size();
	std::vector<number>& integrand = m_aaElem[dynamic_cast<elem_type*>(pElem)];
	if (nIPs != integrand.size())
		UG_THROW("Element attachment vector does not have the required size for integration!");

	// get reference element mapping
	DimReferenceMapping<dim,dim>& mapping = ReferenceMappingProvider::get<dim,dim>(roid);
	mapping.update(&vCornerCoords[0]);

	//	compute det of jacobian at each IP
	std::vector<number> det = std::vector<number>(nIPs);
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
	//	get the sides of the element
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
summarize_err_est_data()
{
	// only called if consider_me()
	for (std::size_t eed = 0; eed < num(); eed++)
		m_vEed[eed]->summarize_err_est_data();
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
