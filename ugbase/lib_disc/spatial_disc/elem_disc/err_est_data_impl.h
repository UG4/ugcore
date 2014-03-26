/*
 * err_est_data_impl.h
 *
 *	Implementation of classes from err_est_data.h
 *
 *  Created on: 26.03.2014
 *      Author: Dmitriy Logashenko
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
	pMG->template attach_to_dv<side_type>(m_aFluxJumpOverSide, 0);
};

/// Called after the computation of the error estimator data in all the elements
template <typename TDomain>
void SideFluxErrEstData<TDomain>::summarize_err_est_data ()
{
	typedef typename domain_traits<dim>::side_type side_type;
	
	const MultiGrid * pMG = m_spSV->subset_handler()->multi_grid();
	
//	Get the access to the flux jumps
	MultiGrid::AttachmentAccessor<side_type, Attachment<number> > aaFluxJump
									(* (MultiGrid *) pMG, m_aFluxJumpOverSide);
	
//	Loop the rim sides and add the jumps
	typedef typename SurfaceView::traits<side_type>::const_iterator t_iterator;
	t_iterator end_rim_side_iter = m_spSV->template end<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
	for (t_iterator rim_side_iter = m_spSV->template begin<side_type> (m_errEstGL, SurfaceView::SHADOW_RIM);
		rim_side_iter != end_rim_side_iter; ++rim_side_iter)
	{
	//	Get the sides on both the levels
		side_type * c_rim_side = *rim_side_iter;
		if (pMG->template num_children<side_type>(c_rim_side) != 1) // we consider _NONCOPY only, no hanging nodes
			UG_THROW ("ConvectionDiffusionFE::summarize_error_estimator:"
					" The error estimator does not accept hanging nodes");
		side_type * f_rim_side = pMG->template get_child<side_type> (c_rim_side, 0);
		
	//	Compute the total jump and save it for both the sides:
		number & c_rim_flux = aaFluxJump[c_rim_side], & f_rim_flux = aaFluxJump[f_rim_side];
		number flux_jump = f_rim_flux - c_rim_flux;
		c_rim_flux = - (f_rim_flux = flux_jump);
	}
};

/// Releases data structures of the error estimator
template <typename TDomain>
void SideFluxErrEstData<TDomain>::release_err_est_data ()
{
//	Release the attachment
	typedef typename domain_traits<dim>::side_type side_type;
	MultiGrid * pMG = (MultiGrid *) (m_spSV->subset_handler()->multi_grid());
	pMG->template detach_from<side_type>(m_aFluxJumpOverSide);
	m_spSV = ConstSmartPtr<SurfaceView> (NULL);
};

} // end of namespace ug

/* End of File */
