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
 * Classes for user data that compute normals or project vectors to low-dimensional
 * subsets. Note that these classes are mainly implemented for visualization purposes
 * and cannot be used for coupling.
 */
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__ONC_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__ONC_USER_DATA__

#include <vector>

// ug4 headers
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_grid/grid_objects/grid_dim_traits.h"
#include "lib_disc/domain_traits.h"
#include "lib_disc/domain_util.h"
#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "lib_disc/spatial_disc/user_data/std_user_data.h"

namespace ug {

/// Projection to the outer normal for a given vector user data
/**
 * This class implements the UserData for evaluation of the normal component of vectors
 * computed by the other user data object. The normal component is evaluation only for
 * WDim-1-dimensional subsets - sides of WDim-dimensional elements.
 *
 * ToDo: Note that the LocalVector provided for the evaluate function is not properly
 * used. Thus, the passed UserData objects should not use it.
 *
 */
template <typename TDomain>
class OutNormCmp
	: public StdUserData<OutNormCmp<TDomain>, MathVector<TDomain::dim>, TDomain::dim, void, UserData<MathVector<TDomain::dim>, TDomain::dim, void> >
{
///	the world dimension
	static const int dim = TDomain::dim;
	
///	the domain type
	typedef TDomain domain_type;
	
///	the grid type
	typedef typename TDomain::grid_type grid_type;
	
///	the original data type
	typedef MathVector<dim> vec_type;
	
///	the original data
	SmartPtr<UserData<vec_type, dim, void> > m_spData;
	
///	the domain
	SmartPtr<domain_type> m_spDomain;
	
///	the subset group for the inner part
	SubsetGroup m_ssGrp;
	
public:

///	constructor
	OutNormCmp
	(
		SmartPtr<domain_type> spDomain, ///< the domain
		SmartPtr<UserData<vec_type, dim, void> > spData, ///< the original data
		const char * ss_names ///< the subsets
	)
	:	m_spData (spData), m_spDomain (spDomain), m_ssGrp (spDomain->subset_handler ())
	{
	// Parse the subset names:
		std::vector<std::string> vssNames;
		try
		{
			TokenizeString (ss_names, vssNames);
			for (size_t k = 0; k < vssNames.size (); k++)
				RemoveWhitespaceFromString (vssNames [k]);
			m_ssGrp.clear ();
			m_ssGrp.add (vssNames);
		} UG_CATCH_THROW ("SubsetIndicatorUserData: Failed to parse subset names.");
	}
	
///	constructor
	OutNormCmp
	(
		SmartPtr<domain_type> spDomain, ///< the domain
		SmartPtr<UserData<vec_type, dim, void> > spData ///< the original data
	)
	:	m_spData (spData), m_spDomain (spDomain), m_ssGrp (spDomain->subset_handler ())
	{}
	
///	Indicator functions are discontinuous
	virtual bool continuous () const {return false;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return m_spData->requires_grid_fct ();}
	// Remark: Note that actually we have not enought data for that. This dependence is neglected.

///	Evaluator
	template <int refDim>
	inline void evaluate
	(
		vec_type vValue[],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		GridObject * elem,
		const MathVector<dim> vCornerCoords [],
		const MathVector<refDim> vLocIP [],
		const size_t nip,
		LocalVector * u,
		const MathMatrix<refDim, dim> * vJT = NULL
	) const
	{
		if (refDim == dim)
		{
			(* m_spData) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			return;
		}
		if (refDim != dim - 1)
		{
			UG_THROW ("OutNormCmp: Wrong dimensionality of the subset.");
		}
		
	//	Get the full-dim. element associated with the given side
		typedef typename grid_dim_traits<dim>::grid_base_object fd_elem_type;
		fd_elem_type * fd_elem = NULL;
		int fd_si = -1;
		
		typedef typename Grid::traits<fd_elem_type>::secure_container fd_elem_secure_container_type;
		fd_elem_secure_container_type fd_elem_list;
		grid_type& grid = * (grid_type*) (m_spDomain->grid().get ());
		grid.associated_elements (fd_elem_list, elem);
		
		if (m_ssGrp.empty ()) // no subsets specified; we assume, elem should be the bnd of the whole domain
		{
			if (fd_elem_list.size () == 1)
			{
				fd_elem = fd_elem_list[0];
				fd_si = m_spDomain->subset_handler()->get_subset_index (fd_elem);
			}
			// otherwise this is no boundary of the domain and we return 0
		}
		else
		{
			for (size_t k = 0; k < fd_elem_list.size (); k++)
			{
				fd_elem_type * e = fd_elem_list[k];
				int e_si = m_ssGrp.subset_handler()->get_subset_index (e);
				if (m_ssGrp.contains (e_si))
				{
					if (fd_elem != NULL) // i.e. e is already the second one
					{ // this is no boundary of the subset; return 0
						fd_elem = NULL;
						break;
					}
					fd_elem = e;
					fd_si = e_si;
				}
			}
		}
		if (fd_elem == NULL)
		{ // this is no bondary, return 0
			for (size_t ip = 0; ip < nip; ip++) vValue[ip] = 0;
			return;
		}
		
	//	Get the coordinates of the corners of the full-dim. element
		std::vector<MathVector<dim> > fd_elem_corner_coords (domain_traits<dim>::MaxNumVerticesOfElem);
		CollectCornerCoordinates (fd_elem->base_object_id (), fd_elem_corner_coords,
			*fd_elem, *m_spDomain, true);
		
	//	Convert the local ip coordinates (use the global coordinates as they are the same)
		const ReferenceObjectID fd_roid = fd_elem->reference_object_id ();
		DimReferenceMapping<dim, dim> & fd_map
			= ReferenceMappingProvider::get<dim, dim> (fd_roid, fd_elem_corner_coords);
		std::vector<MathVector<dim> > fd_loc_ip (nip);
		for (size_t ip = 0; ip < nip; ip++)
			fd_map.global_to_local(fd_loc_ip[ip], vGlobIP[ip]);
	
	//	Call the original UserData, get the vectors
		(* m_spData) (vValue, vGlobIP, time, fd_si, fd_elem, &(fd_elem_corner_coords[0]),
			&(fd_loc_ip[0]), nip, u);
		//TODO: Note that u is here a dummy parameter: We do not have enough data for it.
	
	//	Project the vectors
		const ReferenceObjectID roid = elem->reference_object_id ();
		DimReferenceMapping<refDim, dim> & map
			= ReferenceMappingProvider::get<refDim, dim> (roid, vCornerCoords);
		for (size_t ip = 0; ip < nip; ip++)
		{
			MathMatrix<dim, refDim> J;
			MathVector<dim> p (vValue[ip]);
			
			map.jacobian (J, vLocIP[ip]);
			OrthogProjectVec (p, J);			
			vValue[ip] -= p;
		}
	};
	
///	This function should not be used
	void operator()
	(
		vec_type & vValue,
		const MathVector<dim> & globIP,
		number time,
		int si
	)
	const
	{
		UG_THROW("OutNormCmp: Element required for evaluation, but not passed. Cannot evaluate.");
	}

///	This function should not be used
	void operator()
	(
		vec_type vValue [],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		const size_t nip
	) const
	{
		UG_THROW("OutNormCmp: Element required for evaluation, but not passed. Cannot evaluate.");
	}

};

/// Scaled projection to the outer normal for a given vector user data
/**
 * This class implements the UserData for evaluation of the normal component of vectors
 * computed by the other user data object. The normal component is evaluation only for
 * WDim-1-dimensional subsets - sides of WDim-dimensional elements. It is scaled by
 * a furthe (scalar) UserData.
 *
 * ToDo: Note that the LocalVector provided for the evaluate function is not properly
 * used. Thus, the passed UserData objects should not use it.
 *
 */
template <typename TDomain>
class ScaledOutNormCmp
	: public StdUserData<ScaledOutNormCmp<TDomain>, MathVector<TDomain::dim>, TDomain::dim, void, UserData<MathVector<TDomain::dim>, TDomain::dim, void> >
{
///	the world dimension
	static const int dim = TDomain::dim;
	
///	the domain type
	typedef TDomain domain_type;
	
///	the grid type
	typedef typename TDomain::grid_type grid_type;
	
///	the original data type
	typedef MathVector<dim> vec_type;
	
///	the original vector field data
	SmartPtr<UserData<vec_type, dim, void> > m_spVecData;
	
///	the scaling data
	SmartPtr<UserData<number, dim, void> > m_spScalData;
	
///	the domain
	SmartPtr<domain_type> m_spDomain;
	
///	the subset group for the inner part
	SubsetGroup m_ssGrp;
	
public:

///	constructor
	ScaledOutNormCmp
	(
		SmartPtr<domain_type> spDomain, ///< the domain
		SmartPtr<UserData<number, dim, void> > spScalData, ///< the scaling data
		SmartPtr<UserData<vec_type, dim, void> > spVecData, ///< the original vector field data
		const char * ss_names ///< the subsets
	)
	:	m_spVecData (spVecData), m_spScalData (spScalData), m_spDomain (spDomain), m_ssGrp (spDomain->subset_handler ())
	{
	// Parse the subset names:
		std::vector<std::string> vssNames;
		try
		{
			TokenizeString (ss_names, vssNames);
			for (size_t k = 0; k < vssNames.size (); k++)
				RemoveWhitespaceFromString (vssNames [k]);
			m_ssGrp.clear ();
			m_ssGrp.add (vssNames);
		} UG_CATCH_THROW ("SubsetIndicatorUserData: Failed to parse subset names.");
	}
	
///	constructor
	ScaledOutNormCmp
	(
		SmartPtr<domain_type> spDomain, ///< the domain
		SmartPtr<UserData<number, dim, void> > spScalData, ///< the scaling data
		SmartPtr<UserData<vec_type, dim, void> > spVecData ///< the original vector field data
	)
	:	m_spVecData (spVecData), m_spScalData (spScalData), m_spDomain (spDomain), m_ssGrp (spDomain->subset_handler ())
	{}
	
///	Indicator functions are discontinuous
	virtual bool continuous () const {return false;}

///	Returns true to get the grid element in the evaluation routine
	virtual bool requires_grid_fct () const {return m_spVecData->requires_grid_fct () || m_spScalData->requires_grid_fct ();}
	// Remark: Note that actually we have not enought data for that. This dependence is neglected.

///	Evaluator
	template <int refDim>
	inline void evaluate
	(
		vec_type vValue[],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		GridObject * elem,
		const MathVector<dim> vCornerCoords [],
		const MathVector<refDim> vLocIP [],
		const size_t nip,
		LocalVector * u,
		const MathMatrix<refDim, dim> * vJT = NULL
	) const
	{
		if (refDim == dim)
		{
			std::vector<number> scaling (nip);
			(* m_spVecData) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			(* m_spScalData) (&(scaling[0]), vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
			for (size_t ip = 0; ip < nip; ip++)
				vValue[ip] *= scaling[ip];
			return;
		}
		if (refDim != dim - 1)
		{
			UG_THROW ("ScaledOutNormCmp: Wrong dimensionality of the subset.");
		}
		
	//	Get the full-dim. element associated with the given side
		typedef typename grid_dim_traits<dim>::grid_base_object fd_elem_type;
		fd_elem_type * fd_elem = NULL;
		int fd_si = -1;
		
		typedef typename Grid::traits<fd_elem_type>::secure_container fd_elem_secure_container_type;
		fd_elem_secure_container_type fd_elem_list;
		grid_type& grid = * (grid_type*) (m_spDomain->grid().get ());
		grid.associated_elements (fd_elem_list, elem);
		
		if (m_ssGrp.empty ()) // no subsets specified; we assume, elem should be the bnd of the whole domain
		{
			if (fd_elem_list.size () == 1)
			{
				fd_elem = fd_elem_list[0];
				fd_si = m_spDomain->subset_handler()->get_subset_index (fd_elem);
			}
			// otherwise this is no boundary of the domain and we return 0
		}
		else
		{
			for (size_t k = 0; k < fd_elem_list.size (); k++)
			{
				fd_elem_type * e = fd_elem_list[k];
				int e_si = m_ssGrp.subset_handler()->get_subset_index (e);
				if (m_ssGrp.contains (e_si))
				{
					if (fd_elem != NULL) // i.e. e is already the second one
					{ // this is no boundary of the subset; return 0
						fd_elem = NULL;
						break;
					}
					fd_elem = e;
					fd_si = e_si;
				}
			}
		}
		if (fd_elem == NULL)
		{ // this is no bondary, return 0
			for (size_t ip = 0; ip < nip; ip++) vValue[ip] = 0;
			return;
		}
		
	//	Get the coordinates of the corners of the full-dim. element
		std::vector<MathVector<dim> > fd_elem_corner_coords (domain_traits<dim>::MaxNumVerticesOfElem);
		CollectCornerCoordinates (fd_elem->base_object_id (), fd_elem_corner_coords,
			*fd_elem, *m_spDomain, true);
		
	//	Convert the local ip coordinates (use the global coordinates as they are the same)
		const ReferenceObjectID fd_roid = fd_elem->reference_object_id ();
		DimReferenceMapping<dim, dim> & fd_map
			= ReferenceMappingProvider::get<dim, dim> (fd_roid, fd_elem_corner_coords);
		std::vector<MathVector<dim> > fd_loc_ip (nip);
		for (size_t ip = 0; ip < nip; ip++)
			fd_map.global_to_local(fd_loc_ip[ip], vGlobIP[ip]);
	
	//	Call the original UserData, get the vectors
		std::vector<number> scaling (nip);
		(* m_spVecData) (vValue, vGlobIP, time, fd_si, fd_elem, &(fd_elem_corner_coords[0]),
			&(fd_loc_ip[0]), nip, u);
		(* m_spScalData) (&(scaling[0]), vGlobIP, time, fd_si, fd_elem, &(fd_elem_corner_coords[0]),
			&(fd_loc_ip[0]), nip, u);
		//TODO: Note that u is here a dummy parameter: We do not have enough data for it.
	
	//	Scale and project the vectors
		const ReferenceObjectID roid = elem->reference_object_id ();
		DimReferenceMapping<refDim, dim> & map
			= ReferenceMappingProvider::get<refDim, dim> (roid, vCornerCoords);
		for (size_t ip = 0; ip < nip; ip++)
		{
			vValue[ip] *= scaling[ip];
			
			MathMatrix<dim, refDim> J;
			MathVector<dim> p (vValue[ip]);
			
			map.jacobian (J, vLocIP[ip]);
			OrthogProjectVec (p, J);			
			vValue[ip] -= p;
		}
	};
	
///	This function should not be used
	void operator()
	(
		vec_type & vValue,
		const MathVector<dim> & globIP,
		number time,
		int si
	)
	const
	{
		UG_THROW("ScaledOutNormCmp: Element required for evaluation, but not passed. Cannot evaluate.");
	}

///	This function should not be used
	void operator()
	(
		vec_type vValue [],
		const MathVector<dim> vGlobIP [],
		number time,
		int si,
		const size_t nip
	) const
	{
		UG_THROW("ScaledOutNormCmp: Element required for evaluation, but not passed. Cannot evaluate.");
	}

};

} // end namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__ONC_USER_DATA__

/* End of File */
