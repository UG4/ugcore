/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__
#define __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__

#include "common/common.h"

#include "lib_disc/spatial_disc/user_data/std_user_data.h"
#include "lib_grid/global_attachments.h"
#include "lib_grid/grid_objects/grid_dim_traits.h"


namespace ug{

/**
 * This class provides an interface between data saved in a grid attachment via the
 * coupling interface. Its main purpuse is to access the values saved in the (global,
 * named) attachments in discretizations. The data should be attached to the full-dim.
 * grid objects (grid elements), for ex. to volumens in 3d grids.
 */
template <typename TDomain>
class GlobAttachmentElementUserData
: public StdUserData<GlobAttachmentElementUserData<TDomain>, number, TDomain::dim>
{
	static const int dim = TDomain::dim; ///< the world dimension
	typedef typename grid_dim_traits<dim>::grid_base_object elem_t; ///< type of the full-dim. grid objects

	std::string m_attachment_name; ///< name of the global attachment
	SmartPtr<Grid> m_spGrid; ///< grid of the attachment
	ANumber m_att; ///< the attachment
	Grid::AttachmentAccessor<elem_t, ANumber> m_aatt; ///< the attachment accessor
	
protected:
///	gets the attachment value
	number get_value
	(
		GridObject* elem ///< pointer to the element where evaluate
	) const
	{
		return m_aatt [(elem_t *)elem];
	}
		
public:

/// constructor
	GlobAttachmentElementUserData
	(
		SmartPtr<Grid> grid, ///< grid
		const char* name ///< registered name of the global attachment
	)
	: m_attachment_name(name), m_spGrid(grid)
	{
		m_att = GlobalAttachments::attachment<ANumber>(m_attachment_name);
		m_aatt.access((Grid &) (*m_spGrid), m_att);
	};

	virtual bool continuous() const {return false;} // the data are piecewise constant
	virtual bool requires_grid_fct() const {return true;} // to always get the grid element!
	virtual bool constant() const {return false;} // the data are not globally constant
	virtual bool zero_derivative() const {return true;} // we do not consider any derivatives

/// evaluator for StdUserData interface
	template <int refDim>
	inline void evaluate
	(
		number vValue [],
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
		if (refDim != dim)
			UG_THROW ("GlobAttachmentElementUserData: Only evaluation in full-dim. elements is supported.");
		const number val = get_value(elem);
		for (size_t ip = 0; ip < nip; ++ip) vValue[ip] = val;
	}
	
///	function from the CplUserData interface
	virtual void compute
	(
		LocalVector* u,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		bool bDeriv ///< assumed to be always false
	)
	{
		const number val = get_value(elem);
		for(size_t s = 0; s < this->num_series(); ++s)
			for(size_t ip = 0; ip < this->num_ip(s); ++ip)
				this->value(s,ip) = val;
	}

///	function from the CplUserData interface
	virtual void compute
	(
		LocalVectorTimeSeries* u,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		bool bDeriv ///< assumed to be always false
	)
	{
		const number val = get_value(elem);
		for(size_t s = 0; s < this->num_series(); ++s)
			for(size_t ip = 0; ip < this->num_ip(s); ++ip)
				this->value(s,ip) = val;
	}
	
//	The following two operators may not be used: They do not get the grid element
	virtual void operator()
	(
		number& vValue,
		const MathVector<dim>& globIP,
		number time,
		int si
	) const
	{
		UG_THROW ("GlobAttachmentElementUserData: Cannot provide values basing only on the global coordinates.");
	}
	virtual void operator()
	(
		number vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		const size_t nip
	) const
	{
		UG_THROW ("GlobAttachmentElementUserData: Cannot provide values basing only on the global coordinates.");
	}
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__ */
