/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Shuai Lu
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

template <int WDim, typename TData = number>
class GlobAttachmentElementUserData
: public StdUserData<GlobAttachmentElementUserData<WDim, TData>, TData, WDim>

{
		static constexpr int dim = WDim;
		using data_type = TData;
		using elem_type = typename grid_dim_traits<dim>::grid_base_object;
		using attachment_type = Attachment<data_type>;

	private:
		std::string m_attachment_name;
		SmartPtr<Grid> m_spGrid;
		attachment_type m_att;
		Grid::AttachmentAccessor<elem_type, attachment_type> m_aatt;
		
	///	Evaluation of the attachment in one element
		inline void eval_on_elem
		(
			elem_type * elem, ///< the element to evaluate on
			const size_t nip, ///< number of the values to write
			data_type vValue []  ///< array for the values
		) const
		{
			for (size_t ip = 0; ip < nip; ++ip)
				vValue[ip] = m_aatt [elem];
		}

	public:
	
	/// constructor
		GlobAttachmentElementUserData(SmartPtr<Grid> grid, const char* name)
		:	m_attachment_name(name), m_spGrid(grid)
		{
			if (! GlobalAttachments::is_declared (m_attachment_name))
				UG_THROW ("GlobAttachmentElementUserData: No global attachment '" << m_attachment_name << "' found.");
			m_att = GlobalAttachments::attachment<attachment_type> (m_attachment_name);
			m_aatt.access (*grid, m_att);
		};
	
	//	UserData interface
		
		virtual bool continuous () const
		{
			return false;
		}

		virtual bool requires_grid_fct () const
		{
			return true;
		}

	//	StdUserData interface
		template <int refDim>
		inline void evaluate
		(
			data_type vValue[],
			const MathVector<dim> vGlobIP[],
			number time, int si,
			GridObject* elem,
			const MathVector<dim> vCornerCoords[],
			const MathVector<refDim> vLocIP[],
			const size_t nip,
			LocalVector* u,
			const MathMatrix<refDim, dim>* vJT = nullptr
		) const
		{
			UG_ASSERT (refDim == dim, "GlobAttachmentElementUserData: Dimensionality of the element should be equal to the world dimensionality.");
			
			eval_on_elem ((elem_type *) elem, nip, vValue);
		}
	
	//	StdUserData interface
		
		virtual void compute
		(
			LocalVector* u,
			GridObject* elem,
			const MathVector<dim> vCornerCoords[],
			bool bDeriv = false
		)
		{
			UG_ASSERT (elem->base_object_id() == dim, "GlobAttachmentElementUserData: Dimensionality of the element should be equal to the world dimensionality.");
			
			for (size_t s = 0; s < this->num_series (); ++s)
				eval_on_elem ((elem_type *) elem, this->num_ip (s), this->values (s));
		}

		virtual void compute
		(
			LocalVectorTimeSeries* u,
			GridObject* elem,
			const MathVector<dim> vCornerCoords[],
			bool bDeriv = false
		)
		{
			UG_ASSERT (elem->base_object_id() == dim, "GlobAttachmentElementUserData: Dimensionality of the element should be equal to the world dimensionality.");
			
			for (size_t s = 0; s < this->num_series (); ++s)
				eval_on_elem ((elem_type *) elem, this->num_ip (s), this->values (s));
		}

		virtual void operator () ///< cannot be implemented here
		(
			data_type& value,
			const MathVector<dim>& globIP,
			number time, int si
		) const
		{
			UG_THROW("GlobAttachmentElementUserData: Element required"
					 " for evaluation, but not passed. Cannot evaluate.");
		}

		virtual void operator () ///< cannot be implemented here
		(
			data_type vValue[],
			const MathVector<dim> vGlobIP[],
			number time, int si,
			const size_t nip
		) const
		{
			UG_THROW("GlobAttachmentElementUserData: Element required"
					 " for evaluation, but not passed. Cannot evaluate.");
		}
};


} // end namespace ug

#endif