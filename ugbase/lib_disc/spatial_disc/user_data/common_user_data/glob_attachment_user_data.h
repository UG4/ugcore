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


namespace ug{

template <int dim>
class GlobAttachmentElementUserData
: public StdDependentUserData<GlobAttachmentElementUserData<dim>, number, dim>
{

	private:
		std::string m_attachment_name;
		SmartPtr<Grid> m_spGrid;
		
	public:
		/// constructor
		GlobAttachmentElementUserData(SmartPtr<Grid> grid, const char* name)
		: m_attachment_name(name), m_spGrid(grid)
		{
			// ANumber att = GlobalAttachments::attachment<ANumber>(m_attachment_name);
			// Grid::VolumeAttachmentAccessor<ANumber> aatt((Grid&)*m_spGrid, att);
		};

		virtual bool continuous() const
		{
			return false;
		}

		template <int refDim>
		void eval_and_deriv(number vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<number> > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = NULL) const
		{
				if (refDim != 3)
				{
					UG_THROW ("GlobAttachmentElementUserData: Diamention of the element is not equal to 3.");
				}
			
				ANumber att = GlobalAttachments::attachment<ANumber>(m_attachment_name);
				Grid::VolumeAttachmentAccessor<ANumber> aatt((Grid&)*m_spGrid, att);
		
				//	loop ips
				
				for(size_t ip = 0; ip < nip; ++ip)
				{
					vValue[ip] = aatt[(Volume*)elem];
					
				}

				if(bDeriv){
					UG_THROW ("GlobAttachmentElementUserData: Derivative not implemented.");
					}

		}
	
		virtual void operator() (number& value,
									const MathVector<dim>& globIP,
									number time, int si,
									Vertex* vrt) const
		{
			UG_THROW ("GlobAttachmentElementUserData: Values at vertices of the grid function are not uniquely defined.");
		}
};


} // end namespace ug

#endif /* __H__UG__LIB_DISC__FUNCTION_SPACE__GRID_FUNCTION_USER_DATA__ */
