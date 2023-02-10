/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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
 * A linker that projects a given vector to the manifold of given elements.
 */
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__PROJECT_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__PROJECT_LINKER__

/* ug4 headers */
#include "common/common.h"

#include "lib_disc/reference_element/reference_mapping_provider.h"
#include "linker.h"

namespace ug {

/**
 * User data projecting a given vector to the plane of degenerated elements.
 * The projection is computed only in the degenerated elements. For other elements,
 * 0 is returned.
 */
template <int dim>
class ProjectionLinker
: public StdDataLinker<ProjectionLinker<dim>, MathVector<dim>, dim>
{
///	Base class type
	typedef StdDataLinker<ProjectionLinker<dim>, MathVector<dim>, dim> base_type;
	
public:

///	Constructor
	ProjectionLinker
	(
		SmartPtr<CplUserData<MathVector<dim>, dim> > spVector ///< vector to project
	)
	{
		this->set_num_input (1);
		m_spVector = spVector;
		m_spDVector = spVector.template cast_dynamic<DependentUserData<MathVector<dim>, dim> > ();
		this->set_input (0, spVector, spVector);
	}
	
///	Constructor
	ProjectionLinker
	(
		MathVector<dim> vector ///< vector to project
	)
	{
		this->set_num_input (1);
		m_spVector = make_sp (new ConstUserVector<dim> (vector));
		m_spDVector = m_spVector.template cast_dynamic<DependentUserData<MathVector<dim>, dim> > ();
		this->set_input (0, m_spVector, m_spVector);
	}
	
///	Returns true because without a grid function, we do not get the element to project to!
	virtual bool requires_grid_fct() const {return true;}

///	Evaluation with no element is impossible
	inline void evaluate
	(
		MathVector<dim>& value,
		const MathVector<dim>& glob_ip,
		number time,
		int si
	) const
	{
		UG_THROW ("ProjectionLinker: Cannot evaluate without any specification of the element!");
	}

///	Computation only of the projection
	template <int refDim>
	inline void evaluate
	(
		MathVector<dim> vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
	//	1. Get the vectors themselves
		(*m_spVector) (vValue, vGlobIP, time, si, elem, vCornerCoords, vLocIP, nip, u, vJT);
		
	//	2. Prepare the data of the element
		const ReferenceObjectID roid = elem->reference_object_id ();
		DimReferenceMapping<refDim, dim>& rMapping = ReferenceMappingProvider::get<refDim, dim> (roid);
		rMapping.update (vCornerCoords);
		
	//	3. At every integration point
		for (size_t i = 0; i < nip; i++)
		{
		//	3a. Get the Jacobian
			MathMatrix<dim, refDim> J;
			rMapping.jacobian (J, vLocIP[i]);
		//	3b. Project the vector to the subspace spanned by the columns of the Jacobian
			OrthogProjectVec (vValue[i], J);
		}
	}

///	Computation of the projection and its derivatives
	template <int refDim>
	void eval_and_deriv
	(
		MathVector<dim> vValue[],
		const MathVector<dim> vGlobIP[],
		number time,
		int si,
		GridObject* elem,
		const MathVector<dim> vCornerCoords[],
		const MathVector<refDim> vLocIP[],
		const size_t nip,
		LocalVector* u,
		bool bDeriv,
		int s,
		std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		const MathMatrix<refDim, dim>* vJT = NULL
	) const
	{
	//	1. Get the vectors themselves
		const MathVector<dim>* vVector = m_spVector->values (s);
		
	//	2a. Prepare the data of the element
		const ReferenceObjectID roid = elem->reference_object_id ();
		DimReferenceMapping<refDim, dim>& rMapping = ReferenceMappingProvider::get<refDim, dim> (roid);
		rMapping.update (vCornerCoords);
		
	//	2b. Check if we should compute the derivatives
		if (this->zero_derivative ())
			bDeriv = false;
		else
			this->set_zero (vvvDeriv, nip);
		
	//	3. At every integration point
		for (size_t i = 0; i < nip; i++)
		{
		//	3a. Get the Jacobian
			MathMatrix<dim, refDim> J;
			rMapping.jacobian (J, vLocIP[i]);
		//	3b. Project the vector to the subspace spanned by the columns of the Jacobian
			vValue[i] = vVector[i];
			OrthogProjectVec (vValue[i], J);
		//	3c. Project the derivatives in the same way
			if (! bDeriv) continue;
			for (size_t fct = 0; fct < m_spDVector->num_fct(); fct++)
			{
				const MathVector<dim>* vDVector = m_spDVector->deriv (s, i, fct);
				const size_t c_fct = this->input_common_fct (0, fct);
				for (size_t sh = 0; sh < this->num_sh (c_fct); sh++)
				{
					vvvDeriv[i][c_fct][sh] = vDVector [sh];
					OrthogProjectVec (vvvDeriv[i][c_fct][sh], J);
				}
			}
		}
	}
	
private:

///	data to project
	SmartPtr<CplUserData<MathVector<dim>, dim> > m_spVector;
	SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDVector;
};

} // end namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__PROJECT_LINKER__

/* End of File */
