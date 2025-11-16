/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Susanne Höllbacher
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

#ifndef ELEM_MODIFIER_H_
#define ELEM_MODIFIER_H_

#include "elem_disc_interface.h"

namespace ug{

template <typename TDomain>
class IElemDisc;


template <typename TDomain>
class IElemDiscModifier
{

	protected:
	///	own type
		using this_type = IElemDiscModifier<TDomain>;

	public:
	///	World dimension
		static constexpr int dim = TDomain::dim;

	public:
	///	Constructor (setting default values)
	/// \{
		IElemDiscModifier(): m_pElemDisc(nullptr){};
		IElemDiscModifier(IElemDisc<TDomain>* myElemDisc) : m_pElemDisc(myElemDisc) {};
	 /// \}

	/// Virtual destructor
		virtual ~IElemDiscModifier()= default;


	/// virtual initiates pre-computations before the standard element assembling
		virtual void preprocess(LocalVector& u, LocalVector& d, LocalVector& tmpD, GridObject* elem,
								MathVector<dim> vCornerCoords[], LocalIndices& ind);

	/// virtual initiates pre-computations before the standard element assembling
		virtual void preprocess(LocalVector& u, LocalMatrix& J, GridObject* elem,
								MathVector<dim> vCornerCoords[], LocalIndices& ind);

	/// virtual initiates post-computations after the standard element assembling
		virtual void postprocess(const LocalVector& u, LocalVector& d, LocalIndices& ind);

	/// virtual initiates post-computations after the standard element assembling
		virtual void postprocess(const LocalVector& u, LocalMatrix& J, LocalIndices& ind);

        void set_elem_disc(IElemDisc<TDomain>* myElemDisc){ m_pElemDisc = myElemDisc; }


      protected:
        IElemDisc<TDomain>* m_pElemDisc;
};

/*
 IElemDiscModifier_Local : IElemDiscModifier

 u.acces_by_map(...);
*/
} // end name space ug

#include "elem_modifier_impl.h"


#endif /* ELEM_MODIFIER_H_ */
