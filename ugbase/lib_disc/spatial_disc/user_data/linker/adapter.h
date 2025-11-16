/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Authors: Arne Naegel
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
 * andreasvogel used scale_add_linker ass template
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__LINKER_ADAPTER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__LINKER_ADAPTER__

#include "linker.h"

namespace ug{


template <int dim>
class UserVectorEntryAdapter
	: public StdDataLinker< UserVectorEntryAdapter<dim>, number, dim>
{
public:
	///	Base class type
		using base_type = StdDataLinker< UserVectorEntryAdapter<dim>, number, dim>;

		using data_type = number;
		using user_data_base_type = CplUserData<data_type, dim>;

		using encapsulated_type = MathVector<dim>;
		using input_type = CplUserData<encapsulated_type, dim>;

	public:

		UserVectorEntryAdapter() :  m_index(0), m_spEncaps(nullptr)
		{
			this->set_num_input(_INPUT_+1); //	this linker has one inoput
		}


		inline void evaluate (data_type& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			encapsulated_type dummy;
			(*m_spEncaps)(dummy, globIP, time, si);
			value = dummy[m_index];
		}

		template <int refDim>
		inline void evaluate(data_type vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = nullptr) const
		{
		    std::vector<encapsulated_type> dummy(nip);


		    (*m_spEncaps)(&dummy[0], vGlobIP, time, si,
								elem, vCornerCoords, vLocIP, nip, u, vJT);


			for (size_t ip=0; ip<nip; ++ip)
			{ vValue[ip] = dummy[ip][m_index]; }
		}

		template <int refDim>
		void eval_and_deriv(data_type vValue[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GridObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<data_type > > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = nullptr) const
		{
			//	get the data of the ip series
			const encapsulated_type* vDummy = m_spEncaps->values(s);

			for(size_t ip = 0; ip < nip; ++ip)
			{ vValue[ip] = vDummy[ip][m_index]; }

			//	check if something to do
		 	if(!bDeriv || this->zero_derivative()) return;

		 	//	clear all derivative values
			this->set_zero(vvvDeriv, nip);

			//	Derivatives w.r.t  input
			if( m_spDEncaps.invalid() || m_spDEncaps->zero_derivative()) return;

			// loop integration points
			for(size_t ip = 0; ip < nip; ++ip){

				// loop functions
				for(size_t fct = 0; fct < m_spDEncaps->num_fct(); ++fct)
				{
					//	get derivative of  w.r.t. to all functions
					const encapsulated_type* vDInputFct = m_spDEncaps->deriv(s, ip, fct);

					//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_INPUT_, fct);

					//	loop all shapes and set the derivative
					for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
					{
						UG_ASSERT(commonFct < vvvDeriv[ip].size(), commonFct<<", "<<vvvDeriv[ip].size());
						vvvDeriv[ip][commonFct][sh] += vDInputFct[sh][m_index]; // TODO: double check this!
					}

				}
			}

		}

	public:
	///	set conc import
		void set_vector(SmartPtr<input_type> data, size_t index)
		{
			m_spEncaps = data;  // for evaluation
			m_spDEncaps = data.template cast_dynamic<DependentUserData<encapsulated_type, dim> >(); // for derivatives
			base_type::set_input(_INPUT_, data, data);
			m_index = index;
		}



	protected:
		size_t m_index;

		static constexpr int _INPUT_= 0;

	///	import for concentration
		SmartPtr<input_type> m_spEncaps;
		SmartPtr<DependentUserData<encapsulated_type, dim> > m_spDEncaps;


};

} // namespace ug

#endif /*  __H__UG__LIB_DISC__SPATIAL_DISC__LINKER_ADAPTER__ */
