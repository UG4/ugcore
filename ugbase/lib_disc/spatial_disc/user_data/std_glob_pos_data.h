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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__STD_GLOB_POS_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__STD_GLOB_POS_DATA__

#include "std_user_data.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Base class for Position-Time-UserData
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a base class for all position and time dependent user data.
 * The data thus does not depend on the computed solution.
 * In order to use the interface, the deriving class must implement the method:
 *
 * inline TRet evaluate(TData& D, const MathVector<dim>& x, number time, int si) const
 *
 */
template <typename TImpl, typename TData, int dim, typename TRet = void>
class StdGlobPosData
	: 	public StdUserData<StdGlobPosData<TImpl,TData,dim,TRet>, TData, dim, TRet>
{
	public:
		virtual TRet operator () (TData& value,
								 const MathVector<dim>& globIP,
								 number time, int si) const
		{
			return this->getImpl().evaluate(value, globIP, time, si);
		}

		virtual void operator () (TData vValue[],
								const MathVector<dim> vGlobIP[],
								number time, int si, const size_t nip) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				this->getImpl().evaluate(vValue[ip], vGlobIP[ip], time, si);
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = nullptr) const
		{
			for(size_t ip = 0; ip < nip; ++ip)
				this->getImpl().evaluate(vValue[ip],vGlobIP[ip],time,si);
		}

	///	implement as a UserData
		virtual void compute(LocalVector* u, GridObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
			const number t = this->time();
			const int si = this->subset();

			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					this->getImpl().evaluate(this->value(s,ip), this->ip(s, ip), t, si);
		}

	///	implement as a UserData
		virtual void compute(LocalVectorTimeSeries* u, GridObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
			const int si = this->subset();

			for(size_t s = 0; s < this->num_series(); ++s)
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
					this->getImpl().evaluate(this->value(s,ip), this->ip(s, ip), this->time(s), si);
		}

	///	returns if data is constant
		virtual bool constant() const {return false;}

	///	returns if grid function is needed for evaluation
		virtual bool requires_grid_fct() const {return false;}

	///	returns if provided data is continuous over geometric object boundaries
		virtual bool continuous() const {return true;}

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

} // namespace ug

#endif