/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DATA_LINKER__

#include "../std_user_data.h"
// #include "lib_disc/common/groups_util.h"

namespace ug {

/// combines several UserDatas to a new UserData of a specified type
/**
 * This class provides data at integration points and implements the
 * DependentUserData interface.
 *
 * \tparam 	TData		output Data type
 * \tparam 	dim			World dimension
 */
template <typename TImpl, typename TData, int dim>
class StdDataLinker
	: public StdUserData< StdDataLinker<TImpl, TData, dim>,
						  TData, dim, void,
						  DependentUserData<TData, dim> >
{
	public:
		void operator () (TData& value,
	                  const MathVector<dim>& globIP,
	                  number time, int si) const override;

		void operator () (TData vValue[],
		                  const MathVector<dim> vGlobIP[],
		                  number time, int si, size_t nip) const override;

		template <int refDim>
		inline void evaluate(TData vValue[],
							 const MathVector<dim> vGlobIP[],
							 number time, int si,
							 GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 const MathVector<refDim> vLocIP[],
							 size_t nip,
							 LocalVector* u,
							 const MathMatrix<refDim, dim>* vJT = nullptr) const;

		void compute(LocalVector* u, GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 bool bDeriv = false) override;

		void compute(LocalVectorTimeSeries* u, GridObject* elem,
							 const MathVector<dim> vCornerCoords[],
							 bool bDeriv = false) override;

	protected:
		template <int refDim>
		void eval_deriv(LocalVector* u, GridObject* elem,
		                const MathVector<dim> vCornerCoords[], bool bDeriv = false);
		template <int refDim>
		void eval_deriv(LocalVectorTimeSeries* u, GridObject* elem,
		                const MathVector<dim> vCornerCoords[], bool bDeriv = false);

	public:
	///	returns that a grid function is needed for evaluation
		[[nodiscard]] bool requires_grid_fct() const override;

	///	returns if provided data is continuous over geometric object boundaries
		[[nodiscard]] bool continuous() const override;

	///	returns if derivative is zero
		[[nodiscard]] bool zero_derivative() const override;

	public:
	///	returns if the derivative of the i'th input is zero
		[[nodiscard]] bool zero_derivative(size_t i) const
		{
			if(!m_vspICplUserData[i].valid()) return true;
			return m_vspICplUserData[i]->zero_derivative();
		}

	///	sets the number of inputs
		void set_num_input(size_t num)
		{
			m_vspICplUserData.resize(num);
			m_vspUserDataInfo.resize(num);
		}

	///	sets an input
		virtual void set_input(size_t i,
		                       SmartPtr<ICplUserData<dim> > input,
		                       SmartPtr<UserDataInfo> info)
		{
			UG_ASSERT(i < m_vspICplUserData.size(), "invalid index");
			m_vspICplUserData[i] = input;
			m_vspUserDataInfo[i] = info;
		}

	///	number of inputs
		[[nodiscard]] virtual size_t num_input() const {return num_needed_data();}

	///	number of other Data this data depends on
		[[nodiscard]] size_t num_needed_data() const override {return m_vspICplUserData.size();}

	///	return needed data
		SmartPtr<ICplUserData<dim> > needed_data(size_t i) override
		{
			UG_ASSERT(i < m_vspICplUserData.size(), "Input not needed");
			UG_ASSERT(m_vspICplUserData[i].valid(), "Data input not valid");
			return m_vspICplUserData[i];
		}

	///	returns if data is ok
		void check_setup() const override;

	///	updates the function group
		void set_function_pattern(ConstSmartPtr<FunctionPattern> fctPatt) override;

	protected:
	///	returns number of functions the input depends on
		[[nodiscard]] size_t input_num_fct(size_t i) const
		{
			UG_ASSERT(i < m_vspUserDataInfo.size(), "Input invalid");
			if(!m_vspUserDataInfo[i].valid()) return 0;
			return m_vspUserDataInfo[i]->num_fct();
		}

	///	returns the number in the common FctGrp for a fct of an input
		[[nodiscard]] size_t input_common_fct(size_t i, size_t fct) const
		{
			UG_ASSERT(i < m_vMap.size(), "Input Map invalid");
			UG_ASSERT(fct < m_vMap[i].num_fct(), "Input Map invalid for fct");
			return m_vMap[i][fct];
		}

	///	returns the series id set for the i'th input
		[[nodiscard]] size_t series_id(size_t i, size_t s) const
		{
			UG_ASSERT(i < m_vvSeriesID.size(), "invalid index");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "invalid index");
			return m_vvSeriesID[i][s];
		}

	///	requests series id's from input data
		void local_ip_series_added(size_t seriesID) override;

	///	forwards the local positions to the data inputs
		void local_ips_changed(size_t seriesID, size_t newNumIP) override;

	///	requests cleaning of the ip series in the data inputs
		void local_ip_series_to_be_cleared() override;

	///	forwards the global positions to the data inputs
		void global_ips_changed(size_t seriesID, const MathVector<dim>* vPos, size_t numIP) override;

	protected:
	///	data input
		std::vector<SmartPtr<ICplUserData<dim> > > m_vspICplUserData;

	///	data input casted to IDependend data
		std::vector<SmartPtr<UserDataInfo> > m_vspUserDataInfo;

	///	Function mapping for each input relative to common FunctionGroup
		std::vector<FunctionIndexMapping> m_vMap;

	///	series id the linker uses to get data from input
		std::vector<std::vector<size_t> > m_vvSeriesID;

	protected:
	///	access to implementation
		TImpl& getImpl() {return static_cast<TImpl&>(*this);}

	///	const access to implementation
		const TImpl& getImpl() const {return static_cast<const TImpl&>(*this);}
};

} // end namespace ug

//	include implementation
#include "linker_impl.h"

#endif