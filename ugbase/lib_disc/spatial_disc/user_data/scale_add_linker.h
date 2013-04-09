/*
 * scale_add_linker.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER__

#include "lib_disc/spatial_disc/user_data/std/std_linker_data.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Scaled adding of Data
////////////////////////////////////////////////////////////////////////////////

/**
 * This linker recombines the data like
 *
 * l = s_1 * i_1 + ... + s_k * i_k
 *
 * where, l, i_1, ..., i_k are of the same data type and the s_i are some
 * scaling factor of a (possibly) different type, that is applyable to the
 * data type
 *
 * \tparam		TData		exported and combined Data type
 * \tparam		dim			world dimension
 * \tparam		TDataScale 	type of scaling data
 */
template <typename TData, int dim, typename TDataScale>
class ScaleAddLinker
	: public StdDataLinker<ScaleAddLinker<TData, dim, TDataScale>, TData, dim>
{
	public:
	//	type of base class
		typedef DataLinker<TData, dim> base_type;

	public:
	///	constructor
		ScaleAddLinker() {}

	///	constructor
		ScaleAddLinker(const ScaleAddLinker& linker);

	///	adds an input to the list of summands scaled by a user data factor
	///	\{
		void add(SmartPtr<UserData<TDataScale, dim> > scale,
		         SmartPtr<UserData<TData, dim> > data);
		void add(number scale,
		         SmartPtr<UserData<TData, dim> > data);
		void add(SmartPtr<UserData<TDataScale, dim> > scale,
		         number data);
		void add(number scale,
		         number data);
	/// \}

	///	computes the value
		virtual void compute(LocalVector* u, GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false);

		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const;

		template <int refDim>
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const;

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const;

	protected:
	///	data at ip of input
		const TData& input_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpUserData.size(), "Input not needed");
			UG_ASSERT(m_vpUserData[i].valid(), "Input invalid");
			return m_vpUserData[i]->value(this->series_id(2*i,s), ip);
		}

	///	derivative of data at input at ip
		const TData& input_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");
			UG_ASSERT(m_vpDependData[i].valid(), "Input invalid");
			return m_vpDependData[i]->deriv(this->series_id(2*i,s), ip, fct, dof);
		}

	///	scale at ip of input
		const TDataScale& scale_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpScaleData.size(), "Input not needed");
			UG_ASSERT(m_vpScaleData[i].valid(), "Input invalid");
			return m_vpScaleData[i]->value(this->series_id(2*i+1,s), ip);
		}

	///	derivative of data at input at ip
		const TDataScale& scale_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(i < m_vpScaleDependData.size(), "Input not needed");
			UG_ASSERT(m_vpScaleDependData[i].valid(), "Input invalid");
			return m_vpScaleDependData[i]->deriv(this->series_id(2*i+1,s), ip, fct, dof);
		}

	///	returns number of functions the input depends on
		size_t input_num_fct(size_t i) const {return base_type::input_num_fct(2*i);}

	///	returns the number in the common FctGrp for a fct of an input
		size_t input_common_fct(size_t i, size_t fct) const	{return base_type::input_common_fct(2*i, fct);}

	///	returns number of functions the scaling depends on
		size_t scale_num_fct(size_t i) const {return base_type::input_num_fct(2*i+1);}

	///	returns the number in the common FctGrp for a fct of a scaling
		size_t scale_common_fct(size_t i, size_t fct) const	{return base_type::input_common_fct(2*i+1, fct);}

	protected:
	///	data input
		std::vector<SmartPtr<UserData<TDataScale, dim> > > m_vpScaleData;

	///	data input casted to dependend data
		std::vector<SmartPtr<DependentUserData<TDataScale, dim> > > m_vpScaleDependData;

	///	data input
		std::vector<SmartPtr<UserData<TData, dim> > > m_vpUserData;

	///	data input casted to dependend data
		std::vector<SmartPtr<DependentUserData<TData, dim> > > m_vpDependData;
};

} // end namespace ug

#include "scale_add_linker_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER__ */
