/*
 * scale_add_linker.h
 *
 *  Created on: 04.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER__

#include "lib_disc/spatial_disc/ip_data/std_ip_data.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Linker Traits to allow generic programming
////////////////////////////////////////////////////////////////////////////////

/// Linker Traits
template <typename TData, typename TDataIn>
struct linker_traits
{
///	computes out += s * in1 (with appropriate '*')
	static void mult_add(TData& out, const TData& in1, const TDataIn& s);
};

template <>
struct linker_traits<number, number>
{
	static void mult_add(number& out, const number& in1, const number& s)
	{
		out += in1 * s;
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>, number >
{
	static void mult_add(MathVector<dim>& out,
	                     const MathVector<dim>& in1,
						 const number& s)
	{
		VecScaleAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathVector<dim>, MathMatrix<dim,dim> >
{
	static void mult_add(MathVector<dim>& out,
	                     const MathVector<dim>& in1,
						 const MathMatrix<dim,dim>& s)
	{
		MatVecMultAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathMatrix<dim,dim>, number >
{
	static void mult_add(MathMatrix<dim,dim>& out,
	                     const MathMatrix<dim,dim>& in1,
						 const number& s)
	{
		MatScaleAppend(out, s, in1);
	}
};

template <std::size_t dim>
struct linker_traits< MathTensor<4,dim>, number >
{
	static void mult_add(MathTensor<4,dim>& out,
	                     const MathTensor<4, dim>& in1,
						 const number& s)
	{
		out = in1;
		throw(UGError("linker_traits for MathTensor4 not implemented"));
	}
};

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
	: public StdDataLinker<TData, dim, ScaleAddLinker<TData, dim, TDataScale> >
{
	public:
	//	type of base class
		typedef DataLinker<TData, dim> base_type;

	public:
	///	constructor
		ScaleAddLinker() {}

	///	constructor
		ScaleAddLinker(const ScaleAddLinker& linker)
		{
			if(linker.m_vpIPData.size() != linker.m_vpScaleData.size())
				UG_THROW("ScaleAddLinker: number of scaling factors and data mismatch.");

			for(size_t i = 0; i < linker.m_vpIPData.size(); ++i)
			{
				this->add(linker.m_vpScaleData[i], linker.m_vpIPData[i]);
			}
		}

	///	adds an input to the list of summands scaled by a user data factor
	///	\{
		void add(SmartPtr<IPData<TDataScale, dim> > scale,
		         SmartPtr<IPData<TData, dim> > data);
		void add(number scale,
		         SmartPtr<IPData<TData, dim> > data);
		void add(SmartPtr<IPData<TDataScale, dim> > scale,
		         number data);
		void add(number scale,
		         number data);
	/// \}

	///	computes the value
		virtual void compute(bool bDeriv);

		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			//	reset value
			value = 0.0;

			TData valData;
			TDataScale valScale;

			//	add contribution of each summand
				for(size_t c = 0; c < m_vpIPData.size(); ++c)
				{
					(*m_vpIPData[c])(valData, globIP, time, si);
					(*m_vpScaleData[c])(valScale, globIP, time, si);

					linker_traits<TData, TDataScale>::
					mult_add(value, valData, valScale);
				}
		}

		template <int refDim>
		inline void evaluate (TData& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
			//	reset value
			value = 0.0;

			TData valData;
			TDataScale valScale;

			//	add contribution of each summand
				for(size_t c = 0; c < m_vpIPData.size(); ++c)
				{
					(*m_vpIPData[c])(valData, globIP, time, si, u, elem, vCornerCoords, locIP);
					(*m_vpScaleData[c])(valScale, globIP, time, si, u, elem, vCornerCoords, locIP);

					linker_traits<TData, TDataScale>::
					mult_add(value, valData, valScale);
				}
		}

		template <int refDim>
		inline void evaluate(TData vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			//	reset value
			for(size_t ip = 0; ip < nip; ++ip)
				vValue[ip] = 0.0;

			std::vector<TData> vValData(nip);
			std::vector<TDataScale> vValScale(nip);

			//	add contribution of each summand
				for(size_t c = 0; c < m_vpIPData.size(); ++c)
				{
					(*m_vpIPData[c])(&vValData[0], vGlobIP, time, si, u,
									elem, vCornerCoords, vLocIP, nip, vJT);
					(*m_vpScaleData[c])(&vValScale[0], vGlobIP, time, si, u,
										elem, vCornerCoords, vLocIP, nip, vJT);

					for(size_t ip = 0; ip < nip; ++ip)
						linker_traits<TData, TDataScale>::
						mult_add(vValue[ip], vValData[ip], vValScale[ip]);
				}
		}

	protected:
	///	data at ip of input
		const TData& input_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i].valid(), "Input invalid");
			return m_vpIPData[i]->value(this->series_id(2*i,s), ip);
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
		size_t input_num_fct(size_t i) const
		{
			return base_type::input_num_fct(2*i);
		}

	///	returns the number in the common FctGrp for a fct of an input
		size_t input_common_fct(size_t i, size_t fct) const
		{
			return base_type::input_common_fct(2*i, fct);
		}

	///	returns number of functions the scaling depends on
		size_t scale_num_fct(size_t i) const
		{
			return base_type::input_num_fct(2*i+1);
		}

	///	returns the number in the common FctGrp for a fct of a scaling
		size_t scale_common_fct(size_t i, size_t fct) const
		{
			return base_type::input_common_fct(2*i+1, fct);
		}

	protected:
	///	data input
		std::vector<SmartPtr<IPData<TDataScale, dim> > > m_vpScaleData;

	///	data input casted to dependend data
		std::vector<SmartPtr<DependentIPData<TDataScale, dim> > > m_vpScaleDependData;

	///	data input
		std::vector<SmartPtr<IPData<TData, dim> > > m_vpIPData;

	///	data input casted to dependend data
		std::vector<SmartPtr<DependentIPData<TData, dim> > > m_vpDependData;
};

} // end namespace ug

#include "scale_add_linker_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__SCALE_ADD_LINKER__ */
