/*
 * data_linker.h
 *
 *  Created on: 12.11.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__

#include "ip_data.h"
#include "lib_discretization/common/groups_util.h"

namespace ug{

/// combines several IPDatas to a new IPData of a specified type
/**
 * This class provides data at integration points and implements the
 * DependentIPData interface.
 *
 * \tparam 	TData		output Data type
 * \tparam 	dim			World dimension
 */
template <typename TData, int dim>
class DataLinker
	: public DependentIPData<TData, dim>
{
	public:
	///	Base class type
		typedef DependentIPData<TData, dim> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::time;

	//	explicitly forward methods of IPData
		using base_type::ip;
		using base_type::value;

	public:
	///	constructor
		DataLinker()
		{
			m_vpIIPData.clear();
			m_vpIDependData.clear();
		}

	///	compute method
		virtual void compute(bool compDeriv) = 0;

	///	returns if derivative is zero
		virtual bool zero_derivative() const
		{
			bool bRet = true;

		//	loop inputs
			for(size_t i = 0; i < m_vpIIPData.size(); ++i)
			{
			//	skip unset data ( null as default )
				if(m_vpIIPData[i] == NULL) continue;

			//	flag iff ipdata is dependent
				bRet = bRet && m_vpIIPData[i]->zero_derivative();
			}

			return bRet;
		}

	///	returns if the derivative of the i'th input is zero
		bool zero_derivative(size_t i) const
		{
			if(m_vpIIPData[i] == NULL) return true;
			return m_vpIIPData[i]->zero_derivative();
		}

	///	sets the number of inputs
		void set_num_input(size_t num)
		{
			m_vpIIPData.resize(num, NULL);
			m_vpIDependData.resize(num, NULL);
		}

	///	sets an input
		void set_input(size_t i, IIPDimData<dim>* input)
		{
			UG_ASSERT(i < m_vpIIPData.size(), "invalid index");
			m_vpIIPData[i] = input;
			m_vpIDependData[i] = dynamic_cast<IDependentIPData*>(input);
		}

	///	number of inputs
		size_t num_input() const {return num_needed_data();}

	///	number of other Data this data depends on
		virtual size_t num_needed_data() const {return m_vpIIPData.size();}

	///	return needed data
		virtual IIPData* needed_data(size_t i)
		{
			UG_ASSERT(i < m_vpIIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIIPData[i] != NULL, "Data input not valid");
			return m_vpIIPData[i];
		}

	///	returns if data is ok
		virtual bool make_ready()
		{
		//	check, that all inputs are set
			for(size_t i = 0; i < num_input(); ++i)
				if(m_vpIIPData[i] == NULL)
				{
					UG_LOG("ERROR in 'DataLinker::make_ready': Input number "<<
							i << " missing.\n");
					return false;
				}

		//	if data is dependent, forward function group
			if(!update_function_group())
			{
				UG_LOG("ERROR in 'DataLinker::set_input': Cannot build function"
						" group of linker.\n");
				return false;
			}

		//	everything ok
			return true;
		}


	protected:
	///	returns number of functions the input depends on
		size_t input_num_fct(size_t i) const
		{
			UG_ASSERT(i < m_vpIDependData.size(), "Input invalid");
			if(m_vpIDependData[i] == NULL) return 0;
			return m_vpIDependData[i]->num_fct();
		}

	///	returns the number in the common FctGrp for a fct of an input
		size_t input_common_fct(size_t i, size_t fct) const
		{
			UG_ASSERT(i < m_vMap.size(), "Input Map invalid");
			UG_ASSERT(fct < m_vMap[i].num_fct(), "Input Map invalid for fct");
			return m_vMap[i][fct];
		}

	///	returns the series id set for the i'th input
		size_t series_id(size_t i, size_t s) const
		{
			UG_ASSERT(i < m_vvSeriesID.size(), "invalid index");
			UG_ASSERT(s < m_vvSeriesID[i].size(), "invalid index");
			return m_vvSeriesID[i][s];
		}

	///	updates the function group
		bool update_function_group()
		{
		//	collect all function groups
			std::vector<const FunctionGroup*> vFctGrp(num_input(), NULL);
			for(size_t i = 0; i < m_vpIDependData.size(); ++i)
				if(m_vpIDependData[i] != NULL)
					vFctGrp[i] = &(m_vpIDependData[i]->get_function_group());

		//	create union of all function groups
			if(!CreateUnionOfFunctionGroups(m_commonFctGroup, vFctGrp, true))
			{
				UG_LOG("ERROR in 'DataLinker::update_function_group': Cannot create"
						"common function group.\n");
				return false;
			}

		//	create FunctionIndexMapping for each Disc
			m_vMap.resize(m_vpIDependData.size());
			for(size_t i = 0; i < m_vpIDependData.size(); ++i)
			{
				if(m_vpIDependData[i] != NULL)
				if(!CreateFunctionIndexMapping(m_vMap[i],
											   (m_vpIDependData[i]->get_function_group()),
											   m_commonFctGroup))
				{
					UG_LOG("ERROR in 'DataLinker::update_function_group':"
							"Cannot create Function Index Mapping.\n");
					return false;
				}
			}

		//	set common function group as the function group the data depends on
			this->set_function_group(m_commonFctGroup);

		//	we're done
			return true;
		}

	///	get series id's from input data and resize data
		virtual void adjust_global_ips_and_data(const std::vector<size_t>& vNumIP)
		{
		//	 we need a series id for all inputs
			m_vvSeriesID.resize(m_vpIIPData.size());

		//	loop inputs
			for(size_t i = 0; i < m_vpIIPData.size(); ++i)
			{
			//	resize series ids
				m_vvSeriesID[i].resize(vNumIP.size());

			//	skip unset data
				UG_ASSERT(m_vpIIPData[i] != NULL, "No Input set, but requested.");

			//	request local ips for all series at input data
				for(size_t s = 0; s < vNumIP.size(); ++s)
				{
					switch(this->dim_local_ips())
					{
						case 1:
							m_vvSeriesID[i][s] =
									m_vpIIPData[i]->template register_local_ip_series<1>
											(this->template local_ips<1>(s), num_ip(s));
							break;
						case 2:
							m_vvSeriesID[i][s] =
									m_vpIIPData[i]->template register_local_ip_series<2>
											(this->template local_ips<2>(s), num_ip(s));
							break;
						case 3:
							m_vvSeriesID[i][s] =
									m_vpIIPData[i]->template register_local_ip_series<3>
											(this->template local_ips<3>(s), num_ip(s));
							break;
						default: throw(UGFatalError("Dimension not supported."));
					}
				}
			}

		//	resize data fields
			DependentIPData<TData, dim>::adjust_global_ips_and_data(vNumIP);
		}

	/// implement callback, that is called when global ips changed
		virtual void global_ips_changed(size_t s, const MathVector<dim>* vPos, size_t numIP)
		{
		//	loop inputs
			for(size_t i = 0; i < m_vpIIPData.size(); ++i)
			{
			//	skip unset data
				UG_ASSERT(m_vpIIPData[i] != NULL, "No Input set, but requested.");

			//	request local ips for all series at input data
				for(size_t s = 0; s < m_vvSeriesID[i].size(); ++s)
				{
				//	adjust global ids of imported data
					m_vpIIPData[i]->set_global_ips(m_vvSeriesID[i][s], vPos, numIP);
				}
			}
		}


	protected:
	///	data input
		std::vector<IIPDimData<dim>*> m_vpIIPData;

	///	data input casted to IDependend data
		std::vector<IDependentIPData*> m_vpIDependData;

	///	common functions the data depends on
		FunctionGroup m_commonFctGroup;

	///	Function mapping for each input relative to common FunctionGroup
		std::vector<FunctionIndexMapping> m_vMap;

	///	series id the linker uses to get data from input
		std::vector<std::vector<size_t> > m_vvSeriesID;
};

/// combines several IPDatas of the one data type to a new IPData of a second type
/**
 * This class provides data at integration points and implements the
 * DependentIPData interface. It gets passed several data of one C++-type
 * and combines them to a new value of a second type. In addition also the
 * derivative w.r.t. the primal unknowns are passed.
 *
 * \tparam 	TData		output Data type
 * \tparam 	dim			World dimension
 * \tparam	TDataIn		Input data type
 */
template <typename TData, int dim, typename TDataIn>
class DataLinkerEqualData
	: public DataLinker<TData, dim>
{
	public:
	///	Base class type
		typedef DataLinker<TData, dim> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::time;

	//	explicitly forward methods of IPData
		using base_type::ip;
		using base_type::value;

	//	explicitly forward methods of IDataLinker
		using base_type::num_input;
		using base_type::series_id;

	public:
	///	constructor
		DataLinkerEqualData()
		{
			m_vpIPData.clear();
			m_vpDependData.clear();
		}

	///	compute method
		virtual void compute(bool compDeriv) = 0;

	///	set number of needed inputs
		void set_num_input(size_t num)
		{
		//	resize arrays
			m_vpIPData.resize(num, NULL);
			m_vpDependData.resize(num, NULL);

		//	forward size to base class
			base_type::set_num_input(num);
		}

	///	set input i
		bool set_input(size_t i, IPData<TDataIn, dim>& data)
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");

		//	check input number
			if(i >= num_input())
			{
				UG_LOG("ERROR in 'DataLinker::set_input': Only " << num_input()
				       << " inputs can be set. Use 'set_num_input' to increase"
				       " the number of needed inputs.\n");
				return false;
			}

		//	remember ipdata
			m_vpIPData[i] = &data;

		//	cast to dependent data
			m_vpDependData[i] = dynamic_cast<DependentIPData<TDataIn, dim>*>(&data);

		//	forward to base class
			base_type::set_input(i, &data);

		//	we're done
			return true;
		}

	protected:
	///	data at ip of input
		const TDataIn& input_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i] != NULL, "Input invalid");
			return m_vpIPData[i]->value(series_id(i,s), ip);
		}

	///	data at ip of input
		TDataIn& input_value(size_t i, size_t s, size_t ip)
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i] != NULL, "Input invalid");
			return m_vpIPData[i]->value(series_id(i,s), ip);
		}

	///	derivative of data at input at ip
		const TDataIn& input_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");
			UG_ASSERT(m_vpDependData[i] != NULL, "Input invalid");
			return m_vpDependData[i]->deriv(series_id(i,s), ip, fct, dof);
		}

	///	derivative of data at input at ip
		TDataIn& input_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof)
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");
			UG_ASSERT(m_vpDependData[i] != NULL, "Input invalid");
			return m_vpDependData[i]->deriv(series_id(i,s), ip, fct, dof);
		}

	protected:
	///	data input
		std::vector<IPData<TDataIn, dim>*> m_vpIPData;

	///	data input casted to dependend data
		std::vector<DependentIPData<TDataIn, dim>*> m_vpDependData;
};

/////////////////////////////////////////////////
// Linker Traits to allow generic programming
/////////////////////////////////////////////////

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


/////////////////////////////////////////////////
// Scaled adding of Data
/////////////////////////////////////////////////

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
	: public DataLinker<TData, dim>
{
	public:
	//	type of base class
		typedef DataLinker<TData, dim> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::time;

	//	explicitly forward methods of IPData
		using base_type::ip;
		using base_type::value;

	//	explicitly forward methods of IDependentIPData
		using base_type::num_fct;

	//	explicitly forward methods of DependentIPData
		using base_type::num_sh;
		using base_type::deriv;

	//	explicitly forward methods of DataLinker
		using base_type::series_id;

	public:
	///	constructor
		ScaleAddLinker() {}

	///	adds an input to the list of summands scaled by a user data factor
		bool add(IPData<TDataScale, dim>& scale, IPData<TData, dim>& data)
		{
		//	current number of inputs
			const size_t numInput = base_type::num_input() / 2;

		//	resize scaling
			resize_scaling(numInput+1);

		//	remember ipdata
			m_vpIPData[numInput] = &data;
			UG_ASSERT(m_vpIPData[numInput] != NULL, "Null Pointer as Input set.");
			m_vpDependData[numInput] = dynamic_cast<DependentIPData<TData, dim>*>(&data);

		//	remember ipdata
			m_vpScaleData[numInput] = &scale;
			UG_ASSERT(m_vpScaleData[numInput] != NULL, "Null Pointer as Scale set.");
			m_vpScaleDependData[numInput]
			              = dynamic_cast<DependentIPData<TDataScale, dim>*>(&scale);

		//	increase number of inputs by one
			base_type::set_num_input(2*numInput+2);

		//	add this input
			base_type::set_input(2*numInput, &data);
			base_type::set_input(2*numInput+1, &scale);

		//	done
			return true;
		}

	///	computes the value
		virtual void compute(bool compDeriv)
		{
		//	check that size of Scalings and inputs is equal
			UG_ASSERT(m_vpIPData.size() == m_vpScaleData.size(), "Wrong num Scales.");

		//	compute value
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
				//	reset value
					value(s,ip) = 0.0;

				//	add contribution of each summand
					for(size_t c = 0; c < m_vpIPData.size(); ++c)
					{
						linker_traits<TData, TDataScale>::
						mult_add(value(s, ip),
						         input_value(c, s, ip),
						         scale_value(c, s, ip));
					}
				}

		//	check if derivative is required
			if(!compDeriv || this->zero_derivative()) return;

		//	check sizes
			UG_ASSERT(m_vpDependData.size() == m_vpScaleDependData.size(),
			          	  	  	  	  	  	  	  	  	  "Wrong num Scales.");

		//	clear all derivative values
			this->clear_derivative_values();

		//	loop all inputs
			for(size_t c = 0; c < m_vpIPData.size(); ++c)
			{
			//	check if input has derivative
				if(!m_vpIPData[c]->zero_derivative())
				{
					for(size_t s = 0; s < num_series(); ++s)
						for(size_t ip = 0; ip < num_ip(s); ++ip)
						{
						//	loop functions
							for(size_t fct = 0; fct < input_num_fct(c); ++fct)
							{
							//	get common fct id for this function
								const size_t commonFct = input_common_fct(c, fct);

							//	loop dofs
								for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
								{
									linker_traits<TData, TDataScale>::
									mult_add(deriv(s, ip, commonFct, sh),
									         input_deriv(c, s, ip, fct, sh),
									         scale_value(c, s, ip));
								}
							}
						}
				}

			//	check if scaling has derivative
				if(!m_vpScaleData[c]->zero_derivative())
				{
					for(size_t s = 0; s < num_series(); ++s)
						for(size_t ip = 0; ip < num_ip(s); ++ip)
						{
						//	loop functions
							for(size_t fct = 0; fct < scale_num_fct(c); ++fct)
							{
							//	get common fct id for this function
								const size_t commonFct = scale_common_fct(c, fct);

							//	loop dofs
								for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
								{
									linker_traits<TData, TDataScale>::
									mult_add(deriv(s, ip, commonFct, sh),
											 input_value(c, s, ip),
											 scale_deriv(c, s, ip, fct, sh));
								}
							}
						}
				}
			}
		}

	protected:
	///	resizes the scaling arrays
		void resize_scaling(size_t num)
		{
			m_vpIPData.resize(num, NULL);
			m_vpDependData.resize(num, NULL);
			m_vpScaleData.resize(num, NULL);
			m_vpScaleDependData.resize(num, NULL);
		}

	///	data at ip of input
		const TData& input_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpIPData.size(), "Input not needed");
			UG_ASSERT(m_vpIPData[i] != NULL, "Input invalid");
			return m_vpIPData[i]->value(series_id(2*i,s), ip);
		}

	///	derivative of data at input at ip
		const TData& input_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(i < m_vpDependData.size(), "Input not needed");
			UG_ASSERT(m_vpDependData[i] != NULL, "Input invalid");
			return m_vpDependData[i]->deriv(series_id(2*i,s), ip, fct, dof);
		}

	///	scale at ip of input
		const TDataScale& scale_value(size_t i, size_t s, size_t ip) const
		{
			UG_ASSERT(i < m_vpScaleData.size(), "Input not needed");
			UG_ASSERT(m_vpScaleData[i] != NULL, "Input invalid");
			return m_vpScaleData[i]->value(series_id(2*i+1,s), ip);
		}

	///	derivative of data at input at ip
		const TDataScale& scale_deriv(size_t i, size_t s, size_t ip, size_t fct, size_t dof) const
		{
			UG_ASSERT(i < m_vpScaleDependData.size(), "Input not needed");
			UG_ASSERT(m_vpScaleDependData[i] != NULL, "Input invalid");
			return m_vpScaleDependData[i]->deriv(series_id(2*i+1,s), ip, fct, dof);
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
		std::vector<IPData<TDataScale, dim>*> m_vpScaleData;

	///	data input casted to dependend data
		std::vector<DependentIPData<TDataScale, dim>*> m_vpScaleDependData;

	///	data input
		std::vector<IPData<TData, dim>*> m_vpIPData;

	///	data input casted to dependend data
		std::vector<DependentIPData<TData, dim>*> m_vpDependData;
};


} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DATA_LINKER__ */
