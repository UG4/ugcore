/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Sam Gimbel
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

// C++ lib.
#include <iostream>

// UG4 lib.
#include "common/error.h"
#include "lib_disc/spatial_disc/user_data/linker/linker_traits.h"

#ifdef UG_USE_PYBIND11
#include <Python.h>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace ug
{
////////////////////////////////////////////////////////////////////////////////
// Executing Python Code from Ug4
// Code has been passed as String to UG4
// More Info here: https://docs.python.org/3/extending/embedding.html#
// and also here: https://docs.python.org/3/c-api/ 
// also make sure to have python3-dev installed to use the python header:
// sudo apt install python3-dev
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// PythonUser
// Please Note: This is just a quick fix
// the PythonUserFunction class should be extended like the lua equivalent class
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
// PythonUserFunction
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TDataIn>
PythonUserFunction<TData,dim,TDataIn>::
PythonUserFunction(TFunctionHandle pyCallback, size_t numArgs)
	: m_numArgs(numArgs), m_cbValueRef(pyCallback), m_cbDerivRef(numArgs)//, m_bPosTimeNeed(false)
{
	debug();
	this->set_num_input(numArgs);
}




template <typename TData, int dim, typename TDataIn>
PythonUserFunction<TData,dim,TDataIn>::~PythonUserFunction()
{
	debug();
}


template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::set_deriv(size_t i, TFunctionHandle pyCallback)
{
//	check number of arg
	debug();
	UG_COND_THROW(i >= m_numArgs,
			"PythonUserFunction::set_lua_deriv_callback: Trying " <<
			"to set a derivative for argument " << i <<", that " <<
			"does not exist. Number of arguments is "<<m_numArgs);

//	store name (string) of callback
	m_cbDerivRef[i] = pyCallback;

//	free old reference
//	free_deriv_callback_ref(arg);


}


template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::operator() (TData& out, const std::vector<TDataIn> &input) const
{
	debug();

	UG_ASSERT(input.size() == m_numArgs, "Number of arguments mismatched.");

	//	call python function
	auto pyArgs = py::tuple(py::cast(input));
	py::object result_py = m_cbValueRef(*pyArgs);
	out = result_py.cast<TData>();

}

template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::operator() (TData& out, int numArgs, ...) const
{
	debug();
	UG_ASSERT(numArgs == (int)m_numArgs, "Number of arguments mismatched.");

	// read vector of all arguments.
	std::vector<TDataIn> input(numArgs);

	va_list myArgs;
	va_start(myArgs, numArgs);
	for(int i = 0; i < numArgs; ++i)
	{ input[i] = va_arg(myArgs, TDataIn); }
	va_end(myArgs);					//	end read in of parameters

	// call python function via operator
	(*this)(out, input);
}


template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::eval_value(TData& out, const std::vector<TDataIn>& dataIn,
													const MathVector<dim>& x, number time, int si) const
{
	debug();
	(*this)(out, dataIn);
}


template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::eval_deriv(TData& out, const std::vector<TDataIn>& dataIn,
		 	 	 	 	 	 	 	 	 	 	 	 const MathVector<dim>& x, number time, int si, size_t arg) const
{
	debug();
	UG_ASSERT(dataIn.size() == m_numArgs, "Number of arguments mismatched.");
	UG_ASSERT(arg < m_numArgs, "Argument does not exist.");

	//	push the callback function on the stack

	//	read all arguments and push them to the lua stack
	//	if needed, read additional coordinate, time and subset index arguments and push them to the lua stack
	auto pyArgs = py::tuple(py::cast(dataIn));

	// Call function.
	py::object result_py = m_cbDerivRef[arg](*pyArgs);
	out = result_py.cast<TData>();

}


template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::
evaluate (TData& value,
          const MathVector<dim>& globIP,
          number time, int si) const
{
	debug();
	// PROFILE_CALLBACK();
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t c = 0; c < vDataIn.size(); ++c)
		(*m_vpUserData[c])(vDataIn[c], globIP, time, si);

//	evaluate data at ip
	eval_value(value, vDataIn, globIP, time, si);

	//UG_COND_THROW(IsFiniteAndNotTooBig(value)==false, value);
}

template <typename TData, int dim, typename TDataIn>
template <int refDim>
void PythonUserFunction<TData,dim,TDataIn>::
evaluate(TData vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<refDim> vLocIP[],
         const size_t nip,
         LocalVector* u,
         const MathMatrix<refDim, dim>* vJT) const
{
	debug();
	// PROFILE_CALLBACK();
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t ip = 0; ip < nip; ++ip)
	{
		for(size_t c = 0; c < vDataIn.size(); ++c)
			(*m_vpUserData[c])(vDataIn[c], vGlobIP[ip], time, si, elem, vCornerCoords, vLocIP[ip], u);

	//	evaluate data at ip
		eval_value(vValue[ip], vDataIn, vGlobIP[ip], time, si);
	//	UG_COND_THROW(IsFiniteAndNotTooBig(vValue[ip])==false, vValue[ip]);
	}
}

template <typename TData, int dim, typename TDataIn>
template <int refDim>
void PythonUserFunction<TData,dim,TDataIn>::
eval_and_deriv(TData vValue[],
				const MathVector<dim> vGlobIP[],
				number time, int si,
				GridObject* elem,
				const MathVector<dim> vCornerCoords[],
				const MathVector<refDim> vLocIP[],
				const size_t nip,
				LocalVector* u,
				bool bDeriv,
				int s,
				std::vector<std::vector<TData> > vvvDeriv[],
				const MathMatrix<refDim, dim>* vJT)
{
	debug();
	// PROFILE_CALLBACK();
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	gather all input data for this ip
		for(size_t c = 0; c < vDataIn.size(); ++c)
			vDataIn[c] = m_vpUserData[c]->value(this->series_id(c,s), ip);

	//	evaluate data at ip
		eval_value(vValue[ip], vDataIn, vGlobIP[ip], time, si);
	}

//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return;

//	clear all derivative values
	this->set_zero(vvvDeriv, nip);

//	loop all inputs
	for(size_t c = 0; c < vDataIn.size(); ++c)
	{
	//	check if we have the derivative w.r.t. this input, and the input has derivative
		if(/*m_cbDerivRef[c] == LUA_NOREF || */ m_vpUserData[c]->zero_derivative()) continue;

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	gather all input data for this ip
			for(size_t i = 0; i < vDataIn.size(); ++i)
				vDataIn[i] = m_vpUserData[i]->value(this->series_id(c,s), ip); //< series_id(c,s) or series_id(i,s)

		//	data of derivative w.r.t. one component at ip-values
			TData derivVal;

		//	evaluate data at ip
			eval_deriv(derivVal, vDataIn, vGlobIP[ip], time, si, c);

		//	loop functions
			for(size_t fct = 0; fct < this->input_num_fct(c); ++fct)
			{
			//	get common fct id for this function
				const size_t commonFct = this->input_common_fct(c, fct);

			//	loop dofs
				for(size_t dof = 0; dof < this->num_sh(fct); ++dof)
				{
					linker_traits<TData, TDataIn>::
					mult_add(vvvDeriv[ip][commonFct][dof],
							 derivVal,
							 m_vpDepentData[c]->deriv(this->series_id(c,s), ip, fct, dof));
				//	UG_COND_THROW(IsFiniteAndNotTooBig(vvvDeriv[ip][commonFct][dof])==false, vvvDeriv[ip][commonFct][dof]);
				}
			}
		}
	}
}

/**
 * TODO: Note this is a public (non-virtual) function whose argument
 * should be consistent with the number of the arguments. Should not it also
 * resize the array for the references to the derivatives?
 */
template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::set_num_input(size_t num)
{
	debug();
//	resize arrays
	m_vpUserData.resize(num);
	m_vpDepentData.resize(num);

//	forward size to base class
	base_type::set_num_input(num);
}

template <typename TData, int dim, typename TDataIn>
void PythonUserFunction<TData,dim,TDataIn>::
set_input(size_t i, SmartPtr<CplUserData<TDataIn, dim> > data)
{
	debug();
	UG_ASSERT(i < m_vpUserData.size(), "Input not needed");
	UG_ASSERT(i < m_vpDepentData.size(), "Input not needed");

//	check input number

	UG_COND_THROW(i >= this->num_input(),
			"PythonUserFunction::set_input: Only " << this->num_input() <<
			" inputs can be set. Use 'set_num_input' to increase" <<
			" the number of needed inputs.");

//	remember userdata & cast to dependent data
	m_vpUserData[i] = data;
	m_vpDepentData[i] = data.template cast_dynamic<DependentUserData<TDataIn, dim> >();

//	forward to base class
	base_type::set_input(i, data, data);
}





} // end of namespace ug

#endif // UG_USE_PYBIND11
