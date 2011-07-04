/*
 * elem_disc_interface.cpp
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "elem_disc_interface.h"

namespace ug{

IElemDisc::IElemDisc()
	: 	m_bTimeDependent(false), m_time(0.0),
	  	m_pLocalVectorTimeSeries(NULL), m_id(-1)
{}

bool IElemDisc::set_functions(const char* functions)
{
//	get strings
	std::string fctString = std::string(functions);

//	tokenize string
	TokenizeString(fctString, m_vFct, ',');

//	remove white space
	for(size_t i = 0; i < m_vFct.size(); ++i)
		RemoveWhitespaceFromString(m_vFct[i]);

//	check, that number of symbols is correct
	if(m_vFct.size() != this->num_fct())
	{
		UG_LOG("ERROR in 'IElemDisc::set_functions': Wrong number of symbolic "
				"names of functions passed: Required: "<< this->num_fct() <<
				" Passed: " << m_vFct.size() << " ('"<<functions<<"'). Please "
				" pass correct number of symbolic names separated by ','.\n");
		return false;
	}

//	done
	return true;
}

bool IElemDisc::set_subsets(const char* subsets)
{
//	get strings
	std::string ssString = std::string(subsets);

//	tokenize string
	TokenizeString(ssString, m_vSubset, ',');

//	remove white space
	for(size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

//	done
	return true;
}

void IElemDisc::register_import(IDataImport& Imp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		if(m_vIImport[i] == &Imp)
			throw(UGFatalError("Trying to register import twice."));

//	add it
	m_vIImport.push_back(&Imp);
}

void IElemDisc::register_export(IDataExport& Exp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIExport.size(); ++i)
		if(m_vIExport[i] == &Exp)
			throw(UGFatalError("Trying to register export twice."));

//	add it
	m_vIExport.push_back(&Exp);
}

IDataImport& IElemDisc::get_import(size_t i)
{
	UG_ASSERT(i < num_imports(), "Invalid index");
	return *m_vIImport[i];
}

IDataExport& IElemDisc::get_export(size_t i)
{
	UG_ASSERT(i < num_exports(), "Invalid index");
	return *m_vIExport[i];
}

bool IElemDisc::set_geometric_object_type(int id)
{
	if(function_registered(id))
	{
		m_id = id;
		return true;
	}
	else
	{
		m_id = -1;
		UG_LOG("No or not all functions registered "
				"for object with reference object id " << id << ".\n");
	}
	return false;
};

bool IElemDisc::function_registered(int id)
{
	// loop functions must exist in any case
	if(!prepare_elem_loop_fct_registered(id))
	{
		UG_LOG("prepare_element_loop() function not registered for id " << id <<".\n");
		return false;
	}
	if(!prepare_elem_fct_registered(id))
	{
		UG_LOG("prepare_element(...) function not registered for id " << id <<".\n");
		return false;
	}
	if(!finish_elem_loop_fct_registered(id))
	{
		UG_LOG("finish_element_loop() function not registered for id " << id <<".\n");
		return false;
	}

	return true;
};

bool IElemDisc::prepare_elem_loop_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPrepareElemLoopFct.size())
	{
		if(m_vPrepareElemLoopFct[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::prepare_elem_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPrepareElemFct.size())
	{
		if(m_vPrepareElemFct[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::finish_elem_loop_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vFinishElemLoopFct.size())
	{
		if(m_vFinishElemLoopFct[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_JA_elem_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vElemJAFct.size())
	{
		if(m_vElemJAFct[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_JM_elem_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vElemJMFct.size())
	{
		if(m_vElemJMFct[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_dA_elem_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vElemdAFct.size())
	{
		if(m_vElemdAFct[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_dM_elem_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vElemdMFct.size())
	{
		if(m_vElemdMFct[id] != 0)
			return true;
	}
	return false;
}


bool IElemDisc::ass_rhs_elem_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vElemRHSFct.size())
	{
		if(m_vElemRHSFct[id] != 0)
			return true;
	}
	return false;
}

} // end namespace ug
