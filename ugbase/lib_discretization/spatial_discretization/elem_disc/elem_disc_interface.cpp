/*
 * elem_disc_interface.cpp
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "elem_disc_interface.h"

namespace ug{

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
	if(!prepare_element_loop_fct_registered(id))
	{
		UG_LOG("prepare_element_loop() function not registered for id " << id <<".\n");
		return false;
	}
	if(!prepare_element_fct_registered(id))
	{
		UG_LOG("prepare_element(...) function not registered for id " << id <<".\n");
		return false;
	}
	if(!finish_element_loop_fct_registered(id))
	{
		UG_LOG("finish_element_loop() function not registered for id " << id <<".\n");
		return false;
	}

	return true;
};

bool IElemDisc::prepare_element_loop_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPrepareElemLoopFunc.size())
	{
		if(m_vPrepareElemLoopFunc[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::prepare_element_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vPrepareElemFunc.size())
	{
		if(m_vPrepareElemFunc[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::finish_element_loop_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vFinishElemLoopFunc.size())
	{
		if(m_vFinishElemLoopFunc[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_JA_vol_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssJAFunc.size())
	{
		if(m_vAssJAFunc[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_JM_vol_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssJMFunc.size())
	{
		if(m_vAssJMFunc[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_dA_vol_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssAFunc.size())
	{
		if(m_vAssAFunc[id] != 0)
			return true;
	}
	return false;
}

bool IElemDisc::ass_dM_vol_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssMFunc.size())
	{
		if(m_vAssMFunc[id] != 0)
			return true;
	}
	return false;
}


bool IElemDisc::ass_rhs_vol_fct_registered(int id)
{
	if(id >= 0 && (size_t)id < m_vAssFFunc.size())
	{
		if(m_vAssFFunc[id] != 0)
			return true;
	}
	return false;
}

} // end namespace ug
