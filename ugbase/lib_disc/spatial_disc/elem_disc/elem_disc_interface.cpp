/*
 * elem_disc_interface.cpp
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "elem_disc_interface.h"

namespace ug{

IElemDisc::IElemDisc(int numFct, const char* functions, const char* subsets)
	: 	m_numFct(numFct), m_bTimeDependent(false), m_time(0.0),
	  	m_pLocalVectorTimeSeries(NULL), m_id(ROID_UNKNOWN)
{
	set_functions(functions);
	set_subsets(subsets);
}

void IElemDisc::set_functions(std::string fctString)
{
//	tokenize string
	TokenizeString(fctString, m_vFct, ',');

//	remove white space
	for(size_t i = 0; i < m_vFct.size(); ++i)
		RemoveWhitespaceFromString(m_vFct[i]);

//	check, that number of symbols is correct
	if(m_vFct.size() != m_numFct)
	{
		UG_THROW("IElemDisc: Wrong number of symbolic "
				 "names of functions passed: Required: "<< m_numFct <<
				 " Passed: " << m_vFct.size() << " ('"<<fctString<<"'). Please "
				 " pass correct number of symbolic names separated by ','.\n");
	}
}

void IElemDisc::set_subsets(std::string ssString)
{
//	tokenize string
	TokenizeString(ssString, m_vSubset, ',');

//	remove white space
	for(size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

//	check that at least one subset given
	if(m_vSubset.empty())
		UG_THROW("IElemDisc: At least one Subset must be specified for an element disc.")
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

bool IElemDisc::set_roid(ReferenceObjectID id)
{
	m_id = id;
	return true;

//	\todo: error check
	{
		m_id = ROID_UNKNOWN;
		UG_LOG("No or not all functions registered "
				"for object with reference object id " << id << ".\n");
	}
	return false;
};

} // end namespace ug
