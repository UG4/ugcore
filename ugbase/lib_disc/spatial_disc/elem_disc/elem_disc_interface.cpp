/*
 * elem_disc_interface.cpp
 *
 *  Created on: 03.07.2011
 *      Author: andreasvogel
 */

#include "elem_disc_interface.h"

namespace ug{

IElemDisc::IElemDisc(const char* functions, const char* subsets)
	: 	m_bTimeDependent(false), m_timePoint(0),
	  	m_pLocalVectorTimeSeries(NULL),
	  	m_bFastAssembleEnabled(false), m_id(ROID_UNKNOWN)
{
	m_vFct.clear();
	m_vSubset.clear();
	if(functions) set_functions(functions);
	if(subsets) set_subsets(subsets);
}

IElemDisc::IElemDisc(const std::vector<std::string>& vFct,
                     const std::vector<std::string>& vSubset)
	: 	m_bTimeDependent(false), m_timePoint(0),
		m_pLocalVectorTimeSeries(NULL),
		m_bFastAssembleEnabled(false), m_id(ROID_UNKNOWN)
{
	m_vFct = vFct;
	m_vSubset = vSubset;
}


void IElemDisc::set_functions(std::string fctString)
{
//	tokenize string
	TokenizeString(fctString, m_vFct, ',');

//	remove white space
	for(size_t i = 0; i < m_vFct.size(); ++i)
		RemoveWhitespaceFromString(m_vFct[i]);

//	if no function passed, clear functions
	if(m_vFct.size() == 1 && m_vFct[0].empty()) m_vFct.clear();

//	if functions passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vFct.size(); ++i)
	{
		if(m_vFct.empty())
			UG_THROW("Error while setting functions in an ElemDisc: passed "
							"function string '"<<fctString<<"' lacks a "
							"function specification at position "<<i<<"(of "
							<<m_vFct.size()-1<<")");
	}
}

void IElemDisc::set_subsets(std::string ssString)
{
//	tokenize string
	TokenizeString(ssString, m_vSubset, ',');

//	remove white space
	for(size_t i = 0; i < m_vSubset.size(); ++i)
		RemoveWhitespaceFromString(m_vSubset[i]);

//	if no subset passed, clear subsets
	if(m_vFct.size() == 1 && m_vFct[0].empty()) m_vFct.clear();

//	if subsets passed with separator, but not all tokens filled, throw error
	for(size_t i = 0; i < m_vFct.size(); ++i)
	{
		if(m_vFct.empty())
			UG_THROW("Error while setting subsets in an ElemDisc: passed "
							"subset string '"<<ssString<<"' lacks a "
							"subset specification at position "<<i<<"(of "
							<<m_vFct.size()-1<<")");
	}
}

void IElemDisc::register_import(IDataImport& Imp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIImport.size(); ++i)
		if(m_vIImport[i] == &Imp)
			UG_THROW("Trying to register import twice.");

//	add it
	m_vIImport.push_back(&Imp);
}

void IElemDisc::register_export(SmartPtr<IUserData> Exp)
{
//	check that not already registered
	for(size_t i = 0; i < m_vIExport.size(); ++i)
		if(m_vIExport[i] == Exp)
			UG_THROW("Trying to register export twice.");

//	add it
	m_vIExport.push_back(Exp);
}

IDataImport& IElemDisc::get_import(size_t i)
{
	UG_ASSERT(i < num_imports(), "Invalid index");
	return *m_vIImport[i];
}

SmartPtr<IUserData> IElemDisc::get_export(size_t i)
{
	UG_ASSERT(i < num_exports(), "Invalid index");
	return m_vIExport[i];
}

void IElemDisc::set_roid(ReferenceObjectID id)
{
	m_id = id;

	if(id == ROID_UNKNOWN)
	{
		m_id = ROID_UNKNOWN;
		UG_THROW("Cannot assemble for RefID: "<<id<<".\n");
	}
};

void IElemDisc::set_time_dependent(bool bTimeDependent,
                                   const LocalVectorTimeSeries* locTimeSeries)
{
	m_bTimeDependent = bTimeDependent;
	m_pLocalVectorTimeSeries = locTimeSeries;
}

} // end namespace ug
