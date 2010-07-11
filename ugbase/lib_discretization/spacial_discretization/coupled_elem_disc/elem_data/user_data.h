/*
 * element_data_user_data.h
 *
 *  Created on: 12.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__ELEMENT_DATA_USER_DATA__
#define __H__LIB_DISCRETIZATION__ELEMENT_DATA_USER_DATA__

#include "element_data_import_export.h"

namespace ug {

// Export Possibility for user data.
// It has no slots.
template <typename TDataType, typename TPositionType>
class UserDataExportPossibility : public DataPossibilityItem
{
	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;

	public:
		typedef void (*UserFunction)(const TPositionType& pos, TDataType& data);

	public:
		UserDataExportPossibility(std::string name, UserFunction func) :
			DataPossibilityItem(name, 0, &typeid(TDataType), &typeid(TPositionType))
			{m_vCreatedDataExports.clear();};

	public:
		// create a data export from this possibility
		virtual DataExportItem* create_data_export();

		// Destructor
		virtual ~UserDataExportPossibility();

	protected:
		UserFunction m_userFunc;
};


// Export Possibility of User Data.
// It has no slots.
// TODO: This assumes, that TPositionType for evaluation is the same as the one for Import.
//       But we may have different local and global types
template <typename TDataType, typename TPositionType>
class UserDataExport : public DataExport<TDataType, TPositionType>{
	friend class DataImport<TDataType,TPositionType>;

	public:
		typedef TDataType data_type;
		typedef TPositionType position_type;

	public:
		typedef void (*UserFunction)(const TPositionType& pos, TDataType& data);

	public:
		UserDataExport(std::string name, DataPossibilityItem* possibility, UserFunction func) 	:
			DataExport<TDataType, TPositionType>(name, possibility),
			m_userFunc(func)
			{};

	public:
		// return, if equal to export v
		virtual bool equal(const DataExportItem& v) const;

		// compute
		virtual void compute(bool compute_derivatives)
		{
			for(size_t ip = 0; ip < this->num_ip(); ++ip)
			{
				m_userFunc(this->m_vPosition[ip], this->m_vValues[ip]);
			}
		}

	protected:
		UserFunction m_userFunc;
};

////// IMPLEMENTATION /////

template<typename TDataType, typename TPositionType>
UserDataExportPossibility<TDataType, TPositionType>::
~UserDataExportPossibility()
{
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		UG_ASSERT(delete_data_export(m_vCreatedDataExports[i]),
					"DataClassExportPossibility::~DataClassExportPossibility: Cannot delete Exports.");
	}
}

template<typename TDataType, typename TPositionType>
DataExportItem*
UserDataExportPossibility<TDataType, TPositionType>::
create_data_export()
{
	UserDataExport<TDataType, TPositionType> * exp =
		new UserDataExport<TDataType, TPositionType>(this->name(), m_userFunc);
	m_vCreatedDataExports.push_back(dynamic_cast<DataExportItem*>(exp));

	// A Class export does only depend on one (its) system
	if(!exp->set_num_sys(0)) return NULL;

	return dynamic_cast<DataExportItem*>(exp);
}

template<typename TDataType, typename TPositionType>
bool
UserDataExport<TDataType, TPositionType>::
equal(const DataExportItem& v) const
{
	UserDataExport<TDataType, TPositionType>* exp =
		dynamic_cast<UserDataExport<TDataType, TPositionType>*>(v);

	if(exp == NULL) return false;
	if(this->m_pPossibilityItem != exp->get_possibility_item()) return false;

	// check positions
	if(this->num_ip() != exp->num_ip()) return false;
	for(size_t ip = 0; ip < this->num_ip(); ++ip)
	{
		if(this->m_vPosition[ip] != exp->m_vPosition[ip]) return false;
	}
	return true;
}



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__ELEMENT_DATA_USER_DATA__ */
