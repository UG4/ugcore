/*
 * const_user_data.h
 *
 *  Created on: 12.07.2012
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__CONST_DATA_USER_DATA__
#define __H__LIB_DISCRETIZATION__CONST_DATA_USER_DATA__

#include "import_export.h"

namespace ug {

// Export Possibility for user data.
// It has no slots.
template <typename TDataType>
class ConstUserDataExportPossibility : public DataPossibilityItem
{
	public:
		typedef TDataType data_type;

	public:
		ConstUserDataExportPossibility(std::string name, TDataType& val) :
			DataPossibilityItem(name, 0, &typeid(TDataType)), m_val(val)
			{m_vCreatedDataExports.clear();};

	public:
		// create a data export from this possibility
		virtual DataExportItem* create_data_export();

		// Destructor
		virtual ~ConstUserDataExportPossibility();

	private:
		TDataType m_val;
};


// Export Possibility of User Data.
template <typename TDataType>
class ConstUserDataExport : public DataExport<TDataType>{
	friend class DataImport<TDataType>;

	public:
		typedef TDataType data_type;

	public:
		ConstUserDataExport(std::string name, DataPossibilityItem* possibility, TDataType& val) 	:
			DataExport<TDataType>(name, possibility), m_val(val)
			{};

	public:
		// compute
		virtual void compute(bool compute_derivatives)
		{
			for(size_t ip = 0; ip < this->num_ip(); ++ip)
			{
				this->m_vValue[ip] = m_val;
			}
		}

	protected:
		TDataType m_val;
};

////// IMPLEMENTATION /////

template<typename TDataType>
ConstUserDataExportPossibility<TDataType>::
~ConstUserDataExportPossibility()
{
	for(size_t i = 0; i < m_vCreatedDataExports.size(); ++i)
	{
		UG_ASSERT(delete_data_export(m_vCreatedDataExports[i]),
					"DataClassExportPossibility::~DataClassExportPossibility: Cannot delete Exports.");
	}
}

template<typename TDataType>
DataExportItem*
ConstUserDataExportPossibility<TDataType>::
create_data_export()
{
	ConstUserDataExport<TDataType> * exp =
		new ConstUserDataExport<TDataType>(this->name(), this, m_val);
	m_vCreatedDataExports.push_back(dynamic_cast<DataExportItem*>(exp));

	// A Class export does only depend on one (its) system
	if(!exp->set_num_sys(0)) return NULL;

	return dynamic_cast<DataExportItem*>(exp);
}



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__CONST_DATA_USER_DATA__ */
