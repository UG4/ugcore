/*
 * elementdata_impl.h
 *
 *  Created on: 09.11.2009
 *      Author: andreasvogel
 */
#ifndef __H__LIB_DISCRETIZATION__ELEMENTDATA_IMPL__
#define __H__LIB_DISCRETIZATION__ELEMENTDATA_IMPL__

#include <iostream>
#include <cassert>

namespace ug{

template<typename TDataType, typename TPositionType>
DataExport<TDataType, TPositionType>::DataExport(std::string name, EvalFunction func) : DataExportItem(name)
{
	_DataType = &typeid(TDataType);
	_PositionType = &typeid(TPositionType);
	_evalFunction = func;
};
template<typename TDataType, typename TPositionType>
DataExport<TDataType, TPositionType>::~DataExport()
{
	for(int i = 0; i < _ImportList.size(); ++i)
	{
		delete[] _ImportList[i].valueList;
		delete[] _ImportList[i].positionList;
		(_ImportList[i].import)->_Export = 0;
	}

	_ImportList.clear();
};

template<typename TDataType, typename TPositionType>
void DataExport<TDataType, TPositionType>::compute()
{
	//compute values
	for(int i = 0; i < _ImportList.size(); ++i)
	{
		(*_evalFunction)(_ImportList[i].valueList, _ImportList[i].positionList, _ImportList[i].n);
	}

	//write values to import items
	writeDataToImport();
}

template<typename TDataType, typename TPositionType>
void DataExport<TDataType, TPositionType>::writeDataToImport()
{
	//write values to import items
	for(int i = 0; i < _ImportList.size(); ++i)
	{
		for(int j = 0; j < _ImportList[i].n; ++j)
		{
			((_ImportList[i].import)->_valueList)[j] = (_ImportList[i].valueList)[j];
		}
	}
}

template<typename TDataType, typename TPositionType>
bool DataExport<TDataType, TPositionType>::registerDataImport(DataImport<TDataType, TPositionType>* Import)
{
	int i;
	for(i = 0; i < _ImportList.size(); ++i)
	{
		if(_ImportList[i].import == Import) break;
	}
	if(i < _ImportList.size())
	{
		std::cout << "registerDataImport: DataImport already registered to DataExport. Invalid operation." << std::endl;
		return false;
	}

	DataImportInfo info;

	info.import = Import;
	info.n = Import->_n;
	info.valueList = new TDataType[Import->_n];
	info.positionList = new TPositionType[Import->_n];

	_ImportList.push_back(info);

	std::cout << "registerDataImport: linked DataExport '" << this->name() <<"' with DataImport '" << Import->name() << "'." << std::endl;
	return true;
}

template<typename TDataType, typename TPositionType>
bool DataExport<TDataType, TPositionType>::unregisterDataImport(DataImport<TDataType, TPositionType>* Import)
{
	int i;
	for(i = 0; i < _ImportList.size(); ++i)
	{
		if(_ImportList[i].import == Import) break;
	}
	if(i >= _ImportList.size())
	{
		std::cout << "unregisterDataImport: DataImport not found. Register DataImport first." << std::endl;
		return false;
	}

	delete[] _ImportList[i].valueList;
	delete[] _ImportList[i].positionList;

	_ImportList.erase(_ImportList.begin() + i);

	std::cout << "unregisterDataImport: DataExport '" << this->name() <<"' no longer linked to DataImport '" << Import->name() << "'." << std::endl;
	return true;
}

template<typename TDataType, typename TPositionType>
bool DataExport<TDataType, TPositionType>::updateDataImport(DataImport<TDataType, TPositionType>* Import)
{
	int i;

	for(i = 0; i < _ImportList.size(); ++i)
	{
		if(_ImportList[i].import == Import) break;
	}
	if(i >= _ImportList.size())
	{
		std::cout << "updateDataImport: DataImport not found. Register DataImport first." << std::endl;
		return false;
	}

	// check for reallocation
	if(Import->_n != _ImportList[i].n)
	{
		delete[] _ImportList[i].valueList;
		_ImportList[i].valueList = new TDataType[Import->_n];
		delete[] _ImportList[i].positionList;
		_ImportList[i].positionList = new TPositionType[Import->_n];
		_ImportList[i].n = Import->_n;
	}

	// copy positions
	for(int j = 0; j < _ImportList[i].n; ++j)
	{
		(_ImportList[i].positionList)[j] = (Import->_positionList)[j];
	}

	return true;
}

template<typename TDataType, typename TPositionType>
DataImport<TDataType, TPositionType>::DataImport(std::string name) : DataImportItem(name)
{
	_DataType = &typeid(TDataType);
	_PositionType = &typeid(TPositionType);

	// reset arrays
	_n = 0;
	_valueList = 0;
	_positionList = 0;
	_Export = 0;
};

template<typename TDataType, typename TPositionType>
inline TDataType& DataImport<TDataType, TPositionType>::operator[](size_t i)
{
	assert( i<_n && "ERROR");
	return _valueList[i];
}

template<typename TDataType, typename TPositionType>
inline const TDataType& DataImport<TDataType, TPositionType>::operator[](size_t i) const
{
	assert( i<_n && "ERROR");
	return _valueList[i];
}

template<typename TDataType, typename TPositionType>
void DataImport<TDataType, TPositionType>::setPositions(TPositionType* posList, int n)
{
	//check for memory reallocation
	if(n != _n)
	{
		if(_valueList) delete[] _valueList;
		_valueList = new TDataType[n];
		if(_positionList) delete[] _positionList;
		_positionList = new TPositionType[n];
		_n = n;
	}

	// write new positions
	for(int i = 0; i < _n; ++i)
	{
		_positionList[i] = posList[i];
	}

	// communicate positions to export (if set)
	if(_Export)
	{
		_Export->updateDataImport(this);
	}
}

template<typename TDataType, typename TPositionType>
bool DataImport<TDataType, TPositionType>::registerDataExport(DataExport<TDataType, TPositionType>* Export)
{
	if(_Export)
	{
		std::cout << "setDataExport: Only one DataExport allowed for one DataImport." << std::endl;
		return false;
	}

	// remember export
	_Export = Export;

	// register at export
	_Export->registerDataImport(this);

	// adapt export size
	_Export->updateDataImport(this);

	return true;
}

template<typename TDataType, typename TPositionType>
bool DataImport<TDataType, TPositionType>::registerDataExport(DataExportItem *Export)
{
	typename ug::DataExport<TDataType, TPositionType>* Cast_Export = dynamic_cast< typename ug::DataExport<TDataType, TPositionType>* >(Export);
	return registerDataExport(Cast_Export);
}


template<typename TDataType, typename TPositionType>
bool DataImport<TDataType, TPositionType>::unregisterDataExport()
{
	// adapt export size
	if(_Export)
	{
		_Export->unregisterDataImport(this);
	}

	// reset export
	_Export = 0;

	return true;
}


template<typename TDataType, typename TPositionType>
DataImport<TDataType, TPositionType>::~DataImport()
{
	if(_valueList) delete[] _valueList;
	if(_positionList) delete[] _positionList;

	unregisterDataExport();
}



template <typename TDataType, typename TPositionType>
bool DataContainer::registerItem(DataImport<TDataType, TPositionType>& ImportItem)
{
	ImportItemList::iterator iter;
	iter = find(_importItemList.begin(), _importItemList.end(), dynamic_cast<DataImportItem*>(&ImportItem));
	if(iter != _importItemList.end())
	{
		std::cout << "Container already contains DataImportItem. Can not add." << std::endl;
		return false;
	}

	_importItemList.push_back(&ImportItem);
	return true;
}

template <typename TDataType, typename TPositionType>
bool DataContainer::registerItem(DataExport<TDataType, TPositionType>& ExportItem)
{
	ExportItemList::iterator iter;
	iter = find(_exportItemList.begin(), _exportItemList.end(), dynamic_cast<DataExportItem*>(&ExportItem));
	if(iter != _exportItemList.end())
	{
		std::cout << "Container already contains DataExportItem. Can not add." << std::endl;
		return false;
	}
	_exportItemList.push_back(&ExportItem);
	return true;
}

template <typename TDataType, typename TPositionType>
bool DataContainer::registerItem(DataLinker<TDataType, TPositionType>& LinkerItem)
{
	_linkerItemList.push_back(&LinkerItem);
	return true;
}

template <typename TDataType, typename TPositionType>
bool DataContainer::link(DataImport<TDataType, TPositionType>& ImportItem, DataExport<TDataType, TPositionType>& ExportItem)
{
	std::vector<DataImportItem*>::iterator importIter = _importItemList.begin();
	std::vector<DataExportItem*>::iterator exportIter = _exportItemList.begin();

	importIter = find(_importItemList.begin(), _importItemList.end(), dynamic_cast<DataImportItem*>(&ImportItem));
	if(importIter == _importItemList.end())
	{
		std::cout << "Container does not contain DataImportItem. Can not link." << std::endl;
		return false;
	}
	exportIter = find(_exportItemList.begin(), _exportItemList.end(), dynamic_cast<DataExportItem*>(&ExportItem));
	if(exportIter == _exportItemList.end())
	{
		std::cout << "Container does not contain DataExportItem. Can not link." << std::endl;
		return false;
	}

	//ok, we can link
	//!TODO: Check dependency
	_neddedExportItemList.push_back(&ExportItem);

	ImportItem.registerDataExport(&ExportItem);
	return true;
}


} // End namespace ug

#endif
