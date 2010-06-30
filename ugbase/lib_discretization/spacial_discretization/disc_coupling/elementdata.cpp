/*
 * elementdata.cpp
 *
 *  Created on: 30.06.2009
 *      Author: andreasvogel
 */

#include "elementdata.h"

namespace ug{

DataExportItem::DataExportItem(std::string name)
{
	_name = name;
}

std::string DataExportItem::name()
{
	return _name;
}

const std::type_info* DataExportItem::DataType() const
{
	return _DataType;
}

const std::type_info* DataExportItem::PositionType() const
{
	return _PositionType;
}

DataImportItem::DataImportItem(std::string name)
{
	_name = name;
}

std::string DataImportItem::name()
{
	return _name;
}

const std::type_info* DataImportItem::DataType() const
{
	return _DataType;
}

const std::type_info* DataImportItem::PositionType() const
{
	return _PositionType;
}


bool DataContainer::linkByNumber(int ImportNr, int ExportNr)
{
	// check if linking is possible
	if(ImportNr >= _importItemList.size()) return false;
	if(ExportNr >= _exportItemList.size()) return false;

	if(_importItemList[ImportNr]->DataType() != _exportItemList[ExportNr]->DataType())
	{
		std::cout << "Datatype of Import and Export Parameter does not match" << std::endl;
		return false;
	}
	if(_importItemList[ImportNr]->PositionType() != _exportItemList[ExportNr]->PositionType())
	{
		std::cout << "PositionType of Import and Export Parameter does not match" << std::endl;
		return false;
	}

	//ok, we can link them
	//! TODO: Check dependency
	_neddedExportItemList.push_back(_exportItemList[ExportNr]);

	_importItemList[ImportNr]->registerDataExport(_exportItemList[ExportNr]);
	return true;
}

bool DataContainer::compute()
{
	for(int i=0; i < _neddedExportItemList.size(); ++i)
	{
		_neddedExportItemList[i]->compute();
	}
}


bool DataContainer::printInfo()
{
	for(int i=0; i < _exportItemList.size(); ++i)
	{
		std::cout << "ExportItem " << i << ": " << _exportItemList[i]->name();
		std::cout << "  [DataType: " << _exportItemList[i]->DataType()->name();
		std::cout << ",  PositionType: " << _exportItemList[i]->PositionType()->name() << "]" << std::endl;
	}
	for(int i=0; i < _importItemList.size(); ++i)
	{
		std::cout << "ImportItem " << i << ": " << _importItemList[i]->name();
		std::cout << "  [DataType: " << _importItemList[i]->DataType()->name();
		std::cout << ",  PositionType: " << _importItemList[i]->PositionType()->name() << "]" << std::endl;
	}
}




}
