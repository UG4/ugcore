/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_file_io_asc
#define __H__UG_file_io_asc

#include "common/util/field.h"
#include "common/math/ugmath_types.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/field_util.h"

namespace ug{

class FileReaderASC{
	public:
		FileReaderASC();
		~FileReaderASC();

	///	set an external field in which the data will be loaded
	/**	By default data is loaded into the internal field of the FileReader.
	 * If you want to provide a field into which the data shall be loaded,
	 * e.g. to avoid memory duplication and copying, you may use this method.
	 * In this case you should make sure that the provided field instance outlives
	 * the FileReaders instance.
	 * Call this method with nullptr to indicate that the internal field shall be used again.
	 * \note When an external field-pointer is set or unset, no data will be copied
	 * 		 between the internal and external fields. Unless a new load_file
	 *		 has been executed, the provided data may not be of any use.*/
		void set_field(Field<number>* field);

		void load_file(const char* filename);

		number cell_size() const			{return m_cellSize;}
		number lower_left_corner_x() const	{return m_llcorner.x();}
		number lower_left_corner_y() const	{return m_llcorner.y();}
		number no_data_value() const		{return m_noDataValue;}
		size_t num_rows() const				{return m_field->height();}
		size_t num_columns() const			{return m_field->width();}

		number at(size_t r, size_t c) const	{return m_field->at(r, c);}
		const Field<number>& field() const	{return *m_field;}

	private:
		Field<number>*	m_field;
		number			m_cellSize;
		vector2			m_llcorner;
		number			m_noDataValue;
		SmartPtr<Field<number> >	m_privateField;
};

bool LoadGridFromASC(Grid& grid, const char* filename, AVector3& aPos = aPosition);

}//	end of namespace

#endif