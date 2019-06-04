/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Markus Breit
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

#ifndef __H__UG_file_io_swc
#define __H__UG_file_io_swc

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/common_attachments.h"

namespace ug {

namespace swc_types
{
	enum swc_type
	{
		SWC_UNDF   = 0,
		SWC_SOMA   = 1,
		SWC_AXON   = 2,
		SWC_DEND   = 3,
		SWC_APIC   = 4,
    SWC_FORK   = 5,
    SWC_END    = 6,
    SWC_CUSTOM = 7
	};

	struct SWCPoint
	{
		vector3 coords;
		number radius;
		swc_type type;
		std::vector<size_t> conns;
	};
}


class FileReaderSWC
{
	public:
		typedef swc_types::SWCPoint SWCPoint;

	public:
		FileReaderSWC() {};
		~FileReaderSWC() {};

		bool load_file(const char* fileName);
		bool create_grid(Grid& g, ISubsetHandler* pSH, number scale_length = 1.0);

		const std::vector<swc_types::SWCPoint>& swc_points() const;
		std::vector<swc_types::SWCPoint>& swc_points();

	protected:
		std::vector<swc_types::SWCPoint> m_vPts;
};



class FileWriterSWC
{
	public:
		FileWriterSWC() {};
		~FileWriterSWC() {};

		bool export_grid_to_file(Grid& grid, ISubsetHandler* pSH, const char* filename);
};


bool LoadGridFromSWC(Grid& grid, ISubsetHandler* pSH, const char* filename, AVector3& aPos = aPosition);
bool ExportGridToSWC(Grid& grid, ISubsetHandler* pSH, const char* filename, AVector3& aPos = aPosition);

} // end of namespace ug

#endif // __H__UG_file_io_swc
