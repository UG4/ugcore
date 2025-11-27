/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef __H__TAR_H_
#define __H__TAR_H_


#include <cstdio>
#include <cstring>
#include <ctime>
#include <string>

namespace ug{
struct TarHeader
{
	char filename[100];
	char filemode[8];
	char userID[8];
	char groupID[8];
	char octalFileSize[12];
	char octalModificationTimeStamp[12];
	char checksum[8];
	char linkIndicator;
	char linkedFilename[100];

	char ustarIndicator[6];
	char ustarVersion[2];
	char ustarOwnerUserName[32];
	char ustarOwnerGroupName[32];
	char ustarDeviceMajorNumber[8];
	char ustarDeviceMinorNumber[8];
	char ustarFilenamePrefix[155];
	char padding[12];


	//char padding[512 - (100+8+8+8+12+12+8+1+100)];

	TarHeader()
	{
		memset(this, 0, sizeof(TarHeader));
		strcpy(filemode, "000644 ");
		strcpy(userID, 	"000765 ");
		strcpy(groupID, 	"000024 ");
		strcpy(ustarIndicator, "ustar");
		memcpy(ustarVersion, "00", 2);
		strcpy(ustarOwnerUserName, "");
		strcpy(ustarOwnerGroupName, "");
		strcpy(ustarDeviceMajorNumber, "000000 ");
		strcpy(ustarDeviceMinorNumber, "000000 ");

		sprintf(octalModificationTimeStamp, "%o", static_cast<unsigned int>(time(nullptr)));
//		strncpy(octalModificationTimeStamp, "12253557334 ", 12);
		linkIndicator = '0';
	}
	void set_filename(std::string name)
	{
		strcpy(filename, name.c_str());
	}
	void set_filesize(size_t size)
	{
		char buf[13];
		sprintf(buf, "%011o ", static_cast<unsigned int>(size));
		memcpy(octalFileSize, buf, 12);
	}

	void set_checksum()
	{
		memset(checksum, 32, sizeof(checksum));
		auto p = reinterpret_cast<unsigned char *>(this);
		unsigned int cs = 0;
		for(size_t i=0; i<sizeof(TarHeader); i++)
			cs += p[i];
		sprintf(checksum, "%06o", ((unsigned int)cs));
		//			 std:: cout << "checksum = " << checksum << "\n";
	}
};

}
#endif