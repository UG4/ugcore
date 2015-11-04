/*
 * tar.h
 *
 *  Created on: 16.12.2013
 *      Author: mrupp
 */

#ifndef __H__TAR_H_
#define __H__TAR_H_

#include <stdio.h>
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

		sprintf(octalModificationTimeStamp, "%o", ((unsigned int)time(NULL)));
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
		sprintf(buf, "%011o ", (unsigned int)size);
		memcpy(octalFileSize, buf, 12);
	}

	void set_checksum()
	{
		memset(checksum, 32, sizeof(checksum));
		unsigned char *p = (unsigned char *) this;
		unsigned int cs = 0;
		for(size_t i=0; i<sizeof(TarHeader); i++)
			cs += p[i];
		sprintf(checksum, "%06o", ((unsigned int)cs));
		//			 std:: cout << "checksum = " << checksum << "\n";
	}
};

}
#endif /* TAR_H_ */
