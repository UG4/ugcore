#include <string>
#include <fstream>
#include "pcl/pcl.h"

namespace ug
{

std::string GetParallelName(std::string name, const pcl::ProcessCommunicator &pc, bool bWriteHeader)
{
	char buf[20];
	int rank = pcl::GetProcRank();

	size_t iExtPos = name.find_last_of(".");
	std::string ext = name.substr(iExtPos+1);
	name.resize(iExtPos);
	if(bWriteHeader && rank == pc.get_proc_id(0))
	{
		size_t iSlashPos = name.find_last_of("/");
		if(iSlashPos == std::string::npos) iSlashPos = 0; else iSlashPos++;
		std::string name2 = name.substr(iSlashPos);
		std::fstream file((name+".p"+ext).c_str(), std::ios::out);
		file << pc.size() << "\n";
		for(size_t i=0; i<pc.size(); i++)
		{
			sprintf(buf, "_p%04d.%s", pc.get_proc_id(i), ext.c_str());
			file << name2 << buf << "\n";
		}
	}

	sprintf(buf, "_p%04d.%s", rank, ext.c_str());
	name.append(buf);
	return name;
}


}
