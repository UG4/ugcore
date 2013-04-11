/**
 * \file /pcl_const_rework/ugbase/pcl/pcl_tostring.h
 *
 *  \date	 	05.03.2013
 *  \author 	mrupp
 */

#include <string>
#include <sstream>

#ifndef __H__PCL_TOSTRING_H__
#define __H__PCL_TOSTRING_H__


namespace pcl{

/// \addtogroup pcl
/// \{

inline std::string ToString(const ProcessCommunicator &pc)
{
	if(pc.empty()) return "Empty ProcessCommunicator";
	else if(pc.is_world()) return "MPI_COMM_WORLD ProcessCommunicator";
	else
	{
		std::stringstream out;
		out << "ProcessCommunicator (size = " << pc.size() << "): [";
		for(size_t i=0; i<pc.size(); i++)
			out << pc.get_proc_id(i) << " ";
		out << "] ";
		return out.str();
	}
}

// end group pcl
/// \}


} // namespace ug

#endif /* __UG__PCL_TOSTRING_H__ */
