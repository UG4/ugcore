//	Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m05 d08

#ifndef __H__PCL__
#define __H__PCL__

#include "pcl_base.h"
#include "pcl_methods.h"
#include "pcl_communication_structs.h"
#include "pcl_communicator.h"

////////////////////////////////////////////////////////////////////////
///	this allows us to print messages to the users terminal
/**
 * if an output-processor is specified through \sa pcl::SetOutputProcRank
 * only messages from that processor will be printed to the screen.
 */
#define PCLLOG(msg) if(pcl::IsOutputProc()) {std::cout << msg; std::cout.flush();}

#endif
