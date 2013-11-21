#include "pcl/pcl_process_communicator.h"
#include "algebra_layouts.h"

#include "common/util/string_util.h"


namespace ug
{


std::ostream &operator << (std::ostream &out, const HorizontalAlgebraLayouts &layouts)
{
	out << "HorizontalAlgebraLayouts:\n";
	out << " master: " << OstreamShift(layouts.master()) << "\n";
	out << " slave: " << OstreamShift(layouts.slave()) << "\n";
	out << " process communicator: " << OstreamShift(layouts.proc_comm()) << "\n";
}



std::ostream &operator << (std::ostream &out, const AlgebraLayouts &layouts)
{
	out << "AlgebraLayouts:\n";
	out << " master: " << OstreamShift(layouts.master()) << "\n";
	out << " slave: " << OstreamShift(layouts.slave()) << "\n";
	out << " vertical master: " << OstreamShift(layouts.master()) << "\n";
	out << " vertical slave: " << OstreamShift(layouts.slave()) << "\n";
	out << " process communicator: " << OstreamShift(layouts.proc_comm()) << "\n";
}


}
