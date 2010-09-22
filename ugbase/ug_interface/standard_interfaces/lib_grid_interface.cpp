//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m09 d20

#include "../ug_interface.h"
#include "../registry.h"
#include "lib_grid/lib_grid.h"
#include "lib_grid_interface.h"

using namespace std;

namespace ug{
namespace interface
{

void RegisterLibGridInterface(Registry& reg)
{
	reg.register_object<GridObject>();
	reg.register_object<MultiGridObject>();
	reg.register_object<SubsetHandlerObject>();
	reg.register_object<MGSubsetHandlerObject>();
		
	reg.register_global_function<LoadGridFunc>();
	reg.register_global_function<SaveGridFunc>();
}

}//	end of namespace 
}//	end of namespace 
