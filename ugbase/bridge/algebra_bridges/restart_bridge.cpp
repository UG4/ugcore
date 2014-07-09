/*
 * restart_bridge.cpp
 *
 *  Created on: 05.12.2013
 *      Author: mrupp
 */


// extern headers
#include <iostream>
#include <sstream>
#include <string>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/util_algebra_dependent.h"

// preconditioner
#include "lib_algebra/lib_algebra.h"
#include "common/serialization.h"

#include "pcl/parallel_archive.h"
#include "pcl/parallel_file.h"

#include "../util_overloaded.h"

using namespace std;

namespace ug{


/*
 *

		restartfilename = "uNewtonSolution_T_"..step.."_p_".ug4vector"

  		restartStep = util.GetParamNumber("-restartStep", 0)
		-- apply newton solver
		if step < restartStep then
			ReadFromFile(u, restartfilename)
		else
			-- prepare newton solver
			if newtonSolver:prepare(u) == false then
				print ("Newton solver failed at step "..step.."."); exit();
			end
			if newtonSolver:apply(u) == false then
				print ("Newton solver failed at step "..step.."."); exit();
			end
			SaveToFile(u, restartfilename)
		end
 */

#ifdef UG_PARALLEL
template<typename T, class TOStream>
void Serialize(TOStream &buf, const ParallelVector<T> &v)
{
	uint t= v.get_storage_mask();
	Serialize(buf, t);
	Serialize(buf, *dynamic_cast<const T*>(&v));
}

template<typename T, class TIStream>
void Deserialize(TIStream &buf, ParallelVector<T> &v)
{
	uint t = Deserialize<uint>(buf);
	v.set_storage_type(t);
	Deserialize(buf, *dynamic_cast<T*>(&v));
}
#endif

template<typename T>
void SaveToFile(const T &v, std::string filename)
{
	BinaryBuffer b;
	Serialize(b, v);
	pcl::ParallelFileWrite(b, filename);
}

template<typename T>
void ReadFromFile(T &v, std::string filename)
{
	BinaryBuffer b;
	pcl::ParallelFileRead(b, filename);
	Deserialize(b, v);
}



namespace bridge{
namespace Restart{

/**
 * \defgroup precond_bridge Preconditioner Bridge
 * \ingroup algebra_bridge
 * \{
 */

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */

struct Functionality
{

/**
 * Function called for the registration of Algebra dependent parts.
 * All Functions and Classes depending on Algebra
 * are to be placed here when registering. The method is called for all
 * available Algebra types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

//	typedefs for this algebra
	typedef typename TAlgebra::vector_type vector_type;
	typedef typename TAlgebra::matrix_type matrix_type;


	reg.add_function("SaveToFile", OVERLOADED_FUNCTION_PTR(void, SaveToFile<vector_type>, (const vector_type &, std::string)), grp);
	reg.add_function("ReadFromFile", OVERLOADED_FUNCTION_PTR(void, ReadFromFile<vector_type>, (vector_type &, std::string)), grp);
}

}; // end Functionality

}// end Restart

/// \addtogroup precond_bridge
void RegisterBridge_Restart(Registry& reg, string grp)
{
	grp.append("/Algebra/Restart");
	typedef Restart::Functionality Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
	//reg.add_function("testomato", testomato, grp);
}

} // namespace bridge
} // namespace ug

