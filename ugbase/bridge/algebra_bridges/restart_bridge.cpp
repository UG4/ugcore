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
#include "../util_overloaded.h"
#include "common/util/binary_buffer.h"

using namespace std;

#if UG_PARALLEL
	#include "pcl/parallel_archive.h"
	#include "pcl/parallel_file.h"
#else

namespace pcl {
void WriteCombinedParallelFile(ug::BinaryBuffer &buffer, std::string strFilename)
{
	FILE *f = fopen(strFilename.c_str(), "wb");
	int numProcs = 1;

	int mySize = buffer.write_pos();
	int myNextOffset = 2*sizeof(int) + mySize;

	fwrite(&numProcs, sizeof(numProcs), 1, f);
	fwrite(&myNextOffset, sizeof(myNextOffset), 1, f);
	fwrite(buffer.buffer(), sizeof(char), buffer.write_pos(), f);
	fclose(f);
}

void ReadCombinedParallelFile(ug::BinaryBuffer &buffer, const std::string &strFilename)
{
	FILE *f = fopen(strFilename.c_str(), "rb");
	int numProcs, myNextOffset;
	fread(&numProcs, sizeof(numProcs), 1, f);
	fread(&myNextOffset, sizeof(myNextOffset), 1, f);

	UG_COND_THROW(numProcs != 1, "checkPoint numProcs = " << numProcs << ", but running on 1 cores (and UG_SERIAL!)");
	int firstOffset = 2*sizeof(int);
	int mySize = myNextOffset - firstOffset;

	char *p = new char[mySize];
	fread(p, sizeof(char), mySize, f);

	buffer.clear();
	buffer.reserve(mySize);
	buffer.write(p, mySize);
	delete[] p;

	fclose(f);
}

}
#endif


namespace ug{

/// Serialize for ParallelVector<T>
#ifdef UG_PARALLEL
template<typename T, typename TOStream>
void Serialize(TOStream &buf, const ParallelVector<T> &v)
{
	uint t= v.get_storage_mask();
	Serialize(buf, t);
	Serialize(buf, *dynamic_cast<const T*>(&v));
}

/// Deerialize for ParallelVector<T>
template<typename T, typename TIStream>
void Deserialize(TIStream &buf, ParallelVector<T> &v)
{
	uint t = Deserialize<uint>(buf);
	v.set_storage_type(t);
	Deserialize(buf, *dynamic_cast<T*>(&v));
}
#endif


/*

 Here's an example how to use SaveToFile/ReadFromFile
 Note that they are using MPI I/O to store the BinaryBuffers from all cores into one single file.

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


See also checkpoint_util.lua and time_step_util.lua on how to generate checkpointing/debugging mechanisms with this.
 */

template<typename T>
void SaveToFile(const T &v, std::string filename)
{
	BinaryBuffer b;
	Serialize(b, v);
	pcl::WriteCombinedParallelFile(b, filename);
}

template<typename T>
void ReadFromFile(T &v, std::string filename)
{
	BinaryBuffer b;
	pcl::ReadCombinedParallelFile(b, filename);
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
 * @param grp				group for sorting of functionality
 */
template <typename TAlgebra>
static void Algebra(Registry& reg, string grp)
{
	string suffix = GetAlgebraSuffix<TAlgebra>();
	string tag = GetAlgebraTag<TAlgebra>();

	using vector_type = typename TAlgebra::vector_type;

	reg.add_function("SaveToFile", OVERLOADED_FUNCTION_PTR(void, SaveToFile<vector_type>, (const vector_type &, std::string)), grp);
	reg.add_function("ReadFromFile", OVERLOADED_FUNCTION_PTR(void, ReadFromFile<vector_type>, (vector_type &, std::string)), grp);
}

}; // end Functionality

}// end Restart

/// \addtogroup precond_bridge
void RegisterBridge_Restart(Registry& reg, string grp)
{
	grp.append("/Algebra/Restart");
	using Functionality = Restart::Functionality;

	try{
		RegisterAlgebraDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
	//reg.add_function("testomato", testomato, grp);
}

} // namespace bridge
} // namespace ug

