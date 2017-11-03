/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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


#include "common/util/binary_buffer.h"
#include "common/serialization.h"

namespace pcl {

template <typename TKey, typename TValue, typename Compare>
void MinimalKeyValuePairAcrossAllProcs(TKey& keyInOut, TValue& valInOut, const Compare& cmp)
{
	// in the serial case, input key and value are already output key and value
#ifndef UG_PARALLEL
	return;
#endif
	size_t nProc = NumProcs();
	if (nProc == 1)
		return;

	// write key/value pair into binary buffer
	ug::BinaryBuffer buf;
	ug::Serialize(buf, keyInOut);
	ug::Serialize(buf, valInOut);

	// gather buffers on root
	ProcessCommunicator pc;
	pc.gather(buf, 0);

	// find minimum on root
	size_t rk = ProcRank();
	if (rk == 0)
	{
		TKey key;
		TValue val;

		for (size_t i = 0; i < nProc; ++i)
		{
			ug::Deserialize(buf, key);
			ug::Deserialize(buf, val);

			if (cmp(key, keyInOut))
			{
				keyInOut = key;
				valInOut = val;
			}
		}

		buf.clear();
		ug::Serialize(buf, keyInOut);
		ug::Serialize(buf, valInOut);
	}

	// communicate minimum to all procs
	pc.broadcast(buf, 0);

	// copy to return values
	ug::Deserialize(buf, keyInOut);
	ug::Deserialize(buf, valInOut);
}


} // end of namespace pcl
