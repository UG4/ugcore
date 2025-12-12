/*
The zlib/libpng License

Copyright (c) 2007 Aidin Abedi (www.*)

This software is provided 'as-is', without any express or implied warranty. In no event will
the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial 
applications, and to alter it and redistribute it freely, subject to the following
restrictions:

    1. The origin of this software must not be misrepresented; you must not claim that 
       you wrote the original software. If you use this software in a product, 
       an acknowledgment in the product documentation would be appreciated but is 
       not required.

    2. Altered source versions must be plainly marked as such, and must not be 
       misrepresented as being the original software.

    3. This notice may not be removed or altered from any source distribution.
*/

/*
This file has been altered by Sebastian Reiter (s.b.reiter@googlemail.com).
I replaced some includes and defines, to make it easier to integrate shiny
into ug.
New code is marked with //sreiter
*/

#include "ShinyNodePool.h"
#include "ShinyTools.h"

#include <cstdlib>		//sreiter
#include <memory.h>
//#include <malloc.h>	//sreiter

#if SHINY_PROFILER == TRUE
namespace Shiny {


//-----------------------------------------------------------------------------

	ProfileNodePool* ProfileNodePool::createNodePool(uint32_t a_items) {
		auto* pPool = static_cast<ProfileNodePool*>(
			malloc(sizeof(ProfileNodePool) + sizeof(T) * (a_items - 1)));

		pPool->nextPool = nullptr;
		pPool->_nextItem = &pPool->_items[0];
		pPool->endOfItems = &pPool->_items[a_items];

		memset(&pPool->_items[0], 0, a_items * sizeof(T));
		return pPool;
	}


//-----------------------------------------------------------------------------

	uint32_t ProfileNodePool::memoryUsageChain() {
		uint32_t bytes = ptr32(
			reinterpret_cast<void*>(
				  reinterpret_cast<char*>(endOfItems)
				- reinterpret_cast<char*>(this)));

		if (nextPool) bytes += nextPool->memoryUsageChain();
		return bytes;
	}


//-----------------------------------------------------------------------------

	void ProfileNodePool::destroy() {
		T* pItem = firstItem();

		while (pItem != unusedItem())
			(pItem++)->destroy();

		if (nextPool) nextPool->destroy();
		free(this);
	}

} // namespace Shiny
#endif
