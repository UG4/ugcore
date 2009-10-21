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

#include "ShinyManager.h"

#include <fstream>
#include <memory.h>
#include <stdio.h>

#if SHINY_PROFILER == TRUE
namespace Shiny {


//-----------------------------------------------------------------------------

	ProfileManager ProfileManager::instance = {
		/* _lastTick = */ 0,
		/* _curNode = */ &instance.rootNode,
		/* _tableMask = */ 0,
		/* _nodeTable = */ ProfileManager::_dummyNodeTable,
#if SHINY_PROFILER_LOOKUPRATE == TRUE
		/* _lookupCount = */ 0,
		/* _lookupSuccessCount = */ 0,
#endif
		/* _tableSize = */ 1,
		/* nodeCount = */ 1,
		/* zoneCount = */ 1,
		/* _lastZone = */ &instance.rootZone,
		/* _lastNodePool = */ NULL,
		/* _firstNodePool = */ NULL,
		/* rootNode = */ {
			/* _last = */ { 0, 0 },
			/* zone = */ &instance.rootZone,
			/* parent = */ &instance.rootNode,
			/* nextSibling = */ NULL,
			/* firstChild = */ NULL,
			/* lastChild = */ NULL,
			/* childCount = */ 0,
			/* entryLevel = */ 0,
			/* _cache = */ NULL,
			/* data = */ { { 0, 0 }, { 0, 0 }, { 0, 0 } }
		},
		/* rootZone = */ {
			/* next = */ NULL,
			/* _state = */ ProfileZone::STATE_HIDDEN,
			/* name = */ "<root>",
			/* data = */ { { 0, 0 }, { 0, 0 }, { 0, 0 } }
		},
		/* _initialized = */ false,
		/* _firstUpdate = */ true
	};

	ProfileNode* ProfileManager::_dummyNodeTable[] = { NULL };


//-----------------------------------------------------------------------------

	/* Robert Jenkins' 32 bit integer hash function

	SHINY_INLINE uint32_t hash_value(ProfileNode* a_pParent, ProfileZone* a_pZone) {
		uint32_t a = ptr32(a_pParent) + ptr32(a_pZone);

		a = (a+0x7ed55d16) + (a<<12);
		a = (a^0xc761c23c) ^ (a>>19);
		a = (a+0x165667b1) + (a<<5);
		a = (a+0xd3a2646c) ^ (a<<9);
		a = (a+0xfd7046c5) + (a<<3);
		a = (a^0xb55a4f09) ^ (a>>16);
		return a;
	}
	*/

	/* Old hash function
	
	SHINY_INLINE uint32_t hash_index(ProfileNode* a_pParent, ProfileZone* a_pZone) {
		uint32_t a = ptr32(a_pParent) + ptr32(a_pZone);
		return (a << 8) - (a >> 4);
	}
	*/

	// primary hash function
	SHINY_INLINE uint32_t hash_value(ProfileNode* a_pParent, ProfileZone* a_pZone) {
		uint32_t a = ptr32(a_pParent) + ptr32(a_pZone);

		a = (a+0x7ed55d16) + (a<<12);
		a = (a^0xc761c23c) ^ (a>>19);
		return a;
	}

	// secondary hash used as index offset: force it to be odd
	// so it's relatively prime to the power-of-two table size
	SHINY_INLINE uint32_t hash_offset(uint32_t a) {
		return ((a << 8) + (a >> 4)) | 1;
	}


//-----------------------------------------------------------------------------

	void ProfileManager::preLoad(void) {
		if (!_initialized) {
			_init();

			_createNodeTable(TABLE_SIZE_INIT);
			_createNodePool(TABLE_SIZE_INIT / 2);
		}
	}


//-----------------------------------------------------------------------------

	void ProfileManager::update(float a_damping) {
		_appendTicksToCurNode();

		if (!_firstUpdate) {
			rootZone.preUpdateChain();
			rootNode.updateTree(a_damping);
			rootZone.updateChain(a_damping);

		} else {
			_firstUpdate = false;
			rootZone.preUpdateChain();
			rootNode.updateTree(0);
			rootZone.updateChain(0);
		}
	}


//-----------------------------------------------------------------------------

	void ProfileManager::clear(void) {
		destroy();
		preLoad();
	}


//-----------------------------------------------------------------------------

	void ProfileManager::destroy(void) {
		_resetZones();
		_destroyNodes();
		_uninit();
	}


//-----------------------------------------------------------------------------

	ProfileNode* ProfileManager::_lookupNode(ProfileNodeCache* a_cache, ProfileZone* a_zone) {
		uint32_t nHash = hash_value(_curNode, a_zone);
		uint32_t nIndex = nHash & _tableMask;
		ProfileNode* pNode = _nodeTable[nIndex];

		_incLookup();
		_incLookupSuccess();

		if (pNode) {
			if (pNode->isEqual(_curNode, a_zone)) return pNode; // found it!
			
			// hash collision:

			// compute a secondary hash function for stepping
			uint32_t nStep = hash_offset(nHash);

			for (;;) {
				_incLookup();

				nIndex = (nIndex + nStep) & _tableMask;
				pNode = _nodeTable[nIndex];

				if (!pNode) break;
				else if (pNode->isEqual(_curNode, a_zone)) return pNode;
			}

			// loop is guaranteed to end because the hash table is never full
		}

		if (!a_zone->isInited()) { // zone is not initialized
			a_zone->init(_lastZone);

			_lastZone = a_zone;
			zoneCount++;

			if (_initialized == false) { // first time init
				_init();

				_createNodeTable(TABLE_SIZE_INIT);
				_createNodePool(TABLE_SIZE_INIT / 2);

				// initialization has invalidated nIndex
				// we must compute nIndex again
				return _createNode(a_cache, a_zone);
			}
		}

		// YES nodeCount is not updated
		// but it includes rootNode so it adds up.

		// check if we need to grow the table
		// we keep it at most 1/2 full to be very fast
		if (_tableSize < 2 * nodeCount) {

			_resizeNodeTable(2 * _tableSize);
			_resizeNodePool(nodeCount - 1);

			// expansion has invalidated nIndex
			// we must compute nIndex again
			return _createNode(a_cache, a_zone);
		}
		
		nodeCount++;

		ProfileNode* pNewNode = _lastNodePool->newItem();
		pNewNode->init(_curNode, a_zone, a_cache);

		_nodeTable[nIndex] = pNewNode;
		return pNewNode;
	}


//-----------------------------------------------------------------------------

	ProfileNode* ProfileManager::_createNode(ProfileNodeCache* a_cache, ProfileZone* a_pZone) {
		ProfileNode* pNewNode = _lastNodePool->newItem();
		pNewNode->init(_curNode, a_pZone, a_cache);

		nodeCount++;
		_insertNode(pNewNode);
		return pNewNode;
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_insertNode(ProfileNode* a_pNode) {
		uint32_t nHash = hash_value(a_pNode->parent, a_pNode->zone);
		uint32_t nIndex = nHash & _tableMask;

		if (_nodeTable[nIndex]) {
			uint32_t nStep = hash_offset(nHash);

			while (_nodeTable[nIndex])
				nIndex = (nIndex + nStep) & _tableMask;
		}

		_nodeTable[nIndex] = a_pNode;
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_createNodePool(uint32_t a_nCount) {
		_firstNodePool = ProfileNodePool::createNodePool(a_nCount);
		_lastNodePool = _firstNodePool;
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_resizeNodePool(uint32_t a_nCount) {
		ProfileNodePool* pPool = ProfileNodePool::createNodePool(a_nCount);
		_lastNodePool->nextPool = pPool;
		_lastNodePool = pPool;
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_createNodeTable(uint32_t a_nCount) {
		_tableSize = a_nCount;
		_tableMask = a_nCount - 1;

		_nodeTable = static_cast<ProfileNodeTable*>(
			malloc(sizeof(ProfileNode) * a_nCount));

		memset(_nodeTable, 0, a_nCount * sizeof(ProfileNode*));
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_resizeNodeTable(uint32_t a_nCount) {
		ProfileNodePool* pPool;

		free(_nodeTable);
		_createNodeTable(a_nCount);

		pPool = _firstNodePool;
		while (pPool) {

			ProfileNode *pIter = pPool->firstItem();

			while (pIter != pPool->unusedItem())
				_insertNode(pIter++);

			pPool = pPool->nextPool;
		}
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_resetZones(void) {
		ProfileZone *pZone, *pNextZone;

		pZone = &rootZone;

		for(;;) {
			pZone->uninit();

			pNextZone = pZone->next;
			pZone->next = NULL;
			
			if (!pNextZone) break;
			pZone = pNextZone;
		}

		_lastZone = &rootZone;
		zoneCount = 1;
	}


//-----------------------------------------------------------------------------

	void ProfileManager::_destroyNodes(void) {
		if (_firstNodePool) {
			_firstNodePool->destroy();
			_firstNodePool = NULL;
		}

		if (_nodeTable != instance._dummyNodeTable) {
			free(_nodeTable);

			_nodeTable = instance._dummyNodeTable;
			_tableSize = 1;
			_tableMask = 0;
		}

		_curNode = &rootNode;
		nodeCount = 1;

		_init();
	}


//-----------------------------------------------------------------------------

	bool ProfileManager::output(const char *a_filename) {
		std::ofstream file(a_filename, std::ios_base::out);

		if (!file.is_open()) return false;
		else return output(file);
	}


//-----------------------------------------------------------------------------

	bool ProfileManager::output(std::ostream &a_ostream) {
		a_ostream << outputZonesAsString().c_str()
		          << "\n\n"
		          << outputNodesAsString().c_str()
		          << "\n\n"
				  << std::flush;

		return true;
	}


} // namespace Shiny

#endif // if SHINY_PROFILER == TRUE
