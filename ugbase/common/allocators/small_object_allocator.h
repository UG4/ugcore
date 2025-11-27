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

#ifndef __H__SMALL_OBJECT_ALLOCATOR__
#define __H__SMALL_OBJECT_ALLOCATOR__

// #include <cassert>
#include <vector>

/**	Instances of this class can be used to allocate small objects of the same size
 *	in a highly efficient way.
 */
class FixedAllocator
{
	public:
		FixedAllocator(std::size_t blockSize, unsigned char numBlocksPerChunk);
		void* allocate();
		void deallocate(void* p);
		
	private:
	/**	The Chunk structure has to be used with great care.
	 *	A represents is a block of memory, in which it allocates
	 *	small objects. Be sure that the blockSize used in the calls
	 *	is the same for all calls to an instance.
	 */
		struct Chunk
		{
		///	be careful. numBlocks has to be <= 255.
			void init(std::size_t blockSize, unsigned char numBlocks);

		///	call this method instead of a destructor
			void free() const;
			
		///	returns 0 if no more blocks are available.
			void* allocate(std::size_t blockSize);
		///	deallocates the given pointer.
		/**	Make sure that p was allocated by the same instance on which you
		 *	call deallocate.*/
			void deallocate(void* p, std::size_t blockSize);
			
			unsigned char* m_pData;
			unsigned char m_firstAvailableBlock;
			unsigned char m_numAvailableBlocks;
		};

	private:
		inline bool pointer_is_in_chunk(void* p, Chunk* chunk) const {
			return (p >= chunk->m_pData)
				   && (p < chunk->m_pData + m_blockSize * m_numBlocksPerChunk);
		}
		
	private:
		using Chunks = std::vector<Chunk>;
		
	private:
		std::size_t m_blockSize;
		unsigned char m_numBlocksPerChunk;
		Chunks m_chunks;
		Chunk* m_allocChunk;
		int m_emptyChunkIndex;
		int m_deallocChunkIndex;
		std::size_t m_numFreeBlocks;
};

/**	A singleton that can be used to allocate small objects.*/
template <std::size_t maxObjSize = 64, std::size_t maxChunkSize = 4096>
class SmallObjectAllocator
{
	public:
	///	returns an instance to this singleton
		static SmallObjectAllocator& inst();
		
	///	if numBytes > maxObjSize, allocate will directly call new.
		void* allocate(std::size_t numBytes);
		
	///	make sure that size exactly specifies the number of bytes of the object to which p points.
		void deallocate(void* p, std::size_t size);
		
	private:
		SmallObjectAllocator();
		
	private:
		std::vector<FixedAllocator>	m_allocators;
};


/**	This class implements the operators new and delete, so that they use the
 *	SmallObjectAllocator.
 *	By deriving from this class, your objects will be allocated through the
 *	SmallObjectAllocator too.
 */
template <std::size_t maxObjSize = 64, std::size_t maxChunkSize = 4096>
class SmallObject
{
	public:
		static void* operator new(std::size_t size);
		static void operator delete(void* p, std::size_t size);
		virtual ~SmallObject() = default;
};

////////////////////////////////
//	include implementation
#include "small_object_allocator_impl.h"

#endif
