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

#ifndef __H__SMALL_OBJECT_ALLOCATOR_IMPL__
#define __H__SMALL_OBJECT_ALLOCATOR_IMPL__

#include "small_object_allocator.h"

template <std::size_t maxObjSize, std::size_t maxChunkSize>
SmallObjectAllocator<maxObjSize, maxChunkSize>&
SmallObjectAllocator<maxObjSize, maxChunkSize>::
inst()
{
	static SmallObjectAllocator alloc;
	return alloc;
}

template <std::size_t maxObjSize, std::size_t maxChunkSize>
void* SmallObjectAllocator<maxObjSize, maxChunkSize>::
allocate(std::size_t numBytes)
{
	if(numBytes > maxObjSize)
		return new unsigned char[numBytes];
	
	return m_allocators[numBytes].allocate();
}

template <std::size_t maxObjSize, std::size_t maxChunkSize>
void SmallObjectAllocator<maxObjSize, maxChunkSize>::		
deallocate(void* p, std::size_t size)
{
	if(size > maxObjSize)
		delete[] static_cast<unsigned char*>(p);
	else{
		m_allocators[size].deallocate(p);
	}
}
		
template <std::size_t maxObjSize, std::size_t maxChunkSize>
SmallObjectAllocator<maxObjSize, maxChunkSize>::
SmallObjectAllocator()
{
//	initialize the allocators
	for(std::size_t i = 0; i <= maxObjSize; ++i){
		std::size_t maxNumBlocks = maxChunkSize / maxObjSize;
		if(maxNumBlocks > 255)
			maxNumBlocks = 255;
		m_allocators.push_back(FixedAllocator(maxObjSize,
									static_cast<unsigned char>(maxNumBlocks)));
	}
}


template <std::size_t maxObjSize, std::size_t maxChunkSize>
void* SmallObject<maxObjSize, maxChunkSize>::
operator new(std::size_t size)
{
	return SmallObjectAllocator<maxObjSize, maxChunkSize>::inst().allocate(size);
}

template <std::size_t maxObjSize, std::size_t maxChunkSize>
void SmallObject<maxObjSize, maxChunkSize>::
operator delete(void* p, std::size_t size)
{
	SmallObjectAllocator<maxObjSize, maxChunkSize>::inst().deallocate(p, size);
}

#endif
