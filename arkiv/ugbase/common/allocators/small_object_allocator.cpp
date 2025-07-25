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

#include "small_object_allocator.h"


FixedAllocator::
FixedAllocator(std::size_t blockSize, unsigned char numBlocksPerChunk) :
	m_blockSize(blockSize),
	m_numBlocksPerChunk(numBlocksPerChunk),
	m_allocChunk(0),
	m_emptyChunkIndex(-1),
	m_deallocChunkIndex(-1),
	m_numFreeBlocks(0)
{
}

void* FixedAllocator::
allocate()
{
	if(m_numFreeBlocks == 0){
		//m_chunks.reserve(m_chunks.size() + 1);
		Chunk newChunk;
		newChunk.init(m_blockSize, m_numBlocksPerChunk);
		m_chunks.push_back(newChunk);
		m_allocChunk = &m_chunks.back();
		m_deallocChunkIndex = m_chunks.size() - 1;
		m_numFreeBlocks += m_numBlocksPerChunk;
	}
	else if(m_allocChunk->m_numAvailableBlocks == 0)
	{
		Chunks::iterator iter = m_chunks.begin();
		for(;iter != m_chunks.end(); ++iter){
			if(iter->m_numAvailableBlocks > 0){
				m_allocChunk = &*iter;
				break;
			}
		}
	}
	
	assert(m_allocChunk != 0);
	assert(m_allocChunk->m_numAvailableBlocks > 0);
	
	--m_numFreeBlocks;
	return m_allocChunk->allocate(m_blockSize);;
}

void FixedAllocator::
deallocate(void* p)
{
	assert(m_deallocChunkIndex != -1);
		
	if(!pointer_is_in_chunk(p, &m_chunks[m_deallocChunkIndex]))
	{
	//	find the chunk in which the pointer lies
		int iDown = (int)m_deallocChunkIndex - 1;
		std::size_t iUp = m_deallocChunkIndex + 1;
	
		for(;;){
			if(iDown >= 0){
				if(pointer_is_in_chunk(p, &m_chunks[iDown])){
					m_deallocChunkIndex = (std::size_t)iDown;
					break;
				}
				--iDown;
			}
			
			if(iUp < m_chunks.size()){
				if(pointer_is_in_chunk(p, &m_chunks[iUp])){
					m_deallocChunkIndex = iUp;
					break;
				}
				++iUp;
			}
			
			if(iDown < 0 && iUp == m_chunks.size()){
			//	this can only happen if the pointer was not allocated by
			//	this instance of FixedAllocator.
				m_deallocChunkIndex = -1;
				break;
			}
		}
	}
	
	assert(m_deallocChunkIndex != -1);
	
	Chunk* deallocChunk = &m_chunks[m_deallocChunkIndex];
	deallocChunk->deallocate(p, m_blockSize);
	if(deallocChunk->m_numAvailableBlocks == m_numBlocksPerChunk){
	//	the chunk is empty now. if we already have an empty chunk,
	//	we'll erase it right away and set m_deallocChunk to the middle of
	//	all chunks.
	//	if not, we'll set m_emptyChunkIndex to m_deallocChunkIndex
		if(m_emptyChunkIndex == -1){
			m_emptyChunkIndex = m_deallocChunkIndex;
		}
		else{
		//	swap deallocChunk with the last chunk and erase it afterwards.
			if(m_deallocChunkIndex != (int)m_chunks.size() - 1)
			{
				std::swap(m_chunks[m_deallocChunkIndex], m_chunks.back());
			}
			if(m_chunks.size() > 1){
				m_chunks[m_chunks.size() - 1].free();
				m_chunks.erase(m_chunks.begin() + (m_chunks.size() - 1));
				m_numFreeBlocks -= m_numBlocksPerChunk;
			}
				
		//	find the new index for deallocChunk
			m_deallocChunkIndex = m_chunks.size() / 2;
		}
	}
	++m_numFreeBlocks;
}


void FixedAllocator::Chunk::
init(std::size_t blockSize, unsigned char numBlocks)
{
	m_pData = new unsigned char[blockSize * numBlocks];
	m_firstAvailableBlock = 0;
	m_numAvailableBlocks = numBlocks;
	
	unsigned char* p = m_pData;
	for(unsigned char i = 0; i != numBlocks; p+=blockSize)
		*p = ++i;
}

void FixedAllocator::Chunk::
free()
{
	delete[] m_pData;
}

void* FixedAllocator::Chunk::
allocate(std::size_t blockSize)
{
//	if(!m_numAvailableBlocks)
//		return 0;
		
	unsigned char* pResult = m_pData + (m_firstAvailableBlock * blockSize);
	m_firstAvailableBlock = *pResult;
	--m_numAvailableBlocks;

	return pResult;
}

void FixedAllocator::Chunk::
deallocate(void* p, std::size_t blockSize)
{
	assert(p >= m_pData);
	unsigned char* toRelease = static_cast<unsigned char*>(p);
	assert((toRelease - m_pData) % blockSize == 0);
	*toRelease = m_firstAvailableBlock;
	m_firstAvailableBlock = static_cast<unsigned char>((toRelease - m_pData) / blockSize);
	assert(m_firstAvailableBlock == (toRelease - m_pData) / blockSize);
	++m_numAvailableBlocks;
}
