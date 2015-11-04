//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m12 d03

#ifndef __H__SMALL_OBJECT_ALLOCATOR_IMPL__
#define __H__SMALL_OBJECT_ALLOCATOR_IMPL__

template <std::size_t maxObjSize, std::size_t maxChunkSize>
SmallObjectAllocator<maxObjSize, maxChunkSize>&
SmallObjectAllocator<maxObjSize, maxChunkSize>::
inst()
{
	static SmallObjectAllocator<maxObjSize, maxChunkSize> alloc;
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
