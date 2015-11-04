#ifndef __H__SMALL_OBJECT_ALLOCATOR__
#define __H__SMALL_OBJECT_ALLOCATOR__

#include <cassert>
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
			void free();
			
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
		inline bool pointer_is_in_chunk(void* p, Chunk* chunk)
		{
			return (p >= chunk->m_pData)
				   && (p < chunk->m_pData + m_blockSize * m_numBlocksPerChunk);
		}
		
	private:
		typedef std::vector<Chunk> Chunks;
		
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
		virtual ~SmallObject()	{}
};

////////////////////////////////
//	include implementation
#include "small_object_allocator_impl.h"

#endif
