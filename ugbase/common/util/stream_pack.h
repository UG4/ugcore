// created by Sebastian Reiter
// y09 m11 d04
// s.b.reiter@googlemail.com

#ifndef __H__UG_STREAM_PACK__
#define __H__UG_STREAM_PACK__

#include <map>
#include <iostream>
#include "binary_stream.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class StreamPack;


////////////////////////////////////////////////////////////////////////
//	methods

///	writes a grid-pack to a binary-stream
/**
 * This method is mainly intended to write a grid-pack to a file.
 * If you're using an std::ofstream, be sure to open it with ios::binary.
 */
bool WriteStreamPack(StreamPack& streamPack, std::ostream& out);

///	reads a grid-pack from a binary-stream
/**
 * This method is mainly intended to read a grid-pack from a file.
 * If you're using an std::ifstream, be sure to open it with ios::binary.
 */
bool ReadStreamPack(StreamPack& streamPack, std::istream& in);


////////////////////////////////////////////////////////////////////////
//	StreamPack
///	Container that stores BinaryStream objects with a tag.
/**
 * The StreamPack stores instances of BinaryStream together with a tag.
 * Using this tag you may access a stream.
 * If a stream is requested but doesn't exist yet, it will be created.
 * You may iterate through all existing (tag, stream) pairs.
 * Note that StreamPack::iterator::first holds the tag and
 * StreamPack::iterator::second holds the pointer to the associated stream.
 */
class StreamPack
{
	public:
		typedef std::map<int, BinaryStream*>	StreamMap;
		typedef StreamMap::iterator	iterator;

	public:
		StreamPack()	{}
		~StreamPack();
		
	///	delete all streams.
		void clear();

	///	returns the number of streams
		inline size_t num_streams()					{return m_streamMap.size();}
	///	returns the summed sizes of all streams
		size_t total_size();

	///	returns whether a stream has already been created.
		inline bool has_stream(int tag)				{return m_streamMap.find(tag) != m_streamMap.end();}

	///	returns the stream which is associated with tag.
	/**	If no stream exists, which is associated with tag,
		a new one will be created.*/
		BinaryStream* get_stream(int tag);

	///	resets read and write pointers of all streams
		void reset_streams();
		
	///	erases the stream associated with the given tag.
		void erase_stream(int tag);

	///	returns the begin-iterator to the stored (tag, stream)-pairs
	/**	use iter->first to access the tag and iter->second to access the stream.*/
		inline iterator begin()						{return m_streamMap.begin();}

	///	returns the end-iterator to the stored (tag, stream)-pairs
	/**	use iter->first to access the tag and iter->second to access the stream.*/
		inline iterator end()						{return m_streamMap.end();}

	protected:
		StreamMap	m_streamMap;
};

}

#endif
