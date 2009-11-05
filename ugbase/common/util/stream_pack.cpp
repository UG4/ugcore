// created by Sebastian Reiter
// y09 m11 d04
// s.b.reiter@googlemail.com

#include "stream_pack.h"

using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	methods

//	writes a stream-pack to a binary-stream
bool WriteStreamPack(StreamPack& streamPack, std::ostream& out)
{
//	layout: endianess, magicNumber, version, numStreams, {tag, size, data}
	char endianess = 1;
	int magicNumber = 84197419;
	int version = 1;
	int numStreams = (int)streamPack.num_streams();

//	allows to check whether the file is valid.
	out.write((char*)&endianess, sizeof(char));
	out.write((char*)&magicNumber, sizeof(int));
	out.write((char*)&version, sizeof(int));

	out.write((char*)&numStreams, sizeof(int));

//	write the streams
	StreamPack::iterator iEnd = streamPack.end();
	for(StreamPack::iterator iter = streamPack.begin();
		iter != iEnd; ++iter)
	{
		int tag = (int)iter->first;
		BinaryStream* stream = iter->second;
		int streamSize = (int)stream->size();
		
		out.write((char*)&tag, sizeof(int));
		out.write((char*)&streamSize, sizeof(int));
		out.write((char*)stream->buffer(), streamSize);
	}

//	write the magic number again to allow a correctness check on read.
	out.write((char*)&magicNumber, sizeof(int));

	return true;
}

//	reads a stream-pack from a binary-stream
bool ReadStreamPack(StreamPack& streamPack, std::istream& in)
{
//	clear the pack
	streamPack.clear();

//	layout: endianess, magicNumber, numStreams, {tag, size, data}
	char endianess;
	int magicNumber;
	int version;
	int numStreams;

//	check endianess and magic number
	in.read((char*)&endianess, sizeof(char));
	in.read((char*)&magicNumber, sizeof(int));
	in.read((char*)&version, sizeof(int));

	if(endianess != 1)
		return false;
	if(magicNumber != 84197419)
		return false;
	if(version != 1)
		return false;

	in.read((char*)&numStreams, sizeof(int));

//	write the streams
	for(int i = 0; i < numStreams; ++i)
	{
	//	read the tag
		int tag;
		in.read((char*)&tag, sizeof(int));

	//	read the stream-size
		int streamSize;
		in.read((char*)&streamSize, sizeof(int));

	//	get the stream
		BinaryStream* stream = streamPack.get_stream(tag);
		stream->resize(streamSize);
		in.read((char*)stream->buffer(), streamSize);
	}

//	read the magic number again and check for correctness
	in.read((char*)&magicNumber, sizeof(int));
	if(magicNumber != 84197419)
		return false;

	return true;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	StreamPack - implementation
StreamPack::~StreamPack()
{
	clear();
}

void StreamPack::clear()
{
//	delete all streams
	iterator iEnd = end();
	for(iterator iter = begin(); iter != iEnd; ++iter)
		delete iter->second;
	m_streamMap.clear();
}

size_t StreamPack::total_size()
{
	iterator iEnd = end();
	size_t s = 0;
	for(iterator iter = begin(); iter != iEnd; ++iter)
		s += (size_t)iter->second->size();

	return s;
}

BinaryStream* StreamPack::get_stream(int tag)
{
//	if the entry has already been created, simply return the associated stream.
	iterator iter = m_streamMap.find(tag);
	if(iter != m_streamMap.end())
		return iter->second;

//	we have to create the stream first
	BinaryStream* stream = new BinaryStream;
	m_streamMap[tag] = stream;
	return stream;
}

}
