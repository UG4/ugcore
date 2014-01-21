/*
 * parallel_archive.h
 *
 *  Created on: 13.12.2013
 *      Author: mrupp
 */

#ifndef PARALLEL_ARCHIVE_H_
#define PARALLEL_ARCHIVE_H_

#include "pcl_process_communicator.h"
#include "common/util/binary_stream.h"
#include "common/log.h"
#include <map>
#include <string>
#include <mpi.h>
#include "common/util/string_util.h"
#include "common/util/smart_pointer.h"

namespace pcl{

struct FileBufferDescriptor
{
	FileBufferDescriptor(std::string _name, const char *_buf, size_t _size) :
		name(_name), buf(_buf), size(_size)
	{	}

	FileBufferDescriptor(std::string _name, ug::BinaryBuffer &buf) :
			name(_name), buf(buf.buffer()), size(buf.write_pos())
	{	}

	FileBufferDescriptor(std::string _name, ug::BinaryStream &buf) :
				name(_name), buf((char*)buf.buffer()), size(buf.size())
	{	}

	FileBufferDescriptor(std::string _name, std::stringstream &buf) :
					name(_name), buf(buf.str().c_str()), size(buf.str().length())
	{	}

	std::string name;
	const char *buf;
	size_t size;

};

void WriteParallelArchive(pcl::ProcessCommunicator &pc, std::string strFilename, const std::vector<FileBufferDescriptor> &files);

template<typename TBuffer>
void WriteParallelArchive(pcl::ProcessCommunicator &pc, std::string strFilename, const std::map<std::string, TBuffer> &files)
{
	std::vector<FileBufferDescriptor> fdesc;
	for(typename std::map<std::string, TBuffer>::const_iterator it = files.begin(); it != files.end(); ++it)
		 files.push_back(FileBufferDescriptor(ug::FilenameWithoutPath(it->first), it->second));
	WriteParallelArchive(pc, strFilename, fdesc);
}

/**
 * This class creates one .a archive out of several parallel file writes
 * This has two advantages
 * 1. Instead of writing 1024 files for 1024 cores, this writes 1 file
 * 2. It uses MPI I/O, which will be a LOT faster when number of cores are high
 *
 * For extracting the archive use
 *   ar x thefile.a
 *
 * Example usage:
 * \code
 * 	ParallelArchive pa("test.tar", pcl::ProcessCommunicator(pcl::PCD_WORLD));
	std::stringstream &b = pa.create_stringstream_file(std::string("bla") + ToString(pcl::GetProcRank()));
	b << "Hello World!";
	pa.write();
 *	\endcode
 */
class ParallelArchive
{
private:
	/// internal virtual buffer interface to support different buffers
	struct BufferInterface
	{
		virtual ~BufferInterface() {}
		virtual const char *buffer()=0;
		virtual size_t size()=0;
		virtual void update() {}
	};

	struct BufferBinaryBuffer : public BufferInterface
	{
		ug::BinaryBuffer internal_buffer;
		virtual const char *buffer() { return internal_buffer.buffer(); };
		virtual size_t size() { return internal_buffer.write_pos(); }
	};

	struct BufferBinaryStream : public BufferInterface
	{
		ug::BinaryStream internal_buffer;
		virtual const char *buffer() { return (const char*)internal_buffer.buffer(); };
		virtual size_t size() { return internal_buffer.size(); }
	};

	struct Buffer_stringstream : public BufferInterface
	{
		std::stringstream internal_buffer;
		std::string m_s;
		virtual void update()
		{
			m_s = internal_buffer.str();
		}
		virtual const char *buffer() { return m_s.c_str(); }
		virtual size_t size() { return m_s.length(); }
	};

	struct ConstCharBuffer : public BufferInterface
	{
		ConstCharBuffer(const char *p, size_t s) : m_p(p), m_size(s) {}
		virtual const char *buffer() { return m_p; };
		virtual size_t size() { return m_size; }
		const char *m_p;
		size_t m_size;
	};


public:
	/**
	 * Create a parallel archive.
	 * Note that the file is written on descruction or when you call write() explicitely
	 * @param filename  the name of the archive. add .a for clearness
	 * @param pc		the process communicator used for MPI purposes. default WORLD
	 */
	ParallelArchive(std::string filename, pcl::ProcessCommunicator pc = pcl::ProcessCommunicator(pcl::PCD_WORLD))
		: m_filename(filename), m_pc(pc)
	{
		m_bWritten=false;
		m_bUnsafe = false;
	}

	///
	/**
	 * if the file(s) were not written already, this will write the file
	 * NOTE: Don't do this when using add_raw (unsafe!)
	 * NOTE: Communication happens here.
	 */
	~ParallelArchive()
	{
		if(!m_bWritten)
		{
			UG_COND_THROW(m_bUnsafe, "use ParallelArchive::write when using add_raw.")
			write();
		}
	}

	/**
	 * create a file inside the archive.
	 * @param name the filename inside the archive. can be a path, but only filename is taken
	 * @return a BinaryBuffer to write data to
	 */
	ug::BinaryBuffer &create_BinaryBuffer_file(std::string name)
	{
		SmartPtr<BufferBinaryBuffer> p = new BufferBinaryBuffer;
		files[name] = p;
		return p->internal_buffer;
	}

	/**
	 * create a file inside the archive.
	 * @param name the filename inside the archive. can be a path, but only filename is taken
	 * @return a BinaryStream to write data to
	 */
	ug::BinaryStream &create_BinaryStream_file(std::string name)
	{
		SmartPtr<BufferBinaryStream > p = new BufferBinaryStream;
		files[name] = p;
		return p->internal_buffer;
	}

	/**
	 * create a file inside the archive.
	 * @param name the filename inside the archive. can be a path, but only filename is taken
	 * @return a stringstream to write data to
	 */
	std::stringstream &create_stringstream_file(std::string name)
	{
		SmartPtr<Buffer_stringstream > p = new Buffer_stringstream;
		files[name] = p;
		return p->internal_buffer;
	}

	/**
	 * add raw buffer descriptors (see FileBufferDescriptor)
	 * NOTE: be sure all data/pointers are valid until ::write is called
	 * NOTE: You HAVE to use ParallelArchive::write explicitely
	 * using deconstructors is UNSAFE since data can be deconstructed before ParallelArchive.
	 * @param f
	 */
	void add_raw(FileBufferDescriptor f)
	{
		files[f.name] = new ConstCharBuffer(f.buf, f.size);
		m_bUnsafe = true;
	}

	/**
	 * explicitely writes the data
	 * NOTE: Communication happens here.
	 */
	void write()
	{
		std::vector<FileBufferDescriptor> fdesc;
		for(map_iterator it = files.begin(); it != files.end(); ++it)
		{
			(it->second)->update();
			fdesc.push_back(
					FileBufferDescriptor(ug::FilenameWithoutPath(it->first),
							(it->second)->buffer() , (it->second)->size() ) );
		}
		WriteParallelArchive(m_pc, m_filename, fdesc);
		files.clear();
		m_bWritten = true;
	}

	void create_new_archive(std::string filename)
	{
		if(m_bWritten == false)
			write();
		m_filename = filename;
	}

private:

	typedef std::map<std::string, SmartPtr<BufferInterface> >::iterator map_iterator;

	std::map<std::string, SmartPtr<BufferInterface> > files;
	std::string m_filename;
	pcl::ProcessCommunicator m_pc;
	bool m_bWritten;
	bool m_bUnsafe;
};

}
#endif /* PARALLEL_ARCHIVE_H_ */
