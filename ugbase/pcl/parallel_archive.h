/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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

#ifndef PARALLEL_ARCHIVE_H_
#define PARALLEL_ARCHIVE_H_


#include <map>
#include <string>
#include <utility>
#include <mpi.h>

#include "common/log.h"
#include "common/util/string_util.h"
#include "common/util/smart_pointer.h"
#include "common/util/binary_stream.h"

#include "pcl_process_communicator.h"



namespace pcl{

struct FileBufferDescriptor
{
	FileBufferDescriptor(std::string _name, const char *_buf, size_t _size) :
		name(std::move(_name)), buf(_buf), size(_size)
	{	}

	FileBufferDescriptor(std::string _name, ug::BinaryBuffer& _buf) :
			name(std::move(_name)), buf(_buf.buffer()), size(_buf.write_pos())
	{	}

	FileBufferDescriptor(std::string _name, ug::BinaryStream& _buf) :
				name(std::move(_name)), buf(static_cast<char *>(_buf.buffer())), size(_buf.size())
	{	}

	// Initializing buf with buf.str().c_str() is unsafe, as _buf.str() is only temporary.
	// Therefore, commenting out (not used anywhere anyhow).
	//FileBufferDescriptor(std::string _name, std::stringstream& _buf) :
	//				name(_name), buf(_buf.str().c_str()), size(_buf.str().length())
	//{	}

	std::string name;
	const char *buf;
	size_t size;

};

void WriteParallelArchive(const ProcessCommunicator &pc, const std::string &strFilename, const std::vector<FileBufferDescriptor> &files);

template<typename TBuffer>
void WriteParallelArchive(ProcessCommunicator &pc, std::string strFilename, const std::map<std::string, TBuffer> &files)
{
	std::vector<FileBufferDescriptor> fdesc;
	for(typename std::map<std::string, TBuffer>::const_iterator it = files.begin(); it != files.end(); ++it)
		 files.push_back(FileBufferDescriptor(ug::FilenameWithoutPath(it->first), it->second));
	WriteParallelArchive(pc, std::move(strFilename), fdesc);
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
	std::stringstream &b = pa.create_stringstream_file(std::string("bla") + ToString(pcl::ProcRank()));
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
		virtual ~BufferInterface() = default;
		virtual const char *buffer()=0;
		virtual size_t size()=0;
		virtual void update() {}
	};

	struct BufferBinaryBuffer : BufferInterface
	{
		ug::BinaryBuffer internal_buffer;
		const char *buffer() override { return internal_buffer.buffer(); };
		size_t size() override { return internal_buffer.write_pos(); }
	};

	struct BufferBinaryStream : BufferInterface
	{
		ug::BinaryStream internal_buffer;
		const char *buffer() override { return static_cast<const char *>(internal_buffer.buffer()); };
		size_t size() override { return internal_buffer.size(); }
	};

	struct Buffer_stringstream : BufferInterface
	{
		std::stringstream internal_buffer;
		std::string m_s;
		void update() override { m_s = internal_buffer.str(); }
		const char *buffer() override { return m_s.c_str(); }
		size_t size() override { return m_s.length(); }
	};

	struct ConstCharBuffer : BufferInterface
	{
		ConstCharBuffer(const char *p, size_t s) : m_p(p), m_size(s) {}
		const char *buffer() override { return m_p; };
		size_t size() override { return m_size; }
		const char *m_p;
		size_t m_size;
	};


public:
	/**
	 * Create a parallel archive.
	 * Note that the file is written on destruction or when you call write() explicitly
	 * @param filename  the name of the archive. add .a for clearness
	 * @param pc		the process communicator used for MPI purposes. default WORLD
	 */
	explicit ParallelArchive(std::string filename, ProcessCommunicator pc = ProcessCommunicator(PCD_WORLD))
		: m_filename(std::move(filename)), m_pc(std::move(pc))
	{
		m_bWritten = false;
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
			if (m_bUnsafe)
			{
				UG_LOGN("ParallelArchive: Cannot write from destructor when using add_raw.\n"
					"Use the write() method when using add_raw.");
				return;
			}

			write();
		}
	}

	/**
	 * create a file inside the archive.
	 * @param name the filename inside the archive. can be a path, but only filename is taken
	 * @return a BinaryBuffer to write data to
	 */
	ug::BinaryBuffer &create_BinaryBuffer_file(const std::string &name)
	{
		SmartPtr p(new BufferBinaryBuffer);
		files[name] = p;
		return p->internal_buffer;
	}

	/**
	 * create a file inside the archive.
	 * @param name the filename inside the archive. can be a path, but only filename is taken
	 * @return a BinaryStream to write data to
	 */
	ug::BinaryStream &create_BinaryStream_file(const std::string &name)
	{
		SmartPtr p(new BufferBinaryStream);
		files[name] = p;
		return p->internal_buffer;
	}

	/**
	 * create a file inside the archive.
	 * @param name the filename inside the archive. can be a path, but only filename is taken
	 * @return a stringstream to write data to
	 */
	std::stringstream &create_stringstream_file(const std::string &name)
	{
		SmartPtr p(new Buffer_stringstream);
		files[name] = p;
		return p->internal_buffer;
	}

	/**
	 * add raw buffer descriptors (see FileBufferDescriptor)
	 * NOTE: be sure all data/pointers are valid until ::write is called
	 * NOTE: You HAVE to use ParallelArchive::write explicitly
	 * using deconstructors is UNSAFE since data can be deconstructed before ParallelArchive.
	 * @param f
	 */
	void add_raw(const FileBufferDescriptor &f)
	{
		files[f.name] = make_sp(new ConstCharBuffer(f.buf, f.size));
		m_bUnsafe = true;
	}

	/**
	 * explicitly writes the data
	 * NOTE: Communication happens here.
	 */
	void write()
	{
		std::vector<FileBufferDescriptor> fdesc;
		for(auto it = files.begin(); it != files.end(); ++it)
		{
			(it->second)->update();
			fdesc.emplace_back(ug::FilenameWithoutPath(it->first),
							(it->second)->buffer() , (it->second)->size() );
		}
		WriteParallelArchive(m_pc, m_filename, fdesc);
		files.clear();
		m_bWritten = true;
	}

	void create_new_archive(std::string filename)
	{
		if(m_bWritten == false)
			write();
		m_filename = std::move(filename);
	}

private:
	using map_iterator = std::map<std::string, SmartPtr<BufferInterface> >::iterator;

	std::map<std::string, SmartPtr<BufferInterface> > files;
	std::string m_filename;
	ProcessCommunicator m_pc;
	bool m_bWritten;
	bool m_bUnsafe;
};

}
#endif