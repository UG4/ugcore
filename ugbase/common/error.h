/*
 * error.h
 *
 *  Created on: 22.09.2011
 *      Author: sreiter, avogel
 */

#ifndef __H__UG__COMMON__ERROR__
#define __H__UG__COMMON__ERROR__

#include <string>
#include <vector>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
// UG Throw / Catch
////////////////////////////////////////////////////////////////////////////////

#define UG_THROW(msg)		{std::stringstream ss; ss << msg; \
							throw(ug::UGError(ss.str(),__FILE__,__LINE__));}

#define UG_CATCH_THROW(msg)	catch(ug::UGError& err){std::stringstream ss; ss << msg;\
							  err.push_msg(ss.str(),__FILE__,__LINE__); throw(err);}

////////////////////////////////////////////////////////////////////////////////
// UG Error
////////////////////////////////////////////////////////////////////////////////

namespace ug{

///	Instances of this class or of derived classes are thrown if errors arise.
/**	By default the error-code is 0 and terminate returns false.*/
class UGError
{
	public:
		UGError(const char* msg,
		        const char* file = " -- no file -- ", const unsigned long line = 0)
			{push_msg(msg, file, line);}
		UGError(const std::string& msg,
		        const char* file = " -- no file -- ", const unsigned long line = 0)
			{push_msg(msg, file, line);}

	///	virtual destructor
		virtual ~UGError()	{}

	///	adds a message to the message stack
		void push_msg(const std::string& msg, const char* file = " -- no file -- ",
		              const unsigned long line = 0)
		{
			m_vMsg.push_back(msg);
			m_vFile.push_back(file);
			m_vLine.push_back(line);
		}

	///	adds a message to the message stack
		void push_msg(const char* msg, const char* file = " -- no file -- ",
		              const unsigned long line = 0)
		{
			m_vMsg.push_back(msg);
			m_vFile.push_back(file);
			m_vLine.push_back(line);
		}

	///	returns the initial message
		const std::string& get_msg() const			{return m_vMsg.at(0);}

	///	number of messages in message-stack
		size_t num_msg() const 						{return m_vMsg.size();}

	///	returns a message in the message-stack (innermost is first)
		const std::string& get_msg(size_t i) const	{return m_vMsg.at(i);}

	/// returns the file where a message occured
		const std::string& get_file(size_t i) const	{return m_vFile.at(i);}

	///	returns the line where a message occured
		unsigned long get_line(size_t i) const{return m_vLine.at(i);}

	protected:
		std::vector<std::string> m_vMsg; //< Message stack
		std::vector<std::string> m_vFile; //< File stack
		std::vector<unsigned long> m_vLine; //< Line stack
};

} // end namespace ug

#endif /* __H__UG__COMMON__ERROR__ */
