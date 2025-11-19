/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Andreas Vogel
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

#ifndef __H__UG__COMMON__ERROR__
#define __H__UG__COMMON__ERROR__

#include <string>
#include <vector>
#include <sstream>

#include <exception>

/// \addtogroup ugbase_common
/// \{

void ug_throw_error();

#ifdef __GNUC__
#define PRETTY_FUNCTION __PRETTY_FUNCTION__
#else
#define PRETTY_FUNCTION __FUNCTION__
#endif

////////////////////////////////////////////////////////////////////////////////
// UG Throw / Catch
////////////////////////////////////////////////////////////////////////////////

#define UG_THROW(msg)		{ug_throw_error(); std::stringstream __ss; __ss << msg; \
							throw(ug::UGError(__ss.str(),__FILE__,__LINE__));}

/// UG_COND_THROW(cond, msg) : performs a UG_THROW(msg) if cond == true
#define UG_COND_THROW(cond, msg) { if(cond) { UG_THROW(msg); } }


#define UG_CATCH_THROW(msg)	catch(ug::UGError& err){std::stringstream __ss; __ss << msg;\
							  err.push_msg(__ss.str(),__FILE__,__LINE__); throw(err);} \
	catch(const std::exception& ex)	{	std::stringstream __ss; __ss << msg;\
	  	  	  	  	  	  	  	  throw ug::UGError(__ss.str(), ex,__FILE__,__LINE__); }

#define UG_CATCH_PRINT(msg) \
	catch (ug::UGError& err)\
	{\
		std::stringstream __ss;\
		__ss << msg;\
		err.push_msg(__ss.str(), __FILE__, __LINE__);\
		UG_LOG(err.get_stacktrace());\
	}\
	catch (const std::exception& ex)\
	{\
		std::stringstream __ss;\
		__ss << msg;\
		ug::UGError err(__ss.str(), ex, __FILE__, __LINE__);\
		UG_LOG(err.get_stacktrace());\
	}

#define UG_CATCH_THROW_FUNC()	UG_CATCH_THROW(PRETTY_FUNCTION << "failed. ")

// end group ugbase_common
/// \}

////////////////////////////////////////////////////////////////////////////////
// UG Error
////////////////////////////////////////////////////////////////////////////////

namespace ug{

std::string ErrorStringFromStdException(const std::exception *pex);

/// \addtogroup ugbase_common
/// \{

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

		UGError(const std::string &msg, const std::exception &ex, const char *file, const unsigned long line);

	///	virtual destructor
		virtual ~UGError()	{}

	///	adds a message to the message stack
		void push_msg(const std::string& msg, const char* file = " -- no file -- ",
		              const unsigned long line = 0)
		{
			m_vMsg.push_back(msg);
			m_vFile.emplace_back(file);
			m_vLine.push_back(line);
		}

	///	adds a message to the message stack
		void push_msg(const char* msg, const char* file = " -- no file -- ",
		              const unsigned long line = 0)
		{
			m_vMsg.emplace_back(msg);
			m_vFile.emplace_back(file);
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

		std::string get_stacktrace() const
		{
			std::stringstream ss;
			for(size_t i = 0; i < num_msg(); ++i)
			{
				ss << get_file(i) << ':' << get_line(i) << " : " << get_msg(i) << '\n';
			}
			return ss.str();
		}

	protected:
		std::vector<std::string> m_vMsg; //< Message stack
		std::vector<std::string> m_vFile; //< File stack
		std::vector<unsigned long> m_vLine; //< Line stack
};


////////////////////////////////////////////////////////////////////////////////
///	This special error is used to perform a soft-abort e.g. during script execution.
/**	Note that SoftAbort isn't an error per se. It is just a tool to use the error-handling
 * mechanism to exit through multiple layers of code-execution.*/
class SoftAbort : public UGError
{
	public:
		SoftAbort(std::string msg) : UGError(msg.c_str())	{}
};


////////////////////////////////////////////////////////////////////////////////
// some basic assertions
#define THROW_IF_NOT_EQUAL(s1, s2) { UG_COND_THROW(s1 != s2, "mismatch: " << UG_TO_STRING(s1) << " = " << s1 << "  !=  " << UG_TO_STRING(s2) << " = " << s2 << "."); }
#define THROW_IF_NOT_EQUAL_3(s1, s2, s3) { THROW_IF_NOT_EQUAL(s1, s2); THROW_IF_NOT_EQUAL(s1, s3); }
#define THROW_IF_NOT_EQUAL_4(s1, s2, s3, s4) { THROW_IF_NOT_EQUAL(s1, s2); THROW_IF_NOT_EQUAL(s1, s3); THROW_IF_NOT_EQUAL(s1, s4); }

#define ASSERT_EQUAL(s1, s2) { UG_COND_THROW(s1 != s2, "mismatch: " << UG_TO_STRING(s1) << " = " << s1 << "  !=  " << UG_TO_STRING(s2) << " = " << s2 << "."); }
#define ASSERT_EQUAL_3(s1, s2, s3) { ASSERT_EQUAL(s1, s2); ASSERT_EQUAL(s1, s3); }
#define ASSERT_EQUAL_4(s1, s2, s3, s4) { ASSERT_EQUAL(s1, s2); ASSERT_EQUAL(s1, s3); ASSERT_EQUAL(s1, s4); }


// end group ugbase_common
/// \}

} // end namespace ug

#endif