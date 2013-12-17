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
	catch(std::bad_alloc& ex)	{	std::stringstream __ss; __ss << msg;\
	  	  	  	  	  	  	  	  throw ug::UGError(__ss.str(), ex,__FILE__,__LINE__); } \
	catch(std::bad_cast& ex)	{	std::stringstream __ss; __ss << msg;\
	  	  	  	  	  	  	  	  throw ug::UGError(__ss.str(), ex,__FILE__,__LINE__); } \
	catch(std::exception& ex)	{	std::stringstream __ss; __ss << msg;\
	  	  	  	  	  	  	  	  throw ug::UGError(__ss.str(), ex,__FILE__,__LINE__); }

#define UG_CATCH_THROW_FUNC()	UG_CATCH_THROW(PRETTY_FUNCTION << "failed. ")

// end group ugbase_common
/// \}

////////////////////////////////////////////////////////////////////////////////
// UG Error
////////////////////////////////////////////////////////////////////////////////

namespace ug{

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

		UGError(const std::string &msg, std::bad_alloc &ex, const char *file, const unsigned long line);
		UGError(const std::string &msg, std::bad_cast &ex, const char *file, const unsigned long line);
		UGError(const std::string &msg, std::exception &ex, const char *file, const unsigned long line);

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
// some basic assertions
#define THROW_IF_NOT_EQUAL(s1, s2) { UG_COND_THROW(s1 != s2, "missmatch: " << UG_TO_STRING(s1) << " = " << s1 << "  !=  " << UG_TO_STRING(s2) << " = " << s2 << "."); }
#define THROW_IF_NOT_EQUAL_3(s1, s2, s3) { THROW_IF_NOT_EQUAL(s1, s2); THROW_IF_NOT_EQUAL(s1, s3); }
#define THROW_IF_NOT_EQUAL_4(s1, s2, s3, s4) { THROW_IF_NOT_EQUAL(s1, s2); THROW_IF_NOT_EQUAL(s1, s3); THROW_IF_NOT_EQUAL(s1, s4); }

#define ASSERT_EQUAL(s1, s2) { UG_COND_THROW(s1 != s2, "missmatch: " << UG_TO_STRING(s1) << " = " << s1 << "  !=  " << UG_TO_STRING(s2) << " = " << s2 << "."); }
#define ASSERT_EQUAL_3(s1, s2, s3) { ASSERT_EQUAL(s1, s2); ASSERT_EQUAL(s1, s3); }
#define ASSERT_EQUAL_4(s1, s2, s3, s4) { ASSERT_EQUAL(s1, s2); ASSERT_EQUAL(s1, s3); ASSERT_EQUAL(s1, s4); }


// end group ugbase_common
/// \}

} // end namespace ug

#endif /* __H__UG__COMMON__ERROR__ */
