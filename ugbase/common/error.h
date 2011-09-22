/*
 * error.h
 *
 *  Created on: 22.09.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__COMMON__ERROR__
#define __H__UG__COMMON__ERROR__

#include <string>
#include <vector>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////
// UG Throw / Catch
////////////////////////////////////////////////////////////////////////////////

#define UG_THROW(msg)		{std::stringstream ss; ss << " UG Error in '"\
							<< __FILE__ << "' at line " << __LINE__ << ": "\
							<< msg; throw(ug::UGError(ss.str()));}
#define UG_THROW_FATAL(msg)	{std::stringstream ss; ss << " UG Error in '"\
							<< __FILE__ << "' at line " << __LINE__ << ": "\
							<< msg; throw(UGFatalError(ss.str()));}
#define UG_TRY 				try{
#define UG_THROW_ADD(msg)	}catch(ug::UGError ex) {std::stringstream ss; ss << msg;\
											  ex.push_msg(ss.str()); throw(ex);}

////////////////////////////////////////////////////////////////////////////////
// UG Error
////////////////////////////////////////////////////////////////////////////////

namespace ug{

///	Instances of this class or of derived classes are thrown if errors arise.
/**	By default the error-code is 0 and terminate returns false.*/
class UGError
{
	public:
		UGError(const char* msg) : m_bTerminate(false) {push_msg(msg);}
		UGError(const std::string& msg) : m_bTerminate(false) {push_msg(msg);}

		UGError(bool bExit, const char* msg) : m_bTerminate(bExit) {push_msg(msg);}
		UGError(bool bExit, const std::string& msg) : m_bTerminate(bExit) {push_msg(msg);}

	///	virtual destructor
		virtual ~UGError()	{}

	///	adds a message to the message stack
		void push_msg(const std::string& msg) 		{m_vMsg.push_back(msg);}

	///	adds a message to the message stack
		void push_msg(const char* msg) 		 		{m_vMsg.push_back(msg);}

	///	number of messages in message-stack
		size_t num_msg() const 						{return m_vMsg.size();}

	///	returns a message in the message-stack (innermost is first)
		const std::string& get_msg(size_t i) const	{return m_vMsg.at(i);}

	///	returns the initial message
		const std::string& get_msg() const			{return m_vMsg.at(0);}

	///	returns if program should terminate
		bool terminate() const						{return m_bTerminate;}

	protected:
		std::vector<std::string> m_vMsg; //< Message stack
		bool m_bTerminate;				 //< terminate flag
};

class UGFatalError : public UGError
{
	public:
		UGFatalError(const char* msg) 			  : UGError(true, msg) {}
		UGFatalError(bool bExit, const char* msg) :	UGError(bExit, msg) {}
};

} // end namespace ug

#endif /* __H__UG__COMMON__ERROR__ */
