//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#ifndef __H__COMMON__COMMON__
#define __H__COMMON__COMMON__

#include <iostream>
#include <fstream>
#include <string>

/////////////////////////////////////////////////////////////////
// defines (currently here, should be compiling options)

// Warning and debug logs are disabled when compiling for release
#ifndef NDEBUG
	//#define UG_LOG_TO_FILE
	#define UG_ENABLE_WARNINGS
#else /* NDEBUG */
	//#undef UG_LOG_TO_FILE
	#undef UG_ENABLE_WARNINGS
#endif /* NDEBUG*/

/////////////////////////////////////////////////////////////////
// includes

#include "types.h"
//#include "contract.h"
#include "log.h"
#include "assert.h"
#include "static_assert.h"
#include "util/metaprogramming_util.h"

// depreciated, currently here for backward compatibility
#define LOG(msg) UG_LOG(msg)
#define STATIC_ASSERT(expr, msg) UG_STATIC_ASSERT(expr, msg)

////////////////////////////////////////////////////////////////////////////////////////////////
// save pointer (de-)allocation

#define SAFE_DELETE(a)		{if(a){ delete a; a = NULL;}}
#define SAFE_RELEASE(p)		{if(p) { (p)->Release(); (p)=NULL;}}

////////////////////////////////////////////////////////////////////////
///	Instances of this class or of derived classes are thrown if errors arise.
/**	By default the error-code is 0 and terminate returns false.*/
class UGError
{
	public:
		UGError(const char* msg) :
			m_msg(msg),
			m_terminate(false),
			m_code(0) {}

		UGError(int code, const char* msg) :
			m_msg(msg),
			m_terminate(false),
			m_code(code) {}

		UGError(bool terminate, const char* msg) :
			m_msg(msg),
			m_terminate(terminate),
			m_code(0) {}

		UGError(int code, bool terminate, const char* msg) :
			m_msg(msg),
			m_terminate(terminate),
			m_code(code) {}
		
		virtual ~UGError()	{}
		
		const std::string& get_msg()	{return m_msg;}
		bool terminate()				{return m_terminate;}
		int get_code()					{return m_code;}
		
	protected:
		std::string	m_msg;
		bool m_terminate;
		int m_code;
};

class UGFatalError : public UGError
{
	public:
		UGFatalError(const char* msg) :
			UGError(true, msg) {}
		UGFatalError(int code, const char* msg) :
			UGError(code, true, msg) {}
		UGFatalError(bool terminate, const char* msg) :
			UGError(terminate, msg) {}
		UGFatalError(int code, bool terminate, const char* msg) : 
			UGError(code, terminate, msg) {}
};


#endif
