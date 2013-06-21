/*
 * catch_std.h
 *
 *  Created on: 10.06.2013
 *      Author: mrupp
 */

#ifndef __H_COMMON_CATCH_STD_H_
#define __H_COMMON_CATCH_STD_H_
#include "log.h"
#include "assert.h"

#ifdef NDEBUG
#define UG_LOG_CATCH(expr)\
		UG_LOG_ALL_PROCS(	"\nA std::exception has been thrown:\n"\
				"Description:   " << expr << "\n"\
				"File:        " << __FILE__ << "\n"\
				"Line:        " << __LINE__ << "\n\n"); \
				assert(0);
#else
	#define UG_LOG_CATCH(expr) UG_ASSERT(false, "A std::exception has been thrown:\n" << expr)
#endif

#define CATCH_STD_EXCEPTIONS()\
	catch(std::bad_alloc& ex)	{	UG_LOG_CATCH("bad_alloc caught: " << ex.what() << "\nThis is mostly caused by an OUT OF MEMORY - condition.\n"\
			 	 << "You might have to reduce your problem size, go parallel or increase memory.")	} \
	catch(std::exception& ex)	{	UG_LOG_CATCH("bad_cast caught: " << ex.what() << "\nThis is caused by a Casting error in classes.")	} \



#endif /* __H_COMMON_CATCH_STD_H_ */
