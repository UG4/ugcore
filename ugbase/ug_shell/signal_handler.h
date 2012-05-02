/*
 * signal_handler.h
 *
 *  Created on: 02.05.2012
 *      Author: stephan grein
 */

#ifndef SIGNAL_HANDLER_H_
#define SIGNAL_HANDLER_H_

#include <stdexcept>
using std::runtime_error;

class SignalException : public runtime_error {
	public:
		SignalException(const std::string& err_msg) : std::runtime_error(err_msg) { }
};

class SignalHandler {
public:
    SignalHandler();
    ~SignalHandler();

    static const bool getExitSignal();
    static void setExitSignal();

    void        setupSignalHandlers();
    static void setSignalHandler(int ignored);

	protected:
    	static bool ExitSignal;
};

#endif /* SIGNAL_HANDLER_H_ */
