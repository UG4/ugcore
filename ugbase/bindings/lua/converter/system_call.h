/*
 * system_call.h
 *
 *  Created on: 01.03.2013
 *      Author: mrupp
 */


#ifndef SYSTEM_CALL_H_
#define SYSTEM_CALL_H_

#include <string>

class SystemCall
{
private:
	bool m_bCalled;
	int m_returnCode;
	std::string m_output;

public:
	bool successfull() { return m_bCalled && m_returnCode == 0; }
	bool called() { return m_bCalled; }
	int returnCode() { return m_returnCode; }
	const std::string &output() { return m_output; }

	/**
	 * call a system command (is passed to /usr/bin/sh
	 * @param command
	 */
	SystemCall(const std::string &command);
};


#endif /* SYSTEM_CALL_H_ */
