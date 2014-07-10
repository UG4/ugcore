/*
 * system_call.cpp
 *
 *  Created on: 01.03.2013
 *      Author: mrupp
 */

#include "system_call.h"

#ifndef WIN32
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#else
#include <tchar.h>
#include <windows.h>
#endif

using namespace std;

#if WIN32
SystemCall::SystemCall(const std::string &command);
{
	m_bCalled = false;
	m_returnCode = 0xDEAD;
	m_output = string("Could not call ") + command;

	STARTUPINFO si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);

	PROCESS_INFORMATION pi;
	ZeroMemory(&pi, sizeof(pi));

	//creating new process
	BOOL bResult = CreateProcess(NULL, command.c_str(), NULL, NULL, FALSE, CREATE_NO_WINDOW,
			NULL, NULL, &si, &pi);

	if(bResult)
	{
		//waiting for process termination
		WaitForSingleObject(pi.hProcess, INFINITE);

		m_bCalled = true;
		DWORD exitCode;
		if(GetExitCodeProcess(pi.hProcess, &exitCode))
		{
			m_returnCode = exitCode;
			m_output = string("called ") + command +
					"(todo: get output ( http://msdn.microsoft.com/en-us/library/ms682499(VS.85).aspx )";
		}
		else
			m_returnCode = -1;
	}
	// todo: get output ( http://msdn.microsoft.com/en-us/library/ms682499(VS.85).aspx )
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
}
#else
SystemCall::SystemCall(const std::string &command)
{
	m_bCalled = false;
	m_returnCode = 0xDEAD;
	m_output = string("Could not call ") + command;

	stringstream ss;
	FILE *fp = popen(command.c_str(), "r");
	if(fp == NULL) return;
	m_bCalled = true;
	char buf[255];
	while (fgets(buf, sizeof(buf)-1, fp) != NULL)
	ss << buf;

	m_output = ss.str();
	m_returnCode = pclose(fp);
}
#endif
