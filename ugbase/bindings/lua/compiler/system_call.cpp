/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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

#ifdef WIN32
SystemCall::SystemCall(const std::string &command)
{
	m_bCalled = false;
	m_returnCode = 0xDEAD;
	m_output = string("Could not call ") + command;

	STARTUPINFO si;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);

	PROCESS_INFORMATION pi;
	ZeroMemory(&pi, sizeof(pi));

//	CreateProcess takes a char* string, since it may modify its contents.
//	we thus create a non-const copy of command
	const size_t cmdLen = command.length();
	char* cmd = new char[cmdLen + 1];
	memcpy(cmd, command.c_str(), cmdLen);
	cmd[cmdLen] = 0;

	//creating new process
	BOOL bResult = CreateProcess(NULL, cmd, NULL, NULL, FALSE, CREATE_NO_WINDOW,
			NULL, NULL, &si, &pi);

	delete[] cmd;

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
