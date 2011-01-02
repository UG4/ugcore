/*
 * log.cpp
 *
 *  Created on: 10.03.2010
 *      Author: andreasvogel, sebastianreiter
 */

#include "common/log.h"
using namespace std;

namespace ug{

LogAssistant::
LogAssistant() :
	m_terminalOutputEnabled(true),
	m_fileOutputEnabled(false)
{
	set_debug_levels(-1);
	m_emptyBuf = &m_emptyBufInst;
	m_splitBuf = &m_splitBufInst;
	m_logBuf = clog.rdbuf();
	m_fileBuf = m_fileStream.rdbuf();

	m_splitBufInst.set_buffers(m_logBuf, m_fileBuf);

	update_ostream();
}

LogAssistant& LogAssistant::
instance()
{
	static LogAssistant log;
	return log;
}

bool LogAssistant::
enable_file_output(bool bEnable, const char* filename)
{
	if(bEnable){
		if(!m_fileStream.is_open()){
			m_fileStream.open(filename);
			if(!m_fileStream)
			{
				m_fileOutputEnabled = false;
				return false;
			}
		}
	}

	m_fileOutputEnabled = bEnable;
	update_ostream();
	return true;
}

bool LogAssistant::
enable_terminal_output(bool bEnable)
{
	m_terminalOutputEnabled = bEnable;
	update_ostream();
	return true;
}

void LogAssistant::
update_ostream()
{
	if(m_terminalOutputEnabled){
		if(m_fileOutputEnabled)
			clog.rdbuf(m_splitBuf);
		else
			clog.rdbuf(m_logBuf);
	}
	else{
		if(m_fileOutputEnabled){
			clog.rdbuf(m_fileBuf);
		}
		else{
			clog.rdbuf(m_emptyBuf);
		}
	}
}

bool
LogAssistant::
set_debug_levels(int lev)
{
	for(int i = 0; i < NUM_TAGS; ++i)
	{
		m_TagLevel[i] = lev;
	}
	return true;
}

bool
LogAssistant::
set_debug_level(Tags tags, int lev)
{
	m_TagLevel[tags] = lev;
	return true;
}
/*
std::ostream&
LogAssistant::
file_logger()
{
	static std::ofstream out("uglog.log");
	return out;
}
*/
}
