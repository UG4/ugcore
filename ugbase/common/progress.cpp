/*
 * progress.cpp
 *
 *  Created on: 17.05.2013
 *      Author: mrupp
 */

#include "progress.h"
namespace ug{

int Progress::totalDepth = 0;
int Progress::lastUpdateDepth = -1;

Progress::Progress(int minSecondsUntilProgress)
{
	m_minSecondsUntilProgress = minSecondsUntilProgress;
	m_length=100;
	bStarted=false;
	myDepth = totalDepth++;
}

}
