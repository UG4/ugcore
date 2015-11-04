/*
 * completion.h
 *
 *  Created on: 24.06.2013
 *      Author: mrupp
 */

#ifndef COMPLETION_H_
#define COMPLETION_H_
namespace ug{
namespace bridge{
int CompletionFunction(char *buf, int len, int buflen, int iPrintCompletionList);
void SetOtherCompletions(const char **otherCompletions, int nr);
}
}


#endif /* COMPLETION_H_ */
