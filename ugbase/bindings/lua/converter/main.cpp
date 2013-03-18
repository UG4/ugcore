/*
 * \file	main.cpp
 * \author	Martin Rupp
 *
 * test file for parser
 */

#include <iostream>
#include "pclass.h"

using namespace std;

int main()
{
	pclass parser;
	parser.parse(
	"function bla(x)\n"
	"local s = x*x*25"
	"return math.sin(s) + math.sin(s*x)\n"
	"end\n"	
	);
	const char *a=	"local I_single_star\n"
"local I_single_starw\n"
"if a < 5 and b > 4 then\n"
"	I_single_star = (1+1)*2+1 \n"
"	I_single_starw = 103\n"
"end\n"
"local I_single = 0.48 * I_single_star\n"
"local afff = 32*2+math.sin(32+3)+a\n"
"\n"
"local I_single_star = 10.3\n"
"local I_single = 0.48 * I_single_star\n"
"local H_p = 0.612 * 1e-3\n"
"local H_p_star = 0.34 * 1e-3\n"
"local f_stim = 1 - f_unstim\n"
	"end\n";
	
	
//	parser.reduce();
	if(parser.has_errors())
	{
		cout << parser.err.str() << "\n";
	}
	else
	{
		parser.createC(cout, 1);
		parser.print_globals(cout);
		parser.print_locals(cout);
		//parser.createLUA(cout);
	}
	
	return 0;
}

