/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
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

