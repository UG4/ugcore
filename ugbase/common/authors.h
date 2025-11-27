/*
 * Copyright (c) 2011:  Steinbeis Forschungszentrum (STZ Ölbronn)
 * Author: Michael Hoffer
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

#ifndef AUTHORS_H
#define	AUTHORS_H

#include <string>

namespace ug {

/// \addtogroup ugbase_common
/// \{

    /**
     * Names of the UG Core developer team, seperated by ','.
     */
    static std::string UG_AUTHORS
    	="Sebastian Reiter, Martin Rupp, Andreas Vogel, et al.";

    /**
     * License for binary usage.
     */
    static std::string UG_BINARY_LICENSE =
"ug4 Binary Code License\n"
"\n"
"ug4 (as defined below) is licensed to the licensee only upon \n"
"the condition that the licensee accepts all of the terms \n"
"contained in this Binary Code License Agreement (\"Agreement\"). \n"
"By downloading, installing and/or using ug4, the licensee \n"
"accepts the full terms of this Agreement.\n"
"\n"
"\n"
"Definitions\n"
"1. The binary distribution of the ug4 simulation framework (\"ug4\") \n"
"consists of the following binary components: The ugshell \n"
"application, the VRL-UG4 plugin and the ug4-plugin \n"
"ConvectionDiffusion.\n"
"2. Licensee means the individual that uses ug4 under the \n"
"terms of this License Agreement.\n"
"3. The licensors are the Goethe Center for Scientific Computing (G-CSC),\n"
" Goethe University Frankfurt am Main, Germany and the Steinbeis \n"
"Forschungszentrum, Ölbronn, Germany.\n"
"\n"
"License To Use\n"
"1. Subject to the terms and conditions of this Agreement, the \n"
"licensee is hereby granted a non-exclusive, non-transferable, \n"
"limited license to use the binary distribution of ug4 complete \n"
"and unmodified for the sole purpose of running programs and for \n"
"non-commercial use only.\n"
"2. Usage and/or publication of results obtained with the aid \n"
"of ug4 is only permitted if the following citation is included: \n"
"\"A. Vogel, S. Reiter, M. Rupp, A. Nägel, G. Wittum: UG4 - A Novel \n"
"Flexible Software System for Simulating PDE Based Models on High \n"
"Performance Computers, Comp. Vis. Sci.\"\n"
"3. Redistribution of ug4 is not permitted by this Agreement. \n"
"Additional licenses for developers and/or publishers are not \n"
"covered by this Agreement.\n"
"\n"
"Warranty\n"
"The program is licensed free of charge and there is no warranty \n"
"for the program, to the extent permitted by applicable law. The \n"
"program is provided \"as is\" without warranty of any kind, \n"
"either expressed or implied, including, but not limited to, the \n"
"implied warranties of merchantability and fitness for a particular \n"
"purpose. The entire risk as to the quality and performance of \n"
"the program is with the licensee. Should the program prove \n"
"defective, the licensee assumes the cost of all necessary \n"
"servicing, repair or correction.\n"
"\n"
"In no event unless required by applicable law will the licensor \n"
"be liable to the licensee for damages, including any general, \n"
"special, incidental or consequential damages arising out of the \n"
"use or inability to use the program (including but not limited \n"
"to loss of data or data being rendered inaccurate or losses \n"
"sustained by the licensee or a failure of the program to operate \n"
"with any other programs).\n";

// end group ugbase_common
/// \}

}

#endif