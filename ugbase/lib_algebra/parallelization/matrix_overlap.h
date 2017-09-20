/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_matrix_overlap
#define __H__UG_matrix_overlap

namespace ug{

/* NOTE TO DEVELOPERS:	This method is newer than the one in 'parallel_matrix_overlap_impl.hpp'
 * The main difference is that this method creates separate 'overlap' interfaces,
 * whereas the old implementation tries to create an overlap with only h-master
 * and h-slave interfaces (to my knowledge).
 */

///	Adds coefficients and connections to create an overlapping matrix
/**	The overlap is realized by creating new interfaces: For each H-Master interface
 * a new H-Master-Overlap interface is created and for each H-Slave interface a
 * new H-Slave-Overlap interface is created. Those interfaces are written to
 * masterOverlapOut and slaveOverlapOut.
 *
 * For each H-Slave-DOF we also set a Dirichlet-Row
 *
 * \param matInOut	An additive parallel square matrix without overlap.
 *
 * \returns	A parralel 'consistent' matrix (x=Ay with x,y consistent). H-Slave
 * 			and H-Master-Overlap rows are Dirichlet Rows (a_ii=1, a_ij=0 for i!=j).
 */
template <class TMatrix>
void CreateOverlap (TMatrix& matInOut,
                    IndexLayout& masterOverlapOut,
                    IndexLayout& slaveOverlapOut);

}//	end of namespace

#include "matrix_overlap_impl.h"

#endif	//__H__UG_matrix_overlap
