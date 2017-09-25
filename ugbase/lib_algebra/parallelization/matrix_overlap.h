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
/**	The overlap is realized by first adding new entries to processes with H-Master
 * interfaces so that the full matrix can be stored with each H-Master. New entries
 * are appended at the end of the matrix.
 *
 * For the new entries a new H-Master-Overlap interface is created and for associated
 * entries connected to H-Slave nodes a new H-Slave-Overlap interface is created.
 *
 * For each H-Slave entry we make the matrix row partially consistent, so that all
 * existing connections which connect two slave entries have the same value as
 * the associated connection between respective masters.
 *
 * H-Masters contain the full matrix row once the method is done.
 *
 * Rows of H-Master-Overlap entries are DirichletRows (a_ii=1, a_ij=0 for i!=j).
 *
 * The methods replaces the AlgebraLayout in matInOut with a new one which
 * contains master_overlap and slave_overlap interfaces.
 *
 * \param matInOut	An additive parallel square matrix without overlap.
 */
template <class TMatrix>
void CreateOverlap (TMatrix& matInOut);

}//	end of namespace

#include "matrix_overlap_impl.h"

#endif	//__H__UG_matrix_overlap
