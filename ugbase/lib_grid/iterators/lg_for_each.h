/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_lg_for_each
#define __H__UG_lg_for_each

#include "common/util/vec_for_each.h"	//include end_for

#define lg_for_each(_feType, _feVar, _feCon) \
			for(Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>();\
				_feI != _feCon.end<_feType>(); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_const(_feType, _feVar, _feCon) \
			for(Grid::traits<_feType>::const_iterator _feI = _feCon.begin<_feType>();\
				_feI != _feCon.end<_feType>(); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_template(_feType, _feVar, _feCon) \
			for(typename Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>();\
				_feI != _feCon.end<_feType>(); ++_feI){\
				_feType* _feVar = *_feI;


#define lg_for_each_in_lvl(_feType, _feVar, _feCon, _feLvl) \
			for(Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>(_feLvl);\
				_feI != _feCon.end<_feType>(_feLvl); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_in_lvl_template(_feType, _feVar, _feCon, _feLvl) \
			for(typename Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>(_feLvl);\
				_feI != _feCon.end<_feType>(_feLvl); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_in_subset(_feType, _feVar, _feCon, _feSubset) \
			for(Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>(_feSubset);\
				_feI != _feCon.end<_feType>(_feSubset); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_in_subset_template(_feType, _feVar, _feCon, _feSubset) \
			for(typename Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>(_feSubset);\
				_feI != _feCon.end<_feType>(_feSubset); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_in_subset_lvl(_feType, _feVar, _feCon, _feSubset, _feLvl) \
			for(Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>(_feSubset, _feLvl);\
				_feI != _feCon.end<_feType>(_feSubset, _feLvl); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_in_subset_lvl_template(_feType, _feVar, _feCon, _feSubset, _feLvl) \
			for(typename Grid::traits<_feType>::iterator _feI = _feCon.begin<_feType>(_feSubset, _feLvl);\
				_feI != _feCon.end<_feType>(_feSubset, _feLvl); ++_feI){\
				_feType* _feVar = *_feI;

#define lg_for_each_vertex_in_elem(_feVar, _feElem) \
			for(size_t _feI = 0; _feI < _feElem->num_vertices(); ++_feI){\
				Vertex* _feVar = _feElem->vertex(_feI);

#define lg_end_for	}

#endif	//__H__UG_lg_for_each
