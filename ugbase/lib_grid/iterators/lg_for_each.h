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
