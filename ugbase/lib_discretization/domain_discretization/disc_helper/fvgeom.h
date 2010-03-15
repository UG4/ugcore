/*
 * fvgeom.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef FVGEOM_H_
#define FVGEOM_H_

// extern libraries
#include <cmath>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lib_grid.h"
#include "lib_algebra/lib_algebra.h"

// library intern includes
#include "lib_discretization/reference_element/reference_elements.h"
#include "lib_discretization/local_shape_function_set/local_shape_function_set_factory.h"

namespace ug{

template <typename TElem>
class SD_Values{
	private:
		static const unsigned int _Dim = reference_element_traits<TElem>::dim;
		static const unsigned int _num_co = reference_element_traits<TElem>::num_corners;

	public:
		inline uint num_sh() const {return _num_co;}
		inline number* shapeptr() const {return _shape;}
		inline number shape(int i) const {return _shape[i];}
		inline const MathVector<_Dim>* gradptr_global() const {return _grad_global;}
		inline const MathVector<_Dim>& grad_global(int i) const {return _grad_global[i];}
		inline const MathMatrix<_Dim,_Dim>& J() const {return _J;}
		inline const MathMatrix<_Dim,_Dim>& Jinv() const {return _Jinv;}
		inline number detJ() const {return _detJ;}

	public:
		bool update_local(const MathVector<_Dim>& ip_local)
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			const LocalShapeFunctionSet<ref_elem_type>& TrialSpace = LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);

			for(unsigned int i = 0; i < _num_co; ++i)
			{
				if(TrialSpace.evaluate(i, ip_local, _shape[i]) == false) return false;
				if(TrialSpace.evaluate_grad(i, ip_local, _grad_local[i]) == false) return false;
			}

			return true;
		}

		bool update_global(const MathVector<_Dim>& ip_local, MathVector<_Dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type _refElem;

			// compute transformation inverse and determinate at ip
			if(_refElem.Trafo(corners, ip_local, _Jinv, _detJ) == false) return false;

			// compute _J
			Inverse(_J, _Jinv);

			// compute Shape Gradient
			for(unsigned int i = 0; i < _num_co; ++i)
			{
				MatVecMult(_grad_global[i], _Jinv, _grad_local[i]);
			}
			return true;
		}

	protected:
		number _shape[_num_co];
		MathVector<_Dim> _grad_local[_num_co];

		MathVector<_Dim> _grad_global[_num_co];
		MathMatrix<_Dim,_Dim> _J;
		MathMatrix<_Dim,_Dim> _Jinv;
		number _detJ;
};

// predeclaration
template <typename TElem>
class FVElementGeometry;

template <typename TElem>
class SubControlVolume{
	private:
		static const unsigned int _Dim = reference_element_traits<TElem>::dim;

	public:
		inline unsigned int local_corner_id() const {return _co;}
		inline const MathVector<_Dim>& global_corner_pos() const {return _center;}
		inline number volume() const {return _volume;}

		inline bool update_local(unsigned int corner, const FVElementGeometry<TElem>& geo)
		{
			// corner with id 'corner' is center of the entire control volume
			_co = corner;
			return true;
		}

		inline bool update_global(unsigned int corner, const FVElementGeometry<TElem>& geo, MathVector<_Dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type _refElem;

			// remember center coordinates
			_center = corners[corner];

			// compute size of scv
			if(_Dim == 1)
			{
				_volume = 0.5*(corners[0][0] - corners[1][0]);
			}
			else if(_Dim == 2)
			{
				_volume = qarea(corners[corner],
								geo.obj_midpoint_global(1, _refElem.id(0, corner, 1, 1)),
								geo.obj_midpoint_global(2, _refElem.id(0, corner, 2, 0)),
								geo.obj_midpoint_global(1, _refElem.id(0, corner, 1, 0)) );
			}
			else if(_Dim == 3)
			{
				assert(0 && "Must be implemented");
			}

			return true;
		}

	private:
		number qarea(const MathVector<2>& x0, const MathVector<2>& x1, const MathVector<2>& x2, const MathVector<2>& x3)
		{
			number tmp = (x3[1]-x1[1])*(x2[0]-x0[0])-(x3[0]-x1[0])*(x2[1]-x0[1]);
			return 0.5 * fabs( tmp );
		}

	protected:
		unsigned int _co;
		MathVector<_Dim> _center;
		number _volume;
};

template <typename TElem>
class SubControlVolumeFace{
	private:
		static const unsigned int _Dim = reference_element_traits<TElem>::dim;

	public:
		inline unsigned int from() const {return _from;}
		inline unsigned int to() const {return _to;}
		inline uint num_ip() const {return 1;}
		inline const MathVector<_Dim>& local_ip() const {return _ip_local;}
		inline const MathVector<_Dim>& global_ip() const {return _ip_global;}
		inline const MathVector<_Dim>& normal() const {return _normal;}
		inline const SD_Values<TElem>& sdv() const {return _sdv;}

	public:
		inline bool update_local(int edge, const FVElementGeometry<TElem>& geo)
		{
			static const typename reference_element_traits<TElem>::reference_element_type _refElem;

			// normal of this scvf points from corner '_from' to corner '_to'
			for(unsigned int i = 0; i < _refElem.num_obj_of_obj(1, edge, 0); ++i)
			{
				_from = _refElem.id(1,edge, 0, 0);
				_to = _refElem.id(1,edge, 0, 1);
			}

			// compute ip on scvf
			int scale = 0;
			_ip_local = 0.0;
			for(unsigned int d = _Dim; d >= 1; --d)
			{
				for(unsigned int i = 0; i < _refElem.num_obj_of_obj(1, edge, d); ++i)
				{
					_ip_local += geo.obj_midpoint_local(d, _refElem.id(1, edge, d, i));
					++scale;
				}
			}
			_ip_local *= 1./scale;

			_sdv.update_local(_ip_local);
			return true;
		}

		inline bool update_global(int edge, const FVElementGeometry<TElem>& geo, MathVector<_Dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type _refElem;

			// compute ip on scvf
			int scale = 0;
			_ip_global = 0.0;
			for(unsigned int d = _Dim; d >= 1; --d)
			{
				for(unsigned int i = 0; i < _refElem.num_obj_of_obj(1, edge, d); ++i)
				{
					_ip_global += geo.obj_midpoint_global(d, _refElem.id(1, edge, d, i));
					++scale;
				}
			}
			_ip_global *= 1./scale;

			//compute normal on scvf
			if(_Dim == 1)
			{
				_normal[0] = 1;
			}
			if(_Dim == 2)
			{
				MathVector<_Dim> diff = geo.obj_midpoint_global(2, _refElem.id(1, edge, 2, 0)); // center of element
				diff -= geo.obj_midpoint_global(1, _refElem.id(1, edge, 1, 0)); // edge midpoint

				_normal[0] = diff[1];
				_normal[1] = -diff[0];
			}
			else if(_Dim == 3)
			{
				assert(0 && "Must be implemented");
			}

			_sdv.update_global(_ip_local, corners);
			return true;
		}

	protected:
		unsigned int _from, _to;
		MathVector<_Dim> _ip_local;

		MathVector<_Dim> _ip_global;
		MathVector<_Dim> _normal;
		SD_Values<TElem> _sdv;
};


template <typename TElem>
class FVElementGeometry {
	private:
		static const unsigned int _Dim = reference_element_traits<TElem>::dim;
		static const unsigned int _num_co = reference_element_traits<TElem>::num_corners;
		static const unsigned int _num_ed = reference_element_traits<TElem>::num_edges;
		static const unsigned int _num_fa = reference_element_traits<TElem>::num_faces;
		static const unsigned int _num_vol = reference_element_traits<TElem>::num_volumes;

	public:
		FVElementGeometry()
		{
			static const typename reference_element_traits<TElem>::reference_element_type _refElem;

			// compute local corners
			for(unsigned int i = 0; i < _num_co; ++i)
			{
				_co_local[i] = _refElem.corner(i);
			}

			// compute local midpoints for all geometric objects with  0 < dim <= _Dim
			for(unsigned int d = 1; d <= _Dim; ++d)
			{
				// loop geometric objects of dimension d
				for(unsigned int i = 0; i < _refElem.num_obj(d); ++i)
				{
					_obj_mid_local[d-1][i] = 0.0;
					// add corner coordinates of the corners of the geometric object
					for(unsigned int j = 0; j < _refElem.num_obj_of_obj(d, i, 0); ++j)
					{
						_obj_mid_local[d-1][i] += _co_local[_refElem.id(d, i, 0, j)];
					}
					// scale for correct averaging
					_obj_mid_local[d-1][i] *= 1./(_refElem.num_obj_of_obj(d, i, 0));
				}
			}

			// compute local scvf
			for(unsigned int i = 0; i < num_scvf(); ++i)
			{
				_scvf[i].update_local(i, *this);
			}

			// compute local scv
			for(unsigned int i = 0; i < num_scv(); ++i)
			{
				_scv[i].update_local(i, *this);
			}

		}

		bool update(MathVector<_Dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type _refElem;

			// copy global corners
			for(unsigned int i = 0; i < _num_co; ++i)
			{
				_co_global[i] = corners[i];
			}

			// compute local midpoints for all geometric objects with  0 < dim <= _Dim
			for(unsigned int d = 1; d <= _Dim; ++d)
			{
				// loop geometric objects of dimension d
				for(unsigned int i = 0; i < _refElem.num_obj(d); ++i)
				{
					_obj_mid_global[d-1][i] = 0.0;
					// add corner coordinates of the corners of the geometric object
					for(unsigned int j = 0; j < _refElem.num_obj_of_obj(d, i, 0); ++j)
					{
						_obj_mid_global[d-1][i] += _co_global[_refElem.id(d, i, 0, j)];
					}
					// scale for correct averaging
					_obj_mid_global[d-1][i] *= 1./(_refElem.num_obj_of_obj(d, i, 0));
				}
			}

			// compute scvf
			for(unsigned int i = 0; i < num_scvf(); ++i)
			{
				_scvf[i].update_global(i, *this, corners);
			}

			// compute local scv
			for(unsigned int i = 0; i < num_scv(); ++i)
			{
				_scv[i].update_global(i, *this, corners);
			}

			return true;
		}

	public:
		unsigned int num_scvf() const
		{return _num_ed;};

		unsigned int num_scv() const
		{return _num_co;}

		const SubControlVolumeFace<TElem>& scvf(unsigned int i) const
		{
			return _scvf[i];
		}

		const SubControlVolume<TElem>& scv(unsigned int i) const
		{
			return _scv[i];
		}


		const MathVector<_Dim>& obj_midpoint_local(unsigned int dim_i, unsigned int i) const
		{
			assert(dim_i > 0 && dim_i <= _Dim);
			return _obj_mid_local[dim_i-1][i];
		}

		const MathVector<_Dim>& obj_midpoint_global(unsigned int dim_i, unsigned int i) const
		{
			assert(dim_i > 0 && dim_i <= _Dim);
			return _obj_mid_global[dim_i-1][i];
		}

	public:

		static const unsigned int _max_obj = _num_ed;

		MathVector<_Dim> _co_global[_num_co];
		MathVector<_Dim> _obj_mid_global[_Dim][_max_obj];

		MathVector<_Dim> _co_local[_num_co];
		MathVector<_Dim> _obj_mid_local[_Dim][_max_obj];

		SubControlVolume<TElem> _scv[_num_co];
		SubControlVolumeFace<TElem> _scvf[_num_ed];
};

}

#endif /* FVGEOM_H_ */
