/*
 * fvgeom.h
 *
 *  Created on: 08.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__DISC_HELPER__FVGEOM__
#define __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__DISC_HELPER__FVGEOM__

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

template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::dim>
class SD_Values{
	public:
		static const int dim = reference_element_traits<TElem>::dim;
		static const int world_dim = TWorldDim;
		static const size_t m_num_co = reference_element_traits<TElem>::num_corners;

	public:
		inline size_t num_sh() const {return m_num_co;}
		inline number* shapeptr() const {return m_shape;}
		inline number shape(int i) const {return m_shape[i];}
		inline const MathVector<world_dim>* gradptr_global() const {return m_grad_global;}
		inline const MathVector<world_dim>& grad_global(int i) const {return m_grad_global[i];}
		inline const MathMatrix<world_dim,dim>& J() const {return m_JT;}
		inline const MathMatrix<dim,world_dim>& Jinv() const {return m_JTInv;}
		inline number detJ() const {return m_detJ;}

	public:
		bool update_local(const MathVector<dim>& ip_local)
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			const LocalShapeFunctionSet<ref_elem_type>& TrialSpace = LocalShapeFunctionSetFactory::inst().get_local_shape_function_set<ref_elem_type>(LSFS_LAGRANGEP1);

			for(size_t i = 0; i < m_num_co; ++i)
			{
				if(TrialSpace.evaluate(i, ip_local, m_shape[i]) == false) return false;
				if(TrialSpace.evaluate_grad(i, ip_local, m_grad_local[i]) == false) return false;
			}

			return true;
		}

		bool update_global(const MathVector<dim>& ip_local, MathVector<world_dim> corners[])
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			static ReferenceMapping<ref_elem_type, world_dim> mapping;

			mapping.update(corners);

			// compute transformation inverse and determinate at ip
			if(mapping.jacobian_transposed_inverse(ip_local, m_JTInv) == false) return false;
			if(mapping.jacobian_det(ip_local, m_detJ) == false) return false;

			// compute m_J
			// TODO: Maybe someone needs this. Compute J
			//Inverse(m_J, m_Jinv);

			// compute Shape Gradient
			for(size_t i = 0; i < m_num_co; ++i)
			{
				MatVecMult(m_grad_global[i], m_JTInv, m_grad_local[i]);
			}
			return true;
		}

	protected:
		number m_shape[m_num_co];
		MathVector<dim> m_grad_local[m_num_co];

		MathVector<world_dim> m_grad_global[m_num_co];
		MathMatrix<dim,world_dim> m_JT;
		MathMatrix<world_dim,dim> m_JTInv;
		number m_detJ;
};

// predeclaration
template <typename TElem, int TWorldDim = reference_element_traits<TElem>::dim>
class FVElementGeometry;

template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::dim>
class SubControlVolume{
	public:
		static const int dim = reference_element_traits<TElem>::dim;
		static const int world_dim = TWorldDim;

	public:
		inline size_t local_corner_id() const {return m_co;}
		inline const MathVector<world_dim>& global_corner_pos() const {return m_center;}
		inline number volume() const {return m_volume;}

		inline bool update_local(size_t corner, const FVElementGeometry<TElem, TWorldDim>& geo)
		{
			// corner with id 'corner' is center of the entire control volume
			m_co = corner;
			return true;
		}

		inline bool update_global(const FVElementGeometry<TElem, TWorldDim>& geo, MathVector<world_dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			// remember center coordinates
			m_center = corners[m_co];

			// compute size of scv
			if(world_dim == 1)
			{
				m_volume = 0.5*(corners[0][0] - corners[1][0]);
			}
			else if(world_dim == 2)
			{
				m_volume = qarea(corners[m_co],
								geo.obj_midpoint_global(1, m_refElem.id(0, m_co, 1, 1)),
								geo.obj_midpoint_global(2, m_refElem.id(0, m_co, 2, 0)),
								geo.obj_midpoint_global(1, m_refElem.id(0, m_co, 1, 0)) );
			}
			else if(world_dim == 3)
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
		number qarea(const MathVector<3>& x0, const MathVector<3>& x1, const MathVector<3>& x2, const MathVector<3>& x3){return -1;}

	protected:
		size_t m_co;
		MathVector<world_dim> m_center;
		number m_volume;
};

template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::dim>
class SubControlVolumeFace{
	public:
		static const int dim = reference_element_traits<TElem>::dim;
		static const int world_dim = TWorldDim;

	public:
		inline size_t from() const {return m_from;}
		inline size_t to() const {return m_to;}
		inline size_t num_ip() const {return 1;}
		inline const MathVector<dim>& local_ip() const {return m_ip_local;}
		inline const MathVector<world_dim>& global_ip() const {return m_ip_global;}
		inline const MathVector<world_dim>& normal() const {return m_normal;} // includes area
		inline const SD_Values<TElem, TWorldDim>& sdv() const {return m_sdv;}

	public:
		inline bool update_local(int edge, const FVElementGeometry<TElem, TWorldDim>& geo)
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			m_edge = edge;

			// normal of this scvf points from corner '_from' to corner '_to'
			for(size_t i = 0; i < m_refElem.num_obj_of_obj(1, edge, 0); ++i)
			{
				m_from = m_refElem.id(1,edge, 0, 0);
				m_to = m_refElem.id(1,edge, 0, 1);
			}

			// compute ip on scvf
			int scale = 0;
			m_ip_local = 0.0;
			for(size_t d = dim; d >= 1; --d)
			{
				for(size_t i = 0; i < m_refElem.num_obj_of_obj(1, edge, d); ++i)
				{
					m_ip_local += geo.obj_midpoint_local(d, m_refElem.id(1, edge, d, i));
					++scale;
				}
			}
			m_ip_local *= 1./scale;

			m_sdv.update_local(m_ip_local);
			return true;
		}

		inline bool update_global(const FVElementGeometry<TElem, TWorldDim>& geo, MathVector<world_dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			// compute ip on scvf
			int scale = 0;
			m_ip_global = 0.0;
			for(size_t d = dim; d >= 1; --d)
			{
				for(size_t i = 0; i < m_refElem.num_obj_of_obj(1, m_edge, d); ++i)
				{
					m_ip_global += geo.obj_midpoint_global(d, m_refElem.id(1, m_edge, d, i));
					++scale;
				}
			}
			m_ip_global *= 1./scale;

			//compute normal on scvf
			if(world_dim == 1)
			{
				m_normal[0] = 1;
			}
			if(world_dim == 2)
			{
				MathVector<world_dim> diff = geo.obj_midpoint_global(2, m_refElem.id(1, m_edge, 2, 0)); // center of element
				diff -= geo.obj_midpoint_global(1, m_refElem.id(1, m_edge, 1, 0)); // edge midpoint

				m_normal[0] = diff[1];
				m_normal[1] = -diff[0];
			}
			else if(world_dim == 3)
			{
				assert(0 && "Must be implemented");
			}

			m_sdv.update_global(m_ip_local, corners);
			return true;
		}

	protected:
		size_t m_edge;
		size_t m_from, m_to;
		MathVector<dim> m_ip_local;

		MathVector<world_dim> m_ip_global;
		MathVector<world_dim> m_normal;
		SD_Values<TElem, TWorldDim> m_sdv;
};



template <	typename TElem,
			int TWorldDim = reference_element_traits<TElem>::dim>
class BoundaryFace{
	public:
		static const int dim = reference_element_traits<TElem>::dim;
		static const int world_dim = TWorldDim;

	public:
		inline size_t local_corner_id() const {return m_co;}
		inline size_t face() const {return m_face;}
		inline number size() const {return m_size;}
		inline size_t num_ip() const {return 1;}
		inline const MathVector<dim>& local_ip() const {return m_ip_local;}
		inline const MathVector<world_dim>& global_ip() const {return m_ip_global;}
		inline const MathVector<world_dim>& normal() const {return m_normal;} // includes area
		inline const SD_Values<TElem, TWorldDim>& sdv() const {return m_sdv;}

	public:
		inline bool update_local(int face, int loc_co, const FVElementGeometry<TElem, TWorldDim>& geo)
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			m_face = face;
			m_co = m_refElem.id(1, m_face, 0, loc_co);

			// compute ip on bf
			int scale = 2;
			m_ip_local = 0.0;
			m_ip_local += geo.corner_local(m_co);
			m_ip_local += geo.obj_midpoint_local(dim-1, m_face);

			if(world_dim == 3)
			{
				assert(0 && "3D not supported, average over edge mid points missing");
			}
			m_ip_local *= 1./scale;

			m_sdv.update_local(m_ip_local);
			return true;
		}

		inline bool update_global(const FVElementGeometry<TElem, TWorldDim>& geo, MathVector<world_dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			// compute ip on bf
			int scale = 2;
			m_ip_global = 0.0;
			m_ip_global += geo.corner_global(m_co);
			m_ip_global += geo.obj_midpoint_global(dim-1, m_face);

			if(world_dim == 3)
			{
				assert(0 && "3D not supported, average over edge mid points missing");
			}
			m_ip_global *= 1./scale;

			//compute normal on bf
			MathVector<world_dim> diff ;
			if(world_dim == 1)
			{
				m_normal[0] = 1;
			}
			if(world_dim == 2)
			{
				diff = geo.obj_midpoint_global(1, m_face); // edge midpoint
				diff -= geo.corner_global(m_co); // corner

				m_normal[0] = diff[1];
				m_normal[1] = -diff[0];
			}
			else if(world_dim == 3)
			{
				assert(0 && "Must be implemented");
			}

			m_size = VecTwoNorm(diff);

			m_sdv.update_global(m_ip_local, corners);
			return true;
		}

	protected:
		size_t m_co;
		size_t m_face;

		MathVector<dim> m_ip_local;

		number m_size;
		MathVector<world_dim> m_ip_global;
		MathVector<world_dim> m_normal;
		SD_Values<TElem, TWorldDim> m_sdv;
};


template <	typename TElem,
			int TWorldDim>
class FVElementGeometry {
	public:
		static const int dim = reference_element_traits<TElem>::dim;
		static const int world_dim = TWorldDim;

	private:
		static const size_t m_num_co = reference_element_traits<TElem>::num_corners;
		static const size_t m_num_ed = reference_element_traits<TElem>::num_edges;
		static const size_t m_num_fa = reference_element_traits<TElem>::num_faces;
		static const size_t m_num_vol = reference_element_traits<TElem>::num_volumes;

	public:
		FVElementGeometry()
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			// compute local corners
			for(size_t i = 0; i < m_num_co; ++i)
			{
				m_co_local[i] = m_refElem.corner(i);
			}

			// compute number of bf for each side
			m_num_bnd_sides = m_refElem.num_obj(dim -1);
			for(size_t i = 0; i < m_refElem.num_obj(dim -1); ++i)
			{
				m_num_bf[i] =  m_refElem.num_obj_of_obj(dim - 1, i, 0);
			}

			// compute local midpoints for all geometric objects with  0 < dim <= m_Dim
			for(int d = 1; d <= dim; ++d)
			{
				// loop geometric objects of dimension d
				for(size_t i = 0; i < m_refElem.num_obj(d); ++i)
				{
					m_obj_mid_local[d-1][i] = 0.0;
					// add corner coordinates of the corners of the geometric object
					for(size_t j = 0; j < m_refElem.num_obj_of_obj(d, i, 0); ++j)
					{
						m_obj_mid_local[d-1][i] += m_co_local[m_refElem.id(d, i, 0, j)];
					}
					// scale for correct averaging
					m_obj_mid_local[d-1][i] *= 1./(m_refElem.num_obj_of_obj(d, i, 0));
				}
			}

			// compute local scvf
			for(size_t i = 0; i < num_scvf(); ++i)
			{
				m_scvf[i].update_local(i, *this);
			}

			// compute local scv
			for(size_t i = 0; i < num_scv(); ++i)
			{
				m_scv[i].update_local(i, *this);
			}

			// compute local bf
			for(size_t side = 0; side < num_bnd_sides(); ++side)
			{
				for(size_t i = 0; i < num_bf(side); ++i)
				{
					m_bf[side][i].update_local(side, i, *this);
				}
			}

		}

		bool update(MathVector<world_dim> corners[])
		{
			static const typename reference_element_traits<TElem>::reference_element_type m_refElem;

			// copy global corners
			for(size_t i = 0; i < m_num_co; ++i)
			{
				m_co_global[i] = corners[i];
			}

			// compute local midpoints for all geometric objects with  0 < dim <= m_Dim
			for(int d = 1; d <= dim; ++d)
			{
				// loop geometric objects of dimension d
				for(size_t i = 0; i < m_refElem.num_obj(d); ++i)
				{
					m_obj_mid_global[d-1][i] = 0.0;
					// add corner coordinates of the corners of the geometric object
					for(size_t j = 0; j < m_refElem.num_obj_of_obj(d, i, 0); ++j)
					{
						m_obj_mid_global[d-1][i] += m_co_global[m_refElem.id(d, i, 0, j)];
					}
					// scale for correct averaging
					m_obj_mid_global[d-1][i] *= 1./(m_refElem.num_obj_of_obj(d, i, 0));
				}
			}

			// compute scvf
			for(size_t i = 0; i < num_scvf(); ++i)
			{
				m_scvf[i].update_global(*this, corners);
			}

			// compute local scv
			for(size_t i = 0; i < num_scv(); ++i)
			{
				m_scv[i].update_global(*this, corners);
			}

			// compute local bf
			for(size_t side = 0; side < num_bnd_sides(); ++side)
			{
				for(size_t i = 0; i < num_bf(side); ++i)
				{
					m_bf[side][i].update_global(*this, corners);
				}
			}

			return true;
		}

	public:
		size_t num_scvf() const {return m_num_ed;};
		size_t num_scv() const {return m_num_co;}
		size_t num_bnd_sides() const {return m_num_bnd_sides;}
		size_t num_bf(size_t side) const {return m_num_bf[side];}

		const SubControlVolumeFace<TElem, TWorldDim>& scvf(size_t i) const {return m_scvf[i];}
		const SubControlVolume<TElem, TWorldDim>& scv(size_t i) const {return m_scv[i];}
		const BoundaryFace<TElem, TWorldDim>& bf(size_t side, size_t i) const {return m_bf[side][i];}

		const MathVector<dim>& corner_local(size_t i) const{return m_co_local[i];}
		const MathVector<dim>& obj_midpoint_local(size_t dim_i, size_t i) const{return m_obj_mid_local[dim_i-1][i];}

		const MathVector<world_dim>& corner_global(size_t i) const{return m_co_global[i];}
		const MathVector<world_dim>& obj_midpoint_global(size_t dim_i, size_t i) const{	return m_obj_mid_global[dim_i-1][i];}

	public:

		static const size_t m_max_obj = m_num_ed;

		size_t m_num_bnd_sides;

		// TODO: Handle this max_sides is not 4 in 3D
		size_t m_num_bf[4];

		MathVector<world_dim> m_co_global[m_num_co];
		MathVector<world_dim> m_obj_mid_global[dim][m_max_obj];

		MathVector<dim> m_co_local[m_num_co];
		MathVector<dim> m_obj_mid_local[dim][m_max_obj];

		SubControlVolume<TElem, TWorldDim> m_scv[m_num_co]; // one scv per corner
		SubControlVolumeFace<TElem, TWorldDim> m_scvf[m_num_ed]; // one scvf per edge
		BoundaryFace<TElem, TWorldDim> m_bf[4][4]; // one bf for each corner of each side
};

}

#endif /* __H__LIB_DISCRETIZATION__DOMAIN_DISCRETIZATION__DISC_HELPER__FVGEOM__ */
