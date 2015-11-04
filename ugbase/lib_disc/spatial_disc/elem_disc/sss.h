
#ifndef __H_SINGULAR_SOURCES_AND_SINKS__
#define __H_SINGULAR_SOURCES_AND_SINKS__

#include <vector>
#include <list>

#include "lib_grid/algorithms/bounding_box_util.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/ray_element_intersection_util.h"

#ifdef UG_FOR_LUA
#include "bindings/lua/lua_user_data.h"
#endif

namespace ug {

template<size_t dim, size_t non>
class PointSourceOrSink
{
public:
	MathVector<dim> point;
    size_t marker;
	MathVector<non> data;
	SmartPtr<UserData<MathVector<non>, dim> > spData;

	PointSourceOrSink(const MathVector<dim>& _point,
					  const MathVector<non>& _data)

	    : point(_point), marker(0), data(_data) {}

    PointSourceOrSink(const MathVector<dim>& _point,
					  SmartPtr<UserData<MathVector<non>, dim> > _spData)
	
	    : point(_point), marker(0), spData(_spData) {}
};

template<size_t dim, size_t non>
class LineSourceOrSink
{
public:
    MathVector<dim> point1;
    MathVector<dim> point2;
    MathVector<non> data;
    SmartPtr<UserData<MathVector<non>, dim> > spData;

    LineSourceOrSink(const MathVector<dim>& _point1,
					 const MathVector<dim>& _point2,
					 const MathVector<non>& _data)

		: point1(_point1), point2(_point2), data(_data) {}

	LineSourceOrSink(const MathVector<dim>& _point1,
					 const MathVector<dim>& _point2,
					 SmartPtr<UserData<MathVector<non>, dim> > _spData)

		: point1(_point1), point2(_point2), spData(_spData) {}
};

template<size_t dim, size_t non>
class SingularSourcesAndSinks
{
private:
	typename std::list<PointSourceOrSink<dim, non> > ListP;
	typename std::list<LineSourceOrSink<dim, non>  > ListL;

public:
	void addps(const std::vector<number>& _point,
			   const std::vector<number>& _data)
	{
		if (_point.size() != dim)
			UG_THROW("Expected vector of dimension size (dim: " << dim << ").");

		MathVector<dim> point;
		for(size_t i = 0; i < dim; i++)
			point[i] = _point[i];

		if (_data.size() != non)
			UG_THROW("Expected vector of data size (non: " << non << ").");
		
		MathVector<non> data;
        for(size_t i = 0; i < non; i++)
            data[i] = _data[i];

		ListP.push_back(PointSourceOrSink<dim, non>(point, data));
	}

#ifdef UG_FOR_LUA
	void addps(const std::vector<number>& _point,
			   LuaFunctionHandle func)
	{
		if (_point.size() != dim)
			UG_THROW("Expected vector of dimension size (dim: " << dim << ").");

		MathVector<dim> point;
		for(size_t i = 0; i < dim; i++)
			point[i] = _point[i];

		SmartPtr<UserData<MathVector<non>, dim> > spData = make_sp(new LuaUserData<MathVector<non>, dim>(func));
		ListP.push_back(PointSourceOrSink<dim, non>(point, spData));
	}
#endif

	void addls(const std::vector<number>& _point1,
			   const std::vector<number>& _point2,
			   const std::vector<number>& _data)
	{
		if (_point1.size() != dim || _point2.size() != dim)
			UG_THROW("Expected vector of dimension size (dim: " << dim << ").");

		MathVector<dim> point1, point2;
		for(size_t i = 0; i < dim; i++) {
			point1[i] = _point1[i];
			point2[i] = _point2[i];
		}

		if (_data.size() != non)
            UG_THROW("Expected vector of data size (non: " << non << ").");

        MathVector<non> data;
        for(size_t i = 0; i < non; i++)
            data[i] = _data[i];

		ListL.push_back(LineSourceOrSink<dim, non>(point1, point2, data));
	}

#ifdef UG_FOR_LUA
	void addls(const std::vector<number>& _point1,
			   const std::vector<number>& _point2,
			   LuaFunctionHandle func)
	{
		if (_point1.size() != dim || _point2.size() != dim)
			UG_THROW("Expected vector of dimension size (dim: " << dim << ").");

		MathVector<dim> point1, point2;
		for(size_t i = 0; i < dim; i++) {
			point1[i] = _point1[i];
			point2[i] = _point2[i];
		}

		SmartPtr<UserData<MathVector<non>, dim> > spData = make_sp(new LuaUserData<MathVector<non>, dim>(func));
		ListL.push_back(LineSourceOrSink<dim, non>(point1, point2, spData));
	}
#endif

    void clear_markers()
    {
		typename std::list<PointSourceOrSink<dim, non> >::iterator p;
        for (p = ListP.begin(); p != ListP.end(); p++)
			p->marker = 0;
	}

    template<typename TElem,
	         typename TAAPos,
		     typename TFVGeom>
	number
	get_contrib_of_scv(TElem* elem,
					   Grid& grid,
					   TAAPos aaPos,
					   TFVGeom& geo,
					   size_t scv,
					   number time,
					   MathVector<non>& out)
	{
		typename std::list<PointSourceOrSink<dim, non> >::iterator p;
		AABox<MathVector<dim> > bbox = CalculateBoundingBox(elem, aaPos);
		size_t scv_id = (size_t)elem + scv;
		for (p = ListP.begin(); p != ListP.end(); p++) {
			bool in = true;
			if (p->marker != 0) {
				if (p->marker != scv_id)
					continue;
			}
			else {
				// quick bbox test
				if (!bbox.contains_point(p->point))
					continue;

				// restrict point to element
				if (!ContainsPoint(elem, p->point, aaPos))
					continue;

				// restrict point to scv
				for (size_t i = 0; i < geo.num_scvf(); i++) {
					const typename TFVGeom::SCVF& scvf = geo.scvf(i);
					MathVector<dim> d;
					if (scvf.from() == scv) {
						VecSubtract(d, p->point, scvf.global_ip());
						if (VecDot(d, scvf.normal()) > 0.0) {
							in = false;
							break;
						}
					}
					else if (scvf.to() == scv) {
						VecSubtract(d, p->point, scvf.global_ip());
						if (VecDot(d, scvf.normal()) < 0.0) {
							in = false;
							break;
						}
					}
				}
			}
			if (in) {
				p->marker = scv_id;
				if (p->spData.invalid())
					out = p->data;
				else 
					(*p->spData)(out, p->point, time, 0);
				return 1.0;
			}
		}

		typename std::list<LineSourceOrSink<dim, non> >::iterator l;
        for (l = ListL.begin(); l != ListL.end(); l++)
		{
			// quick bbox test
			if (!bbox.overlaps_line(l->point1, l->point2))
				continue;

			// restrict line segment to element
			MathVector<dim> ls = l->point1;
			MathVector<dim> le = l->point2;
			MathVector<dim> dir;
			number lambda_min, lambda_max;
			VecSubtract(dir, le, ls);
			if (!RayElementIntersection(lambda_min, lambda_max, ls, dir, elem, grid, aaPos))
				continue;
			if (ContainsPoint(elem, ls, aaPos))
				lambda_min = 0.0;
			if (ContainsPoint(elem, le, aaPos))
				lambda_max = 1.0;
			if (lambda_min < 0.0 || lambda_max > 1.0)
				continue;
			VecScaleAdd(le, 1.0, ls, lambda_max, dir);
			VecScaleAdd(ls, 1.0, ls, lambda_min, dir);

			// restrict line segment to scv
			MathVector<dim> p;
			bool in = true;
			for (size_t i = 0; i < geo.num_scvf(); i++) {
				const typename TFVGeom::SCVF& scvf = geo.scvf(i);
				number d_min, d_max, lambda;
				if (scvf.from() == scv) {
					VecSubtract(p, ls, scvf.global_ip());
					d_min = VecDot(p, scvf.normal());
					VecSubtract(p, le, scvf.global_ip());
					d_max = VecDot(p, scvf.normal());
					if (d_min*d_max < 0.0) {
						lambda = fabs(d_min)/(fabs(d_min)+fabs(d_max));
						VecScaleAdd(p, 1.0-lambda, ls, lambda, le);
						if (d_max > 0.0)
							le = p;
						else
							ls = p;
					}
					else if (d_min > 0.0) {
						in = false;
						break;
					}
				}
				else if (scvf.to() == scv) {
					VecSubtract(p, ls, scvf.global_ip());
					d_min = VecDot(p, scvf.normal());
					VecSubtract(p, le, scvf.global_ip());
					d_max = VecDot(p, scvf.normal());
					if (d_min*d_max < 0.0) {
						lambda = fabs(d_min)/(fabs(d_min)+fabs(d_max));
						VecScaleAdd(p, 1.0-lambda, ls, lambda, le);
						if (d_max < 0.0)
							le = p;
						else
							ls = p;
					}
					else if (d_min < 0.0) {
						in = false;
						break;
					}
				}
			}
			if (in) {
				VecSubtract(p, le, ls);
				number len = VecLength(p);
				if (l->spData.invalid())
					out = l->data;
				else
					(*l->spData)(out, ls, time, 0);
				return len;
			}
		}
		return 0.0;
	}
};

template<size_t non>
class SingularSourcesAndSinks<1, non>
{
public:

	void addps(const std::vector<number>& _point,
			   const std::vector<number>& _data)
	{
		UG_THROW("1d unsupported.");
	}

#ifdef UG_FOR_LUA
	void addps(const std::vector<number>& _point,
			   LuaFunctionHandle func)
	{
		UG_THROW("1d unsupported.");
	}
#endif

	void addls(const std::vector<number>& _point1,
			   const std::vector<number>& _point2,
			   const std::vector<number>& _data)
	{
		UG_THROW("1d unsupported.");	
	}

#ifdef UG_FOR_LUA
	void addls(const std::vector<number>& _point1,
			   const std::vector<number>& _point2,
			   LuaFunctionHandle func)
	{
		UG_THROW("1d unsupported.");
	}
#endif

	void clear_markers()
    {
		UG_THROW("1d unsupported.");	
	}

	template<typename TElem,
	         typename TAAPos,
		     typename TFVGeom>
	number
	get_contrib_of_scv(TElem* elem,
					   Grid& grid,
					   TAAPos aaPos,
					   TFVGeom& geo,
					   size_t scv,
					   number time,
					   MathVector<non>& out)
	{
		UG_THROW("1d unsupported.");
		return 0.0;
	}
};

} // end namespace ug

#endif
