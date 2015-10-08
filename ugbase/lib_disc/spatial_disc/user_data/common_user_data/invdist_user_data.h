/*
 * Inverse Distance Weighting (IDW) interpolation for data sets
 *
 * Created on Feb. 15, 2015 by D. Logashenko
 */
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__IDW_USER_DATA__
#define __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__IDW_USER_DATA__

#include <vector>

// ug4 headers
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/spatial_disc/user_data/std_glob_pos_data.h"

namespace ug {

/// Class for inverse distance weighting based on a general data type.
/**
 * The static functions in the class compute the field given by the IDW
 * interpolation of a given order for given interpolation points and values.
 * Type of the values specified by the template parameter (and may be a scalar
 * or a tensor arithmetic type).
 *
 * Having a set of interpolation points \f$ \mathbf{x}_i \f$ with values
 * \f$ u_i \f$ at them, the interpolated value \f$ u \f$ at point
 * \f$ \mathbf{x} \not\in \{ \mathbf{x}_i \} \f$ is computed as
 * \f{eqnarray*}{
 *  u = \frac {\sum_i w_i \cdot u_i} {\sum_i w_i},
 * \f}
 * where
 * \f{eqnarray*}{
 *  w_i := \frac{1}{\| \mathbf{x} - \mathbf{x}_i \|_2^p}
 * \f}
 * (\f$ p \f$ being the order of the interpolation).
 * For \f$ \mathbf{x} = \mathbf{x}_i \f$ (up to some numerical precision), we
 * set \f$ u = u_i \f$.
 *
 * \remark For \f$ p \f$ less or equal the dimension of the geometric space,
 * the interpolated value \f$ u \f$ may be essentially influenced by the values
 * at \f$ \mathbf{x}_i \f$ located far away from \f$ \mathbf{x} \f$. However,
 * this influence may be avoided by specifying the radius \f$ R \f$: Then only
 * those \f$ \mathbf{x}_i \f$ (and therefore \f$ u_i \f$) are considered in the
 * sums, for which \f$ \| \mathbf{x} - \mathbf{x}_i \|_2 \le R \f$ holds.
 *
 * \remark Note that if the radius \f$ R \f$ is specified, then it can happen
 * that the interpolated value depends only on values at points located on one
 * side of \f$ \mathbf{x} \f$, completely ignoring a trend prescribed by points
 * on the other sides. Thus, \f$ R \f$ should be large enough.
 *
 * To loop the interpolation points, iterators of the templated type TPntIterator
 * are used. Every element of the reference elements should have two members:
 * pos and value:
 * TPntIterator ptr;
 * ptr->pos is a MathVector<WDim> object with the coordinates of the interpolation
 * point;
 * ptr->value is a TData object of the value at that point.
 *
 * \tparam WDim			dimensionality of the space
 * \tparam TPntIterator	interpolation point iterator type
 * \tparam TData		type of the values to interpolate
 */
template <int WDim, typename TPntIterator, typename TData = number>
class IDWInterpolation
{
public:

///	dimensionality of the space (i.e. of the coordinate vectors)
	static const int dim = WDim;
	
///	type of the interpolation point iterator
	typedef TPntIterator t_pnt_iter;
	
///	type of the data to extrapolate
	typedef TData data_type;
	
public:

///	computes the interpolation basing on all the interpolation points
	static void compute
	(
		data_type & res, ///< interpolated value
		const MathVector<dim> & pos, ///< geometric position where to interpolate
		t_pnt_iter pnt_beg, ///< the first interpolation point
		t_pnt_iter pnt_end, ///< delimiter of the iterpolation points
		number order, ///< order of the interpolation
		number small_dist = 1e-7 ///< distance at which we do not distinguish the points
	);

///	computes the interpolation basing on the interpolation points in a given ball
	static void compute
	(
		data_type & res, ///< interpolated value
		const MathVector<dim> & pos, ///< geometric position where to interpolate
		number R, ///< radius of the ball (if 0 then the whole space)
		t_pnt_iter pnt_beg, ///< the first interpolation point
		t_pnt_iter pnt_end, ///< delimiter of the iterpolation points
		number order, ///< order of the interpolation
		number small_dist = 1e-7 ///< distance at which we do not distinguish the points
	);
};

/// UserData interface for the IDW interpolation
/**
 * This class implements the UserData interface for the inverse-distance-weighting
 * interpolation.
 *
 * \sa IDWInterpolation
 *
 * Setting the radius to 0 means the unconstrained version of the IDW interpolation.
 *
 * \tparam WDim		dimensionality of the geometric space (the world dimension)
 * \tparam TData	type of the data to interpolate
 */
template <int WDim, typename TData = number>
class IDWUserData
:	public StdGlobPosData<IDWUserData<WDim, TData>, TData, WDim>
{
public:

///	dimensionality of the space (i.e. of the coordinate vectors)
	static const int dim = WDim;
	
///	type of the data to extrapolate
	typedef TData data_type;
	
private:
	
/// type of a interpolation point data item
	struct data_item
	{
		MathVector<dim> pos; ///< (global) geometrical coordinates of the point
		data_type value; ///< value at that point
		
		data_item (const MathVector<dim> & x, const data_type & v) : pos (x), value (v) {};
		data_item (const data_item & dat) : pos (dat.pos), value (dat.value) {};
	};
	
public:

///	class constructor that creates an empty object with default parameters
	IDWUserData ()
	:	m_order (dim + 1), m_R (0)
	{}
	
///	class constructor that creates an empty object with given parameters
	IDWUserData (number order, number R)
	:	m_order (order), m_R (R)
	{}

///	virtual destructor
	~IDWUserData () {}

public:

///	sets the radius of the neighbourhood where the interpolation points are taken from
	void set_radius (number R) {m_R = R;}
	
///	sets the order of the interpolation
	void set_order (number order) {m_order = order;}
	
///	deletes all the interpolation points from the list
	void clear () {m_data.clear ();}
	
///	loads data from a given stream (and appends the loaded points to the current list)
	void load_data_from (std::istream & in);
	
///	loads data from a given file (and appends the loaded points to the current list)
	void load_data_from (const char * file_name);
	
///	appends an interpolation point to the list
	void append (const MathVector<dim> & x, const data_type & val) {m_data.push_back (data_item (x, val));}
	
public:

///	evaluates the data at a given point
	inline void evaluate (data_type & value, const MathVector<dim> & x, number time, int si) const
	{
		typedef typename std::vector<data_item>::const_iterator pnt_iter_type;
		IDWInterpolation<dim, pnt_iter_type, data_type>::compute (value, x, m_R,
										m_data.begin (), m_data.end (), m_order);
	}

private:

	std::vector<data_item> m_data; ///< interpolation points
	number m_order; ///< order of the interpolation
	number m_R; ///< radius of the neighbourhood to look for the interpolation points in (0 == infinite)
};

} // end namespace ug

#include "invdist_user_data_impl.h"

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__USER_DATA__IDW_USER_DATA__

/* End of File */
