/*
 * element_list_traits.h
 * Classes for computing integral properties of geometric objects types in lists.
 *
 * Created by D. Logashenko
 * Sep. 25, 2013
 */
#ifndef __H__UG__LIB_DISC__REFERENCE_ELEMENT__ELEMENT_LIST_TRAITS__
#define __H__UG__LIB_DISC__REFERENCE_ELEMENT__ELEMENT_LIST_TRAITS__

// further lib_disc headers
#include "lib_disc/reference_element/reference_element.h"

// Boost-C++ headers
#include <boost/mpl/transform_view.hpp>
#include <boost/mpl/fold.hpp>
#include <boost/mpl/min_max.hpp>

namespace ug{

/// Class for computing integral properties of geometric objects types in lists
/**
 * This class gets the maximum number of corners, edges etc over all the geometric
 * elements in a given list. The computation is performed at the compile time
 * so that the result (e.g. element_list_traits<ElemList>::max_edges) is a static
 * constant.
 *
 * \tparam ElemList		list of geometric elements
 */
template <typename ElemList>
class element_list_traits
{
///	Metafunction class for counting corners in an element type
	struct mfc_num_corners_of_elem
	{
		template <typename TElem> class apply
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			public: typedef boost::mpl::int_<ref_elem_type::numCorners> type; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting edges in an element type
	struct mfc_num_edges_of_elem
	{
		template <typename TElem> class apply
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			public: typedef boost::mpl::int_<ref_elem_type::numEdges> type; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting faces in an element type
	struct mfc_num_faces_of_elem
	{
		template <typename TElem> class apply
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			public: typedef boost::mpl::int_<ref_elem_type::numFaces> type; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting volumes in an element type
	struct mfc_num_volumes_of_elem
	{
		template <typename TElem> class apply
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			public: typedef boost::mpl::int_<ref_elem_type::numVolumes> type; ///< returned type (i.e. result of the metafunction)
		};
	};
	
///	Metafunction class for counting sides (i.e. edges or faces) in an element type
	struct mfc_num_sides_of_elem
	{
		template <typename TElem> class apply
		{
			typedef typename reference_element_traits<TElem>::reference_element_type ref_elem_type;
			public: typedef boost::mpl::int_<ref_elem_type::numSides> type; ///< returned type (i.e. result of the metafunction)
		};
	};
	
public:
	
/// Max. number of corners of the elements in the element list (as a constant)
	static const int maxCorners
		= boost::mpl::fold
			<
				boost::mpl::transform_view<ElemList, mfc_num_corners_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of edges of the elements in the element list (as a constant)
	static const int maxEdges
		= boost::mpl::fold
			<
				boost::mpl::transform_view<ElemList, mfc_num_edges_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of faces of the elements in the element list (as a constant)
	static const int maxFaces
		= boost::mpl::fold
			<
				boost::mpl::transform_view<ElemList, mfc_num_faces_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of volumes of the elements in the element list (as a constant)
	static const int maxVolumes
		= boost::mpl::fold
			<
				boost::mpl::transform_view<ElemList, mfc_num_volumes_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
	
/// Max. number of sides (edges or faces) of the elements in the element list (as a constant)
	static const int maxSides
		= boost::mpl::fold
			<
				boost::mpl::transform_view<ElemList, mfc_num_sides_of_elem>,
				boost::mpl::int_<0>,
				boost::mpl::max<boost::mpl::_1,boost::mpl::_2>
			>::type::value;
};

} // end namespace ug

#endif // __H__UG__LIB_DISC__REFERENCE_ELEMENT__ELEMENT_LIST_TRAITS__

/* End of File */
