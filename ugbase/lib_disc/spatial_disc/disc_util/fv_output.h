/*
 * fv_output.h
 *
 *  Created on: 06.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_OUTPUT__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_OUTPUT__

// other ug4 modules
#include "common/common.h"
#include "lib_disc/domain_util.h"

// finite volume geometry
#include "fv1_geom.h"
#include "hfv1_geom.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// ConstructGridOfSCVF
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, template <class, int> class TFVGeom, int TWorldDim>
void CreateSCVF(const TElem& elem, TFVGeom<TElem, TWorldDim>& geo, ISubsetHandler& shOut,
				Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut)
{
	// extract dimensions
	static const int refDim = TFVGeom<TElem, TWorldDim>::dim;

	// extract grid
	Grid& grid = *shOut.grid();

	// tmp vector for vertices
	std::vector<VertexBase*> vVert;

	// loop all scv of the element
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		const typename TFVGeom<TElem, TWorldDim>::SCVF& scvf = geo.scvf(i);

		// clear vertices
		vVert.clear();

		// loop corners of scv
		for(size_t co = 0; co < scvf.num_corners(); ++co)
		{
			//	create a new vertex
				Vertex* vrt = *(grid.create<Vertex>());
				vVert.push_back(vrt);

			//	set the coordinates
				aaPosOut[vrt] = scvf.global_corner(co);
		}

		// edge
		if(refDim == 2)
		{
			if (scvf.num_corners() == 2) {
				grid.template create<Edge>(EdgeDescriptor(vVert[0], vVert[1]));
			}
			else {
				UG_THROW("SCVF has a number of nodes, that is not drawable.");
			}
		}
		// face
		else if(refDim == 3)
		{
			if (scvf.num_corners() == 4) {
				grid.template create<Quadrilateral>(QuadrilateralDescriptor(vVert[0], vVert[1],
																			vVert[2], vVert[3]));
			}
			else if (scvf.num_corners() == 3) {
				grid.template create<Triangle>(TriangleDescriptor(vVert[0], vVert[1],
																			vVert[2]));
			}
			else{
				UG_THROW("SCVF has a number of nodes, that is not drawable.");
			}
		}
	}
}

template <typename TElem, template <class, int> class TFVGeom, int TWorldDim>
void ConstructGridOfSCVF(ISubsetHandler& shOut,
                         const SurfaceView& surfView,
                         const Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
                         Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
                         int si)
{
	// Create Geometry
	TFVGeom<TElem, TWorldDim> geo;

	// iterators for primary grid
	typename SurfaceView::traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = surfView.begin<TElem>();
	iterEnd = surfView.end<TElem>();

	// corners of element
	std::vector<MathVector<TWorldDim> > vCornerCoords;

	// iterate over primary grid
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get element
		TElem* elem = *iter;

		// get corner coordinates
		CollectCornerCoordinates(vCornerCoords, *elem, aaPos);

		// update finite volume geometry
		geo.update(elem, &vCornerCoords[0], &(*surfView.subset_handler()));

		// Create dual grid
		CreateSCVF<TElem, TFVGeom, TWorldDim>(*elem, geo, shOut, aaPosOut);
	}
}


template <template <class, int> class TFVGeom, int TWorldDim>
struct ConstructGridOfSCVFWrapper{};

template <template <class, int> class TFVGeom>
struct ConstructGridOfSCVFWrapper<TFVGeom, 1>
{
	static void apply(ISubsetHandler& shOut, const SurfaceView& surfView,
	                  const Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPos,
	                  Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPosOut,
	                  int si, int siDim)
	{
		switch(siDim)
		{
			case 1: ConstructGridOfSCVF<Edge, TFVGeom, 1>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			default: UG_THROW("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 1 <<".");
		}
	}
};

template <template <class, int> class TFVGeom>
struct ConstructGridOfSCVFWrapper<TFVGeom, 2>
{
	static void apply(ISubsetHandler& shOut, const SurfaceView& surfView,
				const Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPos,
				Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPosOut,
				int si, int siDim)
	{
		switch(siDim)
		{
			case 1: ConstructGridOfSCVF<Edge, TFVGeom, 2>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			case 2: ConstructGridOfSCVF<Triangle, TFVGeom, 2>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCVF<Quadrilateral, TFVGeom, 2>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			default: UG_THROW("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 2);
		}
	}
};

template <template <class, int> class TFVGeom>
struct ConstructGridOfSCVFWrapper<TFVGeom, 3>
{
	static void apply(ISubsetHandler& shOut, const SurfaceView& surfView,
				const Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPos,
				Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPosOut,
				int si, int siDim)
	{
		switch(siDim)
		{
			case 1: ConstructGridOfSCVF<Edge, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			case 2: ConstructGridOfSCVF<Triangle, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCVF<Quadrilateral, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			case 3: ConstructGridOfSCVF<Tetrahedron, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCVF<Hexahedron, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCVF<Prism, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCVF<Pyramid, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			default: UG_THROW("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 3);
		}
	}
};

template <template <class, int> class TFVGeom, int TWorldDim>
void ConstructGridOfSCVF(ISubsetHandler& shOut, const SurfaceView& surfView,
						const Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
						Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
						int si, int siDim)
{
	ConstructGridOfSCVFWrapper<TFVGeom, TWorldDim>::apply(shOut, surfView, aaPos, aaPosOut, si, siDim);
}


////////////////////////////////////////////////////////////////////////////////
// ConstructGridOfSCV
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, template <class, int> class TFVGeom, int TWorldDim>
void CreateSCV(const TElem& elem, TFVGeom<TElem, TWorldDim>& geo, ISubsetHandler& shOut,
				Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut)
{
	// extract dimensions
	static const int refDim = TFVGeom<TElem, TWorldDim>::dim;

	// extract grid
	Grid& grid = *shOut.grid();

	// tmp vector for vertices
	std::vector<VertexBase*> vVert;

	// loop all scv of the element
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		const typename TFVGeom<TElem, TWorldDim>::SCV& scv = geo.scv(i);

		// clear vertices
		vVert.clear();

		// loop corners of scv
		for(size_t co = 0; co < scv.num_corners(); ++co)
		{
			//	create a new vertex
				Vertex* vrt = *(grid.create<Vertex>());
				vVert.push_back(vrt);

			//	set the coordinates
				aaPosOut[vrt] = scv.global_corner(co);

			// 	if co == 0, it is a vertex of the primary grid
			//	We assign all of those vertices to subset 0
			//	The other vertices remain in subset -1
				if(co == 0) shOut.assign_subset(vrt, 0);
				else shOut.assign_subset(vrt, -1);

		}

		// edge
		if(refDim == 1)
		{
			grid.template create<Edge>(EdgeDescriptor(vVert[0], vVert[1]));
		}
		// face
		else if(refDim == 2)
		{
			if (scv.num_corners() == 4) {
				grid.template create<Quadrilateral>(QuadrilateralDescriptor(vVert[0], vVert[1],
																		vVert[2], vVert[3]));
			}
			else {
				UG_THROW("SCV has a number of nodes, that is not drawable.");
			}
		}
		// volume
		else if(refDim == 3)
		{
			if (scv.num_corners() == 8) {
				grid.template create<Hexahedron>(HexahedronDescriptor(	vVert[0], vVert[1], vVert[2], vVert[3],
																		vVert[4], vVert[5], vVert[6], vVert[7]));
			}
			else if (scv.num_corners() == 4) {
				grid.template create<Tetrahedron>(TetrahedronDescriptor(	vVert[0], vVert[1], vVert[2], vVert[3]));
			}
			else {
				UG_THROW("SCV has a number of nodes, that is not drawable.");
			}
		}
	}
}

template <typename TElem, template <class, int> class TFVGeom, int TWorldDim>
void ConstructGridOfSCV(ISubsetHandler& shOut, const SurfaceView& surfView,
                        const Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
                        Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
                        int si)
{
	// Create Geometry
	TFVGeom<TElem, TWorldDim> geo;

	// iterators for primary grid
	typename SurfaceView::traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = surfView.begin<TElem>();
	iterEnd = surfView.end<TElem>();

	// corners of element
	std::vector<MathVector<TWorldDim> > vCornerCoords;

	// iterate over primary grid
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
		// get element
		TElem* elem = *iter;

		// get corner coordinates
		CollectCornerCoordinates(vCornerCoords, *elem, aaPos);

		// update finite volume geometry
		geo.update(elem, &vCornerCoords[0], &(*surfView.subset_handler()));

		// Create dual grid
		CreateSCV<TElem, TFVGeom, TWorldDim>(*elem, geo, shOut, aaPosOut);
	}
}


template <template <class, int> class TFVGeom, int TWorldDim>
struct ConstructGridOfSCVWrapper{};

template <template <class TElem, int TWorldDim> class TFVGeom>
struct ConstructGridOfSCVWrapper<TFVGeom, 1>
{
	static void apply(ISubsetHandler& shOut, const SurfaceView& surfView,
	                  const Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPos,
	                  Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPosOut,
	                  int si, int siDim)
	{
		switch(siDim)
		{
			case 1: ConstructGridOfSCV<Edge, TFVGeom, 1>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			default: UG_THROW("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 1);
		}
	}
};

template <template <class, int> class TFVGeom>
struct ConstructGridOfSCVWrapper<TFVGeom, 2>
{
	static void apply(ISubsetHandler& shOut, const SurfaceView& surfView,
	                  const Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPos,
	                  Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPosOut,
	                  int si, int siDim)
	{
		switch(siDim)
		{
			case 1: ConstructGridOfSCV<Edge, TFVGeom, 2>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			case 2: ConstructGridOfSCV<Triangle, TFVGeom, 2>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCV<Quadrilateral, TFVGeom, 2>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			default: UG_THROW("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 2);
		}
	}
};

template <template <class, int> class TFVGeom>
struct ConstructGridOfSCVWrapper<TFVGeom, 3>
{
	static void apply(ISubsetHandler& shOut, const SurfaceView& surfView,
	                  const Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPos,
	                  Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPosOut,
	                  int si, int siDim)
	{
		switch(siDim)
		{
			case 1: ConstructGridOfSCV<Edge, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			case 2: ConstructGridOfSCV<Triangle, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCV<Quadrilateral, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			case 3: ConstructGridOfSCV<Tetrahedron, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCV<Hexahedron, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCV<Prism, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					ConstructGridOfSCV<Pyramid, TFVGeom, 3>(shOut, surfView, aaPos, aaPosOut, si);
					break;
			default: UG_THROW("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 3);
		}
	}
};

template <template <class, int> class TFVGeom, int TWorldDim>
void ConstructGridOfSCV(ISubsetHandler& shOut, const SurfaceView& surfView,
						const Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
						Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
						int si, int siDim)
{
	ConstructGridOfSCVWrapper<TFVGeom, TWorldDim>::apply(shOut, surfView, aaPos, aaPosOut, si, siDim);
}


////////////////////////////////////////////////////////////////////////////////
// Assignement of Subsets
////////////////////////////////////////////////////////////////////////////////

template <typename TElem>
void ColorSubControlVolumeFaces(ISubsetHandler& shOut)
{
	// extract grid
	Grid& grid = *shOut.grid();

	// iterators for primary grid
	typename geometry_traits<TElem>::iterator iter, iterBegin, iterEnd;
	iterBegin = grid.begin<TElem>();
	iterEnd = grid.end<TElem>();

	// iterate over primary grid
	int si = 0;
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
		shOut.assign_subset(*iter, si++);
	}
}

template <typename TElem>
void ColorSubControlVolume(ISubsetHandler& shOut)
{
	// extract grid
	Grid& grid = *shOut.grid();

	// iterators for primary grid
	typename geometry_traits<TElem>::iterator iter, iterBegin, iterEnd;
	iterBegin = grid.begin<TElem>();
	iterEnd = grid.end<TElem>();

	// iterate over primary grid
	int si = 0;
	for(iter = iterBegin; iter != iterEnd; ++iter)
	{
		shOut.assign_subset(*iter, si++);
	}
}


template <int TRefDim>
void ColorControlVolume(ISubsetHandler& shOut)
{
	// extract grid
	Grid& grid = *shOut.grid();

	std::vector<Volume*> vVols;
	std::vector<Face*> vFaces;
	std::vector<EdgeBase*> vEdges;

	int si = 0;
	for(VertexBaseIterator iter = shOut.grid()->begin<VertexBase>();
		iter != shOut.grid()->end<VertexBase>(); ++iter, ++si)
	{
		if(shOut.get_subset_index(*iter) != 0) continue;

		switch(TRefDim)
		{
		case 1:	CollectEdges(vEdges, grid, *iter);
				shOut.assign_subset(vEdges.begin(), vEdges.end(), si);
				break;
		case 2:	CollectFaces(vFaces, grid, *iter);
				shOut.assign_subset(vFaces.begin(), vFaces.end(), si);
				break;
		case 3:	CollectVolumes(vVols, grid, *iter);
				shOut.assign_subset(vVols.begin(), vVols.end(), si);
				break;
		default: UG_THROW("Dimension " << TRefDim << " is not supported.");
		}
	}
}


template <template <class, int> class TFVGeom, typename TAAPosition>
void CreateGridOfSubControlVolumes(ISubsetHandler& shOut, TAAPosition& aaPosOut, const ISubsetHandler& sh, const TAAPosition& aaPos, const SurfaceView& surfView, int si = -1)
{
	static const int dim = TAAPosition::ValueType::Size;

	// Construct dual domain for scv and given subset
	if(si >= 0)
	{
		const int siDim = DimensionOfSubset(sh, si);
		ConstructGridOfSCV<TFVGeom, dim>(shOut, surfView, aaPos, aaPosOut, si, siDim);
	}
	// if no subset selected, construct dual grid for all subsets with dim == worldDim
	else
	{
		for(si = 0; si < sh.num_subsets(); ++si)
		{
			const int siDim = DimensionOfSubset(sh, si);
			if(siDim != dim) continue;

            ConstructGridOfSCV<TFVGeom, dim>(shOut, surfView, aaPos, aaPosOut, si, siDim);
		}
	}

	// Let each SubControlVolume be one subset
	switch(dim)
	{
		case 1: ColorSubControlVolume<EdgeBase>(shOut); break;
		case 2: ColorSubControlVolume<Face>(shOut); break;
		case 3: ColorSubControlVolume<Volume>(shOut); break;
		default: UG_THROW("WriteDualGridToFile: Dimension "<<dim<<" not supported.");
	}
}

template <template <class, int> class TFVGeom, typename TAAPosition, typename TAPosition>
void CreateGridOfControlVolumes(ISubsetHandler& shOut, TAAPosition& aaPosOut, TAPosition& aPosOut, const ISubsetHandler& sh, const TAAPosition& aaPos, const SurfaceView& surfView, int si = -1)
{
	static const int dim = TAAPosition::ValueType::Size;

	// Construct dual domain for scv and given subset
	if(si >= 0)
	{
		const int siDim = DimensionOfSubset(sh, si);
		ConstructGridOfSCV<TFVGeom, dim>(shOut, surfView, aaPos, aaPosOut, si, siDim);
	}
	// if no subset selected, construct dual grid for all subsets with dim == worldDim
	else
	{
		for(si = 0; si < sh.num_subsets(); ++si)
		{
			const int siDim = DimensionOfSubset(sh, si);
			if(siDim != dim) continue;

            ConstructGridOfSCV<TFVGeom, dim>(shOut, surfView, aaPos, aaPosOut, si, siDim);
		}
	}

	// remove doubles
	RemoveDoubles<dim>(*shOut.grid(),
	                   shOut.grid()->begin<VertexBase>(), shOut.grid()->end<VertexBase>(),
	                   aPosOut, 1e-5);

	// Let each SubControlVolume be one subset
	switch(dim)
	{
		case 1: ColorControlVolume<1>(shOut); break;
		case 2: ColorControlVolume<2>(shOut); break;
		case 3: ColorControlVolume<3>(shOut); break;
		default: UG_THROW("WriteDualGridToFile: Dimension "<<dim<<" not supported.");
	}
}

template <template <class, int> class TFVGeom, typename TAAPosition>
void CreateGridOfSubControlVolumeFaces(ISubsetHandler& shOut, TAAPosition& aaPosOut, const ISubsetHandler& sh, const TAAPosition& aaPos, const SurfaceView& surfView, int si = -1)
{
	static const int dim = TAAPosition::ValueType::Size;

	// Construct dual domain for scv and given subset
	if(si >= 0)
	{
		const int siDim = DimensionOfSubset(sh, si);

		ConstructGridOfSCVF<TFVGeom, dim>(shOut, surfView, aaPos, aaPosOut, si, siDim);
	}
	// if no subset selected, construct dual grid for all subsets with dim == worldDim
	else
	{
		for(si = 0; si < sh.num_subsets(); ++si)
		{
			const int siDim = DimensionOfSubset(sh, si);
			if(siDim != dim) continue;

			ConstructGridOfSCVF<TFVGeom, dim>(shOut, surfView, aaPos, aaPosOut, si, siDim);
		}
	}

	// Let each SubControlVolume be one subset
	switch(dim)
	{
		case 1: ColorSubControlVolumeFaces<VertexBase>(shOut); break;
		case 2: ColorSubControlVolumeFaces<EdgeBase>(shOut); break;
		case 3: ColorSubControlVolumeFaces<Face>(shOut); break;
		default: UG_THROW("WriteDualGridToFile: Dimension " << dim << " not supported.");
	}
}

/**
 * Creates a grid consisting of the sub-control-volume faces w.r.t to a given
 * domain.
 *
 * @param domOut	domain to be filled
 * @param domIn		original domain
 * @param si		subset used (-1 for whole domain)
 */
template <template <class, int> class TFVGeom, typename TDomain>
void CreateSubControlVolumeFaceDomain(TDomain& domOut, const TDomain& domIn, const SurfaceView& surfView, int si = -1)
{
	if(&domOut == &domIn)
		UG_THROW("CreateSubControlVolumeFaceDomain: Domains must be different.");

	CreateGridOfSubControlVolumeFaces<TFVGeom, typename TDomain::position_accessor_type>
		( *domOut.subset_handler(),
		  domOut.position_accessor(),
		  *domIn.subset_handler(),
		  domIn.position_accessor(),
		  surfView, si);
}

/**
 * Creates a grid consisting of the sub-control-volume w.r.t to a given
 * domain.
 *
 * @param domOut	domain to be filled
 * @param domIn		original domain
 * @param si		subset used (-1 for whole domain)
 */
template <template <class, int> class TFVGeom, typename TDomain>
void CreateSubControlVolumeDomain(TDomain& domOut, const TDomain& domIn, const SurfaceView& surfView, int si = -1)
{
	if(&domOut == &domIn)
		UG_THROW("CreateSubControlVolumeDomain: Domains must be different.");

	CreateGridOfSubControlVolumes<TFVGeom, typename TDomain::position_accessor_type>
		( *domOut.subset_handler(),
		  domOut.position_accessor(),
		  *domIn.subset_handler(),
		  domIn.position_accessor(),
		  surfView, si);
}

/**
 * Creates a grid consisting of the sub-control-volume faces w.r.t to a given
 * domain.
 *
 * @param domOut	domain to be filled
 * @param domIn		original domain
 * @param si		subset used (-1 for whole domain)
 */
template <template <class, int> class TFVGeom, typename TDomain>
void CreateControlVolumeDomain(TDomain& domOut, const TDomain& domIn, const SurfaceView& surfView, int si = -1)
{
	if(&domOut == &domIn)
		UG_THROW("CreateControlVolumeDomain: Domains must be different.");

	CreateGridOfControlVolumes<TFVGeom, typename TDomain::position_accessor_type, typename TDomain::position_attachment_type>
		( *domOut.subset_handler(),
		  domOut.position_accessor(),
		  domOut.position_attachment(),
		  *domIn.subset_handler(),
		  domIn.position_accessor(),
		  surfView, si);
}

} // end namespace ug


#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_OUTPUT__ */
