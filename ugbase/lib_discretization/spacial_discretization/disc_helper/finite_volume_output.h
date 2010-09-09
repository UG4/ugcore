/*
 * finite_volume_output.h
 *
 *  Created on: 06.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_OUTPUT__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_OUTPUT__

// other ug4 modules
#include "common/common.h"

#include "lib_grid/lib_grid.h"

// finite volume geometry
#include "./finite_volume_geometry.h"
#include "lib_discretization/domain_util.h"

namespace ug{

////////////////////////////////////////////////////
////////////////////////////////////////////////////
// ConstructGridOfSCVF
////////////////////////////////////////////////////
////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
bool CreateSCVF(const TElem& elem, FV1Geometry<TElem, TWorldDim>& geo, SubsetHandler& shOut,
				Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut)
{
	// extract dimensions
	static const int refDim = FV1Geometry<TElem, TWorldDim>::dim;

	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

	// tmp vector for vertices
	std::vector<VertexBase*> vVert;

	// loop all scv of the element
	for(size_t i = 0; i < geo.num_scvf(); ++i)
	{
		const typename FV1Geometry<TElem, TWorldDim>::SCVF& scvf = geo.scvf(i);

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
			grid.template create<Edge>(EdgeDescriptor(vVert[0], vVert[1]));
		}
		// face
		else if(refDim == 3)
		{
			grid.template create<Quadrilateral>(QuadrilateralDescriptor(vVert[0], vVert[1],
																		vVert[2], vVert[3]));
		}
	}
	return true;
}

template <typename TElem, int TWorldDim>
bool ConstructGridOfSCVF(SubsetHandler& shOut, const SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
		int si)
{
	// Create Geometry
	FV1Geometry<TElem, TWorldDim> geo;

	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

	// iterators for primary grid
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = sh.begin<TElem>(si);
	iterEnd = sh.end<TElem>(si);

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
		geo.update(elem, grid, &vCornerCoords[0]);

		// Create dual grid
		CreateSCVF<TElem, TWorldDim>(*elem, geo, shOut, aaPosOut);
	}

	return true;
}


template <int TWorldDim>
bool ConstructGridOfSCVF(SubsetHandler& shOut, const SubsetHandler& sh,
						Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
						Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
						int si);

template <>
bool ConstructGridOfSCVF<1>(SubsetHandler& shOut, const SubsetHandler& sh,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPos,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPosOut,
							int si)
{
	const int siDim = DimensionOfSubset(sh, si);
	switch(siDim)
	{
		case 1: if(!ConstructGridOfSCVF<Edge, 1>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Edges.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 1 <<".\n");
					return false;
	}
	return true;
}

template <>
bool ConstructGridOfSCVF<2>(SubsetHandler& shOut, const SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPos,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPosOut,
		int si)
{
	const int siDim = DimensionOfSubset(sh, si);
	switch(siDim)
	{
		case 1: if(!ConstructGridOfSCVF<Edge, 2>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Edges.\n"); return false;}
				break;
		case 2: if(!ConstructGridOfSCVF<Triangle, 2>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Triangles.\n"); return false;}
				if(!ConstructGridOfSCVF<Quadrilateral, 2>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Quadrilaterals.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 2 <<".\n");
					return false;
	}
	return true;
}

template <>
bool ConstructGridOfSCVF<3>(SubsetHandler& shOut, const SubsetHandler& sh,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPos,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPosOut,
							int si)
{
	const int siDim = DimensionOfSubset(sh, si);
	switch(siDim)
	{
		case 1: if(!ConstructGridOfSCVF<Edge, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Edges.\n"); return false;}
				break;
		case 2: if(!ConstructGridOfSCVF<Triangle, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Triangles.\n"); return false;}
				if(!ConstructGridOfSCVF<Quadrilateral, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Quadrilaterals.\n"); return false;}
				break;
		case 3: if(!ConstructGridOfSCVF<Tetrahedron, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Tetrahedrons.\n"); return false;}
				if(!ConstructGridOfSCVF<Hexahedron, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Hexahedrons.\n"); return false;}
				if(!ConstructGridOfSCVF<Prism, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Prisms.\n"); return false;}
				if(!ConstructGridOfSCVF<Pyramid, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Pyramids.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 3 <<".\n");
					return false;
	}
	return true;
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
// ConstructGridOfSCV
////////////////////////////////////////////////////
////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
bool CreateSCV(const TElem& elem, FV1Geometry<TElem, TWorldDim>& geo, SubsetHandler& shOut,
				Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut)
{
	// extract dimensions
	static const int refDim = FV1Geometry<TElem, TWorldDim>::dim;

	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

	// tmp vector for vertices
	std::vector<VertexBase*> vVert;

	// loop all scv of the element
	for(size_t i = 0; i < geo.num_scv(); ++i)
	{
		const typename FV1Geometry<TElem, TWorldDim>::SCV& scv = geo.scv(i);

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
			grid.template create<Quadrilateral>(QuadrilateralDescriptor(vVert[0], vVert[1],
																		vVert[2], vVert[3]));
		}
		// volume
		else if(refDim == 3)
		{
			grid.template create<Hexahedron>(HexahedronDescriptor(	vVert[0], vVert[1], vVert[2], vVert[3],
																	vVert[4], vVert[5], vVert[6], vVert[7]));
		}
	}
	return true;
}

template <typename TElem, int TWorldDim>
bool ConstructGridOfSCV(SubsetHandler& shOut, const SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
		int si)
{
	// Create Geometry
	FV1Geometry<TElem, TWorldDim> geo;

	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

	// iterators for primary grid
	typename geometry_traits<TElem>::const_iterator iter, iterBegin, iterEnd;
	iterBegin = sh.begin<TElem>(si);
	iterEnd = sh.end<TElem>(si);

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
		geo.update(elem, grid, &vCornerCoords[0]);

		// Create dual grid
		CreateSCV<TElem, TWorldDim>(*elem, geo, shOut, aaPosOut);
	}

	return true;
}


template <int TWorldDim>
bool ConstructGridOfSCV(SubsetHandler& shOut, const SubsetHandler& sh,
						Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPos,
						Grid::VertexAttachmentAccessor<Attachment<MathVector<TWorldDim> > >& aaPosOut,
						int si);

template <>
bool ConstructGridOfSCV<1>(SubsetHandler& shOut, const SubsetHandler& sh,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPos,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<1> > >& aaPosOut,
							int si)
{
	const int siDim = DimensionOfSubset(sh, si);
	switch(siDim)
	{
		case 1: if(!ConstructGridOfSCV<Edge, 1>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Edges.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 1 <<".\n");
					return false;
	}
	return true;
}

template <>
bool ConstructGridOfSCV<2>(SubsetHandler& shOut, const SubsetHandler& sh,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPos,
		Grid::VertexAttachmentAccessor<Attachment<MathVector<2> > >& aaPosOut,
		int si)
{
	const int siDim = DimensionOfSubset(sh, si);
	switch(siDim)
	{
		case 1: if(!ConstructGridOfSCV<Edge, 2>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Edges.\n"); return false;}
				break;
		case 2: if(!ConstructGridOfSCV<Triangle, 2>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Triangles.\n"); return false;}
				if(!ConstructGridOfSCV<Quadrilateral, 2>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Quadrilaterals.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 2 <<".\n");
					return false;
	}
	return true;
}

template <>
bool ConstructGridOfSCV<3>(SubsetHandler& shOut, const SubsetHandler& sh,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPos,
							Grid::VertexAttachmentAccessor<Attachment<MathVector<3> > >& aaPosOut,
							int si)
{
	const int siDim = DimensionOfSubset(sh, si);
	switch(siDim)
	{
		case 1: if(!ConstructGridOfSCV<Edge, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Edges.\n"); return false;}
				break;
		case 2: if(!ConstructGridOfSCV<Triangle, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Triangles.\n"); return false;}
				if(!ConstructGridOfSCV<Quadrilateral, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Quadrilaterals.\n"); return false;}
				break;
		case 3: if(!ConstructGridOfSCV<Tetrahedron, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Tetrahedrons.\n"); return false;}
				if(!ConstructGridOfSCV<Hexahedron, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Hexahedrons.\n"); return false;}
				if(!ConstructGridOfSCV<Prism, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Prisms.\n"); return false;}
				if(!ConstructGridOfSCV<Pyramid, 3>(shOut, sh, aaPos, aaPosOut, si))
					{UG_LOG("CreateDualGrid: Error while processing Pyramids.\n"); return false;}
				break;
		default: UG_LOG("CreateDualGrid: Dimension " << siDim << " not supported. World dimension is " << 3 <<".\n");
					return false;
	}
	return true;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
// Assignement of Subsets
////////////////////////////////////////////////////
////////////////////////////////////////////////////

template <typename TElem>
bool ColorSubControlVolumeFaces(SubsetHandler& shOut)
{
	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

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
	return true;
}

template <typename TElem>
bool ColorSubControlVolume(SubsetHandler& shOut)
{
	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

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
	return true;
}


template <int TRefDim>
bool ColorControlVolume(SubsetHandler& shOut)
{
	// extract grid
	Grid& grid = *shOut.get_assigned_grid();

	std::vector<Volume*> vVols;
	std::vector<Face*> vFaces;
	std::vector<EdgeBase*> vEdges;

	int si = 0;
	for(VertexBaseIterator iter = shOut.begin<VertexBase>(0);
		iter != shOut.end<VertexBase>(0); ++iter, ++si)
	{
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
		default: UG_LOG("Dimension " << 3 << " is not supported.\n");
		}
	}

	return true;
}


template <typename TAPosition>
bool CreateGridOfSubControlVolume(SubsetHandler& shOut, SubsetHandler& sh, TAPosition& aPos, int si = -1)
{
	static const int dim = TAPosition::ValueType::Size;

	// get assigned grid
	Grid& grid = *sh.get_assigned_grid();
	Grid& gridOut = *shOut.get_assigned_grid();

	// create attachment accessor
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<TAPosition> aaPosOut(gridOut, aPos);

	// check ref dim
	int refDim = dim;

	// Construct dual domain for scv and given subset
	if(si >= 0)
	{
		refDim = DimensionOfSubset(sh, si);

		if(!ConstructGridOfSCV<dim>(shOut, sh, aaPos, aaPosOut, si))
			{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
	}
	// if no subset selected, construct dual grid for all subsets with dim == world_dim
	else
	{
		for(si = 0; si < sh.num_subsets(); ++si)
		{
			if(DimensionOfSubset(sh, si) != dim) continue;

			if(!ConstructGridOfSCV<dim>(shOut, sh, aaPos, aaPosOut, si))
				{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
		}
	}

	// Let each SubControlVolume be one subset
	switch(refDim)
	{
		case 1: if(!ColorSubControlVolume<EdgeBase>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		case 2: if(!ColorSubControlVolume<Face>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		case 3: if(!ColorSubControlVolume<Volume>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		default: UG_LOG("WriteDualGridToFile: Dimension " << refDim << " not supported.\n");
					return false;
	}

	return true;
}

template <typename TAPosition>
bool CreateGridOfControlVolume(SubsetHandler& shOut, SubsetHandler& sh, TAPosition& aPos, int si = -1)
{
	static const int dim = TAPosition::ValueType::Size;

	// get assigned grid
	Grid& grid = *sh.get_assigned_grid();
	Grid& gridOut = *shOut.get_assigned_grid();

	// create attachment accessor
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<TAPosition> aaPosOut(gridOut, aPos);

	// check ref dim
	int refDim = dim;

	// Construct dual domain for scv and given subset
	if(si >= 0)
	{
		refDim = DimensionOfSubset(sh, si);

		if(!ConstructGridOfSCV<dim>(shOut, sh, aaPos, aaPosOut, si))
			{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
	}
	// if no subset selected, construct dual grid for all subsets with dim == world_dim
	else
	{
		for(si = 0; si < sh.num_subsets(); ++si)
		{
			if(DimensionOfSubset(sh, si) != dim) continue;

			if(!ConstructGridOfSCV<dim>(shOut, sh, aaPos, aaPosOut, si))
				{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
		}
	}

	// remove doubles
	RemoveDoubles<dim>(gridOut, gridOut.begin<VertexBase>(), gridOut.end<VertexBase>(), aPos, 1e-5);

	// Let each SubControlVolume be one subset
	switch(refDim)
	{
		case 1: if(!ColorControlVolume<1>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		case 2: if(!ColorControlVolume<2>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		case 3: if(!ColorControlVolume<3>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		default: UG_LOG("WriteDualGridToFile: Dimension " << refDim << " not supported.\n");
					return false;
	}

	return true;
}

template <typename TAPosition>
bool CreateGridOfSubControlVolumeFaces(SubsetHandler& shOut, SubsetHandler& sh, TAPosition& aPos, int si = -1)
{
	static const int dim = TAPosition::ValueType::Size;

	// get assigned grid
	Grid& grid = *sh.get_assigned_grid();
	Grid& gridOut = *shOut.get_assigned_grid();

	// create attachment accessor
	Grid::VertexAttachmentAccessor<TAPosition> aaPos(grid, aPos);
	Grid::VertexAttachmentAccessor<TAPosition> aaPosOut(gridOut, aPos);

	// check ref dim
	int refDim = dim;

	// Construct dual domain for scv and given subset
	if(si >= 0)
	{
		refDim = DimensionOfSubset(sh, si);

		if(!ConstructGridOfSCVF<dim>(shOut, sh, aaPos, aaPosOut, si))
			{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
	}
	// if no subset selected, construct dual grid for all subsets with dim == world_dim
	else
	{
		for(si = 0; si < sh.num_subsets(); ++si)
		{
			if(DimensionOfSubset(sh, si) != dim) continue;

			if(!ConstructGridOfSCVF<dim>(shOut, sh, aaPos, aaPosOut, si))
				{UG_LOG("WriteDualGridToFile: Error while writing subset "<<si<<".\n"); return false;}
		}
	}

	// Let each SubControlVolume be one subset
	switch(refDim)
	{
		case 1: if(!ColorSubControlVolumeFaces<VertexBase>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		case 2: if(!ColorSubControlVolumeFaces<EdgeBase>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		case 3: if(!ColorSubControlVolumeFaces<Face>(shOut))
					{UG_LOG("WriteDualGridToFile: Cannot choose subsets for SubVolumes.\n"); return false;}
				break;
		default: UG_LOG("WriteDualGridToFile: Dimension " << refDim << " not supported.\n");
					return false;
	}

	return true;
}

} // end namespace ug


#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_OUTPUT__ */
