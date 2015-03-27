// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include <algorithm>
#include <vector>
#include "horizontal_layers_mesher.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/algorithms/extrusion/extrude.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "lib_grid/algorithms/geom_obj_util/vertex_util.h"
#include "lib_grid/file_io/file_io_asc.h"
#include "lib_grid/iterators/associated_elements_iterator.h"

using namespace std;

namespace ug{


void RasterLayers::
load_from_files(const std::vector<std::string>& filenames)
{
	resize(filenames.size());
	for(size_t i = 0; i < filenames.size(); ++i){
		LoadHeightfieldFromASC(*m_layers[i], filenames[i].c_str());
	}
}

void RasterLayers::
resize(size_t newSize)
{
	size_t oldSize = size();
	m_layers.resize(newSize);
	for(size_t i = oldSize; i < newSize; ++i){
		m_layers[i] = make_sp(new layer_t);
	}
}

void RasterLayers::
invalidate_flat_cells(number minHeight)
{
//	since the top-layer is considered to be the terrain surface, we'll start
//	one layer below the top.
	if(size() <= 1)
		return;
	for(int lvl = (int)size() - 2; lvl >= 0; --lvl){
		Heightfield& curHF = layer(lvl);
		Field<number>& innerField = curHF.field();

	//	iterate over the cells of each heightfield and invalidate values
		for(int iy = 0; iy < (int)innerField.height(); ++iy){
			for(int ix = 0; ix < (int)innerField.width(); ++ix){
				number val = innerField.at(ix, iy);
				if(val != curHF.no_data_value()){
					bool gotDataVal = false;
					vector2 c = curHF.index_to_coordinate(ix, iy);
				//	compare the value against values from higher layers
					for(int ulvl = lvl + 1; ulvl < (int)size(); ++ulvl){
						Heightfield& uHF = layer(ulvl);
						number uval = uHF.interpolate(c);
						if(uval != uHF.no_data_value()){
							gotDataVal = true;
							if(val + minHeight > uval){
							//	the height in the current cell is too small
								innerField.at(ix, iy) = curHF.no_data_value();
							}
							break;
						}
					}

					if(!gotDataVal){
					//	The cell lies outside of the domain which is specified
					//	by the upmost field. The cell thus has to be invalidated
						innerField.at(ix, iy) = curHF.no_data_value();
					}
				}
			}
		}
	}
}


void  RasterLayers::
invalidate_small_lenses(number minArea)
{
	for(size_t lvl = 0; lvl < size(); ++lvl){
		Heightfield& hf = layer(lvl);
		number cellSize = hf.cell_size().x() * hf.cell_size().y();
		if(cellSize > 0){
			InvalidateSmallLenses(hf.field(), minArea / cellSize,
								  hf. no_data_value());
		}
	}
}

void RasterLayers::
remove_small_holes(number maxArea, number minHeight)
{
	using namespace std;

	const int numNbrs = 4;
	const int xadd[numNbrs] = {0, -1, 1, 0};
	const int yadd[numNbrs] = {-1, 0, 0, 1};

	Field<bool>	visited;
	vector<pair<int, int> > cells;

//	this field stores whether we already visited the given cell
	for(int lvl = (int)size() - 2; lvl >= 0; --lvl){
		Heightfield& hf = layer(lvl);
		Field<number>& field = hf.field();
		number noDataValue = hf.no_data_value();

		number cellSize = hf.cell_size().x() * hf.cell_size().y();
		if(cellSize <= 0)
			continue;

		size_t thresholdCellCount(maxArea / cellSize);


		visited.resize_no_copy(field.width(), field.height());
		visited.fill_all(false);

		const int fwidth = (int)field.width();
		const int fheight = (int)field.height();

		for(int outerIy = 0; outerIy < fheight; ++outerIy){
			for(int outerIx = 0; outerIx < fwidth; ++outerIx){
				if(visited.at(outerIx, outerIy)
				   || (field.at(outerIx, outerIy) != noDataValue))
				{
					continue;
				}

				cells.clear();
				cells.push_back(make_pair(outerIx, outerIy));
				size_t curCell = 0;
				while(curCell < cells.size()){
					int ix = cells[curCell].first;
					int iy = cells[curCell].second;

					for(size_t inbr = 0; inbr < numNbrs; ++inbr){
						int nx = ix + xadd[inbr];
						int ny = iy + yadd[inbr];
						if((nx >= 0 && nx < fwidth && ny >= 0 && ny < fheight)
						   &! visited.at(nx, ny))
						{
							visited.at(nx, ny) = true;
							if(field.at(nx, ny) == noDataValue){
								cells.push_back(make_pair(nx, ny));
							}
						}
					}
					++curCell;
				}

				if(cells.size() < thresholdCellCount){
					for(size_t i = 0; i < cells.size(); ++i){
						int ix = cells[i].first;
						int iy = cells[i].second;
						vector2 c = hf.index_to_coordinate(ix, iy);
						pair<int, number> result = trace_line_up(c, (size_t)lvl);
						if(result.first > lvl){
						//	we have to adjust not only the value in the current layer
						//	but also possibly values in lower levels, to avoid
						//	creating new thin layers
							number curVal = result.second;
							field.at(ix, iy) = curVal;
							for(int lowerLvl = lvl - 1; lowerLvl >= 0 ; --lowerLvl){
								Heightfield& lhf = layer(lowerLvl);
								Field<number>& lfield = lhf.field();
								pair<int, int> index = lhf.coordinate_to_index(c.x(), c.y());
								int lx = index.first;
								int ly = index.second;
								if(lx >= 0 && lx < (int)lfield.width()
									&& ly >= 0 && ly < (int)lfield.height())
								{
									number lval = lfield.at(lx, ly);
									if(lval != lhf.no_data_value()){
										if(curVal - lval < minHeight){
											curVal -= minHeight;
											lfield.at(lx, ly) = curVal;
										}
										else
											break;
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


void RasterLayers::
snap_cells_to_higher_layers(number minHeight)
{
	if(size() <= 1)
		return;

	for(int lvl = (int)size() - 2; lvl >= 0; --lvl){
		Heightfield& curHF = layer(lvl);
		Field<number>& curField = curHF.field();
		Heightfield& upperHF = layer(lvl + 1);

		for(int iy = 0; iy < (int)curField.height(); ++iy){
			for(int ix = 0; ix < (int)curField.width(); ++ix){
				number curVal = curField.at(ix, iy);
				vector2 c = curHF.index_to_coordinate(ix, iy);
				number upperVal = upperHF.interpolate(c);
				if(upperVal == upperHF.no_data_value())
					continue;

				if((curVal == curHF.no_data_value()) || (curVal > upperVal)
					|| (upperVal - curVal < minHeight))
				{
					curField.at(ix, iy) = upperVal;
				}
			}
		}
	}
}


void RasterLayers::
eliminate_invalid_cells()
{
	for(size_t lvl = 0; lvl < size(); ++lvl){
		Heightfield& hf = layer(lvl);
		EliminateInvalidCells(hf.field(), hf.no_data_value());
	}
}


void RasterLayers::
blur_layers(number alpha, size_t numIterations)
{
	for(size_t i = 0; i < size(); ++i){
		BlurField(layer(i).field(), alpha, numIterations, layer(i).no_data_value());
	}
}

std::pair<int, number> RasterLayers::
trace_line_down(const vector2& c, size_t firstLayer) const
{
	for(int i = (int)firstLayer; i >= 0; --i){
		number val = layer(i).interpolate(c);
		if(val != layer(i).no_data_value()){
			return make_pair(i, val);
		}
	}

	return make_pair<int, number>(-1, 0);
}

std::pair<int, number> RasterLayers::
trace_line_up(const vector2& c, size_t firstLayer) const
{
	for(size_t i = firstLayer; i < size(); ++i){
		number val = layer(i).interpolate(c);
		if(val != layer(i).no_data_value()){
			return make_pair(i, val);
		}
	}

	return make_pair<int, number>(-1, 0);
}


////////////////////////////////////////////////////////////////////////////////
void MeshLayerBoundaries(Grid& grid, const RasterLayers& layers,
						 Grid::VertexAttachmentAccessor<AVector3> aaPos,
						 ISubsetHandler* pSH)
{
	int defSubInd = -1;
	if(pSH)
		defSubInd = pSH->get_default_subset_index();

	for(size_t lvl = 0; lvl < layers.size(); ++lvl){
		if(pSH)
			pSH->set_default_subset_index((int)lvl);
		CreateGridFromFieldBoundary(grid, layers[lvl], aaPos);
	}

	if(pSH)
		pSH->set_default_subset_index(defSubInd);
}

////////////////////////////////////////////////////////////////////////////////
void MeshLayers(Grid& grid, const RasterLayers& layers,
				Grid::VertexAttachmentAccessor<AVector3> aaPos,
				ISubsetHandler* pSH)
{
	int defSubInd = -1;
	if(pSH)
		defSubInd = pSH->get_default_subset_index();

	for(size_t lvl = 0; lvl < layers.size(); ++lvl){
		if(pSH)
			pSH->set_default_subset_index((int)lvl);
		CreateGridFromField(grid, layers[lvl], aaPos);
	}

	if(pSH)
		pSH->set_default_subset_index(defSubInd);
}

struct ConnectedToOneMarkedVrt{
	ConnectedToOneMarkedVrt(Grid& grid) : m_grid(grid) {}
	bool operator() (Edge* e) const{
		return	(m_grid.is_marked(e->vertex(0)) || m_grid.is_marked(e->vertex(1)))
			&& !(m_grid.is_marked(e->vertex(0)) && m_grid.is_marked(e->vertex(1)));
	}
	Grid& m_grid;
};

void ExtrudeLayers(Grid& grid, const RasterLayers& layers,
				   Grid::VertexAttachmentAccessor<AVector3> aaPos,
				   ISubsetHandler& sh)
{
	UG_COND_THROW(layers.size() < 2, "At least 2 layers are required to perform extrusion!");

	grid.begin_marking();

	vector<Vertex*> curVrts;	// list of vertices that are considered for extrusion
	vector<Vertex*> tmpVrts;	// used to determine the set of vertices that can be extruded
	vector<Vertex*> smoothVrts;
	vector<Face*> curFaces;		// list of faces that are considered for extrusion
	vector<Face*> tmpFaces;		// used to determine the set of faces that can be extruded
	vector<number> vrtHeightVals;	// here we'll record height-values at which new vertices will be placed
	vector<number> volHeightVals;	// here we'll record height-values at the center of new volumes
	vector<int> volSubsetInds;
	vector<Volume*> newVols;
	queue<Volume*> volCandidates;

	Grid::edge_traits::callback cbIsMarked = IsMarked(grid);
	Grid::edge_traits::callback cbConnectedToOneMarkedVrt = ConnectedToOneMarkedVrt(grid);
	AssocElemIter<Vertex, Edge> assocVrtEdgeIterMarkedEdge(cbIsMarked);
	AssocElemIter<Vertex, Edge> assocVrtEdgeIterOneMarked(cbConnectedToOneMarkedVrt);
	AssocElemIter<Face, Edge> assocFaceEdgeIter;
	AssocElemIter<Volume, Edge> assocVolEdgeIter(cbConnectedToOneMarkedVrt);
	
	grid.reserve<Vertex>(grid.num_vertices() * layers.size());
	grid.reserve<Volume>(grid.num_faces() * layers.size() - 1);

//	this accessor is used during smoothing only
	ANumber aHeight;
	grid.attach_to_vertices(aHeight);
	Grid::VertexAttachmentAccessor<ANumber> aaHeight(grid, aHeight);

//	we have to determine the vertices that can be projected onto the top of the
//	given layers-stack. Only those will be used during extrusion
//	all considered vertices will be marked.
	const RasterLayers::layer_t& top = layers.top();
	const int topLayerInd = (int)layers.num_layers() - 1;
	for(VertexIterator i = grid.begin<Vertex>(); i != grid.end<Vertex>(); ++i){
		Vertex* v = *i;
		number val = top.interpolate(vector2(aaPos[v].x(), aaPos[v].y()));
		if(val != top.no_data_value()){
			aaPos[v].z() = val;
			curVrts.push_back(v);
			grid.mark(v);
			sh.assign_subset(v, topLayerInd);
		}
	}

	if(curVrts.size() < 3){
		UG_LOG("Not enough vertices lie in the region of the surface layer\n");
		return;
	}

//	all faces of the initial triangulation that are only connected to marked
//	vertices are considered for extrusion
	for(FaceIterator fi = grid.begin<Face>(); fi != grid.end<Face>(); ++fi){
		Face* f = *fi;
		bool allMarked = true;
		Face::ConstVertexArray vrts = f->vertices();
		const size_t numVrts = f->num_vertices();
		for(size_t i = 0; i < numVrts; ++i){
			if(!grid.is_marked(vrts[i])){
				allMarked = false;
				break;
			}
		}
		if(allMarked)
			curFaces.push_back(f);
	}

	if(curFaces.empty()){
		UG_LOG("Not enough faces lie in the region of the surface layer\n");
		return;
	}


	tmpVrts.reserve(curVrts.size());
	smoothVrts.reserve(curVrts.size());
	vrtHeightVals.reserve(curVrts.size());
	tmpFaces.reserve(curFaces.size());
	newVols.reserve(curFaces.size());
	volHeightVals.reserve(curFaces.size());
	volSubsetInds.reserve(curFaces.size());

	for(int ilayer = (int)layers.size() - 2; ilayer >= 0; --ilayer){

		tmpVrts.clear();
		vrtHeightVals.clear();
		tmpFaces.clear();
		newVols.clear();
		volHeightVals.clear();
		volSubsetInds.clear();
		grid.clear_marks();

	//	trace rays from the current vertices down through the layers until the
	//	next valid entry is found. If none is found, the vertex will be ignored
	//	from then on.
		for(size_t icur = 0; icur < curVrts.size(); ++icur){
			Vertex* v = curVrts[icur];
			vector2 c(aaPos[v].x(), aaPos[v].y());
			pair<int, number> val = layers.trace_line_down(c, ilayer);
			number height;

			if(val.first >= 0){
			//	if val.first == ilayer height will equal val.second. If not,
			//	a linear interpolation is performed, considering the height-val
			//	of the current vertex, the layer distance and the target value.
			//	This height-value will be corrected later on after extrusion
				number ia = 1. / ((number)ilayer - (number)val.first + 1.);
				height = (1. - ia) * aaPos[v].z() + ia * val.second;
				tmpVrts.push_back(v);
				vrtHeightVals.push_back(height);
				sh.assign_subset(v, val.first);
				grid.mark(v);
			}
		}

	//	now find the faces which connect those vertices
		for(size_t iface = 0; iface < curFaces.size(); ++iface){
			Face* f = curFaces[iface];
		//	trace a line from the x-y-center of the face downwards starting at
		//	the current layer to determine the layer in which the face has to
		//	be placed.
			vector3 center = CalculateCenter(f, aaPos);
			vector2 c(center.x(), center.y());
			pair<int, number> val = layers.trace_line_down(c, ilayer);
			if(val.first < 0)
				continue;

		//	now check whether all vertices are marked
			bool allMarked = true;
			Face::ConstVertexArray vrts = f->vertices();
			size_t numVrts = f->num_vertices();
			for(size_t i = 0; i < numVrts; ++i){
				if(!grid.is_marked(vrts[i])){
					allMarked = false;
					break;
				}
			}

			if(allMarked){
				pair<int, number> upVal = layers.trace_line_up(c, ilayer+1);
				if(upVal.first != -1){
					tmpFaces.push_back(f);
					volSubsetInds.push_back(val.first);
					volHeightVals.push_back(upVal.second - val.second);
				}
			}
		}

		if(tmpFaces.empty())
			break;

	//	swap containers and perform extrusion
		curVrts.swap(tmpVrts);
		curFaces.swap(tmpFaces);
		Extrude(grid, &curVrts, NULL, &curFaces, vector3(0, 0, 0), aaPos,
				EO_DEFAULT, &newVols);


	// assign pre-determined subsets
		for(size_t ivol = 0; ivol < newVols.size(); ++ivol)
			sh.assign_subset(newVols[ivol], volSubsetInds[ivol]);

	//	set the precalculated height of new vertices
		for(size_t ivrt = 0; ivrt < curVrts.size(); ++ivrt){
			aaPos[curVrts[ivrt]].z() = vrtHeightVals[ivrt];
		}


	//	finally adjust the height of new vertices which do not belong to the
	//	current layer
		grid.clear_marks();
		grid.mark(curVrts.begin(), curVrts.end());

	//	Consider all current-layer-volumes and apply a minHeight constraint to
	//	those edges which do not lie in the current layer.
	//	We'll also assign subset-indices of vertices connected to volumes of
	//	the current layer to the current layer index
		for(size_t ivol = 0; ivol < newVols.size(); ++ivol){
			Volume* vol = newVols[ivol];
			if(volSubsetInds[ivol] != ilayer)
				continue;

			const number shrinkConst = 0.2;
			const number minHeight = shrinkConst * volHeightVals[ivol];

			for(assocVolEdgeIter.reinit(grid, vol); assocVolEdgeIter.valid();
				++assocVolEdgeIter)
			{
				Edge* e = *assocVolEdgeIter;
				Vertex* from, *to;
				if(grid.is_marked(e->vertex(0))){
					from = e->vertex(1); to = e->vertex(0);
				}
				else{
					from = e->vertex(0); to = e->vertex(1);
				}

				if(sh.get_subset_index(to) != ilayer){
					aaPos[to].z() = max(aaPos[to].z(), aaPos[from].z() - minHeight);
					sh.assign_subset(to, ilayer);
				}
			}
		}

	//	prepare smoothing
	//	Push all new vertices which do not belong to the current layer to smoothVrts
	//	and the vertices directly above them to tmpVrts.
	//	Calculate height for each new vertex on the fly
		smoothVrts.clear();
		tmpVrts.clear();
		for(size_t ivrt = 0; ivrt < curVrts.size(); ++ivrt){
			Vertex* v = curVrts[ivrt];
			Vertex* conVrt = NULL;
			for(assocVrtEdgeIterOneMarked.reinit(grid, v);
				assocVrtEdgeIterOneMarked.valid(); ++assocVrtEdgeIterOneMarked)
			{
				conVrt = GetConnectedVertex(*assocVrtEdgeIterOneMarked, v);
				aaHeight[v] = aaPos[conVrt].z() - aaPos[v].z();
			}
			if((sh.get_subset_index(v) != ilayer) && conVrt){
				smoothVrts.push_back(v);
				tmpVrts.push_back(conVrt);
			}
		}


	//	to perform smoothing we need different marks again
		grid.clear_marks();

	//	we'll mark all smoothing vertices
		grid.mark(smoothVrts.begin(), smoothVrts.end());

	//	and we'll mark all those edges along which we'll smooth
		for(size_t iface = 0; iface < curFaces.size(); ++iface){
			for(assocFaceEdgeIter.reinit(grid, curFaces[iface]); assocFaceEdgeIter.valid();
				++assocFaceEdgeIter)
			{
				grid.mark(*assocFaceEdgeIter);
			}
		}

	//	during smoothing, the height of non-marked vertices will be weighted stronger
	//	than the height of marked ones
		const size_t numSmoothIterations = 100;
		const number markedWeight = 1.0;
		const number smoothAlpha = 0.5;
		for(size_t iSmoothIter = 0; iSmoothIter < numSmoothIterations; ++iSmoothIter){
			for(size_t ivrt = 0; ivrt < smoothVrts.size(); ++ivrt){
				Vertex* v = smoothVrts[ivrt];
				number totalNbrWeight = 0;
				number avHeight = 0;
				for(assocVrtEdgeIterMarkedEdge.reinit(grid, v); assocVrtEdgeIterMarkedEdge.valid();
					++assocVrtEdgeIterMarkedEdge)
				{
					Vertex* conVrt = GetConnectedVertex(*assocVrtEdgeIterMarkedEdge, v);
					if(grid.is_marked(conVrt)){
						avHeight += markedWeight * aaHeight[conVrt];
						totalNbrWeight += markedWeight;
					}
					else{
						avHeight += aaHeight[conVrt];
						totalNbrWeight += 1;
					}
				}

				UG_COND_THROW(totalNbrWeight == 0, "No neighbors found");
				avHeight /= totalNbrWeight;
				
				aaHeight[v] = (1. - smoothAlpha) * aaHeight[v] + smoothAlpha * avHeight;
			}
		}

	//	move vertices upwards only, to avoid invalid volumes in lower layers
		for(size_t ivrt = 0; ivrt < smoothVrts.size(); ++ivrt){
			Vertex* v = smoothVrts[ivrt];
			aaPos[v].z() = max(aaPos[v].z(), aaPos[tmpVrts[ivrt]].z() - aaHeight[v]);
		}
	}

	grid.end_marking();
	grid.detach_from_vertices(aHeight);
}

}//	end of namespace
