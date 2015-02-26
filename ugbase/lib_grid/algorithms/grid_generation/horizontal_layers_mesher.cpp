// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "horizontal_layers_mesher.h"
#include "lib_grid/file_io/file_io_asc.h"

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
		Field<number>& innerField = curHF.field;

	//	iterate over the cells of each heightfield and invalidate values
		for(int iy = 0; iy < (int)innerField.height(); ++iy){
			for(int ix = 0; ix < (int)innerField.width(); ++ix){
				number val = innerField.at(ix, iy);
				if(val != curHF.noDataValue){
					bool gotDataVal = false;
					vector2 c = curHF.index_to_coordinate(ix, iy);
				//	compare the value against values from higher layers
					for(int ulvl = lvl + 1; ulvl < (int)size(); ++ulvl){
						Heightfield& uHF = layer(ulvl);
						number uval = uHF.interpolate(c);
						if(uval != uHF.noDataValue){
							gotDataVal = true;
							if(val + minHeight > uval){
							//	the height in the current cell is too small
								innerField.at(ix, iy) = curHF.noDataValue;
							}
							break;
						}
					}

					if(!gotDataVal){
					//	The cell lies outside of the domain which is specified
					//	by the upmost field. The cell thus has to be invalidated
						innerField.at(ix, iy) = curHF.noDataValue;
					}
				}
			}
		}
	}
}

void RasterLayers::
blur_layers(number alpha, size_t numIterations)
{
	for(size_t i = 0; i < size(); ++i){
		BlurField(layer(i).field, alpha, numIterations, layer(i).noDataValue);
	}
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
}//	end of namespace
