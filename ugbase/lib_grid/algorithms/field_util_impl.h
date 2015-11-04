#ifndef __H__UG_field_util_impl
#define __H__UG_field_util_impl

#include <algorithm>
#include <deque>
#include <queue>

namespace ug{

template <class T>
void BlurField(Field<T>& field, number alpha, size_t numIterations, const T& noDataValue)
{
	using namespace std;
	for(size_t mainIter = 0; mainIter < numIterations; ++mainIter){
		for(int iy = 0; iy < (int)field.height(); ++iy){
			for(int ix = 0; ix < (int)field.width(); ++ix){
				if(field.at(ix, iy) != noDataValue){
					T val = 0;
					number num = 0;
					for(int iny = max<int>(iy - 1, 0); iny < min<int>(iy + 2, (int)field.height()); ++iny){
						for(int inx = max<int>(ix - 1, 0); inx < min<int>(ix + 2, (int)field.width()); ++inx){
							if(!(inx == 0 && iny == 0) && (field.at(inx, iny) != noDataValue)){
								val += field.at(inx, iny);
								++num;
							}
						}
					}

					if(num > 0){
						val *= alpha / num;
						field.at(ix, iy) *= (1.-alpha);
						field.at(ix, iy) += val;
					}
				}
			}
		}
	}
}	

namespace fieldutil{
template <class T>
struct Cell{
	Cell(int _x, int _y, T _val) : x(_x), y(_y), value(_val) {}
	int x;
	int y;
	T value;
};
}

template <class T>
bool EliminateInvalidCells(Field<T>& field, const T& noDataValue)
{
	using namespace std;
	typedef fieldutil::Cell<T>	Cell;

	deque<Cell>	cells;
	number inProgressValue = -noDataValue;

	const int numNbrs = 8;
	const int xadd[numNbrs] = {-1, 0, 1, -1, 1, -1, 0, 1};
	const int yadd[numNbrs] = {-1, -1, -1, 0, 0, 1, 1, 1};

	const int maxNumSteps = 4;
	const int minNumValidNbrsInStep[maxNumSteps] = {4, 3, 2, 1};

//	initially count the number of invalid cells
	size_t numInvalidCells = 0;
	for(int iy = 0; iy < (int)field.height(); ++iy){
		for(int ix = 0; ix < (int)field.width(); ++ix){
			if(field.at(ix, iy) == noDataValue)
				++numInvalidCells;
		}
	}

//	find the initial cells which contain no-data-values and which are neighbors
//	of valid cells
//	We do this in several steps to better smear out the values and to avoid sharp features
//	in the smeared out regions
	for(int istep = 0; (istep < maxNumSteps) && (numInvalidCells > 0); ++istep){
		const int minNumValidNbrs = minNumValidNbrsInStep[istep];

		for(int iy = 0; iy < (int)field.height(); ++iy){
			for(int ix = 0; ix < (int)field.width(); ++ix){
				if(field.at(ix, iy) != noDataValue)
					continue;

				int numValidNbrs = 0;
				for(int i = 0; i < numNbrs; ++i){
					const int inx = ix + xadd[i];
					const int iny = iy + yadd[i];

					if((inx >= 0) && (inx < (int)field.width())
						&& (iny >= 0) && (iny < (int)field.height())
						&& (field.at(inx, iny) != noDataValue)
						&& (field.at(inx, iny) != inProgressValue))
					{
						++numValidNbrs;
					}
				}

				if(numValidNbrs >= minNumValidNbrs){
					cells.push_back(Cell(ix, iy, inProgressValue));
					field.at(ix, iy) = inProgressValue;
					break;
				}
			}
		}

		while(!cells.empty()){
		// iterate over all entries in the queue and calculate their correct values
			for(typename deque<Cell>::iterator cellIter = cells.begin(); cellIter != cells.end(); ++cellIter)
			{
				Cell& cell = *cellIter;
				const int ix = cell.x;
				const int iy = cell.y;

			//	get average value of valid neighbor cells
				T avVal = 0;
				number numValidNbrs = 0;
				for(int i = 0; i < numNbrs; ++i){
					const int inx = ix + xadd[i];
					const int iny = iy + yadd[i];

					if((inx >= 0) && (inx < (int)field.width())
						&& (iny >= 0) && (iny < (int)field.height()))
					{
						if((field.at(inx, iny) != noDataValue)
						   && (field.at(inx, iny) != inProgressValue))
						{
								avVal += field.at(inx, iny);
								++numValidNbrs;
						}
					}
				}

				UG_COND_THROW(numValidNbrs < (number)minNumValidNbrs, "Implementation error!");
				avVal *= 1. / numValidNbrs;
				cell.value = avVal;
				--numInvalidCells;
			}

		//	copy values to the field and collect new candidates
			while(!(cells.empty() || (cells.front().value == inProgressValue)))
			{
				const int ix = cells.front().x;
				const int iy = cells.front().y;
				field.at(ix, iy) = cells.front().value;
				cells.pop_front();

				for(int i = 0; i < numNbrs; ++i){
					const int inx = ix + xadd[i];
					const int iny = iy + yadd[i];

					if((inx >= 0) && (inx < (int)field.width())
						&& (iny >= 0) && (iny < (int)field.height()))
					{
						if(field.at(inx, iny) == noDataValue){
						//	the nbr-cell is a possible new candidate.
						//	count the number of valid nbrs-of that cell and
						//	add it to the queue if there are enough.
							int numValidNbrsOfNbr = 0;
							for(int j = 0; j < numNbrs; ++j){
								const int inxNbr = inx + xadd[j];
								const int inyNbr = iny + yadd[j];

								if((inxNbr >= 0) && (inxNbr < (int)field.width())
									&& (inyNbr >= 0) && (inyNbr < (int)field.height())
									&& (field.at(inxNbr, inyNbr) != noDataValue)
									&& (field.at(inxNbr, inyNbr) != inProgressValue))
								{
									++numValidNbrsOfNbr;
								}
							}
							if(numValidNbrsOfNbr >= minNumValidNbrs){
								cells.push_back(Cell(inx, iny, inProgressValue));
								field.at(inx, iny) = inProgressValue;
							}
						}
					}
				}
			}
		}
	}

	return numInvalidCells == 0;
}


template <class T>
void InvalidateSmallLenses(Field<T>& field, size_t thresholdCellCount,
						   const T& noDataValue)
{
	using namespace std;

	// const int numNbrs = 8;
	// const int xadd[numNbrs] = {-1, 0, 1, -1, 1, -1, 0, 1};
	// const int yadd[numNbrs] = {-1, -1, -1, 0, 0, 1, 1, 1};
	const size_t numNbrs = 4;
	const int xadd[numNbrs] = {0, -1, 1, 0};
	const int yadd[numNbrs] = {-1, 0, 0, 1};

//	this field stores whether we already visited the given cell
	Field<bool>	visited(field.width(), field.height(), false);
	vector<pair<int, int> > cells;

	const int fwidth = (int)field.width();
	const int fheight = (int)field.height();

	for(int outerIy = 0; outerIy < fheight; ++outerIy){
		for(int outerIx = 0; outerIx < fwidth; ++outerIx){
			if(visited.at(outerIx, outerIy)
			   || (field.at(outerIx, outerIy) == noDataValue))
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
					   && !visited.at(nx, ny))
					{
						visited.at(nx, ny) = true;
						if(field.at(nx, ny) != noDataValue){
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
					field.at(ix, iy) = noDataValue;
				}
			}
		}
	}
}

}//	end of namespace

#endif	//__H__UG_field_util_impl
