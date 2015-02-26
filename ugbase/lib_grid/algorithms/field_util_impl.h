// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_field_util_impl
#define __H__UG_field_util_impl

#include <algorithm>

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

}//	end of namespace

#endif	//__H__UG_field_util_impl
