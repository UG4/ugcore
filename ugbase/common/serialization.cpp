// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 24.03.2011 (m,d,y)
 
#include "serialization.h"

namespace ug{

void Serialize(std::ostream& buf, const std::string& str)
{
	size_t len = str.length();

	buf.write((char*)&len, sizeof(size_t));
	if(len > 0)
		buf.write(str.c_str(), sizeof(char) * len);
}

///	deserializes data from a binary stream into a string
void Deserialize(std::istream& buf, std::string& str)
{
//	the buffers allow us to read small strings fast.
//	for bigger ones we have to temporarily reserve memory.
	char staticBuf[64];
	char* flexBuf = NULL;
	char* tBuf = staticBuf;

	size_t len = 0;
	buf.read((char*)&len, sizeof(size_t));

//	check whether we have to allocate memory
//	don't forget that we have to append a zero at the end
	if(len >= 63){
		flexBuf = new char[len + 1];
		tBuf = flexBuf;
	}

	if(len > 0)
		buf.read(tBuf, sizeof(char) * len);
	tBuf[len] = 0;

//	assign data to the out-string
	str = tBuf;

//	clean up
	if(flexBuf)
		delete[] flexBuf;
}

}// end of namespace
