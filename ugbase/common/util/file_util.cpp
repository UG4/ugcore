/**
 * \file file_util.cpp
 * \date 2012-05-15
 * \brief Implementation of OS independent file utility functions
 */

#include "common/util/file_util.h"
#include <vector>


using namespace std;

namespace ug
{

UG_API bool FileExists( const char *filename )
{
  PROFILE_FUNC();

  ifstream in( filename );
  if( in.good() ) {
    in.close();
    return true;
  }
  return false;
}

UG_API size_t FileSize( const char *filename )
{
  PROFILE_FUNC();

  if( !FileExists( filename ) ) {
    UG_THROW( "The file " << filename << " could not be found." );
  }

  ifstream ifs;
  ifs.open( filename, ios_base::in );

  if( !ifs.good() ) {
    UG_THROW( "The file " << filename << " could not be opened." );
  }
  ifs.seekg( 0, ios_base::end );
  size_t length = ifs.tellg();
  ifs.close();
  return length;
}

UG_API bool FileCompare( const char *file1, const char *file2 )
{
  PROFILE_FUNC();

  // Make sure, both files do exist
  if( !FileExists( file1 ) || !FileExists( file2 ) ) {
    UG_THROW( "One or both files could not be found:" << file1 << ", " << file2 );
  }

  // If both files are one and the same, we don't need to compare them any more
  if( strcmp( file1, file2 ) == 0 ) {
    return true;
  }

  // If the files have different size, we know that they are different and can stop
  size_t size1 = FileSize( file1 ), size2 = FileSize( file2 );
  if( size1 != size2 ) {
    return false;
  }
  // If both files have zero size, we can stop as well
  if( size1 == 0 && size2 == 0 ) {
    return true;
  }

  // Open the files
  ifstream f1, f2;
  f1.open( file1, ios_base::in );
  f2.open( file2, ios_base::in );

  // This should be obsolete, as it is checked by FileExists, but we want to be
  // sure.
  if( !f1.good() || !f2.good() ) {
    UG_THROW( "One or both files could not be opened:" << file1 << ", " << file2 );
  }

  string line1 = "", line2 = "";
  bool diff = false;
  while( !f1.eof() || !f2.eof() ) {
    getline( f1, line1 );
    getline( f2, line2 );
    if( ( line1 != line2 ) ||
        ( f1.eof() && !f2.eof() ) ||
        ( !f1.eof() && f2.eof() ) ) {
      diff = true;
      break;
    }
  }

  f1.close();
  f2.close();

  return !diff;
}


bool ReadFile(const char* filename, vector<char> &file, bool bText)
{
	PROFILE_FUNC();
	FILE *f = fopen(filename, bText ? "r" : "rb");
	if(f==NULL)	return false;
	fseek(f, 0, SEEK_END);
	long filesize=ftell(f);
	fseek(f, 0, SEEK_SET);

	long actualFilesize=filesize;
	if(bText) actualFilesize++;
	file.resize(actualFilesize);

	fread(&file[0], 1, filesize, f);
	fclose(f);
	if(bText) file[filesize]=0x00;
	return true;
}

string MakeTmpFile(const string &filename, const string &extension, bool &bSuccess)
{
	bSuccess = true;
	const char *name = (filename+extension).c_str();
	if(!FileExists(name)) return name;
	for(int i=0;i<999999; i++)
	{
		stringstream ss;
		ss << filename << i << extension;
		name = ss.str().c_str();
		if(!FileExists(name)) return name;
	}	
	bSuccess = false;
	return "";
}

} // namespace ug

// EOF
