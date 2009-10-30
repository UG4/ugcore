//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d26

#define SHINY_PROFILER TRUE

#include <iostream>
#include "lib_grid/lib_grid.h"
#include "common/profiler/profiler.h"
#include "lib_grid/advanced/advanced.h"

using namespace std;
using namespace ug;

//	typedefs
typedef Hash<EdgeBase*, EdgeVertices*>	EdgeHash;

//	temporary counters
int numChecks;

////////////////////////////////////////////////////////////////////////
//	FillHash
void FillHash(EdgeHash& hash, Grid& grid)
{
	PROFILE_FUNC();

	for(EdgeBaseIterator iter = grid.edges_begin();
		iter != grid.edges_end(); ++iter)
	{
		hash.add(*iter, *iter);
	}
}

///	vDist can be obtained by a call to Hash::get_distribution
void PrintHashStatistics(std::vector<int>& vDist)
{
	cout << "Hash statistics:\n";
	
	cout << "number of lists with #row elements:\n#elems:\t#lists\n";
//	output statistics
	for(size_t i = 0; i < vDist.size(); ++i)
	{
		cout << "  " << i << ":\t" << vDist[i] << endl;
	}

//	num elements
	int numElems = 0;
	for(size_t i = 1; i < vDist.size(); ++i)
		numElems += vDist[i];

//	average number of lookups for an element
	cout << "average number of lookups per element: ";
	if(numElems > 0)
	{
		float averageNum = 0;
		for(size_t i = 1; i < vDist.size(); ++i)
			averageNum += (i * vDist[i]);
		averageNum /= numElems;
		cout << averageNum << endl;
	}
	else
		cout << "0" << endl;

	cout << endl;
}

EdgeBase* GetEdgeOld(Grid& grid, Face* f, int index)
{
	return grid.get_edge(f, index);
}

EdgeBase* GetEdgeFromHashVec(EdgeHash& hash, Face* f, int index)
{
	EdgeDescriptor ed;
	f->edge(index, ed);

	vector<EdgeBase*> vEdges;
	hash.get_entries(vEdges, &ed);

	for(size_t i = 0; i < vEdges.size(); ++i)
	{
		bool bMatches;
		bMatches = CompareVertices(vEdges[i], &ed);
		if(bMatches)
			return vEdges[i];
	}

	return NULL;
}

EdgeBase* GetEdgeFromHashIter(EdgeHash& hash, Face* f, int index)
{
	EdgeDescriptor ed;
	f->edge(index, ed);

	EdgeHash::Iterator iter;
	EdgeHash::Iterator iterEnd;

	hash.get_iterators(iter, iterEnd, &ed);
	
	for(; iter != iterEnd; ++iter)
	{
		numChecks++;
		bool bMatches;
		bMatches = CompareVertices(*iter, &ed);
		if(bMatches)
			return *iter;
	}

	return NULL;
}

bool PerformGetElemCompare(const char* filename)
{
	PROFILE_FUNC();
	cout << "PerformGetElemCompare\n";
	Grid grid(GRIDOPT_FULL_INTERCONNECTION);

	cout << "options: " << grid.get_options() << endl;

	cout << "loading...\n";
	PROFILE_BEGIN(loading);
	bool loadOK;
	PROFILE_CODE(loadOK = LoadGridFromFile(grid, filename));

	PROFILE_END();

	if(!loadOK)
	{
		cout << "...failed\n";
		return false;
	}
	cout << "...done\n";


//	iterate through all faces and collect check edges
	for(FaceIterator iter = grid.faces_begin(); iter != grid.faces_end(); ++iter)
	{
		Face* f = *iter;
		
		for(int i = 0; i < f->num_edges(); ++i)
		{
			PROFILE_CODE(grid.get_edge(f, i));
		}
	}

	for(VolumeIterator iter = grid.volumes_begin(); iter != grid.volumes_end(); ++iter)
	{
		Volume* v = *iter;
		
		for(int i = 0; i < v->num_edges(); ++i)
		{
			PROFILE_CODE(grid.get_edge(v, i));
		}
	}

	for(VolumeIterator iter = grid.volumes_begin(); iter != grid.volumes_end(); ++iter)
	{
		Volume* v = *iter;
		
		for(int i = 0; i < v->num_faces(); ++i)
		{
			PROFILE_CODE(grid.get_face(v, i));
		}
	}

	cout << "num vertices:\t" << grid.num_vertices() << endl;
	cout << "num edges:\t" << grid.num_edges() << endl;
	cout << "num faces:\t" << grid.num_faces() << endl;
	cout << "num volumes:\t" << grid.num_volumes() << endl;

	return true;
}

////////////////////////////////////////////////////////////////////////
//	main
int main(int argc, char* argv[])
{
//	load a geometry and test the hashing of edges and faces.
	if(argc != 2)
	{
		cout << "please specify a 2d geometry (.txt or .obj).\n";
		return 0;
	}
	
	PerformGetElemCompare(argv[1]);
	
/*
//	begin profiling for the main section.
//	main_section is an arbitrary name.
	cout << "loading grid " << argv[1] << " ..." << endl;
	PROFILE_BEGIN(load_grid);

//	create a grid.
	Grid grid;
	SubsetHandler sh(grid);

	grid.enable_options(VRTOPT_STORE_ASSOCIATED_EDGES |
						FACEOPT_AUTOGENERATE_EDGES |
						FACEOPT_STORE_ASSOCIATED_EDGES |
						VOLOPT_AUTOGENERATE_FACES);

	bool loadOK = LoadGridFromFile(grid, argv[1], sh);
	
	PROFILE_END();

	if(loadOK)
	{
		cout << "...done.\n";

		int num1 = grid.num_vertices();
		int num2 = 2*num1;

		int hashSize = num2;
		cout << "hash size:" << hashSize << endl;
		
	//	the edge hash
		EdgeHash edgeHash(hashSize);

		FillHash(edgeHash, grid);

		vector<int> vDist;
		edgeHash.get_distribution(vDist);
		PrintHashStatistics(vDist);

		cout << "performing edge-access check" << endl;
	//	iterate through all faces and make sure that
	//	both edge-access functions return the same value.
		PROFILE_BEGIN(get_edge_old);
		for(FaceIterator iter = grid.faces_begin();
			iter != grid.faces_end(); ++iter)
		{
			Face* f = *iter;
			for(uint i = 0; i < f->num_edges(); ++i)
				GetEdgeOld(grid, f, i);
		}
		PROFILE_END();
		PROFILE_BEGIN(get_edge_vec);
		for(FaceIterator iter = grid.faces_begin();
			iter != grid.faces_end(); ++iter)
		{
			Face* f = *iter;
			for(uint i = 0; i < f->num_edges(); ++i)
				GetEdgeFromHashVec(edgeHash, f, i);
		}
		PROFILE_END();
		
		numChecks = 0;
		int counter = 0;
		PROFILE_BEGIN(get_edge_iter);
		for(FaceIterator iter = grid.faces_begin();
			iter != grid.faces_end(); ++iter)
		{
			Face* f = *iter;
			for(uint i = 0; i < f->num_edges(); ++i, ++counter)
				GetEdgeFromHashIter(edgeHash, f, i);
		}
		PROFILE_END();
		cout << "average number of checks: ";
		cout << float(numChecks) / float(counter) << endl;
		
	//	check whether all methods return the same edge
		for(FaceIterator iter = grid.faces_begin();
			iter != grid.faces_end(); ++iter)
		{
			Face* f = *iter;
			for(uint i = 0; i < f->num_edges(); ++i)
			{
				EdgeBase* e = GetEdgeOld(grid, f, i);
				if(e != GetEdgeFromHashVec(edgeHash, f, i))
					cout << "edge-vec returns bad edge!" << endl;
				if(e != GetEdgeFromHashIter(edgeHash, f, i))
					cout << "edge-iter returns bad edge!" << endl;
			}
		}
		
	}
	else
	{
		cout << "...failed.\n";
	}
*/

	cout << "done\n\n";

//	call this for output.
	PROFILER_UPDATE();
	PROFILER_OUTPUT();

	return 0;
}

