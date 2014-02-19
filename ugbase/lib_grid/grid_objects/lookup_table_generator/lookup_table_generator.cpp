//	created by Sebastian Reiter
//	s.b.reiter@gmail.com

//	NOTE:	THIS FILE IS NOT PART OF THE UG4-LIBRARY AND SHOULD NOT BE COMPILED ALONG
//			WITH OTHER SOURCES. IT SHOULD ONLY BE COMPILED SEPARATELY IN ORDER
//			TO RE-GENERATE ADDITIONAL LOOKUP-TABLES AFTER CHANGES TO THE ORIGINAL
//			LOOKUP-TABLES in tetrahedron_rules, pyramid_rules, prism_rules or
//			hexahedron_rules HAVE BEEN PERFORMED.
//
//	Compilation hint: Use option -I .../ug4/trunk/ugbase/

//	creates multi-dimensional lookup-tables, which allow to find an edge or
//	a face of an element using their local vertex-indices

#include <iostream>
#include "../tetrahedron_rules.h"
#include "../pyramid_rules.h"
#include "../prism_rules.h"
#include "../hexahedron_rules.h"

using namespace std;
using namespace ug;

void PrintFaceEdgeInds(int numEdges, int numFaces,
						const int edgeVrtInds[][2],
						const int faceVrtInds[][4])
{
	cout << "const int FACE_EDGE_INDS[" << numFaces<< "][" << 4 << "] = {";
	
	for(int i_face = 0; i_face < numFaces; ++i_face){
		const int* f = faceVrtInds[i_face];
		cout << "{";
		for(int i_vrt = 0; i_vrt < 4; ++i_vrt){
		//	since triangles have a -1 in their 4-th component, we have to
		//	make sure to get valid indices.
		//	Note - if i0 == -1, then we won't find a matching edge and the
		//	value -1 is thus written as corresponding edge-index.
			int i0 = f[i_vrt];
				
			int i1 = f[(i_vrt + 1) % 4];
			if(i1 == -1)
				i1 = f[(i_vrt + 2) % 4];
				
		//	now find the edge containing both i0 and i1
			int edgeInd = -1;
			for(int i_edge = 0; i_edge < numEdges; ++i_edge){
				int numMatches = 0;
				const int* e = edgeVrtInds[i_edge];
				if(e[0] == i0 || e[1] == i0)
					++numMatches;
					
				if(e[0] == i1 || e[1] == i1)
					++numMatches;
				
				if(numMatches == 2){
					edgeInd = i_edge;
					break;
				}
			}
			
			cout << edgeInd;
			if(i_vrt < 3)
				cout << ", ";
		}
		cout << "}";
		if(i_face < numFaces - 1)
			cout << ", ";
	}
	
	cout << "};" << endl;
}

void PrintEdgeFromVrts(int numVrts, int numEdges, const int edgeVrtInds[][2])
{
//	print a two-dimensional lookup table, where the edge index for each pair of
//	vertices can be found.
	cout << "const int EDGE_FROM_VRTS[" << numVrts << "][" << numVrts << "] = {";
	
	for(int i_v1 = 0; i_v1 < numVrts; ++i_v1){
		cout << "{";
		for(int i_v2 = 0; i_v2 < numVrts; ++i_v2){
		//	now iterate through the edges and check whether one connects i_v1 and i_v2.
			int eind = -1;
			for(int i = 0; i < numEdges; ++i){
				const int* edge = edgeVrtInds[i];
				if((edge[0] == i_v1 && edge[1] == i_v2)
					|| (edge[0] == i_v2 && edge[1] == i_v1))
				{
					eind = i;
					break;
				}
			}
			
			cout << eind;
			if(i_v2 < numVrts - 1)
				cout << ", ";
		}
		cout << "}";
		if(i_v1 < numVrts - 1)
			cout << ", ";
	}
	cout << "};" << endl;
}

void PrintFaceFromVrts(int numVrts, int numFaces, const int faceVrtInds[][4])
{
//	print a three-dimensional lookup table, where the face index for each triple of
//	vertices can be found. Note that a triple is enough to identify a quadrilateral, too.
	cout << "const int FACE_FROM_VRTS[" << numVrts << "]["
		 << numVrts << "][" << numVrts << "] = {";
	
	for(int i_v1 = 0; i_v1 < numVrts; ++i_v1){
		cout << "{";
		for(int i_v2 = 0; i_v2 < numVrts; ++i_v2){
			cout << "{";
			for(int i_v3 = 0; i_v3 < numVrts; ++i_v3){
			//	now iterate through the faces and check whether one connects i_v1, i_v2 and i_v3.
				int find = -1;
				for(int i_face = 0; i_face < numFaces; ++i_face){
					const int* face = faceVrtInds[i_face];
					int matchCount = 0;
					for(int i = 0; i < 4; ++i){
						if(face[i] == i_v1 || face[i] == i_v2 || face[i] == i_v3)
							++matchCount;
					}
					if(matchCount >= 3){
						find = i_face;
						break;
					}
				}
				cout << find;
				if(i_v3 < numVrts - 1)
					cout << ", ";	
			}
			cout << "}";
			if(i_v2 < numVrts - 1)
				cout << ", ";
		}
		cout << "}";
		if(i_v1 < numVrts - 1)
			cout << ", ";
	}
	cout << "};" << endl;
}

void PrintFaceFromEdges(int numEdges, int numFaces,
						const int edgeVrtInds[][2],
						const int faceVrtInds[][4])
{
	cout << "const int FACE_FROM_EDGES[][" << numEdges << "] = {";
	for(int i_e0 = 0; i_e0 < numEdges; ++i_e0){
		cout << "{";
		const int* e0 = edgeVrtInds[i_e0];
		for(int i_e1 = 0; i_e1 < numEdges; ++i_e1){
			const int* e1 = edgeVrtInds[i_e1];
			
		//	search for the face which contains both e0 and e1
			int faceInd = -1;
			for(int i_face = 0; i_face < numFaces; ++i_face){
				const int* f = faceVrtInds[i_face];
				int match0 = 0;
				int match1 = 0;
				for(int i = 0; i < 4; ++i){
					if(f[i] == e0[0] || f[i] == e0[1])
						++match0;
					if(f[i] == e1[0] || f[i] == e1[1])
						++match1;
				}
				if(match0 == 2 && match1 == 2){
					faceInd = i_face;
					break;
				}
			}
			
			cout << faceInd;
			if(i_e1 < numEdges - 1)
				cout << ", ";
		}
		cout << "}";
		if(i_e0 < numEdges - 1)
			cout << ", ";
	}
	cout << "};" << endl;
}

void PrintFaceContainsEdgeLists(int numEdges, int numFaces,
								const int edgeVrtInds[][2],
								const int faceVrtInds[][4])
{
	cout << "const int FACE_CONTAINS_EDGE[][" << numEdges << "] = {";
	
	for(int i_face = 0; i_face < numFaces; ++i_face){
		cout << "{";
		const int* f = faceVrtInds[i_face];
		for(int i_edge = 0; i_edge < numEdges; ++i_edge){
			const int* e = edgeVrtInds[i_edge];
		//	check whether e is contained in f
			int numMatches = 0;
			for(int i = 0; i < 4; ++i){
				if(f[i] == e[0] || f[i] == e[1])
					++numMatches;
			}
			
			if(numMatches == 2)
				cout << "1";
			else
				cout << "0";
				
			if(i_edge < numEdges - 1)
				cout << ", ";
		}
		cout << "}";
		if(i_face < numFaces - 1)
			cout << ", ";
	}
	cout << "};" << endl;
}

int main(int argc, char* argv[])
{

//	tetrahedron lists
	cout << "\ntetrahedron lists:\n";
	{
		using namespace ug::tet_rules;
		PrintFaceEdgeInds(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintFaceContainsEdgeLists(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintEdgeFromVrts(NUM_VERTICES, NUM_EDGES, EDGE_VRT_INDS);
		PrintFaceFromVrts(NUM_VERTICES, NUM_FACES, FACE_VRT_INDS);
		PrintFaceFromEdges(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
	}
	
//	pyramid lists
	cout << "\npyramid lists:\n";
	{
		using namespace ug::pyra_rules;
		PrintFaceEdgeInds(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintFaceContainsEdgeLists(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintEdgeFromVrts(NUM_VERTICES, NUM_EDGES, EDGE_VRT_INDS);
		PrintFaceFromVrts(NUM_VERTICES, NUM_FACES, FACE_VRT_INDS);
		PrintFaceFromEdges(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
	}
	
//	prism lists
	cout << "\nprism lists:\n";
	{
		using namespace ug::prism_rules;
		PrintFaceEdgeInds(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintFaceContainsEdgeLists(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintEdgeFromVrts(NUM_VERTICES, NUM_EDGES, EDGE_VRT_INDS);
		PrintFaceFromVrts(NUM_VERTICES, NUM_FACES, FACE_VRT_INDS);
		PrintFaceFromEdges(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
	}
	
//	hexahedron lists
	cout << "\nhexahedron lists:\n";
	{
		using namespace ug::hex_rules;
		PrintFaceEdgeInds(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintFaceContainsEdgeLists(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
		PrintEdgeFromVrts(NUM_VERTICES, NUM_EDGES, EDGE_VRT_INDS);
		PrintFaceFromVrts(NUM_VERTICES, NUM_FACES, FACE_VRT_INDS);
		PrintFaceFromEdges(NUM_EDGES, NUM_FACES, EDGE_VRT_INDS, FACE_VRT_INDS);
	}
	return 0;
}
