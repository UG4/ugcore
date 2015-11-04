#include <cassert>
#include "common/math/ugmath.h"
#include "simple_grid.h"

using namespace std;

namespace ug
{
////////////////////////////////////////////////////////////////////////
void SimpleGrid::clear()
{
	vertices.clear();
	vertexNormals.clear();
	triangles.clear();
	triangleNormals.clear();
}

void PrintSimpleGrid(SimpleGrid& sg)
{
	cout << "simple-grid:" << endl;
	cout << "num vrts: " << sg.vertices.size() << endl;
	cout << "num tris: " << sg.triangles.size() << endl;

	cout << "vertices: ";
	for(size_t i = 0; i < sg.vertices.size(); ++i)
		cout << sg.vertices[i] << ", ";
	cout << endl;

	cout << "normals: ";
	for(size_t i = 0; i < sg.vertexNormals.size(); ++i)
		cout << sg.vertexNormals[i] << ", ";
	cout << endl;

	cout << "triangles: ";
	for(size_t i = 0; i < sg.triangles.size(); ++i)
		cout << sg.triangles[i] << ", ";
	cout << endl;

	cout << "triangle normals: ";
	for(size_t i = 0; i < sg.triangleNormals.size(); ++i)
		cout << sg.triangleNormals[i] << ", ";
	cout << endl;
}

////////////////////////////////////////////////////////////////////////
void CalculateTriangleNormal(SimpleGrid& sg, int triIndex)
{
	int i = triIndex * 3;
	assert(((int)sg.triangles.size() > i + 2) && ((int)sg.triangleNormals.size() > triIndex) && "bad triangle index.");

	CalculateTriangleNormal(sg.triangleNormals[triIndex],
							sg.vertices[sg.triangles[i]],
							sg.vertices[sg.triangles[i+1]],
							sg.vertices[sg.triangles[i+2]]);
}

////////////////////////////////////////////////////////////////////////
void CalculateTriangleNormals(SimpleGrid& sg)
{
	size_t numTris = sg.triangles.size() / 3;
	sg.triangleNormals.resize(numTris);
	
	for(size_t i = 0; i < numTris; ++i)
		CalculateTriangleNormal(sg, i);
}

////////////////////////////////////////////////////////////////////////
number GeometricApproximationDegree(SimpleGrid& sg, int triIndex)
{
	size_t i = triIndex * 3;

	return GeometricApproximationDegree(sg.vertexNormals[sg.triangles[i]],
										sg.vertexNormals[sg.triangles[i+1]],
										sg.vertexNormals[sg.triangles[i+2]],
										sg.triangleNormals[triIndex]);
}

////////////////////////////////////////////////////////////////////////
number GeometricApproximationDegree(SimpleGrid& sg)
{
	if(sg.triangles.size() < 3)
		return 0;

	number q = GeometricApproximationDegree(sg, 0);
	size_t numTris = sg.triangles.size() / 3;
	for(size_t i = 1; i < numTris; ++i)
		q = min(q, GeometricApproximationDegree(sg, i));

	return q;
}

////////////////////////////////////////////////////////////////////////
number ShapeQualityDegree(SimpleGrid& sg, int triIndex)
{
	int i = triIndex * 3;
	return TriangleQuality_Area(sg.vertices[sg.triangles[i]],
								sg.vertices[sg.triangles[i+1]],
								sg.vertices[sg.triangles[i+2]]);
}

////////////////////////////////////////////////////////////////////////
number ShapeQualityDegree(SimpleGrid& sg)
{
	if(sg.triangles.size() < 3)
		return 0;

	number q = ShapeQualityDegree(sg, 0);
	size_t numTris = sg.triangles.size() / 3;
	for(size_t i = 1; i < numTris; ++i)
		q = min(q, ShapeQualityDegree(sg, i));

	return q;
}

////////////////////////////////////////////////////////////////////////
bool SwapEdge(SimpleGrid& sg)
{
	if(sg.triangles.size() < 6)
		return false;

//	get the two connected indices
	int ci[2];
	int counter = 0;
	for(size_t i = 0; i < 6; ++i){
		if(sg.triangles[i] > 1){
			if(counter < 2)
				ci[counter] = i;
			++counter;
		}
	}

//	both triangle 0 and 1 have to contain indices 0 and 1.
	if(counter != 2)
		return false;

//	the other indices of the first triangle
	int ind0 = sg.triangles[(ci[0] + 1) % 3];
	int ind1 = sg.triangles[(ci[0] + 2) % 3];

//	assign the vertex indices
	ci[0] = sg.triangles[ci[0]];
	ci[1] = sg.triangles[ci[1]];

//	create the new triangles
	sg.triangles[0] = ci[0];
	sg.triangles[1] = ind0;
	sg.triangles[2] = ci[1];

	sg.triangles[3] = ci[0];
	sg.triangles[4] = ci[1];
	sg.triangles[5] = ind1;

//	calculate the normals of the new triangles
	CalculateTriangleNormal(sg, 0);
	CalculateTriangleNormal(sg, 1);

//	done
	return true;
}

////////////////////////////////////////////////////////////////////////
bool CollapseEdge(SimpleGrid& sg)
{
	if(sg.triangles.size() < 6)
		return false;

//	the index of the new vertex
	int newInd = sg.vertices.size();
	
//	add the new vertex
	vector3 vNew;
	VecAdd(vNew, sg.vertices[0], sg.vertices[1]);
	VecScale(vNew, vNew, 0.5);
	sg.vertices.push_back(vNew);
	
//	and the new normal
	vector3 nNew;
	VecAdd(nNew, sg.vertexNormals[0], sg.vertexNormals[1]);
	VecNormalize(nNew, nNew);
	sg.vertexNormals.push_back(nNew);
	
//	create the new triangle list
//	the new list contains 2 triangles less than the old list
	vector<int> tris(sg.triangles.size() - 6);
	
	for(size_t i = 6; i < sg.triangles.size(); ++i){
		if(sg.triangles[i] < 2)
			tris[i-6] = newInd;
		else
			tris[i-6] = sg.triangles[i];
	}
	
//	swap the lists
	sg.triangles.swap(tris);
	
//	recalculate the normals
	CalculateTriangleNormals(sg);
	
//	done
	return true;
}

}//	end of namespace
