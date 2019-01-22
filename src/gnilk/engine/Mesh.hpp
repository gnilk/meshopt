#pragma once

#include <stdlib.h>
#include <math.h>

#include "vec.h"
#include "gnilk/engine/FaceHelper.hpp"
#include "gnilk/engine/core/Core.hpp"
#include "gnilk/engine/core/Edge.hpp"

#include <boost/unordered/unordered_map.hpp>

#include <vector>

namespace gnilk
{
	namespace engine
	{
		using namespace gnilk::engine::core;

		static float quadCoords[] = {
								-1.0f, 1.0f, 0.0f,
								1.0f, 1.0f, 0.0f,
								1.0f, -1.0f, 0.0f,
								-1.0f, -1.0f, 0.0f,

		};
		static int quadFaces[] = {
				3,2,1,0,
				-1,-1,-1,-1
		};



		static float cubeCoords[] = {
								-1.0f, 1.0f, 1.0f,
								1.0f, 1.0f, 1.0f,
								1.0f, -1.0f, 1.0f,
								-1.0f, -1.0f, 1.0f,

								-1.0f, 1.0f, -1.0f,
								1.0f, 1.0f, -1.0f,
								1.0f, -1.0f, -1.0f,
								-1.0f, -1.0f, -1.0f
		};

		/*     4---5
		 *    /   /|
		 *   0---1 6
		 *   |   |/
		 *   3---2
		 */
		static int cubeFaces[] = {
				3,2,1,0,
	//			-1,-1,-1,-1,
	//			-1,-1,-1,-1,
				4,5,6,7,
				1,5,4,0,
				2,3,7,6,
				5,1,2,6,
				0,4,7,3,
				-1,-1,-1,-1
		};

		static int extrudeCoordFaces[] = {
				 3,4,5,
				 0,3,5, 0,5,2,
				 0,4,3, 0,1,4,
				 5,4,1, 5,1,2,
				 -1, -1, -1,
				 -1, -1, -1, 1,3,0,
				 -1, -1, -1,
				 2,5,4, 2,4,1,
				 -1, -1, -1,

				// *  1,4,6,3 -> 1,4,6 ; 1,6,3
				// *  2,5,4,1 -> 2,5,4 ; 2,4,1
				// *  3,6,5,2 -> 3,6,5 ; 3,5,2

		};

		class Mesh
		{
		private:
			std::vector<Triangle *> triangles;
			std::vector<Edge *> edges;
			std::vector<Vertex3D> vertices;
			std::vector<Normal> faceNormals;
			std::vector<Normal> vertexNormals;
			std::vector<FaceGroup *> groups;


		public:
			float radius;
			Vertex3D mid;

			Mesh(std::vector<Triangle *> triangles, std::vector<Vertex3D> vertices, bool autoInit=true)
			{
				this->triangles = triangles;
				this->vertices = vertices;
				if (autoInit) {
					Initialize();
				}
			}
			virtual ~Mesh()
			{

			}

			void ClearData() {
				edges.clear();
			}

			// static gnilk::engine::Mesh *pMesh;
			static bool edgeLengthComparer(Edge *a, Edge *b) {
				
				return a->Length() < b->Length();
			}

			Mesh *Optimize() {

				std::vector<Triangle *> newTriangles;
				std::vector<Vertex3D> newVertices;
				// Triangles to be eliminated in new mesh
				std::vector<int> triDelete;
				std::vector<int> vtxDelete;

				int edgesBefore = edges.size();
				// remove 10%
				int numToDelete = 2000; // (int)(edges.size() * 0.10);
				int numFacesDeleted = 0;

				// TODO: need to recalc edges after this has been done
				std::sort (edges.begin(), edges.end(), edgeLengthComparer); 

				int f1 = edges[0]->f1;
				int f2 = edges[0]->f2;
				int v1 = edges[0]->v1;
				int v2 = edges[0]->v2;

				vertices.push_back(edges[0]->Mid());
				int idxMid = vertices.size()-1;


				triDelete.push_back(f1);
				triDelete.push_back(f2);		
				//
				// vIni(vertices[v1].data,0,0,0);
				// vIni(vertices[v2].data,0,0,0);

				for(auto tri : triangles) {
					tri->ReplaceVertex(v1, idxMid);
					tri->ReplaceVertex(v2, idxMid);
				}
				// Rebuild new triangle tables
				for(int i=0;i<triangles.size();i++) {
					bool bSkip = false;
					for(int j=0;j<triDelete.size();j++) {
						if (triDelete[j]==i) {
							bSkip = true;
							break;
						}
					}
					if (!bSkip) {
						newTriangles.push_back(triangles[i]);			
					}
				}

				// Create the new mesh
				Mesh *newMesh = new Mesh(newTriangles, vertices);
				return newMesh;
			}
			void Initialize()
			{
				ClearData();
				CalculateEdges();
				CalculateFaceNormals();
				CalculateBoundingSphere();
				// AutoGroupFaces();
				// CalculateVertexNormals();

			}

			void CalculateBoundingSphere()
			{
				radius = FaceHelper::CalculateBoundingSphere(mid, vertices);
			}
			void AutoGroupFaces()
			{
				// Clear groups if already present
				if (groups.size() > 0) {
					groups.clear();
				}

				FaceHelper::CalculateFaceGroups(groups, edges, triangles, faceNormals);
			}
			void CalculateEdges()
			{
				FaceHelper::CalculateEdges(edges, triangles, vertices);
			}

			void CalculateVertexNormals()
			{
				FaceHelper::CalculateVertexNormals(vertices.size(), groups, triangles);
			}

			void CalculateFaceNormals()
			{
				//
				// calculate a vertex -> edge list
				// vertex normal is simply an iteration of all faces normals in the edge list
				//
				//FaceHelper::CalculateFaceNormals(faceNormals, triangles, vertices);
				FaceHelper::CalculateFaceNormals(triangles, vertices);
			}


			// Getters
			std::vector<Vertex3D> &GetVertices()
			{
				return this->vertices;
			}

			std::vector<Triangle *> &GetTriangles()
			{
				return this->triangles;
			}

			std::vector<Edge *> &GetEdges()
			{
				return this->edges;
			}

			// std::vector<Normal> &GetFaceNormals()
			// {
			// 	return this->faceNormals;
			// }

			std::vector<FaceGroup *> &GetGroups()
			{
				return this->groups;
			}

			// Stupid function that allows making of "spheres"
			// 1) Create a cube
			// 2) Subdivide X number of times
			// 3) Normalize the last mesh
			void NormalizeCoords() {
				for(int i=0;i<vertices.size();i++) {
					vNorm(vertices[i].data, vertices[i].data);
				}
				CalculateFaceNormals();
				CalculateBoundingSphere();
				AutoGroupFaces();
				CalculateVertexNormals();
			}

			Mesh *SplitFaces() {
				std::vector<Triangle *> newTris;
				std::vector<Vertex3D> newCoords (this->vertices);
				boost::unordered::unordered_map<int, int> edgemidmap;

				for(int i=0;i<edges.size();i++) {
					int idx = FaceHelper::SplitEdge(newCoords, vertices, edges[i]->v1, edges[i]->v2);
					edgemidmap.insert(std::pair<int, int>(edges[i]->HashCode(), idx));
				}

				for(int i=0;i<triangles.size();i++) {
					Triangle *pTri = triangles[i];

					int m1 = edgemidmap[Edge::HashCode(pTri->v[0],pTri->v[1])];
					int m2 = edgemidmap[Edge::HashCode(pTri->v[1],pTri->v[2])];
					int m3 = edgemidmap[Edge::HashCode(pTri->v[2],pTri->v[0])];

					newTris.push_back(new Triangle(pTri->v[0], m1, m3));
					newTris.push_back(new Triangle(pTri->v[1], m2, m1));
					newTris.push_back(new Triangle(pTri->v[2], m3, m2));
					newTris.push_back(new Triangle(m1, m2, m3));
				}
				return new Mesh(newTris, newCoords);
			}

			void Extrude(std::vector<Triangle *> triangles) {

			}

			void Extrude(int triIdx, float factor) {
				// Note:
				//  - Need to base extrusion on edges rather than triangles
				//  - Must avoid edges sharing two faces (no extrusion of such edges)

				Triangle *triangle = this->triangles[triIdx];
				int v1 = triangle->v[0];
				int v2 = triangle->v[1];
				int v3 = triangle->v[2];

				printf("Extruding\n");
				printf("Tri: %d,%d,%d\n",v1,v2,v3);
				printf("Triangles start: %d\n", this->triangles.size());
				printf("Vertices start: %d\n", this->vertices.size());
				printf("Facenormal: %f, %f, %f\n",triangle->normal.data[0],triangle->normal.data[1],triangle->normal.data[2]);
				printf("Factor: %f\n",factor);

				int v1New = this->vertices.size();
				// Duplicate the three vertices
				this->vertices.push_back(this->vertices[v1]);
				this->vertices.push_back(this->vertices[v2]);
				this->vertices.push_back(this->vertices[v3]);

				// Move vertices
				float tmp[3];
				vAdd(this->vertices[v1New+0].data, this->vertices[v1].data, vMul(tmp, triangle->normal.data, factor));
				vAdd(this->vertices[v1New+1].data, this->vertices[v2].data, vMul(tmp, triangle->normal.data, factor));
				vAdd(this->vertices[v1New+2].data, this->vertices[v3].data, vMul(tmp, triangle->normal.data, factor));

				printf("v1 (%d): %f, %f, %f\n",v1, this->vertices[v1].data[0],this->vertices[v1].data[1],this->vertices[v1].data[2]);
				printf("v2 (%d): %f, %f, %f\n",v2, this->vertices[v2].data[0],this->vertices[v2].data[1],this->vertices[v2].data[2]);
				printf("v3 (%d): %f, %f, %f\n",v3, this->vertices[v3].data[0],this->vertices[v3].data[1],this->vertices[v3].data[2]);

				printf("n1 (%d): %f, %f, %f\n",v1New+0, this->vertices[v1New+0].data[0],this->vertices[v1New+0].data[1],this->vertices[v1New+0].data[2]);
				printf("n2 (%d): %f, %f, %f\n",v1New+1, this->vertices[v1New+1].data[0],this->vertices[v1New+1].data[1],this->vertices[v1New+1].data[2]);
				printf("n3 (%d): %f, %f, %f\n",v1New+2, this->vertices[v1New+2].data[0],this->vertices[v1New+2].data[1],this->vertices[v1New+2].data[2]);

				// Construct new faces
				// 1) replace old face with new vertex information
//				triangle->v[0] = v1New+0;
//				triangle->v[1] = v1New+1;
//				triangle->v[2] = v1New+2;

				printf("Added vertices, new size: %d\n",this->vertices.size());
				printf("Ext.CreateFaces: %d, %d\n",v1, v1New);
				int idx = 0;
				int tbl[]={v1,v2,v3};
				while(extrudeCoordFaces[idx]!=-1) {

					Triangle *pNew = new Triangle();
					for(int i=0;i<3;i++) {
						pNew->v[i] = extrudeCoordFaces[idx+i];
						// Add old vertex index on top of relative vertex index
						if (pNew->v[i] >= 3) pNew->v[i]+=v1New-3;
						else {
							pNew->v[i] = tbl[extrudeCoordFaces[idx+i]];
						}
					}
					printf("t: %d,%d,%d\n",pNew->v[0], pNew->v[1], pNew->v[2]);
					this->triangles.push_back(pNew);
					idx+=3;
				}
				/*
				 *  1,4,6,3 -> 1,4,6 ; 1,6,3
				 *  2,4,4,1 -> 2,5,4 ; 2,4,1
				 *  3,6,5,2 -> 3,6,5 ; 3,5,2
				 *
				 *
				 *              1 --------- 4
				 *            / |          / |
				 *         3    |         6  |
				 *           \  |          \ |
				 *              2 ---------  5
				 *
				 */
			}

			static Mesh *CreateQuad() {
				std::vector<Vertex3D> coords;
				std::vector<Triangle *> triangles;

				for (int i=0;i<4;i++) {
					Vertex3D vertex;
					vDup(vertex.data, &quadCoords[i*3]);
					coords.push_back(vertex);
					
				}
				int i = 0;
				while(quadFaces[i] != -1) {

					Triangle *tri = new Triangle();
					tri->faceData = kFaceData_Vertex;
					tri->v[0] = quadFaces[i+0];
					tri->v[1] = quadFaces[i+1];
					tri->v[2] = quadFaces[i+2];
					triangles.push_back(tri);

					tri = new Triangle();
					tri->faceData = kFaceData_Vertex;
					tri->v[0] = quadFaces[i+0];
					tri->v[1] = quadFaces[i+2];
					tri->v[2] = quadFaces[i+3];

					triangles.push_back(tri);
					i+=4;
				}

				Mesh *pMesh = new Mesh(triangles, coords);
				return pMesh;
			}

			static Mesh *CreateCube() {
				std::vector<Vertex3D> coords;
				std::vector<Triangle *> triangles;
				for(int i=0;i<8;i++) {
					Vertex3D vertex;
					vDup(vertex.data, &cubeCoords[i*3]);
					coords.push_back(vertex);
				}

				int i = 0;
				while(cubeFaces[i] != -1) {

					Triangle *tri = new Triangle();
					tri->faceData = kFaceData_Vertex;
					tri->v[0] = cubeFaces[i+0];
					tri->v[1] = cubeFaces[i+1];
					tri->v[2] = cubeFaces[i+2];
					triangles.push_back(tri);

					tri = new Triangle();
					tri->faceData = kFaceData_Vertex;
					tri->v[0] = cubeFaces[i+0];
					tri->v[1] = cubeFaces[i+2];
					tri->v[2] = cubeFaces[i+3];

					triangles.push_back(tri);
					i+=4;
				}

				Mesh *pMesh = new Mesh(triangles, coords);
				return pMesh;
			}
			static Mesh *CreateSphere(float radius, int segments, int slices) {
				std::vector<Vertex3D> coords;
				std::vector<Triangle *> triangles;

				float mat[4*4];
				int nvps = segments - 1; // vertices per slices
				int nlps = segments * 2 - 1; // lines per slice

				Vertex3D top,bottom;
				vIni(top.data,0,radius,0);
				vIni(bottom.data,0,-radius,0);
				coords.push_back(top);
				coords.push_back(bottom);

				for(int i=0 ; i<slices ; i++)
				{
					vMatIdentity(mat);
					vMatFromAngles(mat, 0,2*M_PI*i/(float)slices,0);
					for(int j=0 ; j<nvps ; j++)
					{
						Vertex3D vertex;
						float v[3];
						v[0] = radius*sin(M_PI*(j+1)/(float)segments);
						v[1] = radius*cos(M_PI*(j+1)/(float)segments);
						v[2] = 0;
						vApplyMat(vertex.data, v, mat);
						coords.push_back(vertex);
					}
				}

				for(int i=0 ; i<slices ; i++)
				{
					Triangle *tri;
					// "vertikala" linjer och polygoner
					for(int j=0 ; j<nvps-1 ; j++)
					{
						tri = new Triangle();
						tri->faceData = kFaceData_Vertex;
						tri->v[0] = i*nvps+j+2+1;
						tri->v[1] = i*nvps+j+2;
						tri->v[2] = ((i+1)%slices)*nvps+j+2;
						triangles.push_back(tri);

						tri = new Triangle();
						tri->faceData = kFaceData_Vertex;
						tri->v[0] = ((i+1)%slices)*nvps+j+2;
						tri->v[1] = ((i+1)%slices)*nvps+j+2+1;
						tri->v[2] = i*nvps+j+2+1;
						triangles.push_back(tri);
					}

					tri = new Triangle();
					tri->faceData = kFaceData_Vertex;
					tri->v[0] = i*nvps+2;
					tri->v[1] = 0;
					tri->v[2] = ((i+1)%slices)*nvps+2;
					triangles.push_back(tri);

					tri = new Triangle();
					tri->faceData = kFaceData_Vertex;
					tri->v[0] = i*nvps+2+nvps-1;
					tri->v[1] = ((i+1)%slices)*nvps+2+nvps-1;
					tri->v[2] = 1;
					triangles.push_back(tri);
				}
				Mesh *pMesh = new Mesh(triangles, coords);
				return pMesh;
			}

		};
	}
}
