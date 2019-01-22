#pragma once

#include <stdlib.h>
#include <math.h>

#include "vec.h"
#include "gnilk/engine/core/Core.hpp"
#include "gnilk/engine/core/Edge.hpp"
#include "gnilk/engine/core/Triangle.hpp"
#include "gnilk/engine/core/FaceGroup.hpp"

#include <boost/unordered/unordered_map.hpp>

#include <vector>

namespace gnilk
{
	namespace engine
	{
		using namespace gnilk::engine::core;

		class FaceHelper
		{
		public:
			//
			// calculates vertex normals per group basis
			// vertex normals are stored directly in to the triangles
			// I did not want unsynchronized buffer sizes nor buffer dependencies on group
			// This is not optimal for rendering - that has to be solved when computing render buffers
			//
			static void CalculateVertexNormals(int vertices, std::vector<int> &faces, std::vector<Triangle *> &triangles) {
				std::vector<Normal> normals;
				normals.reserve(vertices);	// need intermediate storage
				// need to initialize array
				for(int i=0;i<vertices;i++) {
					vIni(normals[i].data,0,0,0);
				}
				// add up normals per vertex - note: normals should already be normalized!
				for(int i=0;i<faces.size();i++) {
					int idxTri = faces[i];
					Triangle *tri = triangles[idxTri];
					for(int j=0;j<3;j++) {
						vAdd(normals[tri->v[j]].data, normals[tri->v[j]].data, tri->normal.data);
					}
				}
				// normalize and assign
				for(int i=0;i<faces.size();i++) {
					Triangle *tri = triangles[faces[i]];
					for(int j=0;j<3;j++) {
						vNorm(tri->vn[j].data,normals[tri->v[j]].data);
					}
				}
			}

			static void CalculateVertexNormals(int vertices, std::vector<FaceGroup*> &groups, std::vector<Triangle *> &triangles)
			{
				for(int i=0;i<groups.size();i++)
				{
					printf("CVN, G: %d, T: %d\n",i, groups[i]->triangles.size());
					CalculateVertexNormals(vertices, groups[i]->triangles, triangles);
				}
			}

			static void CalculateFaceNormal(Normal *faceNormal, std::vector<Vertex3D> &coords, int v1, int v2, int v3) {
				float *pV1 = coords[v1].data;
				float *pV2 = coords[v2].data;
				float *pV3 = coords[v3].data;

				float vE1[3];
				float vE2[3];
				vSub(vE1, pV3, pV1);
				vSub(vE2, pV2, pV1);
				vCross(faceNormal->data, vE1, vE2);
				vNorm(faceNormal->data,faceNormal->data);
			}

			// Calculate normals for a list of faces and assign to face
			static void CalculateFaceNormals(std::vector<Triangle *> &triangles, std::vector<Vertex3D> &coords) {
				for(int i=0;i<triangles.size();i++) {
					Triangle *tri = triangles[i];
					if (!tri->HasNormal()) {
						CalculateFaceNormal(tri->normal.data, coords, tri->v[0], tri->v[1], tri->v[2]);
					}
				}
			}

			// Calculates a single face normal
			static void CalculateFaceNormal(float *normal, std::vector<Vertex3D> &coords, int v1, int v2, int v3) {
				float *pV1 = coords[v1].data;
				float *pV2 = coords[v2].data;
				float *pV3 = coords[v3].data;

				float vE1[3];
				float vE2[3];
				vSub(vE1, pV3, pV1);
				vSub(vE2, pV2, pV1);
				vCross(normal, vE2, vE1);
				vNorm(normal, normal);
			}

			// calculate list of normals for each face
			static void CalculateFaceNormals(std::vector<Normal> &normals, std::vector<Triangle *> &triangles, std::vector<Vertex3D> &coords) {
				normals.reserve(triangles.size());
				for(int i=0;i<triangles.size();i++) {
					CalculateFaceNormal(normals[i].data, coords, triangles[i]->v[0], triangles[i]->v[1], triangles[i]->v[2]);
				}
			}

			// Gets or adds (a new) group for face depending on the deviation between face normal and group normal
			// Not used..
			static FaceGroup *GetOrAddGroup(std::vector<FaceGroup *> &groups, Normal &normal, float threshold) {
				FaceGroup *group = NULL;
				float lastFac = threshold;

				// Find group (if any) which fulfills the the normal threshold
				for(size_t i=0;i<groups.size();i++) {
					float fac = 1.0 - vDot(normal.data, groups[i]->baseNormal.data);
					if ((fac > 0) && (fac < threshold)) {
						group = groups[i];
						lastFac = fac;
					}
				}
				if (group == NULL)
				{
					group = new FaceGroup();
					vDup(group->baseNormal.data, normal.data);
					groups.push_back(group);
				}
				return group;
			}

			// Recursive, assigns faces to groups - based on edges and if two faces have a similar normal they are allowed in the same group
			static void DoCalculateFaceGroups(FaceGroup *pGroup,
					int triangle,
					float threshold,
					std::vector<bool> &visited,
					std::vector<Edge *> &edges,
					std::vector<Triangle *> &triangles,
					std::vector<Normal> &faceNormals,
					int depth)
			{
				if (visited[triangle])
				{
					//printf("%d:%d -> .\n",depth,triangle);
					return;
				}

				pGroup->triangles.push_back(triangle);
				visited[triangle] = true;

				Triangle *pTri = triangles[triangle];
				for(int i=0;i<3;i++)
				{
					Edge *edge = triangles[triangle]->edges[i];
					if (edge == NULL)
					{
						printf("NULL edge %d for triangle %d\n",i, triangle);
						continue;
					}
					int tNext = -1;
					if ((edge->f1 != triangle) && (edge->f1 != -1)) tNext = edge->f1;
					else if ((edge->f2 != triangle) && (edge->f2 != -1)) tNext = edge->f2;

					if (tNext >= 0)
					{
						Triangle *pTriNext = triangles[tNext];
						float faceDiff = 1.0 - vDot(pTri->normal.data, pTriNext->normal.data);
						//float faceDiff = 1.0 - vDot(faceNormals[triangle].data, faceNormals[tNext].data);
						faceDiff = fabs(faceDiff);
						if ((faceDiff >= 0) && (faceDiff <= threshold))
						{
							//printf("%d:%d -> %d - %f\n",depth,triangle, tNext, faceDiff);
							DoCalculateFaceGroups(pGroup, tNext, threshold, visited, edges, triangles, faceNormals,depth+1);
						} else
						{
							// rejected because of threshold
							//printf("%d, r: %f\n",faceDiff);
						}
					}
				} // Foreach edge
				//printf("%d:%d -> x\n",depth,triangle);
			}

			// Automatic group of faces depending on
			static void CalculateFaceGroups(std::vector<FaceGroup *> &groups, std::vector<Edge *> &edges, std::vector<Triangle *> &triangles, std::vector<Normal> &faceNormals)
			{
				std::vector<bool> visited;
				visited.reserve(triangles.size());
				for(size_t i=0;i<triangles.size();i++)
				{
					visited.push_back(false);
				}
				for(size_t i=0;i<triangles.size();i++)
				{
					if (!visited[i])
					{
						FaceGroup *pGroup = new FaceGroup();
						DoCalculateFaceGroups(pGroup, i, 0.5, visited, edges, triangles, faceNormals,0);
						printf("Group %d, triangles: %d\n", groups.size(), pGroup->triangles.size());
						groups.push_back(pGroup);
					}
				}

				for(size_t i=0;i<groups.size();i++)
				{
					FaceGroup *group = groups[i];
					vMul(group->normal.data, group->normal.data, 1.0 / (double)group->triangles.size());
				}
				printf("Groups: %d\n",groups.size());
			}

			// Calculate edges
			static void CalculateEdges(std::vector<Edge *> &edges, std::vector<Triangle *> &triangles, const std::vector<Vertex3D> &vertices)
			{
				boost::unordered::unordered_map<int, Edge *> edgemap;
				for(size_t i=0;i<triangles.size();i++)
				{
					triangles[i]->Reset();
					GetOrAddEdgeMap(edgemap, triangles, vertices, i, 0, 1);
					GetOrAddEdgeMap(edgemap, triangles, vertices, i, 1, 2);
					GetOrAddEdgeMap(edgemap, triangles, vertices, i, 2, 0);
				}
				// I guess there is a better way to generate a vector instead of this...
				int count = 0;
				for (boost::unordered::unordered_map<int, Edge *>::iterator it = edgemap.begin(); it != edgemap.end(); ++it )
				{
					Edge *edge = it->second;
					if (edge->f1 >= 0) triangles[edge->f1]->AddEdge(edge);
					if (edge->f2 >= 0) triangles[edge->f2]->AddEdge(edge);
					edges.push_back(it->second);
					count++;
				}
			} // CalculateEdges

			// Convert faces to triangles - quads are split
			static void FacesToTriangles(
					std::vector<Triangle *> &triangles,
					std::vector<Face *> &faces) {

				for(size_t i=0;i<faces.size();i++) {
					Triangle *tri = new Triangle();
					tri->faceData = faces[i]->faceData;
					for(int j=0;j<3;j++) {
						tri->v[j] = faces[i]->v[j];
						//tri->vn[j] = faces[i]->vn[j];
						tri->uv[j] = faces[i]->uv[j];
					}
					triangles.push_back(tri);
					if (faces[i]->faceType == kFaceType_Quad) {
						tri = new Triangle();
						tri->faceData = faces[i]->faceData;
						for(int j=0;j<3;j++) {
							int idx = ((j+2)>3) ? 0 : j+2;
							tri->v[j] = faces[i]->v[idx];
							//tri->vn[j] = faces[i]->vn[idx];
							tri->uv[j] = faces[i]->uv[idx];
						}
						triangles.push_back(tri);
					}
				}
			} // FacesToTriangles

			static void CollapseEdge(int edge, std::vector<Triangle *> &triangles, std::vector<Edge *> &edges, std::vector<Vertex3D> &coords)
			{
				Edge *pEdge = edges[edge];
				int f1 = pEdge->f1;
				int f2 = pEdge->f2;
				int v1 = pEdge->v1;
				int v2 = pEdge->v2;
				// will this work?
				triangles.erase(triangles.begin()+f1);
				triangles.erase(triangles.begin()+f2);

				edges.erase((edges.begin()+edge));
			}

			static int CountDuplicateVertices(std::vector<Vertex3D> &coords) {
				int duplicates = 0;
				for(int i=0;i<coords.size();i++) {
					float temp[3];
					vNorm(temp,coords[i].data);
					for(int j=0;j<coords.size();j++) {
						if (i==j) continue;
						float temp2[3];
						vNorm(temp2, coords[j].data);
						float fac = vDot(temp,temp2);
						if (fac>0 && fac<0.01) {
							duplicates++;
						}
					}
				}
				return duplicates;
			}

			// Calculates bound sphere radius and mid point
			static double CalculateBoundingSphere(Vertex3D &mid, std::vector<Vertex3D> &coords) {
				Vertex3D cmax,cmin;
				float tmp[3];
				cmax = coords[0];
				cmin = coords[1];
				// Find max/min vertex
				for(int i=0;i<coords.size();i++) {
					cmax.data[0] = fmax(coords[i].data[0], cmax.data[0]);
					cmax.data[1] = fmax(coords[i].data[1], cmax.data[1]);
					cmax.data[2] = fmax(coords[i].data[2], cmax.data[2]);

					cmin.data[0] = fmin(coords[i].data[0], cmin.data[0]);
					cmin.data[1] = fmin(coords[i].data[1], cmin.data[1]);
					cmin.data[2] = fmin(coords[i].data[2], cmin.data[2]);
				}

				// Calculate radius and mid-point
				// Note: Why 2.5???
				double r = vAbs(vSub(tmp,cmax.data, cmin.data)) / 2.5;
				vAdd(tmp, cmax.data, cmin.data);
				vMul(mid.data, tmp, 0.5);
				return r;
			}

			static int SplitEdge(std::vector<Vertex3D> &coordOut, std::vector<Vertex3D> &coords, int v1, int v2) {
				Vertex3D mid;
				vMul(mid.data, vAdd(mid.data,coords[v1].data,coords[v2].data),0.5f);
				coordOut.push_back(mid);
				return (coordOut.size()-1);
			}

			/*      2
			 *     4 5
			 *    1 6 3
			 *
			 * Tris:
			 * 	1,4,6
			 * 	2,5,4
			 * 	3,6,5
			 * 	4,5,6
			 *
			*/
			static void SplitTriangle(std::vector<Triangle *> &triOut, std::vector<Vertex3D> &coordOut, Triangle *triangle, std::vector<Vertex3D> &coords) {
				int v1 = triangle->v[0];
				int v2 = triangle->v[1];
				int v3 = triangle->v[2];

				int idxStart = coordOut.size();
				coordOut.push_back(coords[v1]);
				coordOut.push_back(coords[v2]);
				coordOut.push_back(coords[v3]);
				SplitEdge(coordOut, coords, v1, v2);
				SplitEdge(coordOut, coords, v2, v3);
				SplitEdge(coordOut, coords, v3, v1);
				// Create 4 new triangles

				triOut.push_back(new Triangle(idxStart+0, idxStart+3, idxStart+5));
				triOut.push_back(new Triangle(idxStart+1, idxStart+4, idxStart+3));
				triOut.push_back(new Triangle(idxStart+2, idxStart+5, idxStart+4));
				triOut.push_back(new Triangle(idxStart+3, idxStart+4, idxStart+5));
			}

		private:
			static void GetOrAddEdgeMap(boost::unordered::unordered_map<int, Edge *> &edges, 
				std::vector<Triangle *> &triangles, 
				const std::vector<Vertex3D> &vertices, 
				int face, int v1, int v2) {
				Edge *pEdge = NULL;

				int hashCode = Edge::HashCode(triangles[face]->v[v1], triangles[face]->v[v2]);
				//printf("hc: %d (%d, %d)\n",hashCode,faces[face]->v[v1], faces[face]->v[v2]);
				if (edges.find(hashCode) == edges.end()) {
					pEdge = new Edge(triangles[face]->v[v1], triangles[face]->v[v2], vertices);
					pEdge->AddFace(face);
					edges.insert(std::pair<int, Edge *>(hashCode, pEdge));
				} else {
					if (edges[hashCode]->AddFace(face) != true) {
						//printf("!!!! MULTIPLE FACES SHARING SAME EDGE!!!\n");
						// TODO: Create new face!
					}
				}
			} // GetOrAddEdgeMap
		}; // class FaceHelper
	}
}
