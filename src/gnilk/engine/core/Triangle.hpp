#pragma once

#include <assert.h>
#include <math.h>

#include "gnilk/engine/core/Core.hpp"
#include "gnilk/engine/core/Edge.hpp"

namespace gnilk
{
	namespace engine
	{
		namespace core
		{
			class Triangle
			{
			public:
				kFaceData faceData;
				int v[3];
				//int vn[3];
				Normal vn[3];
				int uv[3];
				Normal normal;
				// edge list is unsorted, use v-table + edge-hash to map correctly
				// Note: Sorting edge list wont give you correct winding since edges are shared
				// between faces, thus face1 may have edge A->B while face2 have same edge B->A
				Edge *edges[3];
			public:
				Triangle(int v1, int v2, int v3) {
					SetIndex(v1,v2,v3);
					Initialize();
				}

				Triangle(int v1, int v2, int v3, int uv1, int uv2, int uv3) {

					SetIndex(v1,v2,v3);
					SetUVIndex(uv1, uv2, uv3);
				}

				Triangle() {
					Initialize();
				}

				bool ShareVertex(int vtx) {
					if (v[0] == vtx) return true;
					if (v[1] == vtx) return true;
					if (v[2] == vtx) return true;
					return false;
				}

				bool ReplaceVertex(int oldVtx, int newVtx) {
					bool bReplaced = false;
					if (v[0] == oldVtx) {
						v[0] = newVtx;
						bReplaced = true;
					} 

					if (v[1] == oldVtx) {
						v[1] = newVtx;
						bReplaced = true;
					} 

					if (v[2] == oldVtx) {
						v[2] = newVtx;
						bReplaced = true;
					}
					return bReplaced;
				}

				void SetIndex(int v1, int v2, int v3) {
					v[0] = v1;
					v[1] = v2;
					v[2] = v3;
					Reset();
				}

				void SetUVIndex(int uv1, int uv2, int uv3) {
					uv[0] = uv1;
					uv[1] = uv2;
					uv[2] = uv3;
				}

				void Reset() {
					Initialize();
				}

				bool HasNormal() {
					return normalAssigned;
				}

				void SetNormal(Normal normal) {
					this->normal = normal;
				}

				void AddEdge(Edge *edge) {

					if (edgeCount > 2)
					{
						printf("triangle is fucked up, too many edges: %d\n",edgeCount);
					} else
					{
						edges[edgeCount++] = edge;
//						if (edgeCount == 3) {
//							SortEdges();
//						}
					}
				}

			private:
				int edgeCount;
				bool normalAssigned;

				void Initialize() {
					edgeCount = 0;
					normalAssigned = false;
				}

//				// Make sure edges are sorted in according to the winding
//				void SortEdges() {
//					Edge *e1 = GetEdge(v[0],v[1]);
//					Edge *e2 = GetEdge(v[1],v[2]);
//					Edge *e3 = GetEdge(v[2],v[0]);
//					edges[0] = e1;
//					edges[1] = e2;
//					edges[2] = e3;
//				}
//				Edge *GetEdge(int v1, int v2) {
//					int hc = Edge::HashCode(v1,v2);
//					for(int i=0;i<edgeCount;i++) {
//						if (edges[i]->HashCode() == hc) return edges[i];
//					}
//					return NULL;
//				}

			};

		}
	}
}
