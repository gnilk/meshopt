#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "vec.h"
#include "gnilk/engine/core/Core.hpp"

namespace gnilk
{
	namespace engine
	{
		namespace core
		{

		class Edge
		{
		public:
			int v1, v2;
			int f1, f2;
			float n1[3], n2[3];
			int sindex;
			const std::vector<Vertex3D> &vertices;			
		public:
			Edge(const std::vector<Vertex3D> &_vertices) : vertices(_vertices)
			{
				v1 = v2 = 0;
				f1 = f2 = -1;	// negative face number
				vIni(n1,0,0,0);
				vIni(n2,0,0,0);
				sindex = 0;
			}

			Edge(int v1, int v2, const std::vector<Vertex3D> &_vertices) : vertices(_vertices)
			{
				this->v1 = v1;
				this->v2 = v2;
				f1 = f2 = -1;	// negative face number
				vIni(n1,0,0,0);
				vIni(n2,0,0,0);
				sindex = 0;
			}
			virtual ~Edge()
			{

			}

			int HashCode() {
				return Edge::HashCode(v1,v2);
			}

			static int HashCode(int v1, int v2)
			{
				return (v1>v2) ? ((v1<<16) | v2) : ((v2<<16) | v1);
			}

			bool Equals(Edge &pOther)
			{
				if (((pOther.GetV1() == v1) && (pOther.GetV2() == v2)) ||
					((pOther.GetV1() == v2) && (pOther.GetV2() == v1))) return true;
				return false;
			}

			std::string ToString()
			{
				std::string foo;
				foo.resize(64);
				snprintf(&foo[0],64,"%d,%d",v1,v2);
				return foo;
			}
			bool IsClosed()
			{
				if ((f1!=-1) && (f2!=-1)) return true;
				return false;
			}
			bool AddFace(int faceIndex)
			{
				if ((f1!=-1) && (f2!=-1)) return false;
				if (f1 == -1) f1 = faceIndex;
				else f2 = faceIndex;
				return true;
			}

			bool AddNormal(float *pNormal)
			{
				if (sindex == 2) return false;
				if (sindex == 0) vDup(n1, pNormal);
				if (sindex == 1) vDup(n2, pNormal);
				sindex++;
				return true;
			}

			Vertex3D Mid() {
				float tmp[3];
				Vertex3D mid;
				vAdd(tmp, (float *)&vertices[v2].data[0], (float *)&vertices[v1].data[0]);
				vMul(mid.data, tmp, 0.5);
				return mid;
				
			}

			float Length() {
				float tmp[3];
				vSub(tmp, (float *)&vertices[v2].data[0], (float *)&vertices[v1].data[0]);
				return vAbs(tmp);
			}

			// properties
			int GetV1() { return v1; }
			int GetV2() { return v2; }
			float *GetN1() {return n1; }
			float *GetN2() { return n2; }
			void SetN1(float *nv) { vDup(n1,nv); }
			void SetN2(float *nv) { vDup(n2,nv); }
		}; // Edge class

		}
	}
}
