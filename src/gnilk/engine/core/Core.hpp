#pragma once

#include <assert.h>
#include <math.h>

namespace gnilk
{
namespace engine
{
	namespace core
	{
		typedef enum {
			kFaceData_None = 0,
			kFaceData_Vertex = 1,
			kFaceData_VertexTexture = 2,
			kFaceData_VertexNormal = 4,
		} kFaceData;

		typedef enum {
			kFaceType_Unknown = 0,
			kFaceType_Tri = 1,
			kFaceType_Quad = 2,
		} kFaceType;

		struct Face {
			kFaceData faceData;
			kFaceType faceType;
			int v[4];
			int vn[4];
			int uv[4];
		};

		struct Normal {
			float data[3];
		};


		struct Vertex3D {	// Vertex coords
			float data[3];
		};
		struct Vertex2D {	// U,V coords
			float data[2];
		};


	}
}

}
