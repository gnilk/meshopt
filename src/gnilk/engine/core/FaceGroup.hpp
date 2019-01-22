#pragma once

#include <assert.h>
#include <math.h>

#include "gnilk/engine/core/Core.hpp"

namespace gnilk
{
	namespace engine
	{
		namespace core
		{
			class FaceGroup
			{
			public:
				std::vector<int> triangles;
				Normal normal;
				Normal baseNormal;
			public:
				FaceGroup()
				{

				}
				virtual ~FaceGroup()
				{

				}
			};

		}
	}
}
