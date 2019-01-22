#pragma once

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <math.h>

#include "vec.h"

#include "gnilk/engine/core/Edge.hpp"
#include "gnilk/engine/core/Core.hpp"
#include "gnilk/engine/FaceHelper.hpp"
#include "gnilk/engine/Mesh.hpp"
#include "gnilk/Tokenizer.hpp"

#include <boost/unordered/unordered_map.hpp>

using namespace gnilk::engine::core;
using namespace gnilk;

namespace gnilk
{
	class Surface
	{

	};

	namespace formats
	{
		namespace wavefront
		{
			typedef enum  {
				kParseResult_Unknown = 0,
				kParseResult_Comment,
				kParseResult_Vertex,
				kParseResult_VertexNormal,
				kParseResult_VertexTexture,
				kParseResult_Face,
			} kParseResult;

			class ObjReader
			{
			public:
				void DumpCoords() {
					int idx =0;
					std::vector<Vertex3D>::iterator it=coords.begin();
					while(it!=coords.end()) {
						Vertex3D vtx = *it;
						printf("%d,%f,%f,%f\n",idx,vtx.data[0],vtx.data[1],vtx.data[2]);
						++it;
						idx++;
					}
				}

				//
				// Parses vertex data, no support for optional and strict length handling
				//
				float *ParseVertex(Tokenizer &tok, float *dst, int maxItems) {
					int idx = 0;
					std::string tmp;

					while(!(tmp = tok.Next()).empty()) {
							dst[idx] = strtod(tmp.c_str(), NULL);
							++idx;
							if (idx >= maxItems) break;
					}

					return dst;
				}

				gnilk::engine::core::Face *ParseFace(Tokenizer &tok, gnilk::engine::core::Face *dstFace) {
					int idx=0;
					int idxComponent = 0;
					int dataFlags = 0;
					bool bExpectOp = false;
					std::string tmp;

					//
					// face format is (index based):
					//  v/vt/vn
					// where v - vertex, vt - vertex texture, vn - vertex normal
					//
					// also possible is: v//vn (i.e. vertex and normal but no texture)
					//
					while(!(tmp = tok.Next()).empty()) {
						if (tmp == "/") {
							bExpectOp = false;
							idxComponent++;	// next component for vertex (v, vt, vn)
						} else {
							if (bExpectOp == true) {
								// two numericals after each-other -> new face vertex!
								idxComponent=0;
								idx++;	// next index in face
							}
							int val = (int)strtol(tmp.c_str(), NULL,10);
							switch (1<<idxComponent) {
								case kFaceData_Vertex:
									dstFace->v[idx] = val - 1;
									break;
								case kFaceData_VertexTexture :
									dstFace->uv[idx] = val - 1;
									break;
								case kFaceData_VertexNormal :
									dstFace->vn[idx] = val - 1;
								default:
									break;
							}
							bExpectOp = true;
							dataFlags |= 1 << (idxComponent);
						}
						if (idx >= 3) break;
					}
					dstFace->faceData = (kFaceData)dataFlags;
					// Zero based.. (0..3=quad, 0..2=tri)
					switch(idx) {
						case 3 :
							dstFace->faceType = kFaceType_Quad;
							break;
						case 2 :
							dstFace->faceType = kFaceType_Tri;
							break;
					}
					return dstFace;
				}

				kParseResult ClassifyLineToken(std::string &tok) {
					static std::string tokens("# v f vn vt");
					kParseResult result = kParseResult_Unknown;
					switch(gnilk::Tokenizer::Case(tok, tokens)) {
						case 0 : // #
							result = kParseResult_Comment;
							break;
						case 1 : // v
							result = kParseResult_Vertex;
							break;
						case 2 : // f
							result = kParseResult_Face;
							break;
						case 3 : // vn
							result = kParseResult_VertexNormal;
							break;
						case 4 : // vt
							result = kParseResult_VertexTexture;
							break;
					}
					return result;
				}

				kParseResult ParseLine(std::string buffer) {
					kParseResult result = kParseResult_Unknown;
					//printf("ParseLine: %s\n",buffer);
					Tokenizer tokenizer(buffer, "/*()");
					std::string tok = tokenizer.Next();
					if (tok.empty()) return result;
					result = ClassifyLineToken(tok);
					switch(result) {
						case kParseResult_Comment :
							// skip
							break;
						case kParseResult_VertexNormal :
						case kParseResult_VertexTexture :
						case kParseResult_Vertex :
							ReadVertex(tokenizer, result);
							break;
						case kParseResult_Face :
							ReadFace(tokenizer);
							break;
						default :
							// Not supported
							break;
					}
					return result;
				}

				void ReadVertex(Tokenizer &tok, kParseResult vertexType) {
					Vertex3D vertex;
					Vertex2D uv;
					switch (vertexType) {
						case kParseResult_Vertex:
							ParseVertex(tok, vertex.data, 3);
							coords.push_back(vertex);
							break;
						case kParseResult_VertexNormal :
							ParseVertex(tok, vertex.data, 3);
							vertexNormals.push_back(vertex);
							break;
						case kParseResult_VertexTexture :
							ParseVertex(tok, uv.data, 2);
							textureCoords.push_back(uv);
							break;
						default:
							break;
					}
				}

				void ReadFace(Tokenizer &tok) {
					Face *pFace = new Face();
					ParseFace(tok, pFace);
					faces.push_back(pFace);
				}

				void ReadStream(FILE *f) {
					char buffer[256];
					while(!feof(f)) {
						ReadLine(f, buffer, 256);
						ParseLine(std::string(buffer));
					}
				}

				void ReadFile(const char *filename) {
					FILE *f = fopen(filename, "rb");
					if (f!=NULL) {
						printf("Reading from file '%s'\n",filename);
						ReadStream(f);
						fclose(f);
					} else {
						fprintf(stderr,"Unable to open file: %s\n",filename);
					}
				}
			protected:
				char *ReadLine(FILE *f, char *buffer, int max) {
					fgets(buffer, max, f);
					return buffer;
				}
			public:
				std::vector<gnilk::engine::core::Face *> faces;
				std::vector<gnilk::engine::core::Vertex3D> coords;
				std::vector<gnilk::engine::core::Vertex3D> vertexNormals;
				std::vector<gnilk::engine::core::Vertex2D> textureCoords;

			};
		}
	}
}

