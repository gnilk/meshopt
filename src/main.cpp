#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <math.h>

#include <GLFW/glfw3.h>
#include <OpenGl/glu.h>
#include "vec.h"

#include "gnilk/engine/core/Edge.hpp"
#include "gnilk/engine/core/Core.hpp"
#include "gnilk/engine/FaceHelper.hpp"
#include "gnilk/engine/Mesh.hpp"
#include "gnilk/Tokenizer.hpp"
#include "gnilk/formats/Wavefront.hpp"

#include <boost/unordered/unordered_map.hpp>


using namespace gnilk::engine::core;
using namespace gnilk::formats::wavefront;
using namespace gnilk;

static void calcFaceNormal(float *faceNormal, ObjReader *pObj, int v1, int v2, int v3);



// ------ end of importer

static void test_ParseFace(ObjReader *pObjReader) {
	int indexList[3];
	Face *pFace = new Face();
	Tokenizer tok("1//1 2//2 3//3");
	pObjReader->ParseFace(tok,pFace);
	assert((pFace->v[0]==1) && (pFace->v[1]=2) && (pFace->v[2]==3));
	assert((pFace->faceData & kFaceData_Vertex) == kFaceData_Vertex);
	assert(pFace->faceType == kFaceType_Tri);
}

static void test_ParseVertex(ObjReader *pObjReader) {
	Vertex3D vertex;
	Tokenizer tok1("1 2.1 3.0 4.0");
	pObjReader->ParseVertex(tok1, vertex.data, 3);
	pObjReader->coords.push_back(vertex);
	Tokenizer tok2("1 2.1 3.0 4.0");
	pObjReader->ParseVertex(tok2, vertex.data, 3);
	pObjReader->coords.push_back(vertex);
	pObjReader->DumpCoords();
}

static void test_ClassifyLine(ObjReader *pObjReader) {
	int res;
	res = pObjReader->ParseLine("# comment");
	assert(res == kParseResult_Comment);
	res = pObjReader->ParseLine("v 1 2 3");
	assert(res == kParseResult_Vertex);
	res = pObjReader->ParseLine("vn 1 2 3");
	assert(res == kParseResult_VertexNormal);
	res = pObjReader->ParseLine("vt 1 2 3");
	assert(res == kParseResult_VertexTexture);
	res = pObjReader->ParseLine("f 2 3 4");
	assert(res == kParseResult_Face);
}

static void test_LineParser(ObjReader *pObjReader) {
	int res;
	res = pObjReader->ParseLine("# comment");
	assert(res == kParseResult_Comment);
	res = pObjReader->ParseLine("v 1 2 3");
	assert(res == kParseResult_Vertex);
	res = pObjReader->ParseLine("f 2 3 4");
	assert(res == kParseResult_Face);
}

static void test_Tokenizer()
{
	std::string tmp;
	gnilk::Tokenizer tokenizer = gnilk::Tokenizer("apa/bpa/cpa","/*()");
	int opc = 0;
	while(!(tmp = tokenizer.Next()).empty()) {
		printf("%s\n",tmp.c_str());
		opc++;
	}
	printf("opc=%d\n",opc);
	assert (opc == 5);

	int idx;
	idx = gnilk::Tokenizer::Case("apa", "apa bpa cpa");
	assert(idx == 0);
	idx = gnilk::Tokenizer::Case("bpa", "apa bpa cpa");
	assert(idx == 1);
	idx = gnilk::Tokenizer::Case("cpa", "apa bpa cpa");
	assert(idx == 2);
	idx = gnilk::Tokenizer::Case("dpa", "apa bpa cpa");
	assert(idx == -1);


}

static void test_GetTok() {

	gnilk::Tokenizer tokenizer = gnilk::Tokenizer("apa/bpa/cpa");
	std::string tNext = tokenizer.Next();
	printf("tok: %s\n",tNext.c_str());
	tNext = tokenizer.Next();
	printf("tok: %s\n",tNext.c_str());
	tNext = tokenizer.Next();
	printf("tok: %s\n",tNext.c_str());
	tNext = tokenizer.Next();
	printf("tok: %s\n",tNext.c_str());
	tNext = tokenizer.Next();
	printf("tok: %s\n",tNext.c_str());
	tNext = tokenizer.Next();
	if (!tNext.empty())
		printf("tok: %s\n",tNext.c_str());
	else
		printf("eos\n");

}

static void test_Edge()
{
	std::vector<Vertex3D> vertices;
	Edge e1(5,6, vertices);
	Edge e2(6,5, vertices);
	Edge e3(2,1, vertices);

	assert(e1.Equals(e2));
	assert(!e1.Equals(e3));
}

static void test_unorderedmap()
{
	std::vector<Vertex3D> vertices;
	Edge e1(5,6,vertices);
	Edge e2(6,5,vertices);
	Edge e3(2,1,vertices);

	boost::unordered::unordered_map<int, Edge> map;

	map.insert(std::pair<int,Edge>(e1.HashCode(), e1));
	assert(map.find(e2.HashCode())!=map.end());
}


// -- end of testers
static void testGLFW(ObjReader *pObj);
int main (int argc, char * const argv[]) {
//	test_Heap();
	// test_Edge();
	//test_unorderedmap();
	//test_GetTok();
//	test_Tokenizer();
//	exit(1);
	char tmp[128];
	getcwd(tmp, 128);
	printf("CWD: %s\n",tmp);
	ObjReader *pObjReader = new ObjReader();
    double tStart, tEnd;
    tStart = glfwGetTime();
	//pObjReader->ReadFile("uh60.obj");
	pObjReader->ReadFile("lillen_256_faces.obj");
	//pObjReader->ReadFile("cube.obj");
	tEnd = glfwGetTime();
	printf("Read file in: %f sec\n",tEnd - tStart);
//	pObjReader->DumpCoords();
	//test_ClassifyLine(pObjReader);
	//test_StringCase();
	//test_GetTok();
//	test_ParseVertex(pObjReader);
//	test_ParseFace(pObjReader);
//	printf("Ok!\n");
	testGLFW(pObjReader);
    return 0;
}
static void calcFaceNormal(float *faceNormal, ObjReader *pObj, int v1, int v2, int v3) {
	float *pV1 = pObj->coords[v1].data;
	float *pV2 = pObj->coords[v2].data;
	float *pV3 = pObj->coords[v3].data;

	float vE1[3];
	float vE2[3];
	vSub(vE1, pV3, pV1);
	vSub(vE2, pV2, pV1);
	vCross(faceNormal, vE1, vE2);
}

static void calcFaceNormal(float *faceNormal, std::vector<gnilk::engine::core::Vertex3D> &coords, int v1, int v2, int v3) {
	float *pV1 = coords[v1].data;
	float *pV2 = coords[v2].data;
	float *pV3 = coords[v3].data;

	float vE1[3];
	float vE2[3];
	vSub(vE1, pV3, pV1);
	vSub(vE2, pV2, pV1);
	vCross(faceNormal, vE1, vE2);
	vNorm(faceNormal,faceNormal);
}

static void renderBoundingSphere(Vertex3D &mid, double radius)
{
	glColor3f(1,0.5,0.5);
	glBegin(GL_LINE_LOOP);
		for(int i=0;i<16;i++)
		{
			double x = radius * sin((double)i * M_PI/8.0);
			double z = radius * cos((double)i * M_PI/8.0);
			glVertex3f(mid.data[0]+x, mid.data[1], mid.data[2]+z);
		}
	glEnd();

	glColor3f(0.5,1,0.5);
	glBegin(GL_LINE_LOOP);
		for(int i=0;i<16;i++)
		{
			double y = radius * sin((double)i * M_PI/8.0);
			double z = radius * cos((double)i * M_PI/8.0);
			glVertex3f(mid.data[0], mid.data[1]+y, mid.data[2]+z);
		}
	glEnd();

	glColor3f(0.5,0.5,1);
	glBegin(GL_LINE_LOOP);
		for(int i=0;i<16;i++)
		{
			double x = radius * sin((double)i * M_PI/8.0);
			double y = radius * cos((double)i * M_PI/8.0);
			glVertex3f(mid.data[0]+x, mid.data[1]+y, mid.data[2]);
		}
	glEnd();
}

// TODO: Calculate edge tables
// 1: Define surface
// 2: Define edge
// 3:
static void renderVertexNormals(Triangle *tri, std::vector<Vertex3D> &coords, float radius)
{
	int v1 = tri->v[0];
	int v2 = tri->v[1];
	int v3 = tri->v[2];

	float tmp[3];
	float dummy[3];
	vIni(dummy,1,0,0);
	glBegin(GL_LINES);
	glVertex3fv(coords[v1].data);
	vAdd(tmp, coords[v1].data, vMul(tmp, tri->vn[0].data,radius));
	glVertex3fv(tmp);

	glVertex3fv(coords[v2].data);
	vAdd(tmp, coords[v2].data, vMul(tmp, tri->vn[1].data,radius));
	glVertex3fv(tmp);

	glVertex3fv(coords[v3].data);
	vAdd(tmp, coords[v3].data, vMul(tmp, tri->vn[2].data,radius));
	glVertex3fv(tmp);
	glEnd();
}
static void renderVertexNormalsForGroup(std::vector<FaceGroup *> &groups, std::vector<Triangle *> &triangles, std::vector<Vertex3D> &coords, float radius)
{
	glColor3f(0.5,0.5,0.5);
	radius = radius / 10.0;
    for(int g=0;g<groups.size();g++) {
    	int nTriangles = groups[g]->triangles.size();
    	for(int i=0;i<nTriangles;i++) {
    		Triangle *tri = triangles[groups[g]->triangles[i]];
    		renderVertexNormals(tri, coords, radius);
    	}
    }
}

static void renderGroups(std::vector<FaceGroup *> &groups, std::vector<Triangle *> &triangles, std::vector<Vertex3D> &coords)
{
	//glDisable(GL_LIGHTING);
    glBegin( GL_TRIANGLES );
    for(int g=0;g<groups.size();g++) {
    	int nTriangles = groups[g]->triangles.size();
    	for(int i=0;i<nTriangles;i++) {
    		Triangle *tri = triangles[groups[g]->triangles[i]];
    			int v1 = tri->v[0];
    			int v2 = tri->v[1];
    			int v3 = tri->v[2];

//    			float faceNormal[3];
//    			calcFaceNormal(faceNormal, coords, v1, v3, v2);
    			//vDup(faceNormal, normals[i].data);

    			double cr = (g&7) / 8.0f;
    			double cg = ((g+1)&7) / 8.0f;
    			double cb = ((g+2)&7) / 8.0f;
    			cr = 1.0 - cr;
    			glColor3f(cr,cg,cb);
    //			if (i&1)
    //			{
    //				glColor3f(0,1,0);
    //			} else
    //			{
    //				glColor3f(0,0,1);
    //			}
    			//glNormal3fv(tri->normal.data);
    			glNormal3fv(tri->vn[0].data);
    			glVertex3fv(coords[v1].data);
    			glNormal3fv(tri->vn[1].data);
    			glVertex3fv(coords[v2].data);
    			glNormal3fv(tri->vn[2].data);
    			glVertex3fv(coords[v3].data);
    	}

    }
	glEnd();
}
static void render(std::vector<Triangle *> &triangles, std::vector<gnilk::engine::core::Vertex3D> &coords)
{
	//glDisable(GL_LIGHTING);
    glBegin( GL_TRIANGLES );
	for(int i=0;i<triangles.size();i++) {
		Triangle *tri = triangles[i];
			int v1 = tri->v[0];
			int v2 = tri->v[1];
			int v3 = tri->v[2];

			float faceNormal[3];
			//calcFaceNormal(faceNormal, coords, v1, v3, v2);
			vDup(faceNormal, tri->normal.data);

//			if (i&1)
//			{
//				glColor3f(0,1,0);
//			} else
//			{
//				glColor3f(0,0,1);
//			}
			glNormal3fv(faceNormal);
			glVertex3fv(coords[v1].data);
			glVertex3fv(coords[v2].data);
			glVertex3fv(coords[v3].data);
	}
	glEnd();
}

static bool isESCPressed = false;
static bool isSPACEPressed = false;

static void glfwOnKey(GLFWwindow *window, int key, int scancode, int action, int mods) {
	switch(key) {
		case GLFW_KEY_ESCAPE:
			//glfwKeyStates[GOA_KEY_ESC] = action;
			if(action == GLFW_PRESS){
				isESCPressed = true;
			}
			return;
		case GLFW_KEY_SPACE:
			//glfwKeyStates[GOA_KEY_ESC] = action;
			if(action == GLFW_PRESS){
				isSPACEPressed = true;
			}
			return;
			break;
	}

	// if ((key >= 0) && (key < 256)) {
	// 	switch(action) {
	// 		case GLFW_PRESS :
	// 			if (glfwKeyStates[key]!=GLFW_PRESS) {
	// 				glfwKeyPress[key]=TRUE;
	// 			}
	// 			break;
	// 		case GLFW_RELEASE :
	// 			if (glfwKeyStates[key]!=GLFW_RELEASE) {
	// 				glfwKeyPress[key]=FALSE;
	// 			}
	// 			break;
	// 		case GLFW_REPEAT :
	// 			if (glfwKeyPress[key]==FALSE) {
	// 				glfwKeyPress[key]=TRUE;
	// 			}
	// 			break;
	// 	}
	// 	glfwKeyStates[key] = action;

	// 	// printf("key: %d, scancode: %d\n",key,scancode);
	// 	// glfwKeyStates[key] = action;
	// }
}


static void render(ObjReader *pObj)
{
        glBegin( GL_QUADS );
		for(int i=0;i<pObj->faces.size();i++) {
			Face *face = pObj->faces[i];
			if (face->faceType == kFaceType_Quad) {
				int v1 = face->v[0];
				int v2 = face->v[1];
				int v3 = face->v[2];
				int v4 = face->v[3];

				float faceNormal[3];
				calcFaceNormal(faceNormal, pObj, v1, v3, v2);

				glNormal3fv(faceNormal);
				glVertex3fv(pObj->coords[v1].data);
				glVertex3fv(pObj->coords[v2].data);
				glVertex3fv(pObj->coords[v3].data);
				glVertex3fv(pObj->coords[v4].data);

			}
		}

		//		glColor3f( 1.0f, 0.0f, 0.0f );
		//		glVertex3f( -5.0f, 0.0f, -4.0f );
		//		glColor3f( 0.0f, 1.0f, 0.0f );
		//		glVertex3f( 5.0f, 0.0f, -4.0f );
		//		glColor3f( 0.0f, 0.0f, 1.0f );
		//		glVertex3f( 0.0f, 0.0f, 6.0f );
        glEnd();

}


static void testGLFW(ObjReader *pObj) {
	int width, height, x;
    double t;

    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        exit( EXIT_FAILURE );
    }

    printf("PI: %f\n",M_PI);

    double tStart, tEnd;
    tStart = glfwGetTime();

    std::vector<gnilk::engine::core::Triangle *> triangles;
    std::vector<Vertex3D> &coords = pObj->coords;

    gnilk::engine::FaceHelper::FacesToTriangles(triangles, pObj->faces);

    gnilk::engine::Mesh *pMesh = NULL;
	gnilk::engine::Mesh *pMesh2 = NULL; 
    gnilk::engine::Mesh mesh_a(triangles, pObj->coords);

	int before = mesh_a.GetTriangles().size();
	pMesh = &mesh_a;
	//pMesh = mesh_a.Optimize();
//	pMesh = gnilk::engine::Mesh::CreateSphere(16, 6, 6);
//	pMesh = pMesh->Optimize();
	// pMesh = pMesh->Optimize();
	// pMesh = pMesh->Optimize();

//	pMesh2 = gnilk::engine::Mesh::CreateSphere(16, 6, 6);

	//pMesh = pMesh->Optimize();
	// for(int i=0;i<2000;i++) {
	// 	pMesh = pMesh->Optimize();
	// }
	printf("Before: %d\n", before);
	printf("After: %d\n", pMesh->GetTriangles().size());
	// for(int i=0;i<100;i++) {
	// 	pMesh = pMesh->Optimize();
	// }
	//pMesh = mesh_a.Optimize();

    //pMesh = &mesh_a;

    tEnd = glfwGetTime();
    printf("took: %f, (%f, %f)\n",tEnd - tStart, tStart, tEnd);
    printf("num edges: %d\n",pMesh->GetEdges().size());
    printf("Bounding sphere Radius: %f\n",pMesh->radius);

	width = 1280;
	height = 720;

    // Open a window and create its OpenGL context
	GLFWwindow *window = glfwCreateWindow(width, height, "meshopt", NULL, NULL);

//    if( !glfwOpenWindow( 640, 480, 8,8,8,0, 24,0, GLFW_WINDOW ) )
	if (window == NULL) {
        fprintf( stderr, "Failed to open GLFW window\n" );

        glfwTerminate();
        exit( EXIT_FAILURE );
    }

	glfwMakeContextCurrent(window);
	glfwSetKeyCallback(window, glfwOnKey);
	//glfwSetCursorPosCallback(window, glfwMousePosCallback);

    // Ensure we can capture the escape key being pressed below
    //glfwEnable( GLFW_STICKY_KEYS );

    // Enable vertical sync (on cards that support it)
    glfwSwapInterval( 1 );
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	  static GLfloat pos[4] = {5.f, 5.f, 10000.f, 0.f};
	glLightfv(GL_LIGHT0, GL_POSITION, pos);
	glEnable(GL_CULL_FACE);
	glEnable(GL_DEPTH_TEST);
	//glDisable(GL_CULL_FACE);

	do
    {	
		if (!isSPACEPressed) {
	        t = glfwGetTime();
		}
        // glfwGetMousePos( &x, NULL );

        // // Get window size (may be different than the requested size)
        // glfwGetWindowSize( &width, &height );

        // Special case: avoid division by zero below
        height = height > 0 ? height : 1;

		double mouse_xpos, mouse_ypos;
		glfwGetCursorPos(window, &mouse_xpos, &mouse_ypos);
        glViewport( 0, 0, width, height );

        // Clear color buffer to black
        glClearColor( 0.0f, 0.0f, 0.0f, 0.0f );
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Select and setup the projection matrix
        glMatrixMode( GL_PROJECTION );
        glLoadIdentity();
        gluPerspective( 65.0f, (GLfloat)width/(GLfloat)height, 1.0f, 1000.0f );

        float dist = pMesh->radius * 2.0f;
        // Select and setup the modelview matrix
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        gluLookAt( 0.0f, 0.0f, dist,    // Eye-position
				  0.0f, 0.0f, 0.0f,   // View-point
				  0.0f, 1.0f, 0.0f );  // Up-vector

        // Draw a rotating colorful triangle
        glTranslatef( 0.0f, 0.0f, 0.0f );
		glRotatef(t*10.0f,0,1,0);
		glRotatef(-t*30.0f,1,0,0);
		glDisable(GL_LIGHTING);
		glColor3f(1.0,0.0,0.0);
//		glBegin(GL_LINES);
//			for(int i=0;i<mesh.GetEdges().size();i++)
//			{
//				int a = mesh.GetEdges()[i]->GetV1();
//				int b = mesh.GetEdges()[i]->GetV2();
//				glVertex3fv(mesh.GetVertices()[a].data);
//				glVertex3fv(mesh.GetVertices()[b].data);
//			}
//		glEnd();
		glEnable(GL_LIGHTING);

  //        glBegin( GL_TRIANGLES );
		// for(int i=0;i<pObj->faces.size();i++) {
		// 	Face *face = pObj->faces[i];
		// 	if (face->faceType == kFaceType_Tri) {
		// 		int v1 = face->v[0]-1;
		// 		int v2 = face->v[1]-1;
		// 		int v3 = face->v[2]-1;


		// 		float faceNormal[3];
		// 		calcFaceNormal(faceNormal, pObj, v1, v3, v2);
		// 		glNormal3fv(faceNormal);
		// 		glVertex3fv(pObj->coords[v1].data);
		// 		glVertex3fv(pObj->coords[v2].data);
		// 		glVertex3fv(pObj->coords[v3].data);

		// 	}
		// }
		// // // glColor3f( 1.0f, 0.0f, 0.0f );
		// // // glVertex3f( -5.0f, 0.0f, -4.0f );
		// // // glColor3f( 0.0f, 1.0f, 0.0f );
		// // // glVertex3f( 5.0f, 0.0f, -4.0f );
		// // // glColor3f( 0.0f, 0.0f, 1.0f );
		// // // glVertex3f( 0.0f, 0.0f, 6.0f );
  //         glEnd();
		render(pMesh->GetTriangles(), pMesh->GetVertices());
		//renderGroups(mesh.GetGroups(), triangles, coords);
		//renderGroups(pMesh->GetGroups(), pMesh->GetTriangles(), pMesh->GetVertices());
		glDisable(GL_LIGHTING);
		glColor3f(1,0,0);
		glBegin(GL_LINES);
		for(auto e : pMesh->GetEdges()) {
			int a = e->GetV1();
			int b = e->GetV2();
			glVertex3fv(pMesh->GetVertices()[a].data);
			glVertex3fv(pMesh->GetVertices()[b].data);
		}
		glEnd();


		//renderBoundingSphere(pMesh->mid, pMesh->radius);


		//renderVertexNormalsForGroup(pMesh->GetGroups(), pMesh->GetTriangles(), pMesh->GetVertices(), pMesh->radius);

        // Swap bufferss
        glfwSwapBuffers(window);
		glfwPollEvents();
//		glfwGetFramebufferSize(window, &pxi, &pyi);


    } // Check if the ESC key was pressed or the window was closed
    while(!isESCPressed);

    // Close OpenGL window and terminate GLFW
    glfwTerminate();

}
