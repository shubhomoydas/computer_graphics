#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "osuGraphics.h"
#include "matlib.h"
#include "ObjLoader.h"
#include "utils.h"
#include "glRadiosity.h"

/******************************************************************************
Draw a polygon.
*****************************************************************************/
void simpleTest();
void blueCube();
void shadowTestScene(ObjModel& sphere, CommandOptions& options);
void drawObj(ObjModel& data, CommandOptions& options);
void objTestAny(CommandOptions& options);
void objTestShadows(CommandOptions& options);
void loadAndDrawObj(const char *fname, CommandOptions& options);
void testTriangle1();
void rayTest(CommandOptions& options);
void radiosityTest(CommandOptions& options);
void main(int argc , char **argv)
{
	CommandOptions options;
	if (!options.parseCommandLine(argc, argv)) {
		options.print_options();
		exit(-1);
	} else {
		options.print_options();
	}

	if (options.testCase >= 8) {
		//simplest tests
		switch (options.testCase) {
		case 8:
			smd::testMakePatch();
			break;
		default:
			break;
		}
		return;
	}

    int xsize = 400;
    int ysize = 400;
	const std::string root = "C:/smd/classes/CS551/project3";
	std::string filenameFlat = root + "/test.obj";
	std::string filenameSmooth = root + "/face.ws.obj";

    /* create a framebuffer */
    osuBeginGraphics (xsize, ysize);

    osuFlush();

   /* initialize the matrix stack */
    osuInitialize();
  
	setGraphicsDebug(options.debug, false);
	//osuColor3f(0.0,1.0,0.0);

	/* go to selected routine */

	switch (options.testCase) {
	case 1:
		simpleTest();
		break;
	case 2:
		if (options.depthTest) osuEnable(OSU_DEPTH_TEST);
		blueCube();
		break;
	case 3:
		if (options.depthTest) osuEnable(OSU_DEPTH_TEST);
		objTestAny(options);
		break;
	case 4:
		if (options.depthTest) osuEnable(OSU_DEPTH_TEST);
		objTestShadows(options);
		break;
	case 5:
		rayTest(options);
		break;
	case 6:
	case 7:
		radiosityTest(options);
		break;
	default:
		fprintf (stderr, "Please use a valid test case.\n");
		exit (-1);
	}

	osuFlush();
	osuWaitOnEscape();
	osuEndGraphics();
}

void objTestShadows(CommandOptions& options) {

	osuPerspective(90.0, 1.0, 1000);
	//osuOrtho(-2, 2, -2, 2, 1.0, 1000);
	osuClear(0, 0, 0);	
	
	std::string spherepath = "C:/smd/classes/CS551/project3/sphere.obj";

	float from[3]={(float)options.e.x, (float)options.e.y, (float)options.e.z};//{2.0,2.0,2.0};
	float at[3] = {(float)options.at.x, (float)options.at.y, (float)options.at.z};//{0.0,0.0,0.0};
	float up[3] = {0.0, 1.0, 0.0};

	osuLookat (from, at, up);

	osuShadeModel(options.shadeModel);

	float lpos[3]={(float)options.l.x, (float)options.l.y, (float)options.l.z}; //{-3.0, 3.0, 3.0};

	if (!smd::VEC_ZERO.equals(options.l)) osuPointLight(lpos,0.7f);
	osuAmbientLight(0.2f);

	getGraphicsManager().setNormMul(options.normMul);

	ObjModel sphere;
	ObjLoader LoaderClass;

	LoaderClass.LoadObj(spherepath.c_str());
	sphere = LoaderClass.ReturnObj();

	if (options.shadows) {
		// First pass to compute shadow map
		int nlights = getGraphicsManager().getNumPointLights();
		getGraphicsManager().setRenderMode(smd::COMPUTE_SHADOW);
		for (int i = 1; i <= nlights; i++) {
			getGraphicsManager().setViewPoint(i); // set to light source view
			shadowTestScene(sphere, options);
		}
		if (options.debug) std::cout << "Completed first pass..." << std::endl;
	}

	// Second pass to actually render image
	getGraphicsManager().setViewPointToEye();
	getGraphicsManager().setRenderMode(smd::COLOR_AND_SHADOW);
	if (options.shadows && options.shadowsOnly) {
		getGraphicsManager().setRenderMode(smd::RENDER_SHADOW);
	}

	shadowTestScene(sphere, options);
	if (options.debug) std::cout << "Completed second pass..." << std::endl;

	smd::GraphicsContext& ctx = getGraphicsContext();
	std::cout << "#ZeroNormals: " << ctx.zeroNormals << "; priorVertexNormals: " << ctx.priorVertexNormals << \
		"; computedVertexNormals: " << ctx.computedVertexNormals << std::endl;
}

void shadowTestScene(ObjModel& sphere, CommandOptions& options) {

	osuDiffuse((float)options.cr.r, (float)options.cr.g, (float)options.cr.b);

	//Specular white color
	osuSpecular((float)options.specular.r, (float)options.specular.g, (float)options.specular.b, (float)options.specular.s);

	if (0) {
	osuDiffuse(1.0, 0.0, 0.0);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(0,0,0);
	osuVertex3f(1,-1,1);
	osuVertex3f(0,-1,0);
	osuEnd();
	}

	if (0) {
	osuDiffuse(0.0, 1.0, 0.0);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(1,-1,-1);
	osuVertex3f(1,-1,1);
	osuVertex3f(-1,-1,1);
	osuEnd();
	}

	testTriangle1();

	if (1) {
	osuPushMatrix();
	osuTranslate(-0.5, 0.0, -0.5);
	drawObj(sphere, options);
	osuPopMatrix();
	}

}

void objTestAny(CommandOptions& options)
{
	osuPerspective(90.0, 1.0, 1000);
	osuClear(0, 0, 0);	
	
	float from[3]={(float)options.e.x, (float)options.e.y, (float)options.e.z};//{2.0,2.0,2.0};
	float at[3] = {(float)options.at.x, (float)options.at.y, (float)options.at.z};//{0.0,0.0,0.0};
	float up[3] = {0.0, 1.0, 0.0};

	osuLookat (from, at, up);

	//Diffuse blue color
	//osuDiffuse(0.0, 0.0 , 1.0 );
	osuDiffuse((float)options.cr.r, (float)options.cr.g, (float)options.cr.b);

	//Specular white color
	osuSpecular((float)options.specular.r, (float)options.specular.g, (float)options.specular.b, (float)options.specular.s);

	osuShadeModel(options.shadeModel);

	float lpos[3]={(float)options.l.x, (float)options.l.y, (float)options.l.z}; //{-3.0, 3.0, 3.0};

	if (!smd::VEC_ZERO.equals(options.l)) osuPointLight(lpos,0.7f);
	osuAmbientLight(0.2f);

	getGraphicsManager().setNormMul(options.normMul);

	loadAndDrawObj(options.file.c_str(), options);

}

void drawObj(ObjModel& data, CommandOptions& options)
{
	bool smooth = (options.shadeModel == OSU_SMOOTH);
	for(int i = 0; i < data.NumTriangle; i++)  {						

		osuBegin(OSU_TRIANGLE);

		//For flat shading, only make one call to osuNormal in the beginning, and use the following
		// to access the faceNormal
		if (!smooth) {
			osuNormal3f(data.TriangleArray[i].faceNormal[0], 
					data.TriangleArray[i].faceNormal[1], 
					data.TriangleArray[i].faceNormal[2]); 
		}
		if(smooth) {
			osuNormal3f(data.NormalArray[data.TriangleArray[i].Vertex[0]].X, 
					data.NormalArray[data.TriangleArray[i].Vertex[0]].Y, 
					data.NormalArray[data.TriangleArray[i].Vertex[0]].Z); 
		}
		osuVertex3f(data.VertexArray[data.TriangleArray[i].Vertex[0]].X, 
				data.VertexArray[data.TriangleArray[i].Vertex[0]].Y, 
				data.VertexArray[data.TriangleArray[i].Vertex[0]].Z);
	
		if(smooth) {
			osuNormal3f(data.NormalArray[data.TriangleArray[i].Vertex[1]].X, 
					data.NormalArray[data.TriangleArray[i].Vertex[1]].Y, 
					data.NormalArray[data.TriangleArray[i].Vertex[1]].Z); 
		}
		osuVertex3f(data.VertexArray[data.TriangleArray[i].Vertex[1]].X, 
				data.VertexArray[data.TriangleArray[i].Vertex[1]].Y, 
				data.VertexArray[data.TriangleArray[i].Vertex[1]].Z);

		if(smooth) {
			osuNormal3f(data.NormalArray[data.TriangleArray[i].Vertex[2]].X, 
					data.NormalArray[data.TriangleArray[i].Vertex[2]].Y, 
					data.NormalArray[data.TriangleArray[i].Vertex[2]].Z); 
		}
		osuVertex3f(data.VertexArray[data.TriangleArray[i].Vertex[2]].X, 
				data.VertexArray[data.TriangleArray[i].Vertex[2]].Y, 
				data.VertexArray[data.TriangleArray[i].Vertex[2]].Z);

		osuEnd();

	}
}

void loadAndDrawObj(const char *fname, CommandOptions& options)
{
	ObjModel data;
	ObjLoader LoaderClass;

	LoaderClass.LoadObj(fname);
	data = LoaderClass.ReturnObj();

	std::cout << "Data loaded..." << std::endl;
	std::cout << "#Triangles: " << data.NumTriangle << 
		"; #Normals: " << data.NumNormal << 
		"; #Vertices: " << data.NumVertex << std::endl;
	//exit(-1);

	osuPushMatrix();

	osuScale(options.scale, options.scale, options.scale);

	drawObj(data, options);

	osuPopMatrix();

	smd::GraphicsContext& ctx = getGraphicsContext();
	std::cout << "#ZeroNormals: " << ctx.zeroNormals << "; priorVertexNormals: " << ctx.priorVertexNormals << \
		"; computedVertexNormals: " << ctx.computedVertexNormals << std::endl;
}

void simpleTest()
{
	osuPerspective(90.0, 1.0, -1000.);

	float from[3]={3.0,0.0,3.0};
	float at[3] = {0.0,0.0,-8.0};
	float up[3] = {0.0, 1.0, 0.0};

	osuLookat (from, at, up);

	osuClear(0,0,0);

	osuDiffuse(0.0, 0.0 , 1.0 );
	osuSpecular(1.0, 1.0, 1.0, 1.0);

	float lpos[3]={0.0, 1.5, 5.0};

	osuPointLight(lpos,1.0);
	osuAmbientLight(0.4f);

	//YOU MUST CONVERT THIS TO TWO TRIANGLES!!!
	osuBegin(OSU_POLYGON);
	osuVertex3f(-4.5, -1.75, -5.5);
	osuVertex3f(-4.5, 1.75, -5.5);
	osuVertex3f(4.5, 1.75, -5.5);
	osuVertex3f(4.5, -1.75, -5.5);
	osuEnd();


	osuDiffuse(1.0, 0.0 , 0.0 );
	osuSpecular(1.0, 1.0, 1.0, 1.0);

	osuColor3f(0.0,0.0,1.0);
	osuBegin(OSU_POLYGON);
	osuVertex3f(0.0, -1.75, -2.5);
	osuVertex3f(0.0, 1.75, -2.5);
	osuVertex3f(0.0, 1.75, -7.5);
	osuVertex3f(0.0, -1.75, -7.5);
	osuEnd();

}

void blueCube()
{
	osuPerspective(40, 6.5, 100);

	float from[3]={5.0,5.0,5.0};
	float at[3] = {0.0,0.0,0.0};
	float up[3] = {0.0, 1.0, 0.0};

	osuLookat (from, at, up);

	osuClear(0, 0, 0);
	osuClearZ();

	osuDiffuse(0.0, 0.0 , 1.0 );
	osuSpecular(1.0, 1.0, 1.0, 1.0);

	float lpos[3]={3.0, 1.5, 5.0};
	float dir[3] = {0.0, -1.0, 0.0};

	osuPointLight(lpos,0.5);
	osuDirectionalLight(dir,0.5);
	osuAmbientLight(0.4f);

	//YOU MUST CONVERT THESE TO USE TRIANGLES!!!
	//back
	
	osuBegin(OSU_POLYGON);
	osuVertex3f( -1, -1, -1);
	osuVertex3f(  1, -1, -1);
	osuVertex3f(  1,  1, -1);
	osuVertex3f( -1,  1, -1);
	osuEnd();
	
	//right
	osuBegin(OSU_POLYGON);
	osuVertex3f(  1, -1, -1);
	osuVertex3f(  1 ,-1,  1);
	osuVertex3f ( 1 , 1,  1);
	osuVertex3f ( 1 , 1, -1);
	osuEnd();
	
	
	//front
	osuBegin(OSU_POLYGON);
	osuVertex3f( -1, -1,  1);
	osuVertex3f( -1,  1,  1);
	osuVertex3f(  1,  1,  1);
	osuVertex3f(  1, -1,  1);
	osuEnd();
	
	//top
	osuDiffuse(1.0, 0.0 , 0.0 );
	osuSpecular(1.0, 1.0, 1.0, 1.0);

	osuBegin(OSU_POLYGON);
	osuVertex3f( -1,  1, -1);
	osuVertex3f(  1,  1, -1);
	osuVertex3f(  1,  1,  1);
	osuVertex3f( -1,  1,  1);
	osuEnd();
	
	
	//bottom
	osuBegin(OSU_POLYGON);
	osuVertex3f( -1, -1, -1);
	osuVertex3f( -1, -1,  1);
	osuVertex3f(  1, -1,  1);
	osuVertex3f(  1 ,-1, -1);
	osuEnd();

	
	//left
	osuDiffuse(0.0, 1.0 , 0.0 );
	osuSpecular(1.0, 1.0, 1.0, 1.0);
	osuBegin(OSU_POLYGON);
	osuVertex3f( -1, -1, -1);
	osuVertex3f( -1,  1 ,-1);
	osuVertex3f( -1,  1,  1);
	osuVertex3f( -1, -1,  1);
	osuEnd();
	
}

