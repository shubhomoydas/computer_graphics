#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "osuGraphics.h"
#include "matlib.h"
#include "ObjLoader.h"
#include "utils.h"

void testTriangle1() {
	osuDiffuse(1.0f, 0.0f, 0.0f);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(3,-1,-3);
	osuVertex3f(1.5,-1,1.5);
	osuVertex3f(-3.5,-1,0);
	osuEnd();
}

void testTriangle2() {
	osuDiffuse(0.1f, 0.1f, 0.1f);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(2,4,-3.5);
	osuVertex3f(2.0,-1.0,1.5);
	osuVertex3f(-1,-3.5,-4);
	osuEnd();
}

void testTriangle3() {
	osuDiffuse(0.0f, 0.2f, 0.2f);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(-3.5,3.0,1.0);
	osuVertex3f(-3.5,-1.0,1.0);
	osuVertex3f(1.0,-1,-5.0);
	osuEnd();
}

void testTriangle4() {
	osuDiffuse(0.0f, 0.2f, 0.2f);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(-1.5,1.0,-1.5);
	osuVertex3f(1.5,1.0,-1.5);
	osuVertex3f(0.0,1.5,2.0);
	osuEnd();
}

void rayTest(CommandOptions& options) {

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

	if (!smd::VEC_ZERO.equals(options.l)) osuPointLight(lpos,options.intensity);
	osuAmbientLight(options.ambient);

	osuSpecular((float)options.specular.r, (float)options.specular.g, (float)options.specular.b, (float)options.specular.s);
	osuDiffuse((float)options.cr.r, (float)options.cr.g, (float)options.cr.b);

	getGraphicsManager().setRendererType(smd::RAY_TRACE);

	osuReflectionDepth(options.reflectionDepth);

	osuReflectance(0.7);
	testTriangle1();

	osuReflectance(0.9);
	testTriangle2();

	osuReflectance(1.0);
	testTriangle3();

	osuReflectance(1.0);
	testTriangle4();

	osuDiffuse(0, 1, 0);
	osuReflectance(0.2);
	osuSphere(-0.5, 0.0, -0.5, 1.0);

	osuDiffuse(0, 0, 1);
	osuReflectance(1.0);
	osuSphere(0.5, -0.5, 1.0, 0.5);

	getGraphicsManager().render();

	if (options.debug) std::cout << "Completed ray tracing..." << std::endl;

	if (0 && options.debug) {
		//std::cout << "M:" << std::endl << getGraphicsContext().getM() << std::endl;
		//std::cout << "Mvp Eye:" << std::endl << getGraphicsContext().getViewPoint().Mvp << std::endl;
		//std::cout << "Mcam Eye:" << std::endl << getGraphicsContext().getViewPoint().Mcam << std::endl;
		//std::cout << "Mvp Light:" << std::endl << getGraphicsContext().pointLights[0]->Mvp << std::endl;
		//std::cout << "Mproj Eye:" << std::endl << getGraphicsContext().getViewPoint().getMprojection() << std::endl;
	}

	smd::GraphicsContext& ctx = getGraphicsContext();
}

