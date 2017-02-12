#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include "osuGraphics.h"
#include "glDatatypes.h"

#ifndef _CS551_UTILS_H_
#define _CS551_UTILS_H_

typedef struct _CommandOptions {
	std::string file;
	int testCase;
	bool depthTest;
	OSUShadeModel shadeModel;
	double scale;
	bool shadows;
	bool shadowsOnly;
	bool debug;
	double normMul;
	float ambient;
	smd::Color cr;
	smd::Color specular;
	smd::Vector e;
	smd::Vector at;
	smd::Vector l;
	float intensity;
	int reflectionDepth;
	int iters; // number of iterations for solution (radiosity)
	int nparts; // Number of patches (radiosity) = nparts x nparts
	int ncontour; // Contour integral (radiosity) = ncontour x ncontour
	bool rndcolor; // assign random colors to patches (radiosity).
	bool smoothcolor; // smoothen colors in patches for each face (radiosity).
	bool occlusion; // check occlusion for each patch (radiosity).
	int scene; // test scene (radiosity).
	bool patchesOnly; // render only patches without computing radiosity.
	bool loadFF; // load pre-computed Form Factor (radiosity)
	std::string FFfile; // pre-computed Form Factor file (radiosity)
	int solver; // Radiosity solver (radiosity).
	bool addPRAmbient; // Add ambient light in Progressive Refinement (Radiosity)
	int argc;
	char **argv;
	_CommandOptions(): \
		file(std::string()), testCase(1), depthTest(true), shadeModel(OSU_FLAT), scale(1.0), \
		shadows(false), shadowsOnly(false), debug(false), normMul(1.0), ambient(0.2f), intensity(0.7f), reflectionDepth(1),
		iters(0), nparts(4), ncontour(4), rndcolor(false), smoothcolor(false), occlusion(false), scene(1), patchesOnly(false),
		loadFF(false), FFfile(std::string()), solver(0), addPRAmbient(false),
		argc(0), argv(0) {}
	void print_options() {
		std::cout << "file=" << file << std::endl;
		std::cout << "testCase=" << testCase << std::endl;
		std::cout << "depthTest=" << (depthTest ? "true" : "false") << std::endl;
		std::cout << "shadeModel=" << shadeModel << std::endl;
		std::cout << "scale=" << scale << std::endl;
		std::cout << "specular=" << specular << std::endl;
		std::cout << "normMul=" << normMul << std::endl;
		std::cout << "ambient=" << ambient << std::endl;
		std::cout << "shadows=" << (shadows ? "true" : "false") << std::endl;
		std::cout << "shadowsOnly=" << (shadowsOnly ? "true" : "false") << std::endl;
		std::cout << "debug=" << (debug ? "true" : "false") << std::endl;
		std::cout << "diffuse (Diffuse color)=" << cr << std::endl;
		std::cout << "e (Eye coords)=" << e << std::endl;
		std::cout << "at (Look at)=" << at << std::endl;
		std::cout << "l (Point light)=" << l << std::endl;
		std::cout << "intensity (Light intensity)=" << intensity << std::endl;
		std::cout << "Reflection depth=" << reflectionDepth << std::endl;
		std::cout << "iters (radiosity)=" << iters << std::endl;
		std::cout << "nparts (radiosity)=" << nparts << std::endl;
		std::cout << "ncontour (radiosity)=" << ncontour << std::endl;
		std::cout << "rndcolor (radiosity)=" << (rndcolor ? "true" : "false") << std::endl;
		std::cout << "smoothcolor (radiosity)=" << (smoothcolor ? "true" : "false") << std::endl;
		std::cout << "occlusion (radiosity)=" << (occlusion ? "true" : "false") << std::endl;
		std::cout << "scene (radiosity)=" << scene << std::endl;
		std::cout << "patchesOnly (radiosity)=" << (patchesOnly ? "true" : "false") << std::endl;
		std::cout << "loadFF (radiosity)=" << (loadFF ? "true" : "false") << std::endl;
		std::cout << "FFfile=" << FFfile << std::endl;
		std::cout << "solver (radiosity)=" << solver << " (0 - Gauss-Seidel; 1 - Progressive Refinement)" << std::endl;
		std::cout << "addPRAmbient (radiosity)=" << (addPRAmbient ? "true" : "false") << std::endl;
	}
	bool parseCommandLine(int argc, char **argv);
} CommandOptions;

template<class T>
void printArray(T* array, int len) {
	for (int i = 0; i < len; i++) {
		printf("%4.3f ",(float)array[i]);
	}
	printf("\n");
}

#endif
