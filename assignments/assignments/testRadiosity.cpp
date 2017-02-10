#include "glRadiosity.h"

void testRadiosityTriangle1();
void testRadiosityTriangle2();
void rBox(float i1, float i2);
void rColumn(double length, double width, double dist=10);
void rSquare();
void scene0(float i);
void scene1(float i);
void scene2(float i);
void scene3(float i);
void scene4(float i);
void scene5(float i);
void renderScene(CommandOptions& options);

void radiosityTest(CommandOptions& options) {

	osuEnable(OSU_DEPTH_TEST);

	osuPerspective(90.0, 1.0, 10000);
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

	//smd::testGaussSeidel();

	if (options.testCase == 7) {
		getGraphicsManager().setRendererType(smd::SCAN_CONVERT);
		renderScene(options);
		return;
	} else {
		getGraphicsManager().setRendererType(smd::RADIOSITY);
		renderScene(options);
		std::list<smd::Shape*>& shapes = getGraphicsManager().getShapes();
		std::vector<smd::Patch*> patches;

		//smd::subDivideByMaxLength(shapes, 10, patches);
		smd::subDivideByNumParts(shapes, options.nparts, patches);

		std::cout << "Made " << patches.size() << " patches with (" << 
			options.nparts << " x " << options.nparts << ") per-face. #faces=" << shapes.size() << std::endl;

		smd::cachePatchVerticesAndInitNormals(patches);

		if(0 && options.nparts == 1) {
			for (unsigned int i = 0; i < patches.size(); i++) {
				for (unsigned int j = i+1; j < patches.size(); j++) {
					if (smd::isPatchVisible(i, j, patches)) {
						std::cout << "Patch " << i << " and " << j << " are visible" << std::endl;
					}
				}
			}
			return;
		}

		if (!options.patchesOnly) {
			smd::Matrix* FF = 0;
			if (options.loadFF) {
				FF = smd::loadCSV(options.FFfile);
				std::cout << "Loaded FF from " << options.FFfile << "; rows=" << FF->rows << "; cols=" << FF->cols << std::endl;
			} else {
				int **HID = 0;
				if (options.occlusion) {
					// compute HID
					HID = smd::getHID(patches);
					std::cout << "Computed HID..." << std::endl;
				}
				FF = new smd::Matrix(patches.size());
				smd::computeSimpleFormFactors(patches, *FF, options.ncontour, HID);
				if (!options.FFfile.empty()) smd::saveData(*FF, options.FFfile);
				std::cout << "Computed Form Factors with " << patches.size() << " patches and each patch (" <<
					options.ncontour << " x " << options.ncontour << ") sub units." << std::endl;
				//std::cout << "Form Factor Mat:" << std::endl << FF;
				if (HID) smd::destroy2DIntArray(HID, patches.size());
			}
			if (options.iters > 0) {
				smd::RadiosityLinearSystemRGB *rad = new smd::RadiosityLinearSystemRGB(patches, *FF, options.solver, options.addPRAmbient);
				std::cout << "Created RGB Linear systems..." << std::endl;
				rad->solve(options.iters);
				std::cout << "Computed solution..." << std::endl;
				rad->loadSolutions(patches);
				std::cout << "Modified patch colrs..." << std::endl;
				rad->printAmibient();
				delete rad;
			}
			if (FF) delete FF;
		}

		if (options.rndcolor) smd::setFaceColorsToRandom(patches);

		smd::setVertexColorsToFace(patches);

		std::map<int, std::vector<smd::Patch*>*>* pmap = smd::groupPatches(patches);
		//smd::printPatches(pmap);
		std::map<int, std::map<int, std::vector<smd::Patch*>* >* >* vmap = smd::groupVertices(pmap);
		std::cout << "Vertexes Grouped..." << std::endl;
		//smd::printVertexGroups(vmap);

		if (options.smoothcolor) smd::averageVertexColorsForGroup(vmap);

		std::vector<smd::Shape*> v;
		std::list<smd::Shape*> s;
		smd::vectorToList<smd::Shape*>(smd::typecastPatchToShapeVector(patches, v), s);
		getGraphicsManager().render(s, getGraphicsContext());

		if (options.debug) std::cout << "Completed radiosity..." << std::endl;
	
		// Cleanup

		smd::release_map_elements(pmap);
		smd::release_map_elements_with_map<int, std::map<int, std::vector<smd::Patch*>* > >(vmap);

		delete pmap;
		delete vmap;
	
		smd::release_vector_elements<smd::Patch*>(&patches);
	}

}

void renderScene(CommandOptions& options) {
	switch (options.scene) {
	case 0:
		scene0(options.intensity);
		break;
	case 1:
		scene1(options.intensity);
		break;
	case 2:
		scene2(options.intensity);
		break;
	case 3:
		scene3(options.intensity);
		break;
	case 4:
		scene4(options.intensity);
		break;
	case 5:
		scene5(options.intensity);
		break;
	default:
		std::cout << "Enter a valid scene" << std::endl;
		std::exit(-1);
	}
}

void scene0(float i) {
	osuPushMatrix();
	osuTranslate(0,0,0);
	rColumn(40.0,20.0);
	osuPopMatrix();

	if(0) {
		// Left face
		osuDiffuse(1.0f, 0.0f, 0.0f);
		osuPushMatrix();
		osuTranslate(-50,0,-100);
		osuScale(50,50,50);
		osuRotate(90,0,1,0);
		rSquare();
		osuPopMatrix();
	}
}

void scene1(float i) {
	rBox(i,i);
}

void scene2(float i) {
	rBox(i,i);
	osuPushMatrix();
	osuTranslate(-20,0,-65);
	rColumn(20,10);
	osuPopMatrix();
}

void scene3(float i) {
	rBox(i,0);
	osuPushMatrix();
	osuTranslate(-20,0,-65);
	rColumn(20,10);
	osuPopMatrix();
}

void scene4(float i) {
	rBox(0,i);
	osuPushMatrix();
	osuTranslate(-20,0,-65);
	rColumn(20,10);
	osuPopMatrix();
}

void scene5(float i) {
	rBox(i,0);
	osuPushMatrix();
	osuTranslate(-20,0,-65);
	rColumn(20,10);
	osuPopMatrix();

	osuPushMatrix();
	osuTranslate(20,0,-75);
	rColumn(5,5);
	osuPopMatrix();
}

void rColumn(double length, double width, double dist) {

	osuEmittance(0,0,0); // faces do not emit light

	// Front face - should be the emitting surface
	osuDiffuse(0.54f, 0.54f, 0.54f);
	//osuDiffuse(1.0f, 0.0f, 0.0f);
	osuPushMatrix();
	osuTranslate(0,0,width/2);
	osuScale(length/2,50,1);
	rSquare();
	osuPopMatrix();

	// Left face
	osuDiffuse(0.54f, 0.54f, 0.54f);
	//osuDiffuse(0.0f, 0.0f, 1.0f);
	osuPushMatrix();
	osuTranslate(-length/2,0,0);
	osuRotate(-90,0,1,0);
	osuScale(width/2,50,1);
	rSquare();
	osuPopMatrix();

	// Back face
	osuDiffuse(0.54f, 0.54f, 0.54f);
	//osuDiffuse(0.0f, 1.0f, 0.0f);
	osuPushMatrix();
	osuTranslate(0,0,-width/2);
	osuRotate(180,0,1,0);
	osuScale(length/2,50,1);
	rSquare();
	osuPopMatrix();

	// Right face
	osuDiffuse(0.54f, 0.54f, 0.54f);
	//osuDiffuse(1.0f, 1.0f, 0.0f);
	osuPushMatrix();
	osuTranslate(length/2,0,0);
	osuRotate(90,0,1,0);
	osuScale(width/2,50,1);
	rSquare();
	osuPopMatrix();

}

void rBox(float i1, float i2) {

	// same configuration as the Cornell box

	// Front face - should be the emitting surface
	osuDiffuse(0.54f, 0.54f, 0.54f);
	/*{
		osuPushMatrix();
		osuTranslate(0,0,-50);
		osuScale(50,50,50);
		osuRotate(180,0,1,0);
		rSquare();
		osuPopMatrix();
	}*/
	{ // comprises of two segments
		osuEmittance(i1,i1,i1);
		osuPushMatrix();
		osuTranslate(-25,0,-50);
		osuScale(25,50,50);
		osuRotate(180,0,1,0);
		rSquare();
		osuPopMatrix();

		osuEmittance(i2,i2,i2);
		osuPushMatrix();
		osuTranslate(25,0,-50);
		osuScale(25,50,50);
		osuRotate(180,0,1,0);
		rSquare();
		osuPopMatrix();
	}

	osuEmittance(0,0,0); // Other faces do not emit light

	// Back face
	osuDiffuse(0.84f, 0.84f, 0.84f);
	osuPushMatrix();
	osuTranslate(0,0,-150);
	osuScale(50,50,50);
	rSquare();
	osuPopMatrix();

	// Left face
	osuDiffuse(1.0f, 0.0f, 0.0f);
	osuPushMatrix();
	osuTranslate(-50,0,-100);
	osuScale(50,50,50);
	osuRotate(90,0,1,0);
	rSquare();
	osuPopMatrix();

	// Right face
	osuDiffuse(0.0f, 0.0f, 1.0f);
	osuPushMatrix();
	osuTranslate(50,0,-100);
	osuScale(50,50,50);
	osuRotate(90,0,-1,0);
	rSquare();
	osuPopMatrix();

	// Bottom face
	osuDiffuse(0.54f, 0.54f, 0.54f);
	osuPushMatrix();
	osuTranslate(0,-50,-100);
	osuScale(50,50,50);
	osuRotate(90,-1,0,0);
	rSquare();
	osuPopMatrix();

	// Top face
	osuDiffuse(0.84f, 0.84f, 0.84f);
	osuPushMatrix();
	osuTranslate(0,50,-100);
	osuScale(50,50,50);
	osuRotate(90,1,0,0);
	rSquare();
	osuPopMatrix();

}

void rSquare() {
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_POLYGON);
	osuVertex3f(1,1,0);
	osuVertex3f(1,-1,0);
	osuVertex3f(-1,-1,0);
	osuVertex3f(-1,1,0);
	osuEnd();
}

void testRadiosityTriangle1() {
	osuDiffuse(1.0f, 0.0f, 0.0f);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(30,-10,-30);
	osuVertex3f(15,-10,15);
	osuVertex3f(-35,-10,0);
	osuEnd();
}

void testRadiosityTriangle2() {
	osuDiffuse(0.0f, 1.0f, 0.0f);
	osuNormal3f(0, 0, 0);
	osuBegin(OSU_TRIANGLE);
	osuVertex3f(30,-10,-30);
	osuVertex3f(-35,-10,0);
	osuVertex3f(-2.5,20,-15);
	osuEnd();
}
