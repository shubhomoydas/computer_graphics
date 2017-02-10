#ifndef OSUGLGRAPHICSMANAGER
#define OSUGLGRAPHICSMANAGER

#include <iostream>
#include <stdio.h>
#include <list>
#include "osuGraphics.h"
#include "glDatatypes.h"

#define round(x) (int)(x+0.5)
#define PI 3.1415926535897
#define floatToDouble3(fl,db) {db[0]=fl[0]; db[1]=fl[1]; db[2]=fl[2];}

#define PRINT_OBJ(msg,obj,out) { out << msg; obj.print(out); out << std::endl; }

namespace smd {

	Vector& matrix_X_vector(const Matrix& m, const Vector& v, Vector& dest);
	double* matrix_X_vecArray(const Matrix& m, const double* v1, double* v2);

	Matrix& translate(double x, double y, double z, Matrix& m, Matrix& temp, Matrix& dest);
	Matrix& scale(double x, double y, double z, Matrix& m, Matrix& temp, Matrix& dest);
	Matrix& lookAt(const Vector& e, const Vector& u, const Vector& v, const Vector& w, Matrix& dest);

	Quaternion& getRotationQuaternion(double theta, double x, double y, double z, Quaternion& q);
	Matrix& getQuaternionToRotationMatrix(Quaternion& q, Matrix& m);

	void setMvp(int w, int h, Matrix& Mvp);
	void setMorth(double l, double r, double b, double t, double n, double f, Matrix& Morth);
	void setMpers(double l, double r, double b, double t, double n, double f, Matrix& Mpers);

	// Assume that p1 -> p2 -> p3 is clockwise
	Vector getNormal(const Point& p1, const Point& p2, const Point& p3);

	template<class T>
	std::vector<T>& listToVector(std::list<T>& l, std::vector<T>& v) {
		for (std::list<T>::iterator it = l.begin(); it != l.end(); it++) {
			v.push_back(*it);
		}
		return v;
	}

	template<class T>
	std::list<T>& vectorToList(std::vector<T>& v, std::list<T>& l) {
		for (unsigned int i=0; i < v.size(); i++) {
			l.push_back(v[i]);
		}
		return l;
	}

	enum RenderMode { COLOR_AND_SHADOW, COLOR, COMPUTE_SHADOW, RENDER_SHADOW, NONE };
	enum RendererType { SCAN_CONVERT, RAY_TRACE, RADIOSITY };

	class Ray {
	public:
		Vector d; // direction
		Vector o; // origin
		Ray() : d(VEC_ZERO), o(VEC_ZERO) {}
		Ray(Vector& d, Vector& o) : d(d), o(o) {
		}
		void getPoint(double t, Point& dest) {
			dest.set(d).mul(t).add(o);
		}
	};

	class HitData {
	public:
		double t[2]; // max two hits supported: 2 - spheres, 1 - triangles
		int nhits;
		double alpha, beta, gamma; // barycentric coordinates for triangles
	};

	class Zbuffer {
	public:
		Matrix* buffer;
		Zbuffer() : buffer(0) {}
		void clearBuffer() {
			if (buffer != 0) {
				buffer->setVal(-1e+300); // some very large number (actually, should be -1 since that is our far-plane)
				//std::cout << "clearBuffer()" << std::endl;
			} else {
				std::cout << "clearBuffer :: Z-buffer not initialized. Call enableDepthTest first." << std::endl;
			}
		}
		void initBuffer() {
			destroyBuffer();
			int w, h;
			osuGetFramebufferSize(&w, &h);
			buffer = new Matrix(h, w);
			clearBuffer();
		}
		void destroyBuffer() {
			if (buffer != 0) {
				delete buffer;
				buffer = 0;
			}
		}
		virtual ~Zbuffer() {
			destroyBuffer();
		}
	};

	class ViewPoint {
	public:
		Vector e; // location
		Vector u; // x-axis for eye coordinates
		Vector v; // up direction
		Vector w; // Gaze direction w == at - e

		double l, r, b, t; // View region - left, right, bottom, top

		Vector at; // Looking at.
		int zmode;
		Zbuffer zbuffer;
		int matrixMode; // 0 - Orthogonal, 1 - Perspective
		int ny; // Window height
		int nx; // Window width
		double nearp;
		double farp;
		Matrix Mcam;
		Matrix Mvp;
		Matrix iMcam;
		Matrix iMvp;
		Matrix iMproj;
		Matrix *Morth, *Mpers;
		Matrix *Mproj[2]; // Two projections being supported
		ViewPoint() : e(VEC_ZERO), u(VEC_ZERO), v(VEC_ZERO), w(VEC_ZERO), at(VEC_ZERO), 
			zmode(-1), zbuffer(Zbuffer()), matrixMode(0), 
			nearp(-1), farp(1),
			Mcam(Matrix(4)), Mvp(Matrix(4)), iMcam(Matrix(4)), iMvp(Matrix(4)), 
			iMproj(Matrix(4)), Morth(0), Mpers(0) {
			Morth = new Matrix(4);
			Mpers = new Matrix(4);
			Mproj[0] = Morth;
			Mproj[1] = Mpers;
			initializeMatrices(false);
		}
		void setNearFar(double nearp, double farp) {
			this->nearp = nearp;
			this->farp = farp;
		}
		void lookAt(double from[3], double at[3], double up[3]);
		void lookAt(const Vector& from, const Vector& at, const Vector& v);
		// e - eye loc, v - up direction, w - gaze direction, v - orthogonal to  u&v
		//void lookAt(const Vector& e, const Vector& u, const Vector& v, const Vector& w);
		void setMorth(double l, double r, double b, double t, double n, double f);
		void setMpers(double fov, double n, double f);
		void worldToEye(const Vector& vWorld, Vector& dest);
		void worldToScreen(const Vector& vWorld, Vector& vEye, Vector& vClip, Vector& dest);
		void screenToWorld(const Vector& vImage, Vector& vClip, Vector& vEye, Vector& dest);
		void screenToWorld(const Vector& vImage, Vector& dest);

		void transformEyeToClip(Vector& v, Vector& dest, bool persDivide=true);
		void transformEyeToScreen(Vector& vEye, Vector& vClip, Vector& dest);
		void transformClipToScreen(Vector& v, Vector& dest);

		bool testAndPopulateZbuffer(int x, int y, double zEye, bool update);
		bool testAndPopulateZbuffer(const Point& vWorld, bool update);

		// get the ray through pixel(i,j)
		Ray getRay(int x, int y);

		void setMatrixMode (int mode) { matrixMode = mode; }
		int getMatrixMode () { return matrixMode; }
		Matrix& getMprojection() { 
			return *Mproj[matrixMode]; 
		}

		bool isDepthEnabled() { return zmode == OSU_DEPTH_TEST; }
		void enableDepthTest() {
			if (zmode != OSU_DEPTH_TEST) {
				zmode = OSU_DEPTH_TEST;
				zbuffer.initBuffer();
			} else {
				std::cout << "enableDepthTest :: depth test already enabled" << std::endl;
			}
		}
		void disableDepthTest() {
			if (zmode == OSU_DEPTH_TEST) {
				zmode = -1;
				zbuffer.destroyBuffer();
			} else {
				std::cout << "disableDepthTest :: depth test already disabled" << std::endl;
			}
		}
		void initMvp() {
			osuGetFramebufferSize(&nx, &ny);
			smd::setMvp(nx, ny, Mvp);
			// in case matrix was singular, set to identity
			if (Mvp.inverse(iMvp) != 0) iMvp.setIdentity();
		}
		void initializeMatrices(bool initMvp=true) {
			Mcam.setIdentity();
			if (initMvp) {
				this->initMvp();
			}
			Morth->setIdentity();
			Mpers->setIdentity();
			iMcam.setIdentity();
			iMproj.setIdentity();
			setMatrixMode(0);
		}
		Matrix* getZbufferMat() {
			return zbuffer.buffer;
		}
		virtual ~ViewPoint() {
			delete Morth;
			delete Mpers;
			zbuffer.destroyBuffer();
		}
	};
	class LightSource : public ViewPoint {
	public:
		double i; // intensity for light sources
		LightSource() : ViewPoint(), i(0) {}
	};

	class GraphicsContext {
	public:
		int zmode;
		RendererType rendererType;
		RenderMode renderMode;
		int reflectionDepth;
		bool debug;
		bool debugFine;

		// debug bookkeeping
		int zeroNormals;
		int priorVertexNormals;
		int computedVertexNormals;

		// Material & color properties
		bool interpolateNormals;
		OSUShadeModel shadeModel;
		Color ambient;
		Color diffuse;
		Color specular;
		Color emittance;
		std::vector<LightSource*> directionalLights;
		std::vector<LightSource*> pointLights;

		GraphicsContext() : eye(ViewPoint()), currViewPoint(0),
			zmode(-1), rendererType(smd::SCAN_CONVERT), renderMode(smd::COLOR_AND_SHADOW), reflectionDepth(1),
			debug(false), debugFine(false), 
			zeroNormals(0), priorVertexNormals(0), computedVertexNormals(0),
			interpolateNormals(true) {
				initMatrices();
				setViewPointToEye();
		}

		ViewPoint& getViewPoint() {
			return *currViewPoint;
		}

		void setViewPoint(ViewPoint* vp) {
			currViewPoint = vp;
		}
		void setViewPoint(unsigned int vpIndex) {
			if (vpIndex == 0) {
				setViewPointToEye();
			} else if (vpIndex > 0 && vpIndex <= pointLights.size()) {
					setViewPoint(pointLights[vpIndex-1]);
			} else {
				std::cout << "Invalid viewpoint index " << vpIndex << std::endl;
				std::exit(-1);
			}
		}
		bool isCurrentViewEye() {
			return currViewPoint == &eye;
		}
		void setViewPointToEye() {
			currViewPoint = &eye;
		}
		
		bool isRenderColor() {
			return (renderMode == smd::COLOR || renderMode == COLOR_AND_SHADOW);
		}

		void transformVectorToWorld(Vector& v, Vector& vWorld, ViewPoint& view);
		void transformWorldToClip(Vector& vWorld, Vector& vEye, Vector& dest, ViewPoint& view, bool persDivide);
		Vector transformVectorToScreen(Vector& v, Vector& wCoords, Vector& eye, Vector& clipCoords, ViewPoint& view, bool persDivide=true);
		void transformVectorToClip(Vector& v, Vector& wCoords, Vector& eye, Vector& dest, ViewPoint& view, bool persDivide);
		void transformVectorToEye(Vector& v, Vector& wCoords, Vector& dest, ViewPoint& view);

		/* Depth test related methods */
		void clearZbuffer();
		Matrix* getZbuffer() {
			return eye.zbuffer.buffer;
		}
		bool isDepthEnabled() { return zmode == OSU_DEPTH_TEST; }
		void enableDepthTest() {
			zmode = OSU_DEPTH_TEST;
			eye.enableDepthTest();
			for (unsigned int i = 0; i < pointLights.size(); i++) {
				pointLights[i]->enableDepthTest();
			}
			// TODO: do the same for directional lights.
		}
		void disableDepthTest() {
			eye.disableDepthTest();
			for (unsigned int i = 0; i < pointLights.size(); i++) {
				pointLights[i]->disableDepthTest();
			}
			// TODO: do the same for directional lights.
		}

		void setInterpolateNormals(bool b) { interpolateNormals = b; }
		void pushMatrix() { matrixStack.push(); }
		void popMatrix() { matrixStack.pop(); }

		Matrix& getM() { 
			return *matrixStack.top; 
		}
		
		void loadIdentity() {
			getM().setIdentity();
		}

		void setMorth(double l, double r, double b, double t, double n, double f);

		void setMpers(double fov, double n, double f);

		void translate(double x, double y, double z);

		void scale(double x, double y, double z);

		void rotate(double angle, double x, double y, double z);

		void lookAt(double from[3], double at[3], double up[3]);

		void initialize();

		// Material and color properties

		void setShadeModel(int model) {
			switch(model) {
			case OSU_SMOOTH:
				shadeModel = OSU_SMOOTH;
				break;
			case OSU_FLAT:
				shadeModel = OSU_FLAT;
				break;
			default:
				std::cout << "Unsupported shade model " << model << std::endl;
				exit(-1);
			}
		}

		void addPointLight(float pos[3], float i) {
			LightSource *ls = new LightSource();
			ls->initializeMatrices(); // sets Mvp, iMvp;
			Vector v, l; v.set(pos[0], pos[1], pos[2], 1);
			matrix_X_vector(getM(), v, l); // point lights are in world coords
			ls->i = i;
			ls->lookAt(l, eye.at, eye.v); // sets Mcam, iMcam
			ls->nearp = eye.nearp;
			ls->farp = eye.farp;
			ls->iMproj.copyFrom(eye.iMproj);
			ls->Morth->copyFrom(*eye.Morth);
			ls->Mpers->copyFrom(*eye.Mpers);
			ls->matrixMode = eye.matrixMode;
			if (0 && debug) {
				std::cout << "Adding light from=" << l << "; at=" << eye.at << "; up=" << eye.v << std::endl;
			}
			// Note: At this point the projection matrix for light source view is orthogonal
			if (isDepthEnabled()) {
				ls->enableDepthTest();
			}
			pointLights.push_back(ls);
		}

		void addDirectionalLight(float dir[3], float i) {
			LightSource *ls = new LightSource();
			Vector v; v.set(dir[0], dir[1], dir[2], 0);
			v.normalizeCoords();
			matrix_X_vector(getM(), v, ls->e); // directional light is in world coords
			ls->i = i;
			if (isDepthEnabled()) ls->enableDepthTest();
			directionalLights.push_back(ls);
		}

		void setAmbientLight(float i) {
			ambient.set(1, 1, 1, i);
		}

		void setDiffuse(float r, float g, float b) {
			diffuse.set(r, g, b);
		}
		void setDiffuse(Color& c) {
			diffuse.set(c);
		}

		void setSpecular(float r, float g, float b, float s) {
			specular.set(r, g, b, s);
		}
		void setSpecular(Color& c) {
			specular.set(c);
		}

		void setEmittance(float r, float g, float b) {
			emittance.set(r, g, b, 0);
		}
		void setEmittance(Color& c) {
			emittance.set(c);
		}

		void setDebug(bool debug, bool debugFine) {
			this->debug = debug;
			this->debugFine = debugFine;
		}
		~GraphicsContext();
	private:
		ViewPoint eye;
		ViewPoint* currViewPoint;
		void initMatrices();
		MatrixStack matrixStack;
		Matrix *tM1; // temporary work matrix
		Matrix *tM2; // temporary work matrix
	};

	class Shape {
	public:
		int id;
		Vector normal;
		double reflectance;
		Color diffuse;
		Color specular;
		Color emittance; // Required for radiosity.
		OSUDrawable getMode() { return mode; }
		virtual void setReflectance(double ks){ reflectance = ks; }
		virtual double getReflectance(){ return reflectance; }
		virtual bool isHit(Ray& ray, HitData& hitdata, double t0, double t1) = 0;
		virtual void addVertex(double x, double y, const Color& color) {
			addVertex(x,y,0,color,VEC_ZERO);
		}
		virtual void addVertex(double x, double y, double z, const Color& color) {
			addVertex(x,y,z,color,VEC_ZERO);
		}
		virtual void addVertex(double x, double y, double z, const Color& color, const Vector& normal, const int id=0) = 0;
		virtual bool canAddVertex() = 0;
		virtual bool readyToRender() = 0;
		virtual void render(GraphicsContext& ctx) = 0;
		virtual void initNormal() { normal.set(smd::VEC_ZERO); };
		virtual const Vector& getNormal() { return normal; };
		virtual std::ostream& print(std::ostream& out) = 0;
		virtual void printColors() = 0;
	protected:
		OSUDrawable mode;
	};

	class ShapeList {
	public:
		std::list<Shape *>& getShapes() { return shapes; }
		void render(GraphicsContext& ctx) {
			if (shapes.size() == 0) {
				printf("No shapes in list to render\n");
			}
			for (std::list<Shape *>::iterator it=shapes.begin(); it != shapes.end(); ++it) {
				(*it)->render(ctx);
			}
		}
		void add(Shape *shape) {
			shapes.push_back(shape);
		}
		Shape *lastShape() {
			if (shapes.empty()) return NULL;
			return shapes.back();
		}
		Shape *popLastShape() {
			if (shapes.empty()) return NULL;
			Shape *shape = shapes.back();
			shapes.pop_back();
			return shape;
		}
		void clear() {
			while (!shapes.empty()) {
				Shape *s = shapes.front();
				shapes.pop_front();
				delete s;
			}
		}
		~ShapeList() {
			clear();
		}
	private:
		std::list<Shape *> shapes;
	};

	class Renderer {
	public:
		virtual void render(ShapeList& shapeList, GraphicsContext& ctx) {
			render(shapeList.getShapes(), ctx);
		}
		virtual void render(std::list<Shape*>& shapes, GraphicsContext& ctx) = 0;
	};

	class RayTraceRenderer : public Renderer {
	public:
		virtual void render(std::list<Shape*>& shapes, GraphicsContext& ctx);
	protected:
		virtual void colorXY(int x, int y, Color& color, std::vector<Shape*>& shapes, GraphicsContext& ctx);
		void colorXY(Ray& ray, Color& color, std::vector<Shape*>& shapes, 
			int maxdepth, double ks, double t0, double t1, GraphicsContext& ctx);
	};

	class RadiosityRenderer : public Renderer {
	public:
		virtual void render(std::list<Shape*>& shapes, GraphicsContext& ctx);
	protected:
	};

	class GraphicsManager : Renderer {
	public:
		GraphicsManager() : 
			reflectance(-1.0),
			shapeCount(0),
				drawMode(OSU_NONE), 
				ctx(GraphicsContext()) {

			currColor.set(0,0,0);
			normal.set(VEC_ZERO);

			ctx.debug = false; // TODO: set false once debugged
			ctx.debugFine = false; // TODO: set false once debugged
			
			normMul = 1.0;

			initialize();

		}
		~GraphicsManager() {
		}

		Shape *constructShape(OSUDrawable mode);

		void setRenderMode (RenderMode renderMode) { ctx.renderMode = renderMode; }

		void setRendererType (RendererType _rendererType) { ctx.rendererType = _rendererType; }
		void setReflectionDepth(int depth) { ctx.reflectionDepth = depth; }
		void setReflectance(double ks) { reflectance = ks; }

		const Color& getColor() {return currColor;}

		void setColor(double r, double g, double b) {
			currColor.set(r,g,b);
		}

		void setNormal(double x, double y, double z) {
			normal.set(normMul*x,normMul*y,normMul*z);
		}
	
		void setNormMul(double val) {
			normMul = val;
		}
	
		void clearZbuffer() {
			if (ctx.debug) std::cout << "clearing zbuffer..." << std::endl;
			ctx.clearZbuffer();
		}
		void clearDepthTest() {
			if (ctx.debug) std::cout << "clearing depth test..." << std::endl;
			ctx.disableDepthTest();
		}
		void setDepthTest(OSUEnable zmode) {
			if (ctx.debug) std::cout << "set depth test to " << zmode << std::endl;
			if (zmode == OSU_DEPTH_TEST) {
				ctx.enableDepthTest();
			} else {
				ctx.disableDepthTest();
			}
		}
		void enable(int depthTestBit) {
			if (depthTestBit == OSU_DEPTH_TEST) 
				setDepthTest(OSU_DEPTH_TEST);
			else
				clearDepthTest();
		}
		void setViewPoint(unsigned int vpIndex) {
			ctx.setViewPoint(vpIndex);
		}
		void setViewPointToEye() {
			ctx.setViewPointToEye();
		}
		int getNumPointLights() {
			return ctx.pointLights.size();
		}

		GraphicsContext& getGraphicsContext() {
			return ctx;
		}

		void setDebug(bool debug, bool debugFine) {
			ctx.debug = debug;
			ctx.debugFine = debugFine;
		}

		void addVertex(double x, double y) {
			addVertex(x,y,0);
		}
		void addVertex(double x, double y, double z) {
			Shape *shape = lastShape();
			Point tv, v; tv.setLoc(x, y, z);
			if (ctx.rendererType == smd::RAY_TRACE || ctx.rendererType == smd::RADIOSITY) {
				// convert to world coordinates if
				// we are in ray trace mode
				matrix_X_vector(ctx.getM(), tv, v);
			} else {
				v.setLoc(x, y, z);
			}
			if (shape != NULL && shape->canAddVertex()) {
				shape->addVertex(v.x,v.y,v.z,ctx.diffuse,normal);
			} else {
				if (shape != NULL && shape->readyToRender()) {
					if (ctx.rendererType == smd::SCAN_CONVERT) {
						popAndRender();
						shape = constructShape(drawMode);
					} else {
						// Ray trace
						shape = constructShape(drawMode); // automatically preserves the last shape
					}
				}
				shape->addVertex(v.x,v.y,v.z,ctx.diffuse,normal);
			}
		}

		Shape *lastShape() {
			return shapeList.lastShape();
		}
		Shape *removeLastShape() {
			return shapeList.popLastShape();
		}
		void begin(OSUDrawable mode) {
			constructShape((OSUDrawable)mode);
		}
		void sphere(double x, double y, double z, double r);
		void end() {
			if (ctx.rendererType == smd::SCAN_CONVERT) {
				popAndRender();
			} else {
				// Do nothing...
				// Render will be called separately later
			}
		}
		void popAndRender() {
			Shape *shape = removeLastShape();
			if (shape != NULL) {
				//if (ctx.debug) PRINT_OBJ("Rendering Shape ",(*shape),std::cout);
				if (ctx.debug && ctx.debugFine) PRINT_OBJ("Rendering ",(*shape),std::cout);
				shape->render(ctx);
				delete shape;
			}
		}
		virtual void render(std::list<Shape*>& shapes, GraphicsContext& ctx) {
			if (ctx.rendererType == smd::RAY_TRACE) {
				rayTracer.render(shapes, ctx);
			} else if (ctx.rendererType == smd::RADIOSITY) {
				radiosityRenderer.render(shapes, ctx);
			} else {
				std::cout << "No renderer defined" << std::endl;
			}
		}
		virtual void render(std::vector<Shape*>& shapes, GraphicsContext& ctx) {
			std::list<Shape*> l;
			vectorToList<Shape*>(shapes, l);
			render(l, ctx);
		}
		void render() {
			render(shapeList.getShapes(), ctx);
		}
		std::list<Shape*> getShapes() {
			return shapeList.getShapes();
		}

		void setDrawMode (OSUDrawable mode) { drawMode = mode; }
		void pushMatrix() { ctx.pushMatrix(); }
		void popMatrix() { ctx.popMatrix(); }

		Matrix& getM() { 
			return ctx.getM(); 
		}

		void initialize() {
			ctx.initialize();
		}

		void loadIdentity() {
			ctx.loadIdentity();
		}

		void setMorth(double l, double r, double b, double t, double n, double f) {
			ctx.setMorth(l, r, b, t, n, f);
		}

		void setMpers(double fov, double n, double f) {
			ctx.setMpers(fov, n, f);
		}

		void translate(double x, double y, double z) {
			ctx.translate(x, y, z);
		}

		void scale(double x, double y, double z) {
			ctx.scale(x, y, z);
		}

		void rotate(double angle, double x, double y, double z) {
			ctx.rotate(angle, x, y, z);
		}

		void lookAt(double from[3], double at[3], double up[3]) {
			ctx.lookAt(from, at, up);
		}

		// Material and color properties

		void setShadeModel(int model) {
			ctx.setShadeModel(model);
		}

		void addPointLight(float pos[3], float i) {
			ctx.addPointLight(pos, i);
		}

		void addDirectionalLight(float dir[3], float i) {
			ctx.addDirectionalLight(dir, i);
		}

		void setAmbientLight(float i) {
			ctx.setAmbientLight(i);
		}

		void setDiffuse(float r, float g, float b) {
			ctx.setDiffuse(r, g, b);
		}

		void setSpecular(float r, float g, float b, float s) {
			ctx.setSpecular(r, g, b, s);
		}

		void setEmittance(float r, float g, float b) {
			ctx.setEmittance(r, g, b);
		}

	private:
		double reflectance;
		int shapeCount;
		OSUDrawable drawMode;
		ShapeList shapeList;
		Color currColor;
		Vector normal;
		GraphicsContext ctx;
		double normMul;
		RayTraceRenderer rayTracer;
		RadiosityRenderer radiosityRenderer;
	}; // GraphicsManager

	void transformVerticesFromModelToWorld(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx);
	bool transformVerticesFromWorldToClip(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx);
	bool transformVerticesToClip(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx);
	void transformVerticesFromClipToScreen(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx);
	void copyToOsuVertex(Vertex& v, osuVertex& dest);
	void copyFromOsuVertex(Vertex& v, osuVertex& src);

} //namespace smd

#endif
