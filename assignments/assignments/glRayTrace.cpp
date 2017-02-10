#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "osuGraphics.h"
#include "glDraw.h"

namespace smd {
	
	/**
	 * col - returned populated with color values
	 * c_l - light intensity
	 * c_r - surface diffuse
	 */
	void addRayColor(Color& col, const double c_l, Color& c_r, const Color& specular, 
		const Vector& n, const Vector& e, const Vector& l, double ks, const bool debug);

	Shape* getHit(Ray& ray, std::vector<Shape*>& shapes, HitData& hd, double t0, double t1, GraphicsContext& ctx);
	void getTriangleColor(Triangle* shape, double alpha, double beta, double gamma, Color& dest);
	void getDiffuseColor(Shape* shape, HitData& hd, Color& color);

	void initPolygonNormals(std::vector<Shape*>& shapes) {
		for (unsigned int i = 0; i < shapes.size(); i++) {
			shapes[i]->initNormal();
		}
	}

	void getNormalAt(Shape* shape, Point& p, Vector& n) {
		switch (shape->getMode()) {
		case OSU_TRIANGLE:
		case OSU_POLYGON:
			n.set(shape->getNormal());
			break;
		case OSU_SPHERE:
			n.set(p).sub(((Sphere*)shape)->c).normalizeCoords();
			break;
		default:
			n.set(smd::VEC_ZERO);
		}
	}

	void RayTraceRenderer::render(std::list<Shape*>& shapeList, GraphicsContext& ctx) {
		std::vector<Shape*> shapes;
		listToVector<Shape*>(shapeList, shapes);
		//std::cout << "# Shapes " << shapes.size() << "; reflectionDepth=" << ctx.reflectionDepth << std::endl;
		//std::cout << "nx=" << ctx.getViewPoint().nx << ", ny=" << ctx.getViewPoint().ny << std::endl;
		initPolygonNormals(shapes);
		int w = ctx.getViewPoint().nx;
		int h = ctx.getViewPoint().ny;
		for (int x = 0; x < w; x++) {
			for (int y = 0; y < h; y++) {
				Color col; col.set(0,0,0);
				colorXY(x, y, col, shapes, ctx);
				osuWritePixel(x,y,round(col.r*255),round(col.g*255),round(col.b*255));
			}
		}
	}
	void RayTraceRenderer::colorXY(int x, int y, Color& color, std::vector<Shape*>& shapes, GraphicsContext& ctx) {
		Ray ray = ctx.getViewPoint().getRay(x,y);
		colorXY(ray, color, shapes, ctx.reflectionDepth, 1, 0, DBL_MAX, ctx);
	}
	void RayTraceRenderer::colorXY(Ray& ray, Color& color, std::vector<Shape*>& shapes, 
		int maxdepth, double ks, double t0, double t1, GraphicsContext& ctx) {
		if (maxdepth <= 0) return;
		HitData hd;
		Shape* hitObj = getHit(ray, shapes, hd, t0, t1, ctx);
		if (hitObj) {
			double c_a = ctx.ambient.s;
			Color c_r; // = ctx.diffuse;
			getDiffuseColor(hitObj, hd, c_r);
			//PRINT_OBJ("Found hit with ",(*hitObj),std::cout);
			//std::cout << "Diffuse color: " << c_r << std::endl;
			Vector p; ray.getPoint(hd.t[0], p);
			color.r = color.r + ks*c_r.r*c_a;
			color.g = color.g + ks*c_r.g*c_a;
			color.b = color.b + ks*c_r.b*c_a;
			for (unsigned int i = 0; i < ctx.pointLights.size(); i++) {
				Vector li = ctx.pointLights[i]->e; li.sub(p).normalizeCoords();
				Ray toLight(li, p);
				if (!getHit(toLight, shapes, hd, 1e-4, DBL_MAX, ctx)) {
					Vector n; getNormalAt(hitObj, p, n);
					Vector e; e.set(ray.o); e.sub(p);
					addRayColor(color, ctx.pointLights[i]->i, c_r, ctx.specular, n, e, li, ks, false);
					if (maxdepth > 1 && hitObj->getReflectance() > 0) {
						// reflection
						Vector tmp; tmp.set(n).mul(2*ray.d.dot(n));
						Vector r; r.set(ray.d).sub(tmp); // d - 2(d.n)n
						Ray reflRay(r, p);
						colorXY(reflRay, color, shapes, maxdepth-1, hitObj->getReflectance()*ks, 1e-4, DBL_MAX, ctx);
					}
				} else {
					// only ambience, no reflection
				}
			}
		}
	}

	void getTriangleColor(Triangle* shape, double alpha, double beta, double gamma, Color& dest) {
		if (shape->vertices.size() != 3) {
			std::cout << "Error in getTriangleColor(): #Triangle vertices=" << shape->vertices.size() << std::endl;
		}
		std::list<Vertex>::iterator iter = ((Triangle*)shape)->vertices.begin();
		Vertex& v1 = *iter; iter++; // a
		Vertex& v2 = *iter; iter++; // b
		Vertex& v3 = *iter;         // c
		combineColors(v1.color, v2.color, v3.color, dest, alpha, beta, gamma);
	}

	void getDiffuseColor(Shape* shape, HitData& hd, Color& color) {
		switch (shape->getMode()) {
		case OSU_TRIANGLE:
			getTriangleColor((Triangle*)shape, hd.alpha, hd.beta, hd.gamma, color);
			break;
		case OSU_SPHERE:
			color.set(((Sphere*)shape)->color);
			break;
		default:
			std::cout << "Shape type " << shape->getMode() << " not supported in getDiffuseColor" << std::endl;
		}
	}

	/**
	 * col - returned populated with color values
	 * c_l - light intensity
	 * c_r - surface diffuse
	 */
	void addRayColor(Color& color, const double c_l, Color& c_r, const Color& specular, 
		const Vector& n, const Vector& e, const Vector& l, double ks, const bool debug) {
		if (debug) {
			std::cout << "e: " << e << "; l: " << l << "; n: " << n << std::endl;
		}
		if (n.dot(e) < 0) return; // vertex not visible to eye
		double n_l = n.dot(l);
		double max_0_nl = std::max(0.0,n_l);
		Vector h; h.set(e).add(l).normalizeCoords();
		double phong = std::pow(std::max(0.0, n.dot(h)),specular.s);
		if (debug) {
			std::cout << "n.l=" << n_l << "; phong=" << phong << std::endl;
		}
		color.r = std::min(1.0, color.r + ks*c_l*(c_r.r*max_0_nl + specular.r*phong));
		color.g = std::min(1.0, color.g + ks*c_l*(c_r.g*max_0_nl + specular.g*phong));
		color.b = std::min(1.0, color.b + ks*c_l*(c_r.b*max_0_nl + specular.b*phong));
	}
	
	Shape* getHit(Ray& ray, std::vector<Shape*>& shapes, HitData& hd, double t0, double t1, GraphicsContext& ctx) {
		Shape* obj = 0;
		double _t1 = t1;
		for (unsigned int i = 0; i < shapes.size(); i++) {
			if (shapes[i]->isHit(ray, hd, t0, _t1)) {
				_t1 = hd.t[0];
				obj = shapes[i];
			}
		}
		return obj;
	}

} // namespace smd
