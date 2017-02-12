#ifndef OSUGLDRAW
#define OSUGLDRAW

#include <iostream>
#include <stdio.h>
#include <list>
#include "osuGraphics.h"
#include "glGraphicsManager.h"

namespace smd {

	class Line : public Shape {
	public:
		union {
			Vertex vertex[2];
			struct {
				Vertex v1, v2;
			};
		};
		Line() {
			vpos = 0; mode = OSU_LINES;
		}
		virtual bool isHit(Ray& ray, HitData& hitdata, double t0, double t1) {
			std::cout << "Line isHit() not supported" << std::endl;
			return false;
		}
		virtual void addVertex(double x, double y, double z, const Color& color, const Vector& normal, const int id=0) {
			if (vpos >= 2) {
				printf("Line: exceeded number of vertices: %d. Ignoring (%f,%f,%f)\n", vpos, x, y, z);
				return;
			}
			vertex[vpos].id = id;
			vertex[vpos].v.setLoc(x,y,z);
			vertex[vpos].wCoords.setLoc(x,y,z);
			vertex[vpos].color.set(color.r, color.g, color.b);
			if (VEC_ZERO.equals(normal)) 
				vertex[vpos].setNormal(normal);
			else
				vertex[vpos].setPriorNormal(normal);
			vpos++;
		}
		virtual bool canAddVertex() {
			return (vpos < 2);
		}
		virtual bool readyToRender() {
			return (vpos == 2);
		}
		virtual void render(GraphicsContext& ctx);
		std::ostream& print(std::ostream& out) {
			out << "<" << id <<  ">: Line <" << v1 << "," << v2 << ">";
			return out;
		}
		virtual void printVertexColors(std::ostream& out) {
			out << "colors <";
			for(int i = 0; i < vpos; i++) {
				out << vertex[i].color;
			}
			out << " >";
		}
		virtual void printColors() {
			printVertexColors(std::cout);
		}
		friend std::ostream& operator<< (std::ostream& stream, Line& line);
	private:
		int vpos;
	};

	class Polygon : public Shape {
	public:
		std::list<Vertex> vertices;
		Polygon() {
			mode = OSU_POLYGON;
		}
		virtual bool isHit(Ray& ray, HitData& hitdata, double t0, double t1) {
			std::cout << "Polygon isHit() not supported" << std::endl;
			return false;
		}
		virtual void addVertex(double x, double y, double z, const Color& color, const Vector& normal, const int id=0) {
			Vertex v;
			v.id = id;
			v.v.setLoc(x,y,z);
			v.wCoords.setLoc(x,y,z);
			v.color.set(color.r, color.g, color.b);
			if (VEC_ZERO.equals(normal)) 
				v.setNormal(normal);
			else
				v.setPriorNormal(normal);
			vertices.push_back(v);
		}
		virtual bool canAddVertex() {
			return true;
		}
		virtual bool readyToRender() {
			return (vertices.size() >= 3);
		}
		virtual void initNormal();
		virtual void initNormal(std::list<Vertex>& vertices);
		virtual void render(GraphicsContext& ctx);
		virtual std::ostream& print(std::ostream& out) {
			out << "<" << id <<  ">: Polygon (" << vertices.size() << ") ";
			printVertices(vertices, out);
			return out;
		}
		virtual void printVertices(std::list<Vertex>& vertices, std::ostream& out) {
			out << "verts <";
			for(std::list<Vertex>::iterator it = vertices.begin(); it != vertices.end(); it++) {
				out << *it;
			}
			out << " >";
		}
		virtual void printVertexColors(std::list<Vertex>& vertices, std::ostream& out) {
			out << "colors <";
			for(std::list<Vertex>::iterator it = vertices.begin(); it != vertices.end(); it++) {
				out << (*it).color;
			}
			out << " >";
		}
		virtual void printColors() {
			printVertexColors(vertices, std::cout);
		}
		friend std::ostream& operator<< (std::ostream& stream, Polygon& polygon);
	protected:
		virtual bool hitTriangle(Point& p1, Point& p2, Point& p3, Ray& ray, HitData& hitdata, double t0, double t1);
		void getClippedPolygonVertices(std::list<Vertex>& vertices, std::list<Vertex>& clippedVertices);
		Vector computeFaceNormal(std::list<Vertex>& vertices);
		void setVertexNormals(std::list<Vertex>& vertices, Vector& n, GraphicsContext& ctx);
		void setVertexColors(std::list<Vertex>& vertices, GraphicsContext& ctx);
	};

	class Triangle : public Polygon {
	public:
		Triangle() {
			mode = OSU_TRIANGLE;
		}
		virtual bool isHit(Ray& ray, HitData& hitdata, double t0, double t1);
		bool canAddVertex() {
			return (vertices.size() < 3);
		}
		bool readyToRender() {
			return (vertices.size() == 3);
		}
		virtual std::ostream& print(std::ostream& out) {
			out << "<" << id <<  ">: Triangle (" << vertices.size() << ") ";
			printVertices(vertices, out);
			return out;
		}
		friend std::ostream& operator<< (std::ostream& stream, Triangle& triangle);
	};

	class Sphere : public Shape {
	public:
		Sphere() {
			mode = OSU_SPHERE; // something for now...
		}
		Sphere(double x, double y, double z, double r) {
			mode = OSU_SPHERE; // something for now...
			c.setLoc(x, y, z);
			this->r = r;
		}
		Point c; // center
		double r; // radius
		Color color;
		virtual bool isHit(Ray& ray, HitData& hitdata, double t0, double t1);
		virtual void addVertex(double x, double y, const Color& color) {
			std::cout << "Sphere::addVertex is unsupported" << std::endl;
		}
		virtual void addVertex(double x, double y, double z, const Color& color) {
			std::cout << "Sphere::addVertex is unsupported" << std::endl;
		}
		virtual void addVertex(double x, double y, double z, const Color& color, const Vector& normal, const int id) {
			std::cout << "Sphere::addVertex is unsupported" << std::endl;
		}
		virtual bool canAddVertex() {
			return false;
		}
		virtual bool readyToRender() {
			return true;
		}
		virtual void render(GraphicsContext& ctx) {
			std::cout << "Sphere::render is unsupported" << std::endl;
		}
		virtual std::ostream& print(std::ostream& out) {
			out << "<" << id <<  ">: Sphere center=" << c << ", r=" << r;
			return out;
		}
		virtual void printColors() {
			std::cout << "colors <" << color << ">";
		}
		friend std::ostream& operator<< (std::ostream& stream, Sphere& triangle);
	};

	struct ImplicitLinear {
		union {
			double array[3];
			struct {
				double a, b, c;
			};
		};
		ImplicitLinear() : a(0), b(0), c(0) {}
		ImplicitLinear(const Point& p1, const Point& p2) {
			set(p1, p2);
		}
		void set(const Point& p1, const Point& p2) {
			a = p1.y - p2.y;
			b = p2.x - p1.x;
			c = p1.x*p2.y - p2.x*p1.y;
		}
		void set(const double x1, const double y1, const Point& p2) {
			a = y1 - p2.y;
			b = p2.x - x1;
			c = x1*p2.y - p2.x*y1;
		}
		void set(const double x1, const double y1, const double x2, const double y2) {
			a = y1 - y2;
			b = x2 - x1;
			c = x1*y2 - x2*y1;
		}
		double operator()(double x, double y) {
			return a*x + b*y + c;
		}
	};

	double max3(double v1, double v2, double v3);
	double min3(double v1, double v2, double v3);

	void combineColors(const Color& c1, const Color& c2, const Color& c3, Color& dest,
					   const double alpha, const double beta, const double gamma, const double scale=1);

} //namespace smd

#endif
