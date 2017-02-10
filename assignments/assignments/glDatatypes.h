#ifndef OSUGLDATATYPES
#define OSUGLDATATYPES

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <vector>
#include "osuGraphics.h"
#include "glMatrix.h"

namespace smd {

	std::string trim(const std::string& str, const std::string& chars = " \t\n");
	void tokenize(const std::string& str, std::vector<std::string>& tokens,
						const std::string& delimiters);

	template<typename T>
	void release_vector_elements(std::vector< T > *v) {
		if (!v) return;
		for (std::vector< T >::iterator it = v->begin(); it != v->end(); it++) {
			delete *it;
		}
	}

	template<typename _K, typename _V>
	void release_map_elements(std::map<_K,_V>* m) {
		if (!m) return;
		typedef std::map<_K, _V> InputType;
		typename InputType::iterator it;
		for (it = m->begin(); it != m->end(); it++) {
			_V val = it->second;
			if (val) delete val;
		}
	}

	template<typename _K, typename _V>
	void release_map_elements_with_map(std::map<_K,_V*>* m) {
		if (!m) return;
		typedef std::map<_K,_V*> InputType;
		typedef _V ContainerType;
		typename InputType::iterator it;
		for (it = m->begin(); it != m->end(); it++) {
			_V* collection = it->second;
			typename ContainerType::iterator vit;
			if (collection) {
				for (vit = collection->begin(); vit != collection->end(); vit++) {
					delete vit->second;
				}
				delete collection;
			}
		}
	}

	template<typename _K, typename _V>
	void release_map_elements_with_container(std::map<_K,_V*>* m) {
		if (!m) return;
		typedef std::map<_K,_V*> InputType;
		typedef _V ContainerType;
		typename InputType::iterator it;
		for (it = m->begin(); it != m->end(); it++) {
			_V* collection = it->second;
			typename ContainerType::iterator vit;
			if (collection) {
				for (vit = collection->begin(); vit != collection->end(); vit++) {
					delete *vit;
				}
				delete collection;
			}
		}
	}

	template<class T>
	void release_objects(T** arr, int len) {
		if (!arr) return;
		for (int i = 0; i < len; i++) {
			delete arr[i];
		}
	}

	class Vector {
	public:
		union {
			double array[4];
			struct {
				double x;
				double y;
				double z;
				double w;
			};
		};
		Vector& set(const Vector& v) {
			x = v.x; y = v.y; z = v.z; w = v.w;
			return *this;
		}
		Vector& set(double _x, double _y, double _z) {
			x = _x; y = _y; z = _z; w = 0;
			return *this;
		}
		Vector& set(double _x, double _y, double _z, double _w) {
			x = _x; y = _y; z = _z; w = _w;
			return *this;
		}
		Vector& setLoc(double _x, double _y, double _z) {
			x = _x; y = _y; z = _z; w = 1;
			return *this;
		}
		double norm() {
			return std::sqrt(x*x + y*y + z*z + w*w);
		}
		Vector& normalize() {
			double s = norm();
			if (s > 0) {
				x = x/s; y = y/s; z = z/s; w = w/s;
			}
			return *this;
		}
		// The *Coords() methods will ignore w
		double normCoords() const {
			return std::sqrt(x*x + y*y + z*z);
		}
		Vector& normalizeCoords() {
			double s = normCoords();
			if (s > 0) {
				x = x/s; y = y/s; z = z/s;
			}
			return *this;
		}
		double dist(const Vector& p) const {
			double x2 = x-p.x, y2 = y-p.y, z2 = z-p.z;
			return std::sqrt(x2*x2 + y2*y2 + z2*z2);
		}
		// this = this + b
		Vector& add(const Vector& b) {
			x = x+b.x;
			y = y+b.y;
			z = z+b.z;
			w = w+b.w;
			return *this;
		}
		// this = this - b
		Vector& sub(const Vector& b) {
			x = x-b.x;
			y = y-b.y;
			z = z-b.z;
			w = w-b.w;
			return *this;
		}
		double dot(const Vector& b) const {
			return x*b.x + y*b.y + z*b.z; // ignoring w here...
		}
		Vector cross(const Vector& b) const {
			Vector cv;
			cv.set(
				y*b.z - z*b.y,
				z*b.x - x*b.z,
				x*b.y - y*b.x
			);
			return cv;
		}
		Vector& mul(double val) {
			x = x*val;
			y = y*val;
			z = z*val;
			w = w*val;
			return *this;
		}
		bool equals(const Vector& v) const {
			return (x==v.x && y==v.y && z==v.z && w==v.w);
		}
		bool isNan() const {
			return x != x || y != y || z != z || w != w;
		}
		bool isInf() const {
			return 
				(x > DBL_MAX || x < -DBL_MAX)
				|| (y > DBL_MAX || y < -DBL_MAX)
				|| (z > DBL_MAX || z < -DBL_MAX)
				|| (w > DBL_MAX || w < -DBL_MAX);
		}
		bool parse(const std::string str);
		static Vector zeroVector() {
			Vector v;
			v.set(0, 0, 0, 0);
			return v;
		}
		std::ostream& print(std::ostream& out) const {
			out << "<" << x << "," << y << "," << z << ";" << w << ">";
			return out;
		}
		friend std::ostream& operator<< (std::ostream& stream, const Vector& vector);
	};

	static const Vector VEC_ZERO = Vector::zeroVector();

	class Color {
	public:
		union {
			double array[4];
			struct {
				double r;
				double g;
				double b;
				double s; // sharpness coefficient for specular; intensity otherwise
			};
		};
		void set(double _r, double _g, double _b) {
			r = _r; g = _g; b = _b; s = 1;
		}
		void set(double _r, double _g, double _b, double _s) {
			r = _r; g = _g; b = _b; s = _s;
		}
		void set(const Color& c) {
			r = c.r; g = c.g; b = c.b; s = c.s;
		}
		bool parse(const std::string str);
		static Color zeroColor() {
			Color c;
			c.set(0, 0, 0, 0);
			return c;
		}
		std::ostream& print(std::ostream& out) const {
			out << "<" << r << "," << g << "," << b << ";" << s << ">";
			return out;
		}
		friend std::ostream& operator<< (std::ostream& stream, const Color& color);
	};

	static const Color COLOR_ZERO = Color::zeroColor();

	typedef Vector Point;
	typedef Vector Quaternion;

	class Vertex {
	public:
		int id;
		Point v;
		Color color;
		Point wCoords; // world coordinates - needs to be preserved for lighting
		Point tEye; // sometimes we need the eye coordinates e.g. in z-buffer
		Vector normal;
		bool priorNormal;
		void set(const Vertex& src) {
			id = src.id;
			v.set(src.v);
			color.set(src.color);
			wCoords.set(src.wCoords);
			tEye.set(src.tEye);
			normal.set(src.normal);
			priorNormal = src.priorNormal;
		}
		void setNormal(const Vector& n) {
			normal.set(n);
			priorNormal = false;
		}
		void setPriorNormal(const Vector& n) {
			normal.set(n);
			priorNormal = true;
		}
		std::ostream& print(std::ostream& out) const {
			out << "[" << id << "]" << v;
			return out;
		}
		friend std::ostream& operator<< (std::ostream& stream, const Vertex& vertex);
	};

	class MatrixStack {
	public:
		Matrix *top;
		MatrixStack() {
			top = 0;
			initialize();
		}
		void initialize();
		void push();
		Matrix *pop();
		int size();
		~MatrixStack();
	private:
		void clear();
		std::list<Matrix *> mats;
	};

	void saveData(const Matrix& mat, const std::string& filepath);
	Matrix* loadCSV(const std::string& filename, const int skiptoplines=0, const int skipleftcols=0, const int skiprightcols=0);

} //namespace smd

#endif
