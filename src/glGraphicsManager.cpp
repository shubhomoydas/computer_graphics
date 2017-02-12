#include <iostream>
#include <stdio.h>
#include <list>
#include "osuGraphics.h"
#include "glGraphicsManager.h"
#include "glDraw.h"

namespace smd {

	// Assume that p1 -> p2 -> p3 is clockwise
	Vector getNormal(const Point& p1, const Point& p2, const Point& p3) {
		Vector u, v;
		u.set(p2); v.set(p3);
		u.sub(p1); v.sub(p2);
		return v.cross(u).normalizeCoords();
	}

	/**
	 *  ViewPoint
	 */
	Ray ViewPoint::getRay(int i, int j) {
		double _u = l + (r-l)*(i+0.5)/nx;
		double _v = b + (t-b)*(j+0.5)/ny;
		Vector d, o;
		if (matrixMode == 0) {
			// orthographic
			o.x = e.x + _u*u.x + _v*v.x;
			o.y = e.y + _u*u.y + _v*v.y;
			o.z = e.z + _u*u.z + _v*v.z;
			o.w = 1;
			d.set(w).mul(-1);
		} else {
			// perspective
			o.set(e);
			d.x = -nearp*w.x + _u*u.x + _v*v.x;
			d.y = -nearp*w.y + _u*u.y + _v*v.y;
			d.z = -nearp*w.z + _u*u.z + _v*v.z;
			d.w = 0;
		}
		return Ray(d, o);
	}

	void ViewPoint::lookAt(double from[3], double at[3], double up[3]) {
		Vector _e; _e.setLoc(from[0],from[1],from[2]);
		Vector _at; _at.setLoc(at[0], at[1], at[2]);
		Vector _up; _up.set(up[0], up[1], up[2], 0);
		lookAt(_e, _at, _up);
	}

	void ViewPoint::lookAt(const Vector& e, const Vector& at, const Vector& up) {
		this->e.set(e);
		this->at.set(at);
		this->v.set(up).normalize();
		this->w.set(this->at).sub(this->e).normalize().mul(-1);
		this->u.set(this->v.cross(this->w).normalize());
		this->v.set(this->w.cross(this->u));

		smd::lookAt(this->e, this->u, this->v, this->w, Mcam);
		Mcam.inverse(iMcam);
		if (0) {
			std::cout << "View from=" << this->e << "; w=" << this->w << 
				"; at=" << this->at << "; up=" << this->v << std::endl
				<< "Mcam:" << std::endl << Mcam << std::endl;
		}
	}

	void ViewPoint::setMorth(double l, double r, double b, double t, double n, double f) {
		smd::setMorth(l, r, b, t, n, f, *Morth);
		setMatrixMode(0);
		setNearFar(n,f);
		this->l = l; this->r = r; this->t = t; this->b = b;
		Morth->inverse(iMproj);
	}

	void ViewPoint::setMpers(double fov, double n, double f) {
		// t = |n|*tan(theta/2) - refer to Shirley pg. 157
		n = std::max(1.0,n);
		double t = abs(n)*std::tan(0.5*fov*PI/180);
		double b = -t;
		double r =  t;
		double l = -r;
		smd::setMpers(l, r, b, t, n, f, *Mpers);
		setMatrixMode(1);
		setNearFar(n,f);
		this->l = l; this->r = r; this->t = t; this->b = b;
		Mpers->inverse(iMproj);
	}

	void ViewPoint::worldToEye(const Vector& vWorld, Vector& dest) {
		matrix_X_vector(Mcam, vWorld, dest);
	}
	void ViewPoint::transformEyeToClip(Vector& vEye, Vector& dest, bool persDivide) {
		matrix_X_vector(getMprojection(), vEye, dest); // afterProj = Mproj*v
		if (persDivide) dest.mul(1/dest.w);
	}
	void ViewPoint::worldToScreen(const Vector& vWorld, Vector& vEye, Vector& vClip, Vector& dest) {
		worldToEye(vWorld, vEye);
		transformEyeToScreen(vEye, vClip, dest);
	}
	void ViewPoint::screenToWorld(const Vector& vImage, Vector& vClip, Vector& vEye, Vector& dest) {
		matrix_X_vector(iMvp, vImage, vClip);
		matrix_X_vector(iMproj, vClip, vEye);
		if (vEye.w != 0.0) vEye.mul(1/vEye.w);
		matrix_X_vector(iMcam, vEye, dest);
	}
	void ViewPoint::screenToWorld(const Vector& vImage, Vector& dest) {
		Vector vClip, vEye;
		screenToWorld(vImage, vClip, vEye, dest);
	}

	void ViewPoint::transformClipToScreen(Vector& v, Vector& dest) {
		matrix_X_vector(Mvp, v, dest); // dest = Mvp*v
	}
	void ViewPoint::transformEyeToScreen(Vector& vEye, Vector& vClip, Vector& dest) {
		transformEyeToClip(vEye, vClip, true);
		transformClipToScreen(vClip, dest);
	}

	bool ViewPoint::testAndPopulateZbuffer(int x, int y, double zEye, bool update) {
		if (zmode != OSU_DEPTH_TEST) return true; // true if zbuffer is disabled
		bool success = false;
		double z = zEye;
		if (matrixMode == 1) { // perspective
			if (z != 0.0) z = z;
		} else {
			z = z;
		}
		if (!update) z = z + 1e-6;
		if (z >= zbuffer.buffer->data[y][x]) {
			if (update) zbuffer.buffer->data[y][x] = z;
			success = true;
		}
		return success;
	}
	bool ViewPoint::testAndPopulateZbuffer(const Point& vWorld, bool update) {
		if (zmode != OSU_DEPTH_TEST) return true; // true if zbuffer is disabled
		Vector vEye, vClip, vImage;
		worldToScreen(vWorld, vEye, vClip, vImage);
		int ix = round(vImage.x), iy = round(vImage.y);
		if (ix >= 0 && iy >= 0 && ix < nx && iy < ny)  {
			return testAndPopulateZbuffer(ix, iy, vEye.z, update);
		}
		return false;
	}

	/**
	 *  GraphicsManager
	 */
	Shape *GraphicsManager::constructShape(OSUDrawable mode) {
		Shape *shape = NULL;
		setDrawMode(mode);
		switch (mode) {
		case OSU_LINES:
			shape = new Line();
			shapeList.add(shape);
			break;
		case OSU_TRIANGLE:
			shape = new Triangle();
			shapeList.add(shape);
			break;
		case OSU_POLYGON:
			shape = new Polygon();
			shapeList.add(shape);
			break;
		case OSU_SPHERE:
			shape = new Sphere();
			shapeList.add(shape);
			break;
		default:
			std::cout << "Unsupported shape type..." << mode << std::endl;
			exit(-1);
		}
		shape->id = shapeCount++;
		shape->setReflectance(reflectance);
		shape->specular.set(ctx.specular);
		shape->diffuse.set(ctx.diffuse);
		shape->emittance.set(ctx.emittance);
		return shape;
	}
	void GraphicsManager::sphere(double x, double y, double z, double r) {
		constructShape(OSU_SPHERE);
		Point c; c.setLoc(x, y, z);
		Sphere* sphere = (Sphere*)shapeList.lastShape();
		matrix_X_vector(ctx.getM(), c, sphere->c);
		sphere->r = r;
		sphere->color.set(ctx.diffuse);
	}

	/**
	 *  GraphicsContext
	 */
	void GraphicsContext::clearZbuffer() {
		eye.zbuffer.clearBuffer();
	}
	void GraphicsContext::transformVectorToWorld(Vector& v, Vector& vWorld, ViewPoint& view) {
		Matrix& M = getM();
		matrix_X_vector(M, v, vWorld); // wCoords = M*v
	}
	void GraphicsContext::transformWorldToClip(Vector& vWorld, Vector& vEye, Vector& dest, ViewPoint& view, bool persDivide) {
		view.worldToEye(vWorld, vEye);
		view.transformEyeToClip(vEye, dest, persDivide);
	}
	void GraphicsContext::transformVectorToEye(Vector& v, Vector& wCoords, Vector& dest, ViewPoint& view) {
		Matrix& M = getM();
		matrix_X_vector(M, v, wCoords); // wCoords = M*v
		view.worldToEye(wCoords, dest); // dest = Mcam*(M*v1)
	}
	Vector GraphicsContext::transformVectorToScreen(Vector& v, Vector& wCoords, Vector& eye, Vector& clipCoords, ViewPoint& view, bool persDivide) {
		Vector v2, v3;
		transformVectorToEye(v, wCoords, eye, view);
		view.transformEyeToClip(eye, v2, persDivide);
		clipCoords.set(v2);
		if (!persDivide) v2.mul(1/v2.w); // v2 was earlier not divided by w; do it now.
		view.transformClipToScreen(v2, v3);
		return v3;
	}
	void GraphicsContext::transformVectorToClip(Vector& v, Vector& wCoords, Vector& eye, Vector& dest, ViewPoint& view, bool persDivide) {
		transformVectorToEye(v, wCoords, eye, view);
		view.transformEyeToClip(eye, dest, persDivide);
	}

	void GraphicsContext::setMorth(double l, double r, double b, double t, double n, double f) {
		if (debug) {
			std::cout << "Orthogonal to l=" << l << ", r=" << r << 
				", b=" << b << ", t=" << t << ", n=" << n << ", f=" << f << std::endl;
		}
		eye.setMorth(l, r, b, t, n, f);
		if (debug && debugFine) {
			std::cout << "Morth=\n" << *(eye.Morth) << std::endl;
		}
	}

	void GraphicsContext::setMpers(double fov, double n, double f) {
		if (debug) {
			std::cout << "Perspecive to fov=" << fov << ", n=" << n << ", f=" << f << std::endl;
		}
		eye.setMpers(fov, n, f);
		if (debug && debugFine) {
			std::cout << "Mpers=\n" << *(eye.Mpers) << std::endl;
		}
	}

	void GraphicsContext::translate(double x, double y, double z) {
		Matrix& M = getM();
		if (debug) {
			std::cout << "Translate to " << x << "," << y << "," << z << std::endl;
		}
		smd::translate(x, y, z, M, *tM1, *tM2);
		M.copyFrom(*tM2);
	}

	void GraphicsContext::scale(double x, double y, double z) {
		Matrix& M = getM();
		if (debug) {
			std::cout << "Scale to " << x << "," << y << "," << z << std::endl;
		}
		smd::scale(x, y, z, M, *tM1, *tM2);
		M.copyFrom(*tM2);
	}

	void GraphicsContext::rotate(double angle, double x, double y, double z) {
		Quaternion q;
		Matrix& M = getM();
		smd::getRotationQuaternion(angle, x, y, z, q);
		smd::getQuaternionToRotationMatrix(q, *tM1);
		M.mul(*tM1, *tM2);
		M.copyFrom(*tM2);
	}

	void GraphicsContext::lookAt(double from[3], double at[3], double up[3]) {
		eye.lookAt(from, at, up);
		if (debug && debugFine) {
			std::cout << "Mcam=\n" << eye.Mcam << std::endl;
		}
	}

	void GraphicsContext::initialize() {
		matrixStack.initialize();
		eye.initializeMatrices();
		if (debug && debugFine) {
			std::cout << "M=\n" << getM() << std::endl;
			std::cout << "Mcam=\n" << eye.Mcam << std::endl;
			std::cout << "Mvp=\n" << eye.Mvp << std::endl;
			std::cout << "Morth=\n" << *(eye.Morth) << std::endl;
			std::cout << "Mpers=\n" << *(eye.Mpers) << std::endl;
		}
		// Material and color
		diffuse.set(0, 0, 0, 0);
		ambient.set(0, 0, 0, 0);
		specular.set(0, 0, 0, 0);
		shadeModel = OSU_FLAT;
	}
	GraphicsContext::~GraphicsContext() {
		release_vector_elements<LightSource*>(&pointLights);
		release_vector_elements<LightSource*>(&directionalLights);
		delete tM1;
		delete tM2;
	}
	void GraphicsContext::initMatrices() {
		eye.initializeMatrices();

		tM1 = new Matrix(4);
		tM2 = new Matrix(4);
	}

	// set up viewport matrix that transforms canonical
	// coordinates to screen coordinates.
	void setMvp(int w, int h, Matrix& Mvp) {
		// Refer to:
		// Shirley pg. 144
		double **a = Mvp.data;
		a[0][0] = (1.0*w)/2.0;   a[0][1] =           0;   a[0][2] = 0;   a[0][3] = (1.0*w-1.0)/2.0;
		a[1][0] =           0;   a[1][1] = (1.0*h)/2.0;   a[1][2] = 0;   a[1][3] = (1.0*h-1.0)/2.0;
		a[2][0] =           0;   a[2][1] =           0;   a[2][2] = 1;   a[2][3] =               0;
		a[3][0] =           0;   a[3][1] =           0;   a[3][2] = 0;   a[3][3] =               1;
	}

	// set up orthogonal projection matrix.
	void setMorth(double l, double r, double b, double t, double n, double f, Matrix& Morth) {
		// Refer to:
		// Shirley pg. 145
		double **a = Morth.data;
		a[0][0] = 2.0/(r-l);   a[0][1] =         0;   a[0][2] =          0;   a[0][3] = -(r+l)/(r-l);
		a[1][0] =         0;   a[1][1] = 2.0/(t-b);   a[1][2] =          0;   a[1][3] = -(t+b)/(t-b);
		a[2][0] =         0;   a[2][1] =         0;   a[2][2] = -2.0/(f-n);   a[2][3] = -(f+n)/(f-n);
		a[3][0] =         0;   a[3][1] =         0;   a[3][2] =          0;   a[3][3] =            1;
	}

	// set up perspective projection matrix.
	void setMpers(double l, double r, double b, double t, double n, double f, Matrix& Mpers) {
		//std::cout << "Pers: t=" << t << ", b=" << b << ", r=" << r << ", l=" << l << std::endl;
		// Refer to:
		// Shirley pg. 155
		double **a = Mpers.data;
		a[0][0] = 2*n/(r-l);   a[0][1] =         0;   a[0][2] =  (r+l)/(r-l);   a[0][3] =            0;
		a[1][0] =         0;   a[1][1] = 2*n/(t-b);   a[1][2] =  (t+b)/(t-b);   a[1][3] =            0;
		a[2][0] =         0;   a[2][1] =         0;   a[2][2] = -(f+n)/(f-n);   a[2][3] = -2*f*n/(f-n);
		a[3][0] =         0;   a[3][1] =         0;   a[3][2] =           -1;   a[3][3] =            0;
	}

	Matrix& translate(double x, double y, double z, Matrix& m, Matrix& temp, Matrix& dest) {
		temp.setIdentity();
		temp.data[0][3] = x;
		temp.data[1][3] = y;
		temp.data[2][3] = z;
		m.mul(temp, dest);
		return dest;
	}

	Matrix& scale(double x, double y, double z, Matrix& m, Matrix& temp, Matrix& dest) {
		temp.setIdentity();
		temp.data[0][0] = x;
		temp.data[1][1] = y;
		temp.data[2][2] = z;
		m.mul(temp, dest);
		return dest;
	}

	// Note: g is a direction vector, NOT the center of focus.
	Matrix& lookAt(const Vector& e, const Vector& u, const Vector& v, const Vector& w, Matrix& dest) {

		// Refer to:
		// Shirley pg. 147

		// Here we have computed the two-matrix product shown in Shirley.
		// The last column becomes the -ve dot products and the block
		// diagonal is a product with identity matrix and hence simple
		// to write out.
		double **a = dest.data;
		a[0][0] = u.x; a[0][1] = u.y; a[0][2] = u.z; a[0][3] = -e.dot(u);
		a[1][0] = v.x; a[1][1] = v.y; a[1][2] = v.z; a[1][3] = -e.dot(v);
		a[2][0] = w.x; a[2][1] = w.y; a[2][2] = w.z; a[2][3] = -e.dot(w);
		a[3][0] =   0; a[3][1] =   0; a[3][2] =   0; a[3][3] =         1;

		return dest;
	}

	Quaternion& getRotationQuaternion(double angle, double x, double y, double z, Quaternion& q) {
		double n = std::sqrt(x*x + y*y + z*z);
		if (n < 1e-18) {
			std::cout << "Unstable norm for quaternion: " << n << std::endl;
			exit(-1);
		}

		// Refer to:
		// http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
		double theta = angle*PI/180;
		double s = std::sin(theta/2);
		q.w = std::cos(theta/2);
		q.x = s*x/n;
		q.y = s*y/n;
		q.z = s*z/n;

		//std::cout << "Quaternion " << q << " has length " << q.norm() << std::endl; // should be 1

		return q;
	}

	Matrix& getQuaternionToRotationMatrix(Quaternion& q, Matrix& m) {
		if (abs(q.norm()-1.0) > 1e-8) {
			std::cout << "Quaternion norm = " << q.norm() << ". Only unit quaternions supported." << std::endl;
			exit(-1);
		}
		double **a = m.data;

		double xx = q.x*q.x, yy = q.y*q.y, zz = q.z*q.z;
		double xy = q.x*q.y, xz = q.x*q.z, xw = q.x*q.w;
		double yz = q.y*q.z, yw = q.y*q.w;
		double zw = q.z*q.w;
		
		// Refer to:
		// http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
		a[0][0] = 1-2*(yy+zz);  a[0][1] =   2*(xy-zw);  a[0][2] =   2*(xz+yw);  a[0][3] = 0;
		a[1][0] =   2*(xy+zw);  a[1][1] = 1-2*(xx+zz);  a[1][2] =   2*(yz-xw);  a[1][3] = 0;
		a[2][0] =   2*(xz-yw);  a[2][1] =   2*(yz+xw);  a[2][2] = 1-2*(xx+yy);  a[2][3] = 0;
		a[3][0] =           0;  a[3][1] =           0;  a[3][2] =           0;  a[3][3] = 1;

		return m;
	}

	// dest = m * v
	Vector& matrix_X_vector(const Matrix& m, const Vector& v, Vector& dest) {
		if (!(m.rows == 4 && m.cols == 4)) {
			std::cout << "Error: Matrix must be 4x4!" << std::endl;
			exit(-1);
		}
		double** a = m.data;
		matrix_X_vecArray(m, v.array, dest.array);
		return dest;
	}

	inline double* matrix_X_vecArray(const Matrix& m, const double* v1, double* v2) {
		if (!(m.rows == 4 && m.cols == 4)) {
			std::cout << "Error: Matrix must be 4x4!" << std::endl;
			exit(-1);
		}
		double** a = m.data;
		v2[0] = a[0][0]*v1[0] + a[0][1]*v1[1] + a[0][2]*v1[2] + a[0][3]*v1[3];
		v2[1] = a[1][0]*v1[0] + a[1][1]*v1[1] + a[1][2]*v1[2] + a[1][3]*v1[3];
		v2[2] = a[2][0]*v1[0] + a[2][1]*v1[1] + a[2][2]*v1[2] + a[2][3]*v1[3];
		v2[3] = a[3][0]*v1[0] + a[3][1]*v1[1] + a[3][2]*v1[2] + a[3][3]*v1[3];
		return v2;
	}

	void transformVerticesToEye(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx) {
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			Vertex v; v.set(*iter);
			ctx.transformVectorToEye(v.v, v.wCoords, v.tEye, ctx.getViewPoint());
			v.v.set(v.tEye);
			dest.push_back(v);
		}
	}

	/* Returns true if all transformed vectors are within finite range */
	void transformVerticesFromModelToWorld(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx) {
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			Vertex v; v.set(*iter);
			ctx.transformVectorToWorld((*iter).v, v.v, ctx.getViewPoint());
			v.wCoords.set(v.v); // world coordinates same as vector
			ctx.transformVectorToWorld((*iter).normal, v.normal, ctx.getViewPoint());
			v.normal.normalizeCoords();
			dest.push_back(v);
		}
	}

	/* Returns true if all transformed vectors are within finite range */
	bool transformVerticesFromWorldToClip(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx) {
		bool degenerate = false;
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			Vertex v; v.set(*iter);
			Vector tv; ctx.transformWorldToClip(v.v, v.tEye, tv, ctx.getViewPoint(), true);
			v.v.set(tv);
			if (v.v.isInf() || v.v.isNan()) {
				degenerate = true;
				//std::cout << "degenerate wCoords=" << v.wCoords << "; tEye=" << v.tEye << "; vClip=" << tv << std::endl;
			}
			dest.push_back(v);
		}
		return !degenerate;
	}

	/* Returns true if all transformed vectors are within finite range */
	bool transformVerticesToClip(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx) {
		bool degenerate = false;
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			Vertex v; v.set(*iter);
			Vector tv; ctx.transformVectorToClip(v.v, v.wCoords, v.tEye, tv, ctx.getViewPoint(), true);
			v.v.set(tv);
			if (v.v.isInf() || v.v.isNan()) {
				degenerate = true;
				//std::cout << "degenerate wCoords=" << v.wCoords << "; tEye=" << v.tEye << "; vClip=" << tv << std::endl;
			}
			dest.push_back(v);
		}
		return !degenerate;
	}

	void transformVerticesFromClipToScreen(std::list<Vertex>& vertices, std::list<Vertex>& dest, GraphicsContext& ctx) {
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			Vertex v; v.set(*iter);
			Vector tv; ctx.getViewPoint().transformClipToScreen(v.v, tv);
			v.v.set(tv);
			dest.push_back(v);
		}
	}

	void copyToOsuVertex(Vertex& v, osuVertex& dest) {
		dest.x = (float)v.v.x;
		dest.y = (float)v.v.y;
		dest.z = (float)v.v.z;
		dest.r = (float)v.color.r;
		dest.g = (float)v.color.g;
		dest.b = (float)v.color.b;
	}

	void copyFromOsuVertex(Vertex& v, osuVertex& src) {
		v.v.x = src.x;
		v.v.y = src.y;
		v.v.z = src.z;
		v.v.w = 1;
		v.color.r = src.r;
		v.color.g = src.g;
		v.color.b = src.b;
		v.color.s = 1;
	}

	double max3(double v1, double v2, double v3) {
		if (v1 > v2)
			return v1 > v3 ? v1 : v3;
		else
			return v2 > v3 ? v2 : v3;
	}

	double min3(double v1, double v2, double v3) {
		if (v1 < v2)
			return v1 < v3 ? v1 : v3;
		else
			return v2 < v3 ? v2 : v3;
	}

} // namespace smd