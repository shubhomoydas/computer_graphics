#ifndef OSUGLRADIOSITY
#define OSUGLRADIOSITY

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include "osuGraphics.h"
#include "ObjLoader.h"
#include "matlib.h"
#include "utils.h"

#define VERTEXID(i,j,lenX) ((i)*(lenX))+(j)

#define POINT_LEN(p1,p2) std::sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) + (p1.z-p2.z)*(p1.z-p2.z))
#define VERTEXPTR_LEN(v1,v2) std::sqrt((v1->v.x-v2->v.x)*(v1->v.x-v2->v.x) + (v1->v.y-v2->v.y)*(v1->v.y-v2->v.y) + (v1->v.z-v2->v.z)*(v1->v.z-v2->v.z))

#define VEC_AVG4(v,v1,v2,v3,v4) { \
	v.x=0.25*(v1.x+v2.x+v3.x+v4.x); \
	v.y=0.25*(v1.y+v2.y+v3.y+v4.y); \
	v.z=0.25*(v1.z+v2.z+v3.z+v4.z); \
	v.w=0.25*(v1.w+v2.w+v3.w+v4.w); \
}
#define VEC_DIFF(v,v1,v2) {v.x=v1.x-v2.x; v.y=v1.y-v2.y; v.z=v1.z-v2.z; v.w=v1.w-v2.w;}
// ignore w in VEC_DOT
#define VEC_DOT(v1,v2) (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z)
// ignore w in VEC_LEN
#define VEC_LEN(v) (std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z))

#define COLOR_INTERPOLATE(c,c1,c2,alpha) { \
	c.r = (1-alpha)*c1.r + alpha*c2.r; \
	c.g = (1-alpha)*c1.g + alpha*c2.g; \
	c.b = (1-alpha)*c1.b + alpha*c2.b; \
	c.s = (1-alpha)*c1.s + alpha*c2.s; \
}

#define VEC_INTERPOLATE(v,v1,v2,alpha) { \
	v.x = (1-alpha)*v1.x + alpha*v2.x; \
	v.y = (1-alpha)*v1.y + alpha*v2.y; \
	v.z = (1-alpha)*v1.z + alpha*v2.z; \
	v.w = (1-alpha)*v1.w + alpha*v2.w; \
}

namespace smd {

	inline double area(const Point& p1, const Point& p2, const Point& p3);

	class Patch : public Polygon {
	public:
		int parentid;
		Color intensity;
		Vertex *varr[4];
		double area;
		// assumes that the vertex array has been initialized for convenient access
		virtual bool isHit(Ray& ray, HitData& hitdata, double t0, double t1) {
			if (vertices.size() != 4) {
				std::cout << "Number of vertices in Patch is not 4" << std::endl;
				std::exit(-1);
			}
			bool hit = false;
			hit = hitTriangle(varr[0]->v, varr[1]->v, varr[2]->v, ray, hitdata, t0, t1);
			if (!hit) {
				hit = hitTriangle(varr[2]->v, varr[3]->v, varr[0]->v, ray, hitdata, t0, t1);
			}
			return hit;
		}
		double computeArea() {
			this->area = smd::area(varr[0]->v, varr[1]->v, varr[2]->v);
			return this->area;
		}
		void cacheVertices() {
			if (vertices.size() != 4) {
				std::cout << "Number of vertices in Patch is not 4" << std::endl;
				std::exit(-1);
			}
			std::list<Vertex>::iterator viter = vertices.begin();
			varr[0] = &(*viter); viter++;
			varr[1] = &(*viter); viter++;
			varr[2] = &(*viter); viter++;
			varr[3] = &(*viter);
		}
	};

	void computePatchAreas(std::vector<Patch*>& patches);

	class PatchDelta {
	public:
		int i; // index into the patch array
		Patch* p;
		double delta;
		int colorband;
		double emittance;
		double ro;
		double Fj;
		PatchDelta(int i, Patch* p, int colorband) : i(i), p(p), delta(0), colorband(colorband), emittance(0), ro(0), Fj(0) {
			emittance = p->emittance.array[colorband];
			ro = p->diffuse.array[colorband];
			delta = emittance;
		}
	};

	class RadiositySolver {
	public:
		Matrix *x;
		RadiositySolver() : x(0) {}
		virtual void solve(int iters) = 0;
		virtual double getAmbient(int i) { return 0; }
		virtual double getAmbient() { return 0; }
	};

	class ProgressiveRadiosity : public RadiositySolver {
	public:
		int colorband;
		Matrix& FF;
		std::vector<Patch*>& patches;
		std::vector<PatchDelta> pds;
		double areasum;
		double R;
		double roavg;
		double ambient;
		ProgressiveRadiosity(std::vector<Patch*>& patches, Matrix &FF, int colorband) : 
			colorband(colorband), patches(patches), FF(FF), pds(std::vector<PatchDelta>()),
			areasum(0), R(0), roavg(0), ambient(0) {
			computePatchAreas(patches);
			setupPatchDeltas(patches, colorband);
			setupMatrices(patches, FF, colorband, pds);
		}
		virtual void solve(int iters) {
			double *d_B = x->data[0];
			for (int i = 0; i < iters; i++) {
				// process patch with max unshot energy
				int bestid = getPatchDeltaWithHighestDeltaArea();
				PatchDelta& pdi = pds[bestid]; // patches are in sorted order of unshot energy
				for (unsigned int j = 0; j < pds.size(); j++) {
					PatchDelta& pdj = pds[j];
					double dRad = pdj.ro * pdi.delta * FF.data[pdi.i][pdj.i] * pdi.p->area / pdj.p->area;
					pdj.delta = pdj.delta + dRad; // update delta B_j
					d_B[pdj.i] = d_B[pdj.i] + dRad; // update B_j
				}
				pdi.delta = 0;
			}
			ambient = computeAmbient(pds, R);
		}
		virtual double getAmbient() {
			return ambient;
		}
		virtual double getAmbient(int i) {
			return pds[i].ro*ambient;
		}
		virtual ~ProgressiveRadiosity() {
			if (x) delete x;
		}
	private:
		int getPatchDeltaWithHighestDeltaArea() {
			int idx = 0;
			double best = -1.0;
			for (unsigned int i = 0; i < pds.size(); i++) {
				PatchDelta& pdi = pds[i];
				if (pdi.delta * pdi.p->area > best) {
					idx = i; best = pdi.delta * pdi.p->area;
				}
			}
			return idx;
		}
		void setupPatchDeltas(std::vector<Patch*>& patches, int colorband) {
			for (unsigned int i = 0; i < patches.size(); i++) {
				pds.push_back(PatchDelta(i, patches[i], colorband));
			}
		}
		double computeAmbient(std::vector<PatchDelta>& pds, double R) {
			double ambient = 0;
			for (unsigned int i = 0; i < pds.size(); i++) {
				ambient = ambient + pds[i].delta*pds[i].Fj;
			}
			return R*ambient;
		}
		void setupMatrices(std::vector<Patch*>& patches, Matrix &FF, int colorband, std::vector<PatchDelta>& pds) {
			if (patches.size() != FF.rows) {
				std::cout << "Formfactor dimensions do not match number of patches" << std::endl;
				std::exit(-1);
			}
			x = new Matrix(1, FF.rows); x->setZeros();
			R = 0;
			roavg = 0;
			areasum = 0;
			double roareasum = 0;
			for (unsigned int i = 0; i < patches.size(); i++) {
				x->data[0][i] = patches[i]->emittance.array[colorband];
				areasum = areasum + patches[i]->area;
				roareasum = roareasum + pds[i].ro*patches[i]->area;
			}
			for (unsigned int i = 0; i < patches.size(); i++) {
				pds[i].Fj = pds[i].p->area / areasum;
			}
			roavg = roareasum / areasum;
			R = 1 / (1-roavg);
			ambient = computeAmbient(pds, R);
			std::cout << "Ambient=" << ambient << "; R=" << R << "; roavg=" << roavg << "; areasum=" << areasum << std::endl;
		}
	};

	// Solves Ax = b
	// x0 is starting point
	class RadiosityLinearSystem : public RadiositySolver {
	public:
		Matrix *A;
		Matrix *b;
		Matrix *x0;
		int colorband;
		RadiosityLinearSystem(std::vector<Patch*>& patches, Matrix &FF, int colorband) : 
			A(0), b(0), x0(0), colorband(colorband) {
			setupMatrices(patches, FF);
		}
		virtual void solve(int iters) {
			gauss_seidel(*A, *b, *x0, *x, iters);
			x0->copyFrom(*x);
		}
		virtual ~RadiosityLinearSystem() {
			if (A) delete A;
			if (b) delete b;
			if (x) delete x;
			if (x0) delete x0;
		}
	private:
		void setupMatrices(std::vector<Patch*>& patches, Matrix &FF) {
			if (patches.size() != FF.rows) {
				std::cout << "Formfactor dimensions do not match number of patches" << std::endl;
				std::exit(-1);
			}
			A = new Matrix(FF.rows);
			b = new Matrix(1, FF.rows);
			x = new Matrix(1, FF.rows); x->setZeros();
			x0 = new Matrix(1, FF.rows); x0->setZeros();
			double **df = FF.data;
			double **da = A->data;
			double *db = b->data[0];
			for (int i = 0; i < FF.rows; i++) {
				Patch *p = patches[i];
				double ro = p->diffuse.array[colorband];
				double *ptr_df = df[i];
				double *ptr_da = da[i];
				for (int j = 0; j < FF.cols; j++, ptr_df++, ptr_da++) {
					*ptr_da = -ro*(*ptr_df);
				}
				da[i][i] = 1-da[i][i];
				db[i] = p->emittance.array[colorband];
			}
		}
	};

	class RadiosityLinearSystemRGB {
	public:
		RadiositySolver *rgb[3];
		bool addPRAmbient;
		// solverType: 0 - Gauss-Seidel; 1 - Progressive
		RadiosityLinearSystemRGB(std::vector<Patch*>& patches, Matrix& FF, const int solverType=0, const bool addPRAmbient=false) {
			this->addPRAmbient = addPRAmbient;
			if (solverType == 1) {
				rgb[0] = new ProgressiveRadiosity(patches, FF, 0);
				rgb[1] = new ProgressiveRadiosity(patches, FF, 1);
				rgb[2] = new ProgressiveRadiosity(patches, FF, 2);
			} else {
				rgb[0] = new RadiosityLinearSystem(patches, FF, 0);
				rgb[1] = new RadiosityLinearSystem(patches, FF, 1);
				rgb[2] = new RadiosityLinearSystem(patches, FF, 2);
			}
		}
		void solve(int iters) {
			rgb[0]->solve(iters);
			rgb[1]->solve(iters);
			rgb[2]->solve(iters);
		}
		void loadSolutions(std::vector<Patch*>& patches) {
			if (patches.size() != rgb[0]->x->cols) {
				std::cout << "Solution dimensions do not match number of patches" << std::endl;
				std::exit(-1);
			}
			unsigned int n = patches.size();
			for (unsigned int i = 0; i < n; i++) {
				Patch* p = patches[i];
				p->diffuse.array[0] = rgb[0]->x->data[0][i] + (addPRAmbient ? rgb[0]->getAmbient(i) : 0);
				p->diffuse.array[1] = rgb[1]->x->data[0][i] + (addPRAmbient ? rgb[1]->getAmbient(i) : 0);
				p->diffuse.array[2] = rgb[2]->x->data[0][i] + (addPRAmbient ? rgb[2]->getAmbient(i) : 0);
			}
		}
		void printAmibient() {
			std::cout << "Ambient R=" << rgb[0]->getAmbient() 
				<< "; Ambient G=" << rgb[1]->getAmbient() 
				<< "; Ambient B=" << rgb[2]->getAmbient() << std::endl;
		}
		~RadiosityLinearSystemRGB() {
			delete rgb[0];
			delete rgb[1];
			delete rgb[2];
		}
	};

	void testMakePatch();
	void cachePatchVerticesAndInitNormals(std::vector<Patch*>& patches); // faster access to vertices
	std::vector<Shape*>& typecastPatchToShapeVector(std::vector<Patch*>& patches, std::vector<Shape*>& shapes);
	void color_mapping(double percentage, double col[3]);
	inline Color& affineInterp(const Color & p1, const Color& p2, const double alpha, Color& dest);
	inline Vector& affineInterp(const Vector & p1, const Vector& p2, const double alpha, Vector& dest);
	std::vector<smd::Patch*>& subDivideByMaxLength(std::list<Shape*>& shapes, double len, std::vector<Patch*>& patches);
	std::vector<smd::Patch*>& subDivideByMaxLength(Polygon& rect, double len, std::vector<Patch*>& patches);
	std::vector<smd::Patch*>& subDivideByMaxLength(const Vertex& p1, const Vertex& p2, const Vertex& p3, const Vertex& p4, 
		double len, std::vector<Patch*>& patches, const Polygon& rect);
	std::vector<smd::Patch*>& subDivideByNumParts(std::list<Shape*>& shapes, int nparts, std::vector<Patch*>& patches);
	std::vector<smd::Patch*>& subDivideByNumParts(Polygon& rect, int nparts, std::vector<Patch*>& patches);
	std::vector<smd::Patch*>& subDivideByNumParts(
		const Vertex& p1, const Vertex& p2, const Vertex& p3, const Vertex& p4, 
		const int cntX, const int cntY, std::vector<Patch*>& patches, const Polygon& rect);

	std::map<int, std::vector<Patch*>*>* groupPatches(std::vector<Patch*>& patches);
	void printPatches(std::map<int, std::vector<Patch*>*> *pmap);
	std::map<int, std::map<int, std::vector<Patch*>* >* >* groupVertices(std::map<int, std::vector<Patch*>*>* pmap);
	void printVertexGroups(std::map<int, std::map<int, std::vector<Patch*>* >* >* vmap);
	void averageVertexColorsForGroup(std::map<int, std::map<int, std::vector<Patch*>* >* >* vmap);

	void setRandomColor(Color& color);
	void setFaceColorsToRandom(std::list<Shape*>& shapes);
	void setFaceColorsToRandom(std::vector<Patch*>& patches);

	// Copies the face color to all vertices.
	void setVertexColorsToFace(Polygon& poly);
	// Copies the face color to all vertices.
	void setVertexColorsToFace(std::vector<Patch*>& patches);

	void computeSimpleFormFactors(std::vector<Patch*>& patches, Matrix& FF, int nparts, int **HID=0);

	int** create2DIntArray(int rows, int cols);
	void destroy2DIntArray(int **arr, int rows);
	bool isPatchVisible(int i1, int i2, std::vector<Patch*>& patches);
	int **getHID(std::vector<Patch*>& patches);

	void set_random_seed(unsigned int seed);
	/**
	 * Generates a uniform random number in range [0.0,1.0]
	 */
	double unif_rand_01();

	/**
	 * Generates a long random integer. The rand() function usually only
	 * returns a value in the range 0-35567. We sometimes need larger values.
	 */
	int randint();

} // namespace smd

#endif // OSUGLRADIOSITY
