#include <map>
#include "glRadiosity.h"

namespace smd {
	
	int** create2DIntArray(int rows, int cols) {
		int **arr = new int*[rows];
		for (int i = 0; i < rows; i++) {
			arr[i] = new int[cols];
		}
		return arr;
	}

	void destroy2DIntArray(int **arr, int rows) {
		for (int i = 0; i < rows; i++) {
			delete[] arr[i];
		}
		delete[] arr;
	}

	void computePatchAreas(std::vector<Patch*>& patches) {
		for (unsigned int i = 0; i < patches.size(); i++) {
			patches[i]->computeArea();
		}
	}

	// Assumes that face normals have been pre computed.
	bool isPatchVisible(int i1, int i2, std::vector<Patch*>& patches) {
		bool visible = false;
		unsigned int n = patches.size();
		Patch* p1 = patches[i1];
		Patch* p2 = patches[i2];
		Point c1, c2;
		VEC_AVG4(c1,p1->varr[0]->v,p1->varr[1]->v,p1->varr[2]->v,p1->varr[3]->v);
		VEC_AVG4(c2,p2->varr[0]->v,p2->varr[1]->v,p2->varr[2]->v,p2->varr[3]->v);
		Vector d; VEC_DIFF(d,c2,c1);
		if (VEC_DOT(p1->normal,d) < 1e-4) {
			// if faces are not facing each other they cannot see
			return false;
		}
		Ray ray(d, c1);
		HitData hd;
		if (p2->isHit(ray, hd, 0, DBL_MAX)) {
			visible = true;
			double t1 = hd.t[0];
			for (unsigned int i = 0; i < n; i++) {
				if (i == i1 || i == i2) continue;
				Point c3;
				Patch* p3 = patches[i];
				VEC_AVG4(c3,p3->varr[0]->v,p3->varr[1]->v,p3->varr[2]->v,p3->varr[3]->v);
				Vector dd; VEC_DIFF(dd,c3,c1);
				if (VEC_DOT(p1->normal,dd) < 1e-4) {
					// if faces are not facing each other they cannot see
					continue;
				}
				if (patches[i]->isHit(ray, hd, 0, t1)) {
					visible = false;
					break;
				}
			}
		}
		return visible;
	}

	int **getHID(std::vector<Patch*>& patches) {
		unsigned int n = patches.size();
		int** HID = create2DIntArray(n,n);
		for (unsigned int i = 0; i < n; i++) {
			HID[i][i] = 0;
			for (unsigned int j = i+1; j < n; j++) {
				// Note: patches with same parent id are on the same plane
				// and therefore cannot see each other
				if ((patches[i]->parentid != patches[j]->parentid) && isPatchVisible(i,j,patches)) {
					HID[i][j] = 1;
				} else {
					HID[i][j] = 0;
					//std::cout << "Patch " << *(patches[i]) << " is not visible from " << *(patches[j]) << std::endl;
				}
				HID[j][i] = HID[i][j];
			}
			if ((i+1) % 50 == 0) {
				std::cout << "Processed " << (i+1) << " rows of HID" << std::endl;
			}
		}
		return HID;
	}

	inline double area(const Point& p1, const Point& p2, const Point& p3) {
		return POINT_LEN(p1,p2)*POINT_LEN(p2,p3);
	}

	inline double area(const Patch* p) {
		Vertex *v1 = p->varr[0], *v2 = p->varr[1], *v3 = p->varr[2];
		double l1 = VERTEXPTR_LEN(v1,v2);
		double l2 = VERTEXPTR_LEN(v2,v3);
		return l1*l2;
	}

	inline double computeSimpleFormFactors(
		const Point& p11, const Point& p12, const Point& p13, const Point&p14, const Vector& n1,
		const Point& p21, const Point& p22, const Point& p23, const Point&p24, const Vector& n2) {

		if (VEC_DOT(n1,n2) > 1e-7) {
			//std::cout << "n1=" << n1 << "; n2=" << n2 << std::endl;
			//std::cout << "Faces do not face each other" << std::endl;
			return 0; // faces do not face each other
		}

		Point c1, c2;
		VEC_AVG4(c1,p11,p12,p13,p14);
		VEC_AVG4(c2,p21,p22,p23,p24);

		Vector d;
		VEC_DIFF(d,c1,c2);

		double r = VEC_LEN(d); //d.normCoords();

		d.normalizeCoords();
		double cos_phi1 = VEC_DOT(d,n1); //d.dot(n1);
		double cos_phi2 = VEC_DOT(d,n2); //d.dot(n2);
		double da1 = area(p11, p12, p13), da2 = area(p21, p22, p23);

		//std::cout << "p11=" << p11 << "; p12=" << p12 << "; p13=" << p13 << "; p14=" << p14 << "; c1=" << c1 << "; c2=" << c2 << std::endl;
		//std::cout << "cos1=" << cos_phi1 << "; cos2=" << cos_phi2 << "; da1=" << da1 << "; da2=" << da2 << "; r=" << r << std::endl;

		return std::fabs(cos_phi1) * std::fabs(cos_phi2) * da1 * da2 / (PI * r *r);
	}

	// Assumes patch normals have been pre-computed
	inline double computeSimpleFormFactors(
		const Point& p11, const Point& p12, const Point& p13, const Point&p14, const Vector& n1,
		Patch* patch2, const int nparts) {
		int cntX = nparts, cntY = nparts;
		double alphaY = 0, dalphaY=1.0/(double)cntY, dalphaX=1.0/(double)cntX;
		//patch2->initNormal();
		Vertex *v1 = patch2->varr[0], *v2 = patch2->varr[1], *v3 = patch2->varr[2], *v4 = patch2->varr[3]; 
		Point pp1, pp2;
		Point pp1i, pp2i;
		pp1.set(v1->v); pp2.set(v2->v);
		double ff = 0;
		int patchid = 0;
		for (int i = 0; i < cntY; i++) {
			alphaY = alphaY+dalphaY;
			VEC_INTERPOLATE(pp1i,v1->v,v4->v,alphaY);
			VEC_INTERPOLATE(pp2i,v2->v,v3->v,alphaY);
			double alphaX = 0;
			Point pp3, pp4, pp3i, pp4i;
			pp3.set(pp1); pp4.set(pp1i);
			for (int j = 0; j < cntX; j++) {
				alphaX = alphaX+dalphaX;
				VEC_INTERPOLATE(pp3i,pp1,pp2,alphaX);
				VEC_INTERPOLATE(pp4i,pp1i,pp2i,alphaX);

				// create patch pp3,pp3i,pp4i,pp4
				ff = ff + computeSimpleFormFactors(p11, p12, p13, p14, n1, pp3, pp3i, pp4i, pp4, patch2->normal);

				pp3.set(pp3i);
				pp4.set(pp4i);
			}
			pp1.set(pp1i);
			pp2.set(pp2i);
		}
		return ff;
	}

	// Compute form factor for emission incoming from patch2 to patch1.
	// Assumes patch normals have been pre-computed
	double computeSimpleFormFactors(Patch* patch1, Patch* patch2, int nparts) {
		int cntX = nparts, cntY = nparts;
		double alphaY = 0, dalphaY=1.0/(double)cntY, dalphaX=1.0/(double)cntX;
		//patch1->initNormal();
		Vertex *v1 = patch1->varr[0], *v2 = patch1->varr[1], *v3 = patch1->varr[2], *v4 = patch1->varr[3]; 
		Point pp1, pp2;
		Point pp1i, pp2i;
		pp1.set(v1->v); pp2.set(v2->v);
		int patchid = 0;
		double ff = 0;
		for (int i = 0; i < cntY; i++) {
			alphaY = alphaY+dalphaY;
			VEC_INTERPOLATE(pp1i,v1->v,v4->v,alphaY);
			VEC_INTERPOLATE(pp2i,v2->v,v3->v,alphaY);
			double alphaX = 0;
			Point pp3, pp4, pp3i, pp4i;
			pp3.set(pp1); pp4.set(pp1i);
			for (int j = 0; j < cntX; j++) {
				alphaX = alphaX+dalphaX;
				VEC_INTERPOLATE(pp3i,pp1,pp2,alphaX);
				VEC_INTERPOLATE(pp4i,pp1i,pp2i,alphaX);

				// compute form factor for patch pp3,pp3i,pp4i,pp4
				ff = ff + computeSimpleFormFactors(pp3, pp3i, pp4i, pp4, patch1->normal, patch2, nparts);

				pp3.set(pp3i);
				pp4.set(pp4i);
			}
			pp1.set(pp1i);
			pp2.set(pp2i);
		}
		return ff;
	}

	void computeSimpleFormFactors(std::vector<Patch*>& patches, Matrix& FF, int nparts, int **HID) {
		if (FF.rows != patches.size() || FF.rows != FF.cols) {
			std::cout << "Matrix dimensions do not match with number of patches (" << patches.size() << ")" << std::endl;
			std::exit(-1);
		}
		double **a = FF.data;
		for (unsigned int i = 0; i < patches.size(); i++) {
			double ai = area(patches[i]);
			a[i][i] = 0;
			for (unsigned int j = i+1; j < patches.size(); j++) {
				if (!HID || (HID && HID[i][j] == 1)) {
					double aj = area(patches[j]);
					double aff = computeSimpleFormFactors(patches[i], patches[j], nparts);
					a[i][j] = (1/ai)*aff;
					a[j][i] = (1/aj)*aff;
				} else {
					a[i][j] = 0;
					a[j][i] = 0;
				}
			}
			if ((i+1) % 50 == 0) {
				std::cout << "Processed " << (i+1) << " rows of FF" << std::endl;
			}
		}
	}

	void averageVertexColorsForGroup(std::vector<Patch*>* patches, int vertexid) {
		Color c; c.set(0,0,0,0);
		for(unsigned int i = 0; i < patches->size(); i++) {
			for (std::list<Vertex>::iterator vit = (*patches)[i]->vertices.begin(); vit != (*patches)[i]->vertices.end(); vit++) {
				if (vit->id != vertexid) continue;
				c.r = c.r+vit->color.r;
				c.g = c.g+vit->color.g;
				c.b = c.b+vit->color.b;
			}
		}
		c.r = c.r/patches->size();
		c.g = c.g/patches->size();
		c.b = c.b/patches->size();
		for(unsigned int i = 0; i < patches->size(); i++) {
			for (std::list<Vertex>::iterator vit = (*patches)[i]->vertices.begin(); vit != (*patches)[i]->vertices.end(); vit++) {
				if (vit->id != vertexid) continue;
				vit->color.set(c);
			}
		}
	}
	void averageVertexColorsForGroup(std::map<int, std::map<int, std::vector<Patch*>* >* >* vmap) {
		for(std::map<int, std::map<int, std::vector<Patch*>* >* >::iterator it = vmap->begin(); it != vmap->end(); it++) {
			int parentid = it->first;
			std::map<int, std::vector<Patch*>* >* vs = it->second;
			for(std::map<int, std::vector<Patch*>*>::iterator viter = vs->begin(); viter != vs->end(); viter++) {
				int vertexid = viter->first;
				std::vector<Patch*>* patches = viter->second;
				averageVertexColorsForGroup(patches, vertexid);
			}
		}
	}

	// Assumes that the vertices have been stored in vertex array in the
	// patch structure for optimized access. Else call groupVerticesUnoptimized()
	std::map<int, std::map<int, std::vector<Patch*>* >* >* groupVertices(std::map<int, std::vector<Patch*>*>* pmap) {
		std::map<int, std::map<int, std::vector<Patch*>* >* > *vmap = new std::map<int, std::map<int, std::vector<Patch*>* >* >();
		for(std::map<int, std::vector<Patch*>*>::iterator iter = pmap->begin(); iter != pmap->end(); iter++) {
			std::map<int, std::vector<Patch*>*>* p = 0;
			std::map<int, std::map<int, std::vector<Patch*>* >* >::iterator it1 = vmap->find(iter->first);
			if (it1 == vmap->end()) {
				p = new std::map<int, std::vector<Patch*>*>();
				vmap->insert(std::pair<int, std::map<int, std::vector<Patch*>* >* >(iter->first, p));
			} else {
				p = it1->second;
			}
			std::vector<Patch*>* patches = iter->second;
			for (unsigned int i = 0; i < patches->size(); i++) {
				for (unsigned int j = 0; j < 4; j++) {
					std::vector<Patch*>* c = 0;
					Vertex* v = (*patches)[i]->varr[j];
					std::map<int, std::vector<Patch*>*>::iterator it = p->find(v->id);
					if (it == p->end()) {
						c = new std::vector<Patch*>();
						p->insert(std::pair<int, std::vector<Patch*>*>(v->id, c));
					} else {
						c = it->second;
					}
					c->push_back((*patches)[i]);
				}
			}
		}
		//release_map_elements_with_map<int, std::vector<Patch*>>(vmap);
		return vmap;
	}

	std::map<int, std::map<int, std::vector<Patch*>* >* >* groupVerticesUnoptimized(std::map<int, std::vector<Patch*>*>* pmap) {
		std::map<int, std::map<int, std::vector<Patch*>* >* > *vmap = new std::map<int, std::map<int, std::vector<Patch*>* >* >();
		for(std::map<int, std::vector<Patch*>*>::iterator iter = pmap->begin(); iter != pmap->end(); iter++) {
			std::map<int, std::vector<Patch*>*>* p = 0;
			std::map<int, std::map<int, std::vector<Patch*>* >* >::iterator it1 = vmap->find(iter->first);
			if (it1 == vmap->end()) {
				p = new std::map<int, std::vector<Patch*>*>();
				vmap->insert(std::pair<int, std::map<int, std::vector<Patch*>* >* >(iter->first, p));
			} else {
				p = it1->second;
			}
			std::vector<Patch*>* patches = iter->second;
			for (unsigned int i = 0; i < patches->size(); i++) {
				for (std::list<Vertex>::iterator vit = (*patches)[i]->vertices.begin(); vit != (*patches)[i]->vertices.end(); vit++) {
					std::vector<Patch*>* c = 0;
					std::map<int, std::vector<Patch*>*>::iterator it = p->find(vit->id);
					if (it == p->end()) {
						c = new std::vector<Patch*>();
						p->insert(std::pair<int, std::vector<Patch*>*>(vit->id, c));
					} else {
						c = it->second;
					}
					c->push_back((*patches)[i]);
				}
			}
		}
		//release_map_elements_with_map<int, std::vector<Patch*>>(vmap);
		return vmap;
	}

	void printVertexGroups(std::map<int, std::map<int, std::vector<Patch*>* >* >* vmap) {
		for(std::map<int, std::map<int, std::vector<Patch*>* >* >::iterator it = vmap->begin(); it != vmap->end(); it++) {
			int parentid = it->first;
			std::map<int, std::vector<Patch*>* >* vs = it->second;
			for(std::map<int, std::vector<Patch*>*>::iterator viter = vs->begin(); viter != vs->end(); viter++) {
				int vertexid = viter->first;
				std::vector<Patch*>* patches = viter->second;
				for(unsigned int i = 0; i < patches->size(); i++) {
					std::cout << "Parentid=" << parentid << "; Vertex=" << vertexid << "; Patchid=" << (*patches)[i]->id << std::endl;
				}
			}
		}
	}

	std::map<int, std::vector<Patch*>*>* groupPatches(std::vector<Patch*>& patches) {
		std::map<int, std::vector<Patch*>*> *pmap = new std::map<int, std::vector<Patch*>*>();
		for (unsigned int i = 0; i < patches.size(); i++) {
			std::vector<Patch*>* c = 0;
			std::map<int, std::vector<Patch*>*>::iterator it = pmap->find(patches[i]->parentid);
			if (it == pmap->end()) {
				c = new std::vector<Patch*>();
				pmap->insert(std::pair<int, std::vector<Patch*>*>(patches[i]->parentid, c));
			} else {
				c = it->second;
			}
			c->push_back(patches[i]);
		}
		//release_map_elements<int, std::vector<Patch*>*>(pmap);
		return pmap;
	}

	void printPatches(std::map<int, std::vector<Patch*>*> *pmap) {
		for(std::map<int, std::vector<Patch*>*>::iterator iter = pmap->begin(); iter != pmap->end(); iter++) {
			std::vector<Patch*>* c = iter->second;
			std::cout << "Parent Patch=" << iter->first << "; #patches=" << c->size() << std::endl;
		}
	}

	void cachePatchVerticesAndInitNormals(std::vector<Patch*>& patches) {
		for(unsigned int i = 0; i < patches.size(); i++) {
			patches[i]->cacheVertices();
			patches[i]->initNormal();
		}
	}
	
	std::vector<Shape*>& typecastPatchToShapeVector(std::vector<Patch*>& patches, std::vector<Shape*>& shapes) {
		for(unsigned int i = 0; i < patches.size(); i++) {
			shapes.push_back(patches[i]);
		}
		return shapes;
	}
	
	// Generate random color
	void setRandomColor(Color& color) {
		double r = unif_rand_01();
		color_mapping(r, color.array);
	}

	void setFaceColorsToRandom(std::vector<Patch*>& patches) {
		for(unsigned int i = 0; i < patches.size(); i++) {
			setRandomColor(patches[i]->diffuse);
		}
	}

	void setFaceColorsToRandom(std::list<Shape*>& shapes) {
		for (std::list<Shape*>::iterator it = shapes.begin(); it != shapes.end(); it++) {
			setRandomColor((*it)->diffuse);
		}
	}

	// Copies the face color to all vertices.
	void setVertexColorsToFace(Polygon& poly) {
		for (std::list<Vertex>::iterator iter = poly.vertices.begin(); iter != poly.vertices.end(); iter++) {
			(*iter).color.set(poly.diffuse);
		}
	}

	// Copies the face color to all vertices.
	void setVertexColorsToFace(std::vector<Patch*>& patches) {
		for(unsigned int i = 0; i < patches.size(); i++) {
			Polygon *p = (Polygon*)patches[i];
			setVertexColorsToFace(*p);
		}
	}

	void RadiosityRenderer::render(std::list<Shape*>& shapes, GraphicsContext& ctx) {
		for (std::list<Shape*>::iterator it = shapes.begin(); it != shapes.end(); it++) {
			Polygon* s = (Polygon*)*it;
			ctx.setDiffuse(s->diffuse);
			ctx.setSpecular(s->specular);
			//setVertexColorsToFace(*s);
			s->render(ctx);
		}
	}

	// Unit test
	void testMakePatch() {
		smd::Shape* s = new smd::Polygon();
		s->addVertex(50,50,0,smd::COLOR_ZERO, smd::VEC_ZERO);
		s->addVertex(50,-50,0,smd::COLOR_ZERO, smd::VEC_ZERO);
		s->addVertex(-50,-50,0,smd::COLOR_ZERO, smd::VEC_ZERO);
		s->addVertex(-50,50,0,smd::COLOR_ZERO, smd::VEC_ZERO);

		smd::Polygon *poly = (smd::Polygon *)s;
		std::vector<smd::Patch*> patches;
		smd::subDivideByMaxLength(*poly, 3, patches);

		std::cout << "Make Patch test ended with " << patches.size() << " patches." << std::endl;

		smd::release_vector_elements<smd::Patch*>(&patches);
		delete s;
	}

	std::vector<smd::Patch*>& subDivideByMaxLength(std::list<Shape*>& shapes, double len, std::vector<Patch*>& patches) {
		for (std::list<Shape*>::iterator it = shapes.begin(); it != shapes.end(); it++) {
			if ((*it)->getMode() != OSU_POLYGON) {
				std::cout << "Only rectangular patches supported..." << std::endl;
				continue;
			}
			Polygon* s = (Polygon*)(*it);
			if (s->vertices.size() != 4) {
				std::cout << "Only rectangular patches supported..." << std::endl;
				continue;
			}
			subDivideByMaxLength(*s, len, patches);
		}
		return patches;
	}
	
	std::vector<smd::Patch*>& subDivideByMaxLength(Polygon& rect, double len, std::vector<Patch*>& patches) {
		std::list<Vertex>::iterator viter = rect.vertices.begin();
		Vertex v1 = *viter; viter++;
		Vertex v2 = *viter; viter++;
		Vertex v3 = *viter; viter++;
		Vertex v4 = *viter; viter++;
		return subDivideByMaxLength(v1, v2, v3, v4, len, patches, rect);
	}

	std::vector<smd::Patch*>& subDivideByMaxLength(
		const Vertex& p1, const Vertex& p2, const Vertex& p3, const Vertex& p4, 
		const double len, std::vector<Patch*>& patches, const Polygon& rect) {
		double minLenX = std::min(p1.v.dist(p2.v),p3.v.dist(p4.v));
		double minLenY = std::min(p1.v.dist(p4.v),p2.v.dist(p3.v));
		if (minLenX < len || minLenY < len) {
			std::cout << "Poly side is smaller than patch size. minLenX=" << minLenX << "; minLenY=" << minLenY << "; len=" << len << std::endl;
		}
		int cntX = (int)std::max(1.0,std::floor((minLenX + len/2) / len));
		int cntY = (int)std::max(1.0,std::floor((minLenY + len/2) / len));
		return subDivideByNumParts(p1, p2, p3, p4, cntX, cntY, patches, rect);
	}

	std::vector<smd::Patch*>& subDivideByNumParts(std::list<Shape*>& shapes, int nparts, std::vector<Patch*>& patches) {
		for (std::list<Shape*>::iterator it = shapes.begin(); it != shapes.end(); it++) {
			if ((*it)->getMode() != OSU_POLYGON) {
				std::cout << "Only rectangular patches supported..." << std::endl;
				continue;
			}
			Polygon* s = (Polygon*)(*it);
			if (s->vertices.size() != 4) {
				std::cout << "Only rectangular patches supported..." << std::endl;
				continue;
			}
			subDivideByNumParts(*s, nparts, patches);
		}
		return patches;
	}
	
	std::vector<smd::Patch*>& subDivideByNumParts(Polygon& rect, int nparts, std::vector<Patch*>& patches) {
		std::list<Vertex>::iterator viter = rect.vertices.begin();
		Vertex v1 = *viter; viter++;
		Vertex v2 = *viter; viter++;
		Vertex v3 = *viter; viter++;
		Vertex v4 = *viter; viter++;
		return subDivideByNumParts(v1, v2, v3, v4, nparts, nparts, patches, rect);
	}

	std::vector<smd::Patch*>& subDivideByNumParts(
		const Vertex& p1, const Vertex& p2, const Vertex& p3, const Vertex& p4, 
		const int cntX, const int cntY, std::vector<Patch*>& patches, const Polygon& rect) {
		double alphaY = 0, dalphaY=1.0/(double)cntY, dalphaX=1.0/(double)cntX;
		Point pp1, pp2;
		Point pp1i, pp2i;
		Color cc1, cc2;
		Color cc1i, cc2i;
		pp1.set(p1.v); pp2.set(p2.v);
		cc1.set(p1.color); cc2.set(p2.color);
		int patchid = 0;
		for (int i = 0; i < cntY; i++) {
			alphaY = alphaY+dalphaY;
			VEC_INTERPOLATE(pp1i,p1.v,p4.v,alphaY);
			VEC_INTERPOLATE(pp2i,p2.v,p3.v,alphaY);
			COLOR_INTERPOLATE(cc1i,p1.color,p4.color,alphaY);
			COLOR_INTERPOLATE(cc2i,p2.color,p3.color,alphaY);
			double alphaX = 0;
			Point pp3, pp4, pp3i, pp4i;
			Color cc3, cc4, cc3i, cc4i;
			pp3.set(pp1); pp4.set(pp1i);
			cc3.set(cc1); cc4.set(cc1i);
			for (int j = 0; j < cntX; j++) {
				alphaX = alphaX+dalphaX;

				VEC_INTERPOLATE(pp3i,pp1,pp2,alphaX);
				VEC_INTERPOLATE(pp4i,pp1i,pp2i,alphaX);
				COLOR_INTERPOLATE(cc3i,cc1,cc2,alphaX);
				COLOR_INTERPOLATE(cc4i,cc1i,cc2i,alphaX);

				// create patch pp3,pp3i,pp4i,pp4
				Patch *p = new Patch();
				p->parentid = rect.id;
				p->id = patchid++;
				p->diffuse.set(rect.diffuse);
				p->specular.set(rect.specular);
				p->emittance.set(rect.emittance);
				p->reflectance = rect.reflectance;
				p->intensity.set(COLOR_ZERO);
				// Note: there are cntX+1 vertices on each row
				p->addVertex(pp3.x, pp3.y, pp3.z, cc3, smd::VEC_ZERO, VERTEXID(i,j,cntX+1));
				p->addVertex(pp3i.x, pp3i.y, pp3i.z, cc3i, smd::VEC_ZERO, VERTEXID(i,j+1,cntX+1));
				p->addVertex(pp4i.x, pp4i.y, pp4i.z, cc4i, smd::VEC_ZERO, VERTEXID(i+1,j+1,cntX+1));
				p->addVertex(pp4.x, pp4.y, pp4.z, cc4, smd::VEC_ZERO, VERTEXID(i+1,j,cntX+1));
				patches.push_back(p);
				//std::cout << "Patch: " << (*p) << std::endl;

				pp3.set(pp3i);
				pp4.set(pp4i);
				cc3.set(cc3i);
				cc4.set(cc4i);
			}
			pp1.set(pp1i);
			pp2.set(pp2i);
			cc1.set(cc1i);
			cc2.set(cc2i);
		}
		return patches;
	}

	inline Color& affineInterp(const Color & p1, const Color& p2, const double alpha, Color& dest) {
		dest.r = (1-alpha)*p1.r + alpha*p2.r;
		dest.g = (1-alpha)*p1.g + alpha*p2.g;
		dest.b = (1-alpha)*p1.b + alpha*p2.b;
		dest.s = (1-alpha)*p1.s + alpha*p2.s;
		return dest;
	}

	inline Vector& affineInterp(const Vector & p1, const Vector& p2, const double alpha, Vector& dest) {
		dest.x = (1-alpha)*p1.x + alpha*p2.x;
		dest.y = (1-alpha)*p1.y + alpha*p2.y;
		dest.z = (1-alpha)*p1.z + alpha*p2.z;
		dest.w = (1-alpha)*p1.w + alpha*p2.w;
		return dest;
	}

	void set_random_seed(unsigned int seed) {
		srand(seed);
	}

	double unif_rand_01()
	{
		return rand() / (double)(RAND_MAX + 1.0);
	}

	int randint() {
		return 65535*(rand()%65535) + rand()%65535;
	}

	void color_mapping(double percentage, double col[3])
	{
		if (percentage == 0.0){
			col[0] = 1.0;
			col[1] = 1.0;
			col[2] = 1.0;
		}
		else if (percentage <= 1.0/3){
			col[0] = 1.0-percentage*3.0;
			col[1] = 1.0;
			col[2] = 1.0-percentage*3.0;
		}
		else if (percentage <= 2.0/3){
			col[0] = percentage*3.0-1.0;
			col[1] = 1.0;
			col[2] = 0.0;
		}
		else if (percentage <= 3.0/3){
			col[0] = 1.0;
			col[1] = 3.0-percentage*3.0;
			col[2] = 0.0;
		}
		else {
			col[0] = 1.0;
			col[1] = 1.0;
			col[2] = 0.0;
		}
	}

} // namespace smd
