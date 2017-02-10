#include "glDraw.h"
#include "lines.h"

#define combineWithWeights(a,b,c,alpha,beta,gamma) (a*alpha+b*beta+c*gamma)

namespace smd {
	
	void vertexListToArray(std::list<Vertex>& vlist, Vertex* varr);
	void computeColorAtPoint(const int x, const int y, const Point& v, const Vector& n, Color& color, GraphicsContext& ctx);
	void combineVectorsWithWeights(const Vector& v1, const Vector& v2, const Vector& v3, 
		Vector& v, const double alpha, const double beta, const double gamma);
	void addLightColor(Color& col, const double c_l, Color& c_r, const Color& specular, 
		const Vector& n, const Vector& e, const Vector& l, const bool debug);
	void lin3x3(Vector& p1, Vector& p2, Vector& p3, Vector& p4, Vector& p5, double& beta, double& gamma, double& t);

	void poly_clip(
		Vertex *verts, int count, 
		Vertex *out_verts, int *out_count,
		double a, double b, double c, double d
	);

	void renderTriangle1(Vertex& v1, Vertex& v2, Vertex& v3, GraphicsContext& ctx);
	void renderTriangle2(Vertex& v1, Vertex& v2, Vertex& v3, GraphicsContext& ctx);
	
	void lin3x3(Vector& p1, Vector& p2, Vector& p3, Vector& p4, Vector& p5, double& beta, double& gamma, double& t) {

		/* Shirley pg. 79
			The intersection is solution of linear system:
			| xa-xb xa-xc xd | | beta  |   | xa-xe |
			| ya-yb ya-yc yd | | gamma | = | ya-ye |
			| za-zb za-zc zd | |   t   |   | za-ze |

		  Assume p4 == d and p5 == e
		  Solution of linear system: (Shirley pg. 79)
			| p1.x-p2.x p1.x-p3.x d.x | | beta  |   | p1.x-e.x |
			| p1.y-p2.y p1.y-p3.y d.y | | gamma | = | p1.y-e.y |
			| p1.z-p2.z p1.z-p3.z d.z | |   t   |   | p1.z-e.z |
			=>
			        | a d g | | beta  |   | j |
			        | b e h | | gamma | = | k |
			        | c f i | |   t   |   | l |
		*/
		double a = p1.x-p2.x, d = p1.x-p3.x, g = p4.x, j = p1.x-p5.x;
		double b = p1.y-p2.y, e = p1.y-p3.y, h = p4.y, k = p1.y-p5.y;
		double c = p1.z-p2.z, f = p1.z-p3.z, i = p4.z, l = p1.z-p5.z;

		// M = a(ei - hf) + b(gf - di) + c(dh - eg)
		double M = a*(e*i - h*f) + b*(g*f - d*i) + c*(d*h - e*g);

		// beta  =  (j(ei - hf) + k(gf - di) + l(dh - eg)) / M
		beta  =  (j*(e*i - h*f) + k*(g*f - d*i) + l*(d*h - e*g)) / M;
		// gamma =  (i(ak - jb) + h(jc - al) + g(bl - kc)) / M
		gamma =  (i*(a*k - j*b) + h*(j*c - a*l) + g*(b*l - k*c)) / M;
		// t     = -(f(ak - jb) + e(jc - al) + d(bl - kc)) / M
		t     = -(f*(a*k - j*b) + e*(j*c - a*l) + d*(b*l - k*c)) / M;

	}

	bool Polygon::hitTriangle(Point& p1, Point& p2, Point& p3, Ray& ray, HitData& hitdata, double t0, double t1) {
		bool hit = false;
		double beta = 0, gamma = 0, t = 0;
		lin3x3(p1, p2, p3, ray.d, ray.o, beta, gamma, t);
		if ((t < t0 || t > t1) || (gamma < 0 || gamma > 1) || (beta < 0 || beta > 1-gamma)) {
			hit = false;
		} else {
			hit = true;
			hitdata.nhits = 1;
			hitdata.t[0] = t;
			hitdata.alpha = 1-beta-gamma;
			hitdata.beta = beta;
			hitdata.gamma = gamma;
		}
		return hit;
	}

	bool Triangle::isHit(Ray& ray, HitData& hitdata, double t0, double t1) {
		if (!readyToRender()) return false;
		std::list<Vertex>::iterator iter = vertices.begin();
		Point& p1 = iter->v; iter++; // a
		Point& p2 = iter->v; iter++; // b
		Point& p3 = iter->v;         // c
		return hitTriangle(p1, p2, p3, ray, hitdata, t0, t1);
	}

	bool Sphere::isHit(Ray& ray, HitData& hitdata, double t0, double t1) {
		hitdata.nhits = 0;
		Vector emc; emc.set(ray.o).sub(c);
		double b = ray.d.dot(emc);
		double dd = ray.d.dot(ray.d);
		double discrim = b*b - dd*(emc.dot(emc) - r*r);
		if (discrim < 0) {
		} else if (discrim == 0)  {
			double t = -b/dd;
			if (t >= t0 && t <= t1) {
				hitdata.nhits = 1;
				hitdata.t[0] = t;
			}
		} else {
			double t = (-b - std::sqrt(discrim))/dd;
			if (t >= t0 && t <= t1) {
				hitdata.t[hitdata.nhits++] = t; // smaller (nearer) one
			}
			t = (-b + std::sqrt(discrim))/dd;
			if (t >= t0 && t <= t1) {
				hitdata.t[hitdata.nhits++] = t; // larger (farther) one
			}
		}
		return hitdata.nhits > 0;
	}

	void Line::render(GraphicsContext& ctx) {
		Point p1, p2;
	
		// get clip coordinates after perspective divide - should actually get before divide.
		Vector beforePersp;
		Vector tv1 = ctx.transformVectorToScreen(v1.v, v1.wCoords, v1.tEye, beforePersp, ctx.getViewPoint());
		Vector tv2 = ctx.transformVectorToScreen(v2.v, v2.wCoords, v2.tEye, beforePersp, ctx.getViewPoint());
		p1.set(tv1); //p1.mul(1/std::abs(p1.w)); // -|w| <= x,y,z <= |w|
		p2.set(tv2); //p2.mul(1/std::abs(p2.w)); // -|w| <= x,y,z <= |w|
		double nearp = -1, farp = 1;
		int bdrawline = near_far_clip(nearp, farp, &p1.x, &p1.y, &p1.z, &p2.x, &p2.y, &p2.z);
		if (bdrawline) {
			draw_line(p1.x, p1.y, p2.x, p2.y);
		} else {
			if (ctx.debug) {
				std::cout << "Line " << *this << " is outside view (" << nearp << "," << farp << ")" << std::endl;
			}
		}
	}
	
	void Polygon::initNormal() {
		initNormal(vertices);
	}
	void Polygon::initNormal(std::list<Vertex>& vertices) {
		normal.set(computeFaceNormal(vertices));
	}
	void Polygon::render(GraphicsContext& ctx) {
		std::list<Vertex> worldVertices;
		std::list<Vertex> tclipVertices;
		std::list<Vertex> clippedVertices;
		std::list<Vertex> screenVertices;

		transformVerticesFromModelToWorld(vertices, worldVertices, ctx);

		initNormal(worldVertices);
		if (VEC_ZERO.equals(normal)) ctx.zeroNormals++;

		setVertexNormals(worldVertices, normal, ctx);
		setVertexColors(worldVertices, ctx);
		
		if (0) {
			//std::cout << "specular=" << ctx.specular << "; diffuse=" << ctx.diffuse << 
			//	"; ambient=" << ctx.ambient << std::endl;
			//printVertexColors(vertices, std::cout); std::cout << std::endl;
			printVertices(worldVertices, std::cout); std::cout << std::endl;
			std::cout << "face normal=" << normal << std::endl;
		}

		// transform the vertices to clip coordinates but retain the
		// z coordinates from the eye(camera) transform. This will be
		// used for z-buffer depth test.
		bool valid = transformVerticesFromWorldToClip(worldVertices, tclipVertices, ctx);
		if (!valid && ctx.debug) {
			std::cout << "Invalid vertex after clipping. (Close to near plane?)" << std::endl;
			//std::cout << "M:" << std::endl << ctx.getViewPoint().getMprojection() << std::endl;
		}
		if (valid) getClippedPolygonVertices(tclipVertices, clippedVertices);
		if (0 && ctx.debug) {
			std::cout << "Polygon original ";
			printVertices(vertices, std::cout);
			std::cout << std::endl;
			std::cout << "Polygon before clipping ";
			printVertices(tclipVertices, std::cout);
			std::cout << std::endl;
			std::cout << "Polygon after clipping ";
			printVertices(clippedVertices, std::cout);
			std::cout << std::endl;
		}
		if (clippedVertices.size() < 3) {
			if (ctx.debug) std::cout << "Polygon outside view region" << std::endl;
			return;
		}
		transformVerticesFromClipToScreen(clippedVertices, screenVertices, ctx);
		//Vertex first = screenVertices.front();
		//if (screenVertices.size() > 3) screenVertices.push_back(first);
		Vertex *varr = new Vertex[screenVertices.size()];
		vertexListToArray(screenVertices, varr);
		// below routine only works for triangles and quadrilaterals.
		Vertex *v[3];
		v[0] = &varr[0]; v[1] = v[2] = 0;
		for (unsigned int i = 0; i < screenVertices.size(); i++) {
			v[1] = v[2]; v[2] = &varr[i];
			if (i > 1) {
				//if (ctx.debug) 
				//	std::cout << "Rendering Triangle: " << v[0]->v << "," << v[1]->v << "," << v[2]->v << std::endl;
				renderTriangle2(*v[0], *v[1], *v[2], ctx);
				//if (ctx.debug) 
				//	std::cout << "Rendering Complete..." << std::endl;
			}
		}
		delete[] varr;
	}

	Vector Polygon::computeFaceNormal(std::list<Vertex>& vertices) {
		if (vertices.size() < 3) return VEC_ZERO;
		Point p[3];
		int i = 0;
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end() && i < 3; iter++, i++) {
			p[i].set((*iter).v);
		}
		return smd::getNormal(p[0], p[1], p[2]);
	}
	
	void Polygon::setVertexColors(std::list<Vertex>& vertices, GraphicsContext& ctx) {
		if (ctx.rendererType == smd::RADIOSITY) return; // no needed for radiosity.
		// Below is for phong shading
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			// Set vertex colors without lighting effects
			computeColorAtPoint(-1, -1, (*iter).v, (*iter).normal, (*iter).color, ctx);
		}
		if (0) { // DEBUG code
			for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
				std::cout << "normal: " << (*iter).normal << std::endl;
			}
		}
	}

	void Polygon::setVertexNormals(std::list<Vertex>& vertices, Vector& n, GraphicsContext& ctx) {
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
			if (!(*iter).priorNormal) {
				(*iter).setNormal(n);
				ctx.computedVertexNormals++;
			} else {
				ctx.priorVertexNormals++;
			}
		}
		if (0) { // DEBUG code
			for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++) {
				std::cout << "normal: " << (*iter).normal << std::endl;
			}
		}
	}

	void Polygon::getClippedPolygonVertices(std::list<Vertex>& vertices, std::list<Vertex>& clippedVertices) {
		Vertex* vertsA = new Vertex[vertices.size()*2];
		Vertex* vertsB = new Vertex[vertices.size()*2];
		Vertex* verts_in, *verts_out, *verts_temp;
		int countA, countB;

		countA = vertices.size(); countB = 0; verts_in = vertsA; verts_out = vertsB;

		int i = 0;
		for (std::list<Vertex>::iterator iter = vertices.begin(); iter != vertices.end(); iter++, i++) {
			//copyToOsuVertex(*iter, verts_in[i]);
			verts_in[i].set(*iter);
		}

		//std::cout << "In count 1: " << countA << std::endl;

		//  x + 1 = 0 -- left plane
		poly_clip(verts_in, countA, verts_out, &countB,  1,  0,  0,  1);
		countA = countB; verts_temp = verts_in; verts_in = verts_out; verts_out = verts_temp;
		//std::cout << "Out count 1: " << countA << std::endl;

		// -x + 1 = 0 -- right plane
		poly_clip(verts_in, countA, verts_out, &countB, -1,  0,  0,  1);
		countA = countB; verts_temp = verts_in; verts_in = verts_out; verts_out = verts_temp;
		//std::cout << "Out count 2: " << countA << std::endl;

		//  y + 1 = 0 -- bottom plane
		poly_clip(verts_in, countA, verts_out, &countB,  0,  1,  0,  1);
		countA = countB; verts_temp = verts_in; verts_in = verts_out; verts_out = verts_temp;
		//std::cout << "Out count 3: " << countA << std::endl;

		// -y + 1 = 0 -- top plane
		poly_clip(verts_in, countA, verts_out, &countB,  0, -1,  0,  1);
		countA = countB; verts_temp = verts_in; verts_in = verts_out; verts_out = verts_temp;
		//std::cout << "Out count 4: " << countA << std::endl;

		// -z - 1 = 0 -- front plane
		poly_clip(verts_in, countA, verts_out, &countB,  0,  0,  1,  1);
		countA = countB; verts_temp = verts_in; verts_in = verts_out; verts_out = verts_temp;
		//std::cout << "Out count 5: " << countA << std::endl;

		//  z - 1 = 0 -- back plane
		poly_clip(verts_in, countA, verts_out, &countB,  0,  0, -1,  1);
		countA = countB; verts_temp = verts_in; verts_in = verts_out; verts_out = verts_temp;
		//std::cout << "Out count 6: " << countA << std::endl;

		// at this point verts_in has the final clipped vertices
		clippedVertices.clear();
		for (i = 0; i < countA; i++) {
			clippedVertices.push_back(verts_in[i]);
		}
		delete[] vertsA;
		delete[] vertsB;
	}

	// DEBUG Code >>>>>>>>
	#define DBG_TOL_X 1e-2
	#define DBG_TOL_Y 1e-2
	#define DBG_TOL_Z 1e-2
	bool insideTestArea(const Vector& vWorld) {
		if (std::abs(vWorld.x-0.30) <= DBG_TOL_X 
			&& std::abs(vWorld.y-(-1)) <= DBG_TOL_Y 
			&& std::abs(vWorld.z-(-0.04)) <= DBG_TOL_Z) {
			return true;
		}
		return false;
	}
	// <<<<<<< DEBUG Code

	/**
	 * Populates the view and light source zbuffers and returns
	 * whether the pixel should be rendered in the view.
	 */
	bool testAndPopulateZbuffer(int x, int y, Point& vWorld, GraphicsContext& ctx) {
		// We will use the eye coordinates for z-buffer test
		Vector vEye; ctx.getViewPoint().worldToEye(vWorld, vEye);
		return ctx.getViewPoint().testAndPopulateZbuffer(x, y, vEye.z, true);
	}

	void addLightColor(Color& col, const double c_l, Color& c_r, const Color& specular, 
		const Vector& n, const Vector& e, const Vector& l, const bool debug) {
		if (debug) {
			std::cout << "e: " << e << "; l: " << l << "; n: " << n << std::endl;
		}
		if (n.dot(e) < 0) return; // vertex not visible to eye
		double n_l = n.dot(l);
		double max_0_nl = std::max(0.0,n_l);
		Vector r = n; r.mul(2*n_l).sub(l);
		double phong = std::pow(std::max(0.0, e.dot(r)),specular.s);
		if (debug) {
			std::cout << "n.l=" << n.dot(e) << "; phong=" << phong << std::endl;
		}
		col.r = col.r + c_l*(c_r.r*max_0_nl + specular.r*phong);
		col.g = col.g + c_l*(c_r.g*max_0_nl + specular.g*phong);
		col.b = col.b + c_l*(c_r.b*max_0_nl + specular.b*phong);
	}
	
	void computeColorAtPoint(const int x, const int y, const Point& vWorld, 
		const Vector& n, Color& color, GraphicsContext& ctx) {

			if (ctx.renderMode == smd::RENDER_SHADOW && ctx.pointLights.size() > 0 && !(x < 0 || y < 0)) {
				Vector vEye, vClip, vImage;
				ctx.getViewPoint().worldToScreen(vWorld, vEye, vClip, vImage);
				//bool visible = ctx.getViewPoint().testAndPopulateZbuffer(x, y, vEye.z, false);
				bool visible = ctx.pointLights[0]->testAndPopulateZbuffer(vWorld, false);
				if (0 && insideTestArea(vWorld)) {
					std::cout << "Coloring at vWorld=" << vWorld << 
						"; visible=" << (visible? "true":"false") << 
						"; x=" << x << "; y=" << y << std::endl << 
						"; vEye=" << vEye << "; vClip=" << vClip << "; vImage=" << vImage << std::endl;
				}
				if (visible) {
					color.set(1,1,1);
				}
				return;
			}

		double c_a = ctx.ambient.s;
		Color c_r = ctx.diffuse;
		Vector e = ctx.getViewPoint().e; e.sub(vWorld).normalizeCoords();
		color.r = c_r.r*c_a;
		color.g = c_r.g*c_a;
		color.b = c_r.b*c_a;

		if (!(x < 0 || y < 0)) {
			for (unsigned int i = 0; i < ctx.directionalLights.size(); i++) {
				double c_l = ctx.directionalLights[i]->i;
				addLightColor(color, c_l, c_r, ctx.specular, n, e, ctx.directionalLights[i]->e, false);
			}

			for (unsigned int i = 0; i < ctx.pointLights.size(); i++) {
				bool visible = true;
				visible = ctx.pointLights[i]->testAndPopulateZbuffer(vWorld, false);
				if (visible) {
					Vector li = ctx.pointLights[i]->e; li.sub(vWorld).normalizeCoords();
					addLightColor(color, ctx.pointLights[i]->i, ctx.diffuse, ctx.specular, n, e, li, false);
				}
			}
		}
	}

	void colorPixelPhong(int x, int y, Vector vWorld, Vector n, GraphicsContext& ctx) {
		Color color;
		computeColorAtPoint(x, y, vWorld, n, color, ctx);
		osuWritePixel(x,y,round(color.r*255),round(color.g*255),round(color.b*255));
	}

	void colorPixel(int x, int y, Vertex& v1, Vertex& v2, Vertex& v3, 
		double alpha, double beta, double gamma, GraphicsContext& ctx) {
		Color color;
		combineColors(v1.color, v2.color, v3.color, color, alpha, beta, gamma, 255.0*ctx.ambient.s);
		osuWritePixel(x,y,round(color.r),round(color.g),round(color.b));
	}

	void renderTriangle2(Vertex& v1, Vertex& v2, Vertex& v3, GraphicsContext& ctx) {

		Point p1; p1.set(v1.v);
		Point p2; p2.set(v2.v);
		Point p3; p3.set(v3.v);

		double xmin=std::floor(min3(p1.x,p2.x,p3.x));
		double xmax =std::ceil(max3(p1.x,p2.x,p3.x));
		double ymin=std::floor(min3(p1.y,p2.y,p3.y));
		double ymax =std::ceil(max3(p1.y,p2.y,p3.y));

		//std::cout << "(" << xmin << "," << ymin << ") - (" << xmax << "," << ymax << ")" << std::endl;
	
		ImplicitLinear f12(p1,p2), f23(p2,p3), f31(p3,p1);
		double f_alpha = f23(p1.x,p1.y);
		double f_beta = f31(p2.x,p2.y);
		double f_gamma = f12(p3.x,p3.y);

		double o_f23 = 0, o_f31 = 0, o_f12 = 0;

		Point op; // Some offscreen point that is not on any triangle edge.
		op.set(-20,-20, 0); // Start searching here
		int i = 0;
		while ((o_f23 == 0 || o_f31 == 0 || o_f12 == 0) && i < 20) {
			// we want to find an offscreen point that does not
			// lie on any of the edges. Try atmost 20 times.
			int m1 = std::rand()%2; m1 = m1 ? m1 : -1;
			int m2 = std::rand()%2; m2 = m2 ? m2 : -1;
			op.x = std::min(-1.0, op.x + m1 * (std::rand()%10));
			op.y = std::min(-1.0, op.y + m2 * (std::rand()%10));
			o_f23 = f23(op.x,op.y);
			o_f31 = f31(op.x,op.y);
			o_f12 = f12(op.x,op.y);
			i++;
		}
		//std::cout << "Offscreen point is " << op << ", found in " << i << " iterations" <<  std::endl;
		for (double y = ymin; y <= ymax; y=y+1) {
			for (double x = xmin; x <= xmax; x=x+1) {
				double alpha = f23(x,y)/f_alpha;
				double beta  = f31(x,y)/f_beta;
				double gamma = f12(x,y)/f_gamma; // gamma = 1-alpha-beta; does not work well
				if (alpha >= 0 && beta >= 0 && gamma >= 0) {
					if ((alpha > 0 || f_alpha * o_f23 > 0) 
							&& (beta > 0 || f_beta * o_f31 > 0) 
							&& (gamma > 0 || f_gamma * o_f12 > 0)) {

						int ix = round(x), iy = round(y);

						// get the world coordinates visible from this pixel
						double z = combineWithWeights(v1.v.z,v2.v.z,v3.v.z,alpha,beta,gamma);
						Vector vWorld, vImage; vImage.setLoc(x, y, z);
						ctx.getViewPoint().screenToWorld(vImage, vWorld);
						if (0 && ctx.debug) {
							Vector vScreen, vEye, vClip;
							ctx.getViewPoint().worldToScreen(vWorld, vEye, vClip, vScreen);
							std::cout << "ix=" << ix << "; iy=" << iy << "; vWorld=" << vWorld << "; vScreen=" << vScreen << std::endl;
						}

						bool draw = false;
						if (ctx.isDepthEnabled()) {
							draw = testAndPopulateZbuffer(ix, iy, vWorld, ctx);
						} else {
							draw = true;
						}

						if (draw) {
							if (ctx.rendererType == smd::RADIOSITY) {
								colorPixel(ix, iy, v1, v2, v3, alpha, beta, gamma, ctx);
							} else if ((ctx.isRenderColor() || ctx.renderMode == smd::RENDER_SHADOW)) {
								Vector n; combineVectorsWithWeights(v1.normal, v2.normal, v3.normal, n, alpha, beta, gamma);
								colorPixelPhong(ix, iy, vWorld, n, ctx);
							}
						}

					}
				}
			}
		}
	}

	void renderTriangle1(Vertex& v1, Vertex& v2, Vertex& v3, GraphicsContext& ctx) {

		Point p1; p1.set(v1.v);
		Point p2; p2.set(v2.v);
		Point p3; p3.set(v3.v);

		//std::cout << "Drawing triangle p1=" << p1 << "; p2=" << p2 << "; p3=" << p3 << std::endl;

		double xmin=std::floor(min3(p1.x,p2.x,p3.x));
		double xmax =std::ceil(max3(p1.x,p2.x,p3.x));
		double ymin=std::floor(min3(p1.y,p2.y,p3.y));
		double ymax =std::ceil(max3(p1.y,p2.y,p3.y));

		ImplicitLinear f12(p1,p2), f23(p2,p3), f31(p3,p1);
		Color color;
		double f_alpha = f23(p1.x,p1.y);
		double f_beta = f31(p2.x,p2.y);
		double f_gamma = f12(p3.x,p3.y);
		for (double y = ymin; y <= ymax; y=y+1) {
			for (double x = xmin; x <= xmax; x=x+1) {
				double alpha = f23(x,y)/f_alpha;
				double beta  = f31(x,y)/f_beta;
				double gamma = f12(x,y)/f_gamma;
				if (alpha > 0 && beta > 0 && gamma > 0) {
					//std::cout << "(" << x << "," << y << ")" << std::endl;
					combineColors(v1.color, v2.color, v3.color, color, alpha, beta, gamma, 255.0);
					osuWritePixel(round(x),round(y),round(color.r),round(color.g),round(color.b));
				}
			}
		}
	}

	/******************************************************************************
	Create a new vertex that is the intersection between a plane and a line
	segment between two given vertices.

	Entry:
	  v0,v1 - two vertex endpoints of the line segment
	  a,b,c,d - coefficients for the plane ax + by + cz + d = 0

	Exit:
	  vnew - the new vertex at the intersection between the line and plane
	******************************************************************************/

	void create_vertex (
		Vertex *v0,
		Vertex *v1,
		Vertex *newv,
		double a,
		double b,
		double c,
		double d
		)
	{
		double t;
		double x0,y0,z0;
		double x1,y1,z1;
		double dx,dy,dz;

		/* shorthands */

		x0 = v0->v.x;
		y0 = v0->v.y;
		z0 = v0->v.z;

		x1 = v1->v.x;
		y1 = v1->v.y;
		z1 = v1->v.z;

		dx = x1 - x0;
		dy = y1 - y0;
		dz = z1 - z0;

		/* find parameter t saying how far between v0 and v1 the intersection is */

		t = -1.0 * (a*x0 + b*y0 + c*z0 + d) / (a*dx + b*dy + c*dz);

		/* interpolate between values in v0 and v1 for location and color */

		newv->v.x = x0 + t * (x1 - x0);
		newv->v.y = y0 + t * (y1 - y0);
		newv->v.z = z0 + t * (z1 - z0);
		newv->v.w = 1; // this is a location

		newv->wCoords.x = v0->wCoords.x + t * (v1->wCoords.x - v0->wCoords.x);
		newv->wCoords.y = v0->wCoords.y + t * (v1->wCoords.y - v0->wCoords.y);
		newv->wCoords.z = v0->wCoords.z + t * (v1->wCoords.z - v0->wCoords.z);
		newv->wCoords.w = 1; // this is a location

		newv->color.r = v0->color.r + t * (v1->color.r - v0->color.r);
		newv->color.g = v0->color.g + t * (v1->color.g - v0->color.g);
		newv->color.b = v0->color.b + t * (v1->color.b - v0->color.b);

		newv->normal.x = v0->normal.x + t * (v1->normal.x - v0->normal.x);
		newv->normal.y = v0->normal.y + t * (v1->normal.y - v0->normal.y);
		newv->normal.z = v0->normal.z + t * (v1->normal.z - v0->normal.z);
		newv->normal.w = 0; // this is a vector
	}


	/******************************************************************************
	Clip a polygon to a plane.

	Entry:
	  verts   - vertices of polygon to clip
	  count   - number of vertices in polygon
	  a,b,c,d - coefficients of plane equation against which to clip:
			positive side described by ax + by + cz + d > 0 are kept

	Exit:
	  out_verts - vertices of clipped polygon
	  out_count - number of vertices in the clipped polygon, or 0 if the entire
			  polygon is on the wrong side of the clipping plane
	******************************************************************************/

	void poly_clip (
		Vertex *verts,
		int count,
		Vertex *out_verts,
		int *out_count,
		double a,
		double b,
		double c,
		double d
		)
	{
		int i; //,ii;
		int new_count = 0;
		Vertex *v0,*v1;
		int in0,in1;  /* are v0 or v1 in the proper half-space */

		v0 = &verts[0];
		in0 = (a * v0->v.x + b * v0->v.y + c * v0->v.z + d > 0);

		for (i = 0; i < count; i++) {

			v0 = &verts[i];
			v1 = &verts[(i+1) % count];
			in1 = (a * v1->v.x + b * v1->v.y + c * v1->v.z + d > 0);

			if (in0 && in1) {
				out_verts[new_count++].set(*v1);
			}
			else if (!in0 && in1) {
				create_vertex (v0, v1, &out_verts[new_count++], a, b, c, d);
				out_verts[new_count++].set(*v1);
			}
			else if (in0 && !in1) {
				create_vertex (v0, v1, &out_verts[new_count++], a, b, c, d);
			}
			else {
				/* both are not in, so we add no vertices to the clipped polygon */
			}

			in0 = in1;
		}

		*out_count = new_count;
	}

	void vertexListToArray(std::list<Vertex>& vlist, Vertex* varr) {
		int i = 0;
		for (std::list<Vertex>::iterator it = vlist.begin(); it != vlist.end(); it++, i++) {
			varr[i].set(*it);
		}
	}

	void combineVectorsWithWeights(const Vector& v1, const Vector& v2, const Vector& v3, 
		Vector& v, const double alpha, const double beta, const double gamma) {
		v.x = (alpha * v1.x + beta * v2.x + gamma * v3.x);
		v.y = (alpha * v1.y + beta * v2.y + gamma * v3.y);
		v.z = (alpha * v1.z + beta * v2.z + gamma * v3.z);
		v.w = (alpha * v1.w + beta * v2.w + gamma * v3.w);
	}

	void combineColors(const Color& c1, const Color& c2, const Color& c3, Color& dest, 
					   const double alpha, const double beta, const double gamma, const double scale) {
		dest.r = scale * (alpha * c1.r + beta * c2.r + gamma * c3.r);
		dest.g = scale * (alpha * c1.g + beta * c2.g + gamma * c3.g);
		dest.b = scale * (alpha * c1.b + beta * c2.b + gamma * c3.b);
		dest.s = scale * (alpha * c1.s + beta * c2.s + gamma * c3.s);
	}

	std::ostream& operator<< (std::ostream& out, Triangle& triangle) {
		return triangle.print(out);
	}

	std::ostream& operator<< (std::ostream& out, Polygon& polygon) {
		return polygon.print(out);
	}

	std::ostream& operator<< (std::ostream& out, Sphere& sphere) {
		return sphere.print(out);
	}

	std::ostream& operator<< (std::ostream& out, Line& line) {
		return line.print(out);
	}

} //namespace smd

