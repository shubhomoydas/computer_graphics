/*

Dummy routines for matrix transformations.

These are for you to write!

*/

#include "glDraw.h"

//-------------------------------------------------

smd::GraphicsManager gfxManager;

smd::GraphicsContext& getGraphicsContext() {
	return gfxManager.getGraphicsContext();
}

smd::GraphicsManager& getGraphicsManager() {
	return gfxManager;
}

void osuOrtho(double left, double right, double bottom, double top, double nearp,double farp)
{ 
	gfxManager.setMorth(left, right, bottom, top, nearp, farp);
}
void osuPerspective(double fovy, double nearp, double farp) 
{  	
	gfxManager.setMpers(fovy, nearp, farp);

	if (0) { // 0 - disable, 1 - enable
		// this is a test code for lookAt() since no other
		// test case was supplied.
		// IMPORTANT:
		// ----------
		// Use this only with commandline argument = 10, else disable
		double from[3] = { 0, 0, -75};
		double at[3]   = { 0.25, 0.05, 74};
		double up[3]   = {    1, 0, 0};
		gfxManager.lookAt(from, at, up);
	}
}

void osuBegin(OSUDrawable mode)
{
	gfxManager.begin((OSUDrawable)mode);
}
void osuColor3f(double red, double green, double blue)
{
	gfxManager.setColor(red, green, blue);
}
void osuVertex2f(double x, double y)
{
	gfxManager.addVertex(x,y);
}
void osuEnd()
{
	gfxManager.end();
}

void osuVertex3f(double x, double y, double z)
{
	gfxManager.addVertex(x,y,z);
}

void osuInitialize() 
{ 
	gfxManager.initialize();
}

void osuPushMatrix() 
{ 
	gfxManager.pushMatrix();
}

void osuPopMatrix() 
{ 
	gfxManager.popMatrix();
}

void osuLoadIdentityMatrix()
{
	gfxManager.loadIdentity();
}

void osuTranslate(double tx, double ty, double tz) 
{ 
	gfxManager.translate(tx, ty, tz);
}

void osuScale(double sx, double sy, double sz) 
{ 
	gfxManager.scale(sx, sy, sz);
}

void osuRotate(double angle, double ax, double ay, double az) { 
	gfxManager.rotate(angle, ax, ay, az);
}

void osuLookat(float from[3], float at[3], float up[3])
{
	double _from[3]; floatToDouble3(from,_from);
	double _at[3]; floatToDouble3(at,_at);
	double _up[3]; floatToDouble3(up,_up);
	gfxManager.lookAt(_from, _at, _up);
}

// Project 3

void setGraphicsDebug(bool debug, bool debugFine) {
	gfxManager.setDebug(debug, debugFine);
}

void osuNormal3f(double x, double y, double z) {
	gfxManager.setNormal(x, y, z);
}

void osuEnable(int depthTestBit) {
	gfxManager.enable(depthTestBit);
}

void osuClearZ() {
	gfxManager.clearZbuffer();
}

void osuShadeModel(int model) {
	gfxManager.setShadeModel(model);
}

void osuPointLight(float pos[3], float i) {
	gfxManager.addPointLight(pos, i);
}

void osuDirectionalLight(float dir[3], float i) {
	gfxManager.addDirectionalLight(dir, i);
}

void osuAmbientLight(float i) {
	gfxManager.setAmbientLight(i);
}

void osuDiffuse(float r, float g, float b) {
	gfxManager.setDiffuse(r, g, b);
}

void osuSpecular(float r, float g, float b, float s) {
	gfxManager.setSpecular(r, g, b, s);
}

void osuEmittance(float r, float g, float b) {
	gfxManager.setEmittance(r, g, b);
}

void osuSphere(double x, double y, double z, double r) {
	gfxManager.sphere(x, y, z, r);
}

void osuReflectionDepth(int depth) {
	gfxManager.setReflectionDepth(depth);
}

void osuReflectance(double ks) {
	gfxManager.setReflectance(ks);
}
