#include "glDraw.h"

/*

Line and matrix header.

*/

void draw_line(double, double, double, double);
int near_far_clip(double,double,double*,double*,double*,double*,double*,double*);
void osuInitialize();
void osuPushMatrix();
void osuPopMatrix();
void osuLoadIdentityMatrix();
void osuTranslate(double, double, double);
void osuScale(double, double, double);
void osuRotate(double angle, double ax, double ay, double az);
void osuOrtho(double, double, double, double, double, double);
void osuPerspective(double, double, double);
void osuLookat(float from[3], float at[3], float up[3]); // modified for project3
void osuBegin(OSUDrawable);
void osuEnd();
void osuVertex3f(double, double, double);

// Project 3
void osuNormal3f(double x, double y, double z);
void osuEnable(int depthTestBit);
void osuClearZ();
void osuShadeModel(int model);
void osuPointLight(float pos[3], float i);
void osuDirectionalLight(float dir[3], float i);
void osuAmbientLight(float i);
void osuDiffuse(float r, float g, float b);
void osuSpecular(float r, float g, float b, float s);

void osuSphere(double x, double y, double z, double r);
void osuReflectionDepth(int depth);
void osuReflectance(double ks);
void osuEmittance(float r, float g, float b);

smd::GraphicsManager& getGraphicsManager();
smd::GraphicsContext& getGraphicsContext();
void setGraphicsDebug(bool debug, bool debugFine);
