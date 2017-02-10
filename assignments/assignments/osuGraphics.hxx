/*==================================================================*/
/* CS551 Graphics Support Code					    */
/* Header File							    */
/* Version 3.0                                                      */
/*==================================================================*/

#ifndef OSUGFXSUPPORT
#define OSUGFXSUPPORT

#include <stdio.h>

/*------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif
  
/*==================================================================*/
/* osuImage Data Type Definition                                     */
/*------------------------------------------------------------------*/

typedef struct {
   /* an image structure holds the image size */
   /* and a dynamic array of pixels */
   int w, h;
   int *values;
} osuImage;

void osuImageInit ( osuImage *I );
void osuImageDestroy ( osuImage *I );

void osuImageGetSize ( osuImage *I, int *w, int *h );
void osuImageSetSize ( osuImage *I, int w, int h );


void osuImageWritePixel ( osuImage *I, int w, int h, int r, int g, int b );
void osuImageReadPixel ( osuImage *I, int w, int h, int *r, int *g, int *b );

/*==================================================================*/
/* A Single Window Interface
/*------------------------------------------------------------------*/

enum OSUWriteMode { OSU_OVERWRITE, OSU_REPLACE, OSU_XOR };
enum OSUColor { OSU_RED, OSU_GREEN, OSU_BLUE };
enum OSUEnable { OSU_DEPTH_TEST };

void osuBeginGraphics ( int w, int h );
void osuEndGraphics ();
void osuClear ( int r, int g, int b );
void osuFlush();
void osuWritePixel ( int x, int y, int r, int g, int b );
void osuSetWriteMode ( int /* OSUWriteMode */ mode );
void osuWaitOnEscape ();
void osuGetFramebufferSize ( int *w, int *h );

/*==================================================================*/
/* Declarations of polygon drawing routines
/*------------------------------------------------------------------*/

enum OSUDrawable { OSU_NONE, OSU_TRIANGLE ,OSU_LINES };

void osuColor3f ( double r, double g, double b );
void osuVertex2f ( double x, double y );
void osuBegin ( int /* OSUDrawable */ );
void osuEnd ();

/*------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif

/*------------------------------------------------------------------*/

#endif

