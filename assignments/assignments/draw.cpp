/*
Line drawing in 3D.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "osuGraphics.h"
#include "matlib.h"
#include "lines.h"
#include <Windows.h>

int xsize = 360;
int ysize = 360;

double rand01()
{
	return ((double)rand()/RAND_MAX);
}

/******************************************************************************
Draw some random lines.
******************************************************************************/

void line_test()
{
  int i;
  float x0,y0,x1,y1;

  /* draw a bunch of random clipped lines*/

  for (i = 0; i < 200; i++) {
    x0 = (float) (600 * rand01()-120);
    y0 = (float) (600 * rand01()-120);
    x1 = (float) (600 * rand01()-120);
    y1 = (float) (600 * rand01()-120);
    draw_line (x0, y0, x1, y1);
  }
}


/******************************************************************************
Test the orthographic projection routine.
******************************************************************************/

void ortho_test()
{
  float nearp = 10.0; 
  float farp = 40.0;

  osuInitialize();

  osuOrtho (-100.0, 100.0, -100.0, 100.0, nearp, farp);

  osuBegin (OSU_LINES);

  /* three axes */
  osuVertex3f (0.0, 0.0, -20.0);
  osuVertex3f (50.0, 0.0, -20.0);

  osuVertex3f (0.0, 0.0, -20.0);
  osuVertex3f (0.0, 50.0, -20.0);

  osuVertex3f (0.0, 0.0, -20.0);
  osuVertex3f (0.0, 0.0, 30.0);

  /* test clipping */
  osuVertex3f (0.0, 25.0, 0.0);
  osuVertex3f (50.0, 25.0, -50.0);

  /* a square */
  osuVertex3f (-10.0, -10.0, -20.0);
  osuVertex3f (-90.0, -10.0, -20.0);

  osuVertex3f (-90.0, -10.0, -20.0);
  osuVertex3f (-90.0, -90.0, -20.0);

  osuVertex3f (-90.0, -90.0, -20.0);
  osuVertex3f (-10.0, -90.0, -20.0);

  osuVertex3f (-10.0, -90.0, -20.0);
  osuVertex3f (-10.0, -10.0, -20.0);

  osuEnd();

}


/******************************************************************************
Draw a cube.
******************************************************************************/

void cube()
{
  osuBegin (OSU_LINES);

  /* top square */

  osuVertex3f (-1.0, -1.0,  1.0);
  osuVertex3f (-1.0,  1.0,  1.0);

  osuVertex3f (-1.0,  1.0,  1.0);
  osuVertex3f ( 1.0,  1.0,  1.0);

  osuVertex3f ( 1.0,  1.0,  1.0);
  osuVertex3f ( 1.0, -1.0,  1.0);

  osuVertex3f ( 1.0, -1.0,  1.0);
  osuVertex3f (-1.0, -1.0,  1.0);

  /* bottom square */

  osuVertex3f (-1.0, -1.0, -1.0);
  osuVertex3f (-1.0,  1.0, -1.0);

  osuVertex3f (-1.0,  1.0, -1.0);
  osuVertex3f ( 1.0,  1.0, -1.0);

  osuVertex3f ( 1.0,  1.0, -1.0);
  osuVertex3f ( 1.0, -1.0, -1.0);

  osuVertex3f ( 1.0, -1.0, -1.0);
  osuVertex3f (-1.0, -1.0, -1.0);

  /* connect top to bottom */

  osuVertex3f (-1.0, -1.0, -1.0);
  osuVertex3f (-1.0, -1.0,  1.0);

  osuVertex3f (-1.0,  1.0, -1.0);
  osuVertex3f (-1.0,  1.0,  1.0);

  osuVertex3f ( 1.0,  1.0, -1.0);
  osuVertex3f ( 1.0,  1.0,  1.0);

  osuVertex3f ( 1.0, -1.0, -1.0);
  osuVertex3f ( 1.0, -1.0,  1.0);

  osuEnd();
}


/******************************************************************************
Orthographic cube.
******************************************************************************/

void ortho_cube()
{
  osuOrtho (-2.0, 2.0, -2.0, 2.0, 0.000, 10000.0);

  osuPushMatrix();
  osuTranslate (0.0, 0.0, -4.0);
  cube();
  osuPopMatrix();
}


/******************************************************************************
Pespective cube, no clipping.
******************************************************************************/

void cube_no_clip()
{
  osuPerspective (60.0, 0.000, 10000.0);

  osuPushMatrix();
  osuTranslate (0.0, 0.0, -4.0);
  cube();
  osuPopMatrix();
}


/******************************************************************************
Pespective cube, clipped at near end.
******************************************************************************/

void cube_clip_near()
{
  osuPerspective (60.0, 3.5, 10000.0);
  osuPushMatrix();
  osuTranslate (0.0, 0.0, -4.0);
  cube();
  osuPopMatrix();
}


/******************************************************************************
Pespective cube, clipped at far end.
******************************************************************************/

void cube_clip_far()
{
  osuPerspective (60.0f, 2.0f, 4.8f);

  osuPushMatrix();
  osuTranslate (0.0, 0.0, -4.0);
  cube();
  osuPopMatrix();
}


/******************************************************************************
Draw X, Y, Z axes.
******************************************************************************/

void axes()
{
  /* x, y, and z axes */

  osuBegin (OSU_LINES);

  osuVertex3f (0.0, 0.0, 0.0);
  osuVertex3f (50.0, 0.0, 0.0);

  osuVertex3f (0.0, 0.0, 0.0);
  osuVertex3f (0.0, 50.0, 0.0);

  osuVertex3f (0.0, 0.0, 0.0);
  osuVertex3f (0.0, 0.0, 50.0);

  osuEnd();

  /* the letter X */

  osuPushMatrix();

  osuTranslate (60.0, 0.0, 0.0);
  osuScale (5.0, 5.0, 5.0);

  osuBegin (OSU_LINES);
  osuVertex3f (1.0, 1.0, 0.0);
  osuVertex3f (-1.0, -1.0, 0.0);
  osuVertex3f (1.0, -1.0, 0.0);
  osuVertex3f (-1.0, 1.0, 0.0);
  osuEnd();

  osuPopMatrix();

  /* the letter Y */

  osuPushMatrix();

  osuTranslate (0.0, 60.0, 0.0);
  osuScale (5.0, 5.0, 5.0);

  osuBegin (OSU_LINES);
  osuVertex3f (1.0, 1.0, 0.0);
  osuVertex3f (0.0, 0.0, 0.0);
  osuVertex3f (-1.0, 1.0, 0.0);
  osuVertex3f (0.0, 0.0, 0.0);
  osuVertex3f (0.0, 0.0, 0.0);
  osuVertex3f (0.0, -1.0, 0.0);
  osuEnd();

  osuPopMatrix();

  /* the letter Z */

  osuPushMatrix();

  osuTranslate (0.0, 0.0, 60.0);
  osuScale (5.0, 5.0, 5.0);

  osuBegin (OSU_LINES);
  osuVertex3f (1.0, 1.0, 0.0);
  osuVertex3f (-1.0, -1.0, 0.0);
  osuVertex3f (-1.0, 1.0, 0.0);
  osuVertex3f (1.0, 1.0, 0.0);
  osuVertex3f (-1.0, -1.0, 0.0);
  osuVertex3f (1.0, -1.0, 0.0);
  osuEnd();

  osuPopMatrix();

  /* three cubes along the -z axis */

  osuPushMatrix();
  osuTranslate (0.0, 0.0, -100.0);
  osuScale (10.0, 10.0, 10.0);
  cube();
  osuPopMatrix();

  osuPushMatrix();
  osuTranslate (0.0, 0.0, -200.0);
  osuScale (10.0, 10.0, 10.0);
  cube();
  osuPopMatrix();

  osuPushMatrix();
  osuTranslate (0.0, 0.0, -400.0);
  osuScale (10.0, 10.0, 10.0);
  cube();
  osuPopMatrix();
}


/******************************************************************************
Draw axes orthographically.
******************************************************************************/

void ortho_axes()
{
 
  osuOrtho (-100.0, 100.0, -100.0, 100.0, 100.0, 1000.0);
  osuPushMatrix();
  osuRotate (45.0, 1.0, -1.0, 0.0);
  osuTranslate (-100.0, -100.0, -150.0);
  axes();
  osuPopMatrix();
}


/******************************************************************************
Draw axes in perspective.
******************************************************************************/

void persp_axes()
{
  osuPerspective (60.0, 0.000, 1000.0);

  osuPushMatrix();
  osuRotate (45.0, 1.0, -1.0, 0.0);
  osuTranslate (-100.0, -100.0, -150.0);
  axes();
  osuPopMatrix();
}


/******************************************************************************
Draw a circle of unit radius.
******************************************************************************/

void circle()
{
  int i;
  float theta;
  float x0,y0,x1,y1;
  float steps = 50;

  osuBegin (OSU_LINES);

  x0 = 1.0;
  y0 = 0.0;
  for (i = 0; i <= steps; i++) {
    theta = 2 * 3.1415926535f * i / steps;
    x1 = (float)cos (theta);
    y1 = (float)sin (theta);
    osuVertex3f (x0, y0, 0.0);
    osuVertex3f (x1, y1, 0.0);
    x0 = x1;
    y0 = y1;
  }

  osuEnd();
}


/******************************************************************************
Draw a face.
******************************************************************************/

void face()
{
  /* head */

  osuPushMatrix();
  osuTranslate (0.5, 0.5, -1.0);
  osuScale (0.4f, 0.4f, 1.0f);
  circle();
  osuPopMatrix();

  /* right eye */

  osuPushMatrix();
  osuTranslate (0.7f, 0.7f, -1.0f);
  osuScale (0.1f, 0.1f, 1.0f);
  circle();
  osuPopMatrix();

  /* left eye */

  osuPushMatrix();
  osuTranslate (0.3f, 0.7f, -1.0f);
  osuScale (0.1f, 0.1f, 1.0f);
  circle();
  osuPopMatrix();

  /* nose */

  osuPushMatrix();
  osuTranslate (0.5f, 0.5f, -1.0f);
  osuScale (0.07f, 0.07f, 1.0f);
  circle();
  osuPopMatrix();

  /* mouth */

  osuPushMatrix();
  osuTranslate (0.5f, 0.25f, -1.0f);
  osuScale (0.2f, 0.1f, 1.0f);
  circle();
  osuPopMatrix();
}


/******************************************************************************
Test the matrix stack by drawing a face.
******************************************************************************/

void face_test()
{
  float nearp = 0.0;
  float farp = 100000.0;

  osuInitialize ();

  osuOrtho (0.0, 1.0, 0.0, 1.0, nearp, farp);

  face();
}


/******************************************************************************
Draw four faces.
******************************************************************************/

void faces()
{
  float nearp = 0.0;
  float farp = 100000.0;

  osuInitialize ();

  osuOrtho (0.0, 1.0, 0.0, 1.0, nearp, farp);

  osuPushMatrix();
  osuTranslate (0.75, 0.25, 0.0);
  osuScale (0.5, 0.5, 1.0);
  osuTranslate (-0.5, -0.5, 0.0);
  face();
  osuPopMatrix();

  osuPushMatrix();
  osuTranslate (0.25, 0.25, 0.0);
  osuScale (0.5, 0.5, 1.0);
  osuTranslate (-0.5, -0.5, 0.0);
  face();
  osuPopMatrix();

  osuPushMatrix();
  osuTranslate (0.75, 0.75, 0.0);
  osuScale (0.5, 0.5, 1.0);
  osuTranslate (-0.5, -0.5, 0.0);
  face();
  osuPopMatrix();

  osuPushMatrix();
  osuTranslate (0.25, 0.75, 0.0);
  osuScale (0.5, 0.5, 1.0);
  osuRotate (30.0, 0.0, 0.0, 1.0);
  osuTranslate (-0.5, -0.5, 0.0);
  face();
  osuPopMatrix();
}

/*****************************************************************

cycle through all of the test cases by hitting the 'n' key.  This
can be invoked at any time.  It will always start over at test 0
and cycle through until the last..then quit.
******************************************************************/
void all()
{
	static int count = 1;
	osuInitialize();
	osuClear(0,0,0);
	switch(count){
		case 1:
			line_test();
			break;
		case 2:
			ortho_test();
			break;	  
		case 3:
			face_test();
			break;
		case 4:
			faces();
			break;
		case 5:
			ortho_cube();
			break;
		case 6:
			cube_no_clip();
			break;
		case 7:
			cube_clip_near();
			break;
		case 8:
			cube_clip_far();
			break;
		case 9:
			ortho_axes();
			break;
		case 10:
			persp_axes();
			break;
		default:
			exit(1);
		}
	count++;
	osuFlush();
	osuWaitOnEscape();
}


/******************************************************************************
Test out drawing routines.
******************************************************************************/
void mainD(int argc , char **argv)
{
  int num;
  num = atoi(argv[1]);

  if((num <0 ) || (num > 11))
	{
		fprintf(stderr, "Please call this program with a number from 1 to 12 \n");
		exit(-1);
	}

  osuBeginGraphics (xsize, ysize);

  /* inialize the matrix stack*/
  osuInitialize();
  /* select routine to execute */
  switch (num) {
    case 1:
      line_test();
      break;
    case 2:
      ortho_test();
      break;
    case 3:
      face_test();
      break;
    case 4:
      faces();
      break;
    case 5:
      ortho_cube();
      break;
    case 6:
      cube_no_clip();
      break;
    case 7:
      cube_clip_near();
      break;
    case 8:
      cube_clip_far();
      break;
    case 9:
      ortho_axes();
      break;
    case 10:
      persp_axes();
      break;
	case 11:
		all();
		break;

    default:
		 fprintf (stderr, "Please use a number from 1 to 11.\n");
         exit (-1);
  } 

  osuFlush();

  printf ("Press 'escape' to exit.\n");
  osuWaitOnEscape();
  osuEndGraphics();
}

