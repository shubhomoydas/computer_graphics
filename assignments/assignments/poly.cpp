/*

Test the polygon scan conversion routines.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "osuGraphics.h"
#include "matlib.h"


//void osuColor3f (double, double, double);
//void osuVertex2f (double, double);
//void osuBegin (int);
//void osuEnd();
void simple_triangle();
void abutting_triangles();
void color_polygons();
void wedges(double,double,double,int);
void overlapTest();

/******************************************************************************
Draw a polygon.
*****************************************************************************/
void mainP(int argc , char **argv)
{

	int num = atoi(argv[1]);

	if((num <0 ) || (num > 5))
	{
		fprintf(stderr, "Please call this program with a number from 1 to 5 \n");
		exit(-1);
	}

  osuBeginGraphics (500, 500);

  /* inialize the matrix stack*/
  osuInitialize(); setGraphicsDebug(true,false);

  /* go to selected routine */

  switch (num) {
      case 1:
      simple_triangle();
      break;
    case 2:
		osuSetWriteMode(OSU_XOR);
      abutting_triangles();
      break;
    case 3:
		osuSetWriteMode(OSU_XOR);
      color_polygons();
      break;
	 case 4:
	   osuColor3f(1.0, 1.0,1.0);
	   osuSetWriteMode(OSU_XOR);
	   wedges(0.501, 0.497, 0.4, 13); 
	   break;
	 case 5:
		//overlap
	   osuSetWriteMode(OSU_XOR);
		overlapTest();
		break;
    default:
      fprintf (stderr, "Please use a number from 1 to 5.\n");
      exit (-1);
  }

  osuFlush();
  osuWaitOnEscape();
  osuEndGraphics();
}


/******************************************************************************
Draw a white triangle.
******************************************************************************/

void simple_triangle()
{
  osuBegin (OSU_TRIANGLE);
  osuColor3f (1.0, 1.0, 1.0);
  osuVertex2f (0.25, 0.25);
  osuVertex2f (0.25, 0.75);
  osuVertex2f (0.75, 0.5);
  osuEnd ();
}


/******************************************************************************
Draw several abutting rectangles.
******************************************************************************/

void abutting_triangles()
{
  double a = 0.1; 
  double b = 0.5;
  double c = 0.9;

  //#ifdef NOTDEF
  osuBegin (OSU_TRIANGLE);
  osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (a, a);
	osuVertex2f (b, b);
	osuVertex2f (a, c);
  osuEnd();
 
  osuBegin(OSU_TRIANGLE);
  osuColor3f(1.0, 0.0, 0.0);
	osuVertex2f(b,b);
	osuVertex2f(c,c);
	osuVertex2f (a,c);
  osuEnd();


  osuBegin(OSU_TRIANGLE);
  osuColor3f(0.0, 1.0, 0.0);
	osuVertex2f(b,b);
	osuVertex2f(c,c);
	osuVertex2f (c,a);
  osuEnd();


  osuBegin(OSU_TRIANGLE);
  osuColor3f(0.0, 0.0, 1.0);
	osuVertex2f(a,a);
	osuVertex2f(c,a);
	osuVertex2f (b,b);
  osuEnd();

}



/******************************************************************************
Draw some triangles that use color interpolation.
******************************************************************************/

void color_polygons()
{
  double x0 = 0.6;
  double y0 = 0.6;
  double x1 = 0.9;
  double y1 = 0.9;
  double x2,y2;

  /* colorful triangle */

  osuBegin (OSU_TRIANGLE);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (0.1, 0.1);
	osuColor3f (0.0, 1.0, 0.0);
	osuVertex2f (0.1, 0.5);
	osuColor3f (0.0, 0.0, 1.0);
	osuVertex2f (0.5, 0.3);
  osuEnd ();

  /* colors for square */
  osuBegin (OSU_TRIANGLE);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (x0, y0);
	osuColor3f (0.0, 1.0, 0.0);
	osuVertex2f (x0, y1);
	osuColor3f (0.0, 0.0, 1.0);
	osuVertex2f (x1, y1);
  osuEnd();

  osuBegin(OSU_TRIANGLE);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (x0, y0);
	osuColor3f(0.0, 0.0, 1.0);
	osuVertex2f(x1,y1);
	osuColor3f(1.0, 1.0, 1.0);
	osuVertex2f(x1, y0);
  osuEnd ();


  x0 = 0.55;
  y0 = 0.15;
  x1 = 0.7;
  y1 = 0.3;
  x2 = 0.85;
  y2 = 0.45;

  osuBegin (OSU_TRIANGLE);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (x0, y1);
	osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (x1, y0);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (x2, y1);
  osuEnd();

  osuBegin(OSU_TRIANGLE);
	osuColor3f(1.0, 0.0, 0.0);
	osuVertex2f(x0,y1);
	osuColor3f(1.0, 0.0, 0.0);
	osuVertex2f(x2,y1);
	osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (x1, y2);
  osuEnd ();

  x0 = 0.15;
  y0 = 0.55;
  x1 = 0.3;
  y1 = 0.7;
  x2 = 0.45;
  y2 = 0.85;

  osuBegin (OSU_TRIANGLE);
	osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (x0, y1);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (x1, y0);
	osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (x2, y1);
  osuEnd();

  osuBegin(OSU_TRIANGLE);
  osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (x0, y1);
	osuColor3f (1.0, 1.0, 1.0);
	osuVertex2f (x2, y1);
	osuColor3f (1.0, 0.0, 0.0);
	osuVertex2f (x1, y2);
  osuEnd ();

}

/******************************************************************************
Draw a polygonal approximation to a circle, but draw it with pie wedges
(thin triangles radiating from the center).

Entry:
  xc,yc  - circle center
  radius - radius of circle
  steps  - number of vertices in approximation
******************************************************************************/

void wedges (double xc, double yc, double radius, int steps)
{
  int i;
  double theta;
  double x,y;
  double xold,yold;

  theta = 2 * 3.1415926535 * (0.5) / (double) steps;
  xold = xc + radius * cos(theta);
  yold = yc + radius * sin(theta);

  for (i = 1; i <= steps; i++) {

    theta = 2 * 3.1415926535 * (i + 0.5) / (double) steps;
    x = xc + radius * cos(theta);
    y = yc + radius * sin(theta);

    osuBegin (OSU_TRIANGLE);
    osuVertex2f (xc, yc);
    osuVertex2f (x, y);
    osuVertex2f (xold, yold);
    osuEnd ();

    xold = x;
    yold = y;
  }
}


void overlapTest()
{
	  osuSetWriteMode(OSU_XOR);
	osuBegin(OSU_TRIANGLE);
	  osuColor3f(1.0,1.0,1.0);
	  osuVertex2f(0.1, 0.1);
	  osuVertex2f(0.6,0.5);
	  osuVertex2f(0.1,0.9);
	osuEnd();



	osuBegin(OSU_TRIANGLE);
		osuVertex2f(0.9, 0.1);
		osuVertex2f(0.4,0.5);
	    osuVertex2f(0.9, 0.9);
	osuEnd();
}
