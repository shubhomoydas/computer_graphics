

#include <stdio.h>
#include "osuGraphics.h"

/* forward declaration */
void draw_rectangle (int, int, int, int, int, int, int);

/******************************************************************************
Draw some rectangles.
******************************************************************************/

void mainR()
{
  osuBeginGraphics (200, 200);

  /* draw some rectangles */
  //draw_rectangle (0,0, 199, 199,  255, 0, 0);
  draw_rectangle (150, 20, 180, 150, 0, 255, 0);
  draw_rectangle (35, 150, 180, 190, 0, 0, 255);
  

  /* write the framebuffer to a file */
  //osuWriteFramebuffer ("rectangles.ppm");

  osuFlush();
  osuWaitOnEscape();
  osuEndGraphics();
}


/******************************************************************************
Draw a rectangle in a given color.

Entry:
  x0,y0 - lower left corner of rectangle
  x1,y1 - upper right corner
  r,g,b - color of rectangle
******************************************************************************/

void draw_rectangle (int x0, int y0, int x1, int y1, int r, int g, int b)
{
  int i,j;

  /* loop through y and x pixel coordinates */

  for (j = y0; j <= y1; j++)
    for (i = x0; i <= x1; i++)
      osuWritePixel (i, j, r, g, b);
}


