/*

Line clip and drawing.

*/


#include <stdio.h>
#include <math.h>
#include "osuGraphics.h"

//typedef struct osuVertex {
//  double x,y,z,w;
//  double r,g,b;
//} osuVertex;


/* the clipping window */
static double xmin = 0.0;
static double xmax = 1.0;
static double ymin = 0.0;
static double ymax = 1.0;

static void set_clip_window(double, double, double, double);
static int clip_line(double *, double *, double *, double *);





/******************************************************************************
Copy the contents of a source vertex to a destination vertex.
******************************************************************************/

void vertex_copy (osuVertex &dest, osuVertex &source)
{
  dest.x = source.x;
  dest.y = source.y;
  dest.z = source.z;

  dest.r = source.r;
  dest.g = source.g;
  dest.b = source.b;
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
  osuVertex &v0,
  osuVertex &v1,
  osuVertex &newv,
  double a,
  double b,
  double c,
  double d
)
{
  double t;
  double x0,y0,z0;
  double x1,y1,z1;
  double r0,g0,b0;
  double r1,g1,b1;
  double dx,dy,dz;

  /* shorthands */

  x0 = v0.x;
  y0 = v0.y;
  z0 = v0.z;

  x1 = v1.x;
  y1 = v1.y;
  z1 = v1.z;

  dx = x1 - x0;
  dy = y1 - y0;
  dz = z1 - z0;

  /* find parameter t saying how far between v0 and v1 the intersection is */

  t = -1.0 * (a*x0 + b*y0 + c*z0 + d) / (a*dx + b*dy + c*dz);

  /* interpolate between values in v0 and v1 for location and color */

  newv.x = x0 + t * (x1 - x0);
  newv.y = y0 + t * (y1 - y0);
  newv.z = z0 + t * (z1 - z0);

  newv.r = v0.r + t * (v1.r - v0.r);
  newv.g = v0.g + t * (v1.g - v0.g);
  newv.b = v0.b + t * (v1.b - v0.b);
}


/******************************************************************************
Draw a white line.

Entry:
  x0,y0 - first endpoint of line
  x1,y1 - second line endpoint
******************************************************************************/

void draw_line(double x0, double y0, double x1, double y1)
{
  int i;
  double x,y;
  double dx,dy;
  double xinc,yinc;
  double lenosuh;
  int result;
  int width,height;

  /* set the clipping window */
  osuGetFramebufferSize (&width, &height);
  set_clip_window (0.0f, 0.0f, width - 0.51f, height - 0.51f);

  /* clip the line in 2D */
  result = clip_line (&x0, &y0, &x1, &y1);

  /* return if line is entirely outside the clip window */
  if (result == 0)
    return;

  /* incremental line drawing */

  dx = x1 - x0;
  dy = y1 - y0;

  /* determine whether horizontal or vertical change is larger */

  if (fabs(dx) > fabs(dy))
    lenosuh = (double)fabs(dx);
  else
    lenosuh = (double)fabs(dy);

  /* special case to avoid dividing by zero */

  if (lenosuh == 0) {
    osuWritePixel ((int) floor(x0+0.5), (int) floor(y0+0.5), 255, 255, 255);
    return;
  }

  xinc = dx / lenosuh;
  yinc = dy / lenosuh;

  x = x0;
  y = y0;

  /* write "lenosuh" number of pixels along the line */

  for (i = 0; i <= lenosuh; i++) {
    osuWritePixel ((int) floor(x+0.5), (int) floor(y+0.5), 255, 255, 255);
    x += xinc;
    y += yinc;
  }
}


/******************************************************************************
Specify a clipping window.

Entry:
  x0,y0 - lower left boundaries of clipping window
  x1,y1 - upper right boundaries
******************************************************************************/

static void set_clip_window(double x0, double y0, double x1, double y1)
{
  xmin = x0;
  ymin = y0;

  xmax = x1;
  ymax = y1;
}


/******************************************************************************
Given a point P outside the window and the rise and run of a line, return
the intersection of line with window that is nearest P.

Entry:
  dx,dy - run and rise of line
  x,y   - the given point P

Exit:
  ix,iy - intersection point
  return 1 if there was a valid intersection, 0 if not
******************************************************************************/

static int clip_helper(
  double dx,
  double dy,
  double x,
  double y,
  double *ix,
  double *iy
)
{
  /* if line not vertical, check against left and right edges of window */

  if (dx != 0) {

    /* check against left edge */
    *iy = dy / dx * (xmin - x) + y;
    if (xmin > x && *iy >= ymin && *iy < ymax) {
      *ix = xmin;
      return (1);
    }

    /* check against right edge */
    *iy = dy / dx * (xmax - x) + y;
    if (xmax < x && *iy > ymin && *iy < ymax) {
      *ix = xmax;
      return (1);
    }
  }

  /* if line not horizontal, check against top and bottom edges of window */

  if (dy != 0) {

    /* check against bottom edge */
    *ix = dx / dy * (ymin - y) + x;
    if (ymin > y && *ix > xmin && *ix < xmax) {
      *iy = ymin;
      return (1);
    }

    /* check against top edge */
    *ix = dx / dy * (ymax - y) + x;
    if (ymax < y && *ix > xmin && *ix < xmax) {
      *iy = ymax;
      return (1);
    }
  }

  /* if we get here, we found no intersection */
  return (0);
}


/******************************************************************************
Clip a line segment to a pre-specified window.

Entry:
  x0,y0 - first line segment endpoint
  x1,y1 - second endpoint

Exit:
  x0,y0,x1,y1 - clipped endpoint positions
  returns 1 if segment is at least partially in window,
  returns 0 if segment is entirely outside window
******************************************************************************/

static int clip_line(double *x0, double *y0, double *x1, double *y1)
{
  int count;
  double dx,dy;
  double xx0 = *x0;
  double yy0 = *y0;
  double xx1 = *x1;
  double yy1 = *y1;

  int code04 = (xx0 < xmin) ? 1 : 0;
  int code03 = (xx0 > xmax) ? 1 : 0;
  int code02 = (yy0 < ymin) ? 1 : 0;
  int code01 = (yy0 > ymax) ? 1 : 0;

  int code14 = (xx1 < xmin) ? 1 : 0;
  int code13 = (xx1 > xmax) ? 1 : 0;
  int code12 = (yy1 < ymin) ? 1 : 0;
  int code11 = (yy1 > ymax) ? 1 : 0;

  int sum0 = code01 + code02 + code03 + code04;
  int sum1 = code11 + code12 + code13 + code14;

  /* completely inside window? */
  if (sum0 == 0 && sum1 == 0)
    return (1);

  /* check for trivial invisibility (both endpoints on wrong side of */
  /* a single side of the window) */

  if (code01 && code11 || code02 && code12 || code03 && code13 ||
      code04 && code14) {
    return (0);
  }

  /* compute run and rise */
  dx = xx1 - xx0;
  dy = yy1 - yy0;

  /* case: only x0,y0 is inside window */
  if (sum0 == 0) {
    int dummy = clip_helper (dx, dy, xx1, yy1, &xx1, &yy1);
    *x0 = xx0;
    *y0 = yy0;
    *x1 = xx1;
    *y1 = yy1;
    return (1);
  }

  /* case: only x1,y1 is inside window */
  if (sum1 == 0) {
    int dummy = clip_helper (dx, dy, xx0, yy0, &xx0, &yy0);
    *x0 = xx0;
    *y0 = yy0;
    *x1 = xx1;
    *y1 = yy1;
    return (1);
  }

  /* neither endpoint is inside the window */

  count = 0;
  count += clip_helper (dx, dy, xx0, yy0, &xx0, &yy0);
  count += clip_helper (dx, dy, xx1, yy1, &xx1, &yy1);

  *x0 = xx0;
  *y0 = yy0;
  *x1 = xx1;
  *y1 = yy1;

  if (count)
    return (1);
  else
    return (0);
}


/******************************************************************************
Clip a line segment to front and back clipping planes.  These clip planes
are along the z-axis.  If your objects are on the negative z portion of
the axis, be sure to specify negative values for "near" and "far".

Entry:
  near,far - clip planes along z-axis
  x0,y0,z0 - first line segment endpoint
  x1,y1,z1 - second endpoint

Exit:
  x0,y0,z0,x1,y1,z1 - clipped endpoint positions
  returns 1 if segment is at least partially in window,
  returns 0 if segment is entirely outside window
******************************************************************************/

int near_far_clip(double nearp, double farp, double *x0, double *y0,
				double *z0,double *x1,double *y1,double *z1)
{
  double temp;
  double fract;
  double xx0 = *x0;
  double yy0 = *y0;
  double zz0 = *z0;
  double xx1 = *x1;
  double yy1 = *y1;
  double zz1 = *z1;
  int code00,code01,code10,code11;

  /* make sure near < far */

  if (nearp > farp) {
    temp = farp;
    farp = nearp;
    nearp = temp;
  }

  /* figure out which endpoints are outside the clipping volume */

  code00 = (zz0 < nearp) ? 1 : 0;
  code01 = (zz0 > farp)  ? 1 : 0;
  code10 = (zz1 < nearp) ? 1 : 0;
  code11 = (zz1 > farp)  ? 1 : 0;

  /* return without clipping if all endpoints are inside clip volume */

  if (code00 + code01 + code10 + code11 == 0)
    return (1);  /* signals inside volume */

  /* if both endpoints are entirely out of clip volume, exit and signal this */

  if (code00 && code10 || code01 && code11)
    return (0);  /* signals outside volume */

  /* clip to near plane if necessary */

  if (code00) {
    fract = (nearp - zz0) / (zz1 - zz0);
    xx0 = xx0 + fract * (xx1 - xx0);
    yy0 = yy0 + fract * (yy1 - yy0);
    zz0 = nearp;
  }
  else if (code10) {
    fract = (nearp - zz1) / (zz0 - zz1);
    xx1 = xx1 + fract * (xx0 - xx1);
    yy1 = yy1 + fract * (yy0 - yy1);
    zz1 = nearp;
  }

  /* clip to far plane if necessary */

  if (code01) {
    fract = (farp - zz0) / (zz1 - zz0);
    xx0 = xx0 + fract * (xx1 - xx0);
    yy0 = yy0 + fract * (yy1 - yy0);
    zz0 = farp;
  }
  else if (code11) {
    fract = (farp - zz1) / (zz0 - zz1);
    xx1 = xx1 + fract * (xx0 - xx1);
    yy1 = yy1 + fract * (yy0 - yy1);
    zz1 = farp;
  }

  /* copy the clipped endpoints */

  *x0 = xx0;
  *y0 = yy0;
  *z0 = zz0;
  *x1 = xx1;
  *y1 = yy1;
  *z1 = zz1;

  /* signal that we're inside the clip volume */
  return (1);
}

/******************************************************************************
Clip a line to a plane.

Entry:
  v0, v1   - vertices of polygon to clip
  a,b,c,d - coefficients of plane equation against which to clip:
	    positive side described by ax + by + cz + d > 0 are kept

Exit:
  int status
	1 partially visible
	0 not visible at all
******************************************************************************/

int lineClip (osuVertex &v0, osuVertex &v1,double a, double b,double c,double d0, double d1)
{
  osuVertex newv0, newv1;
  
  int in0,in1;  /* are v0 or v1 in the proper half-space */
  int status = 1;

  in0 = (a * v0.x + b * v0.y + c * v0.z + d0 > 0);
  in1 = (a * v1.x + b * v1.y + c * v1.z + d1 > 0);

    if (in0 && in1) {
		//do nothing, I'll return v0 and v1
    }
    else if (!in0 && in1) {
	  create_vertex(v0, v1, newv0, a, b, c,d0);
	  vertex_copy(v0, newv0);
	  
    }
    else if (in0 && !in1) {
	  create_vertex(v0, v1, newv1, a, b, c,d1);
	  vertex_copy(v1, newv1);
    }
    else {
      /* both are not in, so we add no vertices to the clipped line*/
		status = 0;
    }
    in0 = in1;
	return(status);
  }
