/*

Routine to perform polygon clipping.

*/


#include <stdio.h>
#include <math.h>
#include "osuGraphics.h"


/******************************************************************************
Copy the contents of a source vertex to a destination vertex.
******************************************************************************/

void vertex_copy (osuVertex *dest, osuVertex *source)
{
  dest->x = source->x;
  dest->y = source->y;
  dest->z = source->z;

  dest->r = source->r;
  dest->g = source->g;
  dest->b = source->b;
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
  osuVertex *v0,
  osuVertex *v1,
  osuVertex *newv,
  float a,
  float b,
  float c,
  float d
)
{
  float t;
  float x0,y0,z0;
  float x1,y1,z1;
  float r0,g0,b0;
  float r1,g1,b1;
  float dx,dy,dz;

  /* shorthands */

  x0 = v0->x;
  y0 = v0->y;
  z0 = v0->z;

  x1 = v1->x;
  y1 = v1->y;
  z1 = v1->z;

  dx = x1 - x0;
  dy = y1 - y0;
  dz = z1 - z0;

  /* find parameter t saying how far between v0 and v1 the intersection is */

  t = -1.0 * (a*x0 + b*y0 + c*z0 + d) / (a*dx + b*dy + c*dz);

  /* interpolate between values in v0 and v1 for location and color */

  newv->x = x0 + t * (x1 - x0);
  newv->y = y0 + t * (y1 - y0);
  newv->z = z0 + t * (z1 - z0);

  newv->r = v0->r + t * (v1->r - v0->r);
  newv->g = v0->g + t * (v1->g - v0->g);
  newv->b = v0->b + t * (v1->b - v0->b);
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
  osuVertex *verts,
  int count,
  osuVertex *out_verts,
  int *out_count,
  float a,
  float b,
  float c,
  float d
)
{
  int i,ii;
  int new_count = 0;
  osuVertex *v0,*v1;
  int in0,in1;  /* are v0 or v1 in the proper half-space */

  v0 = &verts[0];
  in0 = (a * v0->x + b * v0->y + c * v0->z + d > 0);

  for (i = 0; i < count; i++) {

    v0 = &verts[i];
    v1 = &verts[(i+1) % count];
    in1 = (a * v1->x + b * v1->y + c * v1->z + d > 0);

    if (in0 && in1) {
      vertex_copy (&out_verts[new_count++], v1);
    }
    else if (!in0 && in1) {
      create_vertex (v0, v1, &out_verts[new_count++], a, b, c, d);
      vertex_copy (&out_verts[new_count++], v1);
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

