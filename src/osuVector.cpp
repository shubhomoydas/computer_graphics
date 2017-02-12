#include "osuVector.h"
#include <stdlib.h>
#include <math.h>

#define DEBUG

#ifdef DEBUG
#define DB_CHECK( C ) { if ( ! (C) ) { abort(); } } 
#else
#define DB_CHECK( C ) { } 
#endif

double dot(const osuVector4 &a, const osuVector4 &b);

//-----------------------

osuVector4::osuVector4(){
	data[0] = data[1] = data[2] = data[3] = 0.0;
}

osuVector4::osuVector4(double d0, double d1, double d2, double d3){
	data[0] = d0;
	data[1] = d1;
	data[2] = d2;
	data[3] = d3;
}


osuVector4::~osuVector4()
{
}

//---------------------------------------------------------------------
// Inline implementation for osuVector4


double &osuVector4::operator[](unsigned int i) {
	DB_CHECK( i<4 );
	return data[i];
}

double osuVector4::operator[](unsigned int i) const {

	DB_CHECK(i<4);
	return data[i];
}


osuVector4 &osuVector4::operator/=(double d){
	data[0] /= d;
	data[1] /= d;
	data[2] /= d;
	//keep homog coord unchanged!
	return(*this);
}

osuVector4 &osuVector4::operator-(osuVector4 &a){
	data[0] = data[0] - a[0];
	data[1] = data[1] - a[1];
	data[2] = data[2] - a[2];
	data[3] = 1.0;
	return(*this);
	}

	



double dot(const osuVector4 &a, const osuVector4 &b){
	return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]);
}

osuVector4 cross(const osuVector4 &a, const osuVector4 &b) {
	return osuVector4(a[1]*b[2] - b[1]*a[2],
					a[2]*b[0] - b[2]*a[0],
					a[0]*b[1] - b[0]*a[1], 
					1.0);
}

double mag(const osuVector4 &a)
{
	return(sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]));
}
void normalize(osuVector4 &a)
{
	register double m = mag(a);
	if(m!=0) a/=m;
}



////////////////////////////////////////////////////////////////
//-----------------------

osuVector3::osuVector3(){
	data[0] = data[1] = data[2] = 0.0;
}

osuVector3::osuVector3(double d0, double d1, double d2){
	data[0] = d0;
	data[1] = d1;
	data[2] = d2;
}


osuVector3::~osuVector3()
{
}

//---------------------------------------------------------------------
// Inline implementation for osuVector3


double &osuVector3::operator[](unsigned int i) {
	DB_CHECK( i<3 );
	return data[i];
}

double osuVector3::operator[](unsigned int i) const {

	DB_CHECK(i<3);
	return data[i];
}


osuVector3 &osuVector3::operator/=(double d){
	data[0] /= d;
	data[1] /= d;
	//keep homog coord unchanged!
	return(*this);
}

osuVector3 &osuVector3::operator-(osuVector3 &a){
	data[0] = data[0] - a[0];
	data[1] = data[1] - a[1];
	data[2] = 1.0;
	return(*this);
	}


